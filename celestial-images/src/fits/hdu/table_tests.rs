use crate::fits::data::value::TableValue;
use crate::fits::errors::FitsError;
use crate::fits::hdu::Hdu;
use crate::fits::io::reader::FitsFile;
use std::io::Cursor;

fn card(key: &str, value: &str) -> String {
    format!("{key:<8}= {value:<70}")
}

fn block(cards: &[String]) -> Vec<u8> {
    let mut bytes: Vec<u8> = cards.concat().into_bytes();
    bytes.extend_from_slice(format!("{:<80}", "END").as_bytes());
    bytes.resize(bytes.len().div_ceil(2880) * 2880, b' ');
    bytes
}

// A primary HDU with no data, then one extension with these cards and data.
fn file(cards: &[(&str, String)], data: &[u8]) -> Vec<u8> {
    let primary = ["SIMPLE", "BITPIX", "NAXIS", "EXTEND"]
        .iter()
        .zip(["T", "8", "0", "T"])
        .map(|(k, v)| card(k, v))
        .collect::<Vec<_>>();
    let extension = cards.iter().map(|(k, v)| card(k, v)).collect::<Vec<_>>();
    let mut bytes = block(&primary);
    bytes.extend(block(&extension));
    bytes.extend_from_slice(data);
    bytes.resize(bytes.len().div_ceil(2880) * 2880, 0);
    bytes
}

fn quoted(s: &str) -> String {
    format!("'{s:<8}'")
}

// A BINTABLE with the given TFORMs, extra cards and data. NAXIS1 is the
// row length the caller passes, so a header can lie about it.
fn bintable(
    forms: &[&str],
    extra: &[(&str, &str)],
    row_bytes: usize,
    rows: usize,
    data: &[u8],
) -> Vec<u8> {
    let heap = data.len().saturating_sub(row_bytes * rows);
    let mut cards = vec![
        ("XTENSION", quoted("BINTABLE")),
        ("BITPIX", "8".to_string()),
        ("NAXIS", "2".to_string()),
        ("NAXIS1", row_bytes.to_string()),
        ("NAXIS2", rows.to_string()),
        ("PCOUNT", heap.to_string()),
        ("GCOUNT", "1".to_string()),
        ("TFIELDS", forms.len().to_string()),
    ];
    let names: Vec<String> = (1..=forms.len()).map(|n| format!("TFORM{n}")).collect();
    let ttypes: Vec<String> = (1..=forms.len()).map(|n| format!("TTYPE{n}")).collect();
    for (n, form) in forms.iter().enumerate() {
        cards.push((&ttypes[n], quoted(&format!("C{}", n + 1))));
        cards.push((&names[n], quoted(form)));
    }
    cards.extend(extra.iter().map(|(k, v)| (*k, v.to_string())));
    file(&cards, data)
}

fn open(bytes: &[u8], index: usize) -> crate::fits::errors::Result<Hdu> {
    FitsFile::new(Cursor::new(bytes.to_vec()))?.read_hdu(index)
}

fn binary_rows(bytes: &[u8]) -> Vec<Vec<TableValue>> {
    let Hdu::BinaryTable(table) = open(bytes, 1).unwrap() else {
        panic!("not a binary table");
    };
    let mut reader = Cursor::new(bytes.to_vec());
    let rows: Vec<Vec<TableValue>> = table
        .iter_rows(&mut reader)
        .collect::<crate::fits::errors::Result<_>>()
        .unwrap();
    // Every column read gives the same cells as the rows.
    for column in 0..table.column_count().unwrap() {
        let values = table.get_column_values(&mut reader, column).unwrap();
        let from_rows: Vec<TableValue> = rows.iter().map(|r| r[column].clone()).collect();
        assert_eq!(values, from_rows, "column {column}");
    }
    rows
}

fn arr(values: Vec<TableValue>) -> TableValue {
    TableValue::Array(values)
}

fn bits(values: &[u8]) -> TableValue {
    TableValue::Bits(values.iter().map(|&b| b == 1).collect())
}

fn u16v(v: u16) -> TableValue {
    TableValue::U16(v)
}

fn u32v(v: u32) -> TableValue {
    TableValue::U32(v)
}

fn u64v(v: u64) -> TableValue {
    TableValue::U64(v)
}

fn i8v(v: i8) -> TableValue {
    TableValue::I8(v)
}

#[test]
fn bit_column_takes_whole_bytes_and_returns_bits() {
    let data = [0b1010_0000, 0b0000_0001, 0, 0, 0, 7];
    let bytes = bintable(&["16X", "1J"], &[], 6, 1, &data);
    let rows = binary_rows(&bytes);
    let expected = bits(&[1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]);
    assert_eq!(rows, vec![vec![expected, TableValue::I32(7)]]);
}

#[test]
fn bit_column_with_a_partial_byte() {
    let bytes = bintable(&["3X", "1B"], &[], 2, 1, &[0b0110_0000, 9]);
    assert_eq!(
        binary_rows(&bytes),
        vec![vec![bits(&[0, 1, 1]), TableValue::Byte(9)]]
    );
}

#[test]
fn vector_column_returns_every_element() {
    let mut data = Vec::new();
    for v in [1.0f32, 2.0, 3.0] {
        data.extend_from_slice(&v.to_be_bytes());
    }
    let bytes = bintable(&["3E"], &[], 12, 1, &data);
    let f = TableValue::F32;
    assert_eq!(
        binary_rows(&bytes),
        vec![vec![arr(vec![f(1.0), f(2.0), f(3.0)])]]
    );
}

#[test]
fn table19_offsets_give_unsigned_and_signed_byte_values() {
    let mut data = Vec::new();
    for (i, j, k, b) in [
        (-32768i16, i32::MIN, i64::MIN, 0u8),
        (7232, 0, 0, 128),
        (32767, i32::MAX, i64::MAX, 255),
    ] {
        data.extend_from_slice(&i.to_be_bytes());
        data.extend_from_slice(&j.to_be_bytes());
        data.extend_from_slice(&k.to_be_bytes());
        data.push(b);
    }
    let extra = [
        ("TZERO1", "32768"),
        ("TZERO2", "2147483648"),
        ("TZERO3", "9223372036854775808"),
        ("TZERO4", "-128"),
    ];
    let bytes = bintable(&["I", "J", "K", "B"], &extra, 15, 3, &data);
    assert_eq!(
        binary_rows(&bytes),
        vec![
            vec![u16v(0), u32v(0), u64v(0), i8v(-128)],
            vec![u16v(40000), u32v(1 << 31), u64v(1 << 63), i8v(0)],
            vec![u16v(65535), u32v(u32::MAX), u64v(u64::MAX), i8v(127)],
        ]
    );
}

#[test]
fn other_scaling_gives_physical_f64_values() {
    let mut data = Vec::new();
    data.extend_from_slice(&3i16.to_be_bytes());
    data.extend_from_slice(&1.5f32.to_be_bytes());
    let extra = [("TSCAL1", "0.5"), ("TZERO1", "10"), ("TSCAL2", "2.0")];
    let bytes = bintable(&["I", "E"], &extra, 6, 1, &data);
    assert_eq!(
        binary_rows(&bytes),
        vec![vec![TableValue::F64(11.5), TableValue::F64(3.0)]]
    );
}

#[test]
fn integer_tnull_marks_null_cells() {
    let mut data = Vec::new();
    for v in [-999i32, 5] {
        data.extend_from_slice(&v.to_be_bytes());
    }
    let bytes = bintable(&["J"], &[("TNULL1", "-999")], 4, 2, &data);
    assert_eq!(
        binary_rows(&bytes),
        vec![vec![TableValue::Null], vec![TableValue::I32(5)]]
    );
}

#[test]
fn repeated_tzero_and_tscal_use_the_last_card() {
    let mut data = Vec::new();
    for v in [5i32, 6, 7] {
        data.extend_from_slice(&v.to_be_bytes());
    }
    let extra = [
        ("TZERO1", "0"),
        ("TZERO1", "100"),
        ("TSCAL1", "1"),
        ("TSCAL1", "2"),
    ];
    let bytes = bintable(&["J"], &extra, 4, 3, &data);
    let f = TableValue::F64;
    assert_eq!(
        binary_rows(&bytes),
        vec![vec![f(110.0)], vec![f(112.0)], vec![f(114.0)]]
    );
}

#[test]
fn column_info_reports_the_last_tscal_and_tzero() {
    let extra = [
        ("TZERO1", "0"),
        ("TZERO1", "100"),
        ("TSCAL1", "1"),
        ("TSCAL1", "2"),
    ];
    let Hdu::BinaryTable(table) = open(&bintable(&["J"], &extra, 4, 0, &[]), 1).unwrap() else {
        panic!("not a binary table");
    };
    let info = table.column_info(0).unwrap();
    assert_eq!((info.scale, info.zero_offset), (Some(2.0), Some(100.0)));
    let bytes = ascii_table(&[["  42", "   1.500", "ab"]], b' ', &extra);
    let Hdu::AsciiTable(table) = open(&bytes, 1).unwrap() else {
        panic!("not an ASCII table");
    };
    let info = table.column_info(0).unwrap();
    assert_eq!((info.scale, info.zero_offset), (Some(2.0), Some(100.0)));
}

#[test]
fn repeated_tnull_uses_the_last_card() {
    let mut data = Vec::new();
    for v in [5i32, 6, 7] {
        data.extend_from_slice(&v.to_be_bytes());
    }
    let extra = [("TNULL1", "5"), ("TNULL1", "6")];
    let bytes = bintable(&["J"], &extra, 4, 3, &data);
    let i = TableValue::I32;
    assert_eq!(
        binary_rows(&bytes),
        vec![vec![i(5)], vec![TableValue::Null], vec![i(7)]]
    );
}

#[test]
fn tnull_is_compared_before_scaling() {
    let data = (-32768i16).to_be_bytes();
    let extra = [("TZERO1", "32768"), ("TNULL1", "-32768")];
    let bytes = bintable(&["I"], &extra, 2, 1, &data);
    assert_eq!(binary_rows(&bytes), vec![vec![TableValue::Null]]);
}

#[test]
fn byte_vector_is_an_array_of_bytes() {
    let bytes = bintable(&["3B"], &[], 3, 1, &[1, 2, 3]);
    let b = TableValue::Byte;
    assert_eq!(binary_rows(&bytes), vec![vec![arr(vec![b(1), b(2), b(3)])]]);
}

#[test]
fn strings_end_at_the_first_nul_and_drop_trailing_spaces() {
    let data = b"abc\0\0\0\0\0\0\0ab  \0xyz  ";
    let bytes = bintable(&["10A"], &[], 10, 2, data);
    let s = |v: &str| TableValue::String(v.to_string());
    assert_eq!(binary_rows(&bytes), vec![vec![s("abc")], vec![s("ab")]]);
}

#[test]
fn logical_cells_are_true_false_or_null() {
    let bytes = bintable(&["L", "2L"], &[], 3, 2, b"T\0FxT\x01");
    let (l, null) = (TableValue::Logical, TableValue::Null);
    assert_eq!(
        binary_rows(&bytes),
        vec![
            vec![l(true), arr(vec![null.clone(), l(false)])],
            vec![null.clone(), arr(vec![l(true), null])],
        ]
    );
}

// Some writers store heap logicals as bytes 1 and 0. Those elements are
// null, and the other columns of the row still read.
#[test]
fn heap_logicals_stored_as_bytes_leave_the_row_readable() {
    let mut data = Vec::new();
    for (count, offset, j) in [(2i32, 0i32, 7i32), (1, 2, -3)] {
        data.extend_from_slice(&count.to_be_bytes());
        data.extend_from_slice(&offset.to_be_bytes());
        data.extend_from_slice(&j.to_be_bytes());
    }
    data.extend_from_slice(&[1, 0, b'T']);
    let bytes = bintable(&["1PL(2)", "1J"], &[], 12, 2, &data);
    let null = TableValue::Null;
    assert_eq!(
        binary_rows(&bytes),
        vec![
            vec![arr(vec![null.clone(), null]), TableValue::I32(7)],
            vec![arr(vec![TableValue::Logical(true)]), TableValue::I32(-3)],
        ]
    );
}

#[test]
fn columns_wider_than_naxis1_are_rejected_when_opened() {
    let bytes = bintable(&["1J", "1J"], &[], 4, 1, &[0; 4]);
    assert!(matches!(open(&bytes, 1), Err(FitsError::InvalidFormat(_))));
    let bytes = bintable(&["1J"], &[], 8, 1, &[0; 8]);
    assert!(matches!(open(&bytes, 1), Err(FitsError::InvalidFormat(_))));
}

// Rows of a 1PJ column, then a 4-byte gap, then a 12-byte heap. THEAP
// is the end of the gap unless `theap` gives another value.
fn heap_table(theap: Option<&str>, descriptors: &[(i32, i32)]) -> Vec<u8> {
    let gap_end = (8 * descriptors.len() + 4).to_string();
    let theap = theap.unwrap_or(&gap_end);
    let mut data = Vec::new();
    for (count, offset) in descriptors {
        data.extend_from_slice(&count.to_be_bytes());
        data.extend_from_slice(&offset.to_be_bytes());
    }
    data.extend_from_slice(&[0xEE; 4]);
    for v in [10i32, 20, 30] {
        data.extend_from_slice(&v.to_be_bytes());
    }
    bintable(
        &["1PJ(3)"],
        &[("THEAP", theap)],
        8,
        descriptors.len(),
        &data,
    )
}

#[test]
fn variable_length_arrays_are_read_from_the_heap_at_theap() {
    let bytes = heap_table(None, &[(2, 0), (1, 8), (0, 0)]);
    let j = TableValue::I32;
    assert_eq!(
        binary_rows(&bytes),
        vec![
            vec![arr(vec![j(10), j(20)])],
            vec![arr(vec![j(30)])],
            vec![arr(vec![])],
        ]
    );
}

#[test]
fn repeated_theap_uses_the_last_card() {
    let mut data = Vec::new();
    for v in [3i32, 0, 3, 6] {
        data.extend_from_slice(&v.to_be_bytes());
    }
    data.extend_from_slice(&[0x7f; 8]);
    for v in [11i16, 12, 13, 21, 22, 23] {
        data.extend_from_slice(&v.to_be_bytes());
    }
    let extra = [("THEAP", "16"), ("THEAP", "24")];
    let bytes = bintable(&["1PI(3)"], &extra, 8, 2, &data);
    let i = TableValue::I16;
    assert_eq!(
        binary_rows(&bytes),
        vec![
            vec![arr(vec![i(11), i(12), i(13)])],
            vec![arr(vec![i(21), i(22), i(23)])],
        ]
    );
}

#[test]
fn descriptors_outside_the_heap_are_errors() {
    for descriptor in [(4, 0), (1, 12), (-1, 0), (0x7fff_ffff, 0)] {
        let bytes = heap_table(None, &[descriptor]);
        let Hdu::BinaryTable(table) = open(&bytes, 1).unwrap() else {
            panic!("not a binary table");
        };
        let mut reader = Cursor::new(bytes.clone());
        assert!(table.get_row(&mut reader, 0).is_err(), "{descriptor:?}");
    }
}

#[test]
fn theap_outside_the_data_is_rejected_when_opened() {
    for theap in ["4", "40"] {
        let bytes = heap_table(Some(theap), &[(0, 0)]);
        assert!(
            matches!(open(&bytes, 1), Err(FitsError::InvalidFormat(_))),
            "{theap}"
        );
    }
}

#[test]
fn a3dtable_is_read_as_a_binary_table() {
    let mut bytes = bintable(&["1J"], &[], 4, 1, &42i32.to_be_bytes());
    let at = 2880 + 11;
    bytes[at..at + 8].copy_from_slice(b"A3DTABLE");
    let Hdu::BinaryTable(table) = open(&bytes, 1).unwrap() else {
        panic!("A3DTABLE is not read as a binary table");
    };
    let mut reader = Cursor::new(bytes.clone());
    assert_eq!(
        table.get_row(&mut reader, 0).unwrap(),
        vec![TableValue::I32(42)]
    );
}

#[test]
fn binary_table_raw_data_includes_the_heap() {
    use crate::fits::hdu::HduTrait;
    let bytes = heap_table(None, &[(2, 0), (1, 8), (0, 0)]);
    let Hdu::BinaryTable(table) = open(&bytes, 1).unwrap() else {
        panic!("not a binary table");
    };
    let raw = table
        .read_raw_data(&mut Cursor::new(bytes.clone()))
        .unwrap();
    // 3 rows of 8 bytes, a 4-byte gap and a 12-byte heap.
    assert_eq!(raw, bytes[2 * 2880..2 * 2880 + 40]);
}

// An IMAGE extension of two i16 pixels with these PCOUNT and GCOUNT, and a
// data unit of `values`.
fn image_extension(pcount: &str, gcount: &str, values: &[i16]) -> Vec<u8> {
    let cards = [
        ("XTENSION", quoted("IMAGE")),
        ("BITPIX", "16".to_string()),
        ("NAXIS", "1".to_string()),
        ("NAXIS1", "2".to_string()),
        ("PCOUNT", pcount.to_string()),
        ("GCOUNT", gcount.to_string()),
    ];
    let data: Vec<u8> = values.iter().flat_map(|v| v.to_be_bytes()).collect();
    file(&cards, &data)
}

#[test]
fn image_data_is_the_array_without_pcount_or_gcount_values() {
    let bytes = image_extension("4", "2", &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
    let mut fits = FitsFile::new(Cursor::new(bytes)).unwrap();
    assert_eq!(fits.read_hdu_with_data::<i16>(1).unwrap().1, vec![1, 2]);
}

#[test]
fn image_data_with_a_gcount_of_zero_is_an_error() {
    let bytes = image_extension("0", "0", &[]);
    let mut fits = FitsFile::new(Cursor::new(bytes)).unwrap();
    let error = fits
        .read_hdu_with_data::<i16>(1)
        .err()
        .map(|e| e.to_string());
    assert_eq!(
        error.as_deref(),
        Some("Invalid FITS format: HDU 1 needs 4 data bytes, but its data unit holds 0")
    );
}

// An ASCII table whose fields start at TBCOL 1, 9 and 25 in 32-byte rows.
// Each row gives the text of the three fields. Bytes outside the fields are
// `gap`, which the reader must skip.
fn ascii_table(rows: &[[&str; 3]], gap: u8, extra: &[(&str, &str)]) -> Vec<u8> {
    let mut cards = vec![
        ("XTENSION", quoted("TABLE")),
        ("BITPIX", "8".to_string()),
        ("NAXIS", "2".to_string()),
        ("NAXIS1", "32".to_string()),
        ("NAXIS2", rows.len().to_string()),
        ("PCOUNT", "0".to_string()),
        ("GCOUNT", "1".to_string()),
        ("TFIELDS", "3".to_string()),
        ("TTYPE1", quoted("N")),
        ("TBCOL1", "1".to_string()),
        ("TFORM1", quoted("I4")),
        ("TTYPE2", quoted("X")),
        ("TBCOL2", "9".to_string()),
        ("TFORM2", quoted("F8.3")),
        ("TTYPE3", quoted("S")),
        ("TBCOL3", "25".to_string()),
        ("TFORM3", quoted("A2")),
    ];
    cards.extend(extra.iter().map(|(k, v)| (*k, v.to_string())));
    let mut data = Vec::new();
    for fields in rows {
        let mut row = vec![gap; 32];
        for (field, (start, width)) in fields.iter().zip([(0, 4), (8, 8), (24, 2)]) {
            row[start..start + width].copy_from_slice(format!("{field:<width$}").as_bytes());
        }
        data.extend(row);
    }
    file(&cards, &data)
}

fn ascii_rows(bytes: &[u8]) -> Vec<Vec<TableValue>> {
    let Hdu::AsciiTable(table) = open(bytes, 1).unwrap() else {
        panic!("not an ASCII table");
    };
    let mut reader = Cursor::new(bytes.to_vec());
    let rows: Vec<Vec<TableValue>> = table
        .iter_rows(&mut reader)
        .collect::<crate::fits::errors::Result<_>>()
        .unwrap();
    for column in 0..3 {
        let values = table.get_column_values(&mut reader, column).unwrap();
        let from_rows: Vec<TableValue> = rows.iter().map(|r| r[column].clone()).collect();
        assert_eq!(values, from_rows, "column {column}");
    }
    rows
}

#[test]
fn repeated_tbcol_uses_the_last_card() {
    let bytes = ascii_table(&[["  42", "   1.500", "ab"]], b' ', &[("TBCOL1", "9")]);
    assert_eq!(ascii_rows(&bytes)[0][0], TableValue::I64(1));
}

#[test]
fn ascii_fields_are_read_at_tbcol() {
    let bytes = ascii_table(&[["  42", "   1.500", "ab"]], b'x', &[]);
    assert_eq!(
        ascii_rows(&bytes),
        vec![vec![
            TableValue::I64(42),
            TableValue::F64(1.5),
            TableValue::String("ab".to_string()),
        ]]
    );
}

#[test]
fn ascii_tnull_text_marks_null_fields() {
    let bytes = ascii_table(
        &[["-999", "   2.5", "--"], ["   7", "    N/A", "ok"]],
        b' ',
        &[
            ("TNULL1", "'-999'"),
            ("TNULL2", "'    N/A'"),
            ("TNULL3", "'--'"),
        ],
    );
    let s = |v: &str| TableValue::String(v.to_string());
    assert_eq!(
        ascii_rows(&bytes),
        vec![
            vec![TableValue::Null, TableValue::F64(2.5), TableValue::Null],
            vec![TableValue::I64(7), TableValue::Null, s("ok")],
        ]
    );
}

fn groups_file() -> Vec<u8> {
    let cards = [
        ("SIMPLE", "T"),
        ("BITPIX", "16"),
        ("NAXIS", "3"),
        ("NAXIS1", "0"),
        ("NAXIS2", "2"),
        ("NAXIS3", "1"),
        ("GROUPS", "T"),
        ("PCOUNT", "1"),
        ("GCOUNT", "2"),
        ("PSCAL1", "0.5"),
        ("PZERO1", "100.0"),
        ("BZERO", "1000"),
    ];
    let mut bytes = block(&cards.iter().map(|(k, v)| card(k, v)).collect::<Vec<_>>());
    for v in [4i16, 1, 2, -4, 3, 4] {
        bytes.extend_from_slice(&v.to_be_bytes());
    }
    bytes.resize(2 * 2880, 0);
    bytes
}

#[test]
fn random_groups_primary_is_not_read_as_an_image() {
    let bytes = groups_file();
    assert!(matches!(open(&bytes, 0).unwrap(), Hdu::RandomGroups(_)));
    let mut fits = FitsFile::new(Cursor::new(bytes)).unwrap();
    assert!(fits.read_hdu_with_data::<i16>(0).is_err());
    assert!(fits.get_image_info(0).is_err());
}

#[test]
fn random_groups_file_is_not_an_image() {
    let tmp = tempfile::Builder::new().suffix(".fits").tempfile().unwrap();
    std::fs::write(tmp.path(), groups_file()).unwrap();
    assert!(crate::formats::image::Image::open(tmp.path()).is_err());
}

#[test]
fn random_groups_file_error_names_random_groups() {
    let tmp = tempfile::Builder::new().suffix(".fits").tempfile().unwrap();
    std::fs::write(tmp.path(), groups_file()).unwrap();
    let error = crate::formats::image::Image::open(tmp.path())
        .err()
        .map(|e| e.to_string());
    assert_eq!(
        error.as_deref(),
        Some("FITS error: Invalid FITS format: the primary HDU holds random groups, not an image")
    );
}

#[test]
fn random_groups_are_read_with_their_parameter_and_array_scaling() {
    let bytes = groups_file();
    let Hdu::RandomGroups(groups) = open(&bytes, 0).unwrap() else {
        panic!("not random groups");
    };
    let read = groups.read_groups(&mut Cursor::new(bytes)).unwrap();
    assert_eq!(read.len(), 2);
    assert_eq!(read[0].parameters, vec![102.0]);
    assert_eq!(read[0].array, vec![1001.0, 1002.0]);
    assert_eq!(read[1].parameters, vec![98.0]);
    assert_eq!(read[1].array, vec![1003.0, 1004.0]);
}

#[test]
fn random_groups_use_the_last_of_repeated_scaling_and_blank_cards() {
    let cards = [
        ("SIMPLE", "T"),
        ("BITPIX", "16"),
        ("NAXIS", "3"),
        ("NAXIS1", "0"),
        ("NAXIS2", "2"),
        ("NAXIS3", "1"),
        ("GROUPS", "T"),
        ("PCOUNT", "1"),
        ("GCOUNT", "1"),
        ("PZERO1", "0"),
        ("PZERO1", "50"),
        ("BZERO", "0"),
        ("BZERO", "100"),
        ("BLANK", "1"),
        ("BLANK", "3"),
    ];
    let mut bytes = block(&cards.iter().map(|(k, v)| card(k, v)).collect::<Vec<_>>());
    for v in [1i16, 1, 3] {
        bytes.extend_from_slice(&v.to_be_bytes());
    }
    bytes.resize(2 * 2880, 0);
    let Hdu::RandomGroups(groups) = open(&bytes, 0).unwrap() else {
        panic!("not random groups");
    };
    let read = groups.read_groups(&mut Cursor::new(bytes)).unwrap();
    assert_eq!(read[0].parameters, vec![51.0]);
    assert_eq!(read[0].array[0], 101.0);
    assert!(read[0].array[1].is_nan());
}

#[test]
fn random_groups_with_a_gcount_of_zero_read_as_none() {
    let cards = [
        ("SIMPLE", "T"),
        ("BITPIX", "16"),
        ("NAXIS", "2"),
        ("NAXIS1", "0"),
        ("NAXIS2", "1"),
        ("GROUPS", "T"),
        ("PCOUNT", "4294967296"),
        ("GCOUNT", "0"),
    ];
    let bytes = block(&cards.iter().map(|(k, v)| card(k, v)).collect::<Vec<_>>());
    let Hdu::RandomGroups(groups) = open(&bytes, 0).unwrap() else {
        panic!("not random groups");
    };
    assert_eq!(groups.read_groups(&mut Cursor::new(bytes)).unwrap(), vec![]);
}

#[test]
fn random_groups_raw_data_is_the_whole_data_unit() {
    use crate::fits::hdu::HduTrait;
    let bytes = groups_file();
    let Hdu::RandomGroups(groups) = open(&bytes, 0).unwrap() else {
        panic!("not random groups");
    };
    let stored: Vec<i16> = groups.read_data(&mut Cursor::new(bytes)).unwrap();
    assert_eq!(stored, vec![4, 1, 2, -4, 3, 4]);
}
