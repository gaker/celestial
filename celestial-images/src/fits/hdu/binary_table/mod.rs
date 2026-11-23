mod decode;
pub(crate) mod tform;

use super::table::rows::TableRows;
use super::table::TableHdu;
use super::HduTrait;
use tform::BinaryForm;

pub type BinaryTableHdu = TableHdu<BinaryForm>;
pub type BinaryTableRowIterator<'a, R> = TableRows<'a, R, BinaryForm>;

impl BinaryTableHdu {
    pub fn is_compressed_image(&self) -> bool {
        self.header().get("ZIMAGE").as_bool().unwrap_or(false)
    }

    pub fn compression_algorithm(&self) -> Option<&str> {
        if self.is_compressed_image() {
            self.header().get("ZCMPTYPE").as_str()
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::data::value::TableValue;
    use crate::fits::errors::FitsError;
    use crate::fits::hdu::HduType;
    use crate::fits::header::keywords::{set_keyword, Keyword};
    use crate::fits::header::parser::Header;
    use crate::fits::io::reader::HduInfo;
    use std::io::Cursor;

    // One table row per entry of `rows`, with the columns given as
    // (TTYPE, TFORM) pairs and `extra` cards after them.
    // A card in `extra` replaces the card of the same name.
    fn header(columns: &[(&str, &str)], row_bytes: i64, rows: i64, extra: &[Keyword]) -> Header {
        let mut cards = vec![
            Keyword::string("XTENSION", "BINTABLE"),
            Keyword::integer("BITPIX", 8),
            Keyword::integer("NAXIS", 2),
            Keyword::integer("NAXIS1", row_bytes),
            Keyword::integer("NAXIS2", rows),
            Keyword::integer("PCOUNT", 0),
            Keyword::integer("GCOUNT", 1),
            Keyword::integer("TFIELDS", columns.len() as i64),
        ];
        for (i, (name, form)) in columns.iter().enumerate() {
            cards.push(Keyword::string(format!("TTYPE{}", i + 1), name.to_string()));
            cards.push(Keyword::string(format!("TFORM{}", i + 1), form.to_string()));
        }
        for card in extra {
            set_keyword(&mut cards, card.clone());
        }
        let mut h = Header::new();
        for card in cards {
            h.add_keyword(card);
        }
        h
    }

    fn info(data_size: usize) -> HduInfo {
        HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size,
        }
    }

    fn table(
        columns: &[(&str, &str)],
        row_bytes: i64,
        rows: i64,
        extra: &[Keyword],
    ) -> BinaryTableHdu {
        let size = (row_bytes * rows) as usize;
        BinaryTableHdu::new(header(columns, row_bytes, rows, extra), info(size)).unwrap()
    }

    #[test]
    fn metadata_comes_from_the_column_keywords() {
        let extra = [
            Keyword::string("TUNIT1", "m"),
            Keyword::integer("TNULL1", -999),
            Keyword::real("TSCAL1", 2.0),
            Keyword::real("TZERO1", 100.0),
            Keyword::string("TDISP1", "I8"),
            Keyword::string("EXTNAME", "EVENTS"),
            Keyword::integer("EXTVER", 2),
        ];
        let t = table(&[("COL1", "1J"), ("COL2", "2I")], 8, 3, &extra);
        let c = t.column_info(0).unwrap();
        assert_eq!(c.name.as_deref(), Some("COL1"));
        assert_eq!(c.format, "1J");
        assert_eq!(c.unit.as_deref(), Some("m"));
        assert_eq!(c.null_value.as_deref(), Some("-999"));
        assert_eq!(c.scale, Some(2.0));
        assert_eq!(c.zero_offset, Some(100.0));
        assert_eq!(c.display_format.as_deref(), Some("I8"));
        let c = t.column_info(1).unwrap();
        assert_eq!(
            (c.index, c.unit, c.null_value, c.scale),
            (1, None, None, None)
        );
        assert!(t.column_info(2).is_err());
        assert_eq!(t.all_column_info().unwrap().len(), 2);
        assert_eq!(t.column_count().unwrap(), 2);
        assert_eq!(t.number_of_fields(), Some(2));
        assert_eq!(t.number_of_rows(), Some(3));
        assert_eq!(t.row_count(), 3);
        assert_eq!(t.extension_name(), Some("EVENTS"));
        assert_eq!(t.extension_version(), Some(2));
        assert_eq!(t.hdu_type(), HduType::BinaryTable);
    }

    #[test]
    fn names_match_exactly_first_then_without_case() {
        let t = table(
            &[("flux", "1E"), ("FLUX", "1E"), ("Time", "1D")],
            16,
            1,
            &[],
        );
        assert_eq!(t.column_by_name("FLUX").unwrap(), 1);
        assert_eq!(t.column_by_name("flux").unwrap(), 0);
        assert_eq!(t.column_by_name("TIME").unwrap(), 2);
        assert!(matches!(
            t.column_by_name("RATE"),
            Err(FitsError::InvalidFormat(_))
        ));
    }

    #[test]
    fn layout_errors_are_returned_when_the_table_is_built() {
        let bad = [
            header(&[("A", "1J")], 4, 2, &[Keyword::integer("TFIELDS", 1000)]),
            header(&[("A", "1Z")], 4, 1, &[]),
            header(&[("A", "1J")], -4, 1, &[]),
            header(&[("A", "1J")], 4, 1, &[Keyword::real("TSCAL1", f64::NAN)]),
        ];
        for h in bad {
            assert!(BinaryTableHdu::new(h, info(8)).is_err());
        }
        let missing = header(&[("A", "1J")], 4, 1, &[Keyword::integer("TFIELDS", 2)]);
        assert!(BinaryTableHdu::new(missing, info(8)).is_err());
        let short_data = header(&[("A", "1J")], 4, 3, &[]);
        assert!(BinaryTableHdu::new(short_data, info(8)).is_err());
    }

    #[test]
    fn raw_column_bytes_and_values() {
        let t = table(&[("A", "1I"), ("B", "1J")], 6, 2, &[]);
        let data = [0, 1, 0, 0, 0, 2, 0, 3, 0, 0, 0, 4];
        let mut reader = Cursor::new(data.to_vec());
        assert_eq!(
            t.read_column_raw(&mut reader, 1).unwrap(),
            [0, 0, 0, 2, 0, 0, 0, 4]
        );
        assert_eq!(
            t.get_column_by_name(&mut reader, "B").unwrap(),
            vec![TableValue::I32(2), TableValue::I32(4)]
        );
        assert!(t.get_row(&mut reader, 2).is_err());
        let mut rows = t.iter_rows(&mut reader);
        assert_eq!(rows.len(), 2);
        assert_eq!(
            rows.next().unwrap().unwrap(),
            vec![TableValue::I16(1), TableValue::I32(2)]
        );
        assert_eq!(rows.len(), 1);
    }

    #[test]
    fn floats_and_complex_values() {
        let t = table(
            &[("E", "1E"), ("D", "1D"), ("C", "1C"), ("M", "1M")],
            36,
            1,
            &[],
        );
        let mut data = Vec::new();
        data.extend_from_slice(&f32::NAN.to_be_bytes());
        data.extend_from_slice(&2.5f64.to_be_bytes());
        data.extend_from_slice(&1.0f32.to_be_bytes());
        data.extend_from_slice(&(-2.0f32).to_be_bytes());
        data.extend_from_slice(&3.0f64.to_be_bytes());
        data.extend_from_slice(&4.0f64.to_be_bytes());
        let row = t.get_row(&mut Cursor::new(data), 0).unwrap();
        assert!(matches!(row[0], TableValue::F32(v) if v.is_nan()));
        assert_eq!(
            row[1..],
            [
                TableValue::F64(2.5),
                TableValue::Complex32(1.0, -2.0),
                TableValue::Complex64(3.0, 4.0),
            ]
        );
    }

    #[test]
    fn scaled_complex_offsets_the_real_part_only() {
        let extra = [Keyword::real("TSCAL1", 2.0), Keyword::real("TZERO1", 10.0)];
        let t = table(&[("C", "1C")], 8, 1, &extra);
        let mut data = 1.0f32.to_be_bytes().to_vec();
        data.extend_from_slice(&(-3.0f32).to_be_bytes());
        let row = t.get_row(&mut Cursor::new(data), 0).unwrap();
        assert_eq!(row, vec![TableValue::Complex64(12.0, -6.0)]);
    }

    #[test]
    fn logical_bytes_other_than_t_f_are_null() {
        let t = table(&[("L", "1L")], 1, 1, &[]);
        let row = t.get_row(&mut Cursor::new(vec![b'x']), 0).unwrap();
        assert_eq!(row, vec![TableValue::Null]);
    }

    #[test]
    fn empty_and_zero_repeat_columns() {
        let t = table(&[("Z", "0J"), ("P", "0PJ"), ("A", "0A")], 0, 2, &[]);
        let row = t.get_row(&mut Cursor::new(Vec::new()), 1).unwrap();
        let empty = TableValue::Array(Vec::new());
        assert_eq!(
            row,
            vec![empty.clone(), empty, TableValue::String(String::new())]
        );
        let none = table(&[("A", "1J")], 4, 0, &[]);
        assert_eq!(
            none.get_column_values(&mut Cursor::new(Vec::new()), 0)
                .unwrap(),
            vec![]
        );
        assert_eq!(
            none.read_column_raw(&mut Cursor::new(Vec::new()), 0)
                .unwrap(),
            Vec::<u8>::new()
        );
    }

    #[test]
    fn compressed_image_accessors() {
        let zimage = [
            Keyword::logical("ZIMAGE", true),
            Keyword::string("ZCMPTYPE", "GZIP_1"),
        ];
        let t = table(&[("COMPRESSED_DATA", "1PB")], 8, 1, &zimage);
        assert!(t.is_compressed_image());
        assert_eq!(t.compression_algorithm(), Some("GZIP_1"));
        let plain = table(
            &[("A", "1J")],
            4,
            1,
            &[Keyword::string("ZCMPTYPE", "GZIP_1")],
        );
        assert!(!plain.is_compressed_image());
        assert_eq!(plain.compression_algorithm(), None);
        let off = table(&[("A", "1J")], 4, 1, &[Keyword::logical("ZIMAGE", false)]);
        assert!(!off.is_compressed_image());
    }
}
