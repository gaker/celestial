// Tile-compressed images read through `Image::open`. The `.fits.fz` fixtures
// in test_data/compressed were written by fpack 1.7.0. Each `.ref.fit` holds
// the pixels astropy 6.0.1 decodes from the matching `.fits.fz`, written as an
// uncompressed image. funpack decodes the same pixels, except where a test
// says otherwise.
use super::Image;
use crate::core::errors::ImageError;
use crate::fits::compression::CompressionAlgorithm;
use crate::fits::errors::FitsError;
use crate::fits::io::test_support::FitsWriter;
use crate::formats::pixel_data::PixelData;
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("test_data/compressed")
        .join(name)
}

// Floats compare by bit pattern, so NaN matches NaN with the same payload.
fn sample_bits(pixels: &PixelData) -> (&'static str, Vec<u64>) {
    match pixels {
        PixelData::U8(v) => ("U8", v.iter().map(|&x| x as u64).collect()),
        PixelData::I8(v) => ("I8", v.iter().map(|&x| x as u64).collect()),
        PixelData::U16(v) => ("U16", v.iter().map(|&x| x as u64).collect()),
        PixelData::I16(v) => ("I16", v.iter().map(|&x| x as u64).collect()),
        PixelData::U32(v) => ("U32", v.iter().map(|&x| x as u64).collect()),
        PixelData::I32(v) => ("I32", v.iter().map(|&x| x as u64).collect()),
        PixelData::U64(v) => ("U64", v.clone()),
        PixelData::I64(v) => ("I64", v.iter().map(|&x| x as u64).collect()),
        PixelData::F32(v) => ("F32", v.iter().map(|x| x.to_bits() as u64).collect()),
        PixelData::F64(v) => ("F64", v.iter().map(|x| x.to_bits()).collect()),
    }
}

fn assert_matches_reference(stem: &str) {
    let img =
        Image::open(fixture(&format!("{stem}.fits.fz"))).unwrap_or_else(|e| panic!("{stem}: {e}"));
    let reference = Image::open(fixture(&format!("{stem}.ref.fit"))).unwrap();
    assert_eq!(img.dimensions, reference.dimensions, "{stem}");
    assert_eq!(
        sample_bits(&img.pixels),
        sample_bits(&reference.pixels),
        "{stem}"
    );
}

fn open_bytes(bytes: &[u8]) -> crate::core::errors::Result<Image> {
    let tmp = tempfile::Builder::new().suffix(".fits").tempfile().unwrap();
    std::fs::write(tmp.path(), bytes).unwrap();
    Image::open(tmp.path())
}

fn fixture_bytes(name: &str) -> Vec<u8> {
    std::fs::read(fixture(name)).unwrap()
}

// Replaces the first card named `key` in the header of HDU 1. fpack writes a
// one-block primary header, so HDU 1 starts at byte 2880.
fn set_card(bytes: &mut [u8], key: &str, card: &str) {
    let name = format!("{key:<8}");
    let at = (2880..bytes.len())
        .step_by(80)
        .take_while(|&i| !bytes[i..i + 80].starts_with(b"END     "))
        .find(|&i| bytes[i..i + 8] == *name.as_bytes())
        .unwrap_or_else(|| panic!("no {key} card"));
    bytes[at..at + 80].copy_from_slice(format!("{card:<80}").as_bytes());
}

fn card(name: &str, value: &str) -> String {
    format!("{:<80}", format!("{name:<8}= {value:>20}"))
}

fn header_block(cards: &[String]) -> Vec<u8> {
    let mut bytes: Vec<u8> = cards.iter().flat_map(|c| c.bytes()).collect();
    bytes.extend_from_slice(format!("{:<80}", "END").as_bytes());
    bytes.resize(bytes.len().div_ceil(2880) * 2880, b' ');
    bytes
}

// An empty primary HDU and one compressed-image table. `columns` are
// (TTYPE, TFORM) pairs, and `data` is the table rows followed by the heap.
fn compressed_file(
    columns: &[(&str, &str)],
    row_bytes: usize,
    rows: usize,
    z_cards: &[(&str, &str)],
    data: &[u8],
) -> Vec<u8> {
    let primary = [card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "0")];
    let pcount = data.len() - row_bytes * rows;
    let mut cards = vec![
        card("XTENSION", "'BINTABLE'"),
        card("BITPIX", "8"),
        card("NAXIS", "2"),
        card("NAXIS1", &row_bytes.to_string()),
        card("NAXIS2", &rows.to_string()),
        card("PCOUNT", &pcount.to_string()),
        card("GCOUNT", "1"),
        card("TFIELDS", &columns.len().to_string()),
    ];
    for (i, (ttype, tform)) in columns.iter().enumerate() {
        cards.push(card(&format!("TTYPE{}", i + 1), &format!("'{ttype}'")));
        cards.push(card(&format!("TFORM{}", i + 1), &format!("'{tform}'")));
    }
    cards.push(card("ZIMAGE", "T"));
    cards.extend(z_cards.iter().map(|(k, v)| card(k, v)));
    let mut bytes = header_block(&primary);
    bytes.extend(header_block(&cards));
    bytes.extend_from_slice(data);
    bytes.resize(bytes.len().div_ceil(2880) * 2880, 0);
    bytes
}

fn gzip(data: &[u8]) -> Vec<u8> {
    use std::io::Write;
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    enc.write_all(data).unwrap();
    enc.finish().unwrap()
}

// One GZIP_1 tile of BITPIX 8 holding `pixels`, with the image and tile size
// given by the Z cards. The descriptor is (count, offset).
fn gzip_u8_file(z_axes: &[(&str, &str)], pixels: &[u8], descriptor: (i32, i32)) -> Vec<u8> {
    let tile = gzip(pixels);
    let mut data = Vec::new();
    data.extend_from_slice(&descriptor.0.to_be_bytes());
    data.extend_from_slice(&descriptor.1.to_be_bytes());
    data.extend_from_slice(&tile);
    let mut z = vec![("ZCMPTYPE", "'GZIP_1'"), ("ZBITPIX", "8")];
    z.extend_from_slice(z_axes);
    compressed_file(&[("COMPRESSED_DATA", "1PB")], 8, 1, &z, &data)
}

fn assert_invalid(result: crate::core::errors::Result<Image>) {
    match result {
        Err(ImageError::Fits(
            FitsError::InvalidFormat(_)
            | FitsError::InvalidKeywordValue { .. }
            | FitsError::KeywordNotFound { .. },
        )) => {}
        Err(other) => panic!("wrong error: {other}"),
        Ok(img) => panic!("decoded {:?}", img.dimensions),
    }
}

// `text` is part of the message, which tells this reader's rejection of the
// file apart from a rejection of every compressed image.
fn assert_unsupported(result: crate::core::errors::Result<Image>, text: &str) {
    match result {
        Err(ImageError::Fits(FitsError::UnsupportedCompression(m))) if m.contains(text) => {}
        Err(other) => panic!("wrong error for {text:?}: {other}"),
        Ok(img) => panic!("decoded {:?}", img.dimensions),
    }
}

fn assert_invalid_format(result: crate::core::errors::Result<Image>, text: &str) {
    match result {
        Err(ImageError::Fits(FitsError::InvalidFormat(m))) if m.contains(text) => {}
        Err(other) => panic!("wrong error for {text:?}: {other}"),
        Ok(img) => panic!("decoded {:?}", img.dimensions),
    }
}

#[test]
fn open_decodes_rice_tiles_by_ztile() {
    for stem in [
        "rice_i16_t7x5",
        "rice_i32_t7x5",
        "rice_u8_rows",
        "rice_u16_whole",
        "rice_cube_t2x3x2",
    ] {
        assert_matches_reference(stem);
    }
}

#[test]
fn open_decodes_gzip1_tiles_by_ztile() {
    assert_matches_reference("gzip1_i16_t7x5");
}

#[test]
fn open_reads_u16_offset_through_compression() {
    let img = Image::open(fixture("rice_u16_whole.fits.fz")).unwrap();
    assert!(img.pixels.as_u16().is_some());
}

#[test]
fn open_reads_first_image_of_compressed_mef() {
    assert_matches_reference("mef");
    let img = Image::open(fixture("mef.fits.fz")).unwrap();
    assert_eq!(img.get("EXTNAME").as_str(), Some("SCI"));
}

#[test]
fn open_restores_the_image_header() {
    let img = Image::open(fixture("rice_i16_t7x5.fits.fz")).unwrap();
    assert_eq!(img.get("BITPIX").as_i64(), Some(16));
    assert_eq!(img.get("NAXIS").as_i64(), Some(2));
    assert_eq!(img.get("NAXIS1").as_i64(), Some(23));
    assert_eq!(img.get("NAXIS2").as_i64(), Some(17));
    assert_eq!(img.get("OBJECT").as_str(), Some("fixture"));
    for gone in [
        "XTENSION", "ZIMAGE", "ZCMPTYPE", "ZBITPIX", "ZNAXIS1", "ZTILE1", "ZNAME1", "ZVAL1",
        "TTYPE1", "TFORM1", "TFIELDS", "PCOUNT", "EXTNAME", "CHECKSUM",
    ] {
        assert!(img.get_keyword(gone).is_none(), "{gone}");
    }
}

#[test]
fn open_unshuffles_gzip2_tiles() {
    for stem in ["gzip2_i16_t7x5", "gzip2_i32_rows", "gzip2_u8_whole"] {
        assert_matches_reference(stem);
    }
}

#[test]
fn open_uses_the_rice_block_size_keyword() {
    assert_matches_reference("rice_i16_block16");
}

#[test]
fn open_dequantizes_subtractive_dither_1() {
    for stem in ["rice_f32_d1", "rice_f32_seed9999", "rice_f32_wrap"] {
        assert_matches_reference(stem);
    }
}

// The fixtures were written with ZDITHER0 = 1. Without the card the same
// tiles decode to the same pixels, as funpack decodes them.
#[test]
fn open_dithers_with_seed_1_when_zdither0_is_absent() {
    for stem in ["rice_f32_seed1", "rice_f32_d2_seed1"] {
        assert_matches_reference(stem);
        let mut bytes = fixture_bytes(&format!("{stem}.fits.fz"));
        set_card(&mut bytes, "ZDITHER0", &" ".repeat(80));
        let img = open_bytes(&bytes).unwrap_or_else(|e| panic!("{stem}: {e}"));
        let reference = Image::open(fixture(&format!("{stem}.ref.fit"))).unwrap();
        assert_eq!(
            sample_bits(&img.pixels),
            sample_bits(&reference.pixels),
            "{stem}"
        );
    }
}

// funpack 1.7.0 decodes 63 pixels of tile 5 differently: it returns a number
// for the NaN at (7, 5) and shifts the rest of the tile. astropy matches the
// source image to within the quantization step.
#[test]
fn open_dequantizes_float64_dither_with_gzip2() {
    assert_matches_reference("gzip2_f64_d1_t7x5");
}

#[test]
fn open_dequantizes_subtractive_dither_2() {
    for stem in ["gzip1_f32_d2", "rice_f64_d2"] {
        assert_matches_reference(stem);
    }
}

#[test]
fn open_dequantizes_without_dither() {
    assert_matches_reference("rice_f32_nodither");
}

// I * ZSCALE + ZZERO nearly cancels at five pixels of this file. funpack
// rounds once and returns 7.324524e-07 there. astropy rounds the product
// first and returns 0.0. The reference is the funpack output.
#[test]
fn open_dequantizes_without_dither_with_one_rounding() {
    assert_matches_reference("rice_f32_nodither_whole");
}

#[test]
fn open_reads_lossless_gzip_floats() {
    for stem in ["gzip1_f32_lossless", "gzip2_f64_lossless"] {
        assert_matches_reference(stem);
    }
}

#[test]
fn open_reads_gzip_fallback_tiles() {
    for stem in ["rice_f32_fallback", "gzip2_f32_fallback"] {
        assert_matches_reference(stem);
    }
}

// astropy cannot read these files. The references are funpack output.
#[test]
fn open_reads_nocompress_tiles() {
    for stem in ["nocompress_i16", "nocompress_f32"] {
        assert_matches_reference(stem);
    }
}

#[test]
fn open_accepts_rice_one_as_rice_1() {
    let mut bytes = fixture_bytes("rice_i16_t7x5.fits.fz");
    set_card(&mut bytes, "ZCMPTYPE", "ZCMPTYPE= 'RICE_ONE'");
    let img = open_bytes(&bytes).unwrap();
    let reference = Image::open(fixture("rice_i16_t7x5.ref.fit")).unwrap();
    assert_eq!(img.pixels.as_i16(), reference.pixels.as_i16());
}

#[test]
fn open_rejects_compression_names_outside_the_standard() {
    for name in ["RICE", "GZIP", "rice_1", "LZ4"] {
        let mut bytes = fixture_bytes("rice_i16_t7x5.fits.fz");
        set_card(&mut bytes, "ZCMPTYPE", &format!("ZCMPTYPE= '{name}'"));
        assert_unsupported(open_bytes(&bytes), &format!("unknown ZCMPTYPE '{name}'"));
    }
}

#[test]
fn open_rejects_unknown_quantization_method() {
    let mut bytes = fixture_bytes("rice_f32_d1.fits.fz");
    set_card(&mut bytes, "ZQUANTIZ", "ZQUANTIZ= 'SUBTRACTIVE_DITHER_9'");
    assert_unsupported(open_bytes(&bytes), "ZQUANTIZ = 'SUBTRACTIVE_DITHER_9'");
}

#[test]
fn open_requires_a_dither_seed_in_range() {
    for seed in ["0", "10001"] {
        let mut bytes = fixture_bytes("rice_f32_d1.fits.fz");
        set_card(&mut bytes, "ZDITHER0", &format!("ZDITHER0= {seed:>20}"));
        assert_invalid(open_bytes(&bytes));
    }
}

// FITS 4.0 Table 37: BYTEPIX defaults to 4 when ZNAME2/ZVAL2 are absent.
#[test]
fn open_applies_the_bytepix_default() {
    let blank = " ".repeat(80);
    let mut bytes = fixture_bytes("rice_i32_t7x5.fits.fz");
    set_card(&mut bytes, "ZNAME2", &blank);
    set_card(&mut bytes, "ZVAL2", &blank);
    let img = open_bytes(&bytes).unwrap();
    let reference = Image::open(fixture("rice_i32_t7x5.ref.fit")).unwrap();
    assert_eq!(img.pixels.as_i32(), reference.pixels.as_i32());

    let mut bytes = fixture_bytes("rice_i16_t7x5.fits.fz");
    set_card(&mut bytes, "ZNAME2", &blank);
    set_card(&mut bytes, "ZVAL2", &blank);
    assert_invalid(open_bytes(&bytes));
}

#[test]
fn open_reports_hcompress_and_plio_as_unsupported() {
    for (name, algorithm) in [
        ("hcompress_i16.fits.fz", "HCOMPRESS_1"),
        ("plio_mask.fits.fz", "PLIO_1"),
    ] {
        let text = format!("{algorithm} tiles cannot be read");
        assert_unsupported(Image::open(fixture(name)), &text);
    }
}

#[test]
fn write_compressed_image_rejects_hcompress_and_plio() {
    for algorithm in [CompressionAlgorithm::HCompress, CompressionAlgorithm::Plio] {
        let tmp = tempfile::Builder::new().suffix(".fits").tempfile().unwrap();
        let data: Vec<i16> = (0..16).collect();
        let result = FitsWriter::create(tmp.path())
            .unwrap()
            .write_compressed_image(&data, &[4, 4], Some((4, 4)), algorithm, &[]);
        assert!(
            matches!(result, Err(FitsError::UnsupportedCompression(_))),
            "{algorithm:?}: {result:?}"
        );
    }
}

#[test]
fn open_decodes_a_crafted_gzip_tile() {
    let pixels: Vec<u8> = (1..=6).collect();
    let tile_len = gzip(&pixels).len() as i32;
    let axes = [
        ("ZNAXIS", "2"),
        ("ZNAXIS1", "3"),
        ("ZNAXIS2", "2"),
        ("ZTILE1", "3"),
        ("ZTILE2", "2"),
    ];
    let img = open_bytes(&gzip_u8_file(&axes, &pixels, (tile_len, 0))).unwrap();
    assert_eq!(img.pixels.as_u8(), Some(&pixels));
}

#[test]
fn open_rejects_tile_sizes_that_overflow() {
    let pixels = [1u8];
    let len = gzip(&pixels).len() as i32;
    let axes = [
        ("ZNAXIS", "2"),
        ("ZNAXIS1", "4294967296"),
        ("ZNAXIS2", "4294967296"),
        ("ZTILE1", "4294967296"),
        ("ZTILE2", "4294967296"),
    ];
    assert_invalid(open_bytes(&gzip_u8_file(&axes, &pixels, (len, 0))));
}

// A 2^40-pixel tile cannot come from a few compressed bytes, so the reader
// refuses it before allocating the image.
#[test]
fn open_rejects_a_tile_larger_than_its_bytes_can_hold() {
    let pixels = [1u8];
    let len = gzip(&pixels).len() as i32;
    let axes = [("ZNAXIS", "1"), ("ZNAXIS1", "1099511627776")];
    assert_invalid(open_bytes(&gzip_u8_file(&axes, &pixels, (len, 0))));
}

#[test]
fn open_rejects_bad_axis_and_tile_values() {
    let pixels = [1u8, 2];
    let len = gzip(&pixels).len() as i32;
    for axes in [
        [("ZNAXIS", "1"), ("ZNAXIS1", "-2"), ("ZTILE1", "2")],
        [("ZNAXIS", "1"), ("ZNAXIS1", "2"), ("ZTILE1", "0")],
        [("ZNAXIS", "1"), ("ZNAXIS1", "2"), ("ZTILE1", "-1")],
        [("ZNAXIS", "1000"), ("ZNAXIS1", "2"), ("ZTILE1", "2")],
    ] {
        assert_invalid(open_bytes(&gzip_u8_file(&axes, &pixels, (len, 0))));
    }
}

#[test]
fn open_rejects_bad_heap_descriptors() {
    let pixels = [1u8, 2];
    let len = gzip(&pixels).len() as i32;
    let axes = [("ZNAXIS", "1"), ("ZNAXIS1", "2"), ("ZTILE1", "2")];
    for descriptor in [(len, 1), (len + 1, 0), (-1, 0), (len, -1), (len, i32::MAX)] {
        assert_invalid(open_bytes(&gzip_u8_file(&axes, &pixels, descriptor)));
    }
}

#[test]
fn open_rejects_a_tile_that_inflates_past_its_size() {
    let pixels = [1u8, 2, 3];
    let len = gzip(&pixels).len() as i32;
    let axes = [("ZNAXIS", "1"), ("ZNAXIS1", "2"), ("ZTILE1", "2")];
    assert_invalid(open_bytes(&gzip_u8_file(&axes, &pixels, (len, 0))));
}

#[test]
fn open_rejects_a_row_count_that_differs_from_the_tile_count() {
    let pixels = [1u8, 2];
    let len = gzip(&pixels).len() as i32;
    let axes = [("ZNAXIS", "1"), ("ZNAXIS1", "4"), ("ZTILE1", "2")];
    assert_invalid(open_bytes(&gzip_u8_file(&axes, &pixels, (len, 0))));
}

#[test]
fn open_rejects_columns_wider_than_naxis1() {
    let z = [
        ("ZCMPTYPE", "'GZIP_1'"),
        ("ZBITPIX", "8"),
        ("ZNAXIS", "1"),
        ("ZNAXIS1", "1"),
    ];
    let data = vec![0u8; 8];
    let columns = [("COMPRESSED_DATA", "1PB"), ("ZSCALE", "1D")];
    assert_invalid(open_bytes(&compressed_file(&columns, 8, 1, &z, &data)));
}

// FITS 4.0 §7.3.5 gives a P or Q descriptor column a repeat count of 1. A
// count of 0 leaves no bytes in the row for the descriptor.
#[test]
fn open_rejects_descriptor_columns_without_one_descriptor() {
    let z = [
        ("ZCMPTYPE", "'GZIP_1'"),
        ("ZBITPIX", "8"),
        ("ZNAXIS", "2"),
        ("ZNAXIS1", "4"),
        ("ZNAXIS2", "3"),
        ("ZTILE1", "4"),
        ("ZTILE2", "3"),
    ];
    let only = [("COMPRESSED_DATA", "0PB")];
    assert_invalid(open_bytes(&compressed_file(&only, 0, 1, &z, &[])));
    for tform in ["0PB", "0QB", "2PB"] {
        let width = if tform == "2PB" { 16 } else { 0 };
        let columns = [("X", "1J"), ("COMPRESSED_DATA", tform)];
        let data = vec![0u8; 4 + width];
        assert_invalid(open_bytes(&compressed_file(
            &columns,
            4 + width,
            1,
            &z,
            &data,
        )));
    }
}

#[test]
fn open_rejects_uncompressed_tile_of_wrong_length() {
    let z = [
        ("ZCMPTYPE", "'NOCOMPRESS'"),
        ("ZBITPIX", "16"),
        ("ZNAXIS", "1"),
        ("ZNAXIS1", "3"),
    ];
    let mut data = Vec::new();
    data.extend_from_slice(&[0u8; 8]);
    data.extend_from_slice(&2i32.to_be_bytes());
    data.extend_from_slice(&0i32.to_be_bytes());
    data.extend_from_slice(&[0, 1, 0, 2]);
    let columns = [("COMPRESSED_DATA", "1PB"), ("UNCOMPRESSED_DATA", "1PI")];
    assert_invalid(open_bytes(&compressed_file(&columns, 16, 1, &z, &data)));
}

// A file cut inside its data unit fails the size check of the file scan. A
// heap cut short by PCOUNT = 0 leaves every tile descriptor past its end.
#[test]
fn open_rejects_truncated_compressed_files() {
    let bytes = fixture_bytes("rice_f32_d1.fits.fz");
    for cut in [bytes.len() - 2880, 2880 * 3 + 100] {
        assert_invalid_format(open_bytes(&bytes[..cut]), "data bytes, but the file holds");
    }
    let mut bytes = fixture_bytes("rice_f32_d1.fits.fz");
    set_card(&mut bytes, "PCOUNT", &card("PCOUNT", "0"));
    assert_invalid_format(
        open_bytes(&bytes),
        "row 1 of COMPRESSED_DATA points past the heap",
    );
}

// Two RICE_1 tiles of four i16 pixels each, [1, 7, 3, -7] and [7, 2, -7, 4].
// `zblank` adds a ZBLANK column with one value per tile, and `cards` adds
// header cards.
fn rice_i16_file(zblank: Option<[i32; 2]>, cards: &[(&str, &str)]) -> Vec<u8> {
    use crate::ricecomp::RiceCompressible;
    let tiles = [[1i16, 7, 3, -7], [7, 2, -7, 4]].map(|t| i16::compress(&t, 32).unwrap());
    let row_bytes = if zblank.is_some() { 12 } else { 8 };
    let mut data = Vec::new();
    let mut offset = 0;
    for (row, tile) in tiles.iter().enumerate() {
        data.extend_from_slice(&(tile.len() as i32).to_be_bytes());
        data.extend_from_slice(&(offset as i32).to_be_bytes());
        if let Some(values) = zblank {
            data.extend_from_slice(&values[row].to_be_bytes());
        }
        offset += tile.len();
    }
    for tile in &tiles {
        data.extend_from_slice(tile);
    }
    let mut columns = vec![("COMPRESSED_DATA", "1PB")];
    if zblank.is_some() {
        columns.push(("ZBLANK", "1J"));
    }
    let mut z = vec![
        ("ZCMPTYPE", "'RICE_1'"),
        ("ZBITPIX", "16"),
        ("ZNAXIS", "2"),
        ("ZNAXIS1", "4"),
        ("ZNAXIS2", "2"),
        ("ZTILE1", "4"),
        ("ZTILE2", "1"),
        ("ZNAME1", "'BLOCKSIZE'"),
        ("ZVAL1", "32"),
        ("ZNAME2", "'BYTEPIX'"),
        ("ZVAL2", "2"),
    ];
    z.extend_from_slice(cards);
    compressed_file(&columns, row_bytes, 2, &z, &data)
}

fn blank_card(img: &Image) -> Option<i64> {
    img.get("BLANK").as_i64()
}

#[test]
fn open_gives_integer_tiles_the_zblank_keyword_as_blank() {
    let img = open_bytes(&rice_i16_file(None, &[("ZBLANK", "7")])).unwrap();
    assert_eq!(img.pixels.as_i16(), Some(&vec![1, 7, 3, -7, 7, 2, -7, 4]));
    assert_eq!(blank_card(&img), Some(7));

    let scaled = [("ZBLANK", "7"), ("BSCALE", "2.0")];
    let img = open_bytes(&rice_i16_file(None, &scaled)).unwrap();
    let nan = f64::NAN.to_bits();
    let bits = |v: f64| v.to_bits();
    assert_eq!(
        sample_bits(&img.pixels),
        (
            "F64",
            vec![
                bits(2.0),
                nan,
                bits(6.0),
                bits(-14.0),
                nan,
                bits(4.0),
                bits(-14.0),
                bits(8.0)
            ]
        )
    );
}

#[test]
fn open_gives_integer_tiles_the_zblank_column_over_the_keyword() {
    let img = open_bytes(&rice_i16_file(Some([-7, -7]), &[("ZBLANK", "7")])).unwrap();
    assert_eq!(img.pixels.as_i16(), Some(&vec![1, 7, 3, -7, 7, 2, -7, 4]));
    assert_eq!(blank_card(&img), Some(-7));
}

// FITS 4.0 §10.1.3: BLANK gives the undefined value of an integer image.
#[test]
fn open_keeps_blank_over_zblank_on_integer_tiles() {
    let bytes = rice_i16_file(Some([-7, -7]), &[("ZBLANK", "-7"), ("BLANK", "7")]);
    assert_eq!(blank_card(&open_bytes(&bytes).unwrap()), Some(7));
}

#[test]
fn open_rejects_a_zblank_column_that_differs_between_integer_tiles() {
    let bytes = rice_i16_file(Some([7, -7]), &[]);
    assert_unsupported(open_bytes(&bytes), "ZBLANK");
}
