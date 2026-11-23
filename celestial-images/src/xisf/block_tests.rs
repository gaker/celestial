// Data block checksums (XISF 1.0 section 10.5), inline and embedded block
// locations (section 10.3), and the bounds the writer gives a float image
// (section 11.5.1). Expected digests and encodings come from Python's
// hashlib, base64, binascii and zlib.

use crate::xisf::reader::XisfFile;
use crate::xisf::writer::XisfWriter;
use std::io::Cursor;

const HEADER_LEN: usize = 4096;
const DATA_START: usize = 16 + HEADER_LEN;

const SHA1_1234: &str = "12dada1fff4d4787ade3333147202c3b443e376f";
const SHA256_1234: &str = "9f64a747e1b97f131fabb6b447296c9b6f0201e79fb3c5356e6c77e89b6a806a";
const SHA512_1234: &str = "a7c976db1723adb41274178dc82e9b777941ab201c69de61d0f2bc6d27a3598f594fa748e50d88d3c2bf1e2c2e72c3cfef78c3c6d4afa90391f7e33ababca48e";

// zlib.compress(bytes([7] * 16), 9), and the SHA-1 digests of it and of the
// 16 uncompressed bytes.
const ZLIB_SEVENS: [u8; 11] = [120, 218, 99, 103, 71, 5, 0, 3, 200, 0, 113];
const ZLIB_SEVENS_BASE64: &str = "eNpjZ0cFAAPIAHE=";
const SHA1_ZLIB_SEVENS: &str = "b5583aefce62c72ea87ef895826a3597476e6fce";
const SHA1_SEVENS: &str = "7ee8b9b443515808dc6e2fcdd2606cf8677a3b3d";

fn xisf(images: &str, attached: &[u8]) -> Vec<u8> {
    let xml =
        format!(r#"<?xml version="1.0" encoding="UTF-8"?><xisf version="1.0">{images}</xisf>"#);
    assert!(xml.len() <= HEADER_LEN);
    let mut file = b"XISF0100".to_vec();
    file.extend_from_slice(&(HEADER_LEN as u32).to_le_bytes());
    file.extend_from_slice(&[0; 4]);
    file.extend_from_slice(xml.as_bytes());
    file.resize(DATA_START, 0);
    file.extend_from_slice(attached);
    file
}

// An Image element whose block starts `offset` bytes into the attached data.
fn attached(attrs: &str, offset: usize, size: usize) -> String {
    let position = DATA_START + offset;
    format!(r#"<Image {attrs} location="attachment:{position}:{size}"/>"#)
}

fn open(images: &str, attached: &[u8]) -> XisfFile<Cursor<Vec<u8>>> {
    XisfFile::new(Cursor::new(xisf(images, attached))).unwrap()
}

fn read_error(file: &mut XisfFile<Cursor<Vec<u8>>>, index: usize) -> String {
    file.read_image_data_raw(index).unwrap_err().to_string()
}

const GRAY_2X2: &str = r#"geometry="2:2:1" sampleFormat="UInt8""#;

#[test]
fn checksums_that_match_are_accepted_under_both_algorithm_names() {
    for checksum in [
        format!("sha-1:{SHA1_1234}"),
        format!("sha1:{SHA1_1234}"),
        format!("sha-256:{SHA256_1234}"),
        format!("sha256:{SHA256_1234}"),
        format!("sha-512:{SHA512_1234}"),
        format!("sha512:{SHA512_1234}"),
    ] {
        let attrs = format!(r#"{GRAY_2X2} checksum="{checksum}""#);
        let mut file = open(&attached(&attrs, 0, 4), &[1, 2, 3, 4]);
        assert_eq!(file.read_image_data_raw(0).unwrap(), vec![1, 2, 3, 4]);
    }
}

#[test]
fn a_checksum_that_does_not_match_is_an_error() {
    let cases = [
        ("sha-1", SHA1_1234),
        ("sha-256", SHA256_1234),
        ("sha-512", SHA512_1234),
    ];
    for (algorithm, digest) in cases {
        let attrs = format!(r#"{GRAY_2X2} checksum="{algorithm}:{digest}""#);
        let mut file = open(&attached(&attrs, 0, 4), &[1, 2, 3, 5]);
        let message = read_error(&mut file, 0);
        let expected_start =
            format!("Checksum mismatch: the {algorithm} digest of the data block is ");
        assert!(message.starts_with(&expected_start), "{message}");
        assert!(
            message.ends_with(&format!(", the checksum attribute gives {digest}")),
            "{message}"
        );
    }
}

#[test]
fn a_checksum_mismatch_leaves_the_other_images_readable() {
    let bad = format!(r#"{GRAY_2X2} checksum="sha1:{SHA1_SEVENS}""#);
    let images = format!("{}{}", attached(&bad, 0, 4), attached(GRAY_2X2, 4, 4));
    let mut file = open(&images, &[1, 2, 3, 4, 5, 6, 7, 8]);
    assert_eq!(
        read_error(&mut file, 0),
        format!(
            "Checksum mismatch: the sha-1 digest of the data block is {SHA1_1234}, \
             the checksum attribute gives {SHA1_SEVENS}"
        )
    );
    assert_eq!(file.read_image_data_raw(1).unwrap(), vec![5, 6, 7, 8]);
}

// XISF 1.0 section 10.5 makes SHA3-256 and SHA3-512 optional for decoders.
#[test]
fn sha3_and_unknown_algorithms_are_reported_as_unsupported() {
    for algorithm in ["sha3-256", "sha3-512", "md5"] {
        let attrs = format!(r#"{GRAY_2X2} checksum="{algorithm}:00""#);
        let images = format!("{}{}", attached(&attrs, 0, 4), attached(GRAY_2X2, 0, 4));
        let mut file = open(&images, &[1, 2, 3, 4]);
        assert_eq!(
            read_error(&mut file, 0),
            format!("Unsupported checksum algorithm: {algorithm}")
        );
        assert_eq!(file.read_image_data_raw(1).unwrap(), vec![1, 2, 3, 4]);
    }
}

#[test]
fn malformed_checksum_attributes_are_errors() {
    let short = &SHA1_1234[..38];
    let cases = [
        (
            format!("sha1:{short}"),
            "a sha-1 digest is 40 hexadecimal digits",
        ),
        (
            format!("sha1:{SHA1_1234}00"),
            "a sha-1 digest is 40 hexadecimal digits",
        ),
        (
            format!("sha1:{short}zz"),
            "a sha-1 digest is 40 hexadecimal digits",
        ),
        (
            format!("sha1: {short}0"),
            "a sha-1 digest is 40 hexadecimal digits",
        ),
        (
            "sha1:".to_string(),
            "a sha-1 digest is 40 hexadecimal digits",
        ),
        (SHA1_1234.to_string(), "it is not algorithm:digest"),
    ];
    for (checksum, problem) in cases {
        let attrs = format!(r#"{GRAY_2X2} checksum="{checksum}""#);
        let mut file = open(&attached(&attrs, 0, 4), &[1, 2, 3, 4]);
        assert_eq!(
            read_error(&mut file, 0),
            format!("Invalid XISF format: Invalid checksum attribute \"{checksum}\": {problem}")
        );
    }
}

#[test]
fn uppercase_digest_digits_are_accepted() {
    let attrs = format!(r#"{GRAY_2X2} checksum="sha1:{}""#, SHA1_1234.to_uppercase());
    let mut file = open(&attached(&attrs, 0, 4), &[1, 2, 3, 4]);
    assert_eq!(file.read_image_data_raw(0).unwrap(), vec![1, 2, 3, 4]);
}

// XISF 1.0 section 10.6.2: the digest of a compressed block is computed for
// the compressed bytes, and a block that fails the check is not decompressed.
#[test]
fn a_compressed_block_checksum_covers_the_stored_bytes() {
    let image = |digest: &str| {
        let attrs = format!(
            r#"geometry="4:4:1" sampleFormat="UInt8" compression="zlib:16" checksum="sha1:{digest}""#
        );
        attached(&attrs, 0, ZLIB_SEVENS.len())
    };
    let mut file = open(&image(SHA1_ZLIB_SEVENS), &ZLIB_SEVENS);
    assert_eq!(file.read_image_data_raw(0).unwrap(), vec![7; 16]);

    let mut file = open(&image(SHA1_SEVENS), &ZLIB_SEVENS);
    assert!(read_error(&mut file, 0).starts_with("Checksum mismatch"));

    let mut corrupt = ZLIB_SEVENS;
    corrupt[2] ^= 0xFF;
    let mut file = open(&image(SHA1_ZLIB_SEVENS), &corrupt);
    assert!(read_error(&mut file, 0).starts_with("Checksum mismatch"));
}

fn inline(attrs: &str, encoding: &str, text: &str) -> String {
    format!(r#"<Image {attrs} location="inline:{encoding}">{text}</Image>"#)
}

#[test]
fn inline_base64_and_hex_images_are_read() {
    let cases = [
        ("base64", "AQIDBA=="),
        ("base64", "\n  AQID\n  BA==\n"),
        ("hex", "01020304"),
        ("hex", " 01 02\n03\t04 "),
        ("hex", "0102030A"),
        ("hex", "0102030a"),
    ];
    for (encoding, text) in cases {
        let mut file = open(&inline(GRAY_2X2, encoding, text), &[]);
        let expected = if text.ends_with(['A', 'a']) {
            vec![1, 2, 3, 10]
        } else {
            vec![1, 2, 3, 4]
        };
        assert_eq!(file.read_image_data_raw(0).unwrap(), expected, "{text}");
    }
}

#[test]
fn inline_uint16_image_is_little_endian() {
    let attrs = r#"geometry="2:2:1" sampleFormat="UInt16""#;
    let mut file = open(&inline(attrs, "base64", "AgGwoP//AAA="), &[]);
    let samples: Vec<u16> = file.read_image_data_typed(0).unwrap();
    assert_eq!(samples, vec![0x0102, 0xA0B0, 0xFFFF, 0]);
}

// XISF 1.0 section 10.5: the digest of an inline or embedded block is
// computed for the decoded bytes, not for the text.
#[test]
fn inline_block_checksum_covers_the_decoded_bytes() {
    let attrs = format!(r#"{GRAY_2X2} checksum="sha1:{SHA1_1234}""#);
    let mut file = open(&inline(&attrs, "base64", "AQIDBA=="), &[]);
    assert_eq!(file.read_image_data_raw(0).unwrap(), vec![1, 2, 3, 4]);

    // SHA-1 of the eight ASCII characters "AQIDBA==".
    let text_digest = "ff2028b36472f592acd8474a32363f176bf39c52";
    let attrs = format!(r#"{GRAY_2X2} checksum="sha1:{text_digest}""#);
    let mut file = open(&inline(&attrs, "base64", "AQIDBA=="), &[]);
    assert!(read_error(&mut file, 0).starts_with("Checksum mismatch"));
}

#[test]
fn invalid_inline_text_is_an_error_for_that_image_only() {
    let cases = [
        ("base64", "AQID!A==", "Invalid Base64 character '!'"),
        (
            "base64",
            "AQIDBA=",
            "Base64 text has 7 characters, which is not a multiple of 4",
        ),
        (
            "base64",
            "AQID",
            "The stored data block is 3 bytes; the image geometry needs 4",
        ),
        ("base64", "AQ==AQ==", "Invalid Base64 character '='"),
        ("base64", "A===AAAA", "Invalid Base64 character '='"),
        (
            "base64",
            "AAAA====",
            "Base64 text ends with 4 '=' characters; at most 2 are allowed",
        ),
        ("hex", "0102030g", "Invalid hexadecimal character 'g'"),
        (
            "hex",
            "0102030",
            "Hexadecimal text has an odd number of digits (7)",
        ),
        ("hex", "01020304é", "Invalid hexadecimal character 'é'"),
    ];
    for (encoding, text, message) in cases {
        let images = format!(
            "{}{}",
            inline(GRAY_2X2, encoding, text),
            attached(GRAY_2X2, 0, 4)
        );
        let mut file = open(&images, &[1, 2, 3, 4]);
        assert_eq!(
            read_error(&mut file, 0),
            format!("Invalid XISF format: {message}"),
            "{text}"
        );
        assert_eq!(file.read_image_data_raw(1).unwrap(), vec![1, 2, 3, 4]);
    }
}

// XISF 1.0 section 10.3: an element with an inline block has no child
// elements.
#[test]
fn inline_text_with_markup_is_an_error() {
    let text = r#"AQID<FITSKeyword name="A" value="1"/>BA=="#;
    let mut file = open(&inline(GRAY_2X2, "base64", text), &[]);
    assert_eq!(
        read_error(&mut file, 0),
        "Invalid XISF format: The text of an inline data block holds XML markup"
    );
}

#[test]
fn unsupported_and_malformed_locations_are_errors_for_that_image_only() {
    let cases = [
        ("inline:base32", "Unknown data block encoding \"base32\""),
        (
            "url(http://example.com/a.xisb)",
            "Invalid location format: url(http://example.com/a.xisb)",
        ),
        (
            "path(/data/a.xisb)",
            "Invalid location format: path(/data/a.xisb)",
        ),
        ("attachment:x:4", "Invalid offset: x"),
    ];
    for (location, message) in cases {
        let images = format!(
            r#"<Image {GRAY_2X2} location="{location}"/>{}"#,
            attached(GRAY_2X2, 0, 4)
        );
        let mut file = open(&images, &[1, 2, 3, 4]);
        assert_eq!(file.num_images(), 2);
        assert_eq!(
            read_error(&mut file, 0),
            format!("Invalid XISF format: {message}")
        );
        assert_eq!(file.read_image_data_raw(1).unwrap(), vec![1, 2, 3, 4]);
    }
}

fn embedded(attrs: &str, children: &str) -> String {
    format!(r#"<Image {attrs} location="embedded">{children}</Image>"#)
}

#[test]
fn embedded_base64_and_hex_images_are_read() {
    for data in [
        r#"<Data encoding="base64">AQIDBA==</Data>"#,
        r#"<Data encoding="hex">
            01020304
        </Data>"#,
    ] {
        let children = format!(
            r#"<Property id="a" type="Int32" value="1"/>{data}<FITSKeyword name="OBJECT" value="'M31'" comment=""/>"#
        );
        let mut file = open(&embedded(GRAY_2X2, &children), &[]);
        assert_eq!(file.read_image_data_raw(0).unwrap(), vec![1, 2, 3, 4]);
        assert_eq!(file.get("OBJECT").as_str(), Some("M31"));
    }
}

// XISF 1.0 section 10.6: the compression attribute of an embedded block is
// on its Data element.
#[test]
fn embedded_block_compression_is_read_from_the_data_element() {
    let attrs =
        format!(r#"geometry="4:4:1" sampleFormat="UInt8" checksum="sha1:{SHA1_ZLIB_SEVENS}""#);
    let data =
        format!(r#"<Data encoding="base64" compression="zlib:16">{ZLIB_SEVENS_BASE64}</Data>"#);
    let mut file = open(&embedded(&attrs, &data), &[]);
    assert_eq!(file.read_image_data_raw(0).unwrap(), vec![7; 16]);
}

// Only a Data child of the Image element holds the image block. A Data
// element inside a Property holds the property's block.
#[test]
fn a_data_element_inside_a_property_is_not_the_image_block() {
    let children = concat!(
        r#"<Property id="v" type="F64Vector" length="1" location="embedded">"#,
        r#"<Data encoding="base64">AAAAAAAA8D8=</Data></Property>"#,
        r#"<Data encoding="hex">01020304</Data>"#,
    );
    let mut file = open(&embedded(GRAY_2X2, children), &[]);
    assert_eq!(file.read_image_data_raw(0).unwrap(), vec![1, 2, 3, 4]);
}

#[test]
fn embedded_block_problems_are_errors_for_that_image_only() {
    let cases = [
        ("", "An embedded data block has no Data element"),
        (
            r#"<Data>AQIDBA==</Data>"#,
            "The Data element of an embedded data block has no encoding attribute",
        ),
        (
            r#"<Data encoding="base32">AQIDBA==</Data>"#,
            "Unknown data block encoding \"base32\"",
        ),
        (
            r#"<Data encoding="base64">AQ<b/>IDBA==</Data>"#,
            "The text of an embedded data block holds XML markup",
        ),
    ];
    for (children, message) in cases {
        let images = format!(
            "{}{}",
            embedded(GRAY_2X2, children),
            attached(GRAY_2X2, 0, 4)
        );
        let mut file = open(&images, &[1, 2, 3, 4]);
        assert_eq!(
            read_error(&mut file, 0),
            format!("Invalid XISF format: {message}"),
            "{children}"
        );
        assert_eq!(file.read_image_data_raw(1).unwrap(), vec![1, 2, 3, 4]);
    }
}

#[test]
fn an_empty_embedded_image_element_is_an_error() {
    let images = format!(r#"<Image {GRAY_2X2} location="embedded"/>"#);
    let mut file = open(&images, &[]);
    assert_eq!(
        read_error(&mut file, 0),
        "Invalid XISF format: An embedded data block has no Data element"
    );
}

fn written_bounds_f32(data: &[f32]) -> (f64, f64) {
    let mut writer = XisfWriter::new(Cursor::new(Vec::new()));
    writer.add_image(data, data.len(), 1, 1).unwrap();
    let file = XisfFile::new(Cursor::new(writer.write_to_vec().unwrap())).unwrap();
    file.image_info(0).unwrap().bounds
}

fn written_bounds_f64(data: &[f64]) -> (f64, f64) {
    let mut writer = XisfWriter::new(Cursor::new(Vec::new()));
    writer.add_image(data, data.len(), 1, 1).unwrap();
    let file = XisfFile::new(Cursor::new(writer.write_to_vec().unwrap())).unwrap();
    file.image_info(0).unwrap().bounds
}

// XISF 1.0 section 11.5.1 gives bounds as a black point and a white point,
// so a constant image is widened to include [0, 1].
#[test]
fn a_constant_float_image_writes_bounds_that_include_0_and_1() {
    assert_eq!(written_bounds_f32(&[5.0; 4]), (0.0, 5.0));
    assert_eq!(written_bounds_f32(&[0.5; 4]), (0.0, 1.0));
    assert_eq!(written_bounds_f32(&[-3.0, f32::NAN, -3.0]), (-3.0, 1.0));
    assert_eq!(written_bounds_f64(&[1.0; 4]), (0.0, 1.0));
    assert_eq!(written_bounds_f64(&[-0.25; 2]), (-0.25, 1.0));
    assert_eq!(written_bounds_f64(&[f64::NAN, f64::INFINITY]), (0.0, 1.0));
}

#[test]
fn a_float_image_with_a_data_range_keeps_its_range_as_bounds() {
    assert_eq!(written_bounds_f32(&[5.0, 6.0]), (5.0, 6.0));
    assert_eq!(written_bounds_f64(&[-2.0, f64::NAN, 0.5]), (-2.0, 0.5));
}

// An Image element without a location is listed and reported as unreadable,
// so image indices and property indices stay aligned.
#[test]
fn an_image_without_a_location_is_listed_as_unreadable() {
    let unlocated =
        format!(r#"<Image {GRAY_2X2}><Property id="x" type="Int32" value="3"/></Image>"#);
    let located = format!(
        r#"<Image {GRAY_2X2} location="attachment:{DATA_START}:4"><Property id="y" type="Int32" value="4"/></Image>"#
    );
    let images = unlocated + &located;
    let mut file = open(&images, &[1, 2, 3, 4]);
    assert_eq!(file.num_images(), 2);
    assert_eq!(
        read_error(&mut file, 0),
        "Invalid XISF format: The Image element has no location attribute"
    );
    assert_eq!(file.read_image_data_raw(1).unwrap(), vec![1, 2, 3, 4]);
    let ids = |props: Vec<crate::xisf::header::property::XisfProperty>| {
        props.into_iter().map(|p| p.id).collect::<Vec<String>>()
    };
    assert_eq!(ids(file.properties(0).unwrap()), vec!["x"]);
    assert_eq!(ids(file.properties(1).unwrap()), vec!["y"]);
}

// XISF 1.0 section 11.5.1: bounds are two floating point numbers, the black
// point and the white point.
#[test]
fn malformed_bounds_are_an_error_for_that_image_only() {
    for bounds in [
        "abc:xyz",
        "5:invalid",
        "invalid:10",
        "0:1:2",
        "",
        "1",
        "NaN:1",
        "0:inf",
        "-inf:0",
        "2:1",
    ] {
        let attrs = format!(r#"geometry="2:2:1" sampleFormat="Float32" bounds="{bounds}""#);
        let images = format!("{}{}", attached(&attrs, 0, 16), attached(GRAY_2X2, 16, 4));
        let mut data = vec![0u8; 16];
        data.extend_from_slice(&[1, 2, 3, 4]);
        let mut file = open(&images, &data);
        assert_eq!(
            read_error(&mut file, 0),
            format!(
                "Invalid XISF format: Invalid bounds attribute \"{bounds}\": \
                 it is not two finite numbers lower:upper with lower <= upper"
            )
        );
        assert_eq!(file.read_image_data_raw(1).unwrap(), vec![1, 2, 3, 4]);
    }
}

#[test]
fn well_formed_bounds_are_read() {
    for (bounds, expected) in [
        ("0:1", (0.0, 1.0)),
        ("-1.5:1.5", (-1.5, 1.5)),
        ("5:5", (5.0, 5.0)),
        ("1e3:2E3", (1000.0, 2000.0)),
    ] {
        let attrs = format!(r#"geometry="2:2:1" sampleFormat="Float32" bounds="{bounds}""#);
        let mut file = open(&attached(&attrs, 0, 16), &[0; 16]);
        assert_eq!(file.image_info(0).unwrap().bounds, expected);
        assert_eq!(file.read_image_data_raw(0).unwrap(), vec![0; 16]);
    }
}

// XISF 1.0 section 11.5.1: geometry is dim1:...:dimN:channel-count, so two
// fields describe a one-dimensional image. Its dimensions follow a FITS
// NAXIS = 1 image, with channels as the third dimension.
#[test]
fn two_field_geometry_is_one_dimensional() {
    let samples: Vec<u8> = (1..=12).collect();
    for (geometry, dimensions, channels) in [
        ("12:1", vec![12usize], 1),
        ("4:3", vec![4, 1, 3], 3),
        ("4:3:1", vec![4, 3], 1),
        ("2:2:3", vec![2, 2, 3], 3),
    ] {
        let attrs = format!(r#"geometry="{geometry}" sampleFormat="UInt8""#);
        let mut file = open(&attached(&attrs, 0, 12), &samples);
        let info = file.image_info(0).unwrap();
        assert_eq!(info.geometry, dimensions, "{geometry}");
        assert_eq!(info.num_channels(), channels, "{geometry}");
        assert_eq!(file.read_image_data_raw(0).unwrap(), samples);
    }
}

#[test]
fn one_dimensional_normal_storage_is_read_in_planar_order() {
    let attrs = r#"geometry="2:3" sampleFormat="UInt8" pixelStorage="Normal""#;
    let mut file = open(&attached(attrs, 0, 6), &[1, 2, 3, 4, 5, 6]);
    assert_eq!(file.read_image_data_raw(0).unwrap(), vec![1, 4, 2, 5, 3, 6]);
}

#[test]
fn images_with_more_than_two_dimensions_are_rejected() {
    let images = format!(
        r#"<Image geometry="2:2:2:1" sampleFormat="UInt8" location="attachment:{DATA_START}:8"/>"#
    );
    let err = XisfFile::new(Cursor::new(xisf(&images, &[0; 8]))).unwrap_err();
    assert_eq!(
        err.to_string(),
        "Invalid geometry: 2:2:2:1 has 3 dimensions; images with more than two are not supported"
    );
}

// Degenerate bounds are the caller's choice; bounds that are not two finite
// numbers in order are an error.
#[test]
fn caller_bounds_must_be_finite_and_ordered() {
    let write = |bounds: (f64, f64)| {
        let mut writer = XisfWriter::new(Cursor::new(Vec::new()));
        writer.add_image_with_bounds(&[0.5f32; 4], 2, 2, 1, bounds)
    };
    assert!(write((5.0, 5.0)).is_ok());
    assert!(write((-1.0, 2.0)).is_ok());
    for bounds in [
        (f64::NAN, 1.0),
        (0.0, f64::INFINITY),
        (f64::NEG_INFINITY, 0.0),
        (2.0, 1.0),
    ] {
        assert_eq!(
            write(bounds).unwrap_err().to_string(),
            format!(
                "Invalid XISF format: Image bounds ({}, {}) are not two finite numbers \
                 with lower <= upper",
                bounds.0, bounds.1
            )
        );
    }
}

const ZEROS_SHA1: &str = "0000000000000000000000000000000000000000";

// XISF 1.0 section 10.6 moves only the compression and subblocks attributes
// of an embedded block to its Data element. Its checksum and byteOrder stay
// on the Image element (section 10.5), so the Data element's are ignored.
#[test]
fn checksum_and_byte_order_on_the_data_element_are_ignored() {
    let data = |attrs: &str| format!(r#"<Data encoding="base64" {attrs}>AQIDBA==</Data>"#);
    let wrong_on_image = format!(r#"{GRAY_2X2} checksum="sha1:{ZEROS_SHA1}""#);
    let right_on_data = data(&format!(r#"checksum="sha1:{SHA1_1234}""#));
    let mut file = open(&embedded(&wrong_on_image, &right_on_data), &[]);
    assert!(read_error(&mut file, 0).starts_with("Checksum mismatch: "));

    let right_on_image = format!(r#"{GRAY_2X2} checksum="sha1:{SHA1_1234}""#);
    let wrong_on_data = data(&format!(r#"checksum="sha1:{ZEROS_SHA1}""#));
    let mut file = open(&embedded(&right_on_image, &wrong_on_data), &[]);
    assert_eq!(file.read_image_data_raw(0).unwrap(), vec![1, 2, 3, 4]);

    // The samples 1, 2, 3 and 4 as big-endian UInt16.
    let samples = r#"<Data encoding="base64" byteOrder="little">AAEAAgADAAQ=</Data>"#;
    let big = r#"geometry="2:2:1" sampleFormat="UInt16" byteOrder="big""#;
    let mut file = open(&embedded(big, samples), &[]);
    assert_eq!(
        file.read_image_data_raw(0).unwrap(),
        vec![1, 0, 2, 0, 3, 0, 4, 0]
    );

    let samples = r#"<Data encoding="base64" byteOrder="big">AAEAAgADAAQ=</Data>"#;
    let little = r#"geometry="2:2:1" sampleFormat="UInt16""#;
    let mut file = open(&embedded(little, samples), &[]);
    assert_eq!(
        file.read_image_data_raw(0).unwrap(),
        vec![0, 1, 0, 2, 0, 3, 0, 4]
    );
}

// XISF 1.0 section 10.6: subblocks describe the compressed block, so they
// are given on the element that has the compression attribute.
#[test]
fn subblocks_must_be_on_the_element_with_the_compression() {
    let attrs = r#"geometry="4:4:1" sampleFormat="UInt8""#;
    let on_data = format!(
        r#"<Data encoding="base64" compression="zlib:16" subblocks="11,16">{ZLIB_SEVENS_BASE64}</Data>"#
    );
    let mut file = open(&embedded(attrs, &on_data), &[]);
    assert_eq!(file.read_image_data_raw(0).unwrap(), vec![7; 16]);

    let on_image = format!(r#"{attrs} subblocks="11,16""#);
    let data =
        format!(r#"<Data encoding="base64" compression="zlib:16">{ZLIB_SEVENS_BASE64}</Data>"#);
    let message = "Invalid XISF format: A subblocks attribute needs a compression attribute \
                   on the same element";
    let mut file = open(&embedded(&on_image, &data), &[]);
    assert_eq!(read_error(&mut file, 0), message);

    let attached_attrs = format!(r#"{GRAY_2X2} subblocks="4,4""#);
    let mut file = open(&attached(&attached_attrs, 0, 4), &[1, 2, 3, 4]);
    assert_eq!(read_error(&mut file, 0), message);
}

// XISF 1.0 section 11.5.1 makes bounds optional for integer images, and no
// sample decoding uses them, so malformed ones are ignored as if the
// attribute were absent.
#[test]
fn malformed_bounds_of_an_integer_image_are_ignored() {
    for (format, size, high) in [("UInt8", 2, 255.0), ("UInt16", 4, 65535.0)] {
        for bounds in ["abc:xyz", "5:invalid", "", "NaN:1", "0:inf", "2:1"] {
            let attrs = format!(r#"geometry="2:1:1" sampleFormat="{format}" bounds="{bounds}""#);
            let mut file = open(&attached(&attrs, 0, size), &[1, 2, 3, 4]);
            let info = file.image_info(0).unwrap();
            assert_eq!(info.bounds, (0.0, high), "{format} {bounds}");
            assert_eq!(info.block_defect, None, "{format} {bounds}");
            assert_eq!(file.read_image_data_raw(0).unwrap(), [1, 2, 3, 4][..size]);
        }
    }
}

// XISF 1.0 section 11.5.1: without a bounds attribute, an integer image has
// the range of its k-bit sample type, [0, 2^k - 1]. 2^64 - 1 needs 64
// significant bits, so its nearest f64 is 2^64 = 18446744073709551616.
#[test]
fn integer_image_without_bounds_has_the_range_of_its_sample_type() {
    let formats = [
        ("UInt8", 1, 255.0),
        ("UInt16", 2, 65535.0),
        ("UInt32", 4, 4294967295.0),
        ("UInt64", 8, 18446744073709551616.0),
    ];
    for (format, size, high) in formats {
        let attrs = format!(r#"geometry="1:1:1" sampleFormat="{format}""#);
        let file = open(&attached(&attrs, 0, size), &vec![0u8; size]);
        assert_eq!(file.image_info(0).unwrap().bounds, (0.0, high), "{format}");
    }
}

// XML 1.0 section 3.3.3: attribute values are read with character and
// entity references replaced.
#[test]
fn image_attributes_decode_character_references() {
    let attrs = format!(r#"{GRAY_2X2} checksum="sha&#x31;:{SHA1_1234}""#);
    let mut file = open(&attached(&attrs, 0, 4), &[1, 2, 3, 4]);
    assert_eq!(file.read_image_data_raw(0).unwrap(), vec![1, 2, 3, 4]);
}

// A block in the header is described by its location attribute, and its
// decoded bytes are not printed.
#[test]
fn a_header_block_is_not_described_as_an_attachment() {
    let sevens =
        "BwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBwcHBw==";
    let attrs = r#"geometry="8:8:1" sampleFormat="UInt8""#;
    let images = format!(
        "{}{}",
        inline(attrs, "base64", sevens),
        embedded(
            attrs,
            &format!(r#"<Data encoding="base64">{sevens}</Data>"#)
        )
    );
    let file = open(&images, &[]);
    for (index, location) in [(0, "inline:base64"), (1, "embedded")] {
        let info = file.image_info(index).unwrap();
        assert_eq!(info.location.format(), location);
        let debug = format!("{info:?}");
        assert!(!debug.contains("7, 7"), "{debug}");
        assert!(debug.contains("64 bytes"), "{debug}");
    }
}

#[test]
fn two_data_elements_are_an_error_for_that_image_only() {
    let children = r#"<Data encoding="base64">AQIDBA==</Data><Data encoding="base64">AAAA</Data>"#;
    let images = format!(
        "{}{}",
        embedded(GRAY_2X2, children),
        attached(GRAY_2X2, 0, 4)
    );
    let mut file = open(&images, &[1, 2, 3, 4]);
    assert_eq!(
        read_error(&mut file, 0),
        "Invalid XISF format: An embedded data block has more than one Data element"
    );
    assert_eq!(file.read_image_data_raw(1).unwrap(), vec![1, 2, 3, 4]);
}
