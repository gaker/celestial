use super::card;
use crate::core::types::{BitPix, ByteOrder};
use crate::fits::block::{pad_to_block, CARD_SIZE};
use crate::fits::compression::encode::encode_image;
use crate::fits::compression::CompressionAlgorithm;
use crate::fits::data::array::DataArray;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::builder::{primary_cards, primary_header_cards};
use crate::fits::header::keywords::Keyword;
use crate::fits::header::parser::Header;
use crate::fits::util::checksum;
use crate::fits::util::image::shape_text;

const CHECKSUM_PLACEHOLDER: &str = "0000000000000000";

// The bytes of a FITS file whose primary HDU holds `data`.
pub(crate) fn primary_image_bytes<T: DataArray>(
    data: &[T],
    dimensions: &[usize],
    keywords: &[Keyword],
    with_checksum: bool,
) -> Result<Vec<u8>> {
    check_data_length(data.len(), dimensions)?;
    let mut header = build_primary_header(T::BITPIX, dimensions, keywords);
    let mut padded_data = if data.is_empty() {
        Vec::new()
    } else {
        T::to_bytes(data, ByteOrder::BigEndian)?
    };
    pad_to_block(&mut padded_data, 0);

    let mut bytes = if with_checksum {
        checksummed_header_bytes(&mut header, &padded_data)?
    } else {
        header.add_keyword(Keyword::new("END".to_string()));
        header_bytes(&header)?
    };
    bytes.extend(padded_data);
    Ok(bytes)
}

// The bytes of a FITS file with an empty primary HDU and the image
// tile-compressed in the first extension.
pub(crate) fn compressed_image_bytes<T: DataArray>(
    data: &[T],
    dimensions: &[usize],
    tile_size: Option<(usize, usize)>,
    algorithm: CompressionAlgorithm,
    keywords: &[Keyword],
) -> Result<Vec<u8>> {
    check_data_length(data.len(), dimensions)?;
    let pixels = T::to_bytes(data, ByteOrder::BigEndian)?;
    let image = encode_image(
        &pixels,
        T::BITPIX,
        dimensions,
        tile_size,
        algorithm,
        keywords,
    )?;
    let mut header = image.header;
    header.add_keyword(Keyword::new("END".to_string()));
    let mut bytes = header_bytes(&minimal_primary_header())?;
    bytes.extend(header_bytes(&header)?);
    bytes.extend(image.data);
    Ok(bytes)
}

fn build_primary_header(bitpix: BitPix, dimensions: &[usize], keywords: &[Keyword]) -> Header {
    let mut header = Header::new();
    for card in primary_header_cards(bitpix, dimensions, keywords) {
        header.add_keyword(card);
    }
    add_keywords(&mut header, keywords);
    header
}

// Adds DATASUM, CHECKSUM and END to `header` and returns its bytes with the
// CHECKSUM value that makes the HDU sum to zero.
fn checksummed_header_bytes(header: &mut Header, data_bytes: &[u8]) -> Result<Vec<u8>> {
    let data_sum = checksum::calculate_datasum(data_bytes);
    header.add_keyword(Keyword::string(
        "DATASUM",
        &checksum::format_datasum(data_sum),
    ));
    header.add_keyword(Keyword::string("CHECKSUM", CHECKSUM_PLACEHOLDER));
    header.add_keyword(Keyword::new("END".to_string()));

    let with_placeholder = header_bytes_with_checksum(header, CHECKSUM_PLACEHOLDER)?;
    let checksum_value = checksum::create_checksum_card_value(&with_placeholder, data_sum);
    header_bytes_with_checksum(header, &checksum_value)
}

fn header_bytes_with_checksum(header: &Header, checksum_value: &str) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    for keyword in header.keywords() {
        if keyword.name == "CHECKSUM" {
            bytes.extend_from_slice(&format_checksum_card(checksum_value));
        } else {
            bytes.extend(format_card(keyword)?);
        }
    }
    pad_to_block(&mut bytes, b' ');
    Ok(bytes)
}

fn format_checksum_card(checksum_value: &str) -> [u8; CARD_SIZE] {
    let mut card = [b' '; CARD_SIZE];
    card[0..8].copy_from_slice(b"CHECKSUM");
    card[8] = b'=';
    card[9] = b' ';
    let value_str = format!("'{}'", checksum_value);
    let value_bytes = value_str.as_bytes();
    card[10..10 + value_bytes.len()].copy_from_slice(value_bytes);
    card
}

fn header_bytes(header: &Header) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    for keyword in header.keywords() {
        bytes.extend(format_card(keyword)?);
    }
    pad_to_block(&mut bytes, b' ');
    Ok(bytes)
}

fn format_card(keyword: &Keyword) -> Result<Vec<u8>> {
    card::format_keyword(keyword)
}

fn add_keywords(header: &mut Header, keywords: &[Keyword]) {
    for keyword in keywords {
        // Mandatory keywords come from the data being written. An incoming
        // checksum describes the old header and data, so it is dropped and
        // recomputed when a checksum is requested.
        if !keyword.is_mandatory() && !is_checksum_keyword(&keyword.name) {
            header.add_keyword(keyword.clone());
        }
    }
}

// The empty primary HDU before a tile-compressed image extension.
fn minimal_primary_header() -> Header {
    let mut header = Header::new();
    for card in primary_cards(BitPix::U8, &[]) {
        header.add_keyword(card);
    }
    header.add_keyword(Keyword::new("END".to_string()));
    header
}

pub(crate) fn is_checksum_keyword(name: &str) -> bool {
    matches!(
        name.trim().to_ascii_uppercase().as_str(),
        "CHECKSUM" | "DATASUM" | "ZHECKSUM" | "ZDATASUM"
    )
}

// NAXIS = 0 declares no data array, so an empty `dimensions` expects no data.
// NAXISn is written as an i64, so a larger axis cannot be represented.
fn check_data_length(len: usize, dimensions: &[usize]) -> Result<()> {
    if let Some(axis) = dimensions.iter().find(|&&d| i64::try_from(d).is_err()) {
        return Err(FitsError::InvalidFormat(format!(
            "axis length {} exceeds the largest NAXISn value",
            axis
        )));
    }
    let expected = if dimensions.is_empty() {
        Some(0)
    } else {
        dimensions.iter().try_fold(1usize, |n, &d| n.checked_mul(d))
    };
    match expected {
        Some(expected) if expected == len => Ok(()),
        Some(_) if dimensions.is_empty() => Err(FitsError::InvalidFormat(format!(
            "data has {} values but NAXIS = 0 declares no data",
            len
        ))),
        Some(expected) => Err(FitsError::InvalidFormat(format!(
            "data has {} values but dimensions {} need {}",
            len,
            shape_text(dimensions),
            expected
        ))),
        None => Err(FitsError::InvalidFormat(format!(
            "dimensions {} overflow the pixel count",
            shape_text(dimensions)
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::block::FITS_BLOCK_SIZE;
    use crate::fits::errors::FitsError;
    use crate::fits::header::keywords::{Keyword, KeywordValue};
    use crate::fits::io::test_support::FitsWriter;
    use tempfile::NamedTempFile;

    #[test]
    fn fits_writer_create_success() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let writer = FitsWriter::create(path);
        assert!(writer.is_ok());
    }

    #[test]
    fn write_primary_image_i16_data() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let data: Vec<i16> = vec![1, 2, 3, 4, 5, 6];
        let dimensions = vec![2, 3];

        let result = writer.write_primary_image(&data, &dimensions, &[]);
        assert!(result.is_ok());
    }

    #[test]
    fn write_primary_image_with_keywords() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let data: Vec<u8> = vec![0, 1, 2, 3];
        let dimensions = vec![2, 2];

        let keywords = vec![
            Keyword::string("TELESCOP", "Test Telescope"),
            Keyword::real("EXPTIME", 30.5),
            Keyword::logical("ISLIGHT", true),
            Keyword::integer("COUNT", 42),
        ];

        let result = writer.write_primary_image(&data, &dimensions, &keywords);
        assert!(result.is_ok());
    }

    #[test]
    fn write_primary_image_empty_data() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let mut writer = FitsWriter::create(path).unwrap();

        let data: Vec<f32> = vec![];
        let dimensions = vec![0];

        let result = writer.write_primary_image(&data, &dimensions, &[]);
        assert!(result.is_ok());
    }

    #[test]
    fn format_card_end_keyword() {
        let end_keyword = Keyword::new("END".to_string());

        let card = format_card(&end_keyword).unwrap();
        assert_eq!(&card[0..3], b"END");
        assert_eq!(card[3], b' ');
    }

    #[test]
    fn format_card_logical_value() {
        let keyword = Keyword::logical("SIMPLE", true);

        let card = format_card(&keyword).unwrap();
        assert_eq!(&card[0..6], b"SIMPLE");
        assert_eq!(card[8], b'=');
        assert_eq!(card[9], b' ');
        assert!(String::from_utf8_lossy(&card[10..30]).contains('T'));
    }

    #[test]
    fn format_card_integer_value() {
        let keyword = Keyword::integer("NAXIS", 2);

        let card = format_card(&keyword).unwrap();
        assert_eq!(&card[0..5], b"NAXIS");
        assert_eq!(card[8], b'=');
        assert!(String::from_utf8_lossy(&card[10..30]).contains('2'));
    }

    #[test]
    fn format_card_real_value() {
        let keyword = Keyword::real("EXPTIME", 30.5);

        let card = format_card(&keyword).unwrap();
        assert_eq!(&card[0..7], b"EXPTIME");
        assert_eq!(card[8], b'=');
        assert!(String::from_utf8_lossy(&card[10..30]).contains("30.5"));
    }

    #[test]
    fn format_card_string_value() {
        let keyword = Keyword::string("OBJECT", "M31");

        let card = format_card(&keyword).unwrap();
        assert_eq!(&card[0..6], b"OBJECT");
        assert_eq!(card[8], b'=');
        assert!(String::from_utf8_lossy(&card[10..30]).contains("'M31"));
    }

    #[test]
    fn format_card_with_comment() {
        let mut keyword = Keyword::integer("BITPIX", 16);
        keyword.comment = Some("Bits per pixel".to_string());

        let card = format_card(&keyword).unwrap();
        assert_eq!(card[31], b'/');
        assert!(String::from_utf8_lossy(&card[33..]).contains("Bits per pixel"));
    }

    #[test]
    fn format_card_complex_value() {
        let mut keyword = Keyword::new("COMPLEX".to_string());
        keyword.value = Some(KeywordValue::Complex(1.0, 2.0));

        let card = format_card(&keyword).unwrap();
        let value_section = String::from_utf8_lossy(&card[10..30]);
        assert!(value_section.contains("1"));
        assert!(value_section.contains("2"));
    }

    #[test]
    fn header_bytes_pad_to_one_block() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::new("END".to_string()));

        assert_eq!(header_bytes(&header).unwrap().len(), FITS_BLOCK_SIZE);
    }

    #[test]
    fn write_primary_image_with_checksum_creates_valid_file() {
        use crate::fits::io::reader::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let data: Vec<i16> = (0..100).collect();
        let dimensions = vec![10, 10];

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_primary_image_with_checksum(&data, &dimensions, &[])
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        assert_eq!(fits.num_hdus(), 1);

        let header = fits.get_header(0).unwrap();
        assert!(header.get("DATASUM").exists());
        assert!(header.get("CHECKSUM").exists());
    }

    #[test]
    fn write_with_checksum_validates_on_read() {
        use crate::fits::io::reader::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let data: Vec<i16> = (0..64).collect();
        let dimensions = vec![8, 8];

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_primary_image_with_checksum(&data, &dimensions, &[])
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        let result = fits.validate_hdu_checksum(0);
        assert!(
            result.is_ok(),
            "Checksum validation call failed: {:?}",
            result.err()
        );

        let checksum_result = result.unwrap();
        assert!(checksum_result.has_checksums());
        assert!(checksum_result.datasum_valid == Some(true));
    }

    #[test]
    fn write_checksum_empty_data() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let data: Vec<u8> = vec![];
        let dimensions = vec![0];

        let mut writer = FitsWriter::create(path).unwrap();
        let result = writer.write_primary_image_with_checksum(&data, &dimensions, &[]);
        assert!(result.is_ok());
    }

    #[test]
    fn format_checksum_card_correct_format() {
        let checksum_value = "1234567890ABCDEF";

        let card = format_checksum_card(checksum_value);
        assert_eq!(&card[0..8], b"CHECKSUM");
        assert_eq!(card[8], b'=');
        assert_eq!(card[9], b' ');
        assert!(String::from_utf8_lossy(&card).contains(checksum_value));
    }

    #[test]
    fn build_primary_header_includes_required_keywords() {
        let dimensions = vec![10, 20];

        let header = build_primary_header(BitPix::I16, &dimensions, &[]);

        assert!(header.get("SIMPLE").exists());
        assert!(header.get("BITPIX").exists());
        assert!(header.get("NAXIS").exists());
        assert!(header.get("NAXIS1").exists());
        assert!(header.get("NAXIS2").exists());
    }

    #[test]
    fn roundtrip_2d_image_preserves_orientation() {
        use crate::fits::io::reader::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let original: Vec<i16> = (0..100).collect();
        let dimensions = vec![10, 10];

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_primary_image(&original, &dimensions, &[])
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        let (_, read_back): (_, Vec<i16>) = fits.primary_hdu_with_data().unwrap();

        assert_eq!(original, read_back);
    }

    #[test]
    fn roundtrip_gradient_image_correct_orientation() {
        use crate::fits::io::reader::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let mut gradient: Vec<i16> = Vec::with_capacity(25);
        for row in 0..5 {
            for col in 0..5 {
                gradient.push((row * 10 + col) as i16);
            }
        }
        let dimensions = vec![5, 5];

        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_primary_image(&gradient, &dimensions, &[])
                .unwrap();
        }

        let mut fits = FitsFile::open(path).unwrap();
        let (_, read_back): (_, Vec<i16>) = fits.primary_hdu_with_data().unwrap();

        assert_eq!(read_back[0], 0);
        assert_eq!(read_back[4], 4);
        assert_eq!(read_back[20], 40);
        assert_eq!(read_back[24], 44);
        assert_eq!(gradient, read_back);
    }

    fn card_text(keyword: &Keyword) -> String {
        let cards = format_card(keyword).unwrap();
        String::from_utf8(cards[..].to_vec()).unwrap()
    }

    fn padded(cards: &[&str]) -> String {
        cards.iter().map(|c| format!("{:<80}", c)).collect()
    }

    fn real_value_text(card: &str) -> &str {
        let value = &card[10..];
        value.split(" /").next().unwrap().trim()
    }

    #[test]
    fn format_card_string_keeps_values_up_to_68_chars() {
        assert_eq!(
            card_text(&Keyword::string("DATE-OBS", "2024-01-15T12:34:56.789")),
            padded(&["DATE-OBS= '2024-01-15T12:34:56.789'"])
        );
        assert_eq!(
            card_text(&Keyword::string("OBJECT", "NGC 1234 field center")),
            padded(&["OBJECT  = 'NGC 1234 field center'"])
        );
        let longest = "x".repeat(68);
        assert_eq!(
            card_text(&Keyword::string("LONGEST", &longest)),
            padded(&[&format!("LONGEST = '{}'", longest)])
        );
    }

    #[test]
    fn format_card_string_escapes_quotes() {
        assert_eq!(
            card_text(&Keyword::string("OBSERVER", "O'Brien")),
            padded(&["OBSERVER= 'O''Brien'"])
        );
        assert_eq!(
            card_text(&Keyword::string("QUOTE", "'")),
            padded(&["QUOTE   = '''      '"])
        );
    }

    #[test]
    fn format_card_string_pads_short_values_to_eight_chars() {
        assert_eq!(
            card_text(&Keyword::string("OBJECT", "M31").with_comment("target")),
            padded(&["OBJECT  = 'M31     '           / target"])
        );
    }

    #[test]
    fn format_card_long_string_uses_continue_cards() {
        let value: String = (0..100).map(|i| (b'0' + (i % 10) as u8) as char).collect();
        let keyword = Keyword::string("LONGSTR", &value).with_comment("long value");
        assert_eq!(
            card_text(&keyword),
            padded(&[
                &format!("LONGSTR = '{}&'", &value[..67]),
                &format!("CONTINUE  '{}' / long value", &value[67..]),
            ])
        );
    }

    #[test]
    fn format_card_long_string_keeps_escaped_quote_on_one_card() {
        // The 67th escaped character is the first half of a '' pair.
        let value = format!("{}'{}", "a".repeat(66), "b".repeat(10));
        assert_eq!(
            card_text(&Keyword::string("SPLIT", &value)),
            padded(&[
                &format!("SPLIT   = '{}&'", "a".repeat(66)),
                &format!("CONTINUE  '''{}'", "b".repeat(10)),
            ])
        );
    }

    #[test]
    fn format_card_string_rejects_non_ascii_and_control_chars() {
        for value in ["caf\u{e9}", "tab\there", "del\u{7f}"] {
            let result = format_card(&Keyword::string("OBJECT", value));
            assert!(
                matches!(result, Err(FitsError::InvalidKeywordValue { .. })),
                "{:?} was accepted",
                value
            );
        }
    }

    #[test]
    fn format_card_string_with_multibyte_char_returns_invalid_value_error() {
        // Byte 18 falls inside a two-byte character, so cutting the value there would panic.
        let value = format!("a{}", "\u{3b1}".repeat(10));
        let result = format_card(&Keyword::string("OBJECT", &value));
        assert!(matches!(result, Err(FitsError::InvalidKeywordValue { .. })));
    }

    #[test]
    fn format_card_real_round_trips_bits() {
        let values = [
            -2.777_777_777_777_778e-4,
            -1.234_567_890_123_456_8e-10,
            2.345_678_901_234_568e-14,
            4.57e-19,
            -7.5e-19,
            1e25,
            512.0,
            1.0,
            -0.0,
            0.1,
            celestial_core::constants::PI,
            f64::MAX,
            f64::MIN_POSITIVE,
            5e-324,
        ];
        for value in values {
            let card = card_text(&Keyword::real("TEST", value).with_comment("c"));
            let text = real_value_text(&card);
            assert!(
                text.contains('.') && text.parse::<i64>().is_err(),
                "{} written as {:?}",
                value,
                text
            );
            assert_eq!(
                text.parse::<f64>().unwrap().to_bits(),
                value.to_bits(),
                "{}",
                text
            );
        }
    }

    #[test]
    fn format_card_real_layout() {
        assert_eq!(
            card_text(&Keyword::real("CRPIX1", 512.0)),
            padded(&["CRPIX1  =                512.0"])
        );
        assert_eq!(
            card_text(&Keyword::real("A_5_0", 4.57e-19).with_comment("SIP")),
            padded(&["A_5_0   =             4.57E-19 / SIP"])
        );
        assert_eq!(
            card_text(&Keyword::real("TEST", 1e25)),
            padded(&["TEST    =              1.0E+25"])
        );
        assert_eq!(
            card_text(&Keyword::real("CD1_1", -2.777_777_777_777_778e-4).with_comment("CD")),
            padded(&["CD1_1   = -0.0002777777777777778 / CD"])
        );
    }

    #[test]
    fn format_card_real_rejects_non_finite() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let result = format_card(&Keyword::real("BAD", value));
            assert!(
                matches!(result, Err(FitsError::InvalidKeywordValue { .. })),
                "{} was accepted",
                value
            );
        }
    }

    fn written_hdu(path: &std::path::Path) -> (Vec<u8>, Vec<u8>) {
        let bytes = std::fs::read(path).unwrap();
        let end = bytes
            .chunks(80)
            .position(|card| card.starts_with(b"END     "))
            .unwrap();
        let header_len = (end * 80 / FITS_BLOCK_SIZE + 1) * FITS_BLOCK_SIZE;
        let (header, data) = bytes.split_at(header_len);
        (header.to_vec(), data.to_vec())
    }

    #[test]
    fn write_with_checksum_writes_valid_checksum() {
        use crate::fits::io::reader::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();
        let data: Vec<i16> = (0..100).collect();
        {
            let mut writer = FitsWriter::create(path).unwrap();
            writer
                .write_primary_image_with_checksum(&data, &[10, 10], &[])
                .unwrap();
        }

        let (header, data_bytes) = written_hdu(path);
        assert!(checksum::verify_hdu_checksum(&header, &data_bytes));

        let mut fits = FitsFile::open(path).unwrap();
        let result = fits.validate_hdu_checksum(0).unwrap();
        assert_eq!(result.checksum_valid, Some(true));
        assert_eq!(result.datasum_valid, Some(true));
    }

    fn checksum_result_of(bytes: Vec<u8>) -> (Option<bool>, Option<bool>) {
        let mut fits = crate::fits::io::reader::FitsFile::new(std::io::Cursor::new(bytes)).unwrap();
        let result = fits.validate_hdu_checksum(0).unwrap();
        (result.datasum_valid, result.checksum_valid)
    }

    #[test]
    fn checksum_mismatches_are_reported_as_flags() {
        let data: Vec<i16> = (0..100).collect();
        let bytes = primary_image_bytes(&data, &[10, 10], &[], true).unwrap();
        assert_eq!(checksum_result_of(bytes.clone()), (Some(true), Some(true)));

        let mut data_changed = bytes.clone();
        data_changed[FITS_BLOCK_SIZE + 1] ^= 1;
        assert_eq!(checksum_result_of(data_changed), (Some(false), Some(false)));

        let mut header_changed = bytes;
        let comment = header_changed
            .windows(4)
            .position(|w| w == b"file")
            .unwrap();
        header_changed[comment] = b'F';
        assert_eq!(
            checksum_result_of(header_changed),
            (Some(true), Some(false))
        );
    }

    #[test]
    fn compressed_file_primary_header_has_the_image_header_comments() {
        let bytes =
            compressed_image_bytes(&[0i16; 4], &[2, 2], None, CompressionAlgorithm::Gzip, &[])
                .unwrap();
        let header =
            crate::fits::header::parser::HeaderParser::parse_header(&bytes[..2880]).unwrap();
        let cards: Vec<_> = header
            .keywords()
            .iter()
            .map(|k| (k.name.as_str(), k.comment.as_deref()))
            .collect();
        assert_eq!(
            cards,
            [
                ("SIMPLE", Some("file does conform to FITS standard")),
                ("BITPIX", Some("number of bits per data pixel")),
                ("NAXIS", Some("number of data axes")),
                ("EXTEND", Some("FITS dataset may contain extensions")),
            ]
        );
    }

    fn count_cards(header: &[u8], name: &[u8; 8]) -> usize {
        header.chunks(80).filter(|c| c.starts_with(name)).count()
    }

    #[test]
    fn write_drops_incoming_checksum_and_datasum() {
        let stale = [
            Keyword::string("CHECKSUM", "9aa2FaX19aa1CaW1"),
            Keyword::string("DATASUM", "4266173065"),
            Keyword::string("OBJECT", "M31"),
        ];
        let data: Vec<i16> = (0..16).collect();

        let plain = NamedTempFile::new().unwrap();
        FitsWriter::create(plain.path())
            .unwrap()
            .write_primary_image(&data, &[4, 4], &stale)
            .unwrap();
        let (header, _) = written_hdu(plain.path());
        assert_eq!(count_cards(&header, b"CHECKSUM"), 0);
        assert_eq!(count_cards(&header, b"DATASUM "), 0);
        assert_eq!(count_cards(&header, b"OBJECT  "), 1);

        let summed = NamedTempFile::new().unwrap();
        FitsWriter::create(summed.path())
            .unwrap()
            .write_primary_image_with_checksum(&data, &[4, 4], &stale)
            .unwrap();
        let (header, data_bytes) = written_hdu(summed.path());
        assert_eq!(count_cards(&header, b"CHECKSUM"), 1);
        assert_eq!(count_cards(&header, b"DATASUM "), 1);
        assert!(checksum::verify_hdu_checksum(&header, &data_bytes));
    }

    #[test]
    fn write_rejects_data_length_that_does_not_match_dimensions() {
        let temp_file = NamedTempFile::new().unwrap();
        let mut writer = FitsWriter::create(temp_file.path()).unwrap();
        match writer.write_primary_image(&[1i16, 2, 3], &[2, 2], &[]) {
            Err(FitsError::InvalidFormat(m)) => {
                assert_eq!(m, "data has 3 values but dimensions 2x2 need 4")
            }
            other => panic!("{other:?}"),
        }
        assert!(matches!(
            writer.write_primary_image_with_checksum(&[1i16, 2, 3, 4, 5], &[2, 2], &[]),
            Err(FitsError::InvalidFormat(_))
        ));
        match writer.write_primary_image(&[1i16], &[], &[]) {
            Err(FitsError::InvalidFormat(m)) => {
                assert_eq!(m, "data has 1 values but NAXIS = 0 declares no data")
            }
            other => panic!("{other:?}"),
        }
        assert!(matches!(
            writer.write_compressed_image(
                &[1i16, 2, 3],
                &[2, 2],
                Some((2, 2)),
                CompressionAlgorithm::Rice,
                &[]
            ),
            Err(FitsError::InvalidFormat(_))
        ));
    }

    fn count_names(header: &[u8], name: &str) -> usize {
        header
            .chunks(80)
            .filter(|c| {
                String::from_utf8_lossy(&c[..8])
                    .trim()
                    .eq_ignore_ascii_case(name)
            })
            .count()
    }

    #[test]
    fn write_drops_checksum_keywords_in_any_spelling() {
        let stale = [
            Keyword::string("CHECKSUM ", "9aa2FaX19aa1CaW1"),
            Keyword::string("checksum", "9aa2FaX19aa1CaW1"),
            Keyword::string(" DataSum", "4266173065"),
            Keyword::string("ZHECKSUM", "9aa2FaX19aa1CaW1"),
            Keyword::string("ZDATASUM", "4266173065"),
        ];
        let data: Vec<i16> = (0..16).collect();

        let summed = NamedTempFile::new().unwrap();
        FitsWriter::create(summed.path())
            .unwrap()
            .write_primary_image_with_checksum(&data, &[4, 4], &stale)
            .unwrap();
        let (header, data_bytes) = written_hdu(summed.path());
        assert_eq!(count_names(&header, "CHECKSUM"), 1);
        assert_eq!(count_names(&header, "DATASUM"), 1);
        assert_eq!(count_names(&header, "ZHECKSUM"), 0);
        assert_eq!(count_names(&header, "ZDATASUM"), 0);
        assert!(checksum::verify_hdu_checksum(&header, &data_bytes));

        let compressed = NamedTempFile::new().unwrap();
        FitsWriter::create(compressed.path())
            .unwrap()
            .write_compressed_image(
                &data,
                &[4, 4],
                Some((4, 4)),
                CompressionAlgorithm::Rice,
                &stale,
            )
            .unwrap();
        let bytes = std::fs::read(compressed.path()).unwrap();
        for name in ["CHECKSUM", "DATASUM", "ZHECKSUM", "ZDATASUM"] {
            assert_eq!(count_names(&bytes, name), 0, "{}", name);
        }
    }

    #[test]
    fn format_card_long_string_drops_trailing_spaces() {
        let value = format!("{}{}", "d".repeat(60), " ".repeat(20));
        assert_eq!(
            card_text(&Keyword::string("TRAIL", &value)),
            padded(&[&format!("TRAIL   = '{}'", "d".repeat(60))])
        );
        let value = format!("{}{}", "e".repeat(100), " ".repeat(40));
        assert_eq!(
            card_text(&Keyword::string("TRAIL", &value)),
            padded(&[
                &format!("TRAIL   = '{}&'", "e".repeat(67)),
                &format!("CONTINUE  '{}'", "e".repeat(33)),
            ])
        );
    }

    #[test]
    fn format_card_rejects_non_ascii_in_names_comments_and_text() {
        let keywords = [
            Keyword::integer("TEMP", 1).with_comment("deg \u{b0}C"),
            Keyword::string("OBJECT", "M31").with_comment("bell\u{7}"),
            Keyword::real("EXPTIME", 1.0).with_comment("tab\there"),
            Keyword::comment("caf\u{e9}"),
            Keyword::history("line\nbreak"),
            Keyword::integer("T\u{c9}MP", 1),
            Keyword::integer("TE\tMP", 1),
        ];
        for keyword in keywords {
            let result = format_card(&keyword);
            assert!(
                matches!(result, Err(FitsError::InvalidKeywordValue { .. })),
                "{:?} was accepted",
                keyword
            );
        }
    }

    #[test]
    fn failed_write_leaves_existing_file_untouched() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("existing.fits");
        let data: Vec<i16> = (0..4).collect();
        let nan = [Keyword::real("BAD", f64::NAN)];
        let non_ascii = [Keyword::string("OBJECT", "caf\u{e9}")];

        std::fs::write(&path, b"original contents").unwrap();
        let attempts: Vec<Result<()>> = vec![
            FitsWriter::create(&path)
                .unwrap()
                .write_primary_image(&data, &[2, 2], &nan),
            FitsWriter::create(&path)
                .unwrap()
                .write_primary_image_with_checksum(&data, &[2, 2], &non_ascii),
            FitsWriter::create(&path)
                .unwrap()
                .write_primary_image(&data, &[3, 3], &[]),
            FitsWriter::create(&path).unwrap().write_compressed_image(
                &data,
                &[2, 2],
                Some((2, 2)),
                CompressionAlgorithm::Rice,
                &nan,
            ),
        ];
        for (i, attempt) in attempts.iter().enumerate() {
            assert!(attempt.is_err(), "attempt {} succeeded", i);
        }
        let contents = std::fs::read(&path).unwrap();
        assert_eq!(
            String::from_utf8_lossy(&contents[..contents.len().min(40)]),
            "original contents"
        );

        let new_path = dir.path().join("new.fits");
        let result =
            FitsWriter::create(&new_path)
                .unwrap()
                .write_primary_image(&data, &[2, 2], &nan);
        assert!(result.is_err());
        assert!(!new_path.exists());
    }

    #[test]
    fn write_rejects_axis_above_i64_max() {
        let temp_file = NamedTempFile::new().unwrap();
        let mut writer = FitsWriter::create(temp_file.path()).unwrap();
        let empty: [i16; 0] = [];
        assert!(matches!(
            writer.write_primary_image(&empty, &[1usize << 63, 0], &[]),
            Err(FitsError::InvalidFormat(_))
        ));
        match writer.write_primary_image(&empty, &[1usize << 40, 1 << 40], &[]) {
            Err(FitsError::InvalidFormat(m)) => assert_eq!(
                m,
                "dimensions 1099511627776x1099511627776 overflow the pixel count"
            ),
            other => panic!("{other:?}"),
        }
    }

    fn standard_comments() -> [Keyword; 2] {
        [
            Keyword::comment(
                "FITS (Flexible Image Transport System) format is defined in 'Astronomy",
            ),
            Keyword::comment(
                "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H",
            ),
        ]
    }

    #[test]
    fn write_adds_standard_comments_only_when_absent() {
        let data: Vec<i16> = (0..4).collect();
        let mut keywords = standard_comments().to_vec();
        keywords.push(Keyword::comment("user comment"));
        for incoming in [&keywords[..], &[]] {
            let temp_file = NamedTempFile::new().unwrap();
            FitsWriter::create(temp_file.path())
                .unwrap()
                .write_primary_image(&data, &[2, 2], incoming)
                .unwrap();
            let (header, _) = written_hdu(temp_file.path());
            for standard in standard_comments() {
                let card = format!("COMMENT {:<72}", standard.comment.unwrap());
                let copies = header.chunks(80).filter(|c| *c == card.as_bytes()).count();
                assert_eq!(copies, 1, "{}", card);
            }
        }
    }

    #[test]
    fn format_card_writes_long_or_spaced_names_as_hierarch() {
        assert_eq!(
            card_text(&Keyword::real("ESO DET DIT", 1.5).with_comment("c")),
            padded(&["HIERARCH ESO DET DIT = 1.5 / c"])
        );
        assert_eq!(
            card_text(&Keyword::string("ESO OBS NAME", "O'Brien")),
            padded(&["HIERARCH ESO OBS NAME = 'O''Brien'"])
        );
        assert_eq!(
            card_text(&Keyword::integer("LONGKEYWORD", 3)),
            padded(&["HIERARCH LONGKEYWORD = 3"])
        );
        assert_eq!(
            card_text(&Keyword::logical("ESO INS SHUT", true)),
            padded(&["HIERARCH ESO INS SHUT = T"])
        );
    }

    #[test]
    fn format_card_hierarch_cuts_comment_and_rejects_long_name() {
        let comment = "x".repeat(80);
        let card = card_text(&Keyword::integer("ESO DET NDIT", 3).with_comment(&comment));
        assert_eq!(
            card,
            format!("HIERARCH ESO DET NDIT = 3 / {}", &comment[..52])
        );

        let long = Keyword::integer(format!("ESO {}", "Y".repeat(70)), 1);
        assert!(matches!(
            format_card(&long),
            Err(FitsError::InvalidKeywordValue { .. })
        ));
    }

    #[test]
    fn format_card_writes_non_standard_names_as_hierarch() {
        assert_eq!(
            card_text(&Keyword::integer("HIERARCH NAXIS1", 1)),
            padded(&["HIERARCH NAXIS1 = 1"])
        );
        assert_eq!(
            card_text(&Keyword::integer("foo", 1)),
            padded(&["HIERARCH foo = 1"])
        );
        assert_eq!(
            card_text(&Keyword::string("a.b", "x")),
            padded(&["HIERARCH a.b = 'x       '"])
        );
        assert_eq!(
            card_text(&Keyword::integer("DATE-OBS", 1)),
            padded(&["DATE-OBS=                    1"])
        );
    }

    #[test]
    fn format_card_rejects_names_that_cannot_be_hierarch() {
        for name in ["ESO=DET", " ESO", "ESO ", "HIERARCH ", "HIERARCH  ESO"] {
            assert!(
                matches!(
                    format_card(&Keyword::integer(name, 1)),
                    Err(FitsError::InvalidKeyword(_))
                ),
                "{:?}",
                name
            );
        }
    }

    #[test]
    fn format_card_continues_long_hierarch_string() {
        let value: String = ('a'..='z').cycle().take(108).collect();
        let keyword = Keyword::string("HIERARCH ESO OBS TARG NAME", &value).with_comment("c");
        let first = format!("HIERARCH ESO OBS TARG NAME = '{}&'", &value[..48]);
        let second = format!("CONTINUE  '{}' / c", &value[48..]);
        assert_eq!(card_text(&keyword), padded(&[&first, &second]));

        let mut block = card_text(&keyword).into_bytes();
        block.extend(format!("{:<80}", "END").into_bytes());
        block.resize(FITS_BLOCK_SIZE, b' ');
        let header = crate::fits::header::parser::HeaderParser::parse_header(&block).unwrap();
        let read = header.get_keyword("HIERARCH ESO OBS TARG NAME").unwrap();
        assert_eq!(read.value, Some(KeywordValue::String(value)));
        assert_eq!(read.comment.as_deref(), Some("c"));
    }

    #[test]
    fn format_card_writes_undefined_value_with_indicator() {
        assert_eq!(
            card_text(&Keyword::new("EMPTYC".to_string()).with_comment("only comment")),
            padded(&["EMPTYC  =                      / only comment"])
        );
        assert_eq!(
            card_text(&Keyword::new("EMPTY".to_string())),
            padded(&["EMPTY   ="])
        );
    }

    #[test]
    fn hierarch_cards_read_back_under_the_same_name() {
        use crate::fits::header::parser::HeaderParser;

        let keywords = [
            Keyword::real("HIERARCH ESO DET DIT", 1.5),
            Keyword::string("HIERARCH ESO OBS NAME", "a / b"),
            Keyword::integer("HIERARCH NAXIS1", 3),
        ];
        let mut block: Vec<u8> = keywords
            .iter()
            .map(card_text)
            .collect::<String>()
            .into_bytes();
        block.extend(format!("{:<80}", "END").into_bytes());
        block.resize(FITS_BLOCK_SIZE, b' ');
        let header = HeaderParser::parse_header(&block).unwrap();
        for keyword in &keywords {
            assert_eq!(header.get(&keyword.name).value(), keyword.value.as_ref());
        }
    }

    #[test]
    fn format_card_writes_blank_keyword_text_from_column_9() {
        let mut keyword = Keyword::comment("= blank keyword text");
        keyword.name = String::new();
        assert_eq!(
            card_text(&keyword),
            padded(&["        = blank keyword text"])
        );
    }

    fn parsed_keywords(cards: &[&str]) -> Vec<Keyword> {
        let mut block: Vec<u8> = cards
            .iter()
            .chain(&["END"])
            .flat_map(|c| format!("{:<80}", c).into_bytes())
            .collect();
        block.resize(FITS_BLOCK_SIZE, b' ');
        crate::fits::header::parser::HeaderParser::parse_header(&block)
            .unwrap()
            .keywords()
            .to_vec()
    }

    #[test]
    fn cards_without_value_indicator_are_written_back_unchanged() {
        let cards = [
            "KEYWORD  free text",
            "HIERARCHESO=1",
            "HIERARCH",
            "OBJECT  ='M31'",
            "CONTINUE  'orphan'",
            "lower    text",
            "EMPTY   =",
            "EMPTYC  =                      / only comment",
        ];
        let written: String = parsed_keywords(&cards).iter().map(card_text).collect();
        assert_eq!(written, padded(&cards));
    }

    #[test]
    fn format_card_moves_comment_of_long_string_to_final_continue() {
        let value: String = ('a'..='z').cycle().take(108).collect();
        let hierarch =
            Keyword::string("HIERARCH ESO OBS TARG NAME", &value).with_comment("a comment");
        let first = format!("HIERARCH ESO OBS TARG NAME = '{}&'", &value[..48]);
        let second = format!("CONTINUE  '{}&'", &value[48..]);
        assert_eq!(
            card_text(&hierarch),
            padded(&[&first, &second, "CONTINUE  '' / a comment"])
        );

        let long: String = ('a'..='z').cycle().take(127).collect();
        let standard = Keyword::string("LONGKEY", &long).with_comment("a comment");
        let first = format!("LONGKEY = '{}&'", &long[..67]);
        let second = format!("CONTINUE  '{}&'", &long[67..]);
        assert_eq!(
            card_text(&standard),
            padded(&[&first, &second, "CONTINUE  '' / a comment"])
        );

        for keyword in [hierarch, standard] {
            let text = card_text(&keyword);
            let lines: Vec<&str> = (0..text.len())
                .step_by(80)
                .map(|i| &text[i..i + 80])
                .collect();
            let read = parsed_keywords(&lines);
            assert_eq!(read, [keyword]);
        }
    }

    #[test]
    fn format_card_writes_commentary_keyword_with_value_as_value_card() {
        let keyword = Keyword::commentary("KEY", "text").with_value(5);
        assert_eq!(
            card_text(&keyword),
            padded(&["KEY     =                    5 / text"])
        );
    }
}
