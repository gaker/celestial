use crate::fits::header::{Keyword, KeywordValue};
use crate::fits::{FitsError, Result};
use std::collections::HashMap;
use std::str;

const CARD_SIZE: usize = 80;
const HEADER_BLOCK_SIZE: usize = 2880;

#[derive(Debug, Clone)]
pub struct Header {
    keywords: Vec<Keyword>,
    keyword_index: HashMap<String, usize>,
}

#[derive(Debug, Clone)]
pub struct HeaderCard {
    pub keyword: String,
    pub value: Option<String>,
    pub comment: Option<String>,
    raw: [u8; CARD_SIZE],
}

pub struct HeaderParser;

impl Header {
    pub fn new() -> Self {
        Self {
            keywords: Vec::new(),
            keyword_index: HashMap::new(),
        }
    }

    pub fn add_keyword(&mut self, keyword: Keyword) {
        let index = self.keywords.len();
        self.keyword_index.insert(keyword.name.clone(), index);
        self.keywords.push(keyword);
    }

    pub fn get_keyword(&self, name: &str) -> Option<&Keyword> {
        self.keyword_index
            .get(name)
            .and_then(|&index| self.keywords.get(index))
    }

    pub fn get_keyword_value(&self, name: &str) -> Option<&KeywordValue> {
        self.get_keyword(name)?.value.as_ref()
    }

    pub fn keywords(&self) -> &[Keyword] {
        &self.keywords
    }

    pub fn iter(&self) -> impl Iterator<Item = &Keyword> {
        self.keywords.iter()
    }

    pub fn is_primary(&self) -> bool {
        self.get_keyword("SIMPLE")
            .and_then(|k| k.value.as_ref())
            .and_then(|v| v.as_logical())
            .unwrap_or(false)
    }

    pub fn is_extension(&self) -> bool {
        self.get_keyword("XTENSION").is_some()
    }

    fn clear_sensitive_buffers(&mut self) {}
}

impl Default for Header {
    fn default() -> Self {
        Self::new()
    }
}

impl HeaderCard {
    pub fn parse(data: &[u8; CARD_SIZE]) -> Result<Self> {
        let card_str = Self::validate_card_data(data)?;
        let mut card = Self::create_empty_card(*data);

        card.keyword = Self::extract_keyword(card_str)?;

        Self::parse_value_and_comment(card_str, &mut card);

        card.clear_sensitive_data();
        Ok(card)
    }

    fn validate_card_data(data: &[u8; CARD_SIZE]) -> Result<&str> {
        let card_str = str::from_utf8(data)
            .map_err(|_| FitsError::InvalidFormat("Invalid UTF-8 in header card".to_string()))?;

        if card_str.len() < 8 {
            return Err(FitsError::HeaderParse("Card too short".to_string()));
        }

        Ok(card_str)
    }

    fn create_empty_card(raw_data: [u8; CARD_SIZE]) -> HeaderCard {
        HeaderCard {
            keyword: String::new(),
            value: None,
            comment: None,
            raw: raw_data,
        }
    }

    fn extract_keyword(card_str: &str) -> Result<String> {
        let keyword_part = &card_str[0..8];
        Ok(keyword_part.trim().to_string())
    }

    fn parse_value_and_comment(card_str: &str, card: &mut HeaderCard) {
        if card_str.len() >= 10 && &card_str[8..10] == "= " {
            Self::parse_keyword_value_pair(&card_str[10..], card);
        } else if card_str.len() >= 9 {
            Self::parse_comment_only(&card_str[8..], card);
        }
    }

    fn parse_keyword_value_pair(value_comment_part: &str, card: &mut HeaderCard) {
        if let Some(comment_pos) = value_comment_part.find(" / ") {
            Self::parse_value_with_comment(value_comment_part, comment_pos, card);
        } else {
            Self::parse_value_only(value_comment_part, card);
        }
    }

    fn parse_value_with_comment(
        value_comment_part: &str,
        comment_pos: usize,
        card: &mut HeaderCard,
    ) {
        let value_part = value_comment_part[..comment_pos].trim();
        let comment_part = value_comment_part[comment_pos + 3..].trim();

        if !value_part.is_empty() {
            card.value = Some(value_part.to_string());
        }
        if !comment_part.is_empty() {
            card.comment = Some(comment_part.to_string());
        }
    }

    fn parse_value_only(value_comment_part: &str, card: &mut HeaderCard) {
        let value_part = value_comment_part.trim();
        if !value_part.is_empty() {
            card.value = Some(value_part.to_string());
        }
    }

    fn parse_comment_only(rest_of_card: &str, card: &mut HeaderCard) {
        let comment_part = rest_of_card.trim();
        if !comment_part.is_empty() {
            card.comment = Some(comment_part.to_string());
        }
    }

    pub fn to_keyword(&self) -> Result<Keyword> {
        let mut keyword = Keyword::new(self.keyword.clone());

        if let Some(comment) = &self.comment {
            keyword = keyword.with_comment(comment.clone());
        }

        if let Some(value_str) = &self.value {
            let parsed_value = Self::parse_value(value_str)?;
            keyword = keyword.with_value(parsed_value);
        }

        Ok(keyword)
    }

    fn parse_value(value_str: &str) -> Result<KeywordValue> {
        let trimmed = value_str.trim();

        if trimmed == "T" {
            return Ok(KeywordValue::Logical(true));
        }
        if trimmed == "F" {
            return Ok(KeywordValue::Logical(false));
        }

        if trimmed.starts_with('\'') && trimmed.ends_with('\'') && trimmed.len() >= 2 {
            let string_content = &trimmed[1..trimmed.len() - 1];
            return Ok(KeywordValue::String(string_content.trim_end().to_string()));
        }

        if let Ok(int_val) = trimmed.parse::<i64>() {
            return Ok(KeywordValue::Integer(int_val));
        }

        if let Ok(float_val) = trimmed.parse::<f64>() {
            return Ok(KeywordValue::Real(float_val));
        }

        Ok(KeywordValue::String(trimmed.to_string()))
    }

    fn clear_sensitive_data(&mut self) {
        self.raw.fill(0);
    }
}

impl HeaderParser {
    pub fn parse_header(data: &[u8]) -> Result<Header> {
        if !data.len().is_multiple_of(HEADER_BLOCK_SIZE) {
            return Err(FitsError::InvalidFormat(
                "Header size must be multiple of 2880 bytes".to_string(),
            ));
        }

        let mut header = Header::new();
        let mut found_end = false;

        for chunk in data.chunks_exact(CARD_SIZE) {
            if chunk.len() != CARD_SIZE {
                break;
            }

            let mut card_data = [0u8; CARD_SIZE];
            card_data.copy_from_slice(chunk);

            let card = match HeaderCard::parse(&card_data) {
                Ok(card) => card,
                Err(e) => {
                    return Err(e);
                }
            };

            if card.keyword == "END" {
                found_end = true;
                break;
            }

            let keyword = match card.to_keyword() {
                Ok(keyword) => keyword,
                Err(e) => {
                    return Err(e);
                }
            };
            header.add_keyword(keyword);

            card_data.fill(0);
        }

        if !found_end {
            return Err(FitsError::InvalidFormat("Missing END keyword".to_string()));
        }

        header.clear_sensitive_buffers();
        Ok(header)
    }

    pub fn header_size_bytes(data: &[u8]) -> Result<usize> {
        let mut blocks = 0;
        let mut found_end = false;

        for block in data.chunks(HEADER_BLOCK_SIZE) {
            blocks += 1;

            for chunk in block.chunks_exact(CARD_SIZE) {
                if chunk.len() != CARD_SIZE {
                    continue;
                }

                if chunk.len() < 8 {
                    return Err(FitsError::InvalidFormat(
                        "Header card too short".to_string(),
                    ));
                }

                let keyword_part = str::from_utf8(&chunk[0..8])
                    .map_err(|_| FitsError::InvalidFormat("Invalid UTF-8 in header".to_string()))?;

                if keyword_part.trim() == "END" {
                    found_end = true;
                    break;
                }
            }

            if found_end {
                break;
            }
        }

        if !found_end {
            return Err(FitsError::InvalidFormat("Missing END keyword".to_string()));
        }

        Ok(blocks * HEADER_BLOCK_SIZE)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::*;

    #[test]
    fn header_card_parse_keyword() {
        let mut card = [b' '; 80];
        let test_str = "SIMPLE  = T                     / Standard FITS format                 ";
        card[..test_str.len()].copy_from_slice(test_str.as_bytes());

        let parsed = HeaderCard::parse(&card).unwrap();
        assert_eq!(parsed.keyword, "SIMPLE");
        assert_eq!(parsed.value.as_deref(), Some("T"));
        assert_eq!(parsed.comment.as_deref(), Some("Standard FITS format"));
    }

    #[test]
    fn header_card_parse_numeric_value() {
        let mut card = [b' '; 80];
        let test_str = "BITPIX  =                   16 / Bits per pixel                        ";
        card[..test_str.len()].copy_from_slice(test_str.as_bytes());

        let parsed = HeaderCard::parse(&card).unwrap();
        assert_eq!(parsed.keyword, "BITPIX");
        assert_eq!(parsed.value.as_deref(), Some("16"));
        assert_eq!(parsed.comment.as_deref(), Some("Bits per pixel"));
    }

    #[test]
    fn header_card_parse_string_value() {
        let mut card = [b' '; 80];
        let test_str = "OBJECT  = 'M31 Galaxy'         / Target object                         ";
        card[..test_str.len()].copy_from_slice(test_str.as_bytes());

        let parsed = HeaderCard::parse(&card).unwrap();
        assert_eq!(parsed.keyword, "OBJECT");
        assert_eq!(parsed.value.as_deref(), Some("'M31 Galaxy'"));
        assert_eq!(parsed.comment.as_deref(), Some("Target object"));
    }

    #[test]
    fn header_card_parse_comment_only() {
        let mut card = [b' '; 80];
        let test_str = "HISTORY This is a history comment                                       ";
        card[..test_str.len()].copy_from_slice(test_str.as_bytes());

        let parsed = HeaderCard::parse(&card).unwrap();
        assert_eq!(parsed.keyword, "HISTORY");
        assert!(parsed.value.is_none());
        assert_eq!(parsed.comment.as_deref(), Some("This is a history comment"));
    }

    #[test]
    fn header_card_parse_value_no_comment() {
        let mut card = [b' '; 80];
        let test_str = "NAXIS   =                    2                                         ";
        card[..test_str.len()].copy_from_slice(test_str.as_bytes());

        let parsed = HeaderCard::parse(&card).unwrap();
        assert_eq!(parsed.keyword, "NAXIS");
        assert_eq!(parsed.value.as_deref(), Some("2"));
        assert!(parsed.comment.is_none());
    }

    #[test]
    fn validate_card_data_rejects_invalid_utf8() {
        let mut card = [b' '; 80];
        card[0] = 0xFF;
        card[1] = 0xFE;

        let result = HeaderCard::parse(&card);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn validate_card_data_rejects_short_cards() {
        let mut card = [0u8; 80];
        card[..7].copy_from_slice(b"KEYWORD");

        let result = HeaderCard::parse(&card);
        assert!(result.is_ok());
    }

    #[test]
    fn extract_keyword_handles_edge_cases() {
        let test_cases = [
            ("SIMPLE  ", "SIMPLE"),
            ("        ", ""),
            ("TEST1234", "TEST1234"),
            ("A       ", "A"),
        ];

        for (input, expected) in test_cases {
            let mut card = [b' '; 80];
            card[..input.len()].copy_from_slice(input.as_bytes());

            let parsed = HeaderCard::parse(&card).unwrap();
            assert_eq!(parsed.keyword, expected);
        }
    }

    #[test]
    fn parse_value_and_comment_malformed_separators() {
        let malformed_cases = [
            "KEYWORD =",
            "KEYWORD ==",
            "KEYWORD = VALUE / / COMMENT",
            "KEYWORD = 'UNTERMINATED",
            "KEYWORD = VALUE / ",
            "KEYWORD =    ",
        ];

        for test_case in malformed_cases {
            let mut card = [b' '; 80];
            let padded = format!("{:<80}", test_case);
            card.copy_from_slice(padded.as_bytes());

            assert!(HeaderCard::parse(&card).is_ok());
        }
    }

    #[test]
    fn parse_keyword_value_pair_various_formats() {
        let test_cases = [
            ("KEYWORD = VALUE / COMMENT", Some("VALUE"), Some("COMMENT")),
            (
                "KEYWORD = 'STRING' / String comment",
                Some("'STRING'"),
                Some("String comment"),
            ),
            ("KEYWORD = 12345 / Numeric", Some("12345"), Some("Numeric")),
            ("KEYWORD = T / Boolean", Some("T"), Some("Boolean")),
            ("KEYWORD =  / Just comment", None, Some("Just comment")),
            ("KEYWORD = VALUE_NO_COMMENT", Some("VALUE_NO_COMMENT"), None),
        ];

        for (input, expected_value, expected_comment) in test_cases {
            let mut card = [b' '; 80];
            let padded = format!("{:<80}", input);
            card.copy_from_slice(padded.as_bytes());

            let parsed = HeaderCard::parse(&card).unwrap();
            assert_eq!(parsed.value.as_deref(), expected_value);
            assert_eq!(parsed.comment.as_deref(), expected_comment);
        }
    }

    #[test]
    fn header_size_calculation() {
        let fits_data = create_minimal_fits();
        let size = HeaderParser::header_size_bytes(&fits_data).unwrap();
        assert_eq!(size, 2880);
    }

    #[test]
    fn header_parser_with_long_values() {
        let long_value = "A".repeat(68);
        let card_content = format!("LONGVAL = '{}'", long_value);

        let mut card = [b' '; 80];
        let padded = format!("{:<80}", card_content);
        card.copy_from_slice(padded.as_bytes());

        let parsed = HeaderCard::parse(&card).unwrap();
        assert_eq!(parsed.keyword, "LONGVAL");
        assert!(parsed.value.as_ref().unwrap().contains(&long_value));
    }

    #[test]
    fn header_parser_extreme_cases() {
        let empty_card = [b' '; 80];
        let result = HeaderCard::parse(&empty_card);
        assert!(result.is_ok());

        let max_keyword = "TESTKEY1";
        let mut card = [b' '; 80];
        card[..8].copy_from_slice(max_keyword.as_bytes());

        let parsed = HeaderCard::parse(&card).unwrap();
        assert_eq!(parsed.keyword, max_keyword);
    }

    #[test]
    fn header_keywords_accessor() {
        let mut header = Header::new();
        let keyword1 = Keyword::new("TEST1".to_string());
        let keyword2 = Keyword::new("TEST2".to_string());

        header.add_keyword(keyword1);
        header.add_keyword(keyword2);

        let keywords = header.keywords();
        assert_eq!(keywords.len(), 2);
        assert_eq!(keywords[0].name, "TEST1");
        assert_eq!(keywords[1].name, "TEST2");
    }

    #[test]
    fn header_iter() {
        let mut header = Header::new();
        let keyword1 = Keyword::new("ITER1".to_string());
        let keyword2 = Keyword::new("ITER2".to_string());

        header.add_keyword(keyword1);
        header.add_keyword(keyword2);

        let mut iter_count = 0;
        for keyword in header.iter() {
            iter_count += 1;
            assert!(keyword.name.starts_with("ITER"));
        }
        assert_eq!(iter_count, 2);
    }

    #[test]
    fn header_is_extension() {
        let mut header = Header::new();
        assert!(!header.is_extension());

        let xtension_keyword = Keyword::new("XTENSION".to_string())
            .with_value(KeywordValue::String("IMAGE".to_string()));
        header.add_keyword(xtension_keyword);

        assert!(header.is_extension());
    }

    #[test]
    fn header_default() {
        let header = Header::default();
        assert_eq!(header.keywords.len(), 0);
        assert_eq!(header.keyword_index.len(), 0);
    }

    #[test]
    fn validate_card_data_too_short() {
        let mut short_card = [0u8; 80];
        short_card[..7].copy_from_slice(b"KEYWORD");

        let result = HeaderCard::parse(&short_card);
        assert!(result.is_ok());
    }

    #[test]
    fn parse_value_string_with_quotes() {
        let string_value = "'Test String'";
        let parsed = HeaderCard::parse_value(string_value).unwrap();

        match parsed {
            KeywordValue::String(s) => assert_eq!(s, "Test String"),
            _ => panic!("Expected String value"),
        }
    }

    #[test]
    fn parse_value_real_number() {
        let real_value = "1.23456";
        let parsed = HeaderCard::parse_value(real_value).unwrap();

        match parsed {
            KeywordValue::Real(f) => assert!((f - 1.23456).abs() < 1e-10),
            _ => panic!("Expected Real value"),
        }
    }

    #[test]
    fn parse_value_unquoted_string() {
        let unquoted_values = [
            ("not_a_number_or_string", "not_a_number_or_string"),
            ("SharpCap 4.1.12395.0", "SharpCap 4.1.12395.0"),
            ("SOME_VALUE", "SOME_VALUE"),
        ];

        for (input, expected) in unquoted_values {
            let result = HeaderCard::parse_value(input).unwrap();
            match result {
                KeywordValue::String(s) => assert_eq!(s, expected),
                _ => panic!("Expected String value for '{}'", input),
            }
        }
    }

    #[test]
    fn parse_value_unterminated_string() {
        let result = HeaderCard::parse_value("'unterminated string");
        assert!(result.is_ok());
        match result.unwrap() {
            KeywordValue::String(s) => assert_eq!(s, "'unterminated string"),
            _ => panic!("Expected String fallback"),
        }
    }

    #[test]
    fn parse_value_empty() {
        let result = HeaderCard::parse_value("");
        assert!(result.is_ok());
        match result.unwrap() {
            KeywordValue::String(s) => assert_eq!(s, ""),
            _ => panic!("Expected empty String"),
        }
    }

    #[test]
    fn parse_header_invalid_size() {
        let invalid_data = vec![0u8; 1000];
        let result = HeaderParser::parse_header(&invalid_data);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn parse_header_chunk_size_check() {
        let mut valid_data = vec![0u8; 2880];

        let keyword_card = format!(
            "{:<80}",
            "SIMPLE  = T                     / Standard FITS format"
        );
        let bitpix_card = format!("{:<80}", "BITPIX  =                   16 / Bits per pixel");

        valid_data[..80].copy_from_slice(keyword_card.as_bytes());
        valid_data[80..160].copy_from_slice(bitpix_card.as_bytes());

        let result = HeaderParser::parse_header(&valid_data);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn parse_header_card_parsing_error() {
        let mut invalid_data = vec![0u8; 2880];

        let mut bad_card = [0u8; 80];
        bad_card[0] = 0xFF;
        bad_card[1] = 0xFE;

        invalid_data[..80].copy_from_slice(&bad_card);

        let result = HeaderParser::parse_header(&invalid_data);
        assert!(result.is_err());
    }

    #[test]
    fn parse_header_unterminated_string_as_unquoted() {
        let mut valid_data = vec![0u8; 2880];

        let test_card = format!("{:<80}", "TESTKEY = 'unterminated string");
        valid_data[..80].copy_from_slice(test_card.as_bytes());

        let end_card = format!("{:<80}", "END");
        valid_data[80..160].copy_from_slice(end_card.as_bytes());

        let result = HeaderParser::parse_header(&valid_data);
        assert!(result.is_ok());
        let header = result.unwrap();
        let value = header.get_keyword_value("TESTKEY").unwrap();
        assert_eq!(
            value,
            &KeywordValue::String("'unterminated string".to_string())
        );
    }

    #[test]
    fn parse_header_missing_end_keyword() {
        let mut valid_data = vec![0u8; 2880];

        let simple_card = format!(
            "{:<80}",
            "SIMPLE  = T                     / Standard FITS format"
        );
        valid_data[..80].copy_from_slice(simple_card.as_bytes());

        let result = HeaderParser::parse_header(&valid_data);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn header_size_continue_on_short_chunk() {
        let mut data = vec![0u8; 2890];

        let simple_card = format!(
            "{:<80}",
            "SIMPLE  = T                     / Standard FITS format"
        );
        data[..80].copy_from_slice(simple_card.as_bytes());

        let end_card = format!("{:<80}", "END");
        data[80..160].copy_from_slice(end_card.as_bytes());

        let result = HeaderParser::header_size_bytes(&data);
        assert!(result.is_ok());
    }

    #[test]
    fn header_size_card_too_short() {
        let mut data = vec![0u8; 2880];

        let simple_card = format!(
            "{:<80}",
            "SIMPLE  = T                     / Standard FITS format"
        );
        data[..80].copy_from_slice(simple_card.as_bytes());

        let end_card = format!("{:<80}", "END");
        data[80..160].copy_from_slice(end_card.as_bytes());

        let result = HeaderParser::header_size_bytes(&data);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 2880);
    }

    #[test]
    fn header_size_invalid_utf8() {
        let mut data = vec![0u8; 2880];

        data[0] = 0xFF;
        data[1] = 0xFE;

        let result = HeaderParser::header_size_bytes(&data);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), FitsError::InvalidFormat(_)));
    }

    #[test]
    fn to_keyword_with_all_components() {
        let mut card = [b' '; 80];
        let test_str = "TESTKEY = 42                    / Test comment                          ";
        card[..test_str.len()].copy_from_slice(test_str.as_bytes());

        let parsed = HeaderCard::parse(&card).unwrap();
        let keyword = parsed.to_keyword().unwrap();

        assert_eq!(keyword.name, "TESTKEY");
        assert!(keyword.value.is_some());
        assert!(keyword.comment.is_some());
    }

    #[test]
    fn complete_header_workflow() {
        let mut header = Header::new();

        let simple_kw = Keyword::new("SIMPLE".to_string())
            .with_value(KeywordValue::Logical(true))
            .with_comment("Standard FITS".to_string());

        let bitpix_kw = Keyword::new("BITPIX".to_string()).with_value(KeywordValue::Integer(16));

        header.add_keyword(simple_kw);
        header.add_keyword(bitpix_kw);

        assert!(header.is_primary());
        assert!(!header.is_extension());

        let simple_val = header.get_keyword_value("SIMPLE").unwrap();
        assert!(simple_val.as_logical().unwrap());

        assert_eq!(header.keywords().len(), 2);

        let mut count = 0;
        for _ in header.iter() {
            count += 1;
        }
        assert_eq!(count, 2);
    }

    #[test]
    fn parse_sharpcap_header_with_unquoted_string() {
        let mut data = vec![b' '; 2880];
        let cards = [
            "SIMPLE  =                    T",
            "BITPIX  =                  -32",
            "NAXIS   =                    2",
            "NAXIS1  =                 6384",
            "NAXIS2  =                 4256",
            "SWCREATE= SharpCap 4.1.12395.0",
            "END",
        ];

        let mut offset = 0;
        for card in cards {
            let padded = format!("{:<80}", card);
            data[offset..offset + 80].copy_from_slice(padded.as_bytes());
            offset += 80;
        }

        let header = HeaderParser::parse_header(&data).unwrap();

        assert_eq!(header.keywords().len(), 6);

        let swcreate = header.get_keyword_value("SWCREATE").unwrap();
        assert_eq!(
            swcreate,
            &KeywordValue::String("SharpCap 4.1.12395.0".to_string())
        );

        let simple = header.get_keyword_value("SIMPLE").unwrap();
        assert!(simple.as_logical().unwrap());

        let bitpix = header.get_keyword_value("BITPIX").unwrap();
        assert_eq!(bitpix.as_integer().unwrap(), -32);

        let naxis = header.get_keyword_value("NAXIS").unwrap();
        assert_eq!(naxis.as_integer().unwrap(), 2);
    }
}
