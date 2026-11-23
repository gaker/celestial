use crate::fits::block::{CARD_SIZE, FITS_BLOCK_SIZE};
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::keywords::{find_last_keyword, hierarch_alias, Keyword, KeywordValue};
use std::collections::HashMap;
use std::str;

mod card;
mod value;

#[derive(Debug, Clone)]
pub struct Header {
    keywords: Vec<Keyword>,
    keyword_index: HashMap<String, usize>,
}

#[derive(Debug, Clone)]
pub(crate) struct HeaderCard {
    pub keyword: String,
    pub value: Option<String>,
    pub comment: Option<String>,
    commentary: bool,
}

pub(crate) struct HeaderParser;

impl Header {
    pub fn new() -> Self {
        Self {
            keywords: Vec::new(),
            keyword_index: HashMap::new(),
        }
    }

    // A repeated name keeps the index of its first card, which is the card
    // astropy and cfitsio return for that name.
    pub(crate) fn add_keyword(&mut self, keyword: Keyword) {
        let index = self.keywords.len();
        self.keyword_index
            .entry(keyword.name.clone())
            .or_insert(index);
        self.keywords.push(keyword);
    }

    pub fn get_keyword(&self, name: &str) -> Option<&Keyword> {
        let index = self
            .keyword_index
            .get(name)
            .or_else(|| hierarch_alias(name).and_then(|full| self.keyword_index.get(&full)))?;
        self.keywords.get(*index)
    }

    pub fn get<'a>(&'a self, name: &str) -> Field<'a> {
        Field::from_value(self.get_keyword(name).and_then(|k| k.value.as_ref()))
    }

    pub(crate) fn get_last_keyword(&self, name: &str) -> Option<&Keyword> {
        self.keywords.get(find_last_keyword(&self.keywords, name)?)
    }

    // The value of the last card named `name`, which data and WCS readers
    // apply when a scaling, null, table layout or WCS keyword repeats.
    pub fn get_last<'a>(&'a self, name: &str) -> Field<'a> {
        Field::from_value(self.get_last_keyword(name).and_then(|k| k.value.as_ref()))
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
}

impl Default for Header {
    fn default() -> Self {
        Self::new()
    }
}

pub struct Field<'a> {
    value: Option<&'a KeywordValue>,
}

impl<'a> Field<'a> {
    pub(crate) fn from_value(value: Option<&'a KeywordValue>) -> Self {
        Self { value }
    }

    pub fn exists(&self) -> bool {
        self.value.is_some()
    }

    pub fn value(&self) -> Option<&'a KeywordValue> {
        self.value
    }

    pub fn as_bool(&self) -> Option<bool> {
        self.value?.as_logical()
    }

    pub fn as_i64(&self) -> Option<i64> {
        self.value?.as_integer()
    }

    pub fn as_f64(&self) -> Option<f64> {
        self.value?.as_real()
    }

    pub fn as_str(&self) -> Option<&'a str> {
        self.value?.as_string()
    }
}

impl HeaderParser {
    pub(crate) fn parse_header(data: &[u8]) -> Result<Header> {
        if !data.len().is_multiple_of(FITS_BLOCK_SIZE) {
            return Err(FitsError::InvalidFormat(
                "Header size must be multiple of 2880 bytes".to_string(),
            ));
        }

        let mut header = Header::new();
        let mut last_is_quoted = false;

        for chunk in data.as_chunks::<CARD_SIZE>().0 {
            if Self::is_end_card(chunk)? {
                return Ok(header);
            }

            let text = HeaderCard::printable_card(chunk);
            let card = HeaderCard::parse_printable(&text)?;

            if card.keyword == "CONTINUE"
                && last_is_quoted
                && Self::append_continued_string(&mut header, &text)?
            {
                continue;
            }
            last_is_quoted = card
                .value
                .as_deref()
                .is_some_and(HeaderCard::is_quoted_value);
            header.add_keyword(card.to_keyword()?);
        }

        Err(FitsError::InvalidFormat("Missing END keyword".to_string()))
    }

    // FITS 4.0 section 4.4.1 requires spaces in columns 9-80. Other text makes the
    // end of the header ambiguous, so it is an error rather than a guess.
    fn is_end_card(card: &[u8; CARD_SIZE]) -> Result<bool> {
        if card[..8].trim_ascii() != b"END" {
            return Ok(false);
        }
        if card[8..].iter().any(|&b| b != b' ') {
            return Err(FitsError::InvalidFormat(
                "END card holds text after the keyword".to_string(),
            ));
        }
        Ok(true)
    }

    // Long-string convention (FITS 4.0 §4.2.1.2): a quoted string ending in
    // '&' is continued by the next CONTINUE card when that card holds one
    // quoted string and an optional comment. Returns false otherwise, so the
    // '&' stays literal and the CONTINUE card is kept as commentary.
    // The caller checks that the previous value was written in quotes, and
    // passes the card after `HeaderCard::printable_card`.
    fn append_continued_string(header: &mut Header, card: &[u8; CARD_SIZE]) -> Result<bool> {
        let Some(previous) = header.keywords.last_mut() else {
            return Ok(false);
        };
        let Some(KeywordValue::String(text)) = previous.value.as_mut() else {
            return Ok(false);
        };
        if !text.ends_with('&') {
            return Ok(false);
        }
        let Some((next, comment)) = Self::continuation_value(card)? else {
            return Ok(false);
        };
        text.pop();
        text.push_str(&next);
        // Trailing spaces are not significant, as for a single-card string.
        text.truncate(text.trim_end().len());
        if let Some(comment) = comment {
            previous.comment = Some(match previous.comment.take() {
                Some(existing) => format!("{} {}", existing, comment),
                None => comment,
            });
        }
        Ok(true)
    }

    // Returns the quoted string and comment of a CONTINUE card, or `None`
    // when the card does not hold a single quoted string.
    fn continuation_value(card: &[u8; CARD_SIZE]) -> Result<Option<(String, Option<String>)>> {
        let card_str = HeaderCard::card_str(card)?;
        let Some(value_text) = card_str.get(10..) else {
            return Ok(None);
        };
        let Some((quoted, rest)) = HeaderCard::split_quoted(value_text) else {
            return Ok(None);
        };
        let rest = rest.trim();
        let comment = match rest.strip_prefix('/') {
            Some(comment) => Some(comment.trim()).filter(|c| !c.is_empty()),
            None if rest.is_empty() => None,
            None => return Ok(None),
        };
        let KeywordValue::String(next) = HeaderCard::parse_value(quoted)? else {
            return Ok(None);
        };
        Ok(Some((next, comment.map(str::to_string))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn validate_card_data_replaces_invalid_utf8() {
        let mut card = [b' '; 80];
        card[0] = 0xFF;
        card[1] = 0xFE;

        let result = HeaderCard::parse(&card);
        assert_eq!(result.map(|c| c.keyword).ok().as_deref(), Some("??"));
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

    fn typed_header() -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 2));
        header.add_keyword(Keyword::real("EXPTIME", 30.0));
        header.add_keyword(Keyword::string("OBJECT", "M31"));
        header
    }

    #[test]
    fn field_typed_accessors_match_value() {
        let header = typed_header();
        assert_eq!(header.get("SIMPLE").as_bool(), Some(true));
        assert_eq!(header.get("NAXIS").as_i64(), Some(2));
        assert_eq!(header.get("EXPTIME").as_f64(), Some(30.0));
        assert_eq!(header.get("OBJECT").as_str(), Some("M31"));
    }

    #[test]
    fn field_as_f64_coalesces_integer_card() {
        let header = typed_header();
        assert_eq!(header.get("NAXIS").as_f64(), Some(2.0));
    }

    #[test]
    fn field_as_i64_stays_strict_on_real_card() {
        let header = typed_header();
        assert_eq!(header.get("EXPTIME").as_i64(), None);
    }

    #[test]
    fn field_present_but_wrong_type_returns_none() {
        let header = typed_header();
        assert_eq!(header.get("OBJECT").as_i64(), None);
        assert_eq!(header.get("OBJECT").as_f64(), None);
        assert_eq!(header.get("OBJECT").as_bool(), None);
        assert_eq!(header.get("NAXIS").as_str(), None);
        assert_eq!(header.get("NAXIS").as_bool(), None);
    }

    #[test]
    fn field_absent_keyword_returns_none_and_not_exists() {
        let header = typed_header();
        let missing = header.get("MISSING");
        assert!(!missing.exists());
        assert_eq!(missing.as_i64(), None);
        assert_eq!(missing.as_f64(), None);
        assert_eq!(missing.as_str(), None);
        assert_eq!(missing.as_bool(), None);
        assert_eq!(missing.value(), None);
    }

    #[test]
    fn field_exists_distinguishes_present_from_absent() {
        let header = typed_header();
        assert!(header.get("OBJECT").exists());
        assert!(!header.get("MISSING").exists());
    }

    #[test]
    fn field_value_escape_hatch_returns_raw() {
        let header = typed_header();
        assert_eq!(
            header.get("OBJECT").value(),
            Some(&KeywordValue::String("M31".to_string()))
        );
        assert_eq!(header.get("NAXIS").value(), Some(&KeywordValue::Integer(2)));
    }

    #[test]
    fn field_valueless_keyword_does_not_exist() {
        let mut header = Header::new();
        header.add_keyword(Keyword::comment("just a comment"));
        let field = header.get("COMMENT");
        assert!(!field.exists());
        assert_eq!(field.value(), None);
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
            KeywordValue::Real(f) => assert_eq!(f, 1.23456),
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

    fn header_block(cards: &[&str]) -> Vec<u8> {
        let mut block: Vec<u8> = cards
            .iter()
            .chain(std::iter::once(&"END"))
            .flat_map(|c| format!("{:<80}", c).into_bytes())
            .collect();
        block.resize(FITS_BLOCK_SIZE, b' ');
        block
    }

    fn string_value(header: &Header, name: &str) -> Option<String> {
        header.get(name).as_str().map(str::to_string)
    }

    #[test]
    fn parse_header_unescapes_doubled_quotes() {
        let header = HeaderParser::parse_header(&header_block(&[
            "OBSERVER= 'O''Brien'",
            "QUOTES  = ''''''''''",
            "OBJECT  = 'M31 / Andromeda' / slash inside",
        ]))
        .unwrap();
        assert_eq!(
            string_value(&header, "OBSERVER").as_deref(),
            Some("O'Brien")
        );
        assert_eq!(string_value(&header, "QUOTES").as_deref(), Some("''''"));
        assert_eq!(
            string_value(&header, "OBJECT").as_deref(),
            Some("M31 / Andromeda")
        );
        assert_eq!(
            header.get_keyword("OBJECT").unwrap().comment.as_deref(),
            Some("slash inside")
        );
    }

    #[test]
    fn parse_header_joins_continue_cards() {
        let header = HeaderParser::parse_header(&header_block(&[
            "LONGSTR = 'abc  &'",
            "CONTINUE  'd''e&'",
            "CONTINUE  'fgh' / the comment",
            "SHORT   = 'no ampersand'",
            "CONTINUE  'orphan'",
        ]))
        .unwrap();
        assert_eq!(
            string_value(&header, "LONGSTR").as_deref(),
            Some("abc  d'efgh")
        );
        assert_eq!(
            header.get_keyword("LONGSTR").unwrap().comment.as_deref(),
            Some("the comment")
        );
        assert_eq!(
            string_value(&header, "SHORT").as_deref(),
            Some("no ampersand")
        );
        let names: Vec<&str> = header.keywords().iter().map(|k| k.name.as_str()).collect();
        assert_eq!(names, ["LONGSTR", "SHORT", "CONTINUE"]);
    }

    #[test]
    fn parse_header_keeps_ampersand_when_continue_does_not_conform() {
        let header = HeaderParser::parse_header(&header_block(&[
            "LS4     = 'abc&'",
            "CONTINUE  plain text no quotes",
            "LS5     = 'abc&'",
            "CONTINUE  'x' junk after",
            "UNQ     = abc&",
            "CONTINUE  'xyz'",
        ]))
        .unwrap();
        for name in ["LS4", "LS5", "UNQ"] {
            assert_eq!(
                string_value(&header, name).as_deref(),
                Some("abc&"),
                "{}",
                name
            );
        }
        let names: Vec<&str> = header.keywords().iter().map(|k| k.name.as_str()).collect();
        assert_eq!(
            names,
            ["LS4", "CONTINUE", "LS5", "CONTINUE", "UNQ", "CONTINUE"]
        );
    }

    #[test]
    fn parse_header_trims_joined_long_string() {
        let header = HeaderParser::parse_header(&header_block(&[
            "NOTES   = 'abc   &'",
            "CONTINUE  '' / note",
        ]))
        .unwrap();
        assert_eq!(string_value(&header, "NOTES").as_deref(), Some("abc"));
    }

    #[test]
    fn parse_header_rejects_text_after_end_keyword() {
        let mut block = format!("{:<80}", "END                 0").into_bytes();
        block.resize(FITS_BLOCK_SIZE, b' ');
        let result = HeaderParser::parse_header(&block);
        assert!(matches!(result, Err(FitsError::InvalidFormat(_))));
    }

    fn value(header: &Header, name: &str) -> Option<KeywordValue> {
        header.get(name).value().cloned()
    }

    fn comment(header: &Header, name: &str) -> Option<String> {
        header.get_keyword(name).and_then(|k| k.comment.clone())
    }

    #[test]
    fn parse_header_ends_unquoted_value_at_slash() {
        let header = HeaderParser::parse_header(&header_block(&[
            "NOSPACE =                    2/ no space",
            "NOSP2   = 1.5/comment",
            "UNQSL   = 2024/01/15",
        ]))
        .unwrap();
        assert_eq!(value(&header, "NOSPACE"), Some(KeywordValue::Integer(2)));
        assert_eq!(comment(&header, "NOSPACE").as_deref(), Some("no space"));
        assert_eq!(value(&header, "NOSP2"), Some(KeywordValue::Real(1.5)));
        assert_eq!(comment(&header, "NOSP2").as_deref(), Some("comment"));
        assert_eq!(value(&header, "UNQSL"), Some(KeywordValue::Integer(2024)));
        assert_eq!(comment(&header, "UNQSL").as_deref(), Some("01/15"));
    }

    #[test]
    fn parse_header_reads_fortran_d_exponent() {
        let header = HeaderParser::parse_header(&header_block(&[
            "DEXP    =              1.0D+05 / Fortran D exponent",
            "DEXP2   =                 1D5",
            "DEXPL   =              1.0d+05",
            "EEXP    =              -2.5E-3",
        ]))
        .unwrap();
        assert_eq!(value(&header, "DEXP"), Some(KeywordValue::Real(100000.0)));
        assert_eq!(value(&header, "DEXP2"), Some(KeywordValue::Real(100000.0)));
        assert_eq!(value(&header, "DEXPL"), Some(KeywordValue::Real(100000.0)));
        assert_eq!(value(&header, "EEXP"), Some(KeywordValue::Real(-0.0025)));
    }

    #[test]
    fn parse_header_keeps_non_fits_numbers_as_text() {
        let header = HeaderParser::parse_header(&header_block(&[
            "NANV    = NAN",
            "INFV    = inf",
            "HEXV    = 0x10",
        ]))
        .unwrap();
        for (name, text) in [("NANV", "NAN"), ("INFV", "inf"), ("HEXV", "0x10")] {
            assert_eq!(
                value(&header, name),
                Some(KeywordValue::String(text.to_string()))
            );
        }
    }

    #[test]
    fn parse_header_reads_integer_beyond_i64_as_real() {
        let header =
            HeaderParser::parse_header(&header_block(&["BIGINT  = 123456789012345678901"]))
                .unwrap();
        assert_eq!(
            value(&header, "BIGINT"),
            Some(KeywordValue::Real(123456789012345678901.0))
        );
    }

    #[test]
    fn parse_header_reads_hierarch_cards() {
        let header = HeaderParser::parse_header(&header_block(&[
            "HIERARCH ESO DET DIT = 1.5 / hierarch card",
            "HIERARCH ESO OBS NAME = 'abc / def' / c",
            "HIERARCH LONGKEYWORD=3",
        ]))
        .unwrap();
        assert_eq!(
            value(&header, "HIERARCH ESO DET DIT"),
            Some(KeywordValue::Real(1.5))
        );
        assert_eq!(
            comment(&header, "HIERARCH ESO DET DIT").as_deref(),
            Some("hierarch card")
        );
        assert_eq!(
            value(&header, "HIERARCH ESO OBS NAME"),
            Some(KeywordValue::String("abc / def".to_string()))
        );
        assert_eq!(
            comment(&header, "HIERARCH ESO OBS NAME").as_deref(),
            Some("c")
        );
        assert_eq!(
            value(&header, "HIERARCH LONGKEYWORD"),
            Some(KeywordValue::Integer(3))
        );
    }

    #[test]
    fn parse_header_reads_complex_values() {
        let header = HeaderParser::parse_header(&header_block(&[
            "CPLX    =           (1.0, 2.0) / complex",
            "CPLXI   =               (1, 2)",
            "CPLXD   =      (1.0D1, -2.5E0)",
            "CPLXBAD =             (1.0, x)",
        ]))
        .unwrap();
        assert_eq!(
            value(&header, "CPLX"),
            Some(KeywordValue::Complex(1.0, 2.0))
        );
        assert_eq!(comment(&header, "CPLX").as_deref(), Some("complex"));
        assert_eq!(
            value(&header, "CPLXI"),
            Some(KeywordValue::Complex(1.0, 2.0))
        );
        assert_eq!(
            value(&header, "CPLXD"),
            Some(KeywordValue::Complex(10.0, -2.5))
        );
        assert_eq!(
            value(&header, "CPLXBAD"),
            Some(KeywordValue::String("(1.0, x)".to_string()))
        );
    }

    #[test]
    fn parse_header_reads_commentary_cards_as_text() {
        let header = HeaderParser::parse_header(&header_block(&[
            "HISTORY = some history",
            "COMMENT = some comment / not a comment field",
            "        = blank keyword",
        ]))
        .unwrap();
        let cards: Vec<(&str, Option<&KeywordValue>, Option<&str>)> = header
            .keywords()
            .iter()
            .map(|k| (k.name.as_str(), k.value.as_ref(), k.comment.as_deref()))
            .collect();
        assert_eq!(
            cards,
            [
                ("HISTORY", None, Some("= some history")),
                (
                    "COMMENT",
                    None,
                    Some("= some comment / not a comment field")
                ),
                ("", None, Some("= blank keyword")),
            ]
        );
    }

    #[test]
    fn parse_header_keeps_hierarch_prefix_apart_from_standard_names() {
        let header = HeaderParser::parse_header(&header_block(&[
            "HIERARCH NAXIS1 = 1",
            "NAXIS1  =                    4",
            "HIERARCH BZERO = 0",
        ]))
        .unwrap();
        assert_eq!(value(&header, "NAXIS1"), Some(KeywordValue::Integer(4)));
        assert_eq!(
            value(&header, "HIERARCH NAXIS1"),
            Some(KeywordValue::Integer(1))
        );
        assert_eq!(
            value(&header, "HIERARCH BZERO"),
            Some(KeywordValue::Integer(0))
        );
        assert!(!header.get("BZERO").exists());
    }

    #[test]
    fn parse_header_reads_hierarch_without_space_as_commentary() {
        let header = HeaderParser::parse_header(&header_block(&["HIERARCHESO=1"])).unwrap();
        let keyword = &header.keywords()[0];
        assert_eq!(keyword.name, "HIERARCH");
        assert_eq!(keyword.value, None);
        assert_eq!(keyword.comment.as_deref(), Some("ESO=1"));
    }

    fn header_with_cards(cards: &[&[u8]]) -> Vec<u8> {
        let mut block = Vec::new();
        for card in cards.iter().copied().chain([b"END".as_slice()]) {
            let mut card = card.to_vec();
            card.resize(CARD_SIZE, b' ');
            block.extend(card);
        }
        block.resize(FITS_BLOCK_SIZE, b' ');
        block
    }

    #[test]
    fn parse_header_replaces_non_ascii_in_free_text() {
        let header = HeaderParser::parse_header(&header_with_cards(&[
            b"OBJECT  = 'M31'              / caf\xe9",
            b"COUNT   =                    3 / tab\there",
            b"HISTORY   na\xefve text",
            b"COMMENT   nul\0",
            b"        blank \xe9 text",
            b"NOVALUE   free \xff text",
        ]))
        .unwrap();
        let cards: Vec<(&str, Option<&KeywordValue>, Option<&str>)> = header
            .keywords()
            .iter()
            .map(|k| (k.name.as_str(), k.value.as_ref(), k.comment.as_deref()))
            .collect();
        assert_eq!(
            cards,
            [
                (
                    "OBJECT",
                    Some(&KeywordValue::String("M31".to_string())),
                    Some("caf?")
                ),
                ("COUNT", Some(&KeywordValue::Integer(3)), Some("tab?here")),
                ("HISTORY", None, Some("  na?ve text")),
                ("COMMENT", None, Some("  nul?")),
                ("", None, Some("blank ? text")),
                ("NOVALUE", None, Some("  free ? text")),
            ]
        );
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
        let value = header.get("TESTKEY").value().unwrap();
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

        let simple_val = header.get("SIMPLE").value().unwrap();
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

        let swcreate = header.get("SWCREATE").value().unwrap();
        assert_eq!(
            swcreate,
            &KeywordValue::String("SharpCap 4.1.12395.0".to_string())
        );

        let simple = header.get("SIMPLE").value().unwrap();
        assert!(simple.as_logical().unwrap());

        let bitpix = header.get("BITPIX").value().unwrap();
        assert_eq!(bitpix.as_integer().unwrap(), -32);

        let naxis = header.get("NAXIS").value().unwrap();
        assert_eq!(naxis.as_integer().unwrap(), 2);
    }

    type CardParts = (String, Option<KeywordValue>, Option<String>);

    fn card_parts(header: &Header) -> Vec<CardParts> {
        header
            .keywords()
            .iter()
            .map(|k| (k.name.clone(), k.value.clone(), k.comment.clone()))
            .collect()
    }

    fn parts(name: &str, value: Option<KeywordValue>, comment: Option<&str>) -> CardParts {
        (name.to_string(), value, comment.map(str::to_string))
    }

    #[test]
    fn parse_header_replaces_non_ascii_anywhere_in_card() {
        let header = HeaderParser::parse_header(&header_with_cards(&[
            b"OBSERVER= 'Jos\xc3\xa9'",
            b"K\xe9      =                    1",
            b"HIERARCH ESO D\xe9 = 1",
            b"COUNT   = 1\xe9 / comment",
            b"OBJECT  =\xe9 1",
            b"LONG    = 'ab&'",
            b"CONTINUE  'c\xe9'",
        ]))
        .unwrap();
        assert_eq!(
            card_parts(&header),
            [
                parts("OBSERVER", Some(KeywordValue::String("Jos??".into())), None),
                parts("K?", Some(KeywordValue::Integer(1)), None),
                parts("HIERARCH ESO D?", Some(KeywordValue::Integer(1)), None),
                parts(
                    "COUNT",
                    Some(KeywordValue::String("1?".into())),
                    Some("comment")
                ),
                parts("OBJECT", None, Some("=? 1")),
                parts("LONG", Some(KeywordValue::String("abc?".into())), None),
            ]
        );
    }

    #[test]
    fn parse_header_reads_all_nul_card_as_blank() {
        let header = HeaderParser::parse_header(&header_with_cards(&[
            b"OBJECT  = 'M31'",
            &[0u8; CARD_SIZE],
            b"NUL\0    =                    1",
        ]))
        .unwrap();
        assert_eq!(
            card_parts(&header),
            [
                parts("OBJECT", Some(KeywordValue::String("M31".into())), None),
                parts("", None, None),
                parts("NUL?", Some(KeywordValue::Integer(1)), None),
            ]
        );
    }

    #[test]
    fn get_keyword_falls_back_to_hierarch_for_non_standard_names() {
        let header = HeaderParser::parse_header(&header_with_cards(&[
            b"HIERARCH ESO DET DIT = 1.5",
            b"HIERARCH NAXIS1 = 3",
            b"HIERARCH foo = 'bar'",
        ]))
        .unwrap();
        let name = |n: &str| header.get_keyword(n).map(|k| k.name.clone());
        assert_eq!(name("ESO DET DIT").as_deref(), Some("HIERARCH ESO DET DIT"));
        assert_eq!(
            name("HIERARCH ESO DET DIT").as_deref(),
            Some("HIERARCH ESO DET DIT")
        );
        assert_eq!(name("foo").as_deref(), Some("HIERARCH foo"));
        assert_eq!(name("NAXIS1"), None);
        assert_eq!(header.get("ESO DET DIT").as_f64(), Some(1.5));
        assert_eq!(header.get("foo").as_str(), Some("bar"));
        assert!(!header.get("NAXIS1").exists());
    }

    #[test]
    fn get_keyword_prefers_exact_name_over_hierarch() {
        let header = HeaderParser::parse_header(&header_with_cards(&[
            b"HIERARCH ESO DET DIT = 1.5",
            b"HIERARCH HIERARCH ESO DET DIT = 2.5",
        ]))
        .unwrap();
        assert_eq!(header.get("HIERARCH ESO DET DIT").as_f64(), Some(1.5));
    }

    #[test]
    fn parse_header_marks_cards_without_value_indicator_as_commentary() {
        let header = HeaderParser::parse_header(&header_block(&[
            "KEYWORD  free text",
            "HIERARCHESO=1",
            "OBJECT  ='M31'",
            "COMMENT   text",
            "EMPTY   =",
            "EMPTYC  =                      / only comment",
            "COUNT   =                    1",
        ]))
        .unwrap();
        let commentary: Vec<bool> = header
            .keywords()
            .iter()
            .map(|k| k.is_commentary())
            .collect();
        assert_eq!(commentary, [true, true, true, true, false, false, false]);
        assert_eq!(
            card_parts(&header)[..3],
            [
                parts("KEYWORD", None, Some(" free text")),
                parts("HIERARCH", None, Some("ESO=1")),
                parts("OBJECT", None, Some("='M31'")),
            ]
        );
    }

    #[test]
    fn get_returns_first_of_repeated_keyword() {
        let header = HeaderParser::parse_header(&header_block(&[
            "SIMPLE  = T",
            "OBJECT  = 'FIRST'",
            "OBJECT  = 'SECOND'",
            "HIERARCH ESO DET DIT = 1",
            "HIERARCH ESO DET DIT = 2",
        ]))
        .unwrap();
        assert_eq!(header.get("OBJECT").as_str(), Some("FIRST"));
        assert_eq!(header.get("ESO DET DIT").as_i64(), Some(1));
        assert_eq!(header.get("HIERARCH ESO DET DIT").as_i64(), Some(1));
        assert_eq!(header.keywords().len(), 5);
        assert_eq!(header.get_last("OBJECT").as_str(), Some("SECOND"));
        assert_eq!(header.get_last("ESO DET DIT").as_i64(), Some(2));
        assert!(!header.get_last("EXPTIME").exists());
    }

    #[test]
    fn get_keyword_with_hierarch_prefix_finds_bare_name() {
        let mut header = Header::new();
        header.add_keyword(Keyword::integer("ESO DET", 1));
        header.add_keyword(Keyword::integer("EXPTIME", 3));
        assert_eq!(header.get("HIERARCH ESO DET").as_i64(), Some(1));
        assert!(!header.get("HIERARCH EXPTIME").exists());
    }
}
