use super::{HeaderCard, CARD_SIZE};
use crate::fits::errors::{FitsError, Result};
use std::str;

impl HeaderCard {
    #[cfg(test)]
    pub(crate) fn parse(data: &[u8; CARD_SIZE]) -> Result<Self> {
        Self::parse_printable(&Self::printable_card(data))
    }

    pub(super) fn parse_printable(text: &[u8; CARD_SIZE]) -> Result<Self> {
        let card_str = Self::card_str(text)?;
        let mut card = Self::create_empty_card();

        card.keyword = Self::extract_keyword(card_str)?;

        Self::parse_value_and_comment(card_str, &mut card);

        Ok(card)
    }

    pub(super) fn card_str(text: &[u8; CARD_SIZE]) -> Result<&str> {
        str::from_utf8(text)
            .map_err(|_| FitsError::InvalidFormat("Invalid UTF-8 in header card".to_string()))
    }

    // FITS 4.0 §4.1.1 allows only printable ASCII in a card. astropy reads
    // any other byte as '?', so a value such as "José" does not reject the
    // file. A card of NUL bytes reads as a blank card. The result is
    // printable ASCII, so slicing it by index is safe.
    pub(super) fn printable_card(data: &[u8; CARD_SIZE]) -> [u8; CARD_SIZE] {
        if data.iter().all(|&b| b == 0) {
            return [b' '; CARD_SIZE];
        }
        data.map(|b| if (b' '..=b'~').contains(&b) { b } else { b'?' })
    }

    fn create_empty_card() -> Self {
        Self {
            keyword: String::new(),
            value: None,
            comment: None,
            commentary: false,
        }
    }

    fn extract_keyword(card_str: &str) -> Result<String> {
        let keyword_part = &card_str[0..8];
        Ok(keyword_part.trim().to_string())
    }

    fn parse_value_and_comment(card_str: &str, card: &mut Self) {
        let rest = &card_str[8..];
        if Self::is_commentary(&card.keyword) {
            Self::parse_comment_only(rest, card);
        } else if card.keyword == "HIERARCH" && rest.starts_with(' ') {
            Self::parse_hierarch(rest, card);
        } else if let Some(value_comment_part) = rest.strip_prefix("= ") {
            Self::parse_keyword_value_pair(value_comment_part, card);
        } else {
            Self::parse_comment_only(rest, card);
        }
    }

    // COMMENT, HISTORY and blank keywords have no value even when columns
    // 9-10 hold "= " (FITS 4.0 §4.4.2.4).
    fn is_commentary(keyword: &str) -> bool {
        matches!(keyword, "COMMENT" | "HISTORY" | "")
    }

    // HIERARCH convention: the name runs from column 10 to the first '=',
    // and the value field follows it. The name keeps the HIERARCH prefix,
    // so "HIERARCH NAXIS1" never replaces the mandatory NAXIS1.
    fn parse_hierarch(rest: &str, card: &mut Self) {
        match rest.split_once('=') {
            Some((name, value_comment_part)) if !name.trim().is_empty() => {
                card.keyword = format!("HIERARCH {}", name.trim());
                Self::parse_keyword_value_pair(value_comment_part, card);
            }
            _ => Self::parse_comment_only(rest, card),
        }
    }

    fn parse_keyword_value_pair(value_comment_part: &str, card: &mut Self) {
        if let Some((quoted, rest)) = Self::split_quoted(value_comment_part) {
            card.value = Some(quoted.to_string());
            let comment = rest.trim_start().strip_prefix('/').map(str::trim);
            card.comment = comment.filter(|c| !c.is_empty()).map(str::to_string);
            return;
        }
        // An unquoted value cannot hold '/', so the first one starts the
        // comment, with or without spaces around it (FITS 4.0 §4.1.2.3).
        let (value_part, comment_part) = match value_comment_part.split_once('/') {
            Some((value, comment)) => (value.trim(), comment.trim()),
            None => (value_comment_part.trim(), ""),
        };
        if !value_part.is_empty() {
            card.value = Some(value_part.to_string());
        }
        if !comment_part.is_empty() {
            card.comment = Some(comment_part.to_string());
        }
    }

    // Returns the quoted string token and the text after it. A doubled quote
    // is part of the string (FITS 4.0 §4.2.1.1), so it does not end the token.
    pub(super) fn split_quoted(text: &str) -> Option<(&str, &str)> {
        let text = text.trim_start();
        let bytes = text.as_bytes();
        if bytes.first() != Some(&b'\'') {
            return None;
        }
        let mut i = 1;
        while i < bytes.len() {
            if bytes[i] == b'\'' {
                if bytes.get(i + 1) == Some(&b'\'') {
                    i += 2;
                    continue;
                }
                return Some(text.split_at(i + 1));
            }
            i += 1;
        }
        None
    }

    pub(super) fn is_quoted_value(value: &str) -> bool {
        Self::split_quoted(value).is_some_and(|(_, rest)| rest.is_empty())
    }

    // Commentary text keeps its leading spaces, as astropy reads it.
    pub(super) fn parse_comment_only(rest_of_card: &str, card: &mut Self) {
        card.commentary = true;
        let comment_part = rest_of_card.trim_end();
        if !comment_part.is_empty() {
            card.comment = Some(comment_part.to_string());
        }
    }
}
