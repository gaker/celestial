use crate::fits::block::CARD_SIZE;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::keywords::{is_standard_name, Keyword, KeywordValue};
use real::format_real;

mod real;

const VALUE_START: usize = 10;
// Fixed-format numbers and logicals end in column 30 (FITS 4.0 §4.2).
const FIXED_VALUE_WIDTH: usize = 20;
// Space for a string between the quotes on one card: columns 12-79.
const MAX_STRING_ON_CARD: usize = CARD_SIZE - VALUE_START - 2;
// A continued chunk also holds the trailing '&' (FITS 4.0 §4.2.1.2).
const MAX_CONTINUED_CHUNK: usize = MAX_STRING_ON_CARD - 1;
const MIN_STRING_WIDTH: usize = 8;
const NAME_WIDTH: usize = 8;

pub(crate) fn format_keyword(keyword: &Keyword) -> Result<Vec<u8>> {
    if keyword.name == "END" {
        return Ok(card_bytes(b"END", ""));
    }
    check_text(keyword, &keyword.name)?;
    let comment = keyword.comment.as_deref();
    if let Some(text) = comment {
        check_text(keyword, text)?;
    }
    let name = keyword.name.as_str();
    if is_commentary_card(keyword) {
        commentary_card(keyword, comment)
    } else if is_standard_name(name) {
        standard_cards(keyword, comment)
    } else {
        let name = name.strip_prefix("HIERARCH ").unwrap_or(name);
        hierarch_cards(keyword, name, comment)
    }
}

// Commentary cards hold text in columns 9-80 (FITS 4.0 §4.1.2.2 and
// §4.4.2.4), so the name must fit in columns 1-8.
fn commentary_card(keyword: &Keyword, comment: Option<&str>) -> Result<Vec<u8>> {
    let name = keyword.name.as_str();
    if name.len() > NAME_WIDTH {
        return Err(FitsError::InvalidKeyword(keyword.name.clone()));
    }
    Ok(card_bytes(name.as_bytes(), comment.unwrap_or("")))
}

fn standard_cards(keyword: &Keyword, comment: Option<&str>) -> Result<Vec<u8>> {
    let name = keyword.name.as_str();
    let prefix = format!("{:<width$}= ", name, width = NAME_WIDTH);
    match &keyword.value {
        Some(KeywordValue::String(s)) => {
            string_cards(keyword, &prefix, FIXED_VALUE_WIDTH, s, comment)
        }
        Some(value) => Ok(value_card(
            name.as_bytes(),
            &value_text(keyword, value)?,
            comment,
        )),
        // An undefined value keeps the value indicator and a blank value
        // field, so the comment stays a comment (FITS 4.0 §4.1.2.3).
        None => Ok(with_comment(
            text_card(&prefix),
            VALUE_START + FIXED_VALUE_WIDTH,
            comment,
        )),
    }
}

// True when format_keyword rejects the keyword and accepts the same keyword
// under a short standard name, so the name alone makes it unwritable.
pub(crate) fn has_unwritable_name(keyword: &Keyword) -> bool {
    let mut renamed = keyword.clone();
    renamed.name = "NAME".to_string();
    format_keyword(keyword).is_err() && format_keyword(&renamed).is_ok()
}

// A value set on a commentary keyword other than HISTORY, COMMENT and the
// blank name is written as a value card, so it is kept.
pub(crate) fn is_commentary_card(keyword: &Keyword) -> bool {
    let reserved = matches!(keyword.name.as_str(), "HISTORY" | "COMMENT" | "");
    keyword.is_commentary() && (reserved || keyword.value.is_none())
}

// Any other name is written with the HIERARCH convention, as astropy writes
// it. The convention ends the name at the first '=', so a name holding one,
// or with spaces around it, cannot be written.
fn hierarch_cards(keyword: &Keyword, name: &str, comment: Option<&str>) -> Result<Vec<u8>> {
    if name.is_empty() || name.trim() != name || name.contains('=') {
        return Err(FitsError::InvalidKeyword(keyword.name.clone()));
    }
    let prefix = format!("HIERARCH {} = ", name);
    let value = match &keyword.value {
        Some(KeywordValue::String(s)) => return string_cards(keyword, &prefix, 0, s, comment),
        Some(value) => value_text(keyword, value)?,
        None => String::new(),
    };
    let text = format!("{}{}", prefix, value);
    if text.len() > CARD_SIZE {
        return Err(invalid_value(keyword, value));
    }
    Ok(with_comment(text_card(&text), text.len(), comment))
}

// Header text must be printable ASCII, 0x20 to 0x7E.
fn check_text(keyword: &Keyword, text: &str) -> Result<()> {
    if text.bytes().all(|b| (b' '..=b'~').contains(&b)) {
        return Ok(());
    }
    Err(invalid_value(keyword, format!("{:?}", text)))
}

fn invalid_value(keyword: &Keyword, value: String) -> FitsError {
    FitsError::InvalidKeywordValue {
        keyword: keyword.name.clone(),
        value,
    }
}

// The value field as it is written on a single card. A string is quoted
// with no padding; `string_cards` lays out strings with padding.
pub(crate) fn value_text(keyword: &Keyword, value: &KeywordValue) -> Result<String> {
    if let KeywordValue::String(s) = value {
        check_text(keyword, s)?;
    }
    value_field(value, |f| Err(invalid_value(keyword, f.to_string())))
}

// `value_text` without the checks. `non_finite` formats a NaN or infinite
// real, which FITS cannot represent.
pub(crate) fn value_field<E>(
    value: &KeywordValue,
    non_finite: impl Fn(f64) -> std::result::Result<String, E>,
) -> std::result::Result<String, E> {
    let real = |f: f64| format_real(f).map_or_else(|| non_finite(f), Ok);
    Ok(match value {
        KeywordValue::Logical(b) => (if *b { "T" } else { "F" }).to_string(),
        KeywordValue::Integer(i) => i.to_string(),
        KeywordValue::Real(f) => real(*f)?,
        KeywordValue::Complex(re, im) => format!("({}, {})", real(*re)?, real(*im)?),
        KeywordValue::String(s) => format!("'{}'", s.trim_end_matches(' ').replace('\'', "''")),
    })
}

fn card_bytes(name: &[u8], rest: &str) -> Vec<u8> {
    let mut card = vec![b' '; CARD_SIZE];
    card[..name.len()].copy_from_slice(name);
    let rest = &rest.as_bytes()[..rest.len().min(CARD_SIZE - 8)];
    card[8..8 + rest.len()].copy_from_slice(rest);
    card
}

// Callers pass at most CARD_SIZE bytes; anything longer is cut.
fn text_card(text: &str) -> Vec<u8> {
    let mut card = vec![b' '; CARD_SIZE];
    let text = &text.as_bytes()[..text.len().min(CARD_SIZE)];
    card[..text.len()].copy_from_slice(text);
    card
}

fn value_card(name: &[u8], value: &str, comment: Option<&str>) -> Vec<u8> {
    let field = format!("= {:>width$}", value, width = FIXED_VALUE_WIDTH);
    with_comment(
        card_bytes(name, &field),
        VALUE_START + value.len().max(FIXED_VALUE_WIDTH),
        comment,
    )
}

// Comments start at column 32 after a fixed-format value, and right after
// a longer value otherwise. They are cut at column 80.
fn with_comment(mut card: Vec<u8>, value_end: usize, comment: Option<&str>) -> Vec<u8> {
    let Some(comment) = comment else {
        return card;
    };
    let text = format!(" / {}", comment);
    let room = CARD_SIZE.saturating_sub(value_end);
    let text = &text.as_bytes()[..text.len().min(room)];
    card[value_end..value_end + text.len()].copy_from_slice(text);
    card
}

// Trailing spaces in a string are not significant (FITS 4.0 §4.2.1.1).
// Dropping them first keeps a CONTINUE card from holding only spaces.
// `prefix` is the text before the opening quote, and a comment on a single
// card starts no earlier than `min_field` columns after it.
fn string_cards(
    keyword: &Keyword,
    prefix: &str,
    min_field: usize,
    value: &str,
    comment: Option<&str>,
) -> Result<Vec<u8>> {
    check_text(keyword, value)?;
    let value = value.trim_end_matches(' ');
    let escaped = value.replace('\'', "''");
    let quoted = format!("'{:<width$}'", escaped, width = MIN_STRING_WIDTH);
    if prefix.len() + quoted.len() <= CARD_SIZE {
        let card = text_card(&format!("{}{}", prefix, quoted));
        let value_end = prefix.len() + quoted.len().max(min_field);
        return Ok(with_comment(card, value_end, comment));
    }
    continued_string_cards(keyword, prefix, value, comment)
}

// Long-string convention (FITS 4.0 §4.2.1.2). The first chunk fills the
// card after `prefix`; a '' pair needs room for two characters there. A
// comment that does not fit after the last chunk goes on a final
// "CONTINUE  ''" card, as astropy writes it.
fn continued_string_cards(
    keyword: &Keyword,
    prefix: &str,
    value: &str,
    comment: Option<&str>,
) -> Result<Vec<u8>> {
    let first_capacity = CARD_SIZE.saturating_sub(prefix.len() + 3);
    if first_capacity < 2 {
        return Err(invalid_value(keyword, value.to_string()));
    }
    let mut chunks = continued_chunks(value, first_capacity);
    let lead = |i: usize| if i == 0 { prefix } else { "CONTINUE  " };
    if let (Some(text), Some(last)) = (comment, chunks.last()) {
        let last_card_len = lead(chunks.len() - 1).len() + last.len() + 2;
        if last_card_len + " / ".len() + text.len() > CARD_SIZE {
            chunks.push(String::new());
        }
    }
    let mut cards = Vec::with_capacity(chunks.len() * CARD_SIZE);
    for (i, chunk) in chunks.iter().enumerate() {
        let is_last = i + 1 == chunks.len();
        let quoted = format!("'{}{}'", chunk, if is_last { "" } else { "&" });
        let card = text_card(&format!("{}{}", lead(i), quoted));
        let comment = if is_last { comment } else { None };
        cards.extend(with_comment(card, lead(i).len() + quoted.len(), comment));
    }
    Ok(cards)
}

// Splits the escaped value into chunks for the long-string convention.
// A '' pair stays in one chunk, so each chunk reads back on its own.
fn continued_chunks(value: &str, first_capacity: usize) -> Vec<String> {
    let mut chunks = Vec::new();
    let mut current = String::new();
    let mut capacity = first_capacity;
    for c in value.chars() {
        let unit_len = if c == '\'' { 2 } else { 1 };
        if current.len() + unit_len > capacity {
            chunks.push(std::mem::take(&mut current));
            capacity = MAX_CONTINUED_CHUNK;
        }
        if c == '\'' {
            current.push_str("''");
        } else {
            current.push(c);
        }
    }
    chunks.push(current);
    chunks
}
