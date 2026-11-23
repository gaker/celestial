use super::HeaderCard;
use crate::fits::errors::Result;
use crate::fits::header::keywords::{Keyword, KeywordValue};

impl HeaderCard {
    pub(crate) fn to_keyword(&self) -> Result<Keyword> {
        let mut keyword = Keyword::new(self.keyword.clone());
        keyword.commentary = self.commentary;

        if let Some(comment) = &self.comment {
            keyword = keyword.with_comment(comment.clone());
        }

        if let Some(value_str) = &self.value {
            let parsed_value = Self::parse_value(value_str)?;
            keyword = keyword.with_value(parsed_value);
        }

        Ok(keyword)
    }

    pub(crate) fn parse_value(value_str: &str) -> Result<KeywordValue> {
        let trimmed = value_str.trim();

        if trimmed == "T" {
            return Ok(KeywordValue::Logical(true));
        }
        if trimmed == "F" {
            return Ok(KeywordValue::Logical(false));
        }

        if trimmed.starts_with('\'') && trimmed.ends_with('\'') && trimmed.len() >= 2 {
            let string_content = &trimmed[1..trimmed.len() - 1];
            let unescaped = string_content.trim_end().replace("''", "'");
            return Ok(KeywordValue::String(unescaped));
        }

        if let Some(number) = Self::parse_number(trimmed) {
            return Ok(number);
        }

        if let Some((real, imaginary)) = Self::parse_complex(trimmed) {
            return Ok(KeywordValue::Complex(real, imaginary));
        }

        Ok(KeywordValue::String(trimmed.to_string()))
    }

    // An integer that does not fit i64 is kept as the nearest real, so a
    // mandatory keyword that needs an integer rejects it.
    fn parse_number(text: &str) -> Option<KeywordValue> {
        let digits = text.strip_prefix(['+', '-']).unwrap_or(text);
        if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
            return match text.parse::<i64>() {
                Ok(integer) => Some(KeywordValue::Integer(integer)),
                Err(_) => Self::parse_real(text).map(KeywordValue::Real),
            };
        }
        Self::parse_real(text).map(KeywordValue::Real)
    }

    // FITS 4.0 §4.2.4: [sign] digits [. digits] [exponent], where the
    // exponent letter is E or D. Lower-case e and d are read as astropy
    // reads them. NaN, infinities and overflow are not FITS numbers.
    fn parse_real(text: &str) -> Option<f64> {
        let all_digits = |s: &str| s.bytes().all(|b| b.is_ascii_digit());
        let unsigned = text.strip_prefix(['+', '-']).unwrap_or(text);
        let (mantissa, exponent) = match unsigned.split_once(['E', 'e', 'D', 'd']) {
            Some((mantissa, exponent)) => (mantissa, Some(exponent)),
            None => (unsigned, None),
        };
        let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
        if whole.len() + fraction.len() == 0 || !all_digits(whole) || !all_digits(fraction) {
            return None;
        }
        if let Some(exponent) = exponent {
            let exponent = exponent.strip_prefix(['+', '-']).unwrap_or(exponent);
            if exponent.is_empty() || !all_digits(exponent) {
                return None;
            }
        }
        let value: f64 = text.replace(['D', 'd'], "E").parse().ok()?;
        value.is_finite().then_some(value)
    }

    // FITS 4.0 §4.2.5 and §4.2.6: "(real, imaginary)", where each part is
    // an integer or a real.
    fn parse_complex(text: &str) -> Option<(f64, f64)> {
        let inner = text.strip_prefix('(')?.strip_suffix(')')?;
        let (real, imaginary) = inner.split_once(',')?;
        let part = |s: &str| match Self::parse_number(s.trim())? {
            KeywordValue::Integer(i) => Some(i as f64),
            KeywordValue::Real(r) => Some(r),
            _ => None,
        };
        Some((part(real)?, part(imaginary)?))
    }
}
