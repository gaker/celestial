use crate::fits::errors::{FitsError, Result};
use crate::fits::header::keywords::KeywordValue;
use crate::fits::header::parser::Header;

// FITS 4.0 Tables 11 and 19: with a scale of 1, these offsets store i8,
// u16, u32 and u64 in the stored integer type of the same width.
pub(crate) const I8_OFFSET: i128 = -128;
pub(crate) const U16_OFFSET: i128 = 1 << 15;
pub(crate) const U32_OFFSET: i128 = 1 << 31;
pub(crate) const U64_OFFSET: i128 = 1 << 63;

// FITS 4.0 eq. 3 (BSCALE, BZERO) and eq. 7 (TSCALn, TZEROn):
// physical = zero + scale * stored.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Linear {
    scale: f64,
    zero: f64,
    // An integer zero card is compared as an integer, because i64::MAX and
    // 2^63 are the same f64.
    zero_integer: Option<i64>,
}

impl Linear {
    // A missing card takes its default. A card that holds no finite number
    // is an error.
    pub(crate) fn from_header(header: &Header, scale: &str, zero: &str) -> Result<Self> {
        Ok(Self {
            scale: finite_real(header, scale)?.unwrap_or(1.0),
            zero: finite_real(header, zero)?.unwrap_or(0.0),
            zero_integer: header.get_last(zero).as_i64(),
        })
    }

    pub(crate) fn is_identity(&self) -> bool {
        self.scale == 1.0 && self.zero == 0.0
    }

    pub(crate) fn is_offset(&self, zero: i128) -> bool {
        let exact = match self.zero_integer {
            Some(card) => i128::from(card) == zero,
            None => self.zero == zero as f64,
        };
        self.scale == 1.0 && exact
    }

    pub(crate) fn apply(&self, stored: f64) -> f64 {
        self.zero + self.scale * stored
    }

    pub(crate) fn scale(&self) -> f64 {
        self.scale
    }
}

fn finite_real(header: &Header, name: &str) -> Result<Option<f64>> {
    let field = header.get_last(name);
    let Some(value) = field.value() else {
        return Ok(None);
    };
    match field.as_f64() {
        Some(v) if v.is_finite() => Ok(Some(v)),
        _ => Err(FitsError::InvalidKeywordValue {
            keyword: name.to_string(),
            value: value.to_string(),
        }),
    }
}

// Every f64 in this half-open range converts to i64 without saturating.
const I64_REALS: std::ops::Range<f64> = -9223372036854775808.0..9223372036854775808.0;

// An integer card, or a real card that holds an exact i64 such as -32768.0.
pub(crate) fn integer_value(value: &KeywordValue) -> Option<i64> {
    match *value {
        KeywordValue::Integer(i) => Some(i),
        KeywordValue::Real(r) if libm::trunc(r) == r && I64_REALS.contains(&r) => Some(r as i64),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::Keyword;

    fn header(cards: &[Keyword]) -> Header {
        let mut h = Header::new();
        for card in cards {
            h.add_keyword(card.clone());
        }
        h
    }

    #[test]
    fn missing_cards_give_the_identity() {
        let linear = Linear::from_header(&Header::new(), "TSCAL1", "TZERO1").unwrap();
        assert!(linear.is_identity());
        assert!(!linear.is_offset(U16_OFFSET));
        assert_eq!(linear.apply(-3.0), -3.0);
    }

    #[test]
    fn offsets_compare_integer_cards_exactly() {
        let h = header(&[Keyword::integer("TZERO1", i64::MAX)]);
        let linear = Linear::from_header(&h, "TSCAL1", "TZERO1").unwrap();
        assert!(!linear.is_offset(U64_OFFSET));
        let h = header(&[Keyword::real("TZERO1", 9223372036854775808.0)]);
        let linear = Linear::from_header(&h, "TSCAL1", "TZERO1").unwrap();
        assert!(linear.is_offset(U64_OFFSET));
        let h = header(&[
            Keyword::integer("TZERO1", 32768),
            Keyword::real("TSCAL1", 2.0),
        ]);
        let linear = Linear::from_header(&h, "TSCAL1", "TZERO1").unwrap();
        assert!(!linear.is_offset(U16_OFFSET));
        assert_eq!(linear.apply(3.0), 32774.0);
    }

    #[test]
    fn non_numeric_scale_is_an_error() {
        let h = header(&[Keyword::string("TSCAL1", "x")]);
        assert!(Linear::from_header(&h, "TSCAL1", "TZERO1").is_err());
    }

    #[test]
    fn integer_value_accepts_exact_reals_only() {
        assert_eq!(integer_value(&KeywordValue::Integer(-3)), Some(-3));
        assert_eq!(integer_value(&KeywordValue::Real(-32768.0)), Some(-32768));
        assert_eq!(integer_value(&KeywordValue::Real(1.5)), None);
        assert_eq!(integer_value(&KeywordValue::Real(9.3e18)), None);
        assert_eq!(integer_value(&KeywordValue::String("1".into())), None);
    }
}
