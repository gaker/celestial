use crate::fits::header::Keyword;
use celestial_wcs::{WcsKeyword, WcsKeywordValue};

impl From<WcsKeyword> for Keyword {
    fn from(wk: WcsKeyword) -> Self {
        match wk.value {
            WcsKeywordValue::Real(v) => Self::real(wk.name, v),
            WcsKeywordValue::Integer(v) => Self::integer(wk.name, v),
            WcsKeywordValue::String(v) => Self::string(wk.name, v),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::KeywordValue;

    #[test]
    fn real_conversion_preserves_name_and_value() {
        let wk = WcsKeyword::real("CRVAL1", 180.5);
        let k = Keyword::from(wk);
        assert_eq!(k.name, "CRVAL1");
        assert_eq!(k.value, Some(KeywordValue::Real(180.5)));
    }

    #[test]
    fn integer_conversion_preserves_name_and_value() {
        let wk = WcsKeyword::integer("NAXIS", 2);
        let k = Keyword::from(wk);
        assert_eq!(k.name, "NAXIS");
        assert_eq!(k.value, Some(KeywordValue::Integer(2)));
    }

    #[test]
    fn string_conversion_preserves_name_and_value() {
        let wk = WcsKeyword {
            name: "CTYPE1".to_string(),
            value: WcsKeywordValue::String("RA---TAN".to_string()),
        };
        let k = Keyword::from(wk);
        assert_eq!(k.name, "CTYPE1");
        assert_eq!(k.value, Some(KeywordValue::String("RA---TAN".to_string())));
    }

    #[test]
    fn conversion_drops_original() {
        let wk = WcsKeyword::real("CD1_1", 0.0001);
        let _k: Keyword = wk.into();
    }
}
