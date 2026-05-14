#[derive(Debug, Clone)]
pub enum ColorSpace {
    Gray,
    Rgb,
    Unknown(String),
}

impl ColorSpace {
    pub fn parse(s: &str) -> Self {
        match s {
            "Gray" => Self::Gray,
            "RGB" => Self::Rgb,
            other => Self::Unknown(other.to_string()),
        }
    }

    pub fn as_str(&self) -> &str {
        match self {
            Self::Gray => "Gray",
            Self::Rgb => "RGB",
            Self::Unknown(s) => s,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_known_variants() {
        assert!(matches!(ColorSpace::parse("Gray"), ColorSpace::Gray));
        assert!(matches!(ColorSpace::parse("RGB"), ColorSpace::Rgb));
    }

    #[test]
    fn parse_is_case_sensitive() {
        assert!(matches!(ColorSpace::parse("gray"), ColorSpace::Unknown(_)));
        assert!(matches!(ColorSpace::parse("rgb"), ColorSpace::Unknown(_)));
        assert!(matches!(ColorSpace::parse("GRAY"), ColorSpace::Unknown(_)));
    }

    #[test]
    fn parse_unknown_preserves_input() {
        match ColorSpace::parse("CIELab") {
            ColorSpace::Unknown(s) => assert_eq!(s, "CIELab"),
            _ => panic!("expected Unknown variant"),
        }
    }

    #[test]
    fn parse_empty_is_unknown() {
        match ColorSpace::parse("") {
            ColorSpace::Unknown(s) => assert_eq!(s, ""),
            _ => panic!("expected Unknown variant"),
        }
    }

    #[test]
    fn as_str_roundtrip_known() {
        assert_eq!(ColorSpace::parse(ColorSpace::Gray.as_str()).as_str(), "Gray");
        assert_eq!(ColorSpace::parse(ColorSpace::Rgb.as_str()).as_str(), "RGB");
    }

    #[test]
    fn as_str_on_unknown_returns_stored_value() {
        let cs = ColorSpace::Unknown("HSV".to_string());
        assert_eq!(cs.as_str(), "HSV");
    }

    #[test]
    fn parse_then_as_str_roundtrip_unknown() {
        let original = "CMYK";
        let parsed = ColorSpace::parse(original);
        assert_eq!(parsed.as_str(), original);
    }
}
