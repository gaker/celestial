#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PixelStorage {
    #[default]
    Planar,
    Normal,
}

impl PixelStorage {
    pub fn parse(s: &str) -> Self {
        match s {
            "Normal" => Self::Normal,
            _ => Self::Planar,
        }
    }

    pub fn as_str(&self) -> &str {
        match self {
            Self::Planar => "Planar",
            Self::Normal => "Normal",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_normal_returns_normal() {
        assert_eq!(PixelStorage::parse("Normal"), PixelStorage::Normal);
    }

    #[test]
    fn parse_planar_returns_planar() {
        assert_eq!(PixelStorage::parse("Planar"), PixelStorage::Planar);
    }

    #[test]
    fn parse_anything_else_defaults_to_planar() {
        assert_eq!(PixelStorage::parse(""), PixelStorage::Planar);
        assert_eq!(PixelStorage::parse("normal"), PixelStorage::Planar);
        assert_eq!(PixelStorage::parse("NORMAL"), PixelStorage::Planar);
        assert_eq!(PixelStorage::parse("Interleaved"), PixelStorage::Planar);
    }

    #[test]
    fn as_str_matches_canonical_names() {
        assert_eq!(PixelStorage::Planar.as_str(), "Planar");
        assert_eq!(PixelStorage::Normal.as_str(), "Normal");
    }

    #[test]
    fn default_is_planar() {
        assert_eq!(PixelStorage::default(), PixelStorage::Planar);
    }

    #[test]
    fn as_str_then_parse_roundtrip() {
        for variant in [PixelStorage::Planar, PixelStorage::Normal] {
            assert_eq!(PixelStorage::parse(variant.as_str()), variant);
        }
    }
}
