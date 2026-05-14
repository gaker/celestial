#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum XisfCompression {
    #[default]
    None,
    Lz4,
    Lz4Hc,
    Zlib,
    Zstd,
}

impl XisfCompression {
    pub fn parse(s: &str) -> Self {
        let lower = s.to_lowercase();
        if lower.starts_with("lz4+hc") || lower.starts_with("lz4-hc") {
            Self::Lz4Hc
        } else if lower.starts_with("lz4") {
            Self::Lz4
        } else if lower.starts_with("zlib") {
            Self::Zlib
        } else if lower.starts_with("zstd") || lower.starts_with("zstandard") {
            Self::Zstd
        } else {
            Self::None
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::None => "",
            Self::Lz4 => "lz4",
            Self::Lz4Hc => "lz4-hc",
            Self::Zlib => "zlib",
            Self::Zstd => "zstd",
        }
    }

    pub fn is_compressed(&self) -> bool {
        !matches!(self, Self::None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_recognized_lowercase() {
        assert_eq!(XisfCompression::parse("lz4"), XisfCompression::Lz4);
        assert_eq!(XisfCompression::parse("lz4-hc"), XisfCompression::Lz4Hc);
        assert_eq!(XisfCompression::parse("lz4+hc"), XisfCompression::Lz4Hc);
        assert_eq!(XisfCompression::parse("zlib"), XisfCompression::Zlib);
        assert_eq!(XisfCompression::parse("zstd"), XisfCompression::Zstd);
        assert_eq!(XisfCompression::parse("zstandard"), XisfCompression::Zstd);
    }

    #[test]
    fn parse_is_case_insensitive() {
        assert_eq!(XisfCompression::parse("LZ4"), XisfCompression::Lz4);
        assert_eq!(XisfCompression::parse("Lz4-Hc"), XisfCompression::Lz4Hc);
        assert_eq!(XisfCompression::parse("ZLIB"), XisfCompression::Zlib);
        assert_eq!(XisfCompression::parse("ZStd"), XisfCompression::Zstd);
    }

    #[test]
    fn parse_accepts_sublevel_suffixes() {
        assert_eq!(XisfCompression::parse("lz4:5"), XisfCompression::Lz4);
        assert_eq!(XisfCompression::parse("zlib:9"), XisfCompression::Zlib);
        assert_eq!(XisfCompression::parse("zstd:22"), XisfCompression::Zstd);
        assert_eq!(XisfCompression::parse("lz4-hc:12"), XisfCompression::Lz4Hc);
    }

    #[test]
    fn lz4_hc_beats_lz4_prefix() {
        assert_eq!(XisfCompression::parse("lz4-hc"), XisfCompression::Lz4Hc);
        assert_eq!(XisfCompression::parse("lz4+hc"), XisfCompression::Lz4Hc);
    }

    #[test]
    fn parse_unknown_returns_none() {
        assert_eq!(XisfCompression::parse(""), XisfCompression::None);
        assert_eq!(XisfCompression::parse("none"), XisfCompression::None);
        assert_eq!(XisfCompression::parse("gzip"), XisfCompression::None);
        assert_eq!(XisfCompression::parse("rice"), XisfCompression::None);
    }

    #[test]
    fn as_str_maps_all_variants() {
        assert_eq!(XisfCompression::None.as_str(), "");
        assert_eq!(XisfCompression::Lz4.as_str(), "lz4");
        assert_eq!(XisfCompression::Lz4Hc.as_str(), "lz4-hc");
        assert_eq!(XisfCompression::Zlib.as_str(), "zlib");
        assert_eq!(XisfCompression::Zstd.as_str(), "zstd");
    }

    #[test]
    fn is_compressed_is_true_for_every_variant_except_none() {
        assert!(!XisfCompression::None.is_compressed());
        assert!(XisfCompression::Lz4.is_compressed());
        assert!(XisfCompression::Lz4Hc.is_compressed());
        assert!(XisfCompression::Zlib.is_compressed());
        assert!(XisfCompression::Zstd.is_compressed());
    }

    #[test]
    fn parse_as_str_roundtrip_for_non_none() {
        for variant in [
            XisfCompression::Lz4,
            XisfCompression::Lz4Hc,
            XisfCompression::Zlib,
            XisfCompression::Zstd,
        ] {
            assert_eq!(XisfCompression::parse(variant.as_str()), variant);
        }
    }

    #[test]
    fn default_is_none() {
        assert_eq!(XisfCompression::default(), XisfCompression::None);
    }
}
