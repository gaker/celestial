use crate::core::BitPix;
use crate::xisf::{Result, XisfError};

#[derive(Debug, Clone)]
pub enum SampleFormat {
    UInt8,
    UInt16,
    UInt32,
    Float32,
    Float64,
}

impl SampleFormat {
    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "UInt8" => Ok(Self::UInt8),
            "UInt16" => Ok(Self::UInt16),
            "UInt32" => Ok(Self::UInt32),
            "Float32" => Ok(Self::Float32),
            "Float64" => Ok(Self::Float64),
            _ => Err(XisfError::UnsupportedFormat(s.to_string())),
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::UInt8 => "UInt8",
            Self::UInt16 => "UInt16",
            Self::UInt32 => "UInt32",
            Self::Float32 => "Float32",
            Self::Float64 => "Float64",
        }
    }

    pub fn to_bitpix(&self) -> BitPix {
        match self {
            Self::UInt8 => BitPix::U8,
            Self::UInt16 => BitPix::I16,
            Self::UInt32 => BitPix::I32,
            Self::Float32 => BitPix::F32,
            Self::Float64 => BitPix::F64,
        }
    }

    pub fn bytes_per_sample(&self) -> usize {
        match self {
            Self::UInt8 => 1,
            Self::UInt16 => 2,
            Self::UInt32 => 4,
            Self::Float32 => 4,
            Self::Float64 => 8,
        }
    }

    pub fn is_floating_point(&self) -> bool {
        matches!(self, Self::Float32 | Self::Float64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_variants() -> [SampleFormat; 5] {
        [
            SampleFormat::UInt8,
            SampleFormat::UInt16,
            SampleFormat::UInt32,
            SampleFormat::Float32,
            SampleFormat::Float64,
        ]
    }

    #[test]
    fn parse_accepts_every_canonical_name() {
        let cases = [
            ("UInt8", BitPix::U8, 1, false),
            ("UInt16", BitPix::I16, 2, false),
            ("UInt32", BitPix::I32, 4, false),
            ("Float32", BitPix::F32, 4, true),
            ("Float64", BitPix::F64, 8, true),
        ];
        for (name, bitpix, bps, is_fp) in cases {
            let parsed = SampleFormat::parse(name).expect("valid name parses");
            assert_eq!(parsed.as_str(), name);
            assert_eq!(parsed.to_bitpix(), bitpix);
            assert_eq!(parsed.bytes_per_sample(), bps);
            assert_eq!(parsed.is_floating_point(), is_fp);
        }
    }

    #[test]
    fn parse_is_case_sensitive() {
        for bad in ["uint8", "UINT8", "float32", "FLOAT64"] {
            assert!(SampleFormat::parse(bad).is_err(), "'{bad}' should fail");
        }
    }

    #[test]
    fn parse_rejects_unsupported_types() {
        for bad in ["", "Int8", "Int16", "Complex64", "Float16", "Bool", "UInt64"] {
            let err = SampleFormat::parse(bad).expect_err("should fail");
            assert!(
                err.to_string().contains(bad) || bad.is_empty(),
                "error should mention input: {err}"
            );
        }
    }

    #[test]
    fn as_str_then_parse_roundtrip() {
        for variant in all_variants() {
            let s = variant.as_str();
            let reparsed = SampleFormat::parse(s).expect("roundtrip");
            assert_eq!(reparsed.as_str(), s);
        }
    }

    #[test]
    fn integer_variants_are_not_floating_point() {
        assert!(!SampleFormat::UInt8.is_floating_point());
        assert!(!SampleFormat::UInt16.is_floating_point());
        assert!(!SampleFormat::UInt32.is_floating_point());
    }
}
