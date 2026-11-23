use crate::xisf::errors::{Result, XisfError};

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
    // Codec names are the exact lowercase names of XISF 1.0 section 10.6.
    // The flag is true for the byte-shuffled ("+sh") form of the codec.
    pub(crate) fn parse(name: &str) -> Result<(Self, bool)> {
        let (base, shuffled) = match name.strip_suffix("+sh") {
            Some(base) => (base, true),
            None => (name, false),
        };
        let codec = match base {
            "zlib" => Self::Zlib,
            "lz4" => Self::Lz4,
            "lz4hc" => Self::Lz4Hc,
            "zstd" => Self::Zstd,
            _ => return Err(XisfError::UnsupportedCodec(name.to_string())),
        };
        Ok((codec, shuffled))
    }

    pub(crate) fn as_str(&self) -> &'static str {
        match self {
            Self::None => "",
            Self::Lz4 => "lz4",
            Self::Lz4Hc => "lz4hc",
            Self::Zlib => "zlib",
            Self::Zstd => "zstd",
        }
    }

    pub(crate) fn is_compressed(&self) -> bool {
        !matches!(self, Self::None)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CompressionAttr {
    pub(crate) codec: XisfCompression,
    pub(crate) uncompressed_size: u64,
    pub(crate) shuffle_item_size: Option<usize>,
}

impl CompressionAttr {
    // Syntax "codec:uncompressed-size", plus ":item-size" for shuffled codecs.
    pub(crate) fn parse(value: &str) -> Result<Self> {
        let parts: Vec<&str> = value.split(':').collect();
        let (codec, shuffled) = XisfCompression::parse(parts.first().copied().unwrap_or(""))?;
        let field_count = if shuffled { 3 } else { 2 };
        if parts.len() != field_count {
            return Err(invalid_attr("compression", value));
        }
        let uncompressed_size = parse_size(parts[1], "compression", value)?;
        let shuffle_item_size = match parts.get(2) {
            Some(item) => Some(parse_item_size(item, value)?),
            None => None,
        };
        Ok(Self {
            codec,
            uncompressed_size,
            shuffle_item_size,
        })
    }
}

// Syntax "c1,u1:c2,u2:...": compressed and uncompressed lengths of each
// subblock (XISF 1.0 section 10.6).
pub(crate) fn parse_subblocks(value: &str) -> Result<Vec<(u64, u64)>> {
    value
        .split(':')
        .map(|pair| {
            let (compressed, uncompressed) = pair
                .split_once(',')
                .ok_or_else(|| invalid_attr("subblocks", value))?;
            Ok((
                parse_size(compressed, "subblocks", value)?,
                parse_size(uncompressed, "subblocks", value)?,
            ))
        })
        .collect()
}

// Rust's integer parsing accepts a leading '+', which the spec's plain
// decimal sizes do not allow.
fn parse_size(field: &str, attr: &str, value: &str) -> Result<u64> {
    if !field.bytes().all(|b| b.is_ascii_digit()) {
        return Err(invalid_attr(attr, value));
    }
    field.parse().map_err(|_| invalid_attr(attr, value))
}

// The shuffle transform is undefined for an item size of zero.
fn parse_item_size(field: &str, value: &str) -> Result<usize> {
    let size = parse_size(field, "compression", value)?;
    match usize::try_from(size) {
        Ok(size) if size > 0 => Ok(size),
        _ => Err(invalid_attr("compression", value)),
    }
}

fn invalid_attr(attr: &str, value: &str) -> XisfError {
    XisfError::InvalidFormat(format!("Invalid {attr} attribute: \"{value}\""))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_recognizes_spec_codec_names() {
        assert_eq!(
            XisfCompression::parse("zlib").unwrap(),
            (XisfCompression::Zlib, false)
        );
        assert_eq!(
            XisfCompression::parse("lz4").unwrap(),
            (XisfCompression::Lz4, false)
        );
        assert_eq!(
            XisfCompression::parse("lz4hc").unwrap(),
            (XisfCompression::Lz4Hc, false)
        );
        assert_eq!(
            XisfCompression::parse("zstd").unwrap(),
            (XisfCompression::Zstd, false)
        );
    }

    #[test]
    fn parse_recognizes_shuffled_codec_names() {
        assert_eq!(
            XisfCompression::parse("zlib+sh").unwrap(),
            (XisfCompression::Zlib, true)
        );
        assert_eq!(
            XisfCompression::parse("lz4+sh").unwrap(),
            (XisfCompression::Lz4, true)
        );
        assert_eq!(
            XisfCompression::parse("lz4hc+sh").unwrap(),
            (XisfCompression::Lz4Hc, true)
        );
        assert_eq!(
            XisfCompression::parse("zstd+sh").unwrap(),
            (XisfCompression::Zstd, true)
        );
    }

    #[test]
    fn parse_rejects_names_outside_the_spec() {
        for name in [
            "",
            "none",
            "gzip",
            "rice",
            "lz4-hc",
            "lz4+hc",
            "zstandard",
            "LZ4",
            "Zlib",
            "zlib:9",
            "zstd2",
            "lz4hc+sh+sh",
            "+sh",
        ] {
            assert!(
                matches!(
                    XisfCompression::parse(name),
                    Err(XisfError::UnsupportedCodec(_))
                ),
                "{name}"
            );
        }
    }

    #[test]
    fn as_str_maps_all_variants() {
        assert_eq!(XisfCompression::None.as_str(), "");
        assert_eq!(XisfCompression::Lz4.as_str(), "lz4");
        assert_eq!(XisfCompression::Lz4Hc.as_str(), "lz4hc");
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
            assert_eq!(
                XisfCompression::parse(variant.as_str()).unwrap(),
                (variant, false)
            );
        }
    }

    #[test]
    fn default_is_none() {
        assert_eq!(XisfCompression::default(), XisfCompression::None);
    }

    #[test]
    fn compression_attr_reads_size_and_item_size() {
        assert_eq!(
            CompressionAttr::parse("zlib:6220800").unwrap(),
            CompressionAttr {
                codec: XisfCompression::Zlib,
                uncompressed_size: 6220800,
                shuffle_item_size: None,
            }
        );
        assert_eq!(
            CompressionAttr::parse("zstd+sh:201326592:4").unwrap(),
            CompressionAttr {
                codec: XisfCompression::Zstd,
                uncompressed_size: 201326592,
                shuffle_item_size: Some(4),
            }
        );
    }

    #[test]
    fn compression_attr_rejects_malformed_fields() {
        for value in [
            "zlib",
            "zlib:",
            "zlib:-1",
            "zlib:12:2",
            "zlib+sh:12",
            "zlib+sh:12:0",
            "zlib+sh:12:x",
            "lz4:12:2:2",
            "zlib:18446744073709551616",
            "zlib:+126",
            "zlib+sh:12:+2",
        ] {
            assert!(
                matches!(
                    CompressionAttr::parse(value),
                    Err(XisfError::InvalidFormat(_))
                ),
                "{value}"
            );
        }
        assert!(matches!(
            CompressionAttr::parse("brotli:4"),
            Err(XisfError::UnsupportedCodec(_))
        ));
    }

    #[test]
    fn parse_subblocks_reads_pairs_in_order() {
        assert_eq!(parse_subblocks("10,20").unwrap(), vec![(10, 20)]);
        assert_eq!(
            parse_subblocks("10,20:5,5:0,0").unwrap(),
            vec![(10, 20), (5, 5), (0, 0)]
        );
    }

    #[test]
    fn parse_subblocks_rejects_malformed_pairs() {
        for value in [
            "", "10", "10,", ",10", "10,20:", "10;20", "10,20,30", "a,b", "1,-1",
        ] {
            assert!(parse_subblocks(value).is_err(), "{value}");
        }
    }
}
