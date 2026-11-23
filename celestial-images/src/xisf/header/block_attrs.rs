use crate::core::types::ByteOrder;
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::checksum::Checksum;
use crate::xisf::header::compression::{parse_subblocks, CompressionAttr};
use crate::xisf::header::image_info::ImageInfo;

// A data block attribute that could not be parsed. XISF 1.0 section 10.6
// makes such a block unavailable while the rest of the file stays readable,
// so the error is kept here and returned when the block is read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum BlockDefect {
    UnsupportedCodec(String),
    UnsupportedChecksum(String),
    Invalid(String),
}

impl BlockDefect {
    pub(crate) fn from_error(error: XisfError) -> Self {
        match error {
            XisfError::UnsupportedCodec(name) => Self::UnsupportedCodec(name),
            XisfError::UnsupportedChecksum(name) => Self::UnsupportedChecksum(name),
            XisfError::InvalidFormat(message) => Self::Invalid(message),
            other => Self::Invalid(other.to_string()),
        }
    }

    pub(crate) fn to_error(&self) -> XisfError {
        match self {
            Self::UnsupportedCodec(name) => XisfError::UnsupportedCodec(name.clone()),
            Self::UnsupportedChecksum(name) => XisfError::UnsupportedChecksum(name.clone()),
            Self::Invalid(message) => XisfError::InvalidFormat(message.clone()),
        }
    }
}

// Attributes that say how an image data block is encoded. A value that fails
// to parse makes only this block unavailable (XISF 1.0 section 10.6).
#[derive(Default)]
pub(crate) struct BlockAttrs {
    compression: Option<CompressionAttr>,
    subblocks: Vec<(u64, u64)>,
    byte_order: Option<ByteOrder>,
    checksum: Option<Checksum>,
    defect: Option<BlockDefect>,
}

impl BlockAttrs {
    pub(crate) fn apply(&mut self, key: &[u8], value: &str) {
        let parsed = match key {
            b"compression" => CompressionAttr::parse(value).map(|c| self.compression = Some(c)),
            b"subblocks" => parse_subblocks(value).map(|s| self.subblocks = s),
            b"byteOrder" => parse_byte_order(value).map(|b| self.byte_order = Some(b)),
            b"checksum" => Checksum::parse(value).map(|c| self.checksum = Some(c)),
            _ => Ok(()),
        };
        if let Err(error) = parsed {
            self.record(error);
        }
    }

    // Keeps the first problem, which is the one the reader reports.
    pub(crate) fn record(&mut self, error: XisfError) {
        self.defect.get_or_insert(BlockDefect::from_error(error));
    }

    pub(crate) fn into_image_info(self, mut info: ImageInfo) -> ImageInfo {
        self.apply_to(&mut info);
        info
    }

    // True when the attributes give no compression and have no defect.
    pub(crate) fn is_plain(&self) -> bool {
        self.compression.is_none() && self.subblocks.is_empty() && self.defect.is_none()
    }

    // Sets the attributes given here and keeps the others. An embedded
    // block has its compression and subblocks on the Data element and its
    // other attributes on the Image element (XISF 1.0 section 10.6).
    // Subblocks describe the compressed block, so they must be on the
    // element that has the compression attribute.
    pub(crate) fn apply_to(mut self, info: &mut ImageInfo) {
        if !self.subblocks.is_empty() && self.compression.is_none() {
            let message = "A subblocks attribute needs a compression attribute on the same element";
            self.record(XisfError::InvalidFormat(message.to_string()));
        }
        if let Some(attr) = self.compression {
            info.compression = attr.codec;
            info.uncompressed_size = Some(attr.uncompressed_size);
            info.shuffle_item_size = attr.shuffle_item_size;
        }
        if !self.subblocks.is_empty() {
            info.subblocks = std::mem::take(&mut self.subblocks);
        }
        info.byte_order = self.byte_order.unwrap_or(info.byte_order);
        info.checksum = self.checksum.or(info.checksum.take());
        if info.block_defect.is_none() {
            info.block_defect = self.defect;
        }
    }
}

fn parse_byte_order(value: &str) -> Result<ByteOrder> {
    match value {
        "little" => Ok(ByteOrder::LittleEndian),
        "big" => Ok(ByteOrder::BigEndian),
        _ => Err(XisfError::InvalidFormat(format!(
            "Invalid byteOrder attribute: \"{value}\""
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xisf::header::color_space::ColorSpace;
    use crate::xisf::header::compression::XisfCompression;
    use crate::xisf::header::location::DataLocation;
    use crate::xisf::header::pixel_storage::PixelStorage;
    use crate::xisf::header::sample_format::SampleFormat;

    fn plain_info() -> ImageInfo {
        ImageInfo {
            geometry: vec![2, 2],
            sample_format: SampleFormat::UInt16,
            bounds: (0.0, 1.0),
            color_space: ColorSpace::Gray,
            pixel_storage: PixelStorage::Planar,
            location: DataLocation::new(0, 8),
            compression: XisfCompression::None,
            uncompressed_size: None,
            shuffle_item_size: None,
            subblocks: Vec::new(),
            byte_order: ByteOrder::LittleEndian,
            checksum: None,
            block_defect: None,
        }
    }

    fn apply_all(attrs: &[(&str, &str)]) -> ImageInfo {
        let mut block = BlockAttrs::default();
        for (key, value) in attrs {
            block.apply(key.as_bytes(), value);
        }
        block.into_image_info(plain_info())
    }

    #[test]
    fn valid_attributes_fill_the_image_info() {
        let info = apply_all(&[
            ("compression", "lz4hc+sh:8:2"),
            ("subblocks", "5,8"),
            ("byteOrder", "big"),
        ]);
        assert_eq!(info.compression, XisfCompression::Lz4Hc);
        assert_eq!(info.uncompressed_size, Some(8));
        assert_eq!(info.shuffle_item_size, Some(2));
        assert_eq!(info.subblocks, vec![(5, 8)]);
        assert_eq!(info.byte_order, ByteOrder::BigEndian);
        assert_eq!(info.block_defect, None);
    }

    #[test]
    fn unknown_codec_is_kept_as_a_defect() {
        let info = apply_all(&[("compression", "brotli:4")]);
        assert_eq!(
            info.block_defect,
            Some(BlockDefect::UnsupportedCodec("brotli".to_string()))
        );
        assert!(matches!(
            info.block_defect.map(|d| d.to_error()),
            Some(XisfError::UnsupportedCodec(_))
        ));
    }

    #[test]
    fn the_first_malformed_attribute_is_reported() {
        let info = apply_all(&[("byteOrder", "middle"), ("compression", "zlib")]);
        assert_eq!(
            info.block_defect,
            Some(BlockDefect::Invalid(
                "Invalid byteOrder attribute: \"middle\"".to_string()
            ))
        );
    }

    #[test]
    fn subblocks_without_compression_is_a_defect() {
        let info = apply_all(&[("subblocks", "8,8")]);
        assert!(matches!(info.block_defect, Some(BlockDefect::Invalid(_))));
    }

    #[test]
    fn a_checksum_is_parsed_into_the_image_info() {
        let digest = "12dada1fff4d4787ade3333147202c3b443e376f";
        let info = apply_all(&[("checksum", &format!("sha1:{digest}"))]);
        assert_eq!(info.block_defect, None);
        let checksum = info.checksum.unwrap();
        assert!(checksum.verify(&[1, 2, 3, 4]).is_ok());
        assert!(checksum.verify(&[1, 2, 3, 5]).is_err());
    }

    #[test]
    fn a_malformed_or_unsupported_checksum_is_a_defect() {
        let info = apply_all(&[("checksum", "sha-1:00")]);
        assert_eq!(
            info.block_defect,
            Some(BlockDefect::Invalid(
                "Invalid checksum attribute \"sha-1:00\": a sha-1 digest is 40 hexadecimal digits"
                    .to_string()
            ))
        );
        let info = apply_all(&[("checksum", "sha3-512:00")]);
        assert_eq!(
            info.block_defect,
            Some(BlockDefect::UnsupportedChecksum("sha3-512".to_string()))
        );
        assert!(matches!(
            info.block_defect.map(|d| d.to_error()),
            Some(XisfError::UnsupportedChecksum(_))
        ));
    }

    #[test]
    fn other_attributes_are_ignored() {
        let info = apply_all(&[("id", "main"), ("encoding", "base64")]);
        assert_eq!(info.block_defect, None);
        assert_eq!(info.compression, XisfCompression::None);
        assert_eq!(info.checksum, None);
    }

    // The Data element of an embedded block adds its compression to the
    // attributes of the Image element.
    #[test]
    fn apply_to_keeps_attributes_it_does_not_set() {
        let mut info = apply_all(&[("byteOrder", "big"), ("checksum", "sha1:00")]);
        let mut data = BlockAttrs::default();
        data.apply(b"compression", "zlib:8");
        data.apply_to(&mut info);
        assert_eq!(info.compression, XisfCompression::Zlib);
        assert_eq!(info.byte_order, ByteOrder::BigEndian);
        assert!(matches!(info.block_defect, Some(BlockDefect::Invalid(_))));
    }
}
