use crate::core::types::ByteOrder;
use crate::fits::header::keywords::Keyword;
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::block_attrs::BlockDefect;
use crate::xisf::header::checksum::Checksum;
use crate::xisf::header::color_space::ColorSpace;
use crate::xisf::header::compression::XisfCompression;
use crate::xisf::header::location::DataLocation;
use crate::xisf::header::pixel_storage::PixelStorage;
use crate::xisf::header::sample_format::SampleFormat;

#[derive(Debug, Clone)]
pub(crate) struct XisfHeader {
    pub version: String,
    pub images: Vec<ImageInfo>,
    // One list per image, in the order of `images`.
    pub keywords: Vec<Vec<Keyword>>,
    // The FITSKeyword children of the root element.
    pub(crate) root_keywords: Vec<Keyword>,
}

#[derive(Debug, Clone)]
pub struct ImageInfo {
    pub geometry: Vec<usize>,
    pub sample_format: SampleFormat,
    pub bounds: (f64, f64),
    pub color_space: ColorSpace,
    pub pixel_storage: PixelStorage,
    pub(crate) location: DataLocation,
    pub(crate) compression: XisfCompression,
    pub(crate) uncompressed_size: Option<u64>,
    pub(crate) shuffle_item_size: Option<usize>,
    pub(crate) subblocks: Vec<(u64, u64)>,
    pub byte_order: ByteOrder,
    pub(crate) checksum: Option<Checksum>,
    pub(crate) block_defect: Option<BlockDefect>,
}

impl ImageInfo {
    // Byte length of the decoded pixel data. Fails for an empty geometry and
    // for sizes that do not fit in memory.
    pub(crate) fn checked_data_size(&self) -> Result<usize> {
        if self.geometry.is_empty() {
            return Err(XisfError::InvalidGeometry(
                "the image has no dimensions".to_string(),
            ));
        }
        self.geometry
            .iter()
            .try_fold(self.sample_format.bytes_per_sample(), |acc, &d| {
                acc.checked_mul(d)
            })
            .ok_or_else(|| {
                XisfError::InvalidGeometry(format!(
                    "{:?} of {} samples overflows the address space",
                    self.geometry,
                    self.sample_format.as_str()
                ))
            })
    }

    pub(crate) fn format_bounds(&self) -> String {
        format!("{}:{}", self.bounds.0, self.bounds.1)
    }

    // Normal storage interleaves the channels of each pixel (XISF 1.0
    // section 8.5.3.2); one channel is the same in either order.
    pub(crate) fn stores_interleaved(&self) -> bool {
        self.pixel_storage == PixelStorage::Normal && self.num_channels() > 1
    }

    pub(crate) fn num_channels(&self) -> usize {
        if self.geometry.len() >= 3 {
            self.geometry[2]
        } else {
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::{Keyword, KeywordValue};

    fn info(geometry: Vec<usize>, sample_format: SampleFormat) -> ImageInfo {
        ImageInfo {
            geometry,
            sample_format,
            bounds: (0.0, 1.0),
            color_space: ColorSpace::Gray,
            pixel_storage: PixelStorage::Planar,
            location: DataLocation::new(0, 0),
            compression: XisfCompression::None,
            uncompressed_size: None,
            shuffle_item_size: None,
            subblocks: Vec::new(),
            byte_order: ByteOrder::LittleEndian,
            checksum: None,
            block_defect: None,
        }
    }

    #[test]
    fn format_bounds_separates_with_colon() {
        let mut img = info(vec![10, 10], SampleFormat::UInt8);
        img.bounds = (0.0, 65535.0);
        assert_eq!(img.format_bounds(), "0:65535");

        img.bounds = (-1.5, 1.5);
        assert_eq!(img.format_bounds(), "-1.5:1.5");
    }

    #[test]
    fn num_channels_is_one_for_two_dim() {
        let img = info(vec![100, 100], SampleFormat::UInt8);
        assert_eq!(img.num_channels(), 1);
    }

    #[test]
    fn num_channels_reads_third_dimension() {
        let img = info(vec![100, 100, 3], SampleFormat::UInt8);
        assert_eq!(img.num_channels(), 3);

        let img = info(vec![100, 100, 7], SampleFormat::UInt8);
        assert_eq!(img.num_channels(), 7);
    }

    #[test]
    fn num_channels_one_for_short_geometry() {
        let img = info(vec![], SampleFormat::UInt8);
        assert_eq!(img.num_channels(), 1);

        let img = info(vec![10], SampleFormat::UInt8);
        assert_eq!(img.num_channels(), 1);
    }

    #[test]
    fn image_info_preserves_all_fields() {
        let img = ImageInfo {
            geometry: vec![1920, 1080, 3],
            sample_format: SampleFormat::UInt16,
            bounds: (0.0, 65535.0),
            color_space: ColorSpace::Rgb,
            pixel_storage: PixelStorage::Normal,
            location: DataLocation::new(1024, 12441600),
            compression: XisfCompression::Lz4,
            uncompressed_size: Some(12441600),
            shuffle_item_size: Some(2),
            subblocks: vec![(100, 200)],
            byte_order: ByteOrder::BigEndian,
            checksum: None,
            block_defect: None,
        };
        assert_eq!(img.geometry, vec![1920, 1080, 3]);
        assert!(matches!(img.sample_format, SampleFormat::UInt16));
        assert_eq!(img.bounds, (0.0, 65535.0));
        assert!(matches!(img.color_space, ColorSpace::Rgb));
        assert_eq!(img.pixel_storage, PixelStorage::Normal);
        assert_eq!(img.location.offset, 1024);
        assert_eq!(img.location.size, 12441600);
        assert_eq!(img.compression, XisfCompression::Lz4);
        assert_eq!(img.uncompressed_size, Some(12441600));
        assert_eq!(img.shuffle_item_size, Some(2));
        assert_eq!(img.subblocks, vec![(100, 200)]);
        assert_eq!(img.byte_order, ByteOrder::BigEndian);
    }

    #[test]
    fn sizes_do_not_overflow_on_huge_geometry() {
        let huge = info(vec![1 << 32, 1 << 32, 3], SampleFormat::UInt16);
        assert!(matches!(
            huge.checked_data_size(),
            Err(XisfError::InvalidGeometry(_))
        ));
    }

    #[test]
    fn checked_data_size_matches_data_size_and_needs_a_dimension() {
        let rgb = info(vec![1920, 1080, 3], SampleFormat::Float32);
        assert_eq!(rgb.checked_data_size().unwrap(), 1920 * 1080 * 3 * 4);
        assert_eq!(
            info(vec![10], SampleFormat::UInt8)
                .checked_data_size()
                .unwrap(),
            10
        );
        assert!(info(vec![], SampleFormat::UInt8)
            .checked_data_size()
            .is_err());
    }

    #[test]
    fn xisf_header_holds_version_images_and_keywords() {
        let header = XisfHeader {
            version: "1.0".to_string(),
            images: vec![info(vec![10, 10], SampleFormat::UInt8)],
            keywords: vec![vec![
                Keyword::string("TELESCOP", "Hubble").with_comment("name")
            ]],
            root_keywords: Vec::new(),
        };
        assert_eq!(header.version, "1.0");
        assert_eq!(header.images.len(), 1);
        assert_eq!(header.keywords[0].len(), 1);
        assert_eq!(header.keywords[0][0].name, "TELESCOP");
        assert_eq!(
            header.keywords[0][0].value,
            Some(KeywordValue::String("Hubble".to_string()))
        );
    }
}
