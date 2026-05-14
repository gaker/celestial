use crate::fits::header::Keyword;
use crate::xisf::header::{
    ColorSpace, DataLocation, PixelStorage, SampleFormat, XisfCompression,
};

#[derive(Debug, Clone)]
pub struct XisfHeader {
    pub version: String,
    pub images: Vec<ImageInfo>,
    pub keywords: Vec<Keyword>,
}

#[derive(Debug, Clone)]
pub struct ImageInfo {
    pub geometry: Vec<usize>,
    pub sample_format: SampleFormat,
    pub bounds: (f64, f64),
    pub color_space: ColorSpace,
    pub pixel_storage: PixelStorage,
    pub location: DataLocation,
    pub compression: XisfCompression,
    pub uncompressed_size: Option<u64>,
}

impl ImageInfo {
    pub fn data_size(&self) -> u64 {
        if self.geometry.is_empty() {
            return 0;
        }
        let total_pixels: u64 = self.geometry.iter().map(|&d| d as u64).product();
        total_pixels * self.sample_format.bytes_per_sample() as u64
    }

    pub fn format_bounds(&self) -> String {
        format!("{}:{}", self.bounds.0, self.bounds.1)
    }

    pub fn num_channels(&self) -> usize {
        if self.geometry.len() >= 3 {
            self.geometry[2]
        } else {
            1
        }
    }

    pub fn pixels_per_channel(&self) -> usize {
        if self.geometry.len() >= 2 {
            self.geometry[0] * self.geometry[1]
        } else {
            0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Keyword, KeywordValue};

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
        }
    }

    #[test]
    fn data_size_is_pixels_times_bytes_per_sample() {
        let gray = info(vec![1920, 1080], SampleFormat::UInt16);
        assert_eq!(gray.data_size(), 1920 * 1080 * 2);

        let rgb = info(vec![1920, 1080, 3], SampleFormat::UInt8);
        assert_eq!(rgb.data_size(), 1920 * 1080 * 3);

        let f64_cube = info(vec![100, 100, 4], SampleFormat::Float64);
        assert_eq!(f64_cube.data_size(), 100 * 100 * 4 * 8);
    }

    #[test]
    fn data_size_empty_geometry_is_zero() {
        let empty = info(vec![], SampleFormat::UInt8);
        assert_eq!(empty.data_size(), 0);

        let empty_f64 = info(vec![], SampleFormat::Float64);
        assert_eq!(empty_f64.data_size(), 0);
    }

    #[test]
    fn data_size_scales_with_float_formats() {
        let f32_img = info(vec![512, 512], SampleFormat::Float32);
        assert_eq!(f32_img.data_size(), 512 * 512 * 4);

        let f64_img = info(vec![512, 512], SampleFormat::Float64);
        assert_eq!(f64_img.data_size(), 512 * 512 * 8);
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
    fn pixels_per_channel_multiplies_first_two_dims() {
        let img = info(vec![100, 50], SampleFormat::UInt8);
        assert_eq!(img.pixels_per_channel(), 5000);

        let img = info(vec![100, 50, 3], SampleFormat::UInt8);
        assert_eq!(img.pixels_per_channel(), 5000);
    }

    #[test]
    fn pixels_per_channel_zero_for_under_two_dims() {
        assert_eq!(info(vec![], SampleFormat::UInt8).pixels_per_channel(), 0);
        assert_eq!(info(vec![10], SampleFormat::UInt8).pixels_per_channel(), 0);
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
    }

    #[test]
    fn xisf_header_holds_version_images_and_keywords() {
        let header = XisfHeader {
            version: "1.0".to_string(),
            images: vec![info(vec![10, 10], SampleFormat::UInt8)],
            keywords: vec![Keyword::string("TELESCOP", "Hubble").with_comment("name")],
        };
        assert_eq!(header.version, "1.0");
        assert_eq!(header.images.len(), 1);
        assert_eq!(header.keywords.len(), 1);
        assert_eq!(header.keywords[0].name, "TELESCOP");
        assert_eq!(
            header.keywords[0].value,
            Some(KeywordValue::String("Hubble".to_string()))
        );
    }
}
