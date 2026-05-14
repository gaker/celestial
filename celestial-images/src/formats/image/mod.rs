mod fits;
mod metadata;
mod processing;

#[cfg(feature = "standard-formats")]
mod png;
#[cfg(feature = "standard-formats")]
mod tiff;

mod xisf;

use super::format::ImageFormat;
use super::kind::ImageKind;
use super::pixel_data::PixelData;
use crate::core::{ImageError, Result};
use crate::fits::header::Keyword;
use crate::xisf::XisfProperty;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Image {
    pub pixels: PixelData,
    pub dimensions: Vec<usize>,
    pub keywords: Vec<Keyword>,
    pub xisf_properties: Vec<XisfProperty>,
}

impl Image {
    pub fn new(pixels: PixelData, dimensions: impl Into<Vec<usize>>) -> Self {
        Self {
            pixels,
            dimensions: dimensions.into(),
            keywords: Vec::new(),
            xisf_properties: Vec::new(),
        }
    }

    pub fn add_xisf_property(&mut self, property: XisfProperty) {
        self.xisf_properties.push(property);
    }

    pub fn add_xisf_properties(
        &mut self,
        properties: impl IntoIterator<Item = XisfProperty>,
    ) {
        self.xisf_properties.extend(properties);
    }

    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();
        let extension = path_ref
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("");

        let format = match ImageFormat::from_extension(extension) {
            Some(f) => f,
            None => {
                let mut file = std::fs::File::open(path_ref)?;
                ImageFormat::detect(&mut file)?
            }
        };

        match format {
            ImageFormat::Fits => Self::open_fits(path_ref),
            ImageFormat::Xisf => Self::open_xisf(path_ref),
            #[cfg(feature = "standard-formats")]
            ImageFormat::Png => Self::open_png(path_ref),
            #[cfg(feature = "standard-formats")]
            ImageFormat::Tiff => Self::open_tiff(path_ref),
        }
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path_ref = path.as_ref();
        let extension = path_ref
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("");

        match ImageFormat::from_extension(extension) {
            Some(ImageFormat::Fits) => self.save_fits(path_ref),
            Some(ImageFormat::Xisf) => self.save_xisf(path_ref),
            #[cfg(feature = "standard-formats")]
            Some(ImageFormat::Png) => self.save_png(path_ref),
            #[cfg(feature = "standard-formats")]
            Some(ImageFormat::Tiff) => self.save_tiff(path_ref),
            None => Err(ImageError::UnsupportedFormat),
        }
    }

    pub fn width(&self) -> usize {
        self.dimensions.first().copied().unwrap_or(0)
    }

    pub fn height(&self) -> usize {
        self.dimensions.get(1).copied().unwrap_or(1)
    }

    pub fn channels(&self) -> usize {
        self.dimensions.get(2).copied().unwrap_or(1)
    }

    pub fn is_rgb(&self) -> bool {
        self.channels() == 3
    }

    pub fn kind(&self) -> ImageKind {
        ImageKind::from_dimensions(&self.dimensions)
    }

    pub fn get_keyword(&self, name: &str) -> Option<&Keyword> {
        self.keywords.iter().find(|k| k.name == name)
    }

    pub fn set_keyword(&mut self, kw: Keyword) {
        if let Some(existing) = self.keywords.iter_mut().find(|k| k.name == kw.name) {
            *existing = kw;
        } else {
            self.keywords.push(kw);
        }
    }

    pub fn remove_keyword(&mut self, name: &str) {
        self.keywords.retain(|k| k.name != name);
    }

    pub(super) fn extract_dimensions(&self) -> (usize, usize, usize) {
        let width = self.dimensions.first().copied().unwrap_or(1);
        let height = self.dimensions.get(1).copied().unwrap_or(1);
        let channels = self.dimensions.get(2).copied().unwrap_or(1);
        (width, height, channels)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::KeywordValue;
    use crate::xisf::XisfProperty;

    fn mono_u8(width: usize, height: usize) -> Image {
        let pixels = PixelData::U8(vec![0u8; width * height]);
        Image::new(pixels, vec![width, height])
    }

    fn rgb_u8(width: usize, height: usize) -> Image {
        let pixels = PixelData::U8(vec![0u8; width * height * 3]);
        Image::new(pixels, vec![width, height, 3])
    }

    #[test]
    fn new_starts_with_empty_keywords_and_properties() {
        let img = mono_u8(10, 10);
        assert!(img.keywords.is_empty());
        assert!(img.xisf_properties.is_empty());
        assert_eq!(img.dimensions, vec![10usize, 10]);
        assert_eq!(img.pixels.len(), 100);
    }

    #[test]
    fn dimension_accessors_match_shape() {
        let img = mono_u8(320, 240);
        assert_eq!(img.width(), 320);
        assert_eq!(img.height(), 240);
        assert_eq!(img.channels(), 1);
        assert!(!img.is_rgb());

        let rgb = rgb_u8(320, 240);
        assert_eq!(rgb.width(), 320);
        assert_eq!(rgb.height(), 240);
        assert_eq!(rgb.channels(), 3);
        assert!(rgb.is_rgb());
    }

    #[test]
    fn dimension_accessors_handle_missing_axes() {
        let img = Image::new(PixelData::U8(vec![0]), Vec::<usize>::new());
        assert_eq!(img.width(), 0);
        assert_eq!(img.height(), 1);
        assert_eq!(img.channels(), 1);
    }

    #[test]
    fn kind_matches_dimensionality() {
        assert_eq!(mono_u8(10, 10).kind(), ImageKind::Mono);
        assert_eq!(rgb_u8(10, 10).kind(), ImageKind::Rgb);
        let cube = Image::new(PixelData::U8(vec![0; 80]), vec![4usize, 4, 5]);
        assert_eq!(cube.kind(), ImageKind::Cube);
    }

    #[test]
    fn set_keyword_adds_when_absent() {
        let mut img = mono_u8(2, 2);
        img.set_keyword(Keyword::real("EXPTIME", 30.0));
        assert_eq!(
            img.get_keyword("EXPTIME").unwrap().value,
            Some(KeywordValue::Real(30.0))
        );
    }

    #[test]
    fn set_keyword_replaces_when_present() {
        let mut img = mono_u8(2, 2);
        img.set_keyword(Keyword::real("EXPTIME", 10.0));
        img.set_keyword(Keyword::real("EXPTIME", 30.0));
        assert_eq!(img.keywords.len(), 1);
        assert_eq!(
            img.get_keyword("EXPTIME").unwrap().value,
            Some(KeywordValue::Real(30.0))
        );
    }

    #[test]
    fn remove_keyword_drops_matching_name() {
        let mut img = mono_u8(2, 2);
        img.set_keyword(Keyword::real("A", 1.0));
        img.set_keyword(Keyword::real("B", 2.0));
        img.remove_keyword("A");
        assert!(img.get_keyword("A").is_none());
        assert!(img.get_keyword("B").is_some());
    }

    #[test]
    fn get_keyword_returns_none_for_missing_name() {
        let img = mono_u8(2, 2);
        assert!(img.get_keyword("NONEXISTENT").is_none());
    }

    #[test]
    fn add_xisf_property_appends() {
        let mut img = mono_u8(2, 2);
        img.add_xisf_property(XisfProperty::int32("Test:A", 1));
        img.add_xisf_property(XisfProperty::int32("Test:B", 2));
        assert_eq!(img.xisf_properties.len(), 2);
        assert_eq!(img.xisf_properties[0].id, "Test:A");
        assert_eq!(img.xisf_properties[1].id, "Test:B");
    }

    #[test]
    fn add_xisf_properties_extends() {
        let mut img = mono_u8(2, 2);
        img.add_xisf_properties(vec![
            XisfProperty::int32("Test:X", 7),
            XisfProperty::int32("Test:Y", 9),
        ]);
        assert_eq!(img.xisf_properties.len(), 2);
    }

    #[test]
    fn extract_dimensions_handles_mono_and_rgb() {
        assert_eq!(mono_u8(10, 20).extract_dimensions(), (10, 20, 1));
        assert_eq!(rgb_u8(10, 20).extract_dimensions(), (10, 20, 3));
    }

    #[test]
    fn save_rejects_unknown_extension() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("bogus");
        let img = mono_u8(4, 4);
        assert!(img.save(&path).is_err());
    }

    #[test]
    fn open_rejects_nonexistent_file_with_unknown_extension() {
        let result = Image::open("/nonexistent/path/file.unknown");
        assert!(result.is_err());
    }
}
