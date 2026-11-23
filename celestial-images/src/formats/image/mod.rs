mod fits;
#[cfg(test)]
mod fits_compressed_tests;
#[cfg(test)]
mod fits_compressed_write_tests;
mod fits_scaling;
mod metadata;
pub(crate) mod processing;

#[cfg(feature = "standard-formats")]
mod png;
#[cfg(feature = "standard-formats")]
mod tiff;

pub(super) mod xisf;
#[cfg(test)]
mod xisf_wcs_roundtrip_tests;

use super::format::ImageFormat;
use super::kind::ImageKind;
use super::pixel_data::{with_samples, PixelData};
use crate::core::errors::{ImageError, Result};
use crate::core::sample_order::to_planar;
use crate::core::shape::plane_shape;
use crate::fits::header::keywords::{find_keyword, find_last_keyword, set_keyword, Keyword};
use crate::fits::header::parser::Field;
use crate::xisf::header::property::XisfProperty;
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

    // An image read from a file. An Image holds RGB samples interleaved, and
    // the planes of any other image one after another, so samples stored in
    // the other order are reordered.
    pub(super) fn from_parts(
        pixels: PixelData,
        dimensions: Vec<usize>,
        keywords: Vec<Keyword>,
        xisf_properties: Vec<XisfProperty>,
        stored_interleaved: bool,
    ) -> Self {
        let mut img = Self {
            pixels,
            dimensions,
            keywords,
            xisf_properties,
        };
        match (stored_interleaved, img.has_rgb_planes()) {
            (false, true) => img.planar_to_interleaved(),
            (true, false) => {
                let channels = img.channels();
                with_samples!(&mut img.pixels, data => *data = to_planar(data, channels))
            }
            _ => {}
        }
        img
    }

    pub fn add_xisf_property(&mut self, property: XisfProperty) {
        self.xisf_properties.push(property);
    }

    pub fn add_xisf_properties(&mut self, properties: impl IntoIterator<Item = XisfProperty>) {
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
        plane_shape(&self.dimensions).0
    }

    pub fn height(&self) -> usize {
        plane_shape(&self.dimensions).1
    }

    pub fn channels(&self) -> usize {
        plane_shape(&self.dimensions).2
    }

    pub fn is_rgb(&self) -> bool {
        self.channels() == 3
    }

    pub fn kind(&self) -> ImageKind {
        ImageKind::from_dimensions(&self.dimensions)
    }

    pub fn get_keyword(&self, name: &str) -> Option<&Keyword> {
        self.keywords.get(find_keyword(&self.keywords, name)?)
    }

    pub fn get(&self, name: &str) -> Field<'_> {
        Field::from_value(self.get_keyword(name).and_then(|k| k.value.as_ref()))
    }

    // The last card of a repeated name, which data and WCS readers apply.
    pub fn get_last(&self, name: &str) -> Field<'_> {
        let keyword = find_last_keyword(&self.keywords, name).and_then(|i| self.keywords.get(i));
        Field::from_value(keyword.and_then(|k| k.value.as_ref()))
    }

    // Replaces every card of the same name with one card, except commentary
    // keywords, which are appended.
    pub fn set_keyword(&mut self, kw: Keyword) {
        set_keyword(&mut self.keywords, kw);
    }

    pub fn remove_keyword(&mut self, name: &str) {
        let Some(stored) = self.get_keyword(name).map(|k| k.name.clone()) else {
            return;
        };
        self.keywords.retain(|k| k.name != stored);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::compression::CompressionAlgorithm;
    use crate::fits::errors::FitsError;
    use crate::fits::header::keywords::KeywordValue;
    use crate::fits::io::reader::FitsFile;
    use crate::xisf::header::property::XisfProperty;

    fn mono_u8(width: usize, height: usize) -> Image {
        let pixels = PixelData::U8(vec![0u8; width * height]);
        Image::new(pixels, vec![width, height])
    }

    fn rgb_u8(width: usize, height: usize) -> Image {
        let pixels = PixelData::U8(vec![0u8; width * height * 3]);
        Image::new(pixels, vec![width, height, 3])
    }

    fn no_axes(samples: Vec<u8>) -> Image {
        Image::new(PixelData::U8(samples), Vec::<usize>::new())
    }

    #[test]
    fn image_with_no_axes_is_refused_before_any_encoder_runs() {
        let dir = tempfile::tempdir().unwrap();
        let mut extensions = vec!["xisf"];
        if cfg!(feature = "standard-formats") {
            extensions.extend(["png", "tiff"]);
        }
        for samples in [vec![], vec![7]] {
            for ext in &extensions {
                let path = dir.path().join(format!("empty.{ext}"));
                let err = no_axes(samples.clone()).save(&path).unwrap_err();
                assert!(
                    matches!(&err, ImageError::InvalidDimensions(m) if m == "the image has no axes"),
                    "{ext}: {err}"
                );
                assert!(!path.exists(), "{ext}");
            }
        }
    }

    // FITS writes an image with no axes as NAXIS = 0, which declares no
    // data, so only an image without samples has that form.
    #[test]
    fn fits_save_of_an_image_with_no_axes_keeps_the_naxis_0_rules() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.fits");
        no_axes(vec![]).save(&path).unwrap();
        let mut fits = FitsFile::open(&path).unwrap();
        assert_eq!(
            fits.get_header_value(0, "NAXIS").unwrap().as_deref(),
            Some("0")
        );

        let path = dir.path().join("one.fits");
        let err = no_axes(vec![7]).save(&path).unwrap_err();
        assert!(
            matches!(&err, ImageError::Fits(FitsError::InvalidFormat(m))
                if m == "data has 1 values but NAXIS = 0 declares no data"),
            "{err}"
        );
        assert!(!path.exists());

        let path = dir.path().join("tiled.fits");
        let err = no_axes(vec![])
            .save_fits_compressed(&path, CompressionAlgorithm::Rice, None)
            .unwrap_err();
        assert!(
            matches!(&err, ImageError::Fits(FitsError::InvalidFormat(m))
                if m == "a compressed image needs at least one axis"),
            "{err}"
        );
        assert!(!path.exists());
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
    fn get_returns_typed_field() {
        let mut img = mono_u8(2, 2);
        img.set_keyword(Keyword::integer("NAXIS", 2));
        img.set_keyword(Keyword::real("EXPTIME", 30.0));
        img.set_keyword(Keyword::string("OBJECT", "M31"));
        img.set_keyword(Keyword::logical("SIMPLE", true));

        assert_eq!(img.get("NAXIS").as_i64(), Some(2));
        assert_eq!(img.get("NAXIS").as_f64(), Some(2.0));
        assert_eq!(img.get("EXPTIME").as_f64(), Some(30.0));
        assert_eq!(img.get("EXPTIME").as_i64(), None);
        assert_eq!(img.get("OBJECT").as_str(), Some("M31"));
        assert_eq!(img.get("SIMPLE").as_bool(), Some(true));
    }

    #[test]
    fn get_missing_keyword_returns_empty_field() {
        let img = mono_u8(2, 2);
        let field = img.get("NONEXISTENT");
        assert!(!field.exists());
        assert_eq!(field.as_i64(), None);
        assert_eq!(field.value(), None);
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
        assert_eq!(plane_shape(&mono_u8(10, 20).dimensions), (10, 20, 1));
        assert_eq!(plane_shape(&rgb_u8(10, 20).dimensions), (10, 20, 3));
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

    #[test]
    fn keyword_lookups_fall_back_to_hierarch_for_non_standard_names() {
        let mut img = mono_u8(1, 1);
        img.set_keyword(Keyword::real("HIERARCH ESO DET DIT", 1.5));
        img.set_keyword(Keyword::integer("HIERARCH NAXIS1", 3));
        assert_eq!(img.get("ESO DET DIT").as_f64(), Some(1.5));
        assert!(img.get_keyword("NAXIS1").is_none());

        img.set_keyword(Keyword::real("ESO DET DIT", 2.5));
        assert_eq!(img.keywords.len(), 2);
        assert_eq!(img.get("HIERARCH ESO DET DIT").as_f64(), Some(2.5));

        img.remove_keyword("ESO DET DIT");
        assert_eq!(img.keywords.len(), 1);
        img.remove_keyword("NAXIS1");
        assert_eq!(img.keywords.len(), 1);
    }

    #[test]
    fn set_keyword_with_hierarch_prefix_replaces_bare_name() {
        for (bare, full) in [("ESO DET", "HIERARCH ESO DET"), ("foo", "HIERARCH foo")] {
            let mut img = mono_u8(1, 1);
            img.set_keyword(Keyword::integer(bare, 1));
            img.set_keyword(Keyword::integer(full, 2));
            assert_eq!(img.keywords.len(), 1, "{}", full);
            assert_eq!(img.get(bare).as_i64(), Some(2), "{}", full);
            assert_eq!(img.get(full).as_i64(), Some(2), "{}", full);
        }
        let mut img = mono_u8(1, 1);
        img.set_keyword(Keyword::integer("EXPTIME", 1));
        img.set_keyword(Keyword::integer("HIERARCH EXPTIME", 2));
        assert_eq!(img.keywords.len(), 2);
        assert_eq!(img.get("EXPTIME").as_i64(), Some(1));
    }
}
