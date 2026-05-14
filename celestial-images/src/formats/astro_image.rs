use super::format::ImageFormat;
use super::kind::ImageKind;
use crate::core::{ImageError, Result};
use crate::fits::compression::CompressionAlgorithm;
use crate::fits::data::DataArray;
use crate::fits::header::Keyword;
use crate::fits::io::writer::FitsWriter;
use crate::xisf::writer::{XisfDataType, XisfWriter};
use celestial_wcs::Wcs;
use std::path::Path;

const DEFAULT_TILE_SIZE: usize = 32;

pub struct AstroImage<'a, T> {
    data: &'a [T],
    dimensions: Vec<usize>,
    keywords: Vec<Keyword>,
    wcs: Option<&'a Wcs>,
    compressed: bool,
    tile_size: Option<(usize, usize)>,
}

impl<'a, T> AstroImage<'a, T> {
    pub fn new(data: &'a [T], dimensions: impl Into<Vec<usize>>) -> Self {
        Self {
            data,
            dimensions: dimensions.into(),
            keywords: Vec::new(),
            wcs: None,
            compressed: true,
            tile_size: None,
        }
    }

    pub fn wcs(mut self, wcs: &'a Wcs) -> Self {
        self.wcs = Some(wcs);
        self
    }

    pub fn keyword(mut self, kw: Keyword) -> Self {
        self.keywords.push(kw);
        self
    }

    pub fn keywords(mut self, keywords: impl IntoIterator<Item = Keyword>) -> Self {
        self.keywords.extend(keywords);
        self
    }

    pub fn compressed(mut self, compress: bool) -> Self {
        self.compressed = compress;
        self
    }

    pub fn tile_size(mut self, width: usize, height: usize) -> Self {
        self.tile_size = Some((width, height));
        self
    }

    pub fn exposure(self, seconds: f64) -> Self {
        self.keyword(Keyword::real("EXPTIME", seconds))
    }

    pub fn temperature(self, celsius: f64) -> Self {
        self.keyword(Keyword::real("CCD-TEMP", celsius))
    }

    pub fn gain(self, gain: f64) -> Self {
        self.keyword(Keyword::real("GAIN", gain))
    }

    pub fn binning(self, x: u32, y: u32) -> Self {
        self.keyword(Keyword::integer("XBINNING", x as i64))
            .keyword(Keyword::integer("YBINNING", y as i64))
    }

    pub fn filter(self, name: &str) -> Self {
        self.keyword(Keyword::string("FILTER", name))
    }

    pub fn object(self, name: &str) -> Self {
        self.keyword(Keyword::string("OBJECT", name))
    }

    pub fn observer(self, name: &str) -> Self {
        self.keyword(Keyword::string("OBSERVER", name))
    }

    pub fn telescope(self, name: &str) -> Self {
        self.keyword(Keyword::string("TELESCOP", name))
    }

    pub fn instrument(self, name: &str) -> Self {
        self.keyword(Keyword::string("INSTRUME", name))
    }

    pub fn date_obs(self, iso8601: &str) -> Self {
        self.keyword(Keyword::string("DATE-OBS", iso8601))
    }

    pub fn image_kind(&self) -> ImageKind {
        ImageKind::from_dimensions(&self.dimensions)
    }

    pub fn write_fits<P: AsRef<Path>>(&self, path: P) -> crate::fits::Result<()>
    where
        T: DataArray + Clone,
    {
        let all_keywords = self.build_keywords();
        let mut writer = FitsWriter::create(path)?;

        if self.compressed {
            let tile_size = self.compute_tile_size();
            writer.write_compressed_image(
                self.data,
                &self.dimensions,
                tile_size,
                CompressionAlgorithm::Rice,
                &all_keywords,
            )
        } else {
            writer.write_primary_image_with_checksum(self.data, &self.dimensions, &all_keywords)
        }
    }

    pub fn write_xisf<P: AsRef<Path>>(&self, path: P) -> crate::xisf::Result<()>
    where
        T: XisfDataType,
    {
        let mut writer = XisfWriter::create(path)?;

        let (width, height, channels) = self.extract_dimensions();
        writer.add_image(self.data, width, height, channels)?;

        for kw in self.build_keywords() {
            writer.add_keyword(kw);
        }

        writer.write()
    }

    pub fn write_to<P: AsRef<Path>>(&self, path: P) -> Result<()>
    where
        T: DataArray + XisfDataType + Clone,
    {
        let path_ref = path.as_ref();
        let ext = path_ref.extension().and_then(|e| e.to_str()).unwrap_or("");

        match ImageFormat::from_extension(ext) {
            Some(ImageFormat::Fits) => self.write_fits(path).map_err(ImageError::Fits),
            Some(ImageFormat::Xisf) => self.write_xisf(path).map_err(ImageError::Xisf),
            #[cfg(feature = "standard-formats")]
            Some(ImageFormat::Png) | Some(ImageFormat::Tiff) => Err(ImageError::UnsupportedFormat),
            None => Err(ImageError::FormatDetectionFailed(format!(
                "Unknown extension: {}",
                ext
            ))),
        }
    }

    fn extract_dimensions(&self) -> (usize, usize, usize) {
        let width = self.dimensions.first().copied().unwrap_or(1);
        let height = self.dimensions.get(1).copied().unwrap_or(1);
        let channels = self.dimensions.get(2).copied().unwrap_or(1);
        (width, height, channels)
    }

    fn build_keywords(&self) -> Vec<Keyword> {
        let mut keywords = Vec::new();

        if let Some(wcs) = self.wcs {
            keywords.extend(wcs.to_keywords().into_iter().map(Keyword::from));
        }

        keywords.extend(self.keywords.iter().cloned());

        keywords
    }

    fn compute_tile_size(&self) -> (usize, usize) {
        if let Some(size) = self.tile_size {
            return size;
        }

        let width = self.dimensions.first().copied().unwrap_or(1);
        let height = self.dimensions.get(1).copied().unwrap_or(1);

        let tile_w = width.min(DEFAULT_TILE_SIZE);
        let tile_h = height.min(DEFAULT_TILE_SIZE);

        (tile_w, tile_h)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::KeywordValue;
    use celestial_wcs::{Projection, WcsBuilder};

    fn tmp_fits() -> tempfile::NamedTempFile {
        tempfile::Builder::new()
            .suffix(".fits")
            .tempfile()
            .unwrap()
    }

    fn tmp_xisf() -> tempfile::NamedTempFile {
        tempfile::Builder::new()
            .suffix(".xisf")
            .tempfile()
            .unwrap()
    }

    fn sample_wcs() -> Wcs {
        WcsBuilder::new()
            .crpix(256.0, 256.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .build()
            .unwrap()
    }

    #[test]
    fn new_starts_with_defaults() {
        let data: Vec<u8> = vec![0; 16];
        let img = AstroImage::new(&data, [4, 4]);
        assert_eq!(img.dimensions, vec![4usize, 4]);
        assert!(img.keywords.is_empty());
        assert!(img.wcs.is_none());
        assert!(img.compressed);
        assert!(img.tile_size.is_none());
    }

    #[test]
    fn image_kind_reflects_dimensions() {
        let data = vec![0u8; 16];
        assert_eq!(
            AstroImage::new(&data, [4, 4]).image_kind(),
            ImageKind::Mono
        );

        let data = vec![0u8; 48];
        assert_eq!(
            AstroImage::new(&data, [4, 4, 3]).image_kind(),
            ImageKind::Rgb
        );
    }

    #[test]
    fn wcs_builder_sets_wcs() {
        let data = vec![0u8; 16];
        let wcs = sample_wcs();
        let img = AstroImage::new(&data, [4, 4]).wcs(&wcs);
        assert!(img.wcs.is_some());
    }

    #[test]
    fn keyword_appends() {
        let data = vec![0u8; 16];
        let img =
            AstroImage::new(&data, [4, 4]).keyword(Keyword::string("OBJECT", "M31"));
        assert_eq!(img.keywords.len(), 1);
        assert_eq!(img.keywords[0].name, "OBJECT");
    }

    #[test]
    fn keywords_extends() {
        let data = vec![0u8; 16];
        let img = AstroImage::new(&data, [4, 4]).keywords(vec![
            Keyword::real("EXPTIME", 30.0),
            Keyword::string("FILTER", "V"),
        ]);
        assert_eq!(img.keywords.len(), 2);
    }

    #[test]
    fn compressed_and_tile_size_builders() {
        let data = vec![0u8; 16];
        let img = AstroImage::new(&data, [4, 4])
            .compressed(false)
            .tile_size(8, 8);
        assert!(!img.compressed);
        assert_eq!(img.tile_size, Some((8, 8)));
    }

    #[test]
    fn fluent_metadata_setters_emit_keywords() {
        let data = vec![0u8; 16];
        let img = AstroImage::new(&data, [4, 4])
            .exposure(30.0)
            .temperature(-10.0)
            .gain(1.5)
            .binning(2, 2)
            .filter("V")
            .object("M31")
            .observer("GA")
            .telescope("Keck")
            .instrument("WFC3")
            .date_obs("2026-04-21T20:00:00");

        let keywords = &img.keywords;
        let find = |name: &str| keywords.iter().find(|k| k.name == name);
        assert_eq!(find("EXPTIME").unwrap().value, Some(KeywordValue::Real(30.0)));
        assert_eq!(find("CCD-TEMP").unwrap().value, Some(KeywordValue::Real(-10.0)));
        assert_eq!(find("GAIN").unwrap().value, Some(KeywordValue::Real(1.5)));
        assert_eq!(find("XBINNING").unwrap().value, Some(KeywordValue::Integer(2)));
        assert_eq!(find("YBINNING").unwrap().value, Some(KeywordValue::Integer(2)));
        assert_eq!(
            find("FILTER").unwrap().value,
            Some(KeywordValue::String("V".into()))
        );
        assert_eq!(
            find("OBJECT").unwrap().value,
            Some(KeywordValue::String("M31".into()))
        );
        assert_eq!(
            find("OBSERVER").unwrap().value,
            Some(KeywordValue::String("GA".into()))
        );
        assert_eq!(
            find("TELESCOP").unwrap().value,
            Some(KeywordValue::String("Keck".into()))
        );
        assert_eq!(
            find("INSTRUME").unwrap().value,
            Some(KeywordValue::String("WFC3".into()))
        );
        assert_eq!(
            find("DATE-OBS").unwrap().value,
            Some(KeywordValue::String("2026-04-21T20:00:00".into()))
        );
    }

    #[test]
    fn write_fits_compressed_produces_readable_file() {
        let tmp = tmp_fits();
        let data: Vec<u8> = (0u8..=255).collect();
        AstroImage::new(&data, [16, 16])
            .write_fits(tmp.path())
            .unwrap();

        assert!(std::fs::metadata(tmp.path()).unwrap().len() > 0);
    }

    #[test]
    fn write_fits_uncompressed_produces_readable_file() {
        let tmp = tmp_fits();
        let data: Vec<u8> = (0u8..=255).collect();
        AstroImage::new(&data, [16, 16])
            .compressed(false)
            .write_fits(tmp.path())
            .unwrap();

        let restored = crate::formats::Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![16usize, 16]);
    }

    #[test]
    fn write_fits_embeds_keywords() {
        let tmp = tmp_fits();
        let data = vec![0u8; 16];
        AstroImage::new(&data, [4, 4])
            .compressed(false)
            .object("M31")
            .exposure(60.5)
            .write_fits(tmp.path())
            .unwrap();

        let restored = crate::formats::Image::open(tmp.path()).unwrap();
        assert_eq!(
            restored.get_keyword("OBJECT").unwrap().value,
            Some(KeywordValue::String("M31".into()))
        );
        assert_eq!(
            restored.get_keyword("EXPTIME").unwrap().value,
            Some(KeywordValue::Real(60.5))
        );
    }

    #[test]
    fn write_xisf_produces_readable_file() {
        let tmp = tmp_xisf();
        let data: Vec<u8> = (0u8..=255).collect();
        AstroImage::new(&data, [16, 16])
            .write_xisf(tmp.path())
            .unwrap();

        let restored = crate::formats::Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![16usize, 16]);
    }

    #[test]
    fn write_to_dispatches_on_extension() {
        let data: Vec<u8> = (0u8..16).collect();

        let tmp = tmp_fits();
        AstroImage::new(&data, [4, 4])
            .compressed(false)
            .write_to(tmp.path())
            .unwrap();
        assert!(tmp.path().exists());

        let tmp = tmp_xisf();
        AstroImage::new(&data, [4, 4]).write_to(tmp.path()).unwrap();
        assert!(tmp.path().exists());
    }

    #[test]
    fn write_to_rejects_unknown_extension() {
        let tmp = tempfile::Builder::new().suffix(".bogus").tempfile().unwrap();
        let data = vec![0u8; 16];
        assert!(AstroImage::new(&data, [4, 4]).write_to(tmp.path()).is_err());
    }

    #[test]
    fn wcs_keywords_are_included_in_fits_output() {
        let tmp = tmp_fits();
        let data = vec![0u8; 16];
        let wcs = sample_wcs();

        AstroImage::new(&data, [4, 4])
            .compressed(false)
            .wcs(&wcs)
            .write_fits(tmp.path())
            .unwrap();

        let restored = crate::formats::Image::open(tmp.path()).unwrap();
        assert!(restored.get_keyword("CRVAL1").is_some());
        assert!(restored.get_keyword("CRVAL2").is_some());
        assert!(restored.get_keyword("CRPIX1").is_some());
        assert!(restored.get_keyword("CRPIX2").is_some());
    }

    #[test]
    fn compute_tile_size_defaults_to_min_of_dims_and_default() {
        let data = vec![0u8; 64];
        let img = AstroImage::new(&data, [8, 8]);
        assert_eq!(img.compute_tile_size(), (8, 8));

        let data = vec![0u8; 10000];
        let img = AstroImage::new(&data, [100, 100]);
        assert_eq!(img.compute_tile_size(), (32, 32));
    }

    #[test]
    fn compute_tile_size_uses_explicit_override() {
        let data = vec![0u8; 10000];
        let img = AstroImage::new(&data, [100, 100]).tile_size(16, 16);
        assert_eq!(img.compute_tile_size(), (16, 16));
    }

    #[test]
    fn extract_dimensions_handles_defaults() {
        let data = vec![0u8; 16];
        let img = AstroImage::new(&data, [4, 4]);
        assert_eq!(img.extract_dimensions(), (4, 4, 1));

        let rgb = AstroImage::new(&data, vec![4usize, 4, 3]);
        assert_eq!(rgb.extract_dimensions(), (4, 4, 3));

        let empty_dims: Vec<usize> = Vec::new();
        let img_empty = AstroImage::new(&data, empty_dims);
        assert_eq!(img_empty.extract_dimensions(), (1, 1, 1));
    }
}

