use super::atomic_write::write_atomically;
use super::format::ImageFormat;
use super::image::xisf::write_xisf_file;
use super::kind::ImageKind;
use crate::core::errors::{ImageError, Result};
use crate::core::shape::plane_shape_with_axes;
use crate::fits::compression::CompressionAlgorithm;
use crate::fits::data::array::DataArray;
use crate::fits::header::keywords::{set_keyword, Keyword};
use crate::fits::header::metadata;
use crate::fits::io::writer::{compressed_image_bytes, primary_image_bytes};
use crate::xisf::errors::XisfError;
use crate::xisf::writer::data_type::XisfDataType;
use crate::xisf::writer::XisfWriter;
use celestial_wcs::builder::Wcs;
use std::io::Cursor;
use std::path::Path;

pub struct AstroImage<'a, T> {
    data: &'a [T],
    dimensions: Vec<usize>,
    keywords: Vec<Keyword>,
    wcs: Option<&'a Wcs>,
    compression: Option<CompressionAlgorithm>,
    tile_size: Option<(usize, usize)>,
}

impl<'a, T> AstroImage<'a, T> {
    pub fn new(data: &'a [T], dimensions: impl Into<Vec<usize>>) -> Self {
        Self {
            data,
            dimensions: dimensions.into(),
            keywords: Vec::new(),
            wcs: None,
            compression: None,
            tile_size: None,
        }
    }

    pub fn wcs(mut self, wcs: &'a Wcs) -> Self {
        self.wcs = Some(wcs);
        self
    }

    // Replaces every card of the same name with one card, except commentary
    // keywords, which are appended. A keyword replaces the WCS card of its
    // name.
    pub fn keyword(mut self, kw: Keyword) -> Self {
        set_keyword(&mut self.keywords, kw);
        self
    }

    pub fn keywords(self, keywords: impl IntoIterator<Item = Keyword>) -> Self {
        keywords.into_iter().fold(self, Self::keyword)
    }

    pub fn compressed(mut self, algorithm: CompressionAlgorithm) -> Self {
        self.compression = Some(algorithm);
        self
    }

    pub fn tile_size(mut self, width: usize, height: usize) -> Self {
        self.tile_size = Some((width, height));
        self
    }

    pub fn exposure(self, seconds: f64) -> Self {
        self.keyword(metadata::exposure(seconds))
    }

    pub fn temperature(self, celsius: f64) -> Self {
        self.keyword(metadata::temperature(celsius))
    }

    pub fn gain(self, gain: f64) -> Self {
        self.keyword(metadata::gain(gain))
    }

    pub fn binning(self, x: u32, y: u32) -> Self {
        self.keywords(metadata::binning(x, y))
    }

    pub fn filter(self, name: &str) -> Self {
        self.keyword(metadata::filter(name))
    }

    pub fn object(self, name: &str) -> Self {
        self.keyword(metadata::object(name))
    }

    pub fn observer(self, name: &str) -> Self {
        self.keyword(metadata::observer(name))
    }

    pub fn telescope(self, name: &str) -> Self {
        self.keyword(metadata::telescope(name))
    }

    pub fn instrument(self, name: &str) -> Self {
        self.keyword(metadata::instrument(name))
    }

    pub fn date_obs(self, iso8601: &str) -> Self {
        self.keyword(metadata::date_obs(iso8601))
    }

    pub fn image_kind(&self) -> ImageKind {
        ImageKind::from_dimensions(&self.dimensions)
    }

    pub fn write_fits<P: AsRef<Path>>(&self, path: P) -> crate::fits::errors::Result<()>
    where
        T: DataArray + Clone,
    {
        let all_keywords = self.build_keywords();
        let bytes = match self.compression {
            Some(algorithm) => compressed_image_bytes(
                self.data,
                &self.dimensions,
                self.tile_size,
                algorithm,
                &all_keywords,
            ),
            None => primary_image_bytes(self.data, &self.dimensions, &all_keywords, true),
        }?;
        Ok(write_atomically(path.as_ref(), &bytes)?)
    }

    pub fn write_xisf<P: AsRef<Path>>(&self, path: P) -> crate::xisf::errors::Result<()>
    where
        T: XisfDataType,
    {
        let (width, height, channels) =
            plane_shape_with_axes(&self.dimensions, XisfError::InvalidGeometry)?;
        let mut writer = XisfWriter::new(Cursor::new(Vec::new()));
        writer.add_image(self.data, width, height, channels)?;
        write_xisf_file(writer, &self.build_keywords(), &[], path.as_ref())
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

    fn build_keywords(&self) -> Vec<Keyword> {
        let mut keywords = Vec::new();
        if let Some(wcs) = self.wcs {
            keywords.extend(wcs.to_keywords().into_iter().map(Keyword::from));
        }
        for kw in &self.keywords {
            set_keyword(&mut keywords, kw.clone());
        }
        keywords
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::shape::plane_shape;
    use crate::fits::errors::FitsError;
    use crate::fits::header::keywords::KeywordValue;
    use crate::fits::io::reader::FitsFile;
    use celestial_wcs::builder::WcsBuilder;
    use celestial_wcs::spherical::Projection;

    #[test]
    fn xisf_write_of_an_image_with_no_axes_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.xisf");
        let image = AstroImage::new(&[7u8][..], Vec::<usize>::new());
        let err = image.write_xisf(&path).unwrap_err();
        assert!(
            matches!(&err, XisfError::InvalidGeometry(m) if m == "the image has no axes"),
            "{err}"
        );
        let err = image.write_to(&path).unwrap_err();
        assert!(
            matches!(&err, ImageError::Xisf(XisfError::InvalidGeometry(m))
                if m == "the image has no axes"),
            "{err}"
        );
        assert!(!path.exists());
    }

    #[test]
    fn fits_write_of_an_image_with_no_axes_keeps_the_naxis_0_rules() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.fits");
        AstroImage::new(&[][..] as &[u8], Vec::<usize>::new())
            .write_to(&path)
            .unwrap();
        let mut fits = FitsFile::open(&path).unwrap();
        assert_eq!(
            fits.get_header_value(0, "NAXIS").unwrap().as_deref(),
            Some("0")
        );

        let path = dir.path().join("one.fits");
        let err = AstroImage::new(&[7u8][..], Vec::<usize>::new())
            .write_fits(&path)
            .unwrap_err();
        assert!(
            matches!(&err, FitsError::InvalidFormat(m)
                if m == "data has 1 values but NAXIS = 0 declares no data"),
            "{err}"
        );
        assert!(!path.exists());
    }

    fn tmp_fits() -> tempfile::NamedTempFile {
        tempfile::Builder::new().suffix(".fits").tempfile().unwrap()
    }

    fn tmp_xisf() -> tempfile::NamedTempFile {
        tempfile::Builder::new().suffix(".xisf").tempfile().unwrap()
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
        assert_eq!(img.compression, None);
        assert!(img.tile_size.is_none());
    }

    #[test]
    fn image_kind_reflects_dimensions() {
        let data = vec![0u8; 16];
        assert_eq!(AstroImage::new(&data, [4, 4]).image_kind(), ImageKind::Mono);

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
        let img = AstroImage::new(&data, [4, 4]).keyword(Keyword::string("OBJECT", "M31"));
        assert_eq!(img.keywords.len(), 1);
        assert_eq!(img.keywords[0].name, "OBJECT");
    }

    fn written_cards(img: &AstroImage<'_, u8>, name: &str) -> Vec<Option<KeywordValue>> {
        let tmp = tmp_fits();
        img.write_fits(tmp.path()).unwrap();
        let header = crate::fits::io::reader::FitsFile::open(tmp.path())
            .unwrap()
            .get_header(0)
            .unwrap();
        header
            .keywords()
            .iter()
            .filter(|k| k.name == name)
            .map(|k| k.value.clone())
            .collect()
    }

    #[test]
    fn repeated_setters_replace_the_card() {
        let data = vec![0u8; 16];
        let img = AstroImage::new(&data, [4, 4])
            .exposure(1.0)
            .object("A")
            .exposure(2.0)
            .keyword(Keyword::string("OBJECT", "B"));
        assert_eq!(
            written_cards(&img, "EXPTIME"),
            [Some(KeywordValue::Real(2.0))]
        );
        assert_eq!(
            written_cards(&img, "OBJECT"),
            [Some(KeywordValue::String("B".into()))]
        );
    }

    #[test]
    fn keyword_replaces_the_wcs_card_of_the_same_name() {
        let data = vec![0u8; 16];
        let wcs = sample_wcs();
        let img = AstroImage::new(&data, [4, 4])
            .keyword(Keyword::real("CRPIX1", 5.0))
            .wcs(&wcs);
        assert_eq!(
            written_cards(&img, "CRPIX1"),
            [Some(KeywordValue::Real(5.0))]
        );
    }

    #[test]
    fn commentary_keywords_are_appended() {
        let data = vec![0u8; 16];
        let img = AstroImage::new(&data, [4, 4])
            .keyword(Keyword::history("a"))
            .keywords([Keyword::history("b"), Keyword::comment("c")]);
        let texts: Vec<_> = img.keywords.iter().map(|k| k.comment.as_deref()).collect();
        assert_eq!(texts, [Some("a"), Some("b"), Some("c")]);
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
            .compressed(CompressionAlgorithm::Gzip2)
            .tile_size(8, 8);
        assert_eq!(img.compression, Some(CompressionAlgorithm::Gzip2));
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
        assert_eq!(
            find("EXPTIME").unwrap().value,
            Some(KeywordValue::Real(30.0))
        );
        assert_eq!(
            find("CCD-TEMP").unwrap().value,
            Some(KeywordValue::Real(-10.0))
        );
        assert_eq!(find("GAIN").unwrap().value, Some(KeywordValue::Real(1.5)));
        assert_eq!(
            find("XBINNING").unwrap().value,
            Some(KeywordValue::Integer(2))
        );
        assert_eq!(
            find("YBINNING").unwrap().value,
            Some(KeywordValue::Integer(2))
        );
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
            .compressed(CompressionAlgorithm::Rice)
            .tile_size(5, 7)
            .write_fits(tmp.path())
            .unwrap();

        let restored = crate::formats::image::Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![16usize, 16]);
        assert_eq!(restored.pixels.as_u8(), Some(&data));
    }

    #[test]
    fn write_fits_uncompressed_produces_readable_file() {
        let tmp = tmp_fits();
        let data: Vec<u8> = (0u8..=255).collect();
        AstroImage::new(&data, [16, 16])
            .write_fits(tmp.path())
            .unwrap();

        let restored = crate::formats::image::Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![16usize, 16]);
    }

    #[test]
    fn write_fits_embeds_keywords() {
        let tmp = tmp_fits();
        let data = vec![0u8; 16];
        AstroImage::new(&data, [4, 4])
            .object("M31")
            .exposure(60.5)
            .write_fits(tmp.path())
            .unwrap();

        let restored = crate::formats::image::Image::open(tmp.path()).unwrap();
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

        let restored = crate::formats::image::Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![16usize, 16]);
    }

    #[test]
    fn write_xisf_writes_the_astrometric_solution_properties() {
        let tmp = tmp_xisf();
        let data = vec![0u8; 16];
        let wcs = sample_wcs();
        let img = AstroImage::new(&data, [4, 4]).wcs(&wcs);
        img.write_xisf(tmp.path()).unwrap();

        let read = crate::xisf::reader::XisfFile::open(tmp.path())
            .unwrap()
            .properties(0)
            .unwrap();
        let expected = crate::xisf::writer::wcs::solution_properties(&img.build_keywords());
        assert_eq!(expected.len(), 6);
        assert_eq!(read, expected);
    }

    #[test]
    fn write_xisf_drops_structural_keywords() {
        let tmp = tmp_xisf();
        let data = vec![0u8; 16];
        AstroImage::new(&data, [4, 4])
            .keyword(Keyword::real("BZERO", 0.0))
            .object("M31")
            .write_xisf(tmp.path())
            .unwrap();

        let file = crate::xisf::reader::XisfFile::open(tmp.path()).unwrap();
        assert_eq!(file.keywords(), [Keyword::string("OBJECT", "M31")]);
    }

    #[test]
    fn write_to_dispatches_on_extension() {
        let data: Vec<u8> = (0u8..16).collect();

        let tmp = tmp_fits();
        AstroImage::new(&data, [4, 4]).write_to(tmp.path()).unwrap();
        assert!(tmp.path().exists());

        let tmp = tmp_xisf();
        AstroImage::new(&data, [4, 4]).write_to(tmp.path()).unwrap();
        assert!(tmp.path().exists());
    }

    #[test]
    fn write_to_rejects_unknown_extension() {
        let tmp = tempfile::Builder::new()
            .suffix(".bogus")
            .tempfile()
            .unwrap();
        let data = vec![0u8; 16];
        assert!(AstroImage::new(&data, [4, 4]).write_to(tmp.path()).is_err());
    }

    #[test]
    fn wcs_keywords_are_included_in_fits_output() {
        let tmp = tmp_fits();
        let data = vec![0u8; 16];
        let wcs = sample_wcs();

        AstroImage::new(&data, [4, 4])
            .wcs(&wcs)
            .write_fits(tmp.path())
            .unwrap();

        let restored = crate::formats::image::Image::open(tmp.path()).unwrap();
        assert!(restored.get_keyword("CRVAL1").is_some());
        assert!(restored.get_keyword("CRVAL2").is_some());
        assert!(restored.get_keyword("CRPIX1").is_some());
        assert!(restored.get_keyword("CRPIX2").is_some());
    }

    #[test]
    fn extract_dimensions_handles_defaults() {
        let data = vec![0u8; 16];
        let img = AstroImage::new(&data, [4, 4]);
        assert_eq!(plane_shape(&img.dimensions), (4, 4, 1));

        let rgb = AstroImage::new(&data, vec![4usize, 4, 3]);
        assert_eq!(plane_shape(&rgb.dimensions), (4, 4, 3));

        let empty_dims: Vec<usize> = Vec::new();
        let img_empty = AstroImage::new(&data, empty_dims);
        assert_eq!(plane_shape(&img_empty.dimensions), (0, 1, 1));
    }
}
