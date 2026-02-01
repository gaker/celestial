use crate::fits::compression::CompressionAlgorithm;
use crate::fits::data::array::DataArray;
use crate::fits::header::Keyword;
use crate::fits::io::writer::FitsWriter;
use crate::fits::Result;
use celestial_wcs::{Wcs, WcsKeyword, WcsKeywordValue};
use std::path::Path;

const DEFAULT_TILE_SIZE: usize = 32;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImageKind {
    Mono,
    Rgb,
    Cube,
}

impl ImageKind {
    pub fn from_dimensions(dims: &[usize]) -> Self {
        match dims.len() {
            1 | 2 => Self::Mono,
            3 if dims[2] == 3 => Self::Rgb,
            _ => Self::Cube,
        }
    }
}

pub struct FitsImage<'a, T> {
    data: &'a [T],
    dimensions: Vec<usize>,
    keywords: Vec<Keyword>,
    wcs: Option<&'a Wcs>,
    compressed: bool,
    tile_size: Option<(usize, usize)>,
}

impl<'a, T: DataArray + Clone> FitsImage<'a, T> {
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

    pub fn keyword(mut self, keyword: Keyword) -> Self {
        self.keywords.push(keyword);
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

    pub fn image_kind(&self) -> ImageKind {
        ImageKind::from_dimensions(&self.dimensions)
    }

    pub fn write_to<P: AsRef<Path>>(&self, path: P) -> Result<()> {
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

    fn build_keywords(&self) -> Vec<Keyword> {
        let mut keywords = Vec::new();

        if let Some(wcs) = self.wcs {
            keywords.extend(wcs_to_keywords(wcs));
        }

        keywords.extend(self.keywords.clone());
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

fn wcs_to_keywords(wcs: &Wcs) -> Vec<Keyword> {
    wcs.to_keywords()
        .into_iter()
        .map(wcs_keyword_to_fits)
        .collect()
}

fn wcs_keyword_to_fits(wk: WcsKeyword) -> Keyword {
    match wk.value {
        WcsKeywordValue::Real(v) => Keyword::real(wk.name, v),
        WcsKeywordValue::Integer(v) => Keyword::integer(wk.name, v),
        WcsKeywordValue::String(v) => Keyword::string(wk.name, v),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_wcs::{Projection, WcsBuilder};
    use tempfile::NamedTempFile;

    #[test]
    fn image_kind_from_2d_dimensions() {
        assert_eq!(ImageKind::from_dimensions(&[100, 100]), ImageKind::Mono);
        assert_eq!(ImageKind::from_dimensions(&[512, 512]), ImageKind::Mono);
    }

    #[test]
    fn image_kind_from_rgb_dimensions() {
        assert_eq!(ImageKind::from_dimensions(&[100, 100, 3]), ImageKind::Rgb);
    }

    #[test]
    fn image_kind_from_cube_dimensions() {
        assert_eq!(ImageKind::from_dimensions(&[100, 100, 10]), ImageKind::Cube);
        assert_eq!(ImageKind::from_dimensions(&[64, 64, 64]), ImageKind::Cube);
    }

    #[test]
    fn image_kind_from_1d_dimensions() {
        assert_eq!(ImageKind::from_dimensions(&[100]), ImageKind::Mono);
    }

    #[test]
    fn fits_image_builder_defaults() {
        let data: Vec<i16> = vec![1, 2, 3, 4];
        let image = FitsImage::new(&data, [2, 2]);

        assert_eq!(image.image_kind(), ImageKind::Mono);
        assert!(image.compressed);
        assert!(image.wcs.is_none());
        assert!(image.keywords.is_empty());
    }

    #[test]
    fn fits_image_builder_with_keywords() {
        let data: Vec<i16> = vec![1, 2, 3, 4];
        let image = FitsImage::new(&data, [2, 2])
            .keyword(Keyword::string("OBJECT", "M31"))
            .keyword(Keyword::real("EXPTIME", 30.0));

        assert_eq!(image.keywords.len(), 2);
    }

    #[test]
    fn fits_image_builder_uncompressed() {
        let data: Vec<i16> = vec![1, 2, 3, 4];
        let image = FitsImage::new(&data, [2, 2]).compressed(false);

        assert!(!image.compressed);
    }

    #[test]
    fn fits_image_builder_custom_tile_size() {
        let data: Vec<i16> = (0..256).collect();
        let image = FitsImage::new(&data, [16, 16]).tile_size(8, 8);

        assert_eq!(image.tile_size, Some((8, 8)));
    }

    #[test]
    fn fits_image_compute_tile_size_default() {
        let data: Vec<i16> = (0..10000).collect();
        let image = FitsImage::new(&data, [100, 100]);

        let tile_size = image.compute_tile_size();
        assert_eq!(tile_size, (32, 32));
    }

    #[test]
    fn fits_image_compute_tile_size_small_image() {
        let data: Vec<i16> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let image = FitsImage::new(&data, [5, 2]);

        let tile_size = image.compute_tile_size();
        assert_eq!(tile_size, (5, 2));
    }

    #[test]
    fn fits_image_compute_tile_size_custom() {
        let data: Vec<i16> = (0..10000).collect();
        let image = FitsImage::new(&data, [100, 100]).tile_size(64, 64);

        let tile_size = image.compute_tile_size();
        assert_eq!(tile_size, (64, 64));
    }

    #[test]
    fn fits_image_write_compressed() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let data: Vec<i16> = (0..100).collect();
        let result = FitsImage::new(&data, [10, 10]).write_to(path);

        assert!(result.is_ok());
    }

    #[test]
    fn fits_image_write_uncompressed() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let data: Vec<i16> = (0..100).collect();
        let result = FitsImage::new(&data, [10, 10])
            .compressed(false)
            .write_to(path);

        assert!(result.is_ok());
    }

    #[test]
    fn fits_image_write_with_keywords() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let data: Vec<i16> = (0..100).collect();
        let result = FitsImage::new(&data, [10, 10])
            .keyword(Keyword::string("OBJECT", "Test"))
            .keyword(Keyword::real("EXPTIME", 60.0))
            .compressed(false)
            .write_to(path);

        assert!(result.is_ok());
    }

    #[test]
    fn fits_image_write_with_wcs() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let wcs = WcsBuilder::new()
            .crpix(5.0, 5.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .build()
            .unwrap();

        let data: Vec<i16> = (0..100).collect();
        let result = FitsImage::new(&data, [10, 10])
            .wcs(&wcs)
            .compressed(false)
            .write_to(path);

        assert!(result.is_ok());
    }

    #[test]
    fn wcs_keyword_to_fits_real() {
        let wk = WcsKeyword::real("CRPIX1", 512.0);
        let kw = wcs_keyword_to_fits(wk);

        assert_eq!(kw.name, "CRPIX1");
        assert!(kw.value.is_some());
    }

    #[test]
    fn wcs_keyword_to_fits_integer() {
        let wk = WcsKeyword::integer("NAXIS", 2);
        let kw = wcs_keyword_to_fits(wk);

        assert_eq!(kw.name, "NAXIS");
        assert!(kw.value.is_some());
    }

    #[test]
    fn wcs_keyword_to_fits_string() {
        let wk = WcsKeyword::string("CTYPE1", "RA---TAN");
        let kw = wcs_keyword_to_fits(wk);

        assert_eq!(kw.name, "CTYPE1");
        assert!(kw.value.is_some());
    }

    #[test]
    fn fits_image_write_with_wcs_roundtrip() {
        use crate::fits::FitsFile;

        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.path();

        let wcs = WcsBuilder::new()
            .crpix(5.0, 5.0)
            .crval(180.0, 45.0)
            .cd_matrix([[0.001, 0.0], [0.0, 0.001]])
            .projection(Projection::tan())
            .build()
            .unwrap();

        let data: Vec<i16> = (0..100).collect();
        FitsImage::new(&data, [10, 10])
            .wcs(&wcs)
            .compressed(false)
            .write_to(path)
            .unwrap();

        let mut fits = FitsFile::open(path).unwrap();
        let header = fits.get_header(0).unwrap();

        assert!(header.get_keyword_value("CTYPE1").is_some());
        assert!(header.get_keyword_value("CRPIX1").is_some());
        assert!(header.get_keyword_value("CRVAL1").is_some());
        assert!(header.get_keyword_value("CD1_1").is_some());
    }

    #[test]
    fn fits_image_multiple_keywords_method() {
        let data: Vec<i16> = vec![1, 2, 3, 4];
        let kws = vec![
            Keyword::string("OBJECT", "M31"),
            Keyword::real("EXPTIME", 30.0),
            Keyword::integer("GAIN", 100),
        ];

        let image = FitsImage::new(&data, [2, 2]).keywords(kws);

        assert_eq!(image.keywords.len(), 3);
    }
}
