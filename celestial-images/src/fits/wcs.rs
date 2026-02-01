use crate::fits::header::Header;
use crate::fits::{FitsError, Result};
use celestial_wcs::{KeywordProvider, Wcs, WcsBuilder};

struct FitsKeywordAdapter<'a> {
    header: &'a Header,
}

impl KeywordProvider for FitsKeywordAdapter<'_> {
    fn get_string(&self, key: &str) -> Option<String> {
        self.header
            .get_keyword_value(key)?
            .as_string()
            .map(|s| s.to_string())
    }

    fn get_float(&self, key: &str) -> Option<f64> {
        self.header.get_keyword_value(key)?.as_real()
    }

    fn get_int(&self, key: &str) -> Option<i64> {
        self.header.get_keyword_value(key)?.as_integer()
    }
}

pub struct WcsInfo {
    wcs: Wcs,
}

impl WcsInfo {
    pub fn from_header(header: &Header) -> Result<Option<Self>> {
        if header.get_keyword_value("CTYPE1").is_none() {
            return Ok(None);
        }

        let adapter = FitsKeywordAdapter { header };
        let wcs = WcsBuilder::from_header(&adapter)
            .map_err(|e| FitsError::InvalidFormat(format!("WCS: {}", e)))?
            .build()
            .map_err(|e| FitsError::InvalidFormat(format!("WCS: {}", e)))?;

        Ok(Some(Self { wcs }))
    }

    pub fn pix2world(&self, x: f64, y: f64) -> Result<(f64, f64)> {
        self.wcs
            .pix2world(x, y)
            .map_err(|e| FitsError::InvalidFormat(format!("WCS transform: {}", e)))
    }

    pub fn world2pix(&self, ra: f64, dec: f64) -> Result<(f64, f64)> {
        self.wcs
            .world2pix(ra, dec)
            .map_err(|e| FitsError::InvalidFormat(format!("WCS transform: {}", e)))
    }

    #[inline]
    pub fn projection_code(&self) -> &str {
        self.wcs.projection_code()
    }

    #[inline]
    pub fn pixel_scale(&self) -> f64 {
        self.wcs.pixel_scale()
    }

    #[inline]
    pub fn crpix(&self) -> [f64; 2] {
        self.wcs.crpix()
    }

    #[inline]
    pub fn crval(&self) -> (f64, f64) {
        self.wcs.crval()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Keyword, KeywordValue};

    fn create_tan_header() -> Header {
        let mut header = Header::new();
        header.add_keyword(
            Keyword::new("CTYPE1".to_string())
                .with_value(KeywordValue::String("RA---TAN".to_string())),
        );
        header.add_keyword(
            Keyword::new("CTYPE2".to_string())
                .with_value(KeywordValue::String("DEC--TAN".to_string())),
        );
        header
            .add_keyword(Keyword::new("CRPIX1".to_string()).with_value(KeywordValue::Real(512.0)));
        header
            .add_keyword(Keyword::new("CRPIX2".to_string()).with_value(KeywordValue::Real(512.0)));
        header
            .add_keyword(Keyword::new("CRVAL1".to_string()).with_value(KeywordValue::Real(180.0)));
        header.add_keyword(Keyword::new("CRVAL2".to_string()).with_value(KeywordValue::Real(45.0)));
        header.add_keyword(Keyword::new("CD1_1".to_string()).with_value(KeywordValue::Real(-1e-4)));
        header.add_keyword(Keyword::new("CD1_2".to_string()).with_value(KeywordValue::Real(0.0)));
        header.add_keyword(Keyword::new("CD2_1".to_string()).with_value(KeywordValue::Real(0.0)));
        header.add_keyword(Keyword::new("CD2_2".to_string()).with_value(KeywordValue::Real(1e-4)));
        header
    }

    #[test]
    fn test_wcs_from_header() {
        let header = create_tan_header();
        let wcs = WcsInfo::from_header(&header).unwrap().unwrap();

        assert_eq!(wcs.projection_code(), "TAN");
        assert_eq!(wcs.crpix(), [512.0, 512.0]);
        assert_eq!(wcs.crval(), (180.0, 45.0));
    }

    #[test]
    fn test_wcs_missing_ctype_returns_none() {
        let header = Header::new();
        let result = WcsInfo::from_header(&header).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_wcs_roundtrip() {
        let header = create_tan_header();
        let wcs = WcsInfo::from_header(&header).unwrap().unwrap();

        let (ra, dec) = wcs.pix2world(512.0, 512.0).unwrap();
        assert!((ra - 180.0).abs() < 1e-10);
        assert!((dec - 45.0).abs() < 1e-10);

        let (x, y) = wcs.world2pix(ra, dec).unwrap();
        assert!((x - 512.0).abs() < 1e-10);
        assert!((y - 512.0).abs() < 1e-10);
    }

    #[test]
    fn test_wcs_off_center_roundtrip() {
        let header = create_tan_header();
        let wcs = WcsInfo::from_header(&header).unwrap().unwrap();

        let (ra, dec) = wcs.pix2world(300.0, 700.0).unwrap();
        let (x, y) = wcs.world2pix(ra, dec).unwrap();

        assert!((x - 300.0).abs() < 1e-7);
        assert!((y - 700.0).abs() < 1e-7);
    }
}
