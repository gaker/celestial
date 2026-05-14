use celestial_coords::ICRSPosition;
use celestial_core::angle::{parse_dms, parse_hms};
use celestial_images::fits::header::Keyword;
use celestial_images::formats::Image;

use super::error::MetadataError;

/// Reads the position hint from an image's headers.
///
/// Tries three header pairs in order, taking the first one that's present:
///
/// 1. `CRVAL1` / `CRVAL2` — decimal degrees (WCS projection origin)
/// 2. `RA` / `DEC` — either decimal degrees or sexagesimal strings
/// 3. `OBJCTRA` / `OBJCTDEC` — sexagesimal strings (PixInsight, NINA)
///
/// # Errors
///
/// Returns [`MetadataError::NoPositionHint`] if none of the pairs are present,
/// [`MetadataError::InvalidHeaderType`] if a header exists with the wrong value kind,
/// or [`MetadataError::AngleParse`] for malformed sexagesimal strings.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::metadata::hint_from_image;
///
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// let pos = hint_from_image(&img)?;
/// println!("{:.4}° {:.4}°", pos.ra().degrees(), pos.dec().degrees());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn hint_from_image(img: &Image) -> Result<ICRSPosition, MetadataError> {
    if let (Some(cr1), Some(cr2)) = (img.get_keyword("CRVAL1"), img.get_keyword("CRVAL2")) {
        let ra = cr1
            .value
            .as_ref()
            .and_then(|v| v.as_real())
            .ok_or(MetadataError::InvalidHeaderType {
                keyword: "CRVAL1",
                kind: "numeric",
            })?;
        let dec = cr2
            .value
            .as_ref()
            .and_then(|v| v.as_real())
            .ok_or(MetadataError::InvalidHeaderType {
                keyword: "CRVAL2",
                kind: "numeric",
            })?;
        return Ok(ICRSPosition::from_degrees(ra, dec)?);
    }

    if let (Some(ra_kw), Some(dec_kw)) = (img.get_keyword("RA"), img.get_keyword("DEC")) {
        let (ra, dec) = parse_ra_dec_keywords(ra_kw, dec_kw)?;
        return Ok(ICRSPosition::from_degrees(ra, dec)?);
    }

    if let (Some(ra_kw), Some(dec_kw)) =
        (img.get_keyword("OBJCTRA"), img.get_keyword("OBJCTDEC"))
    {
        let ra_str = ra_kw
            .value
            .as_ref()
            .and_then(|v| v.as_string())
            .ok_or(MetadataError::InvalidHeaderType {
                keyword: "OBJCTRA",
                kind: "string",
            })?;
        let dec_str = dec_kw
            .value
            .as_ref()
            .and_then(|v| v.as_string())
            .ok_or(MetadataError::InvalidHeaderType {
                keyword: "OBJCTDEC",
                kind: "string",
            })?;
        let ra = parse_hms(ra_str)
            .map_err(|e| MetadataError::AngleParse {
                value: ra_str.to_string(),
                source: e,
            })?
            .degrees();
        let dec = parse_dms(dec_str)
            .map_err(|e| MetadataError::AngleParse {
                value: dec_str.to_string(),
                source: e,
            })?
            .degrees();
        return Ok(ICRSPosition::from_degrees(ra, dec)?);
    }

    Err(MetadataError::NoPositionHint)
}

fn parse_ra_dec_keywords(
    ra_kw: &Keyword,
    dec_kw: &Keyword,
) -> Result<(f64, f64), MetadataError> {
    let ra_val = ra_kw
        .value
        .as_ref()
        .ok_or(MetadataError::InvalidHeaderType {
            keyword: "RA",
            kind: "any",
        })?;
    let dec_val = dec_kw
        .value
        .as_ref()
        .ok_or(MetadataError::InvalidHeaderType {
            keyword: "DEC",
            kind: "any",
        })?;

    let ra = match ra_val.as_real() {
        Some(v) => v,
        None => {
            let s = ra_val.as_string().ok_or(MetadataError::InvalidHeaderType {
                keyword: "RA",
                kind: "numeric or string",
            })?;
            parse_hms(s)
                .map_err(|e| MetadataError::AngleParse {
                    value: s.to_string(),
                    source: e,
                })?
                .degrees()
        }
    };
    let dec = match dec_val.as_real() {
        Some(v) => v,
        None => {
            let s = dec_val.as_string().ok_or(MetadataError::InvalidHeaderType {
                keyword: "DEC",
                kind: "numeric or string",
            })?;
            parse_dms(s)
                .map_err(|e| MetadataError::AngleParse {
                    value: s.to_string(),
                    source: e,
                })?
                .degrees()
        }
    };
    Ok((ra, dec))
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_images::formats::PixelData;

    fn stub_image() -> Image {
        Image::new(PixelData::F32(vec![0.0; 4]), [2, 2])
    }

    #[test]
    fn crval_pair_is_the_highest_priority_source() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("CRVAL1", 100.0));
        img.set_keyword(Keyword::real("CRVAL2", 20.0));
        img.set_keyword(Keyword::real("RA", 999.0));
        img.set_keyword(Keyword::real("DEC", 999.0));
        img.set_keyword(Keyword::string("OBJCTRA", "00 00 00"));
        img.set_keyword(Keyword::string("OBJCTDEC", "+00 00 00"));
        let pos = hint_from_image(&img).unwrap();
        assert_eq!(pos.ra().degrees(), 100.0);
        assert_eq!(pos.dec().degrees(), 20.0);
    }

    #[test]
    fn crval1_non_numeric_returns_invalid_header_type() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("CRVAL1", "bogus"));
        img.set_keyword(Keyword::real("CRVAL2", 0.0));
        let err = hint_from_image(&img).unwrap_err();
        assert!(
            matches!(err, MetadataError::InvalidHeaderType { keyword: "CRVAL1", .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn crval2_non_numeric_returns_invalid_header_type() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("CRVAL1", 0.0));
        img.set_keyword(Keyword::string("CRVAL2", "bogus"));
        let err = hint_from_image(&img).unwrap_err();
        assert!(
            matches!(err, MetadataError::InvalidHeaderType { keyword: "CRVAL2", .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn ra_dec_numeric_is_treated_as_degrees() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("RA", 83.633));
        img.set_keyword(Keyword::real("DEC", 22.0145));
        let pos = hint_from_image(&img).unwrap();
        assert_eq!(pos.ra().degrees(), 83.633);
        assert_eq!(pos.dec().degrees(), 22.0145);
    }

    #[test]
    fn ra_dec_space_separated_sexagesimal_roundtrips() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("RA", "05 34 31.94"));
        img.set_keyword(Keyword::string("DEC", "+22 00 52.2"));
        let pos = hint_from_image(&img).unwrap();
        let expected_ra = (5.0 + 34.0 / 60.0 + 31.94 / 3600.0) * 15.0;
        let expected_dec = 22.0 + 52.2 / 3600.0;
        assert!((pos.ra().degrees() - expected_ra).abs() < 1e-10);
        assert!((pos.dec().degrees() - expected_dec).abs() < 1e-10);
    }

    #[test]
    fn ra_dec_colon_separated_sexagesimal_handles_negative_dec() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("RA", "05:34:31.94"));
        img.set_keyword(Keyword::string("DEC", "-22:00:52.2"));
        let pos = hint_from_image(&img).unwrap();
        let expected_ra = (5.0 + 34.0 / 60.0 + 31.94 / 3600.0) * 15.0;
        let expected_dec = -(22.0 + 52.2 / 3600.0);
        assert!((pos.ra().degrees() - expected_ra).abs() < 1e-10);
        assert!((pos.dec().degrees() - expected_dec).abs() < 1e-10);
    }

    #[test]
    fn ra_dec_malformed_string_returns_angle_parse() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("RA", "this is not an hms"));
        img.set_keyword(Keyword::real("DEC", 0.0));
        let err = hint_from_image(&img).unwrap_err();
        assert!(matches!(err, MetadataError::AngleParse { .. }), "{err}");
    }

    #[test]
    fn ra_numeric_dec_string_parses_both_paths() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("RA", 83.633));
        img.set_keyword(Keyword::string("DEC", "+22:00:52.2"));
        let pos = hint_from_image(&img).unwrap();
        assert!((pos.ra().degrees() - 83.633).abs() < 1e-10);
        assert!((pos.dec().degrees() - (22.0 + 52.2 / 3600.0)).abs() < 1e-10);
    }

    #[test]
    fn objctra_objctdec_requires_string_values() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("OBJCTRA", "12 30 00.00"));
        img.set_keyword(Keyword::string("OBJCTDEC", "-30 15 30.00"));
        let pos = hint_from_image(&img).unwrap();
        let expected_ra = (12.0 + 30.0 / 60.0) * 15.0;
        let expected_dec = -(30.0 + 15.0 / 60.0 + 30.0 / 3600.0);
        assert!((pos.ra().degrees() - expected_ra).abs() < 1e-10);
        assert!((pos.dec().degrees() - expected_dec).abs() < 1e-10);
    }

    #[test]
    fn objctra_numeric_is_rejected_as_invalid_header_type() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("OBJCTRA", 12.5));
        img.set_keyword(Keyword::string("OBJCTDEC", "+00 00 00"));
        let err = hint_from_image(&img).unwrap_err();
        assert!(
            matches!(err, MetadataError::InvalidHeaderType { keyword: "OBJCTRA", .. }),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn objctra_malformed_returns_angle_parse() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("OBJCTRA", "nonsense"));
        img.set_keyword(Keyword::string("OBJCTDEC", "+00 00 00"));
        let err = hint_from_image(&img).unwrap_err();
        assert!(
            matches!(&err, MetadataError::AngleParse { value, .. } if value == "nonsense"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn objctdec_malformed_returns_angle_parse() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("OBJCTRA", "12 30 00"));
        img.set_keyword(Keyword::string("OBJCTDEC", "bogus"));
        let err = hint_from_image(&img).unwrap_err();
        assert!(
            matches!(&err, MetadataError::AngleParse { value, .. } if value == "bogus"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn no_position_keywords_returns_no_position_hint() {
        let img = stub_image();
        let err = hint_from_image(&img).unwrap_err();
        assert!(matches!(err, MetadataError::NoPositionHint));
    }

    #[test]
    fn dec_malformed_string_with_numeric_ra_returns_angle_parse() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("RA", 83.633));
        img.set_keyword(Keyword::string("DEC", "this is not a dms"));
        let err = hint_from_image(&img).unwrap_err();
        assert!(
            matches!(&err, MetadataError::AngleParse { value, .. } if value == "this is not a dms"),
            "unexpected error: {err}"
        );
    }
}
