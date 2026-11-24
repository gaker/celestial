use celestial_coords::frames::icrs::ICRSPosition;
use celestial_core::angle::{parse_dms, parse_hms, Angle};
use celestial_core::errors::AstroError;
use celestial_images::fits::header::parser::Field;
use celestial_images::formats::image::Image;

use super::error::MetadataError;
use super::header_wcs::{header_wcs, HeaderWcs};

/// Reads the position hint from an image's headers.
///
/// Takes the first of these that the header has:
///
/// 1. The sky position of the image centre under the header's celestial WCS: `CTYPE1`
///    and `CTYPE2` naming right ascension and declination, `CRPIX1/2`, `CRVAL1/2`, and
///    a CD matrix or `CDELT1/2`
/// 2. `CRVAL1` / `CRVAL2`, decimal degrees, when there is no `CTYPE` card or the WCS
///    does not build; skipped when `CTYPE` names other axes
/// 3. `RA` / `DEC`, either decimal degrees or sexagesimal strings
/// 4. `OBJCTRA` / `OBJCTDEC`, sexagesimal strings (PixInsight, NINA)
///
/// A key that repeats is read from its last card.
///
/// # Errors
///
/// Returns [`MetadataError::NoPositionHint`] if none of the sources are present,
/// [`MetadataError::InvalidHeaderType`] if a header exists with the wrong value kind,
/// or [`MetadataError::AngleParse`] for malformed sexagesimal strings.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::metadata::hint::hint_from_image;
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// let pos = hint_from_image(&img)?;
/// println!("{:.4}° {:.4}°", pos.ra().degrees(), pos.dec().degrees());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn hint_from_image(img: &Image) -> Result<ICRSPosition, MetadataError> {
    if let Some(pos) = header_hint(img)? {
        return Ok(pos);
    }

    let (ra_kw, dec_kw) = (img.get_last("RA"), img.get_last("DEC"));
    if ra_kw.exists() && dec_kw.exists() {
        let (ra, dec) = parse_ra_dec_keywords(ra_kw, dec_kw)?;
        return Ok(ICRSPosition::from_degrees(ra, dec)?);
    }

    let (objctra, objctdec) = (img.get_last("OBJCTRA"), img.get_last("OBJCTDEC"));
    if objctra.exists() && objctdec.exists() {
        return objct_hint(objctra, objctdec);
    }

    Err(MetadataError::NoPositionHint)
}

// The image centre under the header WCS, or else CRVAL1/2 unless CTYPE names axes other
// than RA and DEC.
fn header_hint(img: &Image) -> Result<Option<ICRSPosition>, MetadataError> {
    let crval_is_sky = match header_wcs(img) {
        HeaderWcs::Centre(ra, dec) => return Ok(Some(ICRSPosition::from_degrees(ra, dec)?)),
        HeaderWcs::Absent => true,
        HeaderWcs::Unusable(reason) => {
            log::warn!("the header WCS gives no position hint: {reason}");
            true
        }
        HeaderWcs::OtherAxes => {
            log::info!("CTYPE names axes other than RA and DEC, so CRVAL is not a hint");
            false
        }
    };
    if !crval_is_sky {
        return Ok(None);
    }
    crval_hint(img)
}

fn crval_hint(img: &Image) -> Result<Option<ICRSPosition>, MetadataError> {
    let (cr1, cr2) = (img.get_last("CRVAL1"), img.get_last("CRVAL2"));
    if !(cr1.exists() && cr2.exists()) {
        return Ok(None);
    }
    let ra = cr1.as_f64().ok_or(MetadataError::InvalidHeaderType {
        keyword: "CRVAL1",
        kind: "numeric",
    })?;
    let dec = cr2.as_f64().ok_or(MetadataError::InvalidHeaderType {
        keyword: "CRVAL2",
        kind: "numeric",
    })?;
    Ok(Some(ICRSPosition::from_degrees(ra, dec)?))
}

fn objct_hint(objctra: Field<'_>, objctdec: Field<'_>) -> Result<ICRSPosition, MetadataError> {
    let ra_str = objctra.as_str().ok_or(MetadataError::InvalidHeaderType {
        keyword: "OBJCTRA",
        kind: "string",
    })?;
    let dec_str = objctdec.as_str().ok_or(MetadataError::InvalidHeaderType {
        keyword: "OBJCTDEC",
        kind: "string",
    })?;
    let ra = sexagesimal_deg(ra_str, parse_hms)?;
    let dec = sexagesimal_deg(dec_str, parse_dms)?;
    Ok(ICRSPosition::from_degrees(ra, dec)?)
}

fn parse_ra_dec_keywords(ra_kw: Field<'_>, dec_kw: Field<'_>) -> Result<(f64, f64), MetadataError> {
    let ra = degrees_or_sexagesimal(ra_kw, "RA", parse_hms)?;
    let dec = degrees_or_sexagesimal(dec_kw, "DEC", parse_dms)?;
    Ok((ra, dec))
}

type AngleParser = fn(&str) -> Result<Angle, AstroError>;

// The value of `keyword` in decimal degrees, or its sexagesimal string read by `parse`.
fn degrees_or_sexagesimal(
    field: Field<'_>,
    keyword: &'static str,
    parse: AngleParser,
) -> Result<f64, MetadataError> {
    if let Some(v) = field.as_f64() {
        return Ok(v);
    }
    let s = field.as_str().ok_or(MetadataError::InvalidHeaderType {
        keyword,
        kind: "numeric or string",
    })?;
    sexagesimal_deg(s, parse)
}

fn sexagesimal_deg(s: &str, parse: AngleParser) -> Result<f64, MetadataError> {
    let angle = parse(s).map_err(|e| MetadataError::AngleParse {
        value: s.to_string(),
        source: e,
    })?;
    Ok(angle.degrees())
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_images::fits::header::keywords::Keyword;
    use celestial_images::formats::pixel_data::PixelData;
    use celestial_wcs::builder::WcsBuilder;
    use celestial_wcs::coordinate::PixelCoord;

    fn stub_image() -> Image {
        Image::new(PixelData::F32(vec![0.0; 4]), [2, 2])
    }

    // A 100 x 80 image with a TAN WCS whose reference pixel is `crpix`, in FITS pixels.
    fn image_with_tan_wcs(crpix: (f64, f64)) -> Image {
        let mut img = Image::new(PixelData::F32(vec![0.0; 100 * 80]), [100, 80]);
        img.set_keyword(Keyword::string("CTYPE1", "RA---TAN"));
        img.set_keyword(Keyword::string("CTYPE2", "DEC--TAN"));
        img.set_keyword(Keyword::real("CRPIX1", crpix.0));
        img.set_keyword(Keyword::real("CRPIX2", crpix.1));
        img.set_keyword(Keyword::real("CRVAL1", 150.0));
        img.set_keyword(Keyword::real("CRVAL2", 30.0));
        img.set_keyword(Keyword::real("CD1_1", -0.001));
        img.set_keyword(Keyword::real("CD1_2", 0.0));
        img.set_keyword(Keyword::real("CD2_1", 0.0));
        img.set_keyword(Keyword::real("CD2_2", 0.001));
        img
    }

    // (ra, dec) in degrees as an ICRSPosition holds them.
    fn icrs_degrees(ra: f64, dec: f64) -> (f64, f64) {
        let pos = ICRSPosition::from_degrees(ra, dec).unwrap();
        (pos.ra().degrees(), pos.dec().degrees())
    }

    // The sky position of FITS pixel `pixel` under the WCS of `image_with_tan_wcs`.
    fn tan_wcs_sky(crpix: (f64, f64), pixel: (f64, f64)) -> (f64, f64) {
        let wcs = WcsBuilder::new()
            .crpix(crpix.0, crpix.1)
            .crval(150.0, 30.0)
            .cd_matrix([[-0.001, 0.0], [0.0, 0.001]])
            .proj_code("TAN")
            .build()
            .unwrap();
        let sky = wcs
            .pixel_to_celestial(PixelCoord::new(pixel.0, pixel.1))
            .unwrap();
        icrs_degrees(sky.alpha_deg(), sky.delta_deg())
    }

    #[test]
    fn repeated_crval_reads_the_last_card() {
        let mut img = stub_image();
        img.keywords.push(Keyword::real("CRVAL1", 10.0));
        img.keywords.push(Keyword::real("CRVAL1", 11.0));
        img.keywords.push(Keyword::real("CRVAL2", 20.0));
        let pos = hint_from_image(&img).unwrap();
        assert_eq!(
            (pos.ra().degrees(), pos.dec().degrees()),
            icrs_degrees(11.0, 20.0)
        );
    }

    // FITS puts the centre of the first pixel at (1, 1), so the centre of a 100 x 80
    // image is at (50.5, 40.5). CRPIX is 450.5 px and 259.5 px from it, so the hint is
    // about 0.5 deg from CRVAL.
    #[test]
    fn header_wcs_gives_the_sky_position_of_the_image_centre() {
        let crpix = (-400.0, 300.0);
        let pos = hint_from_image(&image_with_tan_wcs(crpix)).unwrap();
        assert_eq!(
            (pos.ra().degrees(), pos.dec().degrees()),
            tan_wcs_sky(crpix, (50.5, 40.5))
        );
    }

    #[test]
    fn header_wcs_reads_the_last_card_of_a_repeated_key() {
        let crpix = (-400.0, 300.0);
        let mut img = image_with_tan_wcs((0.0, 0.0));
        img.keywords.push(Keyword::real("CRPIX1", crpix.0));
        img.keywords.push(Keyword::real("CRPIX2", crpix.1));
        let pos = hint_from_image(&img).unwrap();
        assert_eq!(
            (pos.ra().degrees(), pos.dec().degrees()),
            tan_wcs_sky(crpix, (50.5, 40.5))
        );
    }

    // CRVAL is a position on the sky only when the axes are right ascension and
    // declination.
    #[test]
    fn crval_of_other_axis_types_is_not_a_position_hint() {
        for ctype in [["GLON-TAN", "GLAT-TAN"], ["PIXEL", "PIXEL"]] {
            let mut img = image_with_tan_wcs((50.5, 40.5));
            img.set_keyword(Keyword::string("CTYPE1", ctype[0]));
            img.set_keyword(Keyword::string("CTYPE2", ctype[1]));
            img.set_keyword(Keyword::real("RA", 180.0));
            img.set_keyword(Keyword::real("DEC", 10.0));
            let pos = hint_from_image(&img).unwrap();
            assert_eq!(
                (pos.ra().degrees(), pos.dec().degrees()),
                icrs_degrees(180.0, 10.0),
                "{ctype:?}"
            );
        }
    }

    // A header WCS that does not build leaves CRVAL as the hint: here one without a CD
    // matrix or CDELT, and one with a projection code that does not exist.
    #[test]
    fn header_wcs_that_does_not_build_falls_back_to_crval() {
        let mut no_matrix = image_with_tan_wcs((-400.0, 300.0));
        for key in ["CD1_1", "CD1_2", "CD2_1", "CD2_2"] {
            no_matrix.remove_keyword(key);
        }
        let mut unknown_projection = image_with_tan_wcs((-400.0, 300.0));
        unknown_projection.set_keyword(Keyword::string("CTYPE1", "RA---XYZ"));
        unknown_projection.set_keyword(Keyword::string("CTYPE2", "DEC--XYZ"));
        for img in [no_matrix, unknown_projection] {
            let pos = hint_from_image(&img).unwrap();
            assert_eq!(
                (pos.ra().degrees(), pos.dec().degrees()),
                icrs_degrees(150.0, 30.0)
            );
        }
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
            matches!(
                err,
                MetadataError::InvalidHeaderType {
                    keyword: "CRVAL1",
                    ..
                }
            ),
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
            matches!(
                err,
                MetadataError::InvalidHeaderType {
                    keyword: "CRVAL2",
                    ..
                }
            ),
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
        assert_eq!(
            (pos.ra().degrees(), pos.dec().degrees()),
            (expected_ra, expected_dec)
        );
    }

    #[test]
    fn ra_dec_colon_separated_sexagesimal_handles_negative_dec() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("RA", "05:34:31.94"));
        img.set_keyword(Keyword::string("DEC", "-22:00:52.2"));
        let pos = hint_from_image(&img).unwrap();
        let expected_ra = (5.0 + 34.0 / 60.0 + 31.94 / 3600.0) * 15.0;
        let expected_dec = -(22.0 + 52.2 / 3600.0);
        assert_eq!(
            (pos.ra().degrees(), pos.dec().degrees()),
            (expected_ra, expected_dec)
        );
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
        assert_eq!(
            (pos.ra().degrees(), pos.dec().degrees()),
            (83.633, 22.0 + 52.2 / 3600.0)
        );
    }

    #[test]
    fn objctra_objctdec_requires_string_values() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("OBJCTRA", "12 30 00.00"));
        img.set_keyword(Keyword::string("OBJCTDEC", "-30 15 30.00"));
        let pos = hint_from_image(&img).unwrap();
        let expected_ra = (12.0 + 30.0 / 60.0) * 15.0;
        // -(30 + 15 / 60 + 30 / 3600) is -30.258333333333333; the parsed Dec is 1 ulp
        // (3.6e-15 deg) from it.
        assert_eq!(
            (pos.ra().degrees(), pos.dec().degrees()),
            (expected_ra, -30.25833333333333)
        );
    }

    #[test]
    fn objctra_numeric_is_rejected_as_invalid_header_type() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("OBJCTRA", 12.5));
        img.set_keyword(Keyword::string("OBJCTDEC", "+00 00 00"));
        let err = hint_from_image(&img).unwrap_err();
        assert!(
            matches!(
                err,
                MetadataError::InvalidHeaderType {
                    keyword: "OBJCTRA",
                    ..
                }
            ),
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
