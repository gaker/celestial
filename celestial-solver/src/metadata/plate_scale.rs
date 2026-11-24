use celestial_core::constants::ARCSEC_PER_RAD;
use celestial_images::formats::image::Image;

use super::error::MetadataError;

// The angle in arcseconds of 1 micron seen from 1 mm: a pixel of `pixel_um` behind a
// focal length of `focal_mm` subtends `pixel_um / (1000 * focal_mm)` radians.
pub(crate) const ARCSEC_PER_UM_PER_MM: f64 = ARCSEC_PER_RAD / 1000.0;

/// Plate scale derived from focal length and pixel size.
///
/// `arcsec_per_pixel` follows the small-angle approximation:
/// `206.2648 * pixel_um / focal_mm`.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::metadata::plate_scale::plate_scale_from_image;
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// let scale = plate_scale_from_image(&img, None, None)?;
/// println!("{:.3} arcsec/px", scale.arcsec_per_pixel);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct PlateScale {
    /// Scale in arcseconds per pixel.
    pub arcsec_per_pixel: f64,
    /// Focal length in mm (from override or `FOCALLEN`).
    pub focal_mm: Option<f64>,
    /// Pixel size in microns (from override or `XPIXSZ`).
    pub pixel_um: Option<f64>,
}

/// Reads the plate scale from an image's headers with optional overrides.
///
/// Reads `FOCALLEN` (mm) and `XPIXSZ` (microns) from the image headers, preferring the
/// `focal_override` / `pixel_override` arguments when they're `Some`. Computes the
/// scale as `206.2648 * pixel_um / focal_mm`.
///
/// # Errors
///
/// Returns [`MetadataError::MissingHeader`] if `FOCALLEN` or `XPIXSZ` is missing and
/// no override supplies it, [`MetadataError::NonPositiveValue`] if either is zero
/// or negative, or [`MetadataError::NonFiniteValue`] if either is NaN or infinite.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::metadata::plate_scale::plate_scale_from_image;
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// // Use FOCALLEN and XPIXSZ from headers:
/// let scale = plate_scale_from_image(&img, None, None)?;
///
/// // Override focal length (e.g. scope has a reducer inline):
/// let scale = plate_scale_from_image(&img, Some(500.0), None)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn plate_scale_from_image(
    img: &Image,
    focal_override: Option<f64>,
    pixel_override: Option<f64>,
) -> Result<PlateScale, MetadataError> {
    let focal = focal_override.or_else(|| img.get_last("FOCALLEN").as_f64());
    let pixel = pixel_override.or_else(|| img.get_last("XPIXSZ").as_f64());

    let focal_mm = focal.ok_or(MetadataError::MissingHeader("FOCALLEN"))?;
    let pixel_um = pixel.ok_or(MetadataError::MissingHeader("XPIXSZ"))?;

    positive("focal length", focal_mm)?;
    positive("pixel size", pixel_um)?;

    Ok(PlateScale {
        arcsec_per_pixel: ARCSEC_PER_UM_PER_MM * pixel_um / focal_mm,
        focal_mm: Some(focal_mm),
        pixel_um: Some(pixel_um),
    })
}

// `value` when it is finite and above 0.
pub(super) fn positive(field: &'static str, value: f64) -> Result<f64, MetadataError> {
    if !value.is_finite() {
        return Err(MetadataError::NonFiniteValue { field, value });
    }
    if value <= 0.0 {
        return Err(MetadataError::NonPositiveValue { field, value });
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_images::fits::header::keywords::Keyword;
    use celestial_images::formats::pixel_data::PixelData;

    fn stub_image() -> Image {
        Image::new(PixelData::F32(vec![0.0; 4]), [2, 2])
    }

    #[test]
    fn headers_only_compute_expected_plate_scale() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("FOCALLEN", 1000.0));
        img.set_keyword(Keyword::real("XPIXSZ", 5.0));
        let scale = plate_scale_from_image(&img, None, None).unwrap();
        assert_eq!(scale.focal_mm, Some(1000.0));
        assert_eq!(scale.pixel_um, Some(5.0));
        assert_eq!(scale.arcsec_per_pixel, ARCSEC_PER_UM_PER_MM * 5.0 / 1000.0);
    }

    #[test]
    fn overrides_take_priority_over_headers() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("FOCALLEN", 1000.0));
        img.set_keyword(Keyword::real("XPIXSZ", 5.0));
        let scale = plate_scale_from_image(&img, Some(500.0), Some(3.76)).unwrap();
        assert_eq!(scale.focal_mm, Some(500.0));
        assert_eq!(scale.pixel_um, Some(3.76));
        assert_eq!(scale.arcsec_per_pixel, ARCSEC_PER_UM_PER_MM * 3.76 / 500.0);
    }

    #[test]
    fn overrides_fill_in_when_headers_absent() {
        let img = stub_image();
        let scale = plate_scale_from_image(&img, Some(800.0), Some(2.4)).unwrap();
        assert_eq!(scale.arcsec_per_pixel, ARCSEC_PER_UM_PER_MM * 2.4 / 800.0);
    }

    fn expect_err(result: Result<PlateScale, MetadataError>) -> MetadataError {
        match result {
            Ok(_) => panic!("expected Err, got Ok"),
            Err(e) => e,
        }
    }

    #[test]
    fn missing_focal_length_returns_missing_header() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("XPIXSZ", 5.0));
        let err = expect_err(plate_scale_from_image(&img, None, None));
        assert!(matches!(err, MetadataError::MissingHeader("FOCALLEN")));
    }

    #[test]
    fn missing_pixel_size_returns_missing_header() {
        let mut img = stub_image();
        img.set_keyword(Keyword::real("FOCALLEN", 1000.0));
        let err = expect_err(plate_scale_from_image(&img, None, None));
        assert!(matches!(err, MetadataError::MissingHeader("XPIXSZ")));
    }

    #[test]
    fn zero_or_negative_focal_length_is_rejected() {
        let img = stub_image();
        for bad in [0.0, -1.0] {
            let err = expect_err(plate_scale_from_image(&img, Some(bad), Some(5.0)));
            assert!(
                matches!(err, MetadataError::NonPositiveValue { field: "focal length", value } if value == bad),
                "unexpected error for focal={bad}: {err}"
            );
        }
    }

    #[test]
    fn zero_or_negative_pixel_size_is_rejected() {
        let img = stub_image();
        for bad in [0.0, -0.5] {
            let err = expect_err(plate_scale_from_image(&img, Some(1000.0), Some(bad)));
            assert!(
                matches!(err, MetadataError::NonPositiveValue { field: "pixel size", value } if value == bad),
                "unexpected error for pixel={bad}: {err}"
            );
        }
    }

    #[test]
    fn non_real_header_is_ignored_and_falls_through_to_missing() {
        let mut img = stub_image();
        img.set_keyword(Keyword::string("FOCALLEN", "bogus"));
        img.set_keyword(Keyword::real("XPIXSZ", 5.0));
        let err = expect_err(plate_scale_from_image(&img, None, None));
        assert!(matches!(err, MetadataError::MissingHeader("FOCALLEN")));
    }

    #[test]
    fn non_finite_focal_length_and_pixel_size_are_rejected() {
        let img = stub_image();
        let cases = [
            (f64::NAN, 3.76, "focal length must be finite, got NaN"),
            (f64::INFINITY, 3.76, "focal length must be finite, got inf"),
            (500.0, f64::NAN, "pixel size must be finite, got NaN"),
            (
                500.0,
                f64::NEG_INFINITY,
                "pixel size must be finite, got -inf",
            ),
        ];
        for (focal, pixel, expected) in cases {
            let err = plate_scale_from_image(&img, Some(focal), Some(pixel)).err();
            assert_eq!(err.map(|e| e.to_string()).as_deref(), Some(expected));
        }
    }

    #[test]
    fn repeated_focallen_reads_the_last_card() {
        let mut img = stub_image();
        img.keywords.push(Keyword::real("FOCALLEN", 500.0));
        img.keywords.push(Keyword::real("FOCALLEN", 1000.0));
        img.keywords.push(Keyword::real("XPIXSZ", 5.0));
        let scale = plate_scale_from_image(&img, None, None).unwrap();
        assert_eq!(scale.focal_mm, Some(1000.0));
    }
}
