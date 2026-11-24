use celestial_coords::frames::icrs::ICRSPosition;
use celestial_images::formats::image::Image;
use celestial_time::julian::JulianDate;

use super::error::MetadataError;
use super::hint::hint_from_image;
use super::plate_scale::{plate_scale_from_image, positive};
use super::time::jd_from_image;

/// Everything the solver needs to start a solve: position hint, plate scale, epoch.
///
/// Build this from an [`Image`] via [`metadata_from_image`], or construct it by hand
/// for images without headers and pass it to [`crate::solve::Solver::with_metadata`].
///
/// `focal_mm` and `pixel_um` are optional — they're only used by the WCS fitter as
/// sanity-check priors. `scale_arcsec` is the authoritative plate scale.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::metadata::image_metadata::ImageMetadata;
/// use celestial_coords::frames::icrs::ICRSPosition;
/// use celestial_time::scales::utc::utc_from_calendar;
///
/// let meta = ImageMetadata {
///     hint: ICRSPosition::from_degrees(83.633, 22.014)?,
///     scale_arcsec: 1.5,
///     epoch: utc_from_calendar(2026, 4, 21, 20, 0, 0.0)?.to_julian_date(),
///     focal_mm: None,
///     pixel_um: None,
/// };
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct ImageMetadata {
    /// Approximate center of the field. The catalog cone search is centered here.
    pub hint: ICRSPosition,
    /// Plate scale in arcseconds per pixel.
    pub scale_arcsec: f64,
    /// Observation epoch as a Julian Date in UTC. Used for catalog proper-motion
    /// propagation.
    pub epoch: JulianDate,
    /// Telescope focal length in mm. Optional.
    pub focal_mm: Option<f64>,
    /// Pixel size in microns. Optional.
    pub pixel_um: Option<f64>,
}

impl ImageMetadata {
    // An error when the plate scale, or a focal length or pixel size that is set, is not
    // finite and positive, or the epoch is not finite. The hint is an `ICRSPosition`,
    // which is finite by construction.
    pub(crate) fn validate(&self) -> Result<(), MetadataError> {
        positive("plate scale", self.scale_arcsec)?;
        let jd = self.epoch.jd1() + self.epoch.jd2();
        if !jd.is_finite() {
            return Err(MetadataError::NonFiniteValue {
                field: "epoch",
                value: jd,
            });
        }
        if let Some(focal) = self.focal_mm {
            positive("focal length", focal)?;
        }
        if let Some(pixel) = self.pixel_um {
            positive("pixel size", pixel)?;
        }
        Ok(())
    }
}

/// Selective overrides for header-driven metadata.
///
/// Each `Some` field bypasses its corresponding header lookup. Each `None` field
/// falls through to the header. Constructed directly when using
/// [`metadata_from_image`]; set via builder methods ([`crate::solve::Solver::focal_length_mm`]
/// etc.) when using the [`crate::solve::Solver`].
///
/// `scale_arcsec` takes priority over `focal_mm` / `pixel_um` — set it and the
/// `206.2648 * pixel_um / focal_mm` calculation is skipped.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::metadata::image_metadata::{metadata_from_image, MetadataOverrides};
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// let overrides = MetadataOverrides {
///     focal_mm: Some(500.0),       // override header
///     ..Default::default()         // everything else from headers
/// };
/// let meta = metadata_from_image(&img, &overrides)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Default, Clone)]
pub struct MetadataOverrides {
    /// Position hint override. Bypasses `CRVAL1/2`, `RA/DEC`, `OBJCTRA/DEC`.
    pub hint: Option<ICRSPosition>,
    /// Epoch override. Bypasses `DATE-OBS` and `MJD-OBS`.
    pub epoch: Option<JulianDate>,
    /// Plate scale override (arcseconds per pixel). Takes priority over `focal_mm` and
    /// `pixel_um`.
    pub scale_arcsec: Option<f64>,
    /// Focal length override (mm). Bypasses `FOCALLEN`.
    pub focal_mm: Option<f64>,
    /// Pixel size override (microns). Bypasses `XPIXSZ`.
    pub pixel_um: Option<f64>,
}

/// Builds an [`ImageMetadata`] from an image's headers, with optional per-field overrides.
///
/// For each field, an `overrides.<field>` of `Some(value)` takes priority; `None` falls
/// back to reading the corresponding FITS header.
///
/// # Errors
///
/// Returns [`MetadataError`] when the position hint, the epoch or the plate scale is
/// missing from both the overrides and the headers, when a plate scale, focal length or
/// pixel size is not finite and positive, or when the epoch is not finite.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::metadata::image_metadata::{metadata_from_image, MetadataOverrides};
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// let meta = metadata_from_image(&img, &MetadataOverrides::default())?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn metadata_from_image(
    img: &Image,
    overrides: &MetadataOverrides,
) -> Result<ImageMetadata, MetadataError> {
    let hint = match &overrides.hint {
        Some(h) => h.clone(),
        None => hint_from_image(img)?,
    };
    let epoch = match overrides.epoch {
        Some(e) => e,
        None => jd_from_image(img)?,
    };
    let (scale_arcsec, focal_mm, pixel_um) = plate_scale(img, overrides)?;
    let meta = ImageMetadata {
        hint,
        scale_arcsec,
        epoch,
        focal_mm,
        pixel_um,
    };
    meta.validate()?;
    Ok(meta)
}

// The plate scale in arcseconds per pixel, the focal length and the pixel size: the
// override scale with the override optics, or the scale from the optics and headers.
fn plate_scale(
    img: &Image,
    overrides: &MetadataOverrides,
) -> Result<(f64, Option<f64>, Option<f64>), MetadataError> {
    if let Some(direct) = overrides.scale_arcsec {
        return Ok((direct, overrides.focal_mm, overrides.pixel_um));
    }
    let scale = plate_scale_from_image(img, overrides.focal_mm, overrides.pixel_um)?;
    Ok((scale.arcsec_per_pixel, scale.focal_mm, scale.pixel_um))
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_images::fits::header::keywords::Keyword;
    use celestial_images::formats::pixel_data::PixelData;

    fn full_image() -> Image {
        let mut img = Image::new(PixelData::F32(vec![0.0; 4]), [2, 2]);
        img.set_keyword(Keyword::real("CRVAL1", 180.0));
        img.set_keyword(Keyword::real("CRVAL2", 20.0));
        img.set_keyword(Keyword::string("DATE-OBS", "2024-06-15T03:22:45"));
        img.set_keyword(Keyword::real("FOCALLEN", 1000.0));
        img.set_keyword(Keyword::real("XPIXSZ", 5.0));
        img
    }

    #[test]
    fn populated_image_produces_all_fields() {
        let img = full_image();
        let meta = metadata_from_image(&img, &MetadataOverrides::default()).unwrap();
        assert_eq!(meta.hint.ra().degrees(), 180.0);
        assert_eq!(meta.hint.dec().degrees(), 20.0);
        assert_eq!(meta.focal_mm, Some(1000.0));
        assert_eq!(meta.pixel_um, Some(5.0));
        assert_eq!(
            meta.scale_arcsec,
            crate::metadata::plate_scale::ARCSEC_PER_UM_PER_MM * 5.0 / 1000.0
        );
    }

    #[test]
    fn overrides_take_priority_over_headers() {
        let img = full_image();
        let overrides = MetadataOverrides {
            focal_mm: Some(500.0),
            pixel_um: Some(3.76),
            ..Default::default()
        };
        let meta = metadata_from_image(&img, &overrides).unwrap();
        assert_eq!(meta.focal_mm, Some(500.0));
        assert_eq!(meta.pixel_um, Some(3.76));
    }

    #[test]
    fn missing_hint_propagates_error() {
        let mut img = full_image();
        for name in ["CRVAL1", "CRVAL2", "RA", "DEC", "OBJCTRA", "OBJCTDEC"] {
            img.remove_keyword(name);
        }
        let err = match metadata_from_image(&img, &MetadataOverrides::default()) {
            Err(e) => e,
            Ok(_) => panic!("expected Err"),
        };
        assert!(matches!(err, MetadataError::NoPositionHint));
    }

    #[test]
    fn missing_plate_scale_propagates_error() {
        let mut img = full_image();
        img.remove_keyword("FOCALLEN");
        let err = match metadata_from_image(&img, &MetadataOverrides::default()) {
            Err(e) => e,
            Ok(_) => panic!("expected Err"),
        };
        assert!(matches!(err, MetadataError::MissingHeader("FOCALLEN")));
    }

    #[test]
    fn overrides_fill_in_for_image_without_plate_scale_headers() {
        let mut img = full_image();
        img.remove_keyword("FOCALLEN");
        img.remove_keyword("XPIXSZ");
        let overrides = MetadataOverrides {
            focal_mm: Some(800.0),
            pixel_um: Some(2.4),
            ..Default::default()
        };
        let meta = metadata_from_image(&img, &overrides).unwrap();
        assert_eq!(meta.focal_mm, Some(800.0));
    }

    #[test]
    fn non_finite_or_non_positive_overrides_are_rejected() {
        let img = full_image();
        let scale = |scale_arcsec: f64| MetadataOverrides {
            scale_arcsec: Some(scale_arcsec),
            ..Default::default()
        };
        let cases = [
            (scale(f64::NAN), "plate scale must be finite, got NaN"),
            (scale(-1.5), "plate scale must be positive, got -1.5"),
            (
                MetadataOverrides {
                    focal_mm: Some(f64::INFINITY),
                    ..scale(1.5)
                },
                "focal length must be finite, got inf",
            ),
            (
                MetadataOverrides {
                    epoch: Some(JulianDate::new(f64::NAN, 0.0)),
                    ..Default::default()
                },
                "epoch must be finite, got NaN",
            ),
        ];
        for (overrides, expected) in cases {
            let err = metadata_from_image(&img, &overrides).err();
            assert_eq!(err.map(|e| e.to_string()).as_deref(), Some(expected));
        }
    }
}
