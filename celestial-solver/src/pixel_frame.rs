//! Pixel coordinate conventions, and the conversions between them.
//!
//! The solver works in detection pixels: the centre of the first pixel is at (0, 0), x
//! runs along a row and y runs down the rows. Detection, matching, every fit and the
//! log-odds gate use them. Two outputs use other conventions, with the same axes and
//! directions:
//!
//! - FITS WCS (Calabretta and Greisen 2002, Paper I, section 2.1.4) puts the centre of
//!   the first pixel at (1, 1).
//! - PixInsight puts it at (0.5, 0.5). The FITS CRPIX it writes is its own reference
//!   pixel plus 0.5. On five frames solved by both PixInsight and astrometry.net, the
//!   astrometry.net star positions are our centroids plus 1 px, and its WCS gives the
//!   sky position closest to PixInsight's at its FITS pixel minus 0.5 on four of the
//!   five. On one more frame, PixInsight's spline control points are our centroids plus
//!   about 0.5 px in x. Their y offset changes by 0.35 px with the centroid window, so
//!   y on that frame does not settle the offset.
//!
//! Every conversion between these conventions is in this module.

use celestial_catalog::query::healpix::reduce_ra_deg;
use celestial_wcs::builder::Wcs;
use celestial_wcs::coordinate::{CelestialCoord, PixelCoord};
use celestial_wcs::error::WcsError;

use crate::fit_wcs::solution::WcsSolution;
use crate::solve::SolveResult;

const FITS_OFFSET_PX: f64 = 1.0;
const PIXINSIGHT_OFFSET_PX: f64 = 0.5;

pub(crate) fn detection_to_fits((x, y): (f64, f64)) -> (f64, f64) {
    (x + FITS_OFFSET_PX, y + FITS_OFFSET_PX)
}

pub(crate) fn detection_to_pixinsight((x, y): (f64, f64)) -> (f64, f64) {
    (x + PIXINSIGHT_OFFSET_PX, y + PIXINSIGHT_OFFSET_PX)
}

pub(crate) fn pixinsight_to_detection((x, y): (f64, f64)) -> (f64, f64) {
    (x - PIXINSIGHT_OFFSET_PX, y - PIXINSIGHT_OFFSET_PX)
}

// The sky position, (ra, dec) in degrees with ra in [0, 360), of detection pixel (x, y)
// under `wcs`, a FITS WCS. An error when x or y is not finite or `wcs` does not map the
// pixel to the sky.
pub(crate) fn detection_to_sky(wcs: &Wcs, x: f64, y: f64) -> Result<(f64, f64), WcsError> {
    require_finite(x, y, "pixel")?;
    let (x, y) = detection_to_fits((x, y));
    let sky = wcs.pixel_to_celestial(PixelCoord::new(x, y))?;
    // celestial_wcs returns right ascension in (-180, 180]; the solver reports [0, 360).
    Ok((reduce_ra_deg(sky.alpha_deg()), sky.delta_deg()))
}

// The detection pixel of sky position (ra, dec) in degrees under `wcs`, a FITS WCS. An
// error when ra or dec is not finite or `wcs` does not map the position to a pixel.
pub(crate) fn sky_to_detection(wcs: &Wcs, ra: f64, dec: f64) -> Result<(f64, f64), WcsError> {
    require_finite(ra, dec, "sky position")?;
    let pixel = wcs.celestial_to_pixel(CelestialCoord::new(ra, dec))?;
    Ok((pixel.x() - FITS_OFFSET_PX, pixel.y() - FITS_OFFSET_PX))
}

// PixInsight's linear map from its pixels to the tangent plane, world = CD (pixel -
// reference), as the 2x3 matrix [CD | -CD reference] in row order.
pub(crate) fn pixinsight_linear_transformation(wcs: &WcsSolution) -> [f64; 6] {
    let (rx, ry) = detection_to_pixinsight((wcs.crpix1, wcs.crpix2));
    [
        wcs.cd1_1,
        wcs.cd1_2,
        -(wcs.cd1_1 * rx + wcs.cd1_2 * ry),
        wcs.cd2_1,
        wcs.cd2_2,
        -(wcs.cd2_1 * rx + wcs.cd2_2 * ry),
    ]
}

fn require_finite(a: f64, b: f64, what: &str) -> Result<(), WcsError> {
    if a.is_finite() && b.is_finite() {
        return Ok(());
    }
    Err(WcsError::InvalidParameter {
        message: format!("{what} ({a}, {b}) is not finite"),
    })
}

impl SolveResult {
    /// Transforms a pixel coordinate to sky coordinates (RA, Dec in degrees).
    ///
    /// The pixel is in detection coordinates, with the centre of the first pixel at
    /// (0, 0). Applies the SIP distortion when the solution has one, then the CD matrix
    /// and the TAN de-projection.
    ///
    /// # Errors
    ///
    /// [`WcsError`] when `x` or `y` is not finite or the pixel has no sky position.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// let (ra, dec) = result.pixel_to_sky(100.0, 200.0)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn pixel_to_sky(&self, x: f64, y: f64) -> Result<(f64, f64), WcsError> {
        detection_to_sky(&self.fits_wcs, x, y)
    }

    /// Transforms sky coordinates (RA, Dec in degrees) to pixel coordinates.
    ///
    /// The pixel is in detection coordinates, with the centre of the first pixel at
    /// (0, 0). Applies the SIP distortion when the solution has one.
    ///
    /// # Errors
    ///
    /// [`WcsError`] when `ra` or `dec` is not finite or the position has no pixel, for
    /// example more than 90° from CRVAL.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// let (x, y) = result.sky_to_pixel(180.0, 30.0)?;
    /// println!("(180°, 30°) lands at ({:.1}, {:.1})", x, y);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn sky_to_pixel(&self, ra: f64, dec: f64) -> Result<(f64, f64), WcsError> {
        sky_to_detection(&self.fits_wcs, ra, dec)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::fit_wcs::fits_wcs::build_fits_wcs;
    use crate::fit_wcs::sip::SipSolution;
    use crate::fit_wcs::solution::WcsSolution;

    fn tan() -> WcsSolution {
        WcsSolution {
            crpix1: 1023.0,
            crpix2: 767.0,
            ..WcsSolution::diagonal((180.0, 30.0), (-2.8e-4, 2.8e-4), (2048, 1536))
        }
    }

    // A_2_0 = 2^-14, so a pixel 100 px right of CRPIX moves right by exactly
    // 100^2 * 2^-14 = 0.6103515625 px.
    fn sip(wcs: &WcsSolution) -> SipSolution {
        SipSolution {
            crpix1: wcs.crpix1,
            crpix2: wcs.crpix2,
            crval1: wcs.crval1,
            crval2: wcs.crval2,
            cd1_1: wcs.cd1_1,
            cd1_2: wcs.cd1_2,
            cd2_1: wcs.cd2_1,
            cd2_2: wcs.cd2_2,
            width: wcs.width,
            height: wcs.height,
            a_order: 2,
            b_order: 2,
            a_coeffs: BTreeMap::from([((2, 0), libm::pow(2.0, -14.0))]),
            b_coeffs: BTreeMap::new(),
            n_stars: 0,
            linear_rms_px: 0.0,
            rms_px: 0.0,
            weighted_rms_px: 0.0,
        }
    }

    #[test]
    fn the_reference_pixel_maps_to_the_reference_sky_position() {
        let wcs = tan();
        let fits = build_fits_wcs(&wcs, None).unwrap();
        assert_eq!(fits.crpix(), [1024.0, 768.0]);
        let sky = detection_to_sky(&fits, wcs.crpix1, wcs.crpix2).unwrap();
        assert_eq!(sky, (wcs.crval1, wcs.crval2));
    }

    // The reference pixel of a field at RA 315 degrees, and a pixel 1000 px to its right,
    // where celestial_wcs returns RA as a negative angle.
    #[test]
    fn sky_positions_wrap_into_0_to_360_degrees() {
        let wcs = WcsSolution {
            crval1: 315.0,
            ..tan()
        };
        let fits = build_fits_wcs(&wcs, None).unwrap();
        assert_eq!(
            detection_to_sky(&fits, wcs.crpix1, wcs.crpix2).unwrap(),
            (315.0, 30.0)
        );
        let (ra, dec) = detection_to_sky(&fits, wcs.crpix1 + 1000.0, wcs.crpix2).unwrap();
        let unwrapped = fits
            .pixel_to_celestial(PixelCoord::new(2024.0, 768.0))
            .unwrap();
        assert_eq!(
            (ra, dec),
            (unwrapped.alpha_deg() + 360.0, unwrapped.delta_deg())
        );
        assert_eq!(ra, wcs.pixel_to_sky(wcs.crpix1 + 1000.0, wcs.crpix2).0);
    }

    #[test]
    fn sip_moves_a_pixel_by_its_distortion() {
        let wcs = tan();
        let with_sip = build_fits_wcs(&wcs, Some(&sip(&wcs))).unwrap();
        let without = build_fits_wcs(&wcs, None).unwrap();
        let (x, y) = (wcs.crpix1 + 100.0, wcs.crpix2);
        assert_eq!(
            detection_to_sky(&with_sip, x, y).unwrap(),
            detection_to_sky(&without, x + 0.6103515625, y).unwrap()
        );
    }

    #[test]
    fn non_finite_input_is_an_error() {
        let fits = build_fits_wcs(&tan(), None).unwrap();
        let invalid =
            |r: Result<(f64, f64), WcsError>| matches!(r, Err(WcsError::InvalidParameter { .. }));
        assert!(invalid(detection_to_sky(&fits, f64::NAN, 0.0)));
        assert!(invalid(detection_to_sky(&fits, 0.0, f64::INFINITY)));
        assert!(invalid(sky_to_detection(&fits, f64::NAN, 30.0)));
        assert!(invalid(sky_to_detection(&fits, 180.0, f64::NEG_INFINITY)));
    }
}
