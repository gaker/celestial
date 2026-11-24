//! WCS solution data and pixel↔sky transforms.

use celestial_catalog::query::quad::tan_project_star;
use celestial_core::constants::RAD_TO_DEG;

use crate::metadata::plate_scale::ARCSEC_PER_UM_PER_MM;

/// Per-pair fit residual.
///
/// Computed for each pair a [`WcsSolution`] was fitted on. `err_x` and `err_y` are
/// `detected - predicted` in pixel space. Positive values mean the detected
/// centroid landed to the right of / below where the WCS says the catalog star should
/// be. `err_px` is the Euclidean magnitude.
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
/// let worst = result.wcs.residuals.iter()
///     .max_by(|a, b| a.err_px.partial_cmp(&b.err_px).unwrap()).unwrap();
/// println!("worst residual: {:.2}px at ({:.0}, {:.0})", worst.err_px, worst.px_x, worst.px_y);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct StarResidual {
    /// Detected centroid x in pixel coordinates.
    pub px_x: f64,
    /// Detected centroid y in pixel coordinates.
    pub px_y: f64,
    /// Residual in x (detected minus predicted), in pixels.
    pub err_x: f64,
    /// Residual in y (detected minus predicted), in pixels.
    pub err_y: f64,
    /// Euclidean residual magnitude in pixels.
    pub err_px: f64,
    // Catalog position of the pair, in degrees. With the pixel it identifies the pair,
    // since one detection can be paired with more than one catalog star.
    pub(crate) ra_deg: f64,
    pub(crate) dec_deg: f64,
}

/// Fitted TAN (gnomonic) WCS with residual statistics.
///
/// Standard FITS WCS parameters plus fit diagnostics. Pixel↔sky transforms, with the
/// SIP distortion when there is one, via [`crate::solve::SolveResult::pixel_to_sky`] and
/// [`crate::solve::SolveResult::sky_to_pixel`]. Geometric properties via
/// [`WcsSolution::scale_arcsec`], [`WcsSolution::rotation_deg`],
/// [`WcsSolution::is_flipped`].
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
/// let (ra, dec) = result.pixel_to_sky(1024.0, 768.0)?;
/// println!("center: {:.4}° {:.4}°", ra, dec);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct WcsSolution {
    /// Reference pixel x, with the centre of the first pixel at 0. The FITS CRPIX1 card
    /// is this plus 1. Usually the image center.
    pub crpix1: f64,
    /// Reference pixel y.
    pub crpix2: f64,
    /// RA at CRPIX, in degrees.
    pub crval1: f64,
    /// Dec at CRPIX, in degrees.
    pub crval2: f64,
    /// CD matrix element (1,1) in deg/pixel.
    pub cd1_1: f64,
    /// CD matrix element (1,2).
    pub cd1_2: f64,
    /// CD matrix element (2,1).
    pub cd2_1: f64,
    /// CD matrix element (2,2).
    pub cd2_2: f64,
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Focal length in mm (from metadata — not re-derived).
    pub focal_mm: Option<f64>,
    /// Pixel size in microns (from metadata — not re-derived).
    pub pixel_um: Option<f64>,
    /// Number of inlier pairs in the fit.
    pub n_stars: usize,
    /// Unweighted RMS residual in pixels.
    pub rms_px: f64,
    /// RMS residual in pixels, with each pair weighted by the inverse of its centroid
    /// variance.
    pub weighted_rms_px: f64,
    /// Per-pair residuals for the final fit.
    pub residuals: Vec<StarResidual>,
}

impl WcsSolution {
    /// Transforms a pixel coordinate to sky coordinates (RA, Dec in degrees).
    ///
    /// Applies the CD matrix, then de-projects from the tangent plane. Linear TAN
    /// only: [`crate::solve::SolveResult::pixel_to_sky`] applies SIP as well.
    #[cfg(test)]
    pub(crate) fn pixel_to_sky(&self, px_x: f64, px_y: f64) -> (f64, f64) {
        use celestial_catalog::query::quad::tan_deproject_star;
        use celestial_core::constants::DEG_TO_RAD;

        let u = px_x - self.crpix1;
        let v = px_y - self.crpix2;
        let xi_deg = self.cd1_1 * u + self.cd1_2 * v;
        let eta_deg = self.cd2_1 * u + self.cd2_2 * v;
        tan_deproject_star(
            xi_deg * DEG_TO_RAD,
            eta_deg * DEG_TO_RAD,
            self.crval1,
            self.crval2,
        )
    }

    /// Transforms sky coordinates to pixel coordinates.
    ///
    /// Returns `None` for points that can't be gnomonically projected (e.g. the
    /// anti-center of the TAN projection, more than 90° from CRVAL). Linear TAN only:
    /// [`crate::solve::SolveResult::sky_to_pixel`] applies SIP as well.
    pub(crate) fn sky_to_pixel(&self, ra_deg: f64, dec_deg: f64) -> Option<(f64, f64)> {
        let (xi_rad, eta_rad) = tan_project_star(ra_deg, dec_deg, self.crval1, self.crval2)?;
        let xi = xi_rad * RAD_TO_DEG;
        let eta = eta_rad * RAD_TO_DEG;
        let inverse = InverseCd::new(self.cd1_1, self.cd1_2, self.cd2_1, self.cd2_2);
        let (u, v) = inverse.apply(xi, eta);
        Some((u + self.crpix1, v + self.crpix2))
    }

    /// Plate scale in arcseconds per pixel along each CD-matrix axis.
    ///
    /// Returns `(scale_x, scale_y)` — the column norms of the CD matrix converted
    /// from deg/px to arcsec/px. For a well-behaved image these are nearly equal
    /// (non-square pixels or strong anisotropy widen the gap).
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// let (sx, sy) = result.wcs.scale_arcsec();
    /// println!("{:.3} × {:.3} arcsec/px", sx, sy);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn scale_arcsec(&self) -> (f64, f64) {
        let sx = libm::sqrt(self.cd1_1 * self.cd1_1 + self.cd2_1 * self.cd2_1) * 3600.0;
        let sy = libm::sqrt(self.cd1_2 * self.cd1_2 + self.cd2_2 * self.cd2_2) * 3600.0;
        (sx, sy)
    }

    /// Field rotation in degrees, in (-180, 180].
    ///
    /// The angle from image up to north, positive clockwise, with pixel y growing
    /// downwards as PixInsight displays the image. At zero, north is up, and east is
    /// left when the image is not flipped and right when [`WcsSolution::is_flipped`]
    /// returns `true`. At +90 north points to +x in either parity: the north arrow
    /// turns clockwise on the display, and the camera frame turns counter-clockwise
    /// relative to the sky.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// println!("rotation: {:.2}°", result.wcs.rotation_deg());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn rotation_deg(&self) -> f64 {
        // Pixel y grows downwards, as PixInsight displays it, so at zero rotation north
        // is towards -y. Mirroring the image reverses the sign of the east axis.
        let rotation = if self.is_flipped() {
            libm::atan2(self.cd1_2, -self.cd2_2) * RAD_TO_DEG
        } else {
            libm::atan2(-self.cd1_2, -self.cd2_2) * RAD_TO_DEG
        };
        // atan2(-0.0, x) is -pi for negative x, so a half turn would read -180 for the
        // identity CD and +180 for the same orientation flipped.
        if rotation <= -180.0 {
            rotation + 360.0
        } else {
            rotation
        }
    }

    pub(crate) fn mean_scale_arcsec(&self) -> f64 {
        let (sx, sy) = self.scale_arcsec();
        (sx + sy) * 0.5
    }

    /// `true` if the image is mirrored relative to the sky (CD matrix determinant is
    /// negative).
    ///
    /// Standard after-meridian-flip condition for many mounts, or any setup where
    /// the sensor is oriented to give a visually intuitive left-right matching the
    /// sky through a diagonal.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// if result.wcs.is_flipped() {
    ///     println!("image is flipped");
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn is_flipped(&self) -> bool {
        self.cd1_1 * self.cd2_2 - self.cd1_2 * self.cd2_1 < 0.0
    }

    /// Effective focal length implied by the fitted plate scale and known pixel size.
    ///
    /// Returns `None` if [`WcsSolution::pixel_um`] is unknown or the average scale is
    /// non-positive. Useful as a sanity check against the nominal telescope focal
    /// length — large discrepancies indicate a reducer or Barlow in the light path
    /// that wasn't reported via `FOCALLEN`.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// if let Some(f) = result.wcs.effective_focal_mm() {
    ///     println!("effective focal length: {:.0}mm", f);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn effective_focal_mm(&self) -> Option<f64> {
        let pixel_um = self.pixel_um?;
        let (sx, sy) = self.scale_arcsec();
        let avg_scale = (sx + sy) / 2.0;
        if avg_scale <= 0.0 {
            return None;
        }
        Some(ARCSEC_PER_UM_PER_MM * pixel_um / avg_scale)
    }

    pub(crate) fn nearest_residual(&self, px: f64, py: f64, k: usize) -> (f64, f64) {
        if self.residuals.is_empty() {
            return (0.0, 0.0);
        }
        let mut distances: Vec<(usize, f64)> = self
            .residuals
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let dx = r.px_x - px;
                let dy = r.px_y - py;
                (i, dx * dx + dy * dy)
            })
            .collect();
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let take = k.min(distances.len());
        let (sum_ex, sum_ey) = distances[..take]
            .iter()
            .fold((0.0, 0.0), |(sx, sy), (i, _)| {
                let r = &self.residuals[*i];
                (sx + r.err_x, sy + r.err_y)
            });
        (sum_ex / take as f64, sum_ey / take as f64)
    }
}

// The inverse of a CD matrix: maps a tangent-plane offset (xi, eta) in degrees to a
// pixel offset from CRPIX. Every pixel the solver predicts from a sky position goes
// through it, so the fit, its residuals and the matching round the same way.
#[derive(Debug, Clone, Copy)]
pub(super) struct InverseCd {
    cd: [f64; 4],
    det: f64,
}

impl InverseCd {
    pub(super) fn new(cd1_1: f64, cd1_2: f64, cd2_1: f64, cd2_2: f64) -> Self {
        Self {
            cd: [cd1_1, cd1_2, cd2_1, cd2_2],
            det: cd1_1 * cd2_2 - cd1_2 * cd2_1,
        }
    }

    pub(super) fn apply(&self, xi: f64, eta: f64) -> (f64, f64) {
        let [cd1_1, cd1_2, cd2_1, cd2_2] = self.cd;
        (
            (cd2_2 * xi - cd1_2 * eta) / self.det,
            (-cd2_1 * xi + cd1_1 * eta) / self.det,
        )
    }
}

pub(crate) fn fmt_fov_component(arcmin: f64) -> String {
    if arcmin >= 60.0 {
        let deg = (arcmin / 60.0) as i32;
        let rem = arcmin - deg as f64 * 60.0;
        format!("{}d {:04.1}\"", deg, rem)
    } else {
        let whole = arcmin as i32;
        let arcsec = (arcmin - whole as f64) * 60.0;
        format!("{}' {:04.1}\"", whole, arcsec)
    }
}

// A TAN WCS for tests of a `width` x `height` image: CRVAL `crval`, CRPIX at
// (width / 2, height / 2), the diagonal CD matrix `(cd1_1, cd2_2)` in deg/px, and no fit
// or camera data.
#[cfg(test)]
impl WcsSolution {
    pub(crate) fn diagonal(
        crval: (f64, f64),
        (cd1_1, cd2_2): (f64, f64),
        (width, height): (usize, usize),
    ) -> Self {
        Self {
            crpix1: width as f64 / 2.0,
            crpix2: height as f64 / 2.0,
            crval1: crval.0,
            crval2: crval.1,
            cd1_1,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2,
            width,
            height,
            focal_mm: None,
            pixel_um: None,
            n_stars: 0,
            rms_px: 0.0,
            weighted_rms_px: 0.0,
            residuals: Vec::new(),
        }
    }

    // North up and east left at `scale_arcsec` per pixel.
    pub(crate) fn north_up(crval: (f64, f64), scale_arcsec: f64, size: (usize, usize)) -> Self {
        let s = scale_arcsec / 3600.0;
        Self::diagonal(crval, (-s, s), size)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_wcs() -> WcsSolution {
        WcsSolution {
            focal_mm: Some(1000.0),
            pixel_um: Some(3.76),
            ..WcsSolution::diagonal((180.0, 0.0), (0.001, 0.001), (2048, 1536))
        }
    }

    #[test]
    fn pixel_sky_round_trip_recovers_input() {
        let wcs = sample_wcs();
        let test_points = [
            (1024.0, 768.0),
            (1100.0, 800.0),
            (900.0, 700.0),
            (1024.0, 900.0),
        ];
        let round_trip: Vec<(f64, f64)> = test_points
            .iter()
            .map(|&(px_x, px_y)| {
                let (ra, dec) = wcs.pixel_to_sky(px_x, px_y);
                wcs.sky_to_pixel(ra, dec).unwrap()
            })
            .collect();
        // Two points come back up to 6.6e-12 px from where they started in x.
        assert_eq!(
            round_trip,
            [
                (1024.0, 768.0),
                (1099.9999999999934, 800.0),
                (899.9999999999942, 700.0),
                (1024.0, 900.0)
            ]
        );
    }

    #[test]
    fn pixel_to_sky_at_crpix_returns_crval() {
        let wcs = sample_wcs();
        let (ra, dec) = wcs.pixel_to_sky(wcs.crpix1, wcs.crpix2);
        // crpix is the tangent point; projection of (0, 0) offset is exactly crval
        assert_eq!(ra, wcs.crval1);
        assert_eq!(dec, wcs.crval2);
    }

    #[test]
    fn scale_arcsec_from_diagonal_cd_matches_column_norm_times_3600() {
        let mut wcs = sample_wcs();
        wcs.cd1_1 = 0.001;
        wcs.cd1_2 = 0.0;
        wcs.cd2_1 = 0.0;
        wcs.cd2_2 = 0.002;
        let (sx, sy) = wcs.scale_arcsec();
        // sqrt(cd1_1^2 + cd2_1^2) * 3600 = 3.6
        assert_eq!(sx, 3.6);
        // sqrt(cd1_2^2 + cd2_2^2) * 3600 = 7.2
        assert_eq!(sy, 7.2);
    }

    // Pixel y grows downwards, so an identity CD puts north down and east right: a half
    // turn from north up, east left.
    #[test]
    fn rotation_deg_identity_is_a_half_turn() {
        let wcs = sample_wcs();
        // atan2(-cd1_2, -cd2_2) = atan2(-0, -0.001) = -180, reported as 180
        assert_eq!(wcs.rotation_deg(), 180.0);
    }

    #[test]
    fn rotation_deg_quarter_turn_is_minus_ninety() {
        let mut wcs = sample_wcs();
        wcs.cd1_1 = 0.0;
        wcs.cd1_2 = 0.001;
        wcs.cd2_1 = -0.001;
        wcs.cd2_2 = 0.0;
        // det > 0 so not flipped; atan2(-0.001, -0) * RAD_TO_DEG = -90
        assert_eq!(wcs.rotation_deg(), -90.0);
    }

    fn with_cd(cd: [f64; 4]) -> WcsSolution {
        let mut wcs = sample_wcs();
        (wcs.cd1_1, wcs.cd1_2, wcs.cd2_1, wcs.cd2_2) = (cd[0], cd[1], cd[2], cd[3]);
        wcs
    }

    // Unflipped CD for rotation t, with pixel y down: -s * [[cos t, sin t], [-sin t, cos t]].
    // At t = 0 east is left (cd1_1 < 0) and north is up (cd2_2 < 0).
    #[test]
    fn rotation_deg_unflipped_quarter_turns() {
        let s = 0.001;
        let cases = [
            ([-s, 0.0, 0.0, -s], 0.0),
            ([0.0, -s, s, 0.0], 90.0),
            ([0.0, s, -s, 0.0], -90.0),
            ([s, 0.0, 0.0, s], 180.0),
        ];
        for (cd, rotation) in cases {
            let wcs = with_cd(cd);
            assert!(!wcs.is_flipped());
            assert_eq!(wcs.rotation_deg(), rotation, "CD {cd:?}");
        }
    }

    // Flipped CD for rotation t, with pixel y down: s * [[cos t, sin t], [sin t, -cos t]].
    #[test]
    fn rotation_deg_flipped_quarter_turns() {
        let s = 0.001;
        let cases = [
            ([s, 0.0, 0.0, -s], 0.0),
            ([0.0, s, s, 0.0], 90.0),
            ([0.0, -s, -s, 0.0], -90.0),
            ([-s, 0.0, 0.0, s], 180.0),
        ];
        for (cd, rotation) in cases {
            let wcs = with_cd(cd);
            assert!(wcs.is_flipped());
            assert_eq!(wcs.rotation_deg(), rotation, "CD {cd:?}");
        }
    }

    #[test]
    fn rotation_deg_flipped_uses_negated_cd2_2() {
        let mut wcs = sample_wcs();
        // Flip the x-axis: cd1_1 negative, rest unchanged -> det negative
        wcs.cd1_1 = -0.001;
        assert!(wcs.is_flipped());
        // flipped branch: atan2(cd1_2=0, -cd2_2=-0.001) = atan2(0, -1) = pi -> 180°
        assert_eq!(wcs.rotation_deg(), 180.0);
    }

    #[test]
    fn is_flipped_follows_determinant_sign() {
        let mut wcs = sample_wcs();
        assert!(!wcs.is_flipped());
        wcs.cd1_1 = -0.001;
        assert!(wcs.is_flipped());
    }

    #[test]
    fn effective_focal_mm_none_when_pixel_um_missing() {
        let mut wcs = sample_wcs();
        wcs.pixel_um = None;
        assert!(wcs.effective_focal_mm().is_none());
    }

    #[test]
    fn effective_focal_mm_recovers_input_focal() {
        // with pixel_um = 3.76 and scale = 3.6 arcsec/px (cd=0.001 deg/px = 3.6"/px)
        // effective focal = 206.2648 * 3.76 / 3.6 = 215.43...
        let mut wcs = sample_wcs();
        wcs.cd1_1 = 0.001;
        wcs.cd1_2 = 0.0;
        wcs.cd2_1 = 0.0;
        wcs.cd2_2 = 0.001;
        wcs.pixel_um = Some(3.76);
        assert_eq!(
            wcs.effective_focal_mm(),
            Some(ARCSEC_PER_UM_PER_MM * 3.76 / 3.6)
        );
    }

    #[test]
    fn effective_focal_mm_none_when_scale_zero() {
        let mut wcs = sample_wcs();
        wcs.cd1_1 = 0.0;
        wcs.cd1_2 = 0.0;
        wcs.cd2_1 = 0.0;
        wcs.cd2_2 = 0.0;
        wcs.pixel_um = Some(3.76);
        assert!(wcs.effective_focal_mm().is_none());
    }

    #[test]
    fn nearest_residual_averages_k_nearest_errors() {
        let mut wcs = sample_wcs();
        wcs.residuals = vec![
            StarResidual {
                px_x: 100.0,
                px_y: 100.0,
                err_x: 1.0,
                err_y: 0.0,
                err_px: 1.0,
                ra_deg: 0.0,
                dec_deg: 0.0,
            },
            StarResidual {
                px_x: 200.0,
                px_y: 100.0,
                err_x: 0.0,
                err_y: 2.0,
                err_px: 2.0,
                ra_deg: 0.0,
                dec_deg: 0.0,
            },
            StarResidual {
                px_x: 500.0,
                px_y: 500.0,
                err_x: 100.0,
                err_y: 100.0,
                err_px: 141.4,
                ra_deg: 0.0,
                dec_deg: 0.0,
            },
        ];
        // Query near the first two residuals; k=2 should exclude the far one
        let (ex, ey) = wcs.nearest_residual(150.0, 100.0, 2);
        // mean of (1.0, 0.0) and (0.0, 2.0) = (0.5, 1.0)
        assert_eq!(ex, 0.5);
        assert_eq!(ey, 1.0);
    }

    #[test]
    fn nearest_residual_empty_returns_zero() {
        let wcs = sample_wcs();
        let (ex, ey) = wcs.nearest_residual(0.0, 0.0, 5);
        assert_eq!(ex, 0.0);
        assert_eq!(ey, 0.0);
    }

    #[test]
    fn nearest_residual_clamps_k_to_available_count() {
        let mut wcs = sample_wcs();
        wcs.residuals = vec![StarResidual {
            px_x: 0.0,
            px_y: 0.0,
            err_x: 4.0,
            err_y: 6.0,
            err_px: 7.2,
            ra_deg: 0.0,
            dec_deg: 0.0,
        }];
        // request 10 nearest but only 1 exists
        let (ex, ey) = wcs.nearest_residual(0.0, 0.0, 10);
        assert_eq!(ex, 4.0);
        assert_eq!(ey, 6.0);
    }

    #[test]
    fn fmt_fov_component_small_uses_arcmin_arcsec() {
        // < 60 arcmin: "X' YY.Y""
        let s = fmt_fov_component(30.5);
        assert!(s.contains("30'"));
        assert!(s.contains("30.0\""));
    }

    #[test]
    fn fmt_fov_component_large_uses_degrees() {
        // >= 60 arcmin: "Xd YY.Y""
        let s = fmt_fov_component(120.0);
        assert!(s.contains("2d"));
    }
}
