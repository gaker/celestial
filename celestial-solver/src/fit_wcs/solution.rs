//! WCS solution data and pixel↔sky transforms.

use std::fmt;

use celestial_catalog::query::{tan_deproject_star, tan_project_star};
use celestial_core::constants::{DEG_TO_RAD, RAD_TO_DEG};

/// Per-pair fit residual.
///
/// Computed by [`crate::fit_wcs::fit_wcs`] for each inlier pair. `err_x` and `err_y`
/// are `detected - predicted` in pixel space — positive values mean the detected
/// centroid landed to the right of / below where the WCS says the catalog star should
/// be. `err_px` is the Euclidean magnitude.
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve(&img, &catalog).run()?;
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
}

/// Fitted TAN (gnomonic) WCS with residual statistics.
///
/// Standard FITS WCS parameters plus fit diagnostics. Pixel↔sky transforms via
/// [`WcsSolution::pixel_to_sky`] and [`WcsSolution::sky_to_pixel`]. Geometric
/// properties via [`WcsSolution::scale_arcsec`], [`WcsSolution::rotation_deg`],
/// [`WcsSolution::is_flipped`].
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve(&img, &catalog).run()?;
/// let (ra, dec) = result.wcs.pixel_to_sky(1024.0, 768.0);
/// println!("center: {:.4}° {:.4}°", ra, dec);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct WcsSolution {
    /// Reference pixel x (1-based FITS convention). Usually the image center.
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
    /// SNR-weighted RMS residual in pixels.
    pub weighted_rms_px: f64,
    /// Per-pair residuals for the final fit.
    pub residuals: Vec<StarResidual>,
}

impl WcsSolution {
    /// Transforms a pixel coordinate to sky coordinates (RA, Dec in degrees).
    ///
    /// Applies the inverse CD matrix, then de-projects from the tangent plane. Pure
    /// linear TAN — does not account for SIP distortion. Use [`crate::fit_wcs::SipSolution`]
    /// and its forward evaluator for SIP-corrected transforms.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// let (ra, dec) = result.wcs.pixel_to_sky(100.0, 200.0);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn pixel_to_sky(&self, px_x: f64, px_y: f64) -> (f64, f64) {
        let u = px_x - self.crpix1;
        let v = px_y - self.crpix2;
        let xi_deg = self.cd1_1 * u + self.cd1_2 * v;
        let eta_deg = self.cd2_1 * u + self.cd2_2 * v;
        tan_deproject_star(xi_deg * DEG_TO_RAD, eta_deg * DEG_TO_RAD, self.crval1, self.crval2)
    }

    /// Transforms sky coordinates to pixel coordinates.
    ///
    /// Returns `None` for points that can't be gnomonically projected (e.g. the
    /// anti-center of the TAN projection, more than 90° from CRVAL). Pure linear TAN,
    /// no SIP correction.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// if let Some((x, y)) = result.wcs.sky_to_pixel(180.0, 30.0) {
    ///     println!("(180°, 30°) lands at ({:.1}, {:.1})", x, y);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn sky_to_pixel(&self, ra_deg: f64, dec_deg: f64) -> Option<(f64, f64)> {
        let (xi_rad, eta_rad) = tan_project_star(ra_deg, dec_deg, self.crval1, self.crval2)?;
        let xi = xi_rad * RAD_TO_DEG;
        let eta = eta_rad * RAD_TO_DEG;
        let det = self.cd1_1 * self.cd2_2 - self.cd1_2 * self.cd2_1;
        let px_x = (self.cd2_2 * xi - self.cd1_2 * eta) / det + self.crpix1;
        let px_y = (-self.cd2_1 * xi + self.cd1_1 * eta) / det + self.crpix2;
        Some((px_x, px_y))
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
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// let (sx, sy) = result.wcs.scale_arcsec();
    /// println!("{:.3} × {:.3} arcsec/px", sx, sy);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn scale_arcsec(&self) -> (f64, f64) {
        let sx = libm::sqrt(self.cd1_1 * self.cd1_1 + self.cd2_1 * self.cd2_1) * 3600.0;
        let sy = libm::sqrt(self.cd1_2 * self.cd1_2 + self.cd2_2 * self.cd2_2) * 3600.0;
        (sx, sy)
    }

    /// Field rotation in degrees, accounting for image flip.
    ///
    /// Zero degrees means north-up, east-left (the astronomical convention). Positive
    /// rotation is counter-clockwise.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// println!("rotation: {:.2}°", result.wcs.rotation_deg());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn rotation_deg(&self) -> f64 {
        if self.is_flipped() {
            libm::atan2(self.cd1_2, -self.cd2_2) * RAD_TO_DEG
        } else {
            libm::atan2(self.cd1_2, self.cd2_2) * RAD_TO_DEG
        }
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
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
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
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
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
        Some(206.265 * pixel_um / avg_scale)
    }

    pub(crate) fn nearest_residual(&self, px: f64, py: f64, k: usize) -> (f64, f64) {
        if self.residuals.is_empty() {
            return (0.0, 0.0);
        }
        let mut distances: Vec<(usize, f64)> = self.residuals
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

    pub(crate) fn fmt_bounds(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let w = self.width as f64;
        let h = self.height as f64;
        let corners = [
            ("top-left", 0.0, 0.0),
            ("top-right", w, 0.0),
            ("bottom-left", 0.0, h),
            ("bottom-right", w, h),
        ];
        writeln!(f, "Image bounds:")?;
        for (label, px, py) in &corners {
            let (ra, dec) = self.pixel_to_sky(*px, *py);
            let (ex, ey) = self.nearest_residual(*px, *py, 5);
            writeln!(f, "   {:<17} RA: {}  Dec: {}  ex: {:+.6} px  ey: {:+.6} px",
                label, fmt_ra(ra), fmt_dec(dec), ex, ey)?;
        }
        Ok(())
    }
}

pub(crate) fn fmt_ra(ra_deg: f64) -> String {
    let ra_h = ra_deg / 15.0;
    let hh = ra_h as i32;
    let mm = ((ra_h - hh as f64) * 60.0) as i32;
    let ss = (ra_h - hh as f64 - mm as f64 / 60.0) * 3600.0;
    format!("{:02} {:02} {:06.3}", hh, mm, ss)
}

pub(crate) fn fmt_dec(dec_deg: f64) -> String {
    let sign = if dec_deg < 0.0 { "-" } else { "+" };
    let abs = libm::fabs(dec_deg);
    let dd = abs as i32;
    let mm = ((abs - dd as f64) * 60.0) as i32;
    let ss = (abs - dd as f64 - mm as f64 / 60.0) * 3600.0;
    format!("{}{:02} {:02} {:05.2}", sign, dd, mm, ss)
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

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_wcs() -> WcsSolution {
        WcsSolution {
            crpix1: 1024.0,
            crpix2: 768.0,
            crval1: 180.0,
            crval2: 0.0,
            cd1_1: 0.001,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: 0.001,
            width: 2048,
            height: 1536,
            focal_mm: Some(1000.0),
            pixel_um: Some(3.76),
            n_stars: 0,
            rms_px: 0.0,
            weighted_rms_px: 0.0,
            residuals: Vec::new(),
        }
    }

    #[test]
    fn pixel_sky_round_trip_recovers_input() {
        let wcs = sample_wcs();
        let test_points = [(1024.0, 768.0), (1100.0, 800.0), (900.0, 700.0), (1024.0, 900.0)];
        for (px_x, px_y) in test_points {
            let (ra, dec) = wcs.pixel_to_sky(px_x, px_y);
            let (px_x2, px_y2) = wcs.sky_to_pixel(ra, dec).unwrap();
            // round trip passes through tan projection (f64 trig), so we
            // allow a very tight tolerance rather than exact equality.
            assert!((px_x - px_x2).abs() < 1e-9, "x {} -> {}", px_x, px_x2);
            assert!((px_y - px_y2).abs() < 1e-9, "y {} -> {}", px_y, px_y2);
        }
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

    #[test]
    fn rotation_deg_identity_is_zero() {
        let wcs = sample_wcs();
        // cd1_2 = 0, cd2_2 > 0 -> atan2(0, 1) = 0
        assert_eq!(wcs.rotation_deg(), 0.0);
    }

    #[test]
    fn rotation_deg_quarter_turn_is_ninety() {
        let mut wcs = sample_wcs();
        // 90° rotation: cd1_1=0, cd1_2=scale, cd2_1=-scale, cd2_2=0
        wcs.cd1_1 = 0.0;
        wcs.cd1_2 = 0.001;
        wcs.cd2_1 = -0.001;
        wcs.cd2_2 = 0.0;
        // det > 0 so not flipped; atan2(0.001, 0) * RAD_TO_DEG = 90
        assert!((wcs.rotation_deg() - 90.0).abs() < 1e-12);
    }

    #[test]
    fn rotation_deg_flipped_uses_negated_cd2_2() {
        let mut wcs = sample_wcs();
        // Flip the x-axis: cd1_1 negative, rest unchanged -> det negative
        wcs.cd1_1 = -0.001;
        assert!(wcs.is_flipped());
        // flipped branch: atan2(cd1_2=0, -cd2_2=-0.001) = atan2(0, -1) = pi -> 180°
        assert!((wcs.rotation_deg().abs() - 180.0).abs() < 1e-12);
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
        // effective focal = 206.265 * 3.76 / 3.6 = 215.39...
        let mut wcs = sample_wcs();
        wcs.cd1_1 = 0.001;
        wcs.cd1_2 = 0.0;
        wcs.cd2_1 = 0.0;
        wcs.cd2_2 = 0.001;
        wcs.pixel_um = Some(3.76);
        let f = wcs.effective_focal_mm().unwrap();
        let expected = 206.265 * 3.76 / 3.6;
        assert!((f - expected).abs() < 1e-9, "focal {} vs {}", f, expected);
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
            StarResidual { px_x: 100.0, px_y: 100.0, err_x: 1.0, err_y: 0.0, err_px: 1.0 },
            StarResidual { px_x: 200.0, px_y: 100.0, err_x: 0.0, err_y: 2.0, err_px: 2.0 },
            StarResidual { px_x: 500.0, px_y: 500.0, err_x: 100.0, err_y: 100.0, err_px: 141.4 },
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
            px_x: 0.0, px_y: 0.0, err_x: 4.0, err_y: 6.0, err_px: 7.2,
        }];
        // request 10 nearest but only 1 exists
        let (ex, ey) = wcs.nearest_residual(0.0, 0.0, 10);
        assert_eq!(ex, 4.0);
        assert_eq!(ey, 6.0);
    }

    #[test]
    fn fmt_ra_formats_decimal_degrees_as_hms() {
        // 180° = 12h
        assert_eq!(fmt_ra(180.0), "12 00 00.000");
        // 0° = 00h
        assert_eq!(fmt_ra(0.0), "00 00 00.000");
        // 15° = 1h
        assert_eq!(fmt_ra(15.0), "01 00 00.000");
    }

    #[test]
    fn fmt_dec_shows_sign_and_dms() {
        assert_eq!(fmt_dec(0.0), "+00 00 00.00");
        assert_eq!(fmt_dec(45.5), "+45 30 00.00");
        assert_eq!(fmt_dec(-30.25), "-30 15 00.00");
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
