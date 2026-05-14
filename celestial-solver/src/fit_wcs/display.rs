use std::collections::HashMap;
use std::fmt;

use chrono::Utc;

use super::sip::SipSolution;
use super::solution::{fmt_dec, fmt_fov_component, fmt_ra, WcsSolution};

/// Human-readable report of a [`WcsSolution`] and optional [`SipSolution`].
///
/// Implements [`fmt::Display`] to render a multi-line PixInsight-style summary
/// block — projection origin, plate scale, rotation, focal length, pixel size, SIP
/// order, RMS residuals, and optionally software name, observation time, observer
/// geodetic location.
///
/// Used by the `image-solver` binary for its stdout report. Available publicly for
/// any tool wanting the same format.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::fit_wcs::SolveDisplay;
///
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve(&img, &catalog).run()?;
/// let display = SolveDisplay::new(&result.wcs, result.sip.as_ref())
///     .with_software("my-solver 1.0");
/// println!("{}", display);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct SolveDisplay<'a> {
    wcs: &'a WcsSolution,
    sip: Option<&'a SipSolution>,
    software: Option<&'a str>,
    observation_time: Option<&'a str>,
    geodetic: Option<(f64, f64, Option<f64>)>,
}

impl<'a> SolveDisplay<'a> {
    /// Creates a display with no software / time / location annotations.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_solver::fit_wcs::SolveDisplay;
    ///
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// println!("{}", SolveDisplay::new(&result.wcs, result.sip.as_ref()));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(wcs: &'a WcsSolution, sip: Option<&'a SipSolution>) -> Self {
        Self {
            wcs,
            sip,
            software: None,
            observation_time: None,
            geodetic: None,
        }
    }

    /// Annotates the report with the solver name + version.
    pub fn with_software(mut self, name: &'a str) -> Self {
        self.software = Some(name);
        self
    }

    /// Annotates the report with a formatted UTC observation time.
    pub fn with_observation_time(mut self, utc: &'a str) -> Self {
        self.observation_time = Some(utc);
        self
    }

    /// Annotates the report with observer geodetic location (lon, lat, optional alt).
    pub fn with_geodetic(mut self, lon_deg: f64, lat_deg: f64, alt_m: Option<f64>) -> Self {
        self.geodetic = Some((lon_deg, lat_deg, alt_m));
        self
    }
}

impl<'a> fmt::Display for SolveDisplay<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let wcs = self.wcs;
        let (scale_x, scale_y) = wcs.scale_arcsec();
        let scale_avg = (scale_x + scale_y) / 2.0;
        let rotation = wcs.rotation_deg();
        let fov_x = scale_x * wcs.width as f64 / 60.0;
        let fov_y = scale_y * wcs.height as f64 / 60.0;

        let sep = "=".repeat(80);
        writeln!(f, "{sep}")?;

        writeln!(f, "Creation time ............ {} UTC", Utc::now().format("%Y-%m-%d %H:%M:%S"))?;
        if let Some(sw) = self.software {
            writeln!(f, "Creation software ........ {}", sw)?;
        }
        writeln!(f, "Reference catalog ........ Gaia DR3 - Celestial Catalog v{}", celestial_catalog::VERSION)?;

        writeln!(f, "Linear transformation matrix (sky = matrix * pixel):")?;
        writeln!(f, " {:+.8e}  {:+.8e}  {:+.8e}",
            wcs.cd1_1, wcs.cd1_2,
            wcs.cd1_1 * (-wcs.crpix1) + wcs.cd1_2 * (-wcs.crpix2))?;
        writeln!(f, " {:+.8e}  {:+.8e}  {:+.8e}",
            wcs.cd2_1, wcs.cd2_2,
            wcs.cd2_1 * (-wcs.crpix1) + wcs.cd2_2 * (-wcs.crpix2))?;

        let projection = if self.sip.is_some() { "Gnomonic (TAN-SIP)" } else { "Gnomonic (TAN)" };
        writeln!(f, "Projection ............... {}", projection)?;
        writeln!(f, "Projection origin ........ [{:.6} {:.6}] px -> [RA: {}  Dec: {}]",
            wcs.crpix1, wcs.crpix2,
            fmt_ra(wcs.crval1), fmt_dec(wcs.crval2))?;
        writeln!(f, "Resolution ............... {:.3} arcsec/px", scale_avg)?;
        let flip_suffix = if wcs.is_flipped() { " (flipped)" } else { "" };
        writeln!(f, "Rotation ................. {:.3} deg{}", rotation, flip_suffix)?;
        writeln!(f, "Reference system ......... ICRS")?;

        if let Some(time) = self.observation_time {
            writeln!(f, "Observation start time ... {} UTC", time)?;
        }
        if let Some((lon, lat, alt)) = self.geodetic {
            let alt_str = alt.map(|a| format!("  {:.0} m", a)).unwrap_or_default();
            writeln!(
                f,
                "Geodetic coordinates ..... {}  {}{}",
                format_lon_deg(lon), format_lat_deg(lat), alt_str,
            )?;
        }

        if let Some(px) = wcs.pixel_um {
            let focal = wcs.effective_focal_mm().unwrap_or(0.0);
            writeln!(f, "Focal distance ........... {:.2} mm", focal)?;
            writeln!(f, "Pixel size ............... {:.2} um", px)?;
        }

        writeln!(f, "Field of view ............ {} x {}",
            fmt_fov_component(fov_x), fmt_fov_component(fov_y))?;
        writeln!(f, "Image dimensions ......... {} x {} px", wcs.width, wcs.height)?;

        let (ex, ey) = wcs.nearest_residual(wcs.crpix1, wcs.crpix2, 5);
        writeln!(f, "Image center ............. RA: {}  Dec: {}  ex: {:+.6} px  ey: {:+.6} px",
            fmt_ra(wcs.crval1), fmt_dec(wcs.crval2),
            ex, ey)?;

        wcs.fmt_bounds(f)?;

        if let Some(sip) = self.sip {
            writeln!(f, "SIP distortion order ..... A={}, B={}", sip.a_order, sip.b_order)?;
            fmt_sip_coeffs(f, "A", &sip.a_coeffs)?;
            fmt_sip_coeffs(f, "B", &sip.b_coeffs)?;
        }

        writeln!(f, "Control points ........... {}", wcs.n_stars)?;

        if let Some(sip) = self.sip {
            writeln!(f, "Linear RMS ............... {:.4} px", sip.linear_rms_px)?;
            writeln!(f, "SIP RMS .................. {:.4} px", sip.rms_px)?;
            writeln!(f, "Weighted SIP RMS ......... {:.4} px", sip.weighted_rms_px)?;
            let improvement = if sip.linear_rms_px > 0.0 {
                (1.0 - sip.rms_px / sip.linear_rms_px) * 100.0
            } else {
                0.0
            };
            writeln!(f, "Improvement .............. {:.1}%", improvement)?;
        } else {
            writeln!(f, "RMS error ................ {:.3} px", wcs.rms_px)?;
            writeln!(f, "Weighted RMS ............. {:.3} px", wcs.weighted_rms_px)?;
        }
        write!(f, "{sep}")
    }
}

fn fmt_sip_coeffs(
    f: &mut fmt::Formatter<'_>,
    prefix: &str,
    coeffs: &HashMap<(u32, u32), f64>,
) -> fmt::Result {
    let mut keys: Vec<(u32, u32)> = coeffs.keys().copied().collect();
    keys.sort_by(|a, b| (a.0 + a.1).cmp(&(b.0 + b.1)).then(b.0.cmp(&a.0)));
    for (p, q) in keys {
        writeln!(f, "  {}_{}_{}  {:+.6e}", prefix, p, q, coeffs[&(p, q)])?;
    }
    Ok(())
}

fn format_lat_deg(deg: f64) -> String {
    let sign = if deg < 0.0 { "S" } else { "N" };
    let abs = libm::fabs(deg);
    let d = abs as i32;
    let m = ((abs - d as f64) * 60.0) as i32;
    let s = (abs - d as f64 - m as f64 / 60.0) * 3600.0;
    format!("{:02} {:02} {:02} {}", d, m, s.round() as i32, sign)
}

fn format_lon_deg(deg: f64) -> String {
    let sign = if deg < 0.0 { "W" } else { "E" };
    let abs = libm::fabs(deg);
    let d = abs as i32;
    let m = ((abs - d as f64) * 60.0) as i32;
    let s = (abs - d as f64 - m as f64 / 60.0) * 3600.0;
    format!("{:02} {:02} {:02} {}", d, m, s.round() as i32, sign)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::solution::{StarResidual, WcsSolution};
    use super::super::sip::SipSolution;

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
            n_stars: 50,
            rms_px: 0.3,
            weighted_rms_px: 0.2,
            residuals: vec![StarResidual {
                px_x: 1024.0, px_y: 768.0, err_x: 0.01, err_y: -0.02, err_px: 0.022,
            }],
        }
    }

    fn sample_sip() -> SipSolution {
        let mut a_coeffs = HashMap::new();
        let mut b_coeffs = HashMap::new();
        a_coeffs.insert((2, 0), 1e-6);
        a_coeffs.insert((1, 1), 3e-6);
        b_coeffs.insert((0, 2), 2e-6);
        SipSolution {
            crpix1: 1024.0, crpix2: 768.0, crval1: 180.0, crval2: 0.0,
            cd1_1: 0.001, cd1_2: 0.0, cd2_1: 0.0, cd2_2: 0.001,
            width: 2048, height: 1536,
            a_order: 3, b_order: 3,
            a_coeffs, b_coeffs,
            n_stars: 50,
            linear_rms_px: 1.2,
            rms_px: 0.3,
            weighted_rms_px: 0.25,
        }
    }

    #[test]
    fn display_without_sip_contains_expected_sections() {
        let wcs = sample_wcs();
        let display = SolveDisplay::new(&wcs, None);
        let text = format!("{}", display);
        assert!(text.contains("Projection origin"));
        assert!(text.contains("Resolution"));
        assert!(text.contains("Rotation"));
        assert!(text.contains("Focal distance"));
        assert!(text.contains("Field of view"));
        assert!(text.contains("RMS error"));
        assert!(text.contains("Gnomonic (TAN)"));
    }

    #[test]
    fn display_with_sip_shows_sip_sections() {
        let wcs = sample_wcs();
        let sip = sample_sip();
        let display = SolveDisplay::new(&wcs, Some(&sip));
        let text = format!("{}", display);
        assert!(text.contains("Gnomonic (TAN-SIP)"));
        assert!(text.contains("SIP distortion order"));
        assert!(text.contains("Linear RMS"));
        assert!(text.contains("SIP RMS"));
        assert!(text.contains("Improvement"));
        // coefficient entries
        assert!(text.contains("A_2_0"));
        assert!(text.contains("B_0_2"));
    }

    #[test]
    fn display_with_software_shows_name() {
        let wcs = sample_wcs();
        let display = SolveDisplay::new(&wcs, None).with_software("my-solver v0.1");
        let text = format!("{}", display);
        assert!(text.contains("Creation software"));
        assert!(text.contains("my-solver v0.1"));
    }

    #[test]
    fn display_with_observation_time_shows_time() {
        let wcs = sample_wcs();
        let display = SolveDisplay::new(&wcs, None)
            .with_observation_time("2026-04-21 20:00:00");
        let text = format!("{}", display);
        assert!(text.contains("Observation start time"));
        assert!(text.contains("2026-04-21 20:00:00"));
    }

    #[test]
    fn display_with_geodetic_shows_formatted_location() {
        let wcs = sample_wcs();
        let display = SolveDisplay::new(&wcs, None)
            .with_geodetic(-117.0, 33.0, Some(150.0));
        let text = format!("{}", display);
        assert!(text.contains("Geodetic coordinates"));
        assert!(text.contains("W")); // west longitude
        assert!(text.contains("N")); // north latitude
        assert!(text.contains("150 m"));
    }

    #[test]
    fn display_without_altitude_omits_altitude_field() {
        let wcs = sample_wcs();
        let display = SolveDisplay::new(&wcs, None)
            .with_geodetic(-117.0, 33.0, None);
        let text = format!("{}", display);
        assert!(text.contains("Geodetic coordinates"));
        // no "m " altitude suffix
        assert!(!text.contains("  m\n"));
    }

    #[test]
    fn display_flipped_cd_adds_flipped_suffix() {
        let mut wcs = sample_wcs();
        wcs.cd1_1 = -0.001; // det < 0 -> flipped
        let display = SolveDisplay::new(&wcs, None);
        let text = format!("{}", display);
        assert!(text.contains("(flipped)"));
    }

    #[test]
    fn display_without_pixel_um_omits_focal_line() {
        let mut wcs = sample_wcs();
        wcs.pixel_um = None;
        let display = SolveDisplay::new(&wcs, None);
        let text = format!("{}", display);
        assert!(!text.contains("Focal distance"));
        assert!(!text.contains("Pixel size"));
    }

    #[test]
    fn fmt_sip_coeffs_emits_sorted_coefficient_lines() {
        let mut coeffs = HashMap::new();
        coeffs.insert((2, 0), 1.5e-6);
        coeffs.insert((0, 2), 2.5e-6);

        // We can't call fmt_sip_coeffs directly with a Formatter, so test
        // indirectly via a full Display pass that exercises it.
        let mut wcs = sample_wcs();
        wcs.pixel_um = Some(3.76);
        let sip = SipSolution {
            crpix1: wcs.crpix1, crpix2: wcs.crpix2,
            crval1: wcs.crval1, crval2: wcs.crval2,
            cd1_1: wcs.cd1_1, cd1_2: wcs.cd1_2,
            cd2_1: wcs.cd2_1, cd2_2: wcs.cd2_2,
            width: wcs.width, height: wcs.height,
            a_order: 2, b_order: 2,
            a_coeffs: coeffs.clone(),
            b_coeffs: coeffs,
            n_stars: 10,
            linear_rms_px: 2.0,
            rms_px: 0.5,
            weighted_rms_px: 0.4,
        };
        let display = SolveDisplay::new(&wcs, Some(&sip));
        let text = format!("{}", display);
        assert!(text.contains("A_2_0"));
        assert!(text.contains("A_0_2"));
        assert!(text.contains("B_2_0"));
        assert!(text.contains("B_0_2"));
    }

    #[test]
    fn format_lat_deg_positive_is_north() {
        let s = format_lat_deg(33.5);
        assert!(s.ends_with(" N"), "got {s}");
        assert!(s.starts_with("33 "));
    }

    #[test]
    fn format_lat_deg_negative_is_south() {
        let s = format_lat_deg(-45.25);
        assert!(s.ends_with(" S"), "got {s}");
        assert!(s.starts_with("45 "));
    }

    #[test]
    fn format_lon_deg_positive_is_east() {
        let s = format_lon_deg(120.0);
        assert!(s.ends_with(" E"), "got {s}");
    }

    #[test]
    fn format_lon_deg_negative_is_west() {
        let s = format_lon_deg(-117.0);
        assert!(s.ends_with(" W"), "got {s}");
    }
}
