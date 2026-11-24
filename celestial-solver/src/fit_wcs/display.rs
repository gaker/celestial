mod projection_model;
mod report_time;

use std::collections::BTreeMap;
use std::fmt;

use celestial_time::scales::utc::UTC;

use super::sexagesimal::{fmt_dec, fmt_ra, format_lat_deg, format_lon_deg};
use super::solution::{fmt_fov_component, WcsSolution};
use crate::pixel_frame::{
    detection_to_pixinsight, pixinsight_linear_transformation, pixinsight_to_detection,
};
use crate::solve::SolveResult;
use projection_model::projection_model;
use report_time::to_the_second;

/// Human-readable report of a [`SolveResult`]: its
/// [`WcsSolution`] and optional
/// [`SipSolution`](crate::fit_wcs::sip::SipSolution).
///
/// Implements [`fmt::Display`] to render a multi-line PixInsight-style summary
/// block — projection origin, plate scale, rotation, focal length, pixel size, SIP
/// order, RMS residuals, and optionally creation time, software name, observation time,
/// observer geodetic location. Pixel coordinates are in PixInsight's frame, with the centre of
/// the first pixel at (0.5, 0.5), and sky positions include the SIP distortion.
///
/// Used by the `image-solver` binary for its stdout report. Available publicly for
/// any tool wanting the same format.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::fit_wcs::display::SolveDisplay;
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
/// let display = SolveDisplay::new(&result).with_software("my-solver 1.0");
/// println!("{}", display);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct SolveDisplay<'a> {
    result: &'a SolveResult,
    creation_time: Option<UTC>,
    software: Option<&'a str>,
    observation_time: Option<UTC>,
    geodetic: Option<(f64, f64, Option<f64>)>,
}

impl<'a> SolveDisplay<'a> {
    /// Creates a display with no creation time / software / time / location annotations.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_solver::fit_wcs::display::SolveDisplay;
    ///
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// println!("{}", SolveDisplay::new(&result));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(result: &'a SolveResult) -> Self {
        Self {
            result,
            creation_time: None,
            software: None,
            observation_time: None,
            geodetic: None,
        }
    }

    pub fn with_creation_time(mut self, utc: UTC) -> Self {
        self.creation_time = Some(utc);
        self
    }

    /// Annotates the report with the solver name + version.
    pub fn with_software(mut self, name: &'a str) -> Self {
        self.software = Some(name);
        self
    }

    /// Annotates the report with the observation start time, printed to the second.
    pub fn with_observation_time(mut self, utc: UTC) -> Self {
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
        let model = projection_model(self.result);
        let sep = "=".repeat(80);
        writeln!(f, "{sep}")?;
        self.fmt_provenance(f)?;
        self.fmt_projection(f, &model)?;
        self.fmt_scale_and_rotation(f)?;
        self.fmt_observation(f)?;
        self.fmt_field(f, &model)?;
        self.fmt_sip(f)?;
        self.fmt_residuals(f)?;
        write!(f, "{sep}")
    }
}

// The sections of the report, in the order `fmt` writes them.
impl<'a> SolveDisplay<'a> {
    fn fmt_provenance(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(time) = &self.creation_time {
            writeln!(f, "Creation time ............ {} UTC", to_the_second(time))?;
        }
        if let Some(sw) = self.software {
            writeln!(f, "Creation software ........ {}", sw)?;
        }
        writeln!(
            f,
            "Reference catalog ........ Gaia DR3 - Celestial Catalog v{}",
            celestial_catalog::VERSION
        )
    }

    fn fmt_projection(&self, f: &mut fmt::Formatter<'_>, model: &WcsSolution) -> fmt::Result {
        writeln!(f, "Linear transformation matrix (sky = matrix * pixel):")?;
        let matrix = pixinsight_linear_transformation(model);
        for row in matrix.chunks(3) {
            writeln!(f, " {:+.8e}  {:+.8e}  {:+.8e}", row[0], row[1], row[2])?;
        }

        let projection = if self.result.sip.is_some() {
            "Gnomonic (TAN-SIP)"
        } else {
            "Gnomonic (TAN)"
        };
        writeln!(f, "Projection ............... {}", projection)?;
        let (origin_x, origin_y) = detection_to_pixinsight((model.crpix1, model.crpix2));
        writeln!(
            f,
            "Projection origin ........ [{:.6} {:.6}] px -> [RA: {}  Dec: {}]",
            origin_x,
            origin_y,
            fmt_ra(model.crval1),
            fmt_dec(model.crval2)
        )
    }

    fn fmt_scale_and_rotation(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let wcs = &self.result.wcs;
        let (scale_x, scale_y) = wcs.scale_arcsec();
        let scale_avg = (scale_x + scale_y) / 2.0;
        let rotation = wcs.rotation_deg();
        let linear = if self.result.sip.is_some() {
            " (linear approximation)"
        } else {
            ""
        };
        writeln!(
            f,
            "Resolution ............... {scale_avg:.3} arcsec/px{linear}"
        )?;
        let flip_suffix = if wcs.is_flipped() { " (flipped)" } else { "" };
        writeln!(
            f,
            "Rotation ................. {rotation:.3} deg{flip_suffix}{linear}"
        )?;
        writeln!(f, "Reference system ......... ICRS")
    }

    fn fmt_observation(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(time) = &self.observation_time {
            writeln!(f, "Observation start time ... {} UTC", to_the_second(time))?;
        }
        if let Some((lon, lat, alt)) = self.geodetic {
            let alt_str = alt.map(|a| format!("  {:.0} m", a)).unwrap_or_default();
            writeln!(
                f,
                "Geodetic coordinates ..... {}  {}{}",
                format_lon_deg(lon),
                format_lat_deg(lat),
                alt_str,
            )?;
        }
        let wcs = &self.result.wcs;
        if let Some(px) = wcs.pixel_um {
            let focal = wcs.effective_focal_mm().unwrap_or(0.0);
            writeln!(f, "Focal distance ........... {:.2} mm", focal)?;
            writeln!(f, "Pixel size ............... {:.2} um", px)?;
        }
        Ok(())
    }

    fn fmt_field(&self, f: &mut fmt::Formatter<'_>, model: &WcsSolution) -> fmt::Result {
        let wcs = &self.result.wcs;
        let (scale_x, scale_y) = wcs.scale_arcsec();
        let fov_x = scale_x * wcs.width as f64 / 60.0;
        let fov_y = scale_y * wcs.height as f64 / 60.0;
        writeln!(
            f,
            "Field of view ............ {} x {}",
            fmt_fov_component(fov_x),
            fmt_fov_component(fov_y)
        )?;
        writeln!(
            f,
            "Image dimensions ......... {} x {} px",
            wcs.width, wcs.height
        )?;
        let (w, h) = (wcs.width as f64, wcs.height as f64);
        write!(f, "Image center ............. ")?;
        fmt_sky_at(f, self.result, model, (w / 2.0, h / 2.0))?;
        writeln!(f, "Image bounds:")?;
        let corners = [(0.0, 0.0), (w, 0.0), (0.0, h), (w, h)];
        let labels = ["top-left", "top-right", "bottom-left", "bottom-right"];
        for (label, corner) in labels.iter().zip(corners) {
            write!(f, "   {label:<17} ")?;
            fmt_sky_at(f, self.result, model, corner)?;
        }
        Ok(())
    }

    fn fmt_sip(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Some(sip) = &self.result.sip else {
            return Ok(());
        };
        writeln!(
            f,
            "SIP distortion order ..... A={}, B={}",
            sip.a_order, sip.b_order
        )?;
        fmt_sip_coeffs(f, "A", &sip.a_coeffs)?;
        fmt_sip_coeffs(f, "B", &sip.b_coeffs)
    }

    fn fmt_residuals(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let wcs = &self.result.wcs;
        writeln!(f, "Control points ........... {}", wcs.n_stars)?;
        let Some(sip) = &self.result.sip else {
            writeln!(f, "RMS error ................ {:.3} px", wcs.rms_px)?;
            return writeln!(
                f,
                "Weighted RMS ............. {:.3} px",
                wcs.weighted_rms_px
            );
        };
        writeln!(f, "Linear RMS ............... {:.4} px", sip.linear_rms_px)?;
        writeln!(f, "SIP RMS .................. {:.4} px", sip.rms_px)?;
        writeln!(
            f,
            "Weighted SIP RMS ......... {:.4} px",
            sip.weighted_rms_px
        )?;
        let improvement = if sip.linear_rms_px > 0.0 {
            (1.0 - sip.rms_px / sip.linear_rms_px) * 100.0
        } else {
            0.0
        };
        writeln!(f, "Improvement .............. {:.1}%", improvement)
    }
}

// The sky position of `pixel`, in PixInsight's frame, and the mean residual under
// `model` of the 5 fitted pairs nearest to it.
fn fmt_sky_at(
    f: &mut fmt::Formatter<'_>,
    result: &SolveResult,
    model: &WcsSolution,
    pixel: (f64, f64),
) -> fmt::Result {
    let (x, y) = pixinsight_to_detection(pixel);
    let (ex, ey) = model.nearest_residual(x, y, 5);
    let residual = format!("ex: {ex:+.6} px  ey: {ey:+.6} px");
    match result.pixel_to_sky(x, y) {
        Ok((ra, dec)) => writeln!(f, "RA: {}  Dec: {}  {residual}", fmt_ra(ra), fmt_dec(dec)),
        Err(e) => writeln!(f, "no sky position ({e})  {residual}"),
    }
}

fn fmt_sip_coeffs(
    f: &mut fmt::Formatter<'_>,
    prefix: &str,
    coeffs: &BTreeMap<(u32, u32), f64>,
) -> fmt::Result {
    let mut keys: Vec<(u32, u32)> = coeffs.keys().copied().collect();
    keys.sort_by(|a, b| (a.0 + a.1).cmp(&(b.0 + b.1)).then(b.0.cmp(&a.0)));
    for (p, q) in keys {
        writeln!(f, "  {}_{}_{}  {:+.6e}", prefix, p, q, coeffs[&(p, q)])?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::super::sip::SipSolution;
    use super::super::solution::StarResidual;
    use super::*;

    fn sample_wcs() -> WcsSolution {
        WcsSolution {
            focal_mm: Some(1000.0),
            pixel_um: Some(3.76),
            n_stars: 50,
            rms_px: 0.3,
            weighted_rms_px: 0.2,
            residuals: vec![StarResidual {
                px_x: 1024.0,
                px_y: 768.0,
                err_x: 0.01,
                err_y: -0.02,
                err_px: 0.022,
                ra_deg: 0.0,
                dec_deg: 0.0,
            }],
            ..WcsSolution::diagonal((180.0, 0.0), (0.001, 0.001), (2048, 1536))
        }
    }

    fn sample_sip() -> SipSolution {
        let mut a_coeffs = BTreeMap::new();
        let mut b_coeffs = BTreeMap::new();
        a_coeffs.insert((2, 0), 1e-6);
        a_coeffs.insert((1, 1), 3e-6);
        b_coeffs.insert((0, 2), 2e-6);
        SipSolution {
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
            a_order: 3,
            b_order: 3,
            a_coeffs,
            b_coeffs,
            n_stars: 50,
            linear_rms_px: 1.2,
            rms_px: 0.3,
            weighted_rms_px: 0.25,
        }
    }

    #[test]
    fn display_without_sip_contains_expected_sections() {
        let wcs = sample_wcs();
        let result = SolveResult::of_fit(wcs, None, Vec::new());
        let display = SolveDisplay::new(&result);
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
        let result = SolveResult::of_fit(wcs, Some(sip), Vec::new());
        let display = SolveDisplay::new(&result);
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
        let result = SolveResult::of_fit(wcs, None, Vec::new());
        let display = SolveDisplay::new(&result).with_software("my-solver v0.1");
        let text = format!("{}", display);
        assert!(text.contains("Creation software"));
        assert!(text.contains("my-solver v0.1"));
    }

    #[test]
    fn display_with_observation_time_shows_time() {
        let wcs = sample_wcs();
        let result = SolveResult::of_fit(wcs, None, Vec::new());
        let start = UTC::from_str("2026-04-21T20:00:00").unwrap();
        let display = SolveDisplay::new(&result).with_observation_time(start);
        let text = format!("{}", display);
        assert!(text.contains("Observation start time"));
        assert!(text.contains("2026-04-21 20:00:00"));
    }

    #[test]
    fn display_with_geodetic_shows_formatted_location() {
        let wcs = sample_wcs();
        let result = SolveResult::of_fit(wcs, None, Vec::new());
        let display = SolveDisplay::new(&result).with_geodetic(-117.0, 33.0, Some(150.0));
        let text = format!("{}", display);
        assert!(text.contains("Geodetic coordinates"));
        assert!(text.contains("W")); // west longitude
        assert!(text.contains("N")); // north latitude
        assert!(text.contains("150 m"));
    }

    #[test]
    fn display_without_altitude_omits_altitude_field() {
        let wcs = sample_wcs();
        let result = SolveResult::of_fit(wcs, None, Vec::new());
        let display = SolveDisplay::new(&result).with_geodetic(-117.0, 33.0, None);
        let text = format!("{}", display);
        assert!(text.contains("Geodetic coordinates"));
        // no "m " altitude suffix
        assert!(!text.contains("  m\n"));
    }

    #[test]
    fn display_flipped_cd_adds_flipped_suffix() {
        let mut wcs = sample_wcs();
        wcs.cd1_1 = -0.001; // det < 0 -> flipped
        let result = SolveResult::of_fit(wcs, None, Vec::new());
        let display = SolveDisplay::new(&result);
        let text = format!("{}", display);
        assert!(text.contains("(flipped)"));
    }

    #[test]
    fn display_without_pixel_um_omits_focal_line() {
        let mut wcs = sample_wcs();
        wcs.pixel_um = None;
        let result = SolveResult::of_fit(wcs, None, Vec::new());
        let display = SolveDisplay::new(&result);
        let text = format!("{}", display);
        assert!(!text.contains("Focal distance"));
        assert!(!text.contains("Pixel size"));
    }

    fn square_wcs() -> WcsSolution {
        WcsSolution {
            crpix1: 512.0,
            crpix2: 512.0,
            width: 1024,
            height: 1024,
            ..sample_wcs()
        }
    }

    fn line_starting<'t>(text: &'t str, start: &str) -> &'t str {
        let line = text
            .lines()
            .map(str::trim_start)
            .find(|l| l.starts_with(start));
        line.unwrap_or_else(|| panic!("no line starting {start:?} in\n{text}"))
    }

    // PixInsight puts the centre of the first pixel at (0.5, 0.5), so the centre of a
    // 1024 px image is at 512, and the reference pixel 512 of the solver is at 512.5.
    #[test]
    fn display_gives_pixels_in_the_pixinsight_frame() {
        let wcs = square_wcs();
        let result = SolveResult::of_fit(wcs, None, Vec::new());
        let text = SolveDisplay::new(&result).to_string();
        assert!(
            line_starting(&text, "Projection origin")
                .starts_with("Projection origin ........ [512.500000 512.500000] px -> "),
            "{text}"
        );
        let offset = -(0.001 * 512.5 + 0.0 * 512.5);
        let matrix = text
            .lines()
            .skip_while(|l| !l.starts_with("Linear transformation"));
        assert_eq!(
            matrix.skip(1).take(2).collect::<Vec<_>>(),
            [
                format!(" {:+.8e}  {:+.8e}  {:+.8e}", 0.001, 0.0, offset),
                format!(" {:+.8e}  {:+.8e}  {:+.8e}", 0.0, 0.001, offset),
            ]
        );
    }

    // The SIP polynomials go with the CRVAL and CD matrix fitted with them, so the
    // matrix and origin printed beside them are those of the SIP model. The resolution
    // and rotation are those of the TAN WCS fitted alone, the linear approximation.
    #[test]
    fn display_with_sip_prints_the_origin_and_matrix_of_the_sip_model() {
        let wcs = square_wcs();
        let sip = SipSolution {
            crpix1: 512.0,
            crpix2: 512.0,
            crval1: 180.5,
            crval2: 0.25,
            cd1_1: 0.002,
            cd2_2: 0.002,
            width: 1024,
            height: 1024,
            ..sample_sip()
        };
        let result = SolveResult::of_fit(wcs, Some(sip), Vec::new());
        let text = SolveDisplay::new(&result).to_string();
        let origin = format!(
            "Projection origin ........ [512.500000 512.500000] px -> [RA: {}  Dec: {}]",
            fmt_ra(180.5),
            fmt_dec(0.25)
        );
        assert_eq!(line_starting(&text, "Projection origin"), origin);
        let offset = -(0.002 * 512.5);
        let matrix = text
            .lines()
            .skip_while(|l| !l.starts_with("Linear transformation"));
        assert_eq!(
            matrix.skip(1).take(2).collect::<Vec<_>>(),
            [
                format!(" {:+.8e}  {:+.8e}  {:+.8e}", 0.002, 0.0, offset),
                format!(" {:+.8e}  {:+.8e}  {:+.8e}", 0.0, 0.002, offset),
            ]
        );
        assert_eq!(
            line_starting(&text, "Resolution"),
            "Resolution ............... 3.600 arcsec/px (linear approximation)"
        );
        assert_eq!(
            line_starting(&text, "Rotation"),
            "Rotation ................. 180.000 deg (linear approximation)"
        );
    }

    // Beside positions that include the SIP distortion, the ex and ey columns are the
    // residuals of the SIP model. The TAN residuals here are 1 px; the pairs lie on the
    // SIP model, so their residuals are the rounding of its inverse, printed as -0.
    #[test]
    fn display_with_sip_gives_the_residuals_of_the_sip_model() {
        let sip = SipSolution {
            crpix1: 512.0,
            crpix2: 512.0,
            width: 1024,
            height: 1024,
            ..sample_sip()
        };
        let mut wcs = square_wcs();
        let model = SolveResult::of_fit(wcs.clone(), Some(sip.clone()), Vec::new());
        wcs.residuals = [
            (10.0, 10.0),
            (1000.0, 20.0),
            (30.0, 990.0),
            (1010.0, 1000.0),
        ]
        .iter()
        .map(|&(px_x, px_y)| {
            let (ra_deg, dec_deg) = model.pixel_to_sky(px_x, px_y).unwrap();
            StarResidual {
                px_x,
                px_y,
                err_x: 1.0,
                err_y: 1.0,
                err_px: libm::sqrt(2.0),
                ra_deg,
                dec_deg,
            }
        })
        .collect();
        let result = SolveResult::of_fit(wcs, Some(sip), Vec::new());
        let text = SolveDisplay::new(&result).to_string();
        let top_left = line_starting(&text, "top-left");
        assert!(
            top_left.ends_with("ex: -0.000000 px  ey: -0.000000 px"),
            "{top_left}"
        );
    }

    // The corners are the outer edges of the corner pixels, and their sky positions
    // include the SIP distortion.
    #[test]
    fn display_bounds_are_the_image_corners_with_sip() {
        let wcs = square_wcs();
        let sip = SipSolution {
            crpix1: 512.0,
            crpix2: 512.0,
            width: 1024,
            height: 1024,
            ..sample_sip()
        };
        let result = SolveResult::of_fit(wcs, Some(sip), Vec::new());
        let text = SolveDisplay::new(&result).to_string();
        let (ra, dec) = result.pixel_to_sky(-0.5, -0.5).unwrap();
        let expected = format!(
            "top-left          RA: {}  Dec: {}  ",
            fmt_ra(ra),
            fmt_dec(dec)
        );
        assert!(
            line_starting(&text, "top-left").starts_with(&expected),
            "{expected:?} in\n{text}"
        );
    }

    #[test]
    fn fmt_sip_coeffs_emits_sorted_coefficient_lines() {
        let mut coeffs = BTreeMap::new();
        coeffs.insert((2, 0), 1.5e-6);
        coeffs.insert((0, 2), 2.5e-6);

        // We can't call fmt_sip_coeffs directly with a Formatter, so test
        // indirectly via a full Display pass that exercises it.
        let mut wcs = sample_wcs();
        wcs.pixel_um = Some(3.76);
        let sip = SipSolution {
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
            a_coeffs: coeffs.clone(),
            b_coeffs: coeffs,
            n_stars: 10,
            linear_rms_px: 2.0,
            rms_px: 0.5,
            weighted_rms_px: 0.4,
        };
        let result = SolveResult::of_fit(wcs, Some(sip), Vec::new());
        let display = SolveDisplay::new(&result);
        let text = format!("{}", display);
        assert!(text.contains("A_2_0"));
        assert!(text.contains("A_0_2"));
        assert!(text.contains("B_2_0"));
        assert!(text.contains("B_0_2"));
    }

    // The report holds no clock reading of its own, so the same result gives the same
    // text.
    #[test]
    fn display_without_creation_time_has_no_creation_time_line() {
        let result = SolveResult::of_fit(sample_wcs(), None, Vec::new());
        let text = SolveDisplay::new(&result).to_string();
        assert!(!text.contains("Creation time"), "{text}");
    }

    #[test]
    fn display_prints_the_creation_time_it_is_given() {
        let result = SolveResult::of_fit(sample_wcs(), None, Vec::new());
        let made = UTC::from_str("2026-11-24T12:34:56.789Z").unwrap();
        let text = SolveDisplay::new(&result)
            .with_creation_time(made)
            .to_string();
        assert_eq!(
            line_starting(&text, "Creation time"),
            "Creation time ............ 2026-11-24 12:34:56 UTC"
        );
    }

    // The report truncates to the second. Rounding would print 03:05:00 for the first
    // time and the next day for the second, which is in the leap second at the end of
    // 2016.
    #[test]
    fn display_truncates_times_to_the_second() {
        let result = SolveResult::of_fit(sample_wcs(), None, Vec::new());
        let cases = [
            ("2024-01-15T03:04:59.9996", "2024-01-15 03:04:59"),
            ("2016-12-31T23:59:60.999999", "2016-12-31 23:59:60"),
            ("2016-12-31T23:59:59.9999", "2016-12-31 23:59:59"),
            ("2024-01-15T03:04:05", "2024-01-15 03:04:05"),
        ];
        for (time, shown) in cases {
            let made = UTC::from_str(time).unwrap();
            let text = SolveDisplay::new(&result)
                .with_creation_time(made)
                .with_observation_time(made)
                .to_string();
            assert_eq!(
                line_starting(&text, "Creation time"),
                format!("Creation time ............ {shown} UTC")
            );
            assert_eq!(
                line_starting(&text, "Observation start time"),
                format!("Observation start time ... {shown} UTC")
            );
        }
    }
}
