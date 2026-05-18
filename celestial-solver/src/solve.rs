//! Builder-style entry point for plate solving.
//!
//! [`solve`] returns a [`Solver`], which accepts optional overrides and then runs the
//! full pipeline via [`Solver::run`].

use anyhow::{ensure, Result};
use celestial_catalog::query::Catalog;
use celestial_coords::ICRSPosition;
use celestial_core::math::vincenty_angular_separation;
use celestial_time::JulianDate;

use celestial_images::formats::{Image, PixelData};

use crate::detect::{find_bright_stars, DetectedStar, DetectionParams};
use crate::fit_wcs::{
    fit_sip, fit_wcs, refine_wcs,
    RefineParams, SipSolution, WcsSolution,
};
use crate::match_field::{match_field, FieldMatch, MatchParams, StarPair};

pub(crate) struct RefineOutcome {
    pub wcs: WcsSolution,
    pub pairs: Vec<StarPair>,
    pub refined: bool,
}
use crate::metadata::{metadata_from_image, ImageMetadata, MetadataOverrides};

/// Tuning parameters for the solve pipeline.
///
/// `Default::default()` is the shipped configuration: 4th-order SIP, standard detection
/// thresholds, standard matching weights. Override the individual field you want to tune
/// and leave the rest at `Default::default()`.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::SolveParams;
/// use celestial_solver::detect::DetectionParams;
///
/// let params = SolveParams {
///     detection: DetectionParams {
///         min_snr: 8.0,
///         ..DetectionParams::default()
///     },
///     sip_order: Some(5),
///     ..SolveParams::default()
/// };
/// ```
pub struct SolveParams {
    /// Star detection thresholds and filters.
    pub detection: DetectionParams,
    /// Catalog cone search radius, quad matching vote thresholds.
    pub matching: MatchParams,
    /// Outlier rejection and convergence tuning for WCS refinement.
    pub refine: RefineParams,
    /// Polynomial order for SIP distortion fit. `None` disables SIP. Default is `Some(4)`.
    pub sip_order: Option<u32>,
    /// Reject solves whose fitted plate scale differs from the input hint
    /// by more than this factor in either direction. Catches catastrophic
    /// wrong-WCS fits (matched to a different patch of sky than hinted).
    /// Default 2.0 — i.e. fitted scale must be between hint/2 and hint×2.
    pub max_scale_factor: f64,
    /// Reject solves whose fitted center is more than this many FOV
    /// diameters from the input hint. A real solve lands inside the FOV;
    /// drifting farther means the WCS converged on something unrelated.
    /// Default 1.0 (one FOV diameter).
    pub max_center_offset_fov: f64,
    /// Minimum matched pairs required after refinement. 3 is the
    /// mathematical floor for a TAN fit; below that the WCS isn't even
    /// determined. The scale and center sanity checks catch the
    /// catastrophic wrong-WCS case; this is just the floor.
    pub min_refined_pairs: usize,
}

impl Default for SolveParams {
    fn default() -> Self {
        Self {
            detection: DetectionParams::default(),
            matching: MatchParams::default(),
            refine: RefineParams::default(),
            sip_order: Some(4),
            max_scale_factor: 2.0,
            max_center_offset_fov: 1.0,
            min_refined_pairs: 3,
        }
    }
}

/// Everything the solver produced for one image.
///
/// The primary result is [`wcs`]: the refined gnomonic projection with residual
/// statistics. If SIP fitting succeeded, [`sip`] holds the polynomial coefficients.
/// The other fields expose the raw intermediate data — useful for diagnostics, debug
/// rendering with [`crate::annotate`], and custom downstream analysis.
///
/// [`wcs`]: Self::wcs
/// [`sip`]: Self::sip
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::Image::open("m31.fits")?;
/// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve(&img, &catalog).run()?;
///
/// println!("detected {} stars", result.stars.len());
/// println!("matched {} pairs", result.pairs.len());
/// println!("center: ({:.4}, {:.4})", result.wcs.crval1, result.wcs.crval2);
/// println!("RMS: {:.3} px", result.wcs.rms_px);
/// if let Some(sip) = &result.sip {
///     println!("SIP order: {}", sip.a_order);
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct SolveResult {
    /// Every star the detector found, sorted by descending flux.
    pub stars: Vec<DetectedStar>,
    /// Quad hash matches, image/catalog quad stars, and the initial matched pairs.
    pub field_match: FieldMatch,
    /// WCS fit from the initial matched pairs, before refinement.
    pub initial_wcs: WcsSolution,
    /// Refined WCS. This is what you want for pixel↔sky transforms.
    pub wcs: WcsSolution,
    /// SIP polynomial distortion fit, if the SIP pass converged.
    pub sip: Option<SipSolution>,
    /// `true` if WCS refinement succeeded. `false` means the initial WCS was kept.
    pub refined: bool,
    /// Post-refinement matched pairs. These are the stars that contributed to the
    /// final WCS fit.
    pub pairs: Vec<StarPair>,
}

/// Creates a [`Solver`] bound to `image` and `catalog`.
///
/// The returned builder reads metadata (position hint, epoch, plate scale) from the
/// image's FITS headers when [`Solver::run`] is called. Override individual fields with
/// the setter methods, or replace the whole metadata with [`Solver::with_metadata`] for
/// images that lack headers.
///
/// The `image` parameter accepts any pixel type the FITS/XISF reader produces — `u8`,
/// `u16`, `i16`, `i32`, `f32`, `f64`. No pre-conversion needed.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_images::formats::Image;
/// use celestial_catalog::query::Catalog;
///
/// let img = Image::open("m31.fits")?;
/// let catalog = Catalog::open("celestial-37M.bin")?;
///
/// let result = celestial_solver::solve(&img, &catalog).run()?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn solve<'a>(image: &'a Image, catalog: &'a Catalog) -> Solver<'a> {
    Solver {
        image,
        catalog,
        overrides: MetadataOverrides::default(),
        metadata: None,
        params: SolveParams::default(),
    }
}

/// Builder for a single solve run.
///
/// Construct with [`solve`]. Chain override methods to replace individual pieces of
/// metadata the solver would otherwise read from the image headers. Call [`Solver::run`]
/// to produce a [`SolveResult`].
///
/// The individual setters ([`Solver::hint`], [`Solver::focal_length_mm`], etc.) are
/// overrides — any field you don't set is read from the image headers via
/// [`crate::metadata::metadata_from_image`]. [`Solver::with_metadata`] is the escape
/// hatch for images without FITS headers: pass a fully-populated [`ImageMetadata`] and
/// header parsing is skipped entirely.
pub struct Solver<'a> {
    image: &'a Image,
    catalog: &'a Catalog,
    overrides: MetadataOverrides,
    metadata: Option<ImageMetadata>,
    params: SolveParams,
}

impl<'a> Solver<'a> {
    /// Overrides the position hint. Bypasses `CRVAL1/2`, `RA/DEC`, and `OBJCTRA/DEC`
    /// header lookups.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_coords::ICRSPosition;
    ///
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog)
    ///     .hint(ICRSPosition::from_degrees(83.633, 22.014)?)
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn hint(mut self, pos: ICRSPosition) -> Self {
        self.overrides.hint = Some(pos);
        self
    }

    /// Overrides the observation epoch. Bypasses `DATE-OBS` and `MJD-OBS` header lookups.
    ///
    /// The epoch is used for catalog proper-motion propagation — nearby stars have moved
    /// significantly since J2000.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_time::utc_from_calendar;
    ///
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let jd = utc_from_calendar(2026, 4, 21, 20, 0, 0.0).to_julian_date();
    /// let result = celestial_solver::solve(&img, &catalog).epoch(jd).run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn epoch(mut self, jd: JulianDate) -> Self {
        self.overrides.epoch = Some(jd);
        self
    }

    /// Overrides telescope focal length in mm. Bypasses the `FOCALLEN` header.
    ///
    /// Combined with [`Solver::pixel_size_um`] to compute the plate scale. Use this
    /// when the header reflects the bare OTA but the rig has a focal reducer, flattener,
    /// or Barlow inline.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog)
    ///     .focal_length_mm(500)    // integer literal auto-converts to f64
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn focal_length_mm<T: Into<f64>>(mut self, mm: T) -> Self {
        self.overrides.focal_mm = Some(mm.into());
        self
    }

    /// Overrides pixel size in microns. Bypasses the `XPIXSZ` header.
    ///
    /// Combined with [`Solver::focal_length_mm`] to compute the plate scale. Use this
    /// when the header reports the unbinned pixel size but the frame was binned, or
    /// when the sensor's spec differs from the header.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog)
    ///     .pixel_size_um(3.76)
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn pixel_size_um<T: Into<f64>>(mut self, um: T) -> Self {
        self.overrides.pixel_um = Some(um.into());
        self
    }

    /// Overrides the plate scale directly in arcseconds per pixel.
    ///
    /// Takes priority over [`Solver::focal_length_mm`] / [`Solver::pixel_size_um`] and
    /// skips the `206.265 * pixel_um / focal_mm` calculation. Use this when you know
    /// the plate scale from a prior solve and don't have the optics parameters.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog)
    ///     .plate_scale_arcsec(1.55)
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn plate_scale_arcsec<T: Into<f64>>(mut self, arcsec: T) -> Self {
        self.overrides.scale_arcsec = Some(arcsec.into());
        self
    }

    /// Replaces the full metadata. Skips header parsing entirely.
    ///
    /// The escape hatch for images without FITS headers — PNG, raw TIFF, anything where
    /// the position / epoch / scale live outside the file. When set, the individual
    /// override setters are ignored.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_solver::metadata::ImageMetadata;
    /// use celestial_coords::ICRSPosition;
    /// use celestial_time::utc_from_calendar;
    ///
    /// # let img = celestial_images::formats::Image::open("frame.png")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let meta = ImageMetadata {
    ///     hint: ICRSPosition::from_degrees(83.633, 22.014)?,
    ///     scale_arcsec: 1.5,
    ///     epoch: utc_from_calendar(2026, 4, 21, 20, 0, 0.0).to_julian_date(),
    ///     focal_mm: None,
    ///     pixel_um: None,
    /// };
    /// let result = celestial_solver::solve(&img, &catalog)
    ///     .with_metadata(meta)
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn with_metadata(mut self, meta: ImageMetadata) -> Self {
        self.metadata = Some(meta);
        self
    }

    /// Replaces the solver tuning parameters.
    ///
    /// See [`SolveParams`] for the individual knobs.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_solver::SolveParams;
    ///
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let params = SolveParams { sip_order: None, ..SolveParams::default() };
    /// let result = celestial_solver::solve(&img, &catalog).params(params).run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn params(mut self, params: SolveParams) -> Self {
        self.params = params;
        self
    }

    /// Runs the solve pipeline and returns the [`SolveResult`].
    ///
    /// # Errors
    ///
    /// Returns an error when header parsing fails (missing hint, missing plate scale),
    /// no stars are detected, quad matching fails to find a minimum number of matches,
    /// or the initial WCS fit fails. Refinement and SIP fitting fall back gracefully —
    /// a failed refinement keeps the initial WCS; a failed SIP pass returns `None` in
    /// [`SolveResult::sip`].
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn run(self) -> Result<SolveResult> {
        let meta = match self.metadata {
            Some(m) => m,
            None => metadata_from_image(self.image, &self.overrides)?,
        };
        run_pipeline(self.image, self.catalog, &meta, &self.params)
    }
}

macro_rules! time_stage {
    ($label:literal, $value:ident => $detail:expr, $body:expr) => {{
        let __t = std::time::Instant::now();
        let $value = $body;
        log::info!("[solve] {}: {:?} {}", $label, __t.elapsed(), $detail);
        $value
    }};
}

fn run_pipeline(
    image: &Image,
    catalog: &Catalog,
    meta: &ImageMetadata,
    params: &SolveParams,
) -> Result<SolveResult> {
    let w = image.width();
    let h = image.height();

    let stars = time_stage!("detect", v => format!("({} stars)", v.len()),
        detect_stars(image, &params.detection)?);
    let field_match = time_stage!("match", fm => format!("({} quads, {} pairs)", fm.matches.len(), fm.pairs.len()),
        match_field(&stars, image, catalog, &meta.hint, meta.scale_arcsec, meta.epoch, &params.matching)?);
    let initial_wcs = time_stage!("fit_wcs", _w => "",
        fit_initial_wcs(&field_match.pairs, w, h, meta, params)?);
    let outcome = time_stage!(
        "refine", o => format!("(refined={}, {} pairs)", o.refined, o.pairs.len()),
        refine_or_fallback(&initial_wcs, &stars, catalog, &field_match.pairs, meta.epoch, &params.refine)
    );
    ensure!(
        outcome.pairs.len() >= params.min_refined_pairs,
        "refine produced {} matched pairs (need ≥{}); field has insufficient astrometric content",
        outcome.pairs.len(),
        params.min_refined_pairs,
    );
    check_scale_sanity(&outcome.wcs, meta, params, "refined")?;
    check_center_offset_sanity(&outcome.wcs, meta, w, h, params, "refined")?;
    let sip = time_stage!("sip", s => format!("(fit={})", s.is_some()),
        fit_sip_optional(&outcome.pairs, &outcome.wcs, params.sip_order));

    Ok(SolveResult {
        stars, field_match, initial_wcs,
        wcs: outcome.wcs,
        sip,
        refined: outcome.refined,
        pairs: outcome.pairs,
    })
}

pub(crate) fn detect_stars(
    image: &Image,
    params: &DetectionParams,
) -> Result<Vec<DetectedStar>> {
    let w = image.width();
    let h = image.height();
    let stars = match &image.pixels {
        PixelData::U8(v) => find_bright_stars(v, w, h, params),
        PixelData::U16(v) => find_bright_stars(v, w, h, params),
        PixelData::I16(v) => find_bright_stars(v, w, h, params),
        PixelData::I32(v) => find_bright_stars(v, w, h, params),
        PixelData::F32(v) => find_bright_stars(v, w, h, params),
        PixelData::F64(v) => find_bright_stars(v, w, h, params),
    };
    ensure!(!stars.is_empty(), "no stars detected");
    Ok(stars)
}

pub(crate) fn fit_initial_wcs(
    pairs: &[StarPair],
    width: usize,
    height: usize,
    meta: &ImageMetadata,
    params: &SolveParams,
) -> Result<WcsSolution> {
    let wcs = fit_wcs(
        pairs, width, height,
        meta.hint.ra().degrees(), meta.hint.dec().degrees(),
        meta.focal_mm, meta.pixel_um,
    )?;
    check_scale_sanity(&wcs, meta, params, "initial")?;
    check_center_offset_sanity(&wcs, meta, width, height, params, "initial")?;
    Ok(wcs)
}

pub(crate) fn refine_or_fallback(
    initial_wcs: &WcsSolution,
    stars: &[DetectedStar],
    catalog: &Catalog,
    fallback_pairs: &[StarPair],
    epoch: JulianDate,
    refine_params: &RefineParams,
) -> RefineOutcome {
    match refine_wcs(initial_wcs, stars, catalog, epoch, refine_params) {
        Ok(r) => RefineOutcome { wcs: r.wcs, pairs: r.pairs, refined: true },
        Err(e) => {
            log::warn!("refine_wcs failed: {e:#}");
            RefineOutcome {
                wcs: initial_wcs.clone(),
                pairs: fallback_pairs.to_vec(),
                refined: false,
            }
        }
    }
}

pub(crate) fn fit_sip_optional(
    pairs: &[StarPair],
    wcs: &WcsSolution,
    sip_order: Option<u32>,
) -> Option<SipSolution> {
    sip_order.and_then(|order| match fit_sip(pairs, wcs, order) {
        Ok(s) => Some(s),
        Err(e) => {
            log::warn!("fit_sip failed: {e:#}");
            None
        }
    })
}

pub(crate) fn check_scale_sanity(
    wcs: &WcsSolution,
    meta: &ImageMetadata,
    params: &SolveParams,
    label: &str,
) -> Result<()> {
    let (sx, sy) = wcs.scale_arcsec();
    let fitted_scale = (sx + sy) * 0.5;
    let hint_scale = meta.scale_arcsec;
    let factor = (fitted_scale / hint_scale).max(hint_scale / fitted_scale);
    ensure!(
        factor <= params.max_scale_factor,
        "{label} WCS plate scale {fitted_scale:.3}\"/px diverges from hint {hint_scale:.3}\"/px by factor {factor:.2} (limit {:.2})",
        params.max_scale_factor,
    );
    Ok(())
}

pub(crate) fn check_center_offset_sanity(
    wcs: &WcsSolution,
    meta: &ImageMetadata,
    width: usize,
    height: usize,
    params: &SolveParams,
    label: &str,
) -> Result<()> {
    let sep_deg = wcs_to_hint_separation_deg(wcs, meta);
    let fov_diag_deg = fov_diagonal_deg(wcs, width, height);
    let max_offset_deg = params.max_center_offset_fov * fov_diag_deg;
    ensure!(
        sep_deg <= max_offset_deg,
        "{label} WCS center {sep_deg:.3}° from hint exceeds {:.3}° ({}× FOV diagonal {fov_diag_deg:.3}°)",
        max_offset_deg,
        params.max_center_offset_fov,
    );
    Ok(())
}

fn wcs_to_hint_separation_deg(wcs: &WcsSolution, meta: &ImageMetadata) -> f64 {
    let hint_dec = meta.hint.dec().radians();
    let wcs_dec = wcs.crval2.to_radians();
    vincenty_angular_separation(
        libm::sin(hint_dec), libm::cos(hint_dec),
        libm::sin(wcs_dec), libm::cos(wcs_dec),
        wcs.crval1.to_radians() - meta.hint.ra().radians(),
    ).to_degrees()
}

fn fov_diagonal_deg(wcs: &WcsSolution, width: usize, height: usize) -> f64 {
    let (sx, sy) = wcs.scale_arcsec();
    let fitted_scale = (sx + sy) * 0.5;
    libm::sqrt((width * width + height * height) as f64) * fitted_scale / 3600.0
}

#[cfg(test)]
mod tests {
    use super::*;

    use celestial_images::formats::PixelData;

    use crate::fit_wcs::solution::StarResidual;

    fn meta_at(ra: f64, dec: f64, scale_arcsec: f64) -> ImageMetadata {
        ImageMetadata {
            hint: ICRSPosition::from_degrees(ra, dec).unwrap(),
            scale_arcsec,
            epoch: JulianDate::new(2451545.0, 0.0),
            focal_mm: None,
            pixel_um: None,
        }
    }

    fn wcs_at(crval1: f64, crval2: f64, scale_arcsec_per_px: f64) -> WcsSolution {
        // CD matrix for a north-up, east-left frame: diag(-scale, +scale) deg/px.
        let s = scale_arcsec_per_px / 3600.0;
        WcsSolution {
            crpix1: 512.0,
            crpix2: 512.0,
            crval1,
            crval2,
            cd1_1: -s,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: s,
            width: 1024,
            height: 1024,
            focal_mm: None,
            pixel_um: None,
            n_stars: 0,
            rms_px: 0.0,
            weighted_rms_px: 0.0,
            residuals: Vec::<StarResidual>::new(),
        }
    }

    fn zero_image(w: usize, h: usize) -> Image {
        Image::new(PixelData::from(vec![0u16; w * h]), [w, h])
    }

    fn build_catalog_keep_alive() -> (Catalog, tempfile::NamedTempFile) {
        let stars = vec![
            crate::match_field::test_catalog::SynthStar {
                source_id: 1, ra: 0.0, dec: 0.0, mag: 9.0,
            },
        ];
        let file = crate::match_field::test_catalog::build(4, &stars);
        let cat = Catalog::open(file.path()).unwrap();
        (cat, file)
    }

    #[test]
    fn solve_params_defaults_match_documented_values() {
        let p = SolveParams::default();
        assert_eq!(p.sip_order, Some(4));
        assert_eq!(p.max_scale_factor, 2.0);
        assert_eq!(p.max_center_offset_fov, 1.0);
        assert_eq!(p.min_refined_pairs, 3);
    }

    #[test]
    fn solver_builder_chain_compiles_and_runs_to_no_stars() {
        // End-to-end exercise of every builder setter: hint, epoch, focal_length_mm,
        // pixel_size_um, plate_scale_arcsec, params. The image is zeroed so the
        // detector returns no stars and we exit at the `ensure!(!stars.is_empty())`.
        let (cat, _file) = build_catalog_keep_alive();
        let img = zero_image(1024, 1024);

        let result = solve(&img, &cat)
            .hint(ICRSPosition::from_degrees(10.0, 20.0).unwrap())
            .epoch(JulianDate::new(2451545.0, 0.0))
            .focal_length_mm(500)
            .pixel_size_um(3.76_f64)
            .plate_scale_arcsec(1.55_f64)
            .params(SolveParams::default())
            .run();
        let err = result.err().expect("zero image must error");
        assert!(err.to_string().contains("no stars detected"),
            "expected 'no stars detected', got: {err}");
    }

    #[test]
    fn solver_with_metadata_skips_header_parse_and_reaches_detector() {
        // with_metadata is the escape hatch — if it works, the missing-headers
        // image still reaches the detection stage instead of erroring early.
        let (cat, _file) = build_catalog_keep_alive();
        let img = zero_image(1024, 1024);
        let meta = meta_at(0.0, 0.0, 1.5);

        let err = solve(&img, &cat)
            .with_metadata(meta)
            .run()
            .err()
            .expect("zero image must error");
        assert!(err.to_string().contains("no stars detected"));
    }

    #[test]
    fn run_without_metadata_or_headers_returns_metadata_error() {
        // No overrides, no FITS headers → metadata_from_image fails before
        // we ever reach the detector.
        let (cat, _file) = build_catalog_keep_alive();
        let img = zero_image(1024, 1024);

        let err = solve(&img, &cat)
            .run()
            .err()
            .expect("no headers, no overrides must error");
        // The exact wording comes from MetadataError; just confirm it isn't
        // the detection-stage error — i.e. we failed earlier.
        assert!(!err.to_string().contains("no stars detected"));
    }

    #[test]
    fn check_scale_sanity_accepts_matching_scale() {
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = wcs_at(180.0, 30.0, 1.5);
        let params = SolveParams::default();
        assert!(check_scale_sanity(&wcs, &meta, &params, "test").is_ok());
    }

    #[test]
    fn check_center_offset_sanity_accepts_matching_center() {
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = wcs_at(180.0, 30.0, 1.5);
        let params = SolveParams::default();
        assert!(check_center_offset_sanity(&wcs, &meta, 1024, 1024, &params, "test").is_ok());
    }

    #[test]
    fn check_scale_sanity_rejects_scale_too_large() {
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = wcs_at(180.0, 30.0, 6.0);
        let params = SolveParams::default();
        let err = check_scale_sanity(&wcs, &meta, &params, "label-A").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("label-A"), "label propagated: {msg}");
        assert!(msg.contains("plate scale"), "scale error wording: {msg}");
    }

    #[test]
    fn check_scale_sanity_rejects_scale_too_small() {
        // Exercises the other half of the .max() in the factor computation.
        let meta = meta_at(180.0, 30.0, 6.0);
        let wcs = wcs_at(180.0, 30.0, 1.5);
        let params = SolveParams::default();
        let err = check_scale_sanity(&wcs, &meta, &params, "test").unwrap_err();
        assert!(err.to_string().contains("plate scale"));
    }

    #[test]
    fn check_center_offset_sanity_rejects_center_too_far() {
        // 30° offset on a ~0.43° FOV (1024px × 1.5"/px) far exceeds 1.0 FOV.
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = wcs_at(180.0 + 30.0, 30.0, 1.5);
        let params = SolveParams::default();
        let err = check_center_offset_sanity(&wcs, &meta, 1024, 1024, &params, "centerlabel")
            .unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("centerlabel"));
        assert!(msg.contains("center"), "center error wording: {msg}");
    }

    #[test]
    fn check_center_offset_sanity_respects_relaxed_thresholds() {
        let meta = meta_at(180.0, 30.0, 1.5);
        let wcs = wcs_at(180.0 + 30.0, 30.0, 1.5);
        let params = SolveParams {
            max_center_offset_fov: 1000.0,
            ..SolveParams::default()
        };
        assert!(check_center_offset_sanity(&wcs, &meta, 1024, 1024, &params, "x").is_ok());
    }

    fn gaussian_image_f32(w: usize, h: usize, cx: f64, cy: f64, peak: f32, sigma: f64, bg: f32) -> Vec<f32> {
        let mut out = vec![bg; w * h];
        for y in 0..h {
            for x in 0..w {
                let dx = x as f64 - cx;
                let dy = y as f64 - cy;
                let g = libm::exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
                out[y * w + x] = bg + peak * g as f32;
            }
        }
        out
    }

    fn image_with_one_star(pixels: PixelData, w: usize, h: usize) -> Image {
        Image::new(pixels, [w, h])
    }

    #[test]
    fn detect_stars_zero_image_errors_with_no_stars_detected() {
        let img = zero_image(128, 128);
        let err = detect_stars(&img, &DetectionParams::default())
            .err()
            .expect("zero image must error");
        assert!(err.to_string().contains("no stars detected"));
    }

    #[test]
    fn detect_stars_finds_star_in_each_pixel_type() {
        // Build a synthetic star in f32 once, cast/scale into each PixelData
        // variant, and verify the detector finds it through all six match arms.
        let w = 128;
        let h = 128;
        let base = gaussian_image_f32(w, h, 64.0, 64.0, 5000.0, 2.0, 100.0);

        // (variant_name, pixel_data) — built so each variant carries a detectable star.
        let cases: Vec<(&str, PixelData)> = vec![
            ("U8", PixelData::from(
                base.iter().map(|&v| (v / 25.0).clamp(0.0, 255.0) as u8).collect::<Vec<_>>()
            )),
            ("U16", PixelData::from(
                base.iter().map(|&v| v.clamp(0.0, 65535.0) as u16).collect::<Vec<_>>()
            )),
            ("I16", PixelData::from(
                base.iter().map(|&v| v.clamp(0.0, 32767.0) as i16).collect::<Vec<_>>()
            )),
            ("I32", PixelData::from(
                base.iter().map(|&v| v as i32).collect::<Vec<_>>()
            )),
            ("F32", PixelData::from(base.clone())),
            ("F64", PixelData::from(
                base.iter().map(|&v| v as f64).collect::<Vec<_>>()
            )),
        ];

        for (name, pixels) in cases {
            let img = image_with_one_star(pixels, w, h);
            let stars = detect_stars(&img, &DetectionParams::default())
                .unwrap_or_else(|e| panic!("{name}: detect_stars failed: {e}"));
            assert!(!stars.is_empty(), "{name}: expected ≥1 star");
            // The painted star sits at (64, 64); centroid should land within 1 px.
            let s = &stars[0];
            assert!((s.x - 64.0).abs() < 1.5, "{name}: x={} (expected ~64)", s.x);
            assert!((s.y - 64.0).abs() < 1.5, "{name}: y={} (expected ~64)", s.y);
        }
    }

    fn pairs_on_tangent_plane(center_ra: f64, center_dec: f64, scale_arcsec: f64) -> Vec<StarPair> {
        // Build 8 pairs on a synthetic TAN tangent plane so fit_wcs has data
        // it can actually fit. crpix is at (512, 512); scale = arcsec/px;
        // each grid offset in pixels maps to (xi, eta) in degrees, then we
        // tan-deproject to (ra, dec) at the given center.
        use celestial_catalog::query::tan_deproject_star;
        use celestial_core::constants::DEG_TO_RAD;
        let crpix = 512.0;
        let s_deg = scale_arcsec / 3600.0;
        let mut pairs = Vec::new();
        for &(dx, dy) in &[
            (-200.0, -200.0), (200.0, -200.0), (-200.0, 200.0), (200.0, 200.0),
            (-100.0, 50.0), (100.0, -50.0), (0.0, 0.0), (150.0, 150.0),
        ] {
            let px_x = crpix + dx;
            let px_y = crpix + dy;
            // Match the (-s, 0, 0, +s) CD matrix used by wcs_at().
            let xi_deg = -s_deg * dx;
            let eta_deg = s_deg * dy;
            let (ra, dec) = tan_deproject_star(xi_deg * DEG_TO_RAD, eta_deg * DEG_TO_RAD, center_ra, center_dec);
            pairs.push(StarPair {
                px_x, px_y, ra_deg: ra, dec_deg: dec, votes: 4, snr: 20.0,
            });
        }
        pairs
    }

    #[test]
    fn fit_initial_wcs_recovers_hint_from_clean_pairs() {
        let center_ra = 180.0;
        let center_dec = 30.0;
        let scale = 1.5;
        let pairs = pairs_on_tangent_plane(center_ra, center_dec, scale);
        let meta = meta_at(center_ra, center_dec, scale);
        let params = SolveParams::default();

        let wcs = fit_initial_wcs(&pairs, 1024, 1024, &meta, &params)
            .expect("clean pairs must produce a valid WCS");

        // CRVAL should land within a fraction of a degree of the hint.
        assert!((wcs.crval1 - center_ra).abs() < 0.1, "crval1={} (hint {})", wcs.crval1, center_ra);
        assert!((wcs.crval2 - center_dec).abs() < 0.1, "crval2={} (hint {})", wcs.crval2, center_dec);
        // And the fitted scale should track the hint scale.
        let (sx, sy) = wcs.scale_arcsec();
        let fitted = (sx + sy) * 0.5;
        assert!((fitted - scale).abs() < 0.05, "fitted scale {fitted} vs hint {scale}");
    }

    #[test]
    fn fit_initial_wcs_propagates_fit_wcs_error_on_too_few_pairs() {
        // fit_wcs requires >=3 pairs; pass 2 and check the error surfaces from
        // fit_initial_wcs rather than panicking or returning a bogus WCS.
        let pairs = vec![
            StarPair { px_x: 0.0, px_y: 0.0, ra_deg: 180.0, dec_deg: 30.0, votes: 1, snr: 10.0 },
            StarPair { px_x: 10.0, px_y: 10.0, ra_deg: 180.01, dec_deg: 30.01, votes: 1, snr: 10.0 },
        ];
        let meta = meta_at(180.0, 30.0, 1.5);
        let params = SolveParams::default();
        let err = fit_initial_wcs(&pairs, 1024, 1024, &meta, &params)
            .err()
            .expect("2 pairs must error");
        assert!(err.to_string().contains("at least 3"), "got: {err}");
    }

    #[test]
    fn refine_or_fallback_returns_initial_wcs_and_pairs_on_refine_err() {
        // Empty stars slice → refine_wcs can find no projection matches and
        // returns Err on iteration 0. The fallback path must produce a clone
        // of initial_wcs, a clone of fallback_pairs, and refined=false.
        let (cat, _file) = build_catalog_keep_alive();
        let initial = wcs_at(180.0, 30.0, 1.5);
        let fallback = vec![
            StarPair { px_x: 1.0, px_y: 2.0, ra_deg: 180.0, dec_deg: 30.0, votes: 3, snr: 50.0 },
            StarPair { px_x: 5.0, px_y: 6.0, ra_deg: 180.01, dec_deg: 30.01, votes: 2, snr: 30.0 },
        ];
        let epoch = JulianDate::new(2451545.0, 0.0);
        let refine_params = RefineParams::default();

        let outcome = refine_or_fallback(&initial, &[], &cat, &fallback, epoch, &refine_params);

        assert!(!outcome.refined, "empty stars must take fallback branch");
        assert_eq!(outcome.wcs.crval1, initial.crval1);
        assert_eq!(outcome.wcs.crval2, initial.crval2);
        assert_eq!(outcome.pairs.len(), 2);
        assert_eq!(outcome.pairs[0].px_x, 1.0);
        assert_eq!(outcome.pairs[1].px_x, 5.0);
    }

    #[test]
    fn fit_sip_optional_none_order_returns_none_without_touching_pairs() {
        let wcs = wcs_at(180.0, 30.0, 1.5);
        // Empty pairs would normally trip fit_sip's >=6 check — None bypasses it entirely.
        let out = fit_sip_optional(&[], &wcs, None);
        assert!(out.is_none());
    }

    #[test]
    fn fit_sip_optional_returns_none_on_fit_err() {
        // Some(order) with too few pairs → fit_sip errors → fit_sip_optional
        // logs and returns None. Must not propagate the error.
        let wcs = wcs_at(180.0, 30.0, 1.5);
        let pairs = vec![
            StarPair { px_x: 0.0, px_y: 0.0, ra_deg: 180.0, dec_deg: 30.0, votes: 1, snr: 10.0 },
            StarPair { px_x: 10.0, px_y: 10.0, ra_deg: 180.01, dec_deg: 30.01, votes: 1, snr: 10.0 },
        ];
        let out = fit_sip_optional(&pairs, &wcs, Some(3));
        assert!(out.is_none(), "fit_sip on 2 pairs must yield None");
    }

    #[test]
    fn fit_sip_optional_success_returns_some_solution() {
        // Use the same TAN-projected pairs that fit_initial_wcs accepts; with
        // 8 pairs at SIP order 2, fit_sip should produce a SipSolution.
        let center_ra = 180.0;
        let center_dec = 30.0;
        let scale = 1.5;
        let pairs = pairs_on_tangent_plane(center_ra, center_dec, scale);
        let meta = meta_at(center_ra, center_dec, scale);
        let params = SolveParams::default();
        let wcs = fit_initial_wcs(&pairs, 1024, 1024, &meta, &params)
            .expect("wcs setup");

        let out = fit_sip_optional(&pairs, &wcs, Some(2));
        assert!(out.is_some(), "8 clean pairs at order 2 should yield Some(SipSolution)");
    }
}
