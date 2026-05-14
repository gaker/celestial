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

fn run_pipeline(
    image: &Image,
    catalog: &Catalog,
    meta: &ImageMetadata,
    params: &SolveParams,
) -> Result<SolveResult> {
    let w = image.width();
    let h = image.height();

    let t = std::time::Instant::now();
    let stars = match &image.pixels {
        PixelData::U8(v) => find_bright_stars(v, w, h, &params.detection),
        PixelData::U16(v) => find_bright_stars(v, w, h, &params.detection),
        PixelData::I16(v) => find_bright_stars(v, w, h, &params.detection),
        PixelData::I32(v) => find_bright_stars(v, w, h, &params.detection),
        PixelData::F32(v) => find_bright_stars(v, w, h, &params.detection),
        PixelData::F64(v) => find_bright_stars(v, w, h, &params.detection),
    };
    ensure!(!stars.is_empty(), "no stars detected");
    log::info!("[solve] detect: {:?} ({} stars)", t.elapsed(), stars.len());

    let t = std::time::Instant::now();
    let field_match = match_field(
        &stars, image, catalog, &meta.hint,
        meta.scale_arcsec, meta.epoch, &params.matching,
    )?;
    log::info!("[solve] match: {:?} ({} quads, {} pairs)",
        t.elapsed(), field_match.matches.len(), field_match.pairs.len());

    let t = std::time::Instant::now();
    let initial_wcs = fit_wcs(
        &field_match.pairs, w, h,
        meta.hint.ra().degrees(), meta.hint.dec().degrees(),
        meta.focal_mm, meta.pixel_um,
    )?;
    log::info!("[solve] fit_wcs: {:?}", t.elapsed());
    check_wcs_sanity(&initial_wcs, meta, w, h, params, "initial")?;

    let t = std::time::Instant::now();
    let (wcs, refined_pairs, refined) = match refine_wcs(
        &initial_wcs, &stars, catalog, meta.epoch, &params.refine,
    ) {
        Ok(r) => (r.wcs, r.pairs, true),
        Err(e) => {
            log::warn!("refine_wcs failed: {e:#}");
            (initial_wcs.clone(), field_match.pairs.clone(), false)
        }
    };
    log::info!("[solve] refine: {:?} (refined={refined}, {} pairs)",
        t.elapsed(), refined_pairs.len());
    ensure!(
        refined_pairs.len() >= params.min_refined_pairs,
        "refine produced {} matched pairs (need ≥{}); field has insufficient astrometric content",
        refined_pairs.len(),
        params.min_refined_pairs,
    );
    check_wcs_sanity(&wcs, meta, w, h, params, "refined")?;

    let t = std::time::Instant::now();
    let sip = params.sip_order.and_then(|order| {
        match fit_sip(&refined_pairs, &wcs, order) {
            Ok(s) => Some(s),
            Err(e) => {
                log::warn!("fit_sip failed: {e:#}");
                None
            }
        }
    });
    log::info!("[solve] sip: {:?} (fit={})", t.elapsed(), sip.is_some());

    Ok(SolveResult {
        stars,
        field_match,
        initial_wcs,
        wcs,
        sip,
        refined,
        pairs: refined_pairs,
    })
}

/// Validates a fitted WCS against the input hints. Fails the solve if the
/// fitted plate scale or field center has drifted too far from what the
/// caller asked for — both are strong signals the fit converged on noise
/// or a completely different patch of sky.
fn check_wcs_sanity(
    wcs: &WcsSolution,
    meta: &ImageMetadata,
    width: usize,
    height: usize,
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

    let hint_ra = meta.hint.ra().radians();
    let hint_dec = meta.hint.dec().radians();
    let wcs_ra = wcs.crval1.to_radians();
    let wcs_dec = wcs.crval2.to_radians();
    let sep_rad = vincenty_angular_separation(
        libm::sin(hint_dec), libm::cos(hint_dec),
        libm::sin(wcs_dec), libm::cos(wcs_dec),
        wcs_ra - hint_ra,
    );
    let sep_deg = sep_rad.to_degrees();
    let fov_diag_deg = libm::sqrt((width * width + height * height) as f64)
        * fitted_scale / 3600.0;
    let max_offset_deg = params.max_center_offset_fov * fov_diag_deg;
    ensure!(
        sep_deg <= max_offset_deg,
        "{label} WCS center {sep_deg:.3}° from hint exceeds {:.3}° ({}× FOV diagonal {fov_diag_deg:.3}°)",
        max_offset_deg,
        params.max_center_offset_fov,
    );
    Ok(())
}
