//! Builder-style entry point for plate solving.
//!
//! [`solve`] returns a [`Solver`], which accepts optional overrides and then runs the
//! full pipeline via [`Solver::run`].

mod builder;
mod pipeline;
mod sanity;
#[cfg(test)]
mod test_support;

use celestial_catalog::query::catalog::Catalog;
use celestial_wcs::builder::Wcs;

use celestial_images::formats::image::Image;

use crate::detect::{DetectedStar, DetectionParams};
use crate::error::SolveError;
use crate::fit_wcs::refine::RefineParams;
use crate::fit_wcs::sip::SipSolution;
use crate::fit_wcs::solution::WcsSolution;
use crate::match_field::{FieldMatch, MatchParams, StarPair};
use crate::metadata::image_metadata::{ImageMetadata, MetadataOverrides};

pub(crate) struct RefineOutcome {
    pub(crate) wcs: WcsSolution,
    pub(crate) sip: Option<SipSolution>,
    pub(crate) pairs: Vec<StarPair>,
    pub(crate) refinement: Refinement,
}

// Which WCS a solve returned, and the error behind any fallback. Only `Refined` returns
// everything the refinement fitted.
#[derive(Debug)]
pub enum Refinement {
    // The refined WCS, with the SIP distortion its fit chose.
    Refined,
    // The refined WCS without its SIP distortion: the inverse polynomials AP and BP
    // could not be fitted to it, so its FITS WCS could not map sky to pixels.
    SipLeftOut(SolveError),
    // The initial WCS: the refinement failed with this error.
    Failed(SolveError),
    // The initial WCS: the refined WCS failed verification with this error.
    Rejected(SolveError),
}

impl Refinement {
    // Whether the returned WCS is the refined one.
    pub(crate) fn is_refined(&self) -> bool {
        matches!(self, Self::Refined | Self::SipLeftOut(_))
    }
}

/// Tuning parameters for the solve pipeline.
///
/// `Default::default()` is the shipped configuration: 4th-order SIP, standard detection
/// thresholds, standard matching weights. Override the individual field you want to tune
/// and leave the rest at `Default::default()`.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::solve::SolveParams;
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
    /// Highest polynomial order the SIP distortion fit may choose, from 2 to 9. The
    /// refinement chooses the order, or no SIP, with the least leave-one-out residual
    /// over its pairs. `None` disables SIP. Default is `Some(4)`. Another order makes
    /// `run` return [`SolveError::InvalidParameter`].
    pub sip_order: Option<u32>,
    /// Reject solves whose initial WCS plate scale differs from the input hint
    /// by more than this factor in either direction. The check runs before
    /// refinement, whose catalog search radius grows with the fitted scale.
    /// Default 2.0: the initial scale must be between hint/2 and hint×2. At least 1;
    /// `f64::INFINITY` turns the check off.
    pub max_scale_factor: f64,
    /// Reject solves whose initial WCS center is more than this many FOV
    /// diagonals from the input hint. The check runs before refinement.
    /// Default 1.0 (one FOV diagonal). Above 0; `f64::INFINITY` turns the check off.
    pub max_center_offset_fov: f64,
    /// Minimum pairs the returned WCS must be fitted on. 3 is the
    /// mathematical floor for a TAN fit; below that the WCS isn't even
    /// determined. Whether the match is real is decided by its log-odds
    /// (see [`SolveResult::log_odds`]); this is just the floor.
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
/// statistics. If a SIP order was chosen, [`sip`] holds the polynomial coefficients with
/// their own CRVAL and CD matrix.
/// The other fields expose the raw intermediate data — useful for diagnostics, debug
/// rendering with [`crate::annotate`], and custom downstream analysis.
///
/// [`wcs`]: Self::wcs
/// [`sip`]: Self::sip
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::image::Image::open("m31.fits")?;
/// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
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
    /// Refined TAN WCS fitted alone to `pairs`. [`SolveResult::pixel_to_sky`] and
    /// [`SolveResult::sky_to_pixel`] apply `sip` as well.
    pub wcs: WcsSolution,
    /// SIP distortion with its own CRVAL and CD matrix, fitted together with them to
    /// `pairs`. `None` when SIP is disabled, the refinement chose order 1, the
    /// initial WCS was kept, or the distortion was left out (see `refinement`).
    pub sip: Option<SipSolution>,
    // Whether `wcs` is the refined WCS or the initial one, and the error that made the
    // solve keep the initial WCS or leave the SIP distortion out.
    pub refinement: Refinement,
    /// The matched pairs `wcs` was fitted on. These are the stars that contributed to
    /// the final WCS fit.
    pub pairs: Vec<StarPair>,
    // The solution as a FITS WCS, in FITS pixel coordinates (the centre of the first
    // pixel is (1, 1)): `sip` with its CRVAL and CD matrix and the inverse polynomials
    // AP and BP fitted to it, or `wcs` when there is no SIP. `pixel_to_sky` and
    // `sky_to_pixel` evaluate it in the detection pixels of `stars` and `pairs`, and
    // `to_fits_keywords` writes it.
    pub fits_wcs: Wcs,
    // Natural-log Bayes factor that the returned WCS is real rather than chance (Lang
    // et al. 2010, sec. 2.4); every returned solution exceeds ln(1e9) = 20.72. `wcs`
    // places the catalog stars; each detection in `pairs` is scored against the affine
    // fit to the other pairs.
    pub log_odds: f64,
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
/// use celestial_images::formats::image::Image;
/// use celestial_catalog::query::catalog::Catalog;
///
/// let img = Image::open("m31.fits")?;
/// let catalog = Catalog::open("celestial-37M.bin")?;
///
/// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
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
/// [`crate::metadata::image_metadata::metadata_from_image`]. [`Solver::with_metadata`] is the escape
/// hatch for images without FITS headers: pass a fully-populated [`ImageMetadata`] and
/// header parsing is skipped entirely.
pub struct Solver<'a> {
    image: &'a Image,
    catalog: &'a Catalog,
    overrides: MetadataOverrides,
    metadata: Option<ImageMetadata>,
    params: SolveParams,
}

#[cfg(test)]
impl SolveResult {
    // The result of `wcs` with `sip`, fitted on `pairs`, with no detections or quad
    // matches, for tests of what a result writes and prints.
    pub(crate) fn of_fit(wcs: WcsSolution, sip: Option<SipSolution>, pairs: Vec<StarPair>) -> Self {
        let field_match = FieldMatch {
            image_stars: Vec::new(),
            catalog_stars: Vec::new(),
            matches: Vec::new(),
            pairs: pairs.clone(),
        };
        let fits_wcs = crate::fit_wcs::fits_wcs::build_fits_wcs(&wcs, sip.as_ref())
            .expect("the WCS and SIP of a test result build a FITS WCS");
        SolveResult {
            stars: Vec::new(),
            field_match,
            initial_wcs: wcs.clone(),
            wcs,
            sip,
            refinement: Refinement::Refined,
            pairs,
            fits_wcs,
            log_odds: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::pipeline::fit_initial_wcs;
    use super::test_support::{
        build_catalog_keep_alive, meta_at, pairs_on_tangent_plane, zero_image,
    };

    #[test]
    fn solve_params_defaults_match_documented_values() {
        let p = SolveParams::default();
        assert_eq!(p.sip_order, Some(4));
        assert_eq!(p.max_scale_factor, 2.0);
        assert_eq!(p.max_center_offset_fov, 1.0);
        assert_eq!(p.min_refined_pairs, 3);
    }

    #[test]
    fn solve_rejects_unrelated_star_fields() {
        for seed in 2..=7 {
            let scene = crate::sky_scene::unrelated_field(seed, 100);
            let catalog = Catalog::open(scene.catalog.path()).unwrap();
            let outcome = solve(&scene.image(), &catalog)
                .with_metadata(scene.metadata())
                .run();
            let Err(err) = outcome else {
                panic!("seed {seed}: stars unrelated to the catalog solved");
            };
            assert!(
                err.to_string()
                    .starts_with("no quad match passed verification: det(CD) > 0: "),
                "seed {seed}: {err}"
            );
        }
    }

    #[test]
    fn solve_finds_sparse_field_among_unrelated_stars() {
        for (flip, log_odds) in [(false, 128.71981578078632), (true, 130.10584361683678)] {
            let scene = crate::sky_scene::star_field(3, 40, 30, flip);
            let catalog = Catalog::open(scene.catalog.path()).unwrap();
            let result = solve(&scene.image(), &catalog)
                .with_metadata(scene.metadata())
                .run()
                .unwrap_or_else(|e| panic!("flip={flip}: {e:#}"));
            scene.assert_true_pairs(&result.pairs);
            let w = &result.wcs;
            assert_eq!(w.cd1_1 * w.cd2_2 - w.cd1_2 * w.cd2_1 < 0.0, flip);
            assert_eq!(result.log_odds, log_odds, "flip={flip}");
            assert!(matches!(result.refinement, Refinement::Refined));
        }
    }

    // A caller tells no stars, bad metadata, no match and a failed fit apart by variant.
    #[test]
    fn solve_errors_name_their_cause() {
        let (cat, _file) = build_catalog_keep_alive();
        let blank = zero_image(1024, 1024);
        let no_stars = solve(&blank, &cat)
            .with_metadata(meta_at(0.0, 0.0, 1.5))
            .run();
        let no_headers = solve(&blank, &cat).run();
        let scene = crate::sky_scene::unrelated_field(2, 100);
        let unrelated = Catalog::open(scene.catalog.path()).unwrap();
        let no_match = solve(&scene.image(), &unrelated)
            .with_metadata(scene.metadata())
            .run();
        let pairs = &pairs_on_tangent_plane(180.0, 30.0, 1.5)[..2];
        let meta = meta_at(180.0, 30.0, 1.5);
        let too_few = fit_initial_wcs(pairs, 1024, 1024, &meta, &SolveParams::default());
        let describe = |e: SolveError| e.describe();
        assert_eq!(
            no_stars.err().map(describe).as_deref(),
            Some("NoStars: no stars detected")
        );
        assert_eq!(
            no_headers.err().map(describe).as_deref(),
            Some("Metadata: no position hint: need CRVAL1/2, RA/DEC, or OBJCTRA/OBJCTDEC")
        );
        let no_match = no_match.err().map(describe).unwrap_or_default();
        assert!(
            no_match.starts_with("NoMatch: no quad match passed verification: det(CD) > 0: "),
            "{no_match}"
        );
        assert_eq!(
            too_few.err().map(describe).as_deref(),
            Some("FitFailed: need at least 3 matched pairs, got 2")
        );
    }
}
