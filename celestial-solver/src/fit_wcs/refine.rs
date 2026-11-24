mod catalog;
mod projection_match;
#[cfg(test)]
mod test_support;

use celestial_catalog::query::catalog::Catalog;
use celestial_time::julian::JulianDate;

use crate::detect::DetectedStar;
use crate::error::SolveError;
use crate::input;
use crate::log_odds::POSITION_SIGMA_PX;
use crate::match_field::StarPair;

use super::fits_wcs::build_fits_wcs;
use super::sip::SipSolution;
use super::solution::WcsSolution;
use super::tan_sip_fit::{fit_tan_sip, TanSip};
use catalog::{prep_catalog_for_refine, CatalogEntry};
use projection_match::{match_by_projection, MatchResult};

/// Tuning for the WCS refinement of the solve.
///
/// The refinement loop predicts where each catalog star falls with the current WCS and
/// its SIP distortion, pairs it with the nearest detection within `tolerance_px`, and
/// refits. It stops when the pairs repeat or after `max_iterations`.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::fit_wcs::refine::RefineParams;
///
/// let params = RefineParams {
///     tolerance_px: 2.0,
///     max_iterations: 3,
///     ..RefineParams::default()
/// };
/// ```
pub struct RefineParams {
    /// Maximum detected-to-predicted distance (pixels) to accept a pair.
    pub tolerance_px: f64,
    /// Maximum re-match iterations.
    pub max_iterations: usize,
    /// Target number of catalog stars to pull in the cone search.
    pub target_stars: usize,
}

impl Default for RefineParams {
    fn default() -> Self {
        Self {
            tolerance_px: MATCH_TOLERANCE_SIGMAS * POSITION_SIGMA_PX,
            // The five test images in tests/real-images repeat their pairs by the
            // fourth match pass.
            max_iterations: 5,
            target_stars: 2000,
        }
    }
}

// A correct pair's distance from its prediction, with per-axis error `POSITION_SIGMA_PX`,
// follows a Rayleigh distribution and exceeds 3 sigma with probability exp(-4.5), about
// 1.1%.
const MATCH_TOLERANCE_SIGMAS: f64 = 3.0;

/// Output of [`refine_wcs`].
///
/// Contains the refined [`WcsSolution`], the SIP distortion, and the matched pairs used
/// in the final fit (typically many more than the initial quad matching produced).
pub(crate) struct RefineResult {
    /// TAN WCS fitted alone to `pairs`.
    pub(crate) wcs: WcsSolution,
    // SIP distortion with its own CRVAL and CD matrix, fitted together with them to
    // `pairs`, when an order of 2 or above was chosen.
    pub(crate) sip: Option<SipSolution>,
    /// Final matched pairs after refinement.
    pub(crate) pairs: Vec<StarPair>,
}

/// Re-matches detections against the catalog and refits the WCS with its SIP distortion.
///
/// Takes the initial WCS from [`crate::fit_wcs::fit_wcs`] and uses it to predict
/// catalog positions in pixel space. Matches every detected star to its nearest
/// catalog neighbor within `tolerance_px`, fits the TAN WCS and its SIP distortion
/// together, predicts with both, and iterates until the pairs repeat. The SIP order is
/// the one up to `sip_order` with the least leave-one-out residual over the pairs, among
/// the orders with fewer terms per axis than pairs.
/// `sip_order` `None` fits the TAN WCS alone.
///
/// # Errors
///
/// Returns [`SolveError::InvalidParameter`] when `initial` has a CRPIX, CRVAL or CD
/// matrix that is not finite or a singular CD matrix, a star's position or SNR is not
/// finite, `epoch` is not finite, `params.tolerance_px` is not finite and positive, or
/// `sip_order` is not from 2 to 9; [`SolveError::Catalog`] when the cone search fails;
/// and [`SolveError::FitFailed`] or [`SolveError::Wcs`] when too few matches are found
/// or the fit fails.
pub(crate) fn refine_wcs(
    initial: &WcsSolution,
    stars: &[DetectedStar],
    catalog: &Catalog,
    epoch: JulianDate,
    params: &RefineParams,
    sip_order: Option<u32>,
) -> Result<RefineResult, SolveError> {
    check_inputs(initial, stars, epoch, params, sip_order)?;
    let (catalog_stars, _mag) =
        prep_catalog_for_refine(initial, catalog, epoch, params.target_stars)?;
    let mut model = TanSip {
        wcs: initial.clone(),
        sip: None,
    };
    let mut fitted: Option<MatchResult> = None;
    for iter in 0..params.max_iterations {
        let m = rematch(&model, stars, &catalog_stars, params.tolerance_px)?;
        if stops_at(iter, &m, fitted.as_ref())? {
            break;
        }
        model = fit_tan_sip(&m.pairs, &m.weights, &model.wcs, sip_order)?;
        log_iteration(iter, &model);
        fitted = Some(m);
    }
    let pairs = fitted.map(|f| f.pairs).ok_or_else(|| {
        SolveError::fit_failed(format!(
            "refine produced no matched pairs over {} iterations",
            params.max_iterations
        ))
    })?;
    let TanSip { wcs, sip } = model;
    Ok(RefineResult { wcs, sip, pairs })
}

fn check_inputs(
    initial: &WcsSolution,
    stars: &[DetectedStar],
    epoch: JulianDate,
    params: &RefineParams,
    sip_order: Option<u32>,
) -> Result<(), SolveError> {
    input::wcs("initial", initial)?;
    input::stars(stars)?;
    input::epoch(epoch)?;
    input::positive("tolerance_px", params.tolerance_px)?;
    input::sip_order("sip_order", sip_order)
}

// Whether refinement stops at pass `iter` with the pairs `m`: they are too few to fit,
// an error on the first pass, or they repeat the pairs `fitted` was fitted on.
fn stops_at(
    iter: usize,
    m: &MatchResult,
    fitted: Option<&MatchResult>,
) -> Result<bool, SolveError> {
    if m.pairs.len() < 3 {
        if iter == 0 {
            return Err(SolveError::fit_failed(format!(
                "refine could not extend match: only {} projection pairs within tolerance",
                m.pairs.len(),
            )));
        }
        log::debug!("refine iter {iter}: only {} pairs, stopping", m.pairs.len());
        return Ok(true);
    }
    let repeat = fitted.is_some_and(|f| f.links == m.links);
    if repeat {
        log::debug!("refine iter {iter}: the pairs repeat");
    }
    Ok(repeat)
}

// The pairs within `tolerance_px` of where `model`, with its SIP distortion, places the
// catalog stars.
fn rematch(
    model: &TanSip,
    stars: &[DetectedStar],
    catalog_stars: &[CatalogEntry],
    tolerance_px: f64,
) -> Result<MatchResult, SolveError> {
    let predictor = build_fits_wcs(&model.wcs, model.sip.as_ref())?;
    Ok(match_by_projection(
        &predictor,
        &model.wcs,
        stars,
        catalog_stars,
        tolerance_px,
    ))
}

fn log_iteration(iter: usize, model: &TanSip) {
    log::debug!(
        "refine iter {iter}: {} stars, TAN rms={:.4} px, SIP order {:?} rms={:?} px",
        model.wcs.n_stars,
        model.wcs.rms_px,
        model.sip.as_ref().map(|s| s.a_order),
        model.sip.as_ref().map(|s| s.rms_px),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    use celestial_catalog::query::catalog::Catalog;
    use celestial_catalog::query::quad::tan_project_star;
    use celestial_core::constants::RAD_TO_DEG;
    use celestial_wcs::builder::Wcs;

    use crate::detect::DetectedStar;
    use crate::sky_scene::SynthStar;

    use test_support::{det, synth_catalog};

    #[test]
    fn refine_params_defaults_match_documented_values() {
        let p = RefineParams::default();
        assert_eq!(p.tolerance_px, 3.0);
        assert_eq!(p.max_iterations, 5);
        assert_eq!(p.target_stars, 2000);
    }

    // A 3000 x 2000 px field at 3"/px in the SIP convention of Shupe et al. (2005,
    // eq. 1): x = u + A(u, v), y = v + B(u, v), with a cubic radial term that moves the
    // corners 8 px and two quadratic terms. Beyond about 1300 px from CRPIX the
    // distortion exceeds the 3 px match tolerance under the TAN WCS alone. One star sits
    // in each 100 px cell, its detection up to 0.1 px from its true pixel.
    struct DistortedField {
        truth: Wcs,
        initial: WcsSolution,
        stars: Vec<DetectedStar>,
        catalog: Catalog,
        _file: tempfile::NamedTempFile,
    }

    fn distorted_field() -> DistortedField {
        let (w, h) = (3000, 2000);
        let initial = WcsSolution::north_up((150.0, 30.0), 3.0, (w, h));
        let radial = 8.0 / libm::pow(libm::hypot(1500.0, 1000.0), 3.0);
        let mut a_coeffs = std::collections::BTreeMap::new();
        let mut b_coeffs = std::collections::BTreeMap::new();
        for (p, q, a, b) in [
            (3, 0, 1.0, 0.0),
            (1, 2, 1.0, 0.0),
            (2, 1, 0.0, 1.0),
            (0, 3, 0.0, 1.0),
        ] {
            a_coeffs.insert((p, q), a * radial);
            b_coeffs.insert((p, q), b * radial);
        }
        a_coeffs.insert((2, 0), 2.0e-7);
        b_coeffs.insert((1, 1), 1.5e-7);
        let sip = SipSolution {
            crpix1: initial.crpix1,
            crpix2: initial.crpix2,
            crval1: initial.crval1,
            crval2: initial.crval2,
            cd1_1: initial.cd1_1,
            cd1_2: initial.cd1_2,
            cd2_1: initial.cd2_1,
            cd2_2: initial.cd2_2,
            width: w,
            height: h,
            a_order: 3,
            b_order: 3,
            a_coeffs,
            b_coeffs,
            n_stars: 0,
            linear_rms_px: 0.0,
            rms_px: 0.0,
            weighted_rms_px: 0.0,
        };
        let truth = build_fits_wcs(&initial, Some(&sip)).unwrap();
        let mut next = crate::sky_scene::uniform(11);
        let mut stars = Vec::new();
        let mut synth = Vec::new();
        for k in 0..600 {
            let x = f64::from(k % 30) * 100.0 + 10.0 + 80.0 * next();
            let y = f64::from(k / 30) * 100.0 + 10.0 + 80.0 * next();
            let (ra, dec) = crate::pixel_frame::detection_to_sky(&truth, x, y).unwrap();
            let (dx, dy) = (0.2 * next() - 0.1, 0.2 * next() - 0.1);
            stars.push(det(x + dx, y + dy, 50.0));
            let mag = 8.0 + f64::from(k) as f32 * 0.01;
            synth.push(SynthStar {
                source_id: i64::from(k) + 1,
                ra,
                dec,
                mag,
            });
        }
        let (catalog, _file) = synth_catalog(synth);
        DistortedField {
            truth,
            initial,
            stars,
            catalog,
            _file,
        }
    }

    // An 8000 x 600 px field at 3"/px with a cubic radial distortion of 8 px at the
    // corners and a sixth-order term of 4 px at the ends of the long axis, one star in
    // each 100 px cell. Sky positions are the TAN deprojection of the SIP-corrected
    // pixel, so the fixture fits no inverse polynomials.
    fn elongated_field() -> DistortedField {
        let (w, h) = (8000, 600);
        let initial = WcsSolution::north_up((150.0, 30.0), 3.0, (w, h));
        let radial = 8.0 / libm::pow(libm::hypot(4000.0, 300.0), 3.0);
        let a_6_0 = 4.0 / libm::pow(4000.0, 6.0);
        let distort = |u: f64, v: f64| {
            let r2 = u * u + v * v;
            (
                u + radial * u * r2 + a_6_0 * libm::pow(u, 6.0),
                v + radial * v * r2,
            )
        };
        let mut next = crate::sky_scene::uniform(13);
        let mut stars = Vec::new();
        let mut synth = Vec::new();
        for k in 0..480 {
            let x = f64::from(k % 80) * 100.0 + 10.0 + 80.0 * next();
            let y = f64::from(k / 80) * 100.0 + 10.0 + 80.0 * next();
            let (u, v) = distort(x - initial.crpix1, y - initial.crpix2);
            let (ra, dec) = initial.pixel_to_sky(u + initial.crpix1, v + initial.crpix2);
            let (dx, dy) = (0.2 * next() - 0.1, 0.2 * next() - 0.1);
            stars.push(det(x + dx, y + dy, 50.0));
            let mag = 8.0 + f64::from(k) as f32 * 0.01;
            synth.push(SynthStar {
                source_id: i64::from(k) + 1,
                ra,
                dec,
                mag,
            });
        }
        let (catalog, _file) = synth_catalog(synth);
        DistortedField {
            truth: build_fits_wcs(&initial, None).unwrap(),
            initial,
            stars,
            catalog,
            _file,
        }
    }

    // The long-axis monomials of a field 13 times wider than tall are scaled by the long
    // extent; fits of order 6 and its order 7 inverse still determine every term.
    #[test]
    fn refine_fits_a_high_order_sip_on_an_elongated_field() {
        let field = elongated_field();
        let epoch = JulianDate::new(2451545.0, 0.0);
        let params = RefineParams::default();
        let r = refine_wcs(
            &field.initial,
            &field.stars,
            &field.catalog,
            epoch,
            &params,
            Some(6),
        )
        .unwrap();
        let sip = r.sip.as_ref().map(|s| (s.a_order, s.n_stars));
        assert_eq!((sip, r.pairs.len()), (Some((6, 480)), 480));
    }

    // Separation in arcsec between the sky positions of `pixel` under `a` and `b`.
    fn separation_arcsec(a: &Wcs, b: &Wcs, pixel: (f64, f64)) -> f64 {
        let (ra_a, dec_a) = crate::pixel_frame::detection_to_sky(a, pixel.0, pixel.1).unwrap();
        let (ra_b, dec_b) = crate::pixel_frame::detection_to_sky(b, pixel.0, pixel.1).unwrap();
        let (xi, eta) = tan_project_star(ra_b, dec_b, ra_a, dec_a).unwrap();
        libm::hypot(xi, eta) * RAD_TO_DEG * 3600.0
    }

    fn nearest(points: impl Iterator<Item = (f64, f64)>, to: (f64, f64)) -> f64 {
        points
            .map(|p| libm::hypot(p.0 - to.0, p.1 - to.1))
            .fold(f64::INFINITY, f64::min)
    }

    // The refined pairs reach every corner: the nearest fitted pair to each corner is
    // its nearest detection. The fitted model then places each corner closer to the
    // truth than the standard deviation of one detection offset along one axis,
    // 0.2 / sqrt(12) px for offsets uniform in [-0.1, 0.1] px, or 0.17" at 3"/px.
    #[test]
    fn refine_matches_the_corners_of_a_distorted_field_and_recovers_its_sip() {
        let field = distorted_field();
        let epoch = JulianDate::new(2451545.0, 0.0);
        let params = RefineParams::default();
        let r = refine_wcs(
            &field.initial,
            &field.stars,
            &field.catalog,
            epoch,
            &params,
            Some(4),
        )
        .unwrap();
        let model = build_fits_wcs(&r.wcs, r.sip.as_ref()).unwrap();
        let corners = [(0.0, 0.0), (2999.0, 0.0), (0.0, 1999.0), (2999.0, 1999.0)];
        let fitted: Vec<f64> = corners
            .iter()
            .map(|&c| nearest(r.pairs.iter().map(|p| (p.px_x, p.px_y)), c))
            .collect();
        let detected: Vec<f64> = corners
            .iter()
            .map(|&c| nearest(field.stars.iter().map(|s| (s.x, s.y)), c))
            .collect();
        assert_eq!(fitted, detected);
        let errors: Vec<f64> = corners
            .iter()
            .map(|&c| separation_arcsec(&field.truth, &model, c))
            .collect();
        // The detections are up to 0.1 px from their true pixels, a standard deviation
        // of 0.2 / sqrt(12) px, 0.173" at 3"/px. Every corner error is below that.
        assert_eq!(
            errors,
            [
                0.06857467454275454,
                0.06615867652590171,
                0.1034837574089729,
                0.08963094933175193
            ]
        );
        let sip = r.sip.as_ref().map(|s| (s.a_order, s.n_stars));
        assert_eq!((sip, r.pairs.len()), (Some((3, 600)), 600));
    }

    // A NaN in the initial WCS reached the catalog cone search, whose error named
    // neither the WCS nor the field.
    #[test]
    fn refine_wcs_rejects_non_finite_inputs() {
        let star = SynthStar {
            source_id: 1,
            ra: 180.0,
            dec: 30.0,
            mag: 9.0,
        };
        let (catalog, _file) = synth_catalog(vec![star]);
        let good = WcsSolution::north_up((180.0, 30.0), 1.5, (1024, 1024));
        let stars = vec![det(512.0, 512.0, 20.0)];
        let epoch = JulianDate::new(2451545.0, 0.0);
        let refine = |initial: &WcsSolution,
                      stars: &[DetectedStar],
                      epoch: JulianDate,
                      tolerance_px: f64,
                      sip_order: Option<u32>| {
            let params = RefineParams {
                tolerance_px,
                ..RefineParams::default()
            };
            let out = refine_wcs(initial, stars, &catalog, epoch, &params, sip_order);
            out.err().map(|e| e.describe())
        };
        let nan_cd = WcsSolution {
            cd1_1: f64::NAN,
            ..good.clone()
        };
        let nan_crpix = WcsSolution {
            crpix2: f64::NAN,
            ..good.clone()
        };
        let zero_cd = WcsSolution {
            cd1_1: 0.0,
            cd2_2: 0.0,
            ..good.clone()
        };
        let nan_snr = vec![det(512.0, 512.0, f64::NAN)];
        let tol = 3.0;
        let cases = [
            (
                refine(&nan_cd, &stars, epoch, tol, Some(4)),
                "InvalidParameter: invalid parameter initial: cd1_1 must be finite, got NaN",
            ),
            (
                refine(&nan_crpix, &stars, epoch, tol, Some(4)),
                "InvalidParameter: invalid parameter initial: crpix2 must be finite, got NaN",
            ),
            (
                refine(&zero_cd, &stars, epoch, tol, Some(4)),
                "InvalidParameter: invalid parameter initial: the CD matrix must be invertible, got determinant 0",
            ),
            (
                refine(&good, &stars, epoch, f64::NAN, Some(4)),
                "InvalidParameter: invalid parameter tolerance_px: must be finite and positive, got NaN",
            ),
            (
                refine(&good, &nan_snr, epoch, tol, Some(4)),
                "InvalidParameter: invalid parameter stars: star 0 is not finite: pixel (512, 512), snr NaN",
            ),
            (
                refine(&good, &stars, JulianDate::new(f64::NAN, 0.0), tol, Some(4)),
                "InvalidParameter: invalid parameter epoch: must be finite, got (NaN, 0)",
            ),
            (
                refine(&good, &stars, epoch, tol, Some(10)),
                "InvalidParameter: invalid parameter sip_order: must be from 2 to 9, got 10",
            ),
        ];
        let (got, want): (Vec<_>, Vec<_>) = cases
            .iter()
            .map(|(err, expected)| (err.as_deref(), Some(*expected)))
            .unzip();
        assert_eq!(got, want);
    }
}
