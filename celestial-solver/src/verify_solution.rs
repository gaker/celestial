//! The log-odds gate on the WCS that `solve` returns.
//!
//! A solution is scored as a match hypothesis (see `crate::log_odds`): its TAN WCS places
//! the brightest catalog stars about its CRVAL in the image, and the affine fit to
//! exactly the pairs the WCS was fitted on gives the variance of each placement. Each
//! fitted detection is scored against the fit to the other pairs. The refined WCS is
//! returned when it passes; otherwise the initial WCS is returned when it passes. A
//! refined solution keeps the SIP distortion fitted with its WCS for its output; the
//! initial WCS has none.
//!
//! The SIP distortion does not place the references. Leaving a pair out of the score
//! removes it from the affine fit only, and SIP terms fitted on that pair would still
//! bend the WCS onto it: with SIP of order 4, a CD matrix turned by 5 degrees passes.

use celestial_catalog::query::catalog::Catalog;
use celestial_time::julian::JulianDate;
use celestial_wcs::builder::Wcs;

use crate::detect::DetectedStar;
use crate::error::SolveError;
use crate::fit_wcs::fits_wcs::build_fits_wcs;
use crate::fit_wcs::sip::SipSolution;
use crate::fit_wcs::solution::WcsSolution;
use crate::log_odds::field::Field;
use crate::log_odds::{FitStars, Hypothesis, ACCEPT_LOG_ODDS};
use crate::match_field::quad_stars::{reference_count, reference_stars, search_radius};
use crate::match_field::StarPair;
use crate::solve::{RefineOutcome, Refinement};

// What a solution is verified against, and the fewest pairs its WCS may be fitted on.
pub(crate) struct Evidence<'a> {
    pub(crate) catalog: &'a Catalog,
    pub(crate) stars: &'a [DetectedStar],
    pub(crate) epoch: JulianDate,
    pub(crate) min_pairs: usize,
}

// A solution that passed verification: the outcome with only the pairs its WCS was
// fitted on and the SIP distortion its FITS WCS includes, that FITS WCS, and the
// log-odds of its TAN WCS.
pub(crate) struct Verified {
    pub(crate) outcome: RefineOutcome,
    pub(crate) fits_wcs: Wcs,
    pub(crate) log_odds: f64,
}

// The verified solution. When `outcome` was refined but fails verification, or its TAN
// WCS does not build, the initial WCS with `initial_pairs` is verified in its place and
// returned with `Refinement::Rejected`, as `refine_or_fallback` returns it when refine
// fails. Other failures, such as a catalog error, are operational rather than a
// rejection of the solution, and are returned at once.
pub(crate) fn verify_outcome(
    evidence: &Evidence<'_>,
    outcome: RefineOutcome,
    initial_wcs: &WcsSolution,
    initial_pairs: &[StarPair],
) -> Result<Verified, SolveError> {
    if !outcome.refinement.is_refined() {
        return verify_initial(evidence, outcome, None);
    }
    let refined = verify(evidence, outcome);
    verify_after_refined(evidence, refined, initial_wcs, initial_pairs)
}

// The refined solution when `refined`, the result of verifying it, passed, and otherwise
// as `verify_outcome` describes.
fn verify_after_refined(
    evidence: &Evidence<'_>,
    refined: Result<Result<Verified, String>, SolveError>,
    initial_wcs: &WcsSolution,
    initial_pairs: &[StarPair],
) -> Result<Verified, SolveError> {
    let rejection = match refined {
        Ok(Ok(verified)) => return Ok(verified),
        Ok(Err(reason)) => SolveError::no_match(reason),
        Err(e @ SolveError::Wcs(_)) => in_solution("refined", e),
        // A catalog error is an operational failure, not a rejection of the refined WCS,
        // so it is returned rather than answered with the initial WCS.
        Err(e) => return Err(in_solution("refined", e)),
    };
    let first = rejection.chain();
    let initial = RefineOutcome {
        wcs: initial_wcs.clone(),
        sip: None,
        pairs: initial_pairs.to_vec(),
        refinement: Refinement::Rejected(rejection),
    };
    let verified = verify_initial(evidence, initial, Some(&first))?;
    log::warn!(
        "{first}; keeping the initial solution, log-odds {:.2}",
        verified.log_odds
    );
    Ok(verified)
}

// `candidate`, the initial WCS, verified. `first` is why the refined WCS was rejected,
// when it was.
fn verify_initial(
    evidence: &Evidence<'_>,
    candidate: RefineOutcome,
    first: Option<&str>,
) -> Result<Verified, SolveError> {
    match verify(evidence, candidate) {
        Ok(Ok(verified)) => Ok(verified),
        Ok(Err(reason)) => match first {
            Some(first) => Err(not_verified(&format!("{first}; {reason}"))),
            None => Err(not_verified(&reason)),
        },
        Err(e) => {
            if let Some(first) = first {
                log::warn!("{first}; the initial solution could not be verified");
            }
            Err(in_solution("initial", e))
        }
    }
}

fn in_solution(solution: &'static str, error: SolveError) -> SolveError {
    SolveError::Verification {
        solution,
        source: Box::new(error),
    }
}

fn not_verified(rejections: &str) -> SolveError {
    SolveError::no_match(format!(
        "{rejections} (acceptance needs above {ACCEPT_LOG_ODDS:.2}, odds of 1e9)"
    ))
}

// `candidate` with only the pairs its WCS was fitted on, or why it is rejected. The
// outer error is a failure to verify at all.
fn verify(
    evidence: &Evidence<'_>,
    candidate: RefineOutcome,
) -> Result<Result<Verified, String>, SolveError> {
    let pairs = fitted_pairs(&candidate.wcs, &candidate.pairs);
    let log_odds = match passing_log_odds(evidence, &candidate, &pairs)? {
        Ok(log_odds) => log_odds,
        Err(rejection) => return Ok(Err(rejection)),
    };
    let (fits_wcs, sip, sip_error) = solution_wcs(&candidate.wcs, candidate.sip)?;
    let refinement = match sip_error {
        Some(e) => Refinement::SipLeftOut(e),
        None => candidate.refinement,
    };
    let outcome = RefineOutcome {
        wcs: candidate.wcs,
        sip,
        pairs,
        refinement,
    };
    Ok(Ok(Verified {
        outcome,
        fits_wcs,
        log_odds,
    }))
}

// The log-odds of the WCS of `candidate` fitted on `pairs` when it passes, or why it is
// rejected: too few such pairs, pairs that determine no affine model, or a log-odds at
// or below the threshold.
fn passing_log_odds(
    evidence: &Evidence<'_>,
    candidate: &RefineOutcome,
    pairs: &[StarPair],
) -> Result<Result<f64, String>, SolveError> {
    let label = if candidate.refinement.is_refined() {
        "refined"
    } else {
        "initial"
    };
    let n = pairs.len();
    if n < evidence.min_pairs {
        let need = evidence.min_pairs;
        return Ok(Err(format!(
            "{label} solution is fitted on {n} pairs (need at least {need})"
        )));
    }
    let Some(log_odds) = solution_log_odds(evidence, &candidate.wcs, pairs)? else {
        return Ok(Err(format!(
            "{label} solution: the {n} pairs the WCS was fitted on do not determine an \
             affine model"
        )));
    };
    if log_odds > ACCEPT_LOG_ODDS {
        return Ok(Ok(log_odds));
    }
    Ok(Err(format!(
        "{label} solution failed verification: log-odds {log_odds:.2} with {n} fitted pairs"
    )))
}

// The FITS WCS of `wcs` with `sip`, and that SIP. Without SIP when there is none or its
// inverse fit fails, with the error of that fit.
fn solution_wcs(
    wcs: &WcsSolution,
    sip: Option<SipSolution>,
) -> Result<(Wcs, Option<SipSolution>, Option<SolveError>), SolveError> {
    let Some(sip) = sip else {
        return Ok((build_fits_wcs(wcs, None)?, None, None));
    };
    match build_fits_wcs(wcs, Some(&sip)) {
        Ok(fits_wcs) => Ok((fits_wcs, Some(sip), None)),
        Err(e) => {
            log::warn!("SIP distortion left out of the solution: {}", e.chain());
            Ok((build_fits_wcs(wcs, None)?, None, Some(e)))
        }
    }
}

// Log-odds of `wcs`, without distortion, fitted to `fitted`. `None` when those pairs do
// not determine an affine model.
fn solution_log_odds(
    evidence: &Evidence<'_>,
    wcs: &WcsSolution,
    fitted: &[StarPair],
) -> Result<Option<f64>, SolveError> {
    let references = solution_references(evidence, wcs)?;
    let centre = (wcs.crval1, wcs.crval2);
    let field = Field::new(evidence.stars, (wcs.width, wcs.height), centre, &references);
    let hypothesis = Hypothesis::of_wcs(&field, fitted, wcs);
    Ok(hypothesis.map(|h| h.log_odds(FitStars::LeftOut, false)))
}

// The reference stars for `wcs`, by the rule `match_field` uses about its hint: the
// brightest stars of a cone about CRVAL sized from the image and the fitted scale.
fn solution_references(
    evidence: &Evidence<'_>,
    wcs: &WcsSolution,
) -> Result<Vec<(f64, f64)>, SolveError> {
    let (w, h) = (wcs.width, wcs.height);
    let scale = wcs.mean_scale_arcsec();
    let radius = search_radius(w, h, scale, None);
    let n_references = reference_count(w, h, scale, radius);
    let centre = (wcs.crval1, wcs.crval2);
    reference_stars(
        evidence.catalog,
        centre,
        radius,
        n_references,
        evidence.epoch,
    )
}

// The pairs whose detected pixel and catalog position are those of one of the WCS fit's
// residuals. Residuals carry both unchanged from their pairs, so the comparison is exact.
fn fitted_pairs(wcs: &WcsSolution, pairs: &[StarPair]) -> Vec<StarPair> {
    let mut keys: Vec<[u64; 4]> = wcs
        .residuals
        .iter()
        .map(|r| pair_key(r.px_x, r.px_y, r.ra_deg, r.dec_deg))
        .collect();
    keys.sort_unstable();
    pairs
        .iter()
        .filter(|p| {
            let key = pair_key(p.px_x, p.px_y, p.ra_deg, p.dec_deg);
            keys.binary_search(&key).is_ok()
        })
        .cloned()
        .collect()
}

fn pair_key(px_x: f64, px_y: f64, ra_deg: f64, dec_deg: f64) -> [u64; 4] {
    [px_x, px_y, ra_deg, dec_deg].map(f64::to_bits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fit_wcs::fit_wcs;
    use crate::fit_wcs::sip::fit_sip_to_tan;
    use crate::sky_scene::{epoch, star_field, Scene, CENTER, HEIGHT, WIDTH};
    use celestial_core::constants::DEG_TO_RAD;
    use celestial_wcs::error::WcsError;

    fn true_pairs(scene: &Scene) -> Vec<StarPair> {
        let pairs = scene
            .truth
            .iter()
            .map(|&(px_x, px_y, ra_deg, dec_deg)| StarPair {
                px_x,
                px_y,
                ra_deg,
                dec_deg,
                votes: 1,
                snr: 50.0,
            });
        pairs.collect()
    }

    fn wcs_of(pairs: &[StarPair]) -> WcsSolution {
        fit_wcs(pairs, WIDTH, HEIGHT, CENTER.0, CENTER.1, None, None).unwrap()
    }

    // The scene's pairs with every pixel moved 20 px along x: a WCS fitted to them
    // places every reference 20 px from its detection.
    fn shifted(pairs: &[StarPair]) -> Vec<StarPair> {
        let moved = pairs.iter().map(|p| StarPair {
            px_x: p.px_x + 20.0,
            ..p.clone()
        });
        moved.collect()
    }

    fn refined(pairs: Vec<StarPair>) -> RefineOutcome {
        RefineOutcome {
            wcs: wcs_of(&pairs),
            sip: None,
            pairs,
            refinement: Refinement::Refined,
        }
    }

    #[test]
    fn a_refined_solution_that_fails_falls_back_to_the_initial_one() {
        let scene = star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let evidence = Evidence {
            catalog: &catalog,
            stars: &stars,
            epoch: epoch(),
            min_pairs: 3,
        };
        let pairs = true_pairs(&scene);
        let initial = wcs_of(&pairs);
        let outcome = refined(shifted(&pairs));
        let verified = verify_outcome(&evidence, outcome, &initial, &pairs).unwrap();
        let kept = verified.outcome;
        assert!(matches!(
            kept.refinement,
            Refinement::Rejected(SolveError::NoMatch { .. })
        ));
        assert_eq!(
            (kept.wcs.crval1, kept.wcs.crval2),
            (initial.crval1, initial.crval2)
        );
        assert_eq!(verified.log_odds, 130.43308940175476);
    }

    #[test]
    fn a_solution_is_an_error_when_neither_fit_passes() {
        let scene = star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let evidence = Evidence {
            catalog: &catalog,
            stars: &stars,
            epoch: epoch(),
            min_pairs: 3,
        };
        let bad = shifted(&true_pairs(&scene));
        let outcome = refined(bad.clone());
        let Err(err) = verify_outcome(&evidence, outcome, &wcs_of(&bad), &bad) else {
            panic!("a WCS 20 px off every star passed verification");
        };
        assert!(matches!(err, SolveError::NoMatch { .. }));
        assert_eq!(
            err.to_string(),
            "refined solution failed verification: log-odds -1.39 with 18 fitted pairs; \
             initial solution failed verification: log-odds -1.39 with 18 fitted pairs \
             (acceptance needs above 20.72, odds of 1e9)"
        );
    }

    // The CD matrix turned 2 degrees about CRPIX: references 300 px from the centre move
    // about 10 px from their detections, while the pairs are those the WCS was fitted on.
    fn verify_turned_wcs(sip_order: Option<u32>) -> Result<Verified, SolveError> {
        let scene = star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let evidence = Evidence {
            catalog: &catalog,
            stars: &stars,
            epoch: epoch(),
            min_pairs: 3,
        };
        let pairs = true_pairs(&scene);
        let mut wcs = wcs_of(&pairs);
        let (s, c) = (libm::sin(2.0 * DEG_TO_RAD), libm::cos(2.0 * DEG_TO_RAD));
        let (a, b, d, e) = (wcs.cd1_1, wcs.cd1_2, wcs.cd2_1, wcs.cd2_2);
        (wcs.cd1_1, wcs.cd1_2) = (c * a - s * d, c * b - s * e);
        (wcs.cd2_1, wcs.cd2_2) = (s * a + c * d, s * b + c * e);
        let sip = sip_order.map(|order| fit_sip_to_tan(&pairs, &wcs, order));
        let outcome = RefineOutcome {
            wcs: wcs.clone(),
            sip,
            pairs,
            refinement: Refinement::Failed(SolveError::fit_failed("not refined")),
        };
        verify_outcome(&evidence, outcome, &wcs, &[])
    }

    const TURNED_WCS_ERROR: &str = "initial solution failed verification: log-odds 2.69 \
        with 18 fitted pairs (acceptance needs above 20.72, odds of 1e9)";

    #[test]
    fn a_wcs_that_misplaces_its_own_stars_fails() {
        let Err(err) = verify_turned_wcs(None) else {
            panic!("a WCS turned 2 degrees passed verification");
        };
        assert_eq!(err.to_string(), TURNED_WCS_ERROR);
    }

    // A SIP polynomial fitted to the pairs bends the turned WCS back onto them, so a gate
    // that placed references with it would pass this WCS.
    #[test]
    fn a_wcs_that_misplaces_its_own_stars_fails_with_sip_at_every_order() {
        for order in 2..=4 {
            let Err(err) = verify_turned_wcs(Some(order)) else {
                panic!("a WCS turned 2 degrees passed verification with SIP order {order}");
            };
            assert_eq!(err.to_string(), TURNED_WCS_ERROR, "SIP order {order}");
        }
    }

    #[test]
    fn pairs_the_wcs_was_not_fitted_on_are_neither_scored_nor_returned() {
        let scene = star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let evidence = Evidence {
            catalog: &catalog,
            stars: &stars,
            epoch: epoch(),
            min_pairs: 3,
        };
        let pairs = true_pairs(&scene);
        let wcs = wcs_of(&pairs);
        // A junk detection paired with a catalog star 1 degree away.
        let stray = StarPair {
            px_x: scene.junk[0].0,
            px_y: scene.junk[0].1,
            ra_deg: CENTER.0 + 1.0,
            dec_deg: CENTER.1,
            votes: 1,
            snr: 50.0,
        };
        let mut with_stray = pairs.clone();
        with_stray.push(stray);
        let verify = |pairs: Vec<StarPair>| {
            let outcome = RefineOutcome {
                wcs: wcs.clone(),
                sip: None,
                pairs,
                refinement: Refinement::Refined,
            };
            let verified = verify_outcome(&evidence, outcome, &wcs, &[]).unwrap();
            let kept = verified.outcome.pairs;
            let pixels: Vec<(f64, f64)> = kept.iter().map(|p| (p.px_x, p.px_y)).collect();
            (pixels, verified.log_odds)
        };
        assert_eq!(verify(with_stray), verify(pairs));
    }

    // A second catalog star paired with a fitted detection: the WCS was fitted on the
    // detection's other pair only.
    #[test]
    fn a_pair_at_a_fitted_pixel_with_another_sky_position_is_not_fitted() {
        let scene = star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let evidence = Evidence {
            catalog: &catalog,
            stars: &stars,
            epoch: epoch(),
            min_pairs: 3,
        };
        let pairs = true_pairs(&scene);
        let wcs = wcs_of(&pairs);
        let second = StarPair {
            ra_deg: pairs[0].ra_deg + 1.0,
            ..pairs[0].clone()
        };
        let mut with_second = pairs.clone();
        with_second.push(second);
        let verify = |pairs: Vec<StarPair>| {
            let outcome = RefineOutcome {
                wcs: wcs.clone(),
                sip: None,
                pairs,
                refinement: Refinement::Refined,
            };
            let verified = verify_outcome(&evidence, outcome, &wcs, &[]).unwrap();
            let kept = verified.outcome.pairs;
            let sky: Vec<(f64, f64)> = kept.iter().map(|p| (p.ra_deg, p.dec_deg)).collect();
            (sky, verified.log_odds)
        };
        assert_eq!(verify(with_second), verify(pairs));
    }

    // The FITS WCS of a SIP distortion needs the inverse polynomials AP and BP, and a
    // term of f64::MAX overflows the grid they are fitted on. The distortion is left out
    // and the result records why.
    #[test]
    fn a_sip_distortion_without_an_inverse_is_left_out_and_recorded() {
        let scene = star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let evidence = Evidence {
            catalog: &catalog,
            stars: &stars,
            epoch: epoch(),
            min_pairs: 3,
        };
        let pairs = true_pairs(&scene);
        let wcs = wcs_of(&pairs);
        let mut sip = fit_sip_to_tan(&pairs, &wcs, 2);
        sip.a_coeffs.insert((2, 0), f64::MAX);
        let outcome = RefineOutcome {
            wcs: wcs.clone(),
            sip: Some(sip),
            pairs: pairs.clone(),
            refinement: Refinement::Refined,
        };
        let verified = verify_outcome(&evidence, outcome, &wcs, &pairs).unwrap();
        assert!(verified.outcome.sip.is_none());
        let Refinement::SipLeftOut(err) = verified.outcome.refinement else {
            panic!("the SIP distortion was left out with no record");
        };
        assert_eq!(
            err.describe(),
            "FitFailed: inverse SIP polynomials of order 3: singular normal equations at col 4"
        );
    }

    // The refined WCS fails, and the initial WCS passes its log-odds on 18 pairs but is
    // held to the same floor on fitted pairs.
    #[test]
    fn the_initial_solution_needs_the_minimum_fitted_pairs() {
        let scene = star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let evidence = Evidence {
            catalog: &catalog,
            stars: &stars,
            epoch: epoch(),
            min_pairs: 19,
        };
        let pairs = true_pairs(&scene);
        let outcome = refined(shifted(&pairs));
        let Err(err) = verify_outcome(&evidence, outcome, &wcs_of(&pairs), &pairs) else {
            panic!("a solution fitted on 18 pairs passed a floor of 19");
        };
        assert_eq!(
            err.to_string(),
            "refined solution is fitted on 18 pairs (need at least 19); \
             initial solution is fitted on 18 pairs (need at least 19) \
             (acceptance needs above 20.72, odds of 1e9)"
        );
    }

    // A CRVAL2 outside -90 to 90 fails the cone search for the references, an operational
    // failure that is returned with the solution it was verifying.
    #[test]
    fn a_catalog_error_while_verifying_names_the_solution() {
        let scene = star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let evidence = Evidence {
            catalog: &catalog,
            stars: &stars,
            epoch: epoch(),
            min_pairs: 3,
        };
        let pairs = true_pairs(&scene);
        let mut off_sphere = wcs_of(&pairs);
        off_sphere.crval2 = 91.0;
        let run = |refinement: Refinement| {
            let outcome = RefineOutcome {
                wcs: off_sphere.clone(),
                sip: None,
                pairs: pairs.clone(),
                refinement,
            };
            let out = verify_outcome(&evidence, outcome, &wcs_of(&pairs), &pairs);
            out.err().map(|e| e.describe())
        };
        let got = [
            run(Refinement::Refined),
            run(Refinement::Failed(SolveError::fit_failed("not refined"))),
        ];
        let want = [
            "Verification: refined solution: Invalid input: cone search: Dec 91 is outside \
             [-90, 90] degrees",
            "Verification: initial solution: Invalid input: cone search: Dec 91 is outside \
             [-90, 90] degrees",
        ];
        assert_eq!(got.each_ref().map(|g| g.as_deref()), want.map(Some));
    }

    // A refined TAN WCS that celestial_wcs does not build. A WCS that passes the log-odds
    // has the finite CRPIX and invertible CD matrix the builder checks, so the error is
    // given directly.
    #[test]
    fn a_refined_wcs_that_does_not_build_falls_back_to_the_initial_one() {
        let scene = star_field(3, 40, 30, false);
        let catalog = Catalog::open(scene.catalog.path()).unwrap();
        let stars = scene.detections();
        let evidence = Evidence {
            catalog: &catalog,
            stars: &stars,
            epoch: epoch(),
            min_pairs: 3,
        };
        let pairs = true_pairs(&scene);
        let initial = wcs_of(&pairs);
        let unbuilt = Err(SolveError::Wcs(WcsError::NonInvertibleMatrix {
            determinant: 0.0,
        }));
        let verified = verify_after_refined(&evidence, unbuilt, &initial, &pairs).unwrap();
        let Refinement::Rejected(rejection) = &verified.outcome.refinement else {
            panic!("the refined WCS was not recorded as rejected");
        };
        assert_eq!(
            rejection.describe(),
            "Verification: refined solution: Non-invertible matrix (determinant = 0)"
        );
        assert_eq!(verified.outcome.wcs.crval1, initial.crval1);
        assert_eq!(verified.log_odds, 130.43308940175476);
    }
}
