//! WCS fitting: gnomonic projection, CD matrix, SIP distortion, residual refinement.
//!
//! The solve fits a TAN (gnomonic) WCS, a [`solution::WcsSolution`], to the matched star
//! pairs, then re-matches against the catalog and refits the TAN WCS with its SIP
//! distortion, a [`sip::SipSolution`].

mod cd_matrix;
pub mod display;
pub(crate) mod fits_wcs;
mod median_clip;
mod normal_equations;
mod polynomial_fit;
mod ransac;
pub mod refine;
mod sexagesimal;
pub mod sip;
mod sip_inverse;
pub mod solution;
mod tan_sip_fit;
mod tangent_projection;

use celestial_catalog::query::quad::tan_deproject_star;
use celestial_core::constants::{DEG_TO_RAD, PI};

use crate::detect::image_centre;
use crate::error::SolveError;
use crate::input;
use crate::log_odds::POSITION_SIGMA_PX;
use crate::match_field::StarPair;
use solution::WcsSolution;

use cd_matrix::{
    compute_pixel_errors, compute_residuals, fit_cd_weighted, median_clip_loop, rms_error,
    weighted_rms_error, CdParams, WcsParams,
};
use ransac::ransac_fit;
use tangent_projection::{extract_kept_pairs, project_pairs_to_tangent, TangentProjection};

/// Fits a TAN (gnomonic) WCS from matched star pairs.
///
/// RANSAC inlier selection, then weighted least-squares CD-matrix fit, then
/// sigma-clipping loop, then a final weighted fit. The CRVAL is refined from the
/// initial `crval_ra_deg` / `crval_dec_deg` prior using the tangent-plane offsets.
/// Each pair is weighted by the inverse of its centroid variance, which falls with its
/// SNR towards a 1 px floor.
///
/// # Errors
///
/// Returns [`SolveError::InvalidParameter`] when a pair or the prior CRVAL is not
/// finite, the prior declination is outside -90 to 90, or a set focal length or pixel
/// size is not finite and positive; and [`SolveError::FitFailed`] when fewer than 3
/// pairs are supplied, the tangent-plane projection fails on too many of them, RANSAC
/// fails to find a consistent inlier set, or the fit is singular.
pub(crate) fn fit_wcs(
    pairs: &[StarPair],
    width: usize,
    height: usize,
    crval_ra_deg: f64,
    crval_dec_deg: f64,
    focal_mm: Option<f64>,
    pixel_um: Option<f64>,
) -> Result<WcsSolution, SolveError> {
    check_inputs(pairs, (crval_ra_deg, crval_dec_deg), focal_mm, pixel_um)?;
    check_pair_count(pairs.len())?;
    let (crpix1, crpix2) = image_centre(width, height);
    let (cd, (crval1, crval2), (kept_pairs, kept_weights)) =
        fit_stable_cd(pairs, (crval_ra_deg, crval_dec_deg), (crpix1, crpix2))?;

    let wcs_params = WcsParams {
        crpix1,
        crpix2,
        crval_ra: crval1,
        crval_dec: crval2,
        cd1_1: cd.cd1_1,
        cd1_2: cd.cd1_2,
        cd2_1: cd.cd2_1,
        cd2_2: cd.cd2_2,
    };
    let residuals = compute_residuals(&kept_pairs, wcs_params);
    let sizes = || residuals.iter().map(|r| r.err_px);
    Ok(WcsSolution {
        crpix1,
        crpix2,
        crval1,
        crval2,
        cd1_1: cd.cd1_1,
        cd1_2: cd.cd1_2,
        cd2_1: cd.cd2_1,
        cd2_2: cd.cd2_2,
        width,
        height,
        focal_mm,
        pixel_um,
        n_stars: residuals.len(),
        rms_px: rms_error(sizes()),
        weighted_rms_px: weighted_rms_error(sizes(), &kept_weights),
        residuals,
    })
}

// The CD matrix and CRVAL fitted to the inliers among `pairs`, projected about `prior`
// and then refitted until the inliers repeat, with the pairs and weights of that fit.
fn fit_stable_cd<'a>(
    pairs: &'a [StarPair],
    prior: (f64, f64),
    crpix: (f64, f64),
) -> Result<(CdParams, (f64, f64), Kept<'a>), SolveError> {
    let weights: Vec<f64> = pairs.iter().map(|p| pair_weight(p.snr)).collect();
    let tp = project_pairs_to_tangent(pairs, &weights, prior.0, prior.1, crpix.0, crpix.1);
    check_projected_count(tp.rows.len())?;
    let (prior_cd, mask) = fit_inlier_cd(&tp)?;
    let kept = extract_kept_pairs(&tp, &mask);
    let centre = fitted_centre(&prior_cd, prior);
    refit_until_stable(pairs, &weights, kept, centre, crpix)
}

fn check_inputs(
    pairs: &[StarPair],
    (crval_ra_deg, crval_dec_deg): (f64, f64),
    focal_mm: Option<f64>,
    pixel_um: Option<f64>,
) -> Result<(), SolveError> {
    input::finite("crval_ra_deg", crval_ra_deg)?;
    input::declination("crval_dec_deg", crval_dec_deg)?;
    input::positive_if_set("focal_mm", focal_mm)?;
    input::positive_if_set("pixel_um", pixel_um)?;
    input::pairs(pairs)
}

pub(super) fn check_pair_count(n: usize) -> Result<(), SolveError> {
    if n < 3 {
        return Err(SolveError::fit_failed(format!(
            "need at least 3 matched pairs, got {n}"
        )));
    }
    Ok(())
}

pub(super) fn check_projected_count(n: usize) -> Result<(), SolveError> {
    if n < 3 {
        return Err(SolveError::fit_failed(
            "too few stars projected onto tangent plane",
        ));
    }
    Ok(())
}

// RANSAC inlier selection, then a weighted CD fit, sigma clipping and a final weighted
// fit. Returns the CD parameters and the mask of pairs they were fitted to.
fn fit_inlier_cd(tp: &TangentProjection<'_>) -> Result<(CdParams, Vec<bool>), SolveError> {
    // Five times the per-axis positional error. Against the exact pixel map, a correct
    // pair's residual exceeds it with probability exp(-12.5), about 4e-6, under that
    // Gaussian error. RANSAC tests pairs against a model fitted to 3 of them: the
    // residual then has per-axis variance σ²(1 + h), where σ is `POSITION_SIGMA_PX` and
    // h is the pair's leverage under the 3-point fit, and lens distortion adds an
    // offset that no affine model removes. Correct pairs far from the sample can
    // therefore miss the tolerance, and the least-squares refit on the inliers
    // recovers only those it brings back within it.
    let mut mask = ransac_fit(&tp.rows, &tp.pairs, INLIER_PX)?;
    let mut cd = fit_cd_weighted(&tp.rows, &mask, &tp.weights)?;
    median_clip_loop(&tp.rows, &mut mask, &mut cd, &tp.weights)?;
    cd = fit_cd_weighted(&tp.rows, &mask, &tp.weights)?;
    Ok((cd, mask))
}

const INLIER_PX: f64 = 5.0 * POSITION_SIGMA_PX;

// Least-squares weight of a pair: the inverse of the variance of its centroid on one
// axis, the sum of two terms.
// - The Cramér-Rao bound on the centroid of a Gaussian star in background-limited noise
//   is sqrt(8π)·s²·σ_bg/F px per axis, for PSF sigma s px, background noise σ_bg per
//   pixel and flux F. With F = 2π·s²·A for peak A and SNR = A/σ_bg, this is
//   sqrt(2/π)/SNR px, whatever the PSF width.
// - `POSITION_SIGMA_PX`, the 1 px centroid sigma of astrometry.net 0.97
//   (DEFAULT_VERIFY_PIX), covers the errors the bound leaves out: undersampling, blends
//   and lens distortion that the affine fit cannot model. With it, no pair outweighs a
//   pair of SNR s by more than a factor of 1 + 2/(π·s²).
fn pair_weight(snr: f64) -> f64 {
    let centroid_sigma = libm::sqrt(2.0 / PI) / snr;
    1.0 / (centroid_sigma * centroid_sigma + POSITION_SIGMA_PX * POSITION_SIGMA_PX)
}

// The pairs within `INLIER_PX` of `model`, the test RANSAC applies, clipped by
// `median_clip_loop`, and the CD fitted to them. The tolerance admits a pair a few
// pixels off, from a blend or a saturated core. The clip's first pass measures
// residuals against `model`, which was fitted before such pairs were readmitted, so a
// group of them cannot pull the fit through themselves before the clip sees them.
fn select_inliers(
    tp: &TangentProjection<'_>,
    model: CdParams,
) -> Result<(CdParams, Vec<bool>), SolveError> {
    let all = vec![true; tp.rows.len()];
    let errors = compute_pixel_errors(&tp.rows, &all, model);
    let mut mask: Vec<bool> = errors.iter().map(|&e| e <= INLIER_PX).collect();
    let mut cd = model;
    median_clip_loop(&tp.rows, &mut mask, &mut cd, &tp.weights)?;
    cd = fit_cd_weighted(&tp.rows, &mask, &tp.weights)?;
    Ok((cd, mask))
}

// The sky position of the fitted tangent-plane offset (`xi0`, `eta0`) about `tangent`.
fn fitted_centre(cd: &CdParams, tangent: (f64, f64)) -> (f64, f64) {
    let (xi0, eta0) = (cd.xi0 * DEG_TO_RAD, cd.eta0 * DEG_TO_RAD);
    tan_deproject_star(xi0, eta0, tangent.0, tangent.1)
}

// The selection can alternate between two sets when a pair sits at the inlier threshold;
// with exact pairs it repeats on the second pass.
const MAX_RESELECTIONS: usize = 4;

type Kept<'a> = (Vec<&'a StarPair>, Vec<f64>);

// Inliers chosen about a tangent point far from the field centre fit a model that is
// wrong there, so they exclude correct pairs and can include wrong ones. The CD matrix
// is refitted about `centre` and the inliers are chosen again from all `pairs` against
// it, moving the centre each time, until the inliers repeat. The inliers can repeat on
// the first pass while the refit's CRVAL is still away from `centre`, by 0.7" for a
// prior 1.2 degrees off, so the kept pairs are fitted once more about that CRVAL.
// Returns the last refit, its CRVAL and the pairs it was fitted on.
fn refit_until_stable<'a>(
    pairs: &'a [StarPair],
    weights: &[f64],
    mut kept: Kept<'a>,
    mut centre: (f64, f64),
    crpix: (f64, f64),
) -> Result<(CdParams, (f64, f64), Kept<'a>), SolveError> {
    let (mut cd, mut crval) = refit_at(&kept.0, &kept.1, centre, crpix)?;
    for _ in 0..MAX_RESELECTIONS {
        let tp = project_pairs_to_tangent(pairs, weights, centre.0, centre.1, crpix.0, crpix.1);
        let Ok((selected_cd, mask)) = select_inliers(&tp, cd) else {
            break;
        };
        let selected = extract_kept_pairs(&tp, &mask);
        if same_pairs(&selected.0, &kept.0) {
            break;
        }
        kept = selected;
        centre = fitted_centre(&selected_cd, centre);
        (cd, crval) = refit_at(&kept.0, &kept.1, centre, crpix)?;
    }
    (cd, crval) = refit_at(&kept.0, &kept.1, crval, crpix)?;
    Ok((cd, crval, kept))
}

// Both lists borrow from the same slice of pairs, in its order.
fn same_pairs(a: &[&StarPair], b: &[&StarPair]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(&x, &y)| std::ptr::eq(x, y))
}

// A CD matrix fitted about one tangent point does not describe the projection about
// another, and the difference grows with their separation and the declination. The
// kept pairs are reprojected about `centre` and the CD matrix is fitted once more; the
// CRVAL moves by that fit's offset, which is second order in the first.
fn refit_at(
    kept: &[&StarPair],
    weights: &[f64],
    centre: (f64, f64),
    (crpix1, crpix2): (f64, f64),
) -> Result<(CdParams, (f64, f64)), SolveError> {
    let kept: Vec<StarPair> = kept.iter().map(|&p| p.clone()).collect();
    let tp = project_pairs_to_tangent(&kept, weights, centre.0, centre.1, crpix1, crpix2);
    let all = vec![true; tp.rows.len()];
    let cd = fit_cd_weighted(&tp.rows, &all, &tp.weights)?;
    let crval = fitted_centre(&cd, centre);
    Ok((cd, crval))
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_catalog::query::quad::tan_deproject_star;
    use celestial_core::constants::TWOPI;

    struct SynthGrid {
        crval_ra: f64,
        crval_dec: f64,
        cd1_1: f64,
        cd2_2: f64,
        crpix1: f64,
        crpix2: f64,
        grid_step: f64,
        grid_half: i32,
    }

    fn synthesize_pairs(g: &SynthGrid) -> Vec<StarPair> {
        let mut pairs = Vec::new();
        for gy in -g.grid_half..=g.grid_half {
            for gx in -g.grid_half..=g.grid_half {
                let px_x = g.crpix1 + gx as f64 * g.grid_step;
                let px_y = g.crpix2 + gy as f64 * g.grid_step;
                let u = px_x - g.crpix1;
                let v = px_y - g.crpix2;
                let xi_deg = g.cd1_1 * u;
                let eta_deg = g.cd2_2 * v;
                let (ra, dec) = tan_deproject_star(
                    xi_deg * DEG_TO_RAD,
                    eta_deg * DEG_TO_RAD,
                    g.crval_ra,
                    g.crval_dec,
                );
                pairs.push(StarPair {
                    px_x,
                    px_y,
                    ra_deg: ra,
                    dec_deg: dec,
                    votes: 4,
                    snr: 50.0,
                });
            }
        }
        pairs
    }

    #[test]
    fn fit_wcs_recovers_known_cd_matrix_and_crval() {
        let crval_ra = 180.0;
        let crval_dec = 0.0;
        let cd1_1 = 0.001;
        let cd2_2 = 0.001;
        let width = 1000;
        let height = 1000;
        let pairs = synthesize_pairs(&SynthGrid {
            crval_ra,
            crval_dec,
            cd1_1,
            cd2_2,
            crpix1: width as f64 / 2.0,
            crpix2: height as f64 / 2.0,
            grid_step: 20.0,
            grid_half: 5,
        });
        let wcs = fit_wcs(&pairs, width, height, crval_ra, crval_dec, None, None).unwrap();

        // crpix stored as width/2, height/2
        assert_eq!(wcs.crpix1, 500.0);
        assert_eq!(wcs.crpix2, 500.0);
        // The pairs were synthesized at CRVAL (180, 0) with CD = diag(0.001, 0.001) and
        // no noise. The fitted Dec is 2.8e-19 deg from it, the CD terms at most 4.4e-17
        // deg/px (4.4e-14 of the scale), and the RMS is 6.3e-12 px.
        assert_eq!((wcs.crval1, wcs.crval2), (180.0, 2.7808532181887886e-19));
        assert_eq!(
            (wcs.cd1_1, wcs.cd1_2, wcs.cd2_1, wcs.cd2_2),
            (
                0.0010000000000000438,
                -7.606505436419194e-21,
                -1.1277598793763638e-20,
                0.0009999999999999998
            )
        );
        assert_eq!(wcs.rms_px, 6.263205042073299e-12);
    }

    #[test]
    fn fit_wcs_preserves_focal_and_pixel_metadata() {
        let pairs = synthesize_pairs(&SynthGrid {
            crval_ra: 180.0,
            crval_dec: 0.0,
            cd1_1: 0.001,
            cd2_2: 0.001,
            crpix1: 500.0,
            crpix2: 500.0,
            grid_step: 20.0,
            grid_half: 5,
        });
        let wcs = fit_wcs(&pairs, 1000, 1000, 180.0, 0.0, Some(500.0), Some(3.76)).unwrap();
        assert_eq!(wcs.focal_mm, Some(500.0));
        assert_eq!(wcs.pixel_um, Some(3.76));
    }

    // A field at Dec 60 with the prior 1.2 degrees off its centre along both axes. The
    // pairs are exact. A CD matrix fitted about the prior leaves residuals of 28 px RMS,
    // and sigma clipping about the prior drops 2 of the 361 pairs; the selection is
    // repeated about the fitted centre, which keeps them all.
    #[test]
    fn fit_wcs_refits_the_cd_matrix_about_the_fitted_centre() {
        let pairs = synthesize_pairs(&SynthGrid {
            crval_ra: 180.0,
            crval_dec: 60.0,
            cd1_1: -0.001,
            cd2_2: 0.001,
            crpix1: 1500.0,
            crpix2: 1000.0,
            grid_step: 100.0,
            grid_half: 9,
        });
        let wcs = fit_wcs(&pairs, 3000, 2000, 182.4, 61.2, None, None).unwrap();
        assert_eq!((pairs.len(), wcs.n_stars), (361, 361));
        assert_eq!(wcs.rms_px, 2.2615635925947653e-7);
    }

    // Deterministic uniform and Gaussian samples for the pulled-fit scenes.
    struct Lcg(u64);

    impl Lcg {
        fn uniform(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
        }

        fn gaussian(&mut self) -> f64 {
            let (a, b) = (self.uniform(), self.uniform());
            libm::sqrt(-2.0 * libm::log(a)) * libm::cos(TWOPI * b)
        }
    }

    // 150 correct pairs on a 3000 x 2000 field at 1.5"/px, with 0.2 px centroid noise
    // and SNR 20 to 220. The first `n_bad` pairs instead sit `bad_px` from their true
    // pixel at SNR `bad_snr`, as a blend or a saturated core would. Returns the pairs
    // and the true pixel of each.
    fn pulled_scene(n_bad: usize, bad_px: f64, bad_snr: f64) -> (Vec<StarPair>, Vec<(f64, f64)>) {
        let s = 1.5 / 3600.0;
        let mut rng = Lcg(42);
        let mut pairs = Vec::new();
        let mut truth = Vec::new();
        for i in 0..150 {
            let x = 60.0 + rng.uniform() * 2880.0;
            let y = 60.0 + rng.uniform() * 1880.0;
            let (xi, eta) = (-s * (x - 1500.0), -s * (y - 1000.0));
            let (ra, dec) = tan_deproject_star(xi * DEG_TO_RAD, eta * DEG_TO_RAD, 150.0, 50.0);
            truth.push((x, y));
            let snr = 20.0 + 200.0 * rng.uniform();
            let mut p = StarPair {
                px_x: x + 0.2 * rng.gaussian(),
                px_y: y + 0.2 * rng.gaussian(),
                ra_deg: ra,
                dec_deg: dec,
                votes: 4,
                snr,
            };
            if i < n_bad {
                let a = rng.uniform() * TWOPI;
                p.px_x = x + bad_px * libm::cos(a);
                p.px_y = y + bad_px * libm::sin(a);
                p.snr = bad_snr;
            }
            pairs.push(p);
        }
        (pairs, truth)
    }

    // Largest distance between a pair's true pixel and the fitted WCS's pixel for its
    // catalog position.
    fn max_error_px(wcs: &WcsSolution, pairs: &[StarPair], truth: &[(f64, f64)]) -> f64 {
        let errors = pairs.iter().zip(truth).filter_map(|(p, t)| {
            let (x, y) = wcs.sky_to_pixel(p.ra_deg, p.dec_deg)?;
            Some(libm::hypot(x - t.0, y - t.1))
        });
        errors.fold(0.0, f64::max)
    }

    #[test]
    fn fit_wcs_is_not_pulled_by_one_bright_pair_inside_the_inlier_tolerance() {
        let (pairs, truth) = pulled_scene(1, 4.5, 2000.0);
        let wcs = fit_wcs(&pairs, 3000, 2000, 150.3, 49.8, None, None).unwrap();
        assert_eq!(
            (wcs.n_stars, max_error_px(&wcs, &pairs, &truth)),
            (149, 0.059001181086606645)
        );
    }

    #[test]
    fn fit_wcs_keeps_correct_pairs_beside_three_bright_pairs_2_px_off() {
        let (pairs, truth) = pulled_scene(3, 2.0, 5000.0);
        let wcs = fit_wcs(&pairs, 3000, 2000, 150.3, 49.8, None, None).unwrap();
        assert_eq!(
            (wcs.n_stars, max_error_px(&wcs, &pairs, &truth)),
            (147, 0.08856145248223467)
        );
    }

    #[test]
    fn fit_wcs_is_not_pulled_by_twenty_bright_pairs_2_px_off() {
        let (pairs, truth) = pulled_scene(20, 2.0, 2000.0);
        let wcs = fit_wcs(&pairs, 3000, 2000, 150.3, 49.8, None, None).unwrap();
        assert_eq!(
            (wcs.n_stars, max_error_px(&wcs, &pairs, &truth)),
            (130, 0.13073350933207176)
        );
    }

    #[test]
    fn fit_wcs_is_not_pulled_by_thirty_bright_pairs_2_px_off() {
        let (pairs, truth) = pulled_scene(30, 2.0, 5000.0);
        let wcs = fit_wcs(&pairs, 3000, 2000, 150.3, 49.8, None, None).unwrap();
        assert_eq!(
            (wcs.n_stars, max_error_px(&wcs, &pairs, &truth)),
            (120, 0.09723393143248982)
        );
    }

    #[test]
    fn fit_wcs_is_not_pulled_by_fifty_bright_pairs_3_px_off() {
        let (pairs, truth) = pulled_scene(50, 3.0, 2000.0);
        let wcs = fit_wcs(&pairs, 3000, 2000, 150.3, 49.8, None, None).unwrap();
        assert_eq!(
            (wcs.n_stars, max_error_px(&wcs, &pairs, &truth)),
            (100, 0.06679714876678405)
        );
    }

    // A pair more than 90 degrees from the tangent point has no gnomonic projection.
    fn behind_tangent_plane() -> StarPair {
        StarPair {
            px_x: 1.0,
            px_y: 1.0,
            ra_deg: 0.0,
            dec_deg: 0.0,
            votes: 4,
            snr: 50.0,
        }
    }

    #[test]
    fn fit_wcs_ignores_pairs_behind_the_tangent_plane() {
        let grid = synthesize_pairs(&SynthGrid {
            crval_ra: 180.0,
            crval_dec: 0.0,
            cd1_1: 0.001,
            cd2_2: 0.001,
            crpix1: 500.0,
            crpix2: 500.0,
            grid_step: 20.0,
            grid_half: 5,
        });
        let mut pairs = vec![behind_tangent_plane()];
        pairs.extend(grid.iter().cloned());
        let wcs = fit_wcs(&pairs, 1000, 1000, 180.0, 0.0, None, None).unwrap();
        let fitted: Vec<(f64, f64)> = wcs.residuals.iter().map(|r| (r.px_x, r.px_y)).collect();
        let expected: Vec<(f64, f64)> = grid.iter().map(|p| (p.px_x, p.px_y)).collect();
        assert_eq!(fitted, expected);
    }

    #[test]
    fn fit_wcs_rejects_fewer_than_three_pairs() {
        let pairs = vec![
            StarPair {
                px_x: 0.0,
                px_y: 0.0,
                ra_deg: 0.0,
                dec_deg: 0.0,
                votes: 4,
                snr: 10.0,
            },
            StarPair {
                px_x: 1.0,
                px_y: 0.0,
                ra_deg: 0.1,
                dec_deg: 0.0,
                votes: 4,
                snr: 10.0,
            },
        ];
        let err = fit_wcs(&pairs, 100, 100, 0.0, 0.0, None, None).unwrap_err();
        assert!(err.to_string().contains("at least 3 matched pairs"));
    }

    #[test]
    fn fit_wcs_rejects_non_finite_inputs() {
        let pairs = synthesize_pairs(&SynthGrid {
            crval_ra: 180.0,
            crval_dec: 0.0,
            cd1_1: 0.001,
            cd2_2: 0.001,
            crpix1: 500.0,
            crpix2: 500.0,
            grid_step: 20.0,
            grid_half: 2,
        });
        let mut bad_pair = pairs.clone();
        bad_pair[2] = StarPair {
            px_x: 1.0,
            px_y: 2.0,
            ra_deg: f64::NAN,
            dec_deg: 4.0,
            votes: 1,
            snr: 50.0,
        };
        let fit = |pairs: &[StarPair], ra: f64, dec: f64, focal: Option<f64>| {
            let out = fit_wcs(pairs, 1000, 1000, ra, dec, focal, None);
            out.err().map(|e| e.describe())
        };
        let cases = [
            (
                fit(&pairs, f64::NAN, 0.0, None),
                "InvalidParameter: invalid parameter crval_ra_deg: must be finite, got NaN",
            ),
            (
                fit(&pairs, 180.0, 95.0, None),
                "InvalidParameter: invalid parameter crval_dec_deg: must be from -90 to 90, got 95",
            ),
            (
                fit(&pairs, 180.0, 0.0, Some(f64::NAN)),
                "InvalidParameter: invalid parameter focal_mm: must be finite and positive, got NaN",
            ),
            (
                fit(&bad_pair, 180.0, 0.0, None),
                "InvalidParameter: invalid parameter pairs: pair 2 is not finite: pixel (1, 2), sky (NaN, 4), snr 50",
            ),
        ];
        let (got, want): (Vec<_>, Vec<_>) = cases
            .iter()
            .map(|(err, expected)| (err.as_deref(), Some(*expected)))
            .unzip();
        assert_eq!(got, want);
    }
}
