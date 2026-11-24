use celestial_catalog::query::quad::tan_project_star;
use celestial_core::constants::RAD_TO_DEG;

use crate::error::SolveError;
use crate::match_field::StarPair;

use super::median_clip::median_clip;
use super::normal_equations::solve_normal_equations;
use super::solution::{InverseCd, StarResidual};

#[derive(Debug, Clone, Copy)]
pub(super) struct CdParams {
    pub(super) cd1_1: f64,
    pub(super) cd1_2: f64,
    pub(super) xi0: f64,
    pub(super) cd2_1: f64,
    pub(super) cd2_2: f64,
    pub(super) eta0: f64,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct WcsParams {
    pub(super) crpix1: f64,
    pub(super) crpix2: f64,
    pub(super) crval_ra: f64,
    pub(super) crval_dec: f64,
    pub(super) cd1_1: f64,
    pub(super) cd1_2: f64,
    pub(super) cd2_1: f64,
    pub(super) cd2_2: f64,
}

// The inverse of an affine CD fit: maps a tangent-plane position back to a pixel offset.
pub(super) struct InverseAffine {
    cd: InverseCd,
    xi0: f64,
    eta0: f64,
}

impl InverseAffine {
    pub(super) fn new(p: CdParams) -> Self {
        Self {
            cd: InverseCd::new(p.cd1_1, p.cd1_2, p.cd2_1, p.cd2_2),
            xi0: p.xi0,
            eta0: p.eta0,
        }
    }

    // Distance in pixels between the row's pixel offset and the one the model predicts
    // from the row's tangent-plane position.
    pub(super) fn pixel_error(&self, row: &(f64, f64, f64, f64)) -> f64 {
        let (pred_u, pred_v) = self.cd.apply(row.2 - self.xi0, row.3 - self.eta0);
        let ex = row.0 - pred_u;
        let ey = row.1 - pred_v;
        libm::sqrt(ex * ex + ey * ey)
    }
}

pub(super) fn compute_pixel_errors(
    projected: &[(f64, f64, f64, f64)],
    mask: &[bool],
    p: CdParams,
) -> Vec<f64> {
    let model = InverseAffine::new(p);
    projected
        .iter()
        .zip(mask)
        .map(|(row, &on)| if on { model.pixel_error(row) } else { 0.0 })
        .collect()
}

pub(super) fn fit_cd_weighted(
    projected: &[(f64, f64, f64, f64)],
    mask: &[bool],
    weights: &[f64],
) -> Result<CdParams, SolveError> {
    let (cd1_1, cd1_2, xi0) = solve_3x3_weighted(projected, mask, weights, |t| t.2)?;
    let (cd2_1, cd2_2, eta0) = solve_3x3_weighted(projected, mask, weights, |t| t.3)?;
    check_cd_determinant(cd1_1 * cd2_2 - cd1_2 * cd2_1)?;
    Ok(CdParams {
        cd1_1,
        cd1_2,
        xi0,
        cd2_1,
        cd2_2,
        eta0,
    })
}

// A CD matrix whose determinant is NaN or at most 1e-30 in size is singular.
pub(super) fn check_cd_determinant(det: f64) -> Result<(), SolveError> {
    let size = libm::fabs(det);
    if size.is_nan() || size <= 1e-30 {
        return Err(SolveError::fit_failed("singular CD matrix"));
    }
    Ok(())
}

// Clips the pairs whose unweighted residual exceeds `clip_threshold` and refits the CD
// matrix with `weights` after each pass, until a pass clips nothing or fewer than 3
// pairs remain. The median ignores the weights, so a group of heavily weighted pairs
// that pulls the fit towards itself does not raise the threshold.
pub(super) fn median_clip_loop(
    projected: &[(f64, f64, f64, f64)],
    mask: &mut [bool],
    cd: &mut CdParams,
    weights: &[f64],
) -> Result<(), SolveError> {
    median_clip(
        mask,
        cd,
        |cd, mask| compute_pixel_errors(projected, mask, *cd),
        |mask| fit_cd_weighted(projected, mask, weights),
    )
}

fn solve_3x3_weighted(
    data: &[(f64, f64, f64, f64)],
    mask: &[bool],
    weights: &[f64],
    target: impl Fn(&(f64, f64, f64, f64)) -> f64,
) -> Result<(f64, f64, f64), SolveError> {
    let mut s_uu = 0.0_f64;
    let mut s_uv = 0.0_f64;
    let mut s_u = 0.0_f64;
    let mut s_vv = 0.0_f64;
    let mut s_v = 0.0_f64;
    let mut s_1 = 0.0_f64;
    let mut s_ub = 0.0_f64;
    let mut s_vb = 0.0_f64;
    let mut s_b = 0.0_f64;

    for (i, row) in data.iter().enumerate() {
        if !mask[i] {
            continue;
        }
        let w = weights[i];
        let (u, v) = (row.0, row.1);
        let b = target(row);
        s_uu += w * u * u;
        s_uv += w * u * v;
        s_u += w * u;
        s_vv += w * v * v;
        s_v += w * v;
        s_1 += w;
        s_ub += w * u * b;
        s_vb += w * v * b;
        s_b += w * b;
    }

    let mut a = [s_uu, s_uv, s_u, s_uv, s_vv, s_v, s_u, s_v, s_1];
    let mut x = [s_ub, s_vb, s_b];
    solve_normal_equations(&mut a, &mut x, 3)?;
    Ok((x[0], x[1], x[2]))
}

pub(super) fn compute_residuals(pairs: &[&StarPair], wcs: WcsParams) -> Vec<StarResidual> {
    let det = wcs.cd1_1 * wcs.cd2_2 - wcs.cd1_2 * wcs.cd2_1;
    if libm::fabs(det) < 1e-30 {
        return Vec::new();
    }
    let inverse = InverseCd::new(wcs.cd1_1, wcs.cd1_2, wcs.cd2_1, wcs.cd2_2);
    pairs
        .iter()
        .filter_map(|p| pair_residual(p, &wcs, &inverse))
        .collect()
}

// The residual of `p` under `wcs`, whose CD matrix inverts to `inverse`. `None` when the
// pair has no gnomonic projection about the CRVAL.
fn pair_residual(p: &StarPair, wcs: &WcsParams, inverse: &InverseCd) -> Option<StarResidual> {
    let (xi_rad, eta_rad) = tan_project_star(p.ra_deg, p.dec_deg, wcs.crval_ra, wcs.crval_dec)?;
    let xi = xi_rad * RAD_TO_DEG;
    let eta = eta_rad * RAD_TO_DEG;
    let (u, v) = inverse.apply(xi, eta);
    let (pred_x, pred_y) = (u + wcs.crpix1, v + wcs.crpix2);
    let err_x = p.px_x - pred_x;
    let err_y = p.px_y - pred_y;
    let err_px = libm::sqrt(err_x * err_x + err_y * err_y);
    Some(StarResidual {
        px_x: p.px_x,
        px_y: p.px_y,
        err_x,
        err_y,
        err_px,
        ra_deg: p.ra_deg,
        dec_deg: p.dec_deg,
    })
}

// The RMS of the residual sizes `errors` in pixels. 0 when there are none.
pub(super) fn rms_error(errors: impl IntoIterator<Item = f64>) -> f64 {
    let (sum_sq, n) = errors
        .into_iter()
        .fold((0.0_f64, 0_usize), |(s, n), e| (s + e * e, n + 1));
    if n == 0 {
        return 0.0;
    }
    libm::sqrt(sum_sq / n as f64)
}

// The RMS of the residual sizes `errors`, each weighted by its entry of `weights`. 0 when
// the weights sum to 0 or less.
pub(super) fn weighted_rms_error(errors: impl IntoIterator<Item = f64>, weights: &[f64]) -> f64 {
    let (sum_w, sum_we2) = errors
        .into_iter()
        .zip(weights)
        .fold((0.0_f64, 0.0_f64), |(sw, swe), (e, &w)| {
            (sw + w, swe + w * e * e)
        });
    if sum_w > 0.0 {
        libm::sqrt(sum_we2 / sum_w)
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn residual(err_x: f64, err_y: f64) -> StarResidual {
        let err_px = libm::sqrt(err_x * err_x + err_y * err_y);
        StarResidual {
            px_x: 0.0,
            px_y: 0.0,
            err_x,
            err_y,
            err_px,
            ra_deg: 0.0,
            dec_deg: 0.0,
        }
    }

    fn pair(px_x: f64, px_y: f64, ra_deg: f64, dec_deg: f64) -> StarPair {
        StarPair {
            px_x,
            px_y,
            ra_deg,
            dec_deg,
            votes: 4,
            snr: 20.0,
        }
    }

    #[test]
    fn rms_error_empty_returns_zero() {
        assert_eq!(rms_error([]), 0.0);
    }

    #[test]
    fn rms_error_matches_hand_calculation() {
        // residuals with err_px 3, 4 -> sqrt((9+16)/2) = sqrt(12.5)
        let residuals = [residual(3.0, 0.0), residual(0.0, 4.0)];
        let rms = rms_error(residuals.iter().map(|r| r.err_px));
        assert_eq!(rms, libm::sqrt(12.5));
    }

    #[test]
    fn weighted_rms_error_empty_returns_zero() {
        assert_eq!(weighted_rms_error([], &[]), 0.0);
    }

    #[test]
    fn weighted_rms_error_zero_weight_returns_zero() {
        let residuals = [residual(3.0, 0.0)];
        let weights = vec![0.0];
        assert_eq!(weighted_rms_error([residuals[0].err_px], &weights), 0.0);
    }

    #[test]
    fn weighted_rms_error_shifts_toward_heavily_weighted_residuals() {
        // equal weights: rms = sqrt((9+16)/2) = sqrt(12.5)
        let residuals = [residual(3.0, 0.0), residual(0.0, 4.0)];
        let sizes = || residuals.iter().map(|r| r.err_px);
        let equal = weighted_rms_error(sizes(), &[1.0, 1.0]);
        assert_eq!(equal, libm::sqrt(12.5));

        // heavier on the bigger residual -> rms pulled up
        let skewed_big = weighted_rms_error(sizes(), &[1.0, 100.0]);
        assert!(skewed_big > equal);

        // heavier on the smaller residual -> rms pulled down
        let skewed_small = weighted_rms_error(sizes(), &[100.0, 1.0]);
        assert!(skewed_small < equal);
    }

    // Pairs on one line in the pixel plane fix no CD matrix. Rounding leaves pivots far
    // above 1e-30 but below the working precision of the sums, so an absolute floor
    // returned a CD matrix of about 1e12 deg/px for some of these lines.
    #[test]
    fn fit_cd_weighted_rejects_pairs_on_a_line() {
        for (step, slope) in [(0.1, 3.0), (0.7, 0.3), (1.1, 1.7), (13.3, 0.1)] {
            let rows: Vec<(f64, f64, f64, f64)> = (1..=5)
                .map(|k| {
                    let u = step * k as f64;
                    let v = u * slope;
                    (u, v, 0.001 * u, 0.002 * v)
                })
                .collect();
            for n in [3, 5] {
                let fit = fit_cd_weighted(&rows[..n], &vec![true; n], &vec![1.0; n]);
                assert!(
                    matches!(fit, Err(SolveError::FitFailed { .. })),
                    "step {step}, slope {slope}, {n} rows: {:?}",
                    fit.map(|c| (c.cd1_1, c.cd1_2, c.cd2_1, c.cd2_2))
                );
            }
        }
    }

    #[test]
    fn compute_residuals_zero_when_pair_lies_on_wcs() {
        // With identity CD at 1 deg/px and crval at origin, a pair at pixel
        // (crpix + 10, crpix) should correspond to RA exactly 10 deg away.
        let wcs = WcsParams {
            crpix1: 0.0,
            crpix2: 0.0,
            crval_ra: 0.0,
            crval_dec: 0.0,
            cd1_1: 1.0,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: 1.0,
        };

        // Construct a catalog position, project it through the WCS forward to get
        // the pixel, then feed back in. That round-trip is exact.
        let ra = 0.0;
        let dec = 0.0;
        let p = pair(0.0, 0.0, ra, dec);
        let refs: Vec<&StarPair> = vec![&p];
        let residuals = compute_residuals(&refs, wcs);
        assert_eq!(residuals.len(), 1);
        // At the tangent point, predicted = crpix exactly.
        assert_eq!(residuals[0].err_px, 0.0);
    }

    #[test]
    fn compute_residuals_returns_empty_for_singular_cd() {
        let wcs = WcsParams {
            crpix1: 0.0,
            crpix2: 0.0,
            crval_ra: 0.0,
            crval_dec: 0.0,
            cd1_1: 0.0,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: 0.0,
        };
        let p = pair(100.0, 100.0, 10.0, 10.0);
        let refs: Vec<&StarPair> = vec![&p];
        let residuals = compute_residuals(&refs, wcs);
        assert!(residuals.is_empty());
    }

    fn synthesize_pairs(
        crval_ra: f64,
        crval_dec: f64,
        cd: &CdParams,
        samples: &[(f64, f64)],
    ) -> Vec<StarPair> {
        // Forward-project (u, v) -> (xi, eta) via CD -> (ra, dec) via tan deprojection
        use celestial_catalog::query::quad::tan_deproject_star;
        use celestial_core::constants::DEG_TO_RAD;

        samples
            .iter()
            .map(|&(u, v)| {
                let xi_deg = cd.cd1_1 * u + cd.cd1_2 * v + cd.xi0;
                let eta_deg = cd.cd2_1 * u + cd.cd2_2 * v + cd.eta0;
                let (ra, dec) = tan_deproject_star(
                    xi_deg * DEG_TO_RAD,
                    eta_deg * DEG_TO_RAD,
                    crval_ra,
                    crval_dec,
                );
                pair(u, v, ra, dec)
            })
            .collect()
    }

    #[test]
    fn fit_cd_weighted_recovers_known_cd_matrix() {
        let crval_ra = 180.0;
        let crval_dec = 0.0;
        let cd = CdParams {
            cd1_1: 0.001,
            cd1_2: 0.0,
            xi0: 0.0,
            cd2_1: 0.0,
            cd2_2: 0.001,
            eta0: 0.0,
        };
        let samples = [
            (-100.0, -100.0),
            (100.0, -100.0),
            (0.0, 0.0),
            (-100.0, 100.0),
            (100.0, 100.0),
            (50.0, -50.0),
        ];
        let pairs = synthesize_pairs(crval_ra, crval_dec, &cd, &samples);

        // Project through the "known" WCS to get (u, v, xi, eta)
        use celestial_catalog::query::quad::tan_project_star;
        use celestial_core::constants::RAD_TO_DEG;
        let projected: Vec<(f64, f64, f64, f64)> = pairs
            .iter()
            .filter_map(|p| {
                tan_project_star(p.ra_deg, p.dec_deg, crval_ra, crval_dec)
                    .map(|(xi, eta)| (p.px_x, p.px_y, xi * RAD_TO_DEG, eta * RAD_TO_DEG))
            })
            .collect();

        let mask = vec![true; projected.len()];
        let weights = vec![1.0; projected.len()];
        let fit = fit_cd_weighted(&projected, &mask, &weights).unwrap();

        // The pairs were synthesized with CD = diag(0.001, 0.001). The fitted terms
        // differ from it by at most 4.3e-17 deg/px, 4.3e-14 of the scale.
        assert_eq!(
            (fit.cd1_1, fit.cd1_2, fit.cd2_1, fit.cd2_2),
            (
                0.0010000000000000425,
                5.16672658535462e-19,
                -1.3862385157639302e-21,
                0.0010000000000000002
            )
        );
    }

    #[test]
    fn median_clip_loop_removes_obvious_outlier() {
        let crval_ra = 180.0;
        let crval_dec = 0.0;
        let cd = CdParams {
            cd1_1: 0.001,
            cd1_2: 0.0,
            xi0: 0.0,
            cd2_1: 0.0,
            cd2_2: 0.001,
            eta0: 0.0,
        };
        // 20 clean points in a grid so the fit is well-determined even with the
        // outlier present.
        let mut samples: Vec<(f64, f64)> = Vec::new();
        for gy in -2..=2 {
            for gx in -2..=2 {
                samples.push((gx as f64 * 20.0, gy as f64 * 20.0));
            }
        }
        // one outrageous outlier
        samples.push((100.0, 100.0));
        let pairs = synthesize_pairs(crval_ra, crval_dec, &cd, &samples);

        use celestial_catalog::query::quad::tan_project_star;
        use celestial_core::constants::RAD_TO_DEG;
        let mut projected: Vec<(f64, f64, f64, f64)> = pairs
            .iter()
            .filter_map(|p| {
                tan_project_star(p.ra_deg, p.dec_deg, crval_ra, crval_dec)
                    .map(|(xi, eta)| (p.px_x, p.px_y, xi * RAD_TO_DEG, eta * RAD_TO_DEG))
            })
            .collect();

        // Fit an initial CD matrix on the clean points only — this is the normal
        // bootstrap for sigma-clipping: you start from a reasonable model.
        let last = projected.len() - 1;
        let mut mask = vec![true; projected.len()];
        mask[last] = false;
        let weights = vec![1.0; projected.len()];
        let mut cd_fit = fit_cd_weighted(&projected, &mask, &weights).unwrap();

        // Now enable the outlier in the mask and corrupt its pixel coord
        mask[last] = true;
        projected[last].0 += 500.0;
        projected[last].1 -= 500.0;

        median_clip_loop(&projected, &mut mask, &mut cd_fit, &weights).unwrap();
        // outlier must be flagged
        assert!(!mask[last], "expected outlier to be rejected");
        // clean points must be kept
        assert!(
            mask[..last].iter().all(|&m| m),
            "clean points should survive"
        );
    }
}
