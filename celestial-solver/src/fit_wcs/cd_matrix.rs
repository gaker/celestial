use anyhow::{ensure, Result};
use celestial_catalog::query::tan_project_star;
use celestial_core::constants::RAD_TO_DEG;

use crate::match_field::StarPair;

use super::solution::StarResidual;

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

pub(super) fn ransac_fit(
    projected: &[(f64, f64, f64, f64)],
    pairs: &[StarPair],
    inlier_threshold: f64,
) -> Result<Vec<bool>> {
    let n = projected.len();
    let max_trials = 1000;
    let mut best_inliers = 0usize;
    let mut best_mask = vec![false; n];

    let max_v = pairs.iter().map(|p| p.votes).max().unwrap_or(0);
    let high_vote: Vec<usize> = pairs.iter().enumerate()
        .filter(|(_, p)| p.votes >= max_v / 3)
        .map(|(i, _)| i)
        .filter(|&i| i < n)
        .collect();

    let sample_pool = if high_vote.len() >= 3 { &high_vote } else {
        return Ok(vec![true; n]);
    };

    let mut rng_state: u64 = 0xDEAD_BEEF_CAFE_1234;
    let mut next_rand = || -> u64 {
        rng_state ^= rng_state << 13;
        rng_state ^= rng_state >> 7;
        rng_state ^= rng_state << 17;
        rng_state
    };

    let all_on = vec![true; n];

    for _ in 0..max_trials {
        let a = sample_pool[next_rand() as usize % sample_pool.len()];
        let b = sample_pool[next_rand() as usize % sample_pool.len()];
        let c = sample_pool[next_rand() as usize % sample_pool.len()];
        if a == b || b == c || a == c {
            continue;
        }

        let mut sample_mask = vec![false; n];
        sample_mask[a] = true;
        sample_mask[b] = true;
        sample_mask[c] = true;

        let Ok(cd_sample) = fit_cd_matrix_6(projected, &sample_mask) else {
            continue;
        };
        if libm::fabs(cd_sample.cd1_1 * cd_sample.cd2_2 - cd_sample.cd1_2 * cd_sample.cd2_1) < 1e-30 {
            continue;
        }

        let errs = compute_pixel_errors(projected, &all_on, cd_sample);
        let inlier_count = errs.iter().filter(|&&e| e <= inlier_threshold).count();

        if inlier_count > best_inliers {
            best_inliers = inlier_count;
            best_mask = errs.iter().map(|&e| e <= inlier_threshold).collect();

            if let Ok(cd_refined) = fit_cd_matrix_6(projected, &best_mask) {
                let errs2 = compute_pixel_errors(projected, &all_on, cd_refined);
                let refined = errs2.iter().filter(|&&e| e <= inlier_threshold).count();
                if refined >= inlier_count {
                    best_inliers = refined;
                    best_mask = errs2.iter().map(|&e| e <= inlier_threshold).collect();
                }
            }
        }
    }

    ensure!(best_inliers >= 3, "RANSAC found no consistent model ({} pairs)", n);
    log::debug!("RANSAC: {} inliers from {} pairs", best_inliers, n);
    Ok(best_mask)
}

fn fit_cd_matrix_6(
    projected: &[(f64, f64, f64, f64)],
    mask: &[bool],
) -> Result<CdParams> {
    let (cd1_1, cd1_2, xi0) = solve_3x3(projected, mask, |t| t.2)?;
    let (cd2_1, cd2_2, eta0) = solve_3x3(projected, mask, |t| t.3)?;
    let det = cd1_1 * cd2_2 - cd1_2 * cd2_1;
    ensure!(libm::fabs(det) > 1e-30, "singular CD matrix");
    Ok(CdParams { cd1_1, cd1_2, xi0, cd2_1, cd2_2, eta0 })
}

pub(super) fn compute_pixel_errors(
    projected: &[(f64, f64, f64, f64)],
    mask: &[bool],
    p: CdParams,
) -> Vec<f64> {
    let det = p.cd1_1 * p.cd2_2 - p.cd1_2 * p.cd2_1;
    let inv_00 = p.cd2_2 / det;
    let inv_01 = -p.cd1_2 / det;
    let inv_10 = -p.cd2_1 / det;
    let inv_11 = p.cd1_1 / det;

    projected.iter().enumerate().map(|(i, row)| {
        if !mask[i] {
            return 0.0;
        }
        let dxi = row.2 - p.xi0;
        let deta = row.3 - p.eta0;
        let pred_u = inv_00 * dxi + inv_01 * deta;
        let pred_v = inv_10 * dxi + inv_11 * deta;
        let ex = row.0 - pred_u;
        let ey = row.1 - pred_v;
        libm::sqrt(ex * ex + ey * ey)
    }).collect()
}

pub(super) fn fit_cd_weighted(
    projected: &[(f64, f64, f64, f64)],
    mask: &[bool],
    weights: &[f64],
) -> Result<CdParams> {
    let (cd1_1, cd1_2, xi0) = solve_3x3_weighted(projected, mask, weights, |t| t.2)?;
    let (cd2_1, cd2_2, eta0) = solve_3x3_weighted(projected, mask, weights, |t| t.3)?;
    let det = cd1_1 * cd2_2 - cd1_2 * cd2_1;
    ensure!(libm::fabs(det) > 1e-30, "singular CD matrix");
    Ok(CdParams { cd1_1, cd1_2, xi0, cd2_1, cd2_2, eta0 })
}

pub(super) fn sigma_clip_loop_weighted(
    projected: &[(f64, f64, f64, f64)],
    mask: &mut [bool],
    cd: &mut CdParams,
    weights: &[f64],
) -> Result<()> {
    for iter in 0..10 {
        let errs = compute_pixel_errors(projected, mask, *cd);
        let mut sum_sr2 = 0.0_f64;
        let mut n_active = 0_usize;
        for (i, &e) in errs.iter().enumerate() {
            if !mask[i] { continue; }
            sum_sr2 += weights[i] * e * e;
            n_active += 1;
        }
        if n_active < 3 { break; }

        let sr_rms = libm::sqrt(sum_sr2 / n_active as f64);
        let threshold = 3.0 * sr_rms;

        let mut rejected = 0usize;
        for (i, &e) in errs.iter().enumerate() {
            if mask[i] && e * libm::sqrt(weights[i]) > threshold {
                mask[i] = false;
                rejected += 1;
            }
        }

        let n_active = mask.iter().filter(|&&m| m).count();
        let rms = libm::sqrt(
            errs.iter().enumerate()
                .filter(|(i, _)| mask[*i])
                .map(|(_, &e)| e * e)
                .sum::<f64>() / n_active.max(1) as f64
        );
        log::debug!("fit_wcs iter {}: {} stars, rms={:.3}px, rejected {}", iter, n_active, rms, rejected);

        if rejected == 0 || n_active < 3 {
            break;
        }

        *cd = fit_cd_weighted(projected, mask, weights)?;
    }
    Ok(())
}

fn solve_3x3_weighted(
    data: &[(f64, f64, f64, f64)],
    mask: &[bool],
    weights: &[f64],
    target: impl Fn(&(f64, f64, f64, f64)) -> f64,
) -> Result<(f64, f64, f64)> {
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

    let a = [
        [s_uu, s_uv, s_u],
        [s_uv, s_vv, s_v],
        [s_u, s_v, s_1],
    ];
    let rhs = [s_ub, s_vb, s_b];
    solve_3x3_system(&a, &rhs)
}

fn solve_3x3(
    data: &[(f64, f64, f64, f64)],
    mask: &[bool],
    target: impl Fn(&(f64, f64, f64, f64)) -> f64,
) -> Result<(f64, f64, f64)> {
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
        let (u, v) = (row.0, row.1);
        let b = target(row);
        s_uu += u * u;
        s_uv += u * v;
        s_u += u;
        s_vv += v * v;
        s_v += v;
        s_1 += 1.0;
        s_ub += u * b;
        s_vb += v * b;
        s_b += b;
    }

    let a = [
        [s_uu, s_uv, s_u],
        [s_uv, s_vv, s_v],
        [s_u, s_v, s_1],
    ];
    let rhs = [s_ub, s_vb, s_b];
    solve_3x3_system(&a, &rhs)
}

fn solve_3x3_system(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Result<(f64, f64, f64)> {
    let mut m = [[0.0_f64; 4]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = a[i][j];
        }
        m[i][3] = b[i];
    }

    for col in 0..3 {
        let mut pivot_row = col;
        let mut pivot_val = libm::fabs(m[col][col]);
        for (row, mr) in m.iter().enumerate().skip(col + 1) {
            let v = libm::fabs(mr[col]);
            if v > pivot_val {
                pivot_val = v;
                pivot_row = row;
            }
        }
        ensure!(pivot_val > 1e-30, "singular normal equations in 3x3");
        if pivot_row != col {
            m.swap(col, pivot_row);
        }
        let diag = m[col][col];
        let col_snapshot = m[col];
        for mr in m.iter_mut().skip(col + 1) {
            let factor = mr[col] / diag;
            for (x, &mc) in mr[col..].iter_mut().zip(col_snapshot[col..].iter()) {
                *x -= factor * mc;
            }
        }
    }

    let z = m[2][3] / m[2][2];
    let y = (m[1][3] - m[1][2] * z) / m[1][1];
    let x = (m[0][3] - m[0][1] * y - m[0][2] * z) / m[0][0];
    Ok((x, y, z))
}

pub(super) fn compute_residuals(pairs: &[&StarPair], wcs: WcsParams) -> Vec<StarResidual> {
    let det = wcs.cd1_1 * wcs.cd2_2 - wcs.cd1_2 * wcs.cd2_1;
    if libm::fabs(det) < 1e-30 {
        return Vec::new();
    }
    let inv_00 = wcs.cd2_2 / det;
    let inv_01 = -wcs.cd1_2 / det;
    let inv_10 = -wcs.cd2_1 / det;
    let inv_11 = wcs.cd1_1 / det;

    pairs
        .iter()
        .filter_map(|p| {
            let (xi_rad, eta_rad) = tan_project_star(
                p.ra_deg, p.dec_deg, wcs.crval_ra, wcs.crval_dec,
            )?;
            let xi = xi_rad * RAD_TO_DEG;
            let eta = eta_rad * RAD_TO_DEG;
            let pred_x = inv_00 * xi + inv_01 * eta + wcs.crpix1;
            let pred_y = inv_10 * xi + inv_11 * eta + wcs.crpix2;
            let err_x = p.px_x - pred_x;
            let err_y = p.px_y - pred_y;
            let err_px = libm::sqrt(err_x * err_x + err_y * err_y);
            Some(StarResidual {
                px_x: p.px_x, px_y: p.px_y,
                err_x, err_y, err_px,
            })
        })
        .collect()
}

pub(super) fn rms_error(residuals: &[StarResidual]) -> f64 {
    if residuals.is_empty() {
        return 0.0;
    }
    let sum_sq: f64 = residuals.iter().map(|r| r.err_px * r.err_px).sum();
    libm::sqrt(sum_sq / residuals.len() as f64)
}

pub(super) fn weighted_rms_error(residuals: &[StarResidual], weights: &[f64]) -> f64 {
    if residuals.is_empty() {
        return 0.0;
    }
    let mut sum_w = 0.0_f64;
    let mut sum_we2 = 0.0_f64;
    for (r, &w) in residuals.iter().zip(weights.iter()) {
        sum_w += w;
        sum_we2 += w * r.err_px * r.err_px;
    }
    if sum_w > 0.0 { libm::sqrt(sum_we2 / sum_w) } else { 0.0 }
}

pub(super) fn solve_normal_equations(
    a: &mut [f64],
    b: &mut [f64],
    n: usize,
) -> Result<Vec<f64>> {
    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_val = libm::fabs(a[col * n + col]);
        for row in (col + 1)..n {
            let v = libm::fabs(a[row * n + col]);
            if v > pivot_val {
                pivot_val = v;
                pivot_row = row;
            }
        }
        ensure!(pivot_val > 1e-30, "singular SIP normal equations at col {}", col);
        if pivot_row != col {
            for k in 0..n {
                a.swap(col * n + k, pivot_row * n + k);
            }
            b.swap(col, pivot_row);
        }
        let diag = a[col * n + col];
        for row in (col + 1)..n {
            let factor = a[row * n + col] / diag;
            for k in col..n {
                let above = a[col * n + k];
                a[row * n + k] -= factor * above;
            }
            b[row] -= factor * b[col];
        }
    }

    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in (i + 1)..n {
            sum -= a[i * n + j] * x[j];
        }
        x[i] = sum / a[i * n + i];
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn residual(err_x: f64, err_y: f64) -> StarResidual {
        let err_px = libm::sqrt(err_x * err_x + err_y * err_y);
        StarResidual { px_x: 0.0, px_y: 0.0, err_x, err_y, err_px }
    }

    fn pair(px_x: f64, px_y: f64, ra_deg: f64, dec_deg: f64) -> StarPair {
        StarPair { px_x, px_y, ra_deg, dec_deg, votes: 4, snr: 20.0 }
    }

    #[test]
    fn rms_error_empty_returns_zero() {
        assert_eq!(rms_error(&[]), 0.0);
    }

    #[test]
    fn rms_error_matches_hand_calculation() {
        // residuals with err_px 3, 4 -> sqrt((9+16)/2) = sqrt(12.5)
        let residuals = vec![residual(3.0, 0.0), residual(0.0, 4.0)];
        let rms = rms_error(&residuals);
        assert_eq!(rms, libm::sqrt(12.5));
    }

    #[test]
    fn weighted_rms_error_empty_returns_zero() {
        assert_eq!(weighted_rms_error(&[], &[]), 0.0);
    }

    #[test]
    fn weighted_rms_error_zero_weight_returns_zero() {
        let residuals = vec![residual(3.0, 0.0)];
        let weights = vec![0.0];
        assert_eq!(weighted_rms_error(&residuals, &weights), 0.0);
    }

    #[test]
    fn weighted_rms_error_shifts_toward_heavily_weighted_residuals() {
        // equal weights: rms = sqrt((9+16)/2) = sqrt(12.5)
        let residuals = vec![residual(3.0, 0.0), residual(0.0, 4.0)];
        let equal = weighted_rms_error(&residuals, &[1.0, 1.0]);
        assert_eq!(equal, libm::sqrt(12.5));

        // heavier on the bigger residual -> rms pulled up
        let skewed_big = weighted_rms_error(&residuals, &[1.0, 100.0]);
        assert!(skewed_big > equal);

        // heavier on the smaller residual -> rms pulled down
        let skewed_small = weighted_rms_error(&residuals, &[100.0, 1.0]);
        assert!(skewed_small < equal);
    }

    #[test]
    fn solve_normal_equations_identity_returns_b() {
        // 3x3 identity matrix, b = [1, 2, 3]
        let mut a = vec![
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
        ];
        let mut b = vec![1.0, 2.0, 3.0];
        let x = solve_normal_equations(&mut a, &mut b, 3).unwrap();
        assert_eq!(x, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn solve_normal_equations_diagonal_returns_componentwise_divide() {
        // diag(2, 4, 8), b = [2, 8, 24] -> x = [1, 2, 3]
        let mut a = vec![
            2.0, 0.0, 0.0,
            0.0, 4.0, 0.0,
            0.0, 0.0, 8.0,
        ];
        let mut b = vec![2.0, 8.0, 24.0];
        let x = solve_normal_equations(&mut a, &mut b, 3).unwrap();
        assert_eq!(x, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn solve_normal_equations_needs_pivoting_still_works() {
        // first-column pivot is zero, so row swap is required
        let mut a = vec![
            0.0, 1.0,
            2.0, 0.0,
        ];
        let mut b = vec![3.0, 4.0];
        // x_0 = 4/2 = 2, x_1 = 3/1 = 3
        let x = solve_normal_equations(&mut a, &mut b, 2).unwrap();
        assert_eq!(x, vec![2.0, 3.0]);
    }

    #[test]
    fn solve_normal_equations_singular_errors() {
        // rank-1 matrix: row 1 = 2x row 0
        let mut a = vec![
            1.0, 2.0,
            2.0, 4.0,
        ];
        let mut b = vec![1.0, 2.0];
        assert!(solve_normal_equations(&mut a, &mut b, 2).is_err());
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
        assert!(residuals[0].err_px < 1e-12);
    }

    #[test]
    fn compute_residuals_returns_empty_for_singular_cd() {
        let wcs = WcsParams {
            crpix1: 0.0, crpix2: 0.0, crval_ra: 0.0, crval_dec: 0.0,
            cd1_1: 0.0, cd1_2: 0.0, cd2_1: 0.0, cd2_2: 0.0,
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
        use celestial_catalog::query::tan_deproject_star;
        use celestial_core::constants::DEG_TO_RAD;

        samples.iter().map(|&(u, v)| {
            let xi_deg = cd.cd1_1 * u + cd.cd1_2 * v + cd.xi0;
            let eta_deg = cd.cd2_1 * u + cd.cd2_2 * v + cd.eta0;
            let (ra, dec) = tan_deproject_star(
                xi_deg * DEG_TO_RAD, eta_deg * DEG_TO_RAD, crval_ra, crval_dec,
            );
            pair(u, v, ra, dec)
        }).collect()
    }

    #[test]
    fn fit_cd_weighted_recovers_known_cd_matrix() {
        let crval_ra = 180.0;
        let crval_dec = 0.0;
        let cd = CdParams {
            cd1_1: 0.001, cd1_2: 0.0, xi0: 0.0,
            cd2_1: 0.0, cd2_2: 0.001, eta0: 0.0,
        };
        let samples = [
            (-100.0, -100.0), (100.0, -100.0), (0.0, 0.0),
            (-100.0, 100.0), (100.0, 100.0), (50.0, -50.0),
        ];
        let pairs = synthesize_pairs(crval_ra, crval_dec, &cd, &samples);

        // Project through the "known" WCS to get (u, v, xi, eta)
        use celestial_catalog::query::tan_project_star;
        use celestial_core::constants::RAD_TO_DEG;
        let projected: Vec<(f64, f64, f64, f64)> = pairs.iter().filter_map(|p| {
            tan_project_star(p.ra_deg, p.dec_deg, crval_ra, crval_dec)
                .map(|(xi, eta)| (p.px_x, p.px_y, xi * RAD_TO_DEG, eta * RAD_TO_DEG))
        }).collect();

        let mask = vec![true; projected.len()];
        let weights = vec![1.0; projected.len()];
        let fit = fit_cd_weighted(&projected, &mask, &weights).unwrap();

        // numerical round trip through tan projection; should recover within
        // a very tight tolerance
        assert!((fit.cd1_1 - 0.001).abs() < 1e-12);
        assert!((fit.cd2_2 - 0.001).abs() < 1e-12);
        assert!(fit.cd1_2.abs() < 1e-12);
        assert!(fit.cd2_1.abs() < 1e-12);
    }

    #[test]
    fn sigma_clip_loop_removes_obvious_outlier() {
        let crval_ra = 180.0;
        let crval_dec = 0.0;
        let cd = CdParams {
            cd1_1: 0.001, cd1_2: 0.0, xi0: 0.0,
            cd2_1: 0.0, cd2_2: 0.001, eta0: 0.0,
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

        use celestial_catalog::query::tan_project_star;
        use celestial_core::constants::RAD_TO_DEG;
        let mut projected: Vec<(f64, f64, f64, f64)> = pairs.iter().filter_map(|p| {
            tan_project_star(p.ra_deg, p.dec_deg, crval_ra, crval_dec)
                .map(|(xi, eta)| (p.px_x, p.px_y, xi * RAD_TO_DEG, eta * RAD_TO_DEG))
        }).collect();

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

        sigma_clip_loop_weighted(&projected, &mut mask, &mut cd_fit, &weights).unwrap();
        // outlier must be flagged
        assert!(!mask[last], "expected outlier to be rejected");
        // clean points must be kept
        assert!(mask[..last].iter().all(|&m| m), "clean points should survive");
    }
}
