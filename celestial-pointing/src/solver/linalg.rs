use crate::error::{Error, Result};
use nalgebra::{DMatrix, DVector};

pub(super) struct SolveOutcome {
    pub(super) x: DVector<f64>,
    pub(super) cov_diag: Vec<f64>,
    pub(super) leverage: Vec<f64>,
    pub(super) rank: usize,
    pub(super) sigma_max: f64,
    pub(super) sigma_min: f64,
}

pub(super) fn solve_damped(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    lambda: f64,
    fit_tol: f64,
) -> Result<SolveOutcome> {
    if lambda <= 0.0 {
        return solve_and_covariance(a, b, fit_tol);
    }
    let m = a.ncols();
    let rows = a.nrows();
    let mut a_aug = DMatrix::zeros(rows + m, m);
    let mut b_aug = DVector::zeros(rows + m);
    for r in 0..rows {
        for c in 0..m {
            a_aug[(r, c)] = a[(r, c)];
        }
        b_aug[r] = b[r];
    }
    let sqrt_lambda = lambda.sqrt();
    for c in 0..m {
        a_aug[(rows + c, c)] = sqrt_lambda;
    }
    solve_and_covariance(&a_aug, &b_aug, fit_tol)
}

pub(super) fn solve_and_covariance(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    fit_tol: f64,
) -> Result<SolveOutcome> {
    let m = a.ncols();
    let mut col_scale = vec![1.0_f64; m];
    let mut a_scaled = a.clone();
    for j in 0..m {
        let norm = a.column(j).norm();
        if norm > 0.0 {
            col_scale[j] = 1.0 / norm;
            for i in 0..a_scaled.nrows() {
                a_scaled[(i, j)] *= col_scale[j];
            }
        }
    }

    let svd = a_scaled.svd(true, true);
    let u = svd
        .u
        .as_ref()
        .ok_or_else(|| Error::Fit("SVD did not produce U".into()))?;
    let v_t = svd
        .v_t
        .as_ref()
        .ok_or_else(|| Error::Fit("SVD did not produce V^T".into()))?;
    let singular = &svd.singular_values;
    let k = singular.len();
    let sigma_max = if k > 0 { singular[0] } else { 0.0 };
    let sigma_min = if k > 0 { singular[k - 1] } else { 0.0 };

    let threshold = if fit_tol > 0.0 {
        fit_tol * sigma_max
    } else {
        0.0
    };

    let mut inv_s = vec![0.0_f64; k];
    let mut rank = 0;
    for (i, &s) in singular.iter().enumerate() {
        if s > threshold {
            inv_s[i] = 1.0 / s;
            rank += 1;
        }
    }

    let u_t_b = u.transpose() * b;
    let mut z = DVector::zeros(k);
    for i in 0..k {
        z[i] = inv_s[i] * u_t_b[i];
    }
    let v = v_t.transpose();
    let x_scaled = &v * z;
    let mut x = DVector::zeros(m);
    for j in 0..m {
        x[j] = x_scaled[j] * col_scale[j];
    }

    let mut cov_diag = vec![0.0_f64; m];
    for j in 0..m {
        let mut acc = 0.0_f64;
        for i in 0..k {
            let inv = inv_s[i];
            if inv == 0.0 {
                continue;
            }
            let v_ji = v[(j, i)];
            acc += (v_ji * inv) * (v_ji * inv);
        }
        cov_diag[j] = acc * col_scale[j] * col_scale[j];
    }

    let n_rows = a.nrows();
    let mut leverage = vec![0.0_f64; n_rows];
    for r in 0..n_rows {
        let mut acc = 0.0_f64;
        for c in 0..k {
            if inv_s[c] == 0.0 {
                continue;
            }
            let u_rc = u[(r, c)];
            acc += u_rc * u_rc;
        }
        leverage[r] = acc;
    }

    Ok(SolveOutcome {
        x,
        cov_diag,
        leverage,
        rank,
        sigma_max,
        sigma_min,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // Identity system: A = I, b = [1, 2, 3] → x = [1, 2, 3], rank = 3.
    #[test]
    fn solve_and_covariance_identity_system_recovers_b() {
        let a = DMatrix::identity(3, 3);
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let out = solve_and_covariance(&a, &b, 1.0e-9).unwrap();
        assert_eq!(out.rank, 3);
        for i in 0..3 {
            assert!((out.x[i] - b[i]).abs() < 1e-9);
        }
        // Leverage = 1 per row for a square full-rank identity.
        for r in 0..3 {
            assert!((out.leverage[r] - 1.0).abs() < 1e-9);
        }
        assert!((out.sigma_max - 1.0).abs() < 1e-9);
        assert!((out.sigma_min - 1.0).abs() < 1e-9);
    }

    // Diagonal scaling system: a non-trivial column-scaling case. The
    // function normalizes columns to unit norm before SVD, then unscales.
    #[test]
    fn solve_and_covariance_handles_disparate_column_scales() {
        let mut a = DMatrix::zeros(3, 2);
        a[(0, 0)] = 1.0;
        a[(1, 0)] = 1.0;
        a[(2, 0)] = 1.0;
        a[(0, 1)] = 1000.0;
        a[(1, 1)] = 2000.0;
        a[(2, 1)] = 3000.0;
        // True solution: x0 = 1, x1 = 0.001 → b = [2, 3, 4]
        let b = DVector::from_vec(vec![2.0, 3.0, 4.0]);
        let out = solve_and_covariance(&a, &b, 1.0e-9).unwrap();
        assert!((out.x[0] - 1.0).abs() < 1e-6);
        assert!((out.x[1] - 0.001).abs() < 1e-9);
    }

    // Rank deficiency: A has two identical columns. With reasonable tolerance,
    // rank should be 1, not 2 — the SVD truncation guards against the
    // degeneracy.
    #[test]
    fn solve_and_covariance_truncates_rank_for_collinear_columns() {
        let mut a = DMatrix::zeros(3, 2);
        for r in 0..3 {
            a[(r, 0)] = 1.0;
            a[(r, 1)] = 1.0; // identical
        }
        let b = DVector::from_vec(vec![1.0, 1.0, 1.0]);
        let out = solve_and_covariance(&a, &b, 1.0e-6).unwrap();
        assert_eq!(out.rank, 1, "collinear columns should produce rank 1");
    }

    #[test]
    fn solve_and_covariance_zero_tolerance_keeps_full_rank() {
        // With fit_tol = 0, threshold = 0, all positive singular values count.
        let a = DMatrix::identity(3, 3);
        let b = DVector::from_vec(vec![1.0, 0.0, 0.0]);
        let out = solve_and_covariance(&a, &b, 0.0).unwrap();
        assert_eq!(out.rank, 3);
    }

    #[test]
    fn solve_and_covariance_zero_column_does_not_panic() {
        // Column 1 is all zeros. col_scale stays at 1.0 for that column;
        // the resulting singular value is 0 and gets truncated.
        let mut a = DMatrix::zeros(3, 2);
        a[(0, 0)] = 1.0;
        a[(1, 0)] = 1.0;
        a[(2, 0)] = 1.0;
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let out = solve_and_covariance(&a, &b, 1.0e-6).unwrap();
        assert_eq!(out.rank, 1);
        // x for the zero-column term should be exactly 0.
        assert_eq!(out.x[1], 0.0);
    }

    // Covariance diagonal is non-negative (variance interpretation).
    #[test]
    fn solve_and_covariance_cov_diag_is_nonnegative() {
        let mut a = DMatrix::zeros(4, 2);
        a[(0, 0)] = 1.0;
        a[(1, 0)] = 2.0;
        a[(2, 1)] = 1.0;
        a[(3, 1)] = 3.0;
        let b = DVector::from_vec(vec![1.0, 1.0, 1.0, 1.0]);
        let out = solve_and_covariance(&a, &b, 1.0e-9).unwrap();
        for &c in &out.cov_diag {
            assert!(c >= 0.0, "covariance diagonal must be non-negative");
        }
    }

    #[test]
    fn solve_and_covariance_leverage_sums_bounded_by_rank() {
        // For a full-rank n×m system with n >= m, leverages should sum to rank.
        let a = DMatrix::identity(5, 5);
        let b = DVector::from_vec(vec![1.0; 5]);
        let out = solve_and_covariance(&a, &b, 1.0e-9).unwrap();
        let total: f64 = out.leverage.iter().sum();
        assert!((total - out.rank as f64).abs() < 1e-6);
    }

    // --- solve_damped ------------------------------------------------

    #[test]
    fn solve_damped_zero_lambda_delegates_to_undamped_solve() {
        let a = DMatrix::identity(3, 3);
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let undamped = solve_and_covariance(&a, &b, 1.0e-9).unwrap();
        let damped = solve_damped(&a, &b, 0.0, 1.0e-9).unwrap();
        for i in 0..3 {
            assert!((undamped.x[i] - damped.x[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn solve_damped_negative_lambda_also_delegates_undamped() {
        let a = DMatrix::identity(2, 2);
        let b = DVector::from_vec(vec![1.0, 1.0]);
        let damped = solve_damped(&a, &b, -1.0, 1.0e-9).unwrap();
        assert!((damped.x[0] - 1.0).abs() < 1e-12);
        assert!((damped.x[1] - 1.0).abs() < 1e-12);
    }

    // Heavy damping pulls the solution toward zero — large lambda should
    // shrink the recovered coefficients.
    #[test]
    fn solve_damped_large_lambda_shrinks_solution_toward_zero() {
        let a = DMatrix::identity(2, 2);
        let b = DVector::from_vec(vec![1.0, 1.0]);
        let lightly_damped = solve_damped(&a, &b, 0.01, 1.0e-9).unwrap();
        let heavily_damped = solve_damped(&a, &b, 100.0, 1.0e-9).unwrap();
        assert!(heavily_damped.x[0].abs() < lightly_damped.x[0].abs());
        assert!(heavily_damped.x[0].abs() < 0.5);
    }

    // --- SolveOutcome accessors --------------------------------------

    #[test]
    fn solve_outcome_fields_are_visible_within_super() {
        // Compile-time check that pub(super) fields are accessible.
        let a = DMatrix::identity(2, 2);
        let b = DVector::from_vec(vec![1.0, 1.0]);
        let out = solve_and_covariance(&a, &b, 1.0e-9).unwrap();
        assert_eq!(out.x.len(), 2);
        assert_eq!(out.cov_diag.len(), 2);
        assert_eq!(out.leverage.len(), 2);
        assert!(out.rank > 0);
        assert!(out.sigma_max >= out.sigma_min);
    }
}
