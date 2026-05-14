use super::linalg::SolveOutcome;
use super::ObsDiagnostic;
use crate::observation::Observation;
use nalgebra::DVector;

pub(crate) fn compute_sky_rms(residuals: &DVector<f64>, observations: &[&Observation]) -> f64 {
    let n = observations.len();
    if n == 0 {
        return 0.0;
    }
    libm::sqrt(sum_sq_sky(residuals, observations) / n as f64)
}

pub(super) fn compute_popn_sd(
    residuals: &DVector<f64>,
    observations: &[&Observation],
    n_free_params: usize,
) -> f64 {
    let n = observations.len();
    if n == 0 || n <= n_free_params {
        return 0.0;
    }
    let dof = n - n_free_params;
    libm::sqrt(sum_sq_sky(residuals, observations) / dof as f64)
}

fn sum_sq_sky(residuals: &DVector<f64>, observations: &[&Observation]) -> f64 {
    let mut sum_sq = 0.0;
    for i in 0..observations.len() {
        let dh = residuals[2 * i];
        let dd = residuals[2 * i + 1];
        let cos_dec = libm::cos(observations[i].catalog_dec.radians());
        let dx = dh * cos_dec;
        sum_sq += dx * dx + dd * dd;
    }
    sum_sq
}

pub(super) fn compute_diagnostics(
    free_residuals: &DVector<f64>,
    leverage: &[f64],
    rank: usize,
) -> Vec<ObsDiagnostic> {
    let n_rows = free_residuals.len();
    let n_obs = n_rows / 2;
    let rss = free_residuals.dot(free_residuals);
    let dof = n_rows.saturating_sub(rank).max(1);
    let s2 = rss / dof as f64;
    let s = libm::sqrt(s2);

    let mut out = Vec::with_capacity(n_obs);
    for i in 0..n_obs {
        let r_ha = free_residuals[2 * i];
        let r_dec = free_residuals[2 * i + 1];
        let h_ha = leverage[2 * i].clamp(0.0, 1.0);
        let h_dec = leverage[2 * i + 1].clamp(0.0, 1.0);
        let residual_sky = libm::sqrt(r_ha * r_ha + r_dec * r_dec);
        let t_ha = if s > 0.0 && h_ha < 1.0 {
            r_ha / (s * libm::sqrt(1.0 - h_ha))
        } else {
            0.0
        };
        let t_dec = if s > 0.0 && h_dec < 1.0 {
            r_dec / (s * libm::sqrt(1.0 - h_dec))
        } else {
            0.0
        };
        let studentized = libm::sqrt((t_ha * t_ha + t_dec * t_dec) / 2.0);
        out.push(ObsDiagnostic {
            residual_sky,
            leverage: h_ha + h_dec,
            studentized,
        });
    }
    out
}

pub(super) fn compute_sigma_from_covariance(
    solved: &SolveOutcome,
    residuals: &DVector<f64>,
    free_indices: &[usize],
    total_terms: usize,
    n_observations: usize,
) -> Vec<f64> {
    let dof = n_observations.saturating_sub(solved.rank).max(1);
    let s2 = residuals.dot(residuals) / dof as f64;
    let mut sigma = vec![0.0_f64; total_terms];
    for (fi, &idx) in free_indices.iter().enumerate() {
        sigma[idx] = libm::sqrt((s2 * solved.cov_diag[fi]).abs());
    }
    sigma
}

#[cfg(test)]
mod tests {
    use super::super::linalg::{solve_and_covariance, SolveOutcome};
    use super::*;
    use crate::test_support::ObsBuilder;
    use nalgebra::DMatrix;

    fn obs(dec_deg: f64) -> Observation {
        ObsBuilder::new()
            .catalog_dec_deg(dec_deg)
            .observed_dec_deg(dec_deg)
            .build()
    }

    // --- compute_sky_rms ----------------------------------------------

    #[test]
    fn compute_sky_rms_empty_observations_returns_zero() {
        let residuals = DVector::from_vec(vec![]);
        let observations: Vec<&Observation> = Vec::new();
        assert_eq!(compute_sky_rms(&residuals, &observations), 0.0);
    }

    #[test]
    fn compute_sky_rms_uses_cos_dec_weighting_on_ha_residuals() {
        // Single obs at dec=60° (cos=0.5). HA residual = 100", Dec residual = 0".
        // dx = 100 * 0.5 = 50; sky_rms = sqrt(50^2 / 1) = 50.
        let o = obs(60.0);
        let observations: Vec<&Observation> = vec![&o];
        let residuals = DVector::from_vec(vec![100.0, 0.0]);
        let rms = compute_sky_rms(&residuals, &observations);
        assert!((rms - 50.0).abs() < 1e-6);
    }

    #[test]
    fn compute_sky_rms_at_equator_no_cos_dec_attenuation() {
        let o = obs(0.0); // cos(0) = 1
        let observations: Vec<&Observation> = vec![&o];
        let residuals = DVector::from_vec(vec![3.0, 4.0]);
        // dx=3, dd=4 → sqrt((9+16)/1) = 5
        let rms = compute_sky_rms(&residuals, &observations);
        assert!((rms - 5.0).abs() < 1e-9);
    }

    #[test]
    fn compute_sky_rms_averages_across_observations() {
        let o1 = obs(0.0);
        let o2 = obs(0.0);
        let observations: Vec<&Observation> = vec![&o1, &o2];
        // Two obs: dx=3,dd=4 (norm 5); dx=0,dd=0. RMS = sqrt(25 / 2) ≈ 3.535
        let residuals = DVector::from_vec(vec![3.0, 4.0, 0.0, 0.0]);
        let rms = compute_sky_rms(&residuals, &observations);
        assert!((rms - libm::sqrt(12.5)).abs() < 1e-9);
    }

    // --- compute_popn_sd ----------------------------------------------

    #[test]
    fn compute_popn_sd_empty_returns_zero() {
        let observations: Vec<&Observation> = Vec::new();
        let residuals = DVector::from_vec(vec![]);
        assert_eq!(compute_popn_sd(&residuals, &observations, 0), 0.0);
    }

    #[test]
    fn compute_popn_sd_returns_zero_when_n_le_n_free_params() {
        let o = obs(0.0);
        let observations: Vec<&Observation> = vec![&o, &o];
        let residuals = DVector::from_vec(vec![1.0, 1.0, 1.0, 1.0]);
        // 2 obs, 2 free params → dof would be 0 → returns 0
        assert_eq!(compute_popn_sd(&residuals, &observations, 2), 0.0);
        // 2 obs, 3 free params → returns 0
        assert_eq!(compute_popn_sd(&residuals, &observations, 3), 0.0);
    }

    // The denominator is (n - n_free_params), so popn_sd > sky_rms whenever
    // n_free_params > 0.
    #[test]
    fn compute_popn_sd_is_larger_than_sky_rms_when_dof_consumed() {
        let o1 = obs(0.0);
        let o2 = obs(0.0);
        let o3 = obs(0.0);
        let observations: Vec<&Observation> = vec![&o1, &o2, &o3];
        let residuals = DVector::from_vec(vec![3.0, 4.0, 3.0, 4.0, 3.0, 4.0]);
        let rms = compute_sky_rms(&residuals, &observations);
        let sd = compute_popn_sd(&residuals, &observations, 1);
        assert!(sd > rms, "popn_sd ({}) should exceed sky_rms ({})", sd, rms);
    }

    #[test]
    fn compute_popn_sd_with_zero_free_params_equals_sky_rms() {
        let o = obs(0.0);
        let observations: Vec<&Observation> = vec![&o, &o];
        let residuals = DVector::from_vec(vec![3.0, 4.0, 3.0, 4.0]);
        let rms = compute_sky_rms(&residuals, &observations);
        let sd = compute_popn_sd(&residuals, &observations, 0);
        assert!((sd - rms).abs() < 1e-9);
    }

    // --- compute_diagnostics ------------------------------------------

    #[test]
    fn compute_diagnostics_produces_one_entry_per_observation() {
        // 4 free residuals = 2 observations.
        let resid = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let leverage = vec![0.1, 0.2, 0.3, 0.4];
        let out = compute_diagnostics(&resid, &leverage, 1);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn compute_diagnostics_residual_sky_is_l2_norm_of_pair() {
        let resid = DVector::from_vec(vec![3.0, 4.0]); // single obs
        let leverage = vec![0.0, 0.0];
        let out = compute_diagnostics(&resid, &leverage, 1);
        assert!((out[0].residual_sky - 5.0).abs() < 1e-9);
    }

    #[test]
    fn compute_diagnostics_leverage_sum_is_clamped_per_row() {
        let resid = DVector::from_vec(vec![0.0, 0.0]);
        // Out-of-range leverages clamp to [0, 1].
        let leverage = vec![-0.5, 1.5];
        let out = compute_diagnostics(&resid, &leverage, 1);
        assert_eq!(out[0].leverage, 0.0 + 1.0);
    }

    // Leverage at 1.0 makes the denominator (1 - h_ha) = 0 → studentized
    // term defaults to 0 to avoid div-by-zero.
    #[test]
    fn compute_diagnostics_leverage_one_gives_zero_studentized() {
        let resid = DVector::from_vec(vec![100.0, 100.0]);
        let leverage = vec![1.0, 1.0];
        let out = compute_diagnostics(&resid, &leverage, 1);
        assert_eq!(out[0].studentized, 0.0);
    }

    // Zero residuals across the board produce zero RSS → s = 0 →
    // studentized = 0 (the s > 0 guard).
    #[test]
    fn compute_diagnostics_zero_residuals_produce_zero_studentized() {
        let resid = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0]);
        let leverage = vec![0.1, 0.1, 0.1, 0.1];
        let out = compute_diagnostics(&resid, &leverage, 2);
        for d in &out {
            assert_eq!(d.studentized, 0.0);
            assert_eq!(d.residual_sky, 0.0);
        }
    }

    #[test]
    fn compute_diagnostics_full_rank_dof_falls_back_to_one() {
        // n_rows = 2, rank = 2 → dof would be 0; saturating_sub guards it.
        let resid = DVector::from_vec(vec![1.0, 1.0]);
        let leverage = vec![0.0, 0.0];
        let out = compute_diagnostics(&resid, &leverage, 10); // rank > n_rows
        // Should not panic, should still produce one diag.
        assert_eq!(out.len(), 1);
    }

    // --- compute_sigma_from_covariance --------------------------------

    fn make_outcome() -> SolveOutcome {
        let a = DMatrix::identity(3, 3);
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        solve_and_covariance(&a, &b, 1.0e-9).unwrap()
    }

    #[test]
    fn compute_sigma_zero_indices_yields_zeros() {
        let solved = make_outcome();
        let residuals = DVector::from_vec(vec![1.0, 1.0, 1.0]);
        let sigma = compute_sigma_from_covariance(&solved, &residuals, &[], 5, 3);
        assert_eq!(sigma, vec![0.0; 5]);
    }

    #[test]
    fn compute_sigma_maps_free_indices_into_total_terms_layout() {
        let solved = make_outcome();
        let residuals = DVector::from_vec(vec![1.0, 1.0, 1.0]);
        // free indices [0, 2] in a 4-term model → sigma[1] and sigma[3] should be 0.
        let sigma = compute_sigma_from_covariance(&solved, &residuals, &[0, 2], 4, 3);
        assert_eq!(sigma.len(), 4);
        assert!(sigma[0] > 0.0);
        assert_eq!(sigma[1], 0.0);
        assert!(sigma[2] > 0.0);
        assert_eq!(sigma[3], 0.0);
    }

    // The dof denominator saturates at 1 if rank >= n_observations, so it
    // never panics. (This is an unusual but reachable code path — e.g.,
    // when the problem is exactly determined.)
    #[test]
    fn compute_sigma_saturates_dof_at_one() {
        let solved = make_outcome();
        let residuals = DVector::from_vec(vec![0.1, 0.1, 0.1]);
        // n_observations = solved.rank = 3 → dof would be 0, saturates to 1.
        let sigma = compute_sigma_from_covariance(&solved, &residuals, &[0], 1, 3);
        assert!(sigma[0].is_finite());
    }
}
