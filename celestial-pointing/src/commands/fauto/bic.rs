use crate::error::Result;
use crate::observation::Observation;
use crate::solver::{fit_model, FitResult};
use crate::terms::create_term;

pub(super) const MAX_COND_NUMBER: f64 = 1.0e4;

pub(super) fn compute_bic(n_obs: usize, n_terms: usize, sky_rms: f64) -> f64 {
    let n = n_obs as f64;
    let k = n_terms as f64;
    let weighted_rss = sky_rms * sky_rms * n;
    n * libm::log(weighted_rss / n) + k * libm::log(n)
}

pub(super) fn try_fit(
    observations: &[&Observation],
    term_names: &[String],
    latitude: f64,
    fit_tol: f64,
) -> Result<FitResult> {
    let terms: Vec<_> = term_names
        .iter()
        .map(|n| create_term(n))
        .collect::<Result<Vec<_>>>()?;
    let fixed = vec![false; terms.len()];
    let coeffs = vec![0.0; terms.len()];
    fit_model(observations, &terms, &fixed, &coeffs, latitude, fit_tol)
}

pub(super) fn fit_is_well_conditioned(fit: &FitResult) -> bool {
    !fit.rank_info.is_deficient() && fit.rank_info.condition_number() <= MAX_COND_NUMBER
}

pub(super) fn bic_of(
    observations: &[&Observation],
    active: &[String],
    latitude: f64,
    fit_tol: f64,
) -> Result<f64> {
    let fit = try_fit(observations, active, latitude, fit_tol)?;
    Ok(compute_bic(observations.len(), active.len(), fit.sky_rms))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::RankInfo;
    use crate::test_support::FitResultBuilder;

    fn fit_with_rank(n_free: usize, rank: usize, sigma_max: f64, sigma_min: f64) -> FitResult {
        FitResultBuilder::new()
            .rank_info(RankInfo {
                n_free,
                rank,
                sigma_max,
                sigma_min,
                tol: 0.0,
            })
            .build()
    }

    // BIC reference values: for N observations and one term with RMS=R,
    // BIC = N * ln(R^2) + ln(N).
    #[test]
    fn compute_bic_matches_hand_calculation() {
        let n = 100_usize;
        let rms = 2.0_f64;
        let expected = (n as f64) * libm::log(rms * rms) + libm::log(n as f64);
        let got = compute_bic(n, 1, rms);
        assert!((got - expected).abs() < 1e-9, "got {}, expected {}", got, expected);
    }

    #[test]
    fn compute_bic_increases_with_term_count() {
        let n = 50;
        let rms = 1.5;
        let bic_1 = compute_bic(n, 1, rms);
        let bic_5 = compute_bic(n, 5, rms);
        assert!(bic_5 > bic_1, "more terms should cost more BIC");
    }

    #[test]
    fn compute_bic_decreases_with_lower_rms() {
        let n = 50;
        let lower = compute_bic(n, 3, 0.5);
        let higher = compute_bic(n, 3, 2.0);
        assert!(lower < higher, "smaller RMS should yield smaller BIC");
    }

    #[test]
    fn fit_is_well_conditioned_accepts_full_rank_small_cond() {
        // cond = sigma_max / sigma_min = 10 / 1 = 10 << 1e4
        let fit = fit_with_rank(3, 3, 10.0, 1.0);
        assert!(fit_is_well_conditioned(&fit));
    }

    #[test]
    fn fit_is_well_conditioned_rejects_rank_deficient() {
        // rank < n_free → deficient
        let fit = fit_with_rank(3, 2, 10.0, 1.0);
        assert!(!fit_is_well_conditioned(&fit));
    }

    #[test]
    fn fit_is_well_conditioned_rejects_huge_condition_number() {
        // cond = 1e5 / 1.0 = 1e5 > MAX_COND_NUMBER (1e4)
        let fit = fit_with_rank(3, 3, 1.0e5, 1.0);
        assert!(!fit_is_well_conditioned(&fit));
    }

    #[test]
    fn fit_is_well_conditioned_accepts_at_max_cond_number() {
        // cond exactly at the limit
        let fit = fit_with_rank(3, 3, MAX_COND_NUMBER, 1.0);
        assert!(fit_is_well_conditioned(&fit));
    }
}
