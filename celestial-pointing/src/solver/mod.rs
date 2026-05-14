mod design;
mod iterated;
mod linalg;
mod single_step;
mod stats;
mod validate;

use crate::error::Result;
use crate::observation::Observation;
use crate::terms::Term;

pub(crate) use design::build_residuals;
pub(crate) use iterated::fit_model_iterated;
pub(crate) use stats::compute_sky_rms;

use single_step::fit_model_with_inputs;

#[derive(Clone)]
pub struct FitResult {
    pub coefficients: Vec<f64>,
    pub change: Vec<f64>,
    pub sigma: Vec<f64>,
    pub sky_rms: f64,
    pub popn_sd: f64,
    pub term_names: Vec<String>,
    pub rank_info: RankInfo,
    pub diagnostics: Vec<ObsDiagnostic>,
    pub iter_report: IterReport,
}

#[derive(Clone, Debug, Default)]
pub struct IterReport {
    pub iterations: usize,
    pub converged: bool,
    pub final_lambda: f64,
    pub final_step_norm: f64,
}

#[derive(Clone, Debug)]
pub struct IterOptions {
    pub max_iter: usize,
    pub tol_rad: f64,
    pub initial_lambda: f64,
    pub lambda_up: f64,
    pub lambda_down: f64,
    pub robust: bool,
    pub huber_k: f64,
}

impl Default for IterOptions {
    fn default() -> Self {
        Self {
            max_iter: 50,
            tol_rad: 1.0e-6,
            initial_lambda: 1.0e-3,
            lambda_up: 10.0,
            lambda_down: 10.0,
            robust: false,
            huber_k: 1.345,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct ObsDiagnostic {
    pub residual_sky: f64,
    pub leverage: f64,
    pub studentized: f64,
}

#[derive(Clone, Debug, Default)]
pub struct RankInfo {
    pub n_free: usize,
    pub rank: usize,
    pub sigma_max: f64,
    pub sigma_min: f64,
    pub tol: f64,
}

impl RankInfo {
    pub fn is_deficient(&self) -> bool {
        self.rank < self.n_free
    }

    pub fn condition_number(&self) -> f64 {
        if self.sigma_min > 0.0 {
            self.sigma_max / self.sigma_min
        } else {
            f64::INFINITY
        }
    }
}

/// All inputs needed to run a fit. Bundles arguments that were previously
/// passed individually so call signatures stay readable.
pub struct FitInputs<'a> {
    pub observations: &'a [&'a Observation],
    pub terms: &'a [Box<dyn Term>],
    pub fixed: &'a [bool],
    pub parallel: &'a [bool],
    pub coefficients: &'a [f64],
    pub latitude: f64,
    pub fit_tol: f64,
}

pub fn fit_model(
    observations: &[&Observation],
    terms: &[Box<dyn Term>],
    fixed: &[bool],
    coefficients: &[f64],
    latitude: f64,
    fit_tol: f64,
) -> Result<FitResult> {
    let parallel = vec![true; terms.len()];
    let inputs = FitInputs {
        observations,
        terms,
        fixed,
        parallel: &parallel,
        coefficients,
        latitude,
        fit_tol,
    };
    fit_model_with_inputs(&inputs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::terms::create_term;
    use crate::test_support::obs_with_ha_offset as obs;

    // --- IterOptions::default ----------------------------------------

    #[test]
    fn iter_options_default_values() {
        let o = IterOptions::default();
        assert_eq!(o.max_iter, 50);
        assert_eq!(o.tol_rad, 1.0e-6);
        assert_eq!(o.initial_lambda, 1.0e-3);
        assert_eq!(o.lambda_up, 10.0);
        assert_eq!(o.lambda_down, 10.0);
        assert!(!o.robust);
        assert_eq!(o.huber_k, 1.345);
    }

    // --- RankInfo ----------------------------------------------------

    #[test]
    fn rank_info_is_deficient_when_rank_below_n_free() {
        let info = RankInfo {
            n_free: 5,
            rank: 3,
            sigma_max: 10.0,
            sigma_min: 1.0,
            tol: 0.0,
        };
        assert!(info.is_deficient());
    }

    #[test]
    fn rank_info_is_not_deficient_when_full_rank() {
        let info = RankInfo {
            n_free: 5,
            rank: 5,
            sigma_max: 10.0,
            sigma_min: 1.0,
            tol: 0.0,
        };
        assert!(!info.is_deficient());
    }

    #[test]
    fn rank_info_condition_number_is_ratio_when_sigma_min_positive() {
        let info = RankInfo {
            n_free: 3,
            rank: 3,
            sigma_max: 100.0,
            sigma_min: 5.0,
            tol: 0.0,
        };
        assert!((info.condition_number() - 20.0).abs() < 1e-12);
    }

    #[test]
    fn rank_info_condition_number_is_infinite_when_sigma_min_zero() {
        let info = RankInfo {
            n_free: 3,
            rank: 2,
            sigma_max: 10.0,
            sigma_min: 0.0,
            tol: 0.0,
        };
        assert!(info.condition_number().is_infinite());
    }

    #[test]
    fn rank_info_default_is_deficient_no_and_inf_condition() {
        // RankInfo::default() → all zeros; sigma_min=0 → infinity.
        let info = RankInfo::default();
        // 0 < 0 is false, so not deficient.
        assert!(!info.is_deficient());
        assert!(info.condition_number().is_infinite());
    }

    // --- fit_model wrapper -------------------------------------------

    #[test]
    fn fit_model_wraps_fit_model_with_inputs_using_all_parallel() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap()];
        let result =
            fit_model(&observations, &terms, &[false], &[0.0], 0.0, 1.0e-9).unwrap();
        assert!((result.coefficients[0] - (-100.0)).abs() < 1e-6);
    }

    #[test]
    fn fit_model_propagates_validation_error() {
        let observations: Vec<&Observation> = Vec::new();
        let terms: Vec<Box<dyn Term>> = Vec::new();
        let result = fit_model(&observations, &terms, &[], &[], 0.0, 1.0e-9);
        assert!(result.is_err());
    }

    // --- ObsDiagnostic / IterReport default ---------------------------

    #[test]
    fn obs_diagnostic_default_is_all_zeros() {
        let d = ObsDiagnostic::default();
        assert_eq!(d.residual_sky, 0.0);
        assert_eq!(d.leverage, 0.0);
        assert_eq!(d.studentized, 0.0);
    }

    #[test]
    fn iter_report_default_is_unconverged_at_zero() {
        let r = IterReport::default();
        assert_eq!(r.iterations, 0);
        assert!(!r.converged);
        assert_eq!(r.final_lambda, 0.0);
        assert_eq!(r.final_step_norm, 0.0);
    }
}
