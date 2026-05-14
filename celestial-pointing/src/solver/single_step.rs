use super::design::{
    build_design_matrix, build_residuals, cos_dec_per_obs, extract_columns,
    scale_a_ha_rows_by_cos_dec, scale_ha_rows_by_cos_dec, subtract_all_contributions,
};
use super::linalg::{solve_and_covariance, SolveOutcome};
use super::stats::{compute_diagnostics, compute_sigma_from_covariance, compute_sky_rms, compute_popn_sd};
use super::validate::{collect_indices, validate_fit_inputs};
use super::{FitInputs, FitResult, IterReport, RankInfo};
use crate::error::Result;
use nalgebra::DVector;

pub(super) fn fit_model_with_inputs(inputs: &FitInputs<'_>) -> Result<FitResult> {
    validate_fit_inputs(inputs)?;
    let free_indices = collect_indices(inputs.fixed, false);
    let fixed_indices = collect_indices(inputs.fixed, true);
    let cos_dec = cos_dec_per_obs(inputs.observations);

    let outcome = single_step_fit(inputs, &cos_dec, &free_indices, &fixed_indices)?;

    let sigma = compute_sigma_from_covariance(
        &outcome.solved,
        &outcome.free_residuals,
        &free_indices,
        inputs.terms.len(),
        inputs.observations.len(),
    );
    let sky_rms = compute_sky_rms(&outcome.actual_residuals_ha, inputs.observations);
    let popn_sd = compute_popn_sd(
        &outcome.actual_residuals_ha,
        inputs.observations,
        free_indices.len(),
    );
    let diagnostics = compute_diagnostics(
        &outcome.free_residuals,
        &outcome.solved.leverage,
        outcome.solved.rank,
    );
    let term_names = inputs.terms.iter().map(|t| t.name().to_string()).collect();
    let rank_info = RankInfo {
        n_free: free_indices.len(),
        rank: outcome.solved.rank,
        sigma_max: outcome.solved.sigma_max,
        sigma_min: outcome.solved.sigma_min,
        tol: inputs.fit_tol,
    };

    Ok(FitResult {
        coefficients: outcome.coefficients,
        change: outcome.change,
        sigma,
        sky_rms,
        popn_sd,
        term_names,
        rank_info,
        diagnostics,
        iter_report: IterReport {
            iterations: 1,
            converged: true,
            final_lambda: 0.0,
            final_step_norm: 0.0,
        },
    })
}

struct FitOutcome {
    coefficients: Vec<f64>,
    change: Vec<f64>,
    solved: SolveOutcome,
    free_residuals: DVector<f64>,
    actual_residuals_ha: DVector<f64>,
}

fn single_step_fit(
    inputs: &FitInputs<'_>,
    cos_dec: &[f64],
    free_indices: &[usize],
    fixed_indices: &[usize],
) -> Result<FitOutcome> {
    let mut b = build_residuals(inputs.observations);
    let a_full = build_design_matrix(inputs.observations, inputs.terms, inputs.latitude);

    scale_ha_rows_by_cos_dec(&mut b, cos_dec);
    let a_full_sky = scale_a_ha_rows_by_cos_dec(&a_full, cos_dec);

    subtract_all_contributions(&mut b, &a_full_sky, inputs.coefficients);

    let a_free = extract_columns(&a_full_sky, free_indices);
    let solved = solve_and_covariance(&a_free, &b, inputs.fit_tol)?;
    let free_residuals = &b - &a_free * &solved.x;

    let mut change = vec![0.0; inputs.terms.len()];
    for (fi, &idx) in free_indices.iter().enumerate() {
        change[idx] = solved.x[fi];
    }
    for &idx in fixed_indices {
        change[idx] = 0.0;
    }

    let new_coeffs: Vec<f64> = inputs
        .coefficients
        .iter()
        .zip(change.iter())
        .map(|(c, d)| c + d)
        .collect();

    let full_residuals_ha = build_residuals(inputs.observations);
    let new_coeffs_dv = DVector::from_vec(new_coeffs.clone());
    let actual_residuals_ha = &full_residuals_ha - &a_full * &new_coeffs_dv;

    Ok(FitOutcome {
        coefficients: new_coeffs,
        change,
        solved,
        free_residuals,
        actual_residuals_ha,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observation::Observation;
    use crate::terms::{create_term, Term};
    use crate::test_support::obs_with_ha_offset as obs;

    fn run_fit<'a>(
        observations: &'a [&'a Observation],
        terms: &'a [Box<dyn Term>],
        fixed: &'a [bool],
        coefficients: &'a [f64],
    ) -> FitResult {
        let parallel = vec![true; terms.len()];
        let inputs = FitInputs {
            observations,
            terms,
            fixed,
            parallel: &parallel,
            coefficients,
            latitude: 0.0,
            fit_tol: 1.0e-9,
        };
        fit_model_with_inputs(&inputs).unwrap()
    }

    // Three observations all showing actual_ha = 100" above commanded → IH = -100".
    #[test]
    fn fit_recovers_ih_offset_from_three_observations() {
        let o1 = obs(100.0, 30.0);
        let o2 = obs(100.0, 45.0);
        let o3 = obs(100.0, 60.0);
        let observations: Vec<&Observation> = vec![&o1, &o2, &o3];
        let terms = vec![create_term("IH").unwrap()];
        let result = run_fit(&observations, &terms, &[false], &[0.0]);
        assert_eq!(result.term_names, vec!["IH"]);
        assert!(
            (result.coefficients[0] - (-100.0)).abs() < 1e-6,
            "expected -100, got {}",
            result.coefficients[0],
        );
    }

    // single-step marks itself converged with iterations=1.
    #[test]
    fn fit_iter_report_records_single_iteration() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap()];
        let result = run_fit(&observations, &terms, &[false], &[0.0]);
        assert_eq!(result.iter_report.iterations, 1);
        assert!(result.iter_report.converged);
        assert_eq!(result.iter_report.final_lambda, 0.0);
    }

    #[test]
    fn fit_starting_from_non_zero_coefficient_still_converges() {
        let o1 = obs(100.0, 30.0);
        let o2 = obs(100.0, 45.0);
        let observations: Vec<&Observation> = vec![&o1, &o2];
        let terms = vec![create_term("IH").unwrap()];
        // Starting from -50, change should be -50 to land at -100.
        let result = run_fit(&observations, &terms, &[false], &[-50.0]);
        assert!((result.coefficients[0] - (-100.0)).abs() < 1e-6);
        assert!((result.change[0] - (-50.0)).abs() < 1e-6);
    }

    // Fixed term contributes via subtract_all_contributions but is not solved
    // for. Its coefficient should not change.
    #[test]
    fn fit_fixed_term_keeps_its_coefficient_and_zero_change() {
        let o1 = obs(100.0, 30.0);
        let o2 = obs(100.0, 45.0);
        let observations: Vec<&Observation> = vec![&o1, &o2];
        let terms = vec![
            create_term("IH").unwrap(),
            create_term("ID").unwrap(),
        ];
        // Fix ID at 50, leave IH free.
        let result = run_fit(&observations, &terms, &[false, true], &[0.0, 50.0]);
        assert_eq!(result.coefficients[1], 50.0);
        assert_eq!(result.change[1], 0.0);
    }

    #[test]
    fn fit_change_vector_equals_new_coefficients_minus_initial() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap()];
        let result = run_fit(&observations, &terms, &[false], &[0.0]);
        assert!((result.coefficients[0] - result.change[0]).abs() < 1e-12);
    }

    #[test]
    fn fit_sigma_vector_length_matches_terms() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap(), create_term("ID").unwrap()];
        let result = run_fit(&observations, &terms, &[false, false], &[0.0, 0.0]);
        assert_eq!(result.sigma.len(), 2);
    }

    #[test]
    fn fit_diagnostics_length_matches_observations() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap()];
        let result = run_fit(&observations, &terms, &[false], &[0.0]);
        assert_eq!(result.diagnostics.len(), 3);
    }

    #[test]
    fn fit_invalid_inputs_propagate_validation_error() {
        // Empty terms → "no terms to fit".
        let observations: Vec<&Observation> = Vec::new();
        let terms: Vec<Box<dyn Term>> = Vec::new();
        let parallel: Vec<bool> = Vec::new();
        let inputs = FitInputs {
            observations: &observations,
            terms: &terms,
            fixed: &[],
            parallel: &parallel,
            coefficients: &[],
            latitude: 0.0,
            fit_tol: 1.0e-9,
        };
        assert!(fit_model_with_inputs(&inputs).is_err());
    }

    // With perfect data (zero residual after fit), sky_rms should be tiny.
    #[test]
    fn fit_perfect_data_has_near_zero_sky_rms() {
        let o1 = obs(100.0, 30.0);
        let o2 = obs(100.0, 45.0);
        let observations: Vec<&Observation> = vec![&o1, &o2];
        let terms = vec![create_term("IH").unwrap()];
        let result = run_fit(&observations, &terms, &[false], &[0.0]);
        assert!(result.sky_rms.abs() < 1e-6, "got sky_rms = {}", result.sky_rms);
    }
}
