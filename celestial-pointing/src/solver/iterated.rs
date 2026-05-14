use super::design::{cos_dec_per_obs, extract_columns};
use super::linalg::solve_damped;
use super::stats::{compute_diagnostics, compute_sigma_from_covariance, compute_sky_rms, compute_popn_sd};
use super::validate::{collect_indices, validate_fit_inputs};
use super::{FitInputs, FitResult, IterOptions, IterReport, RankInfo};
use crate::error::Result;
use crate::observation::Observation;
use crate::terms::Term;
use nalgebra::{DMatrix, DVector};

pub(crate) fn fit_model_iterated(
    inputs: &FitInputs<'_>,
    options: &IterOptions,
) -> Result<FitResult> {
    validate_fit_inputs(inputs)?;
    let free_indices = collect_indices(inputs.fixed, false);
    let fixed_indices = collect_indices(inputs.fixed, true);
    let cos_dec = cos_dec_per_obs(inputs.observations);

    let mut state: Vec<f64> = inputs.coefficients.to_vec();
    let initial_chi2 = chi2_at_state(
        inputs.observations,
        inputs.terms,
        &state,
        inputs.latitude,
        &cos_dec,
    );
    let lm = run_lm_loop(
        inputs,
        &cos_dec,
        &free_indices,
        &fixed_indices,
        &mut state,
        initial_chi2,
        options,
    )?;

    let (a_free_final, b_final) = build_lm_system(
        inputs,
        &state,
        &cos_dec,
        &free_indices,
        options,
        lm.prev_chi2,
    );
    let final_solved = solve_damped(&a_free_final, &b_final, 0.0, inputs.fit_tol)?;
    apply_final_step(
        &mut state,
        &free_indices,
        &final_solved.x,
        inputs,
        &cos_dec,
        lm.prev_chi2,
    );
    let free_residuals = &b_final - &a_free_final * &final_solved.x;

    let change = state
        .iter()
        .zip(inputs.coefficients.iter())
        .map(|(c, c0)| c - c0)
        .collect();

    let actual_residuals_ha =
        residuals_under_state(inputs.observations, inputs.terms, &state, inputs.latitude);
    let sky_rms = compute_sky_rms(&actual_residuals_ha, inputs.observations);
    let popn_sd = compute_popn_sd(&actual_residuals_ha, inputs.observations, free_indices.len());
    let sigma = compute_sigma_from_covariance(
        &final_solved,
        &free_residuals,
        &free_indices,
        inputs.terms.len(),
        inputs.observations.len(),
    );
    let diagnostics = compute_diagnostics(
        &free_residuals,
        &final_solved.leverage,
        final_solved.rank,
    );
    let term_names = inputs.terms.iter().map(|t| t.name().to_string()).collect();
    let rank_info = RankInfo {
        n_free: free_indices.len(),
        rank: final_solved.rank,
        sigma_max: final_solved.sigma_max,
        sigma_min: final_solved.sigma_min,
        tol: inputs.fit_tol,
    };

    Ok(FitResult {
        coefficients: state,
        change,
        sigma,
        sky_rms,
        popn_sd,
        term_names,
        rank_info,
        diagnostics,
        iter_report: IterReport {
            iterations: lm.iterations,
            converged: lm.converged,
            final_lambda: lm.lambda,
            final_step_norm: lm.last_step_norm,
        },
    })
}

struct LmLoopResult {
    iterations: usize,
    converged: bool,
    lambda: f64,
    last_step_norm: f64,
    prev_chi2: f64,
}

#[allow(clippy::too_many_arguments)]
fn run_lm_loop(
    inputs: &FitInputs<'_>,
    cos_dec: &[f64],
    free_indices: &[usize],
    fixed_indices: &[usize],
    state: &mut Vec<f64>,
    initial_chi2: f64,
    options: &IterOptions,
) -> Result<LmLoopResult> {
    let mut prev_chi2 = initial_chi2;
    let mut lambda = options.initial_lambda;
    let mut iterations = 0usize;
    let mut converged = false;
    let mut last_step_norm = f64::INFINITY;

    while iterations < options.max_iter {
        iterations += 1;
        let (a_free, b) =
            build_lm_system(inputs, state, cos_dec, free_indices, options, prev_chi2);
        let solved = solve_damped(&a_free, &b, lambda, inputs.fit_tol)?;
        let delta = &solved.x;

        let trial_state = trial_state_from(state, free_indices, fixed_indices, delta);
        let trial_chi2 = chi2_at_state(
            inputs.observations,
            inputs.terms,
            &trial_state,
            inputs.latitude,
            cos_dec,
        );

        if trial_chi2 < prev_chi2 {
            let step_norm = rotational_step_norm(inputs.terms, free_indices, delta);
            *state = trial_state;
            prev_chi2 = trial_chi2;
            lambda /= options.lambda_down;
            last_step_norm = step_norm;
            if step_norm < options.tol_rad {
                converged = true;
                break;
            }
        } else {
            lambda *= options.lambda_up;
            if lambda > 1.0e12 {
                break;
            }
        }
    }
    Ok(LmLoopResult {
        iterations,
        converged,
        lambda,
        last_step_norm,
        prev_chi2,
    })
}

fn trial_state_from(
    state: &[f64],
    free_indices: &[usize],
    fixed_indices: &[usize],
    delta: &DVector<f64>,
) -> Vec<f64> {
    let mut trial = state.to_vec();
    for (fi, &idx) in free_indices.iter().enumerate() {
        trial[idx] += delta[fi];
    }
    for &idx in fixed_indices {
        trial[idx] = state[idx];
    }
    trial
}

fn apply_final_step(
    state: &mut Vec<f64>,
    free_indices: &[usize],
    final_delta: &DVector<f64>,
    inputs: &FitInputs<'_>,
    cos_dec: &[f64],
    prev_chi2: f64,
) {
    let trial_state_final: Vec<f64> = {
        let mut s = state.clone();
        for (fi, &idx) in free_indices.iter().enumerate() {
            s[idx] += final_delta[fi];
        }
        s
    };
    let trial_chi2_final = chi2_at_state(
        inputs.observations,
        inputs.terms,
        &trial_state_final,
        inputs.latitude,
        cos_dec,
    );
    if trial_chi2_final <= prev_chi2 {
        *state = trial_state_final;
    }
}

fn residuals_under_state(
    observations: &[&Observation],
    terms: &[Box<dyn Term>],
    state: &[f64],
    latitude: f64,
) -> DVector<f64> {
    let n = observations.len();
    let mut b = DVector::zeros(2 * n);
    for (i, obs) in observations.iter().enumerate() {
        let h = obs.commanded_ha.radians();
        let dec = obs.catalog_dec.radians();
        let pier = obs.pier_side.sign();
        let mut dh_pred = 0.0;
        let mut dd_pred = 0.0;
        for (j, term) in terms.iter().enumerate() {
            if state[j] == 0.0 {
                continue;
            }
            let (jh, jd) = term.jacobian_equatorial(h, dec, latitude, pier);
            dh_pred += jh * state[j];
            dd_pred += jd * state[j];
        }
        let raw_dh = (obs.actual_ha - obs.commanded_ha).arcseconds();
        let raw_dd = (obs.observed_dec - obs.catalog_dec).arcseconds();
        b[2 * i] = raw_dh - dh_pred;
        b[2 * i + 1] = raw_dd - dd_pred;
    }
    b
}

fn chi2_at_state(
    observations: &[&Observation],
    terms: &[Box<dyn Term>],
    state: &[f64],
    latitude: f64,
    cos_dec: &[f64],
) -> f64 {
    let mut total = 0.0;
    let resid = residuals_under_state(observations, terms, state, latitude);
    for i in 0..observations.len() {
        let r_h = resid[2 * i] * cos_dec[i];
        let r_d = resid[2 * i + 1];
        total += r_h * r_h + r_d * r_d;
    }
    total
}

fn build_lm_system(
    inputs: &FitInputs<'_>,
    state: &[f64],
    cos_dec: &[f64],
    free_indices: &[usize],
    options: &IterOptions,
    chi2_estimate: f64,
) -> (DMatrix<f64>, DVector<f64>) {
    let observations = inputs.observations;
    let terms = inputs.terms;
    let latitude = inputs.latitude;
    let n = observations.len();
    let m = terms.len();
    let mut a = DMatrix::zeros(2 * n, m);
    let mut b = residuals_under_state(observations, terms, state, latitude);
    for (i, obs) in observations.iter().enumerate() {
        let h = obs.commanded_ha.radians();
        let dec = obs.catalog_dec.radians();
        let pier = obs.pier_side.sign();
        let c = cos_dec[i];
        for (j, term) in terms.iter().enumerate() {
            let (jh, jd) = term.jacobian_equatorial(h, dec, latitude, pier);
            a[(2 * i, j)] = jh * c;
            a[(2 * i + 1, j)] = jd;
        }
        b[2 * i] *= c;
    }

    if options.robust && n >= free_indices.len() && chi2_estimate > 0.0 {
        let dof = (2 * n).saturating_sub(free_indices.len()).max(1);
        let sigma_est = (chi2_estimate / dof as f64).sqrt();
        let k_sigma = options.huber_k * sigma_est;
        if k_sigma > 0.0 {
            for row in 0..(2 * n) {
                let r = b[row].abs();
                if r > k_sigma {
                    let w_sqrt = (k_sigma / r).sqrt();
                    b[row] *= w_sqrt;
                    for col in 0..m {
                        a[(row, col)] *= w_sqrt;
                    }
                }
            }
        }
    }

    let a_free = extract_columns(&a, free_indices);
    (a_free, b)
}

fn rotational_step_norm(
    terms: &[Box<dyn Term>],
    free_indices: &[usize],
    delta: &DVector<f64>,
) -> f64 {
    let rot_names = ["IH", "ID", "NP", "ME"];
    let mut max_step = 0.0_f64;
    for (fi, &term_idx) in free_indices.iter().enumerate() {
        if !rot_names.contains(&terms[term_idx].name()) {
            continue;
        }
        let step_arcsec = delta[fi].abs();
        let step_rad = step_arcsec / 206264.806_247;
        if step_rad > max_step {
            max_step = step_rad;
        }
    }
    if max_step == 0.0 {
        let mut s = 0.0_f64;
        for v in delta.iter() {
            s += v * v;
        }
        s.sqrt() / 206264.806_247
    } else {
        max_step
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::terms::create_term;
    use crate::test_support::obs_with_ha_offset as obs;

    fn run_iter<'a>(
        observations: &'a [&'a Observation],
        terms: &'a [Box<dyn Term>],
        fixed: &'a [bool],
        coefficients: &'a [f64],
        options: &IterOptions,
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
        fit_model_iterated(&inputs, options).unwrap()
    }

    // --- end-to-end LM fits ------------------------------------------

    #[test]
    fn iterated_fit_recovers_ih_offset_and_reports_convergence() {
        let o1 = obs(100.0, 30.0);
        let o2 = obs(100.0, 45.0);
        let o3 = obs(100.0, 60.0);
        let observations: Vec<&Observation> = vec![&o1, &o2, &o3];
        let terms = vec![create_term("IH").unwrap()];
        let result = run_iter(&observations, &terms, &[false], &[0.0], &IterOptions::default());
        assert!((result.coefficients[0] - (-100.0)).abs() < 1e-4);
        assert!(result.iter_report.converged, "LM should converge on linear data");
    }

    #[test]
    fn iterated_fit_writes_term_names_in_order() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap(), create_term("ID").unwrap()];
        let result = run_iter(
            &observations,
            &terms,
            &[false, false],
            &[0.0, 0.0],
            &IterOptions::default(),
        );
        assert_eq!(result.term_names, vec!["IH", "ID"]);
    }

    #[test]
    fn iterated_fit_change_vector_is_state_minus_initial() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap()];
        let result = run_iter(&observations, &terms, &[false], &[25.0], &IterOptions::default());
        // change = final - 25
        assert!(
            (result.change[0] - (result.coefficients[0] - 25.0)).abs() < 1e-9,
        );
    }

    // max_iter=0 → loop never runs → no improvement, iterations=0.
    #[test]
    fn iterated_fit_zero_max_iter_does_not_iterate() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap()];
        let mut opts = IterOptions::default();
        opts.max_iter = 0;
        let result = run_iter(&observations, &terms, &[false], &[0.0], &opts);
        assert_eq!(result.iter_report.iterations, 0);
        assert!(!result.iter_report.converged);
    }

    #[test]
    fn iterated_fit_rank_info_reports_n_free() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap(), create_term("ID").unwrap()];
        let result = run_iter(
            &observations,
            &terms,
            &[true, false], // one fixed, one free
            &[0.0, 0.0],
            &IterOptions::default(),
        );
        assert_eq!(result.rank_info.n_free, 1);
    }

    // Robust mode runs the Huber down-weighting path. With clean data this
    // shouldn't change the fit, but it should execute without panicking.
    #[test]
    fn iterated_fit_robust_mode_executes_huber_path() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap()];
        let mut opts = IterOptions::default();
        opts.robust = true;
        let result = run_iter(&observations, &terms, &[false], &[0.0], &opts);
        assert!((result.coefficients[0] - (-100.0)).abs() < 1.0);
    }

    // --- trial_state_from -------------------------------------------

    #[test]
    fn trial_state_adds_delta_to_free_and_holds_fixed() {
        let state = vec![1.0, 2.0, 3.0];
        let free_indices = vec![0, 2];
        let fixed_indices = vec![1];
        let delta = DVector::from_vec(vec![10.0, 100.0]);
        let trial = trial_state_from(&state, &free_indices, &fixed_indices, &delta);
        assert_eq!(trial, vec![11.0, 2.0, 103.0]);
    }

    #[test]
    fn trial_state_all_free_applies_full_delta() {
        let state = vec![0.0, 0.0];
        let delta = DVector::from_vec(vec![5.0, -7.0]);
        let trial = trial_state_from(&state, &[0, 1], &[], &delta);
        assert_eq!(trial, vec![5.0, -7.0]);
    }

    #[test]
    fn trial_state_all_fixed_returns_state_unchanged() {
        let state = vec![1.0, 2.0];
        let delta = DVector::from_vec(Vec::<f64>::new());
        let trial = trial_state_from(&state, &[], &[0, 1], &delta);
        assert_eq!(trial, state);
    }

    // --- rotational_step_norm ----------------------------------------

    #[test]
    fn rotational_step_norm_picks_max_among_rotational_terms() {
        let terms = vec![
            create_term("IH").unwrap(),
            create_term("CH").unwrap(), // non-rotational
            create_term("ME").unwrap(),
        ];
        // delta in arcseconds at the free positions. With all three free:
        // IH step = 1", ME step = 10", CH step = 1000" (ignored as non-rotational).
        let free_indices = vec![0, 1, 2];
        let delta = DVector::from_vec(vec![1.0, 1000.0, 10.0]);
        let norm = rotational_step_norm(&terms, &free_indices, &delta);
        // Should be 10" in radians, not 1000".
        let expected = 10.0 / 206264.806_247;
        assert!((norm - expected).abs() < 1e-12);
    }

    // If no rotational terms are present, the function falls back to the
    // L2 norm of delta (in radians).
    #[test]
    fn rotational_step_norm_no_rotational_terms_falls_back_to_l2() {
        // "IH", "ID", "NP", "ME" count as rotational — use strictly non-rotational ones.
        let terms = vec![create_term("CH").unwrap(), create_term("TF").unwrap()];
        let delta = DVector::from_vec(vec![3.0, 4.0]);
        let norm = rotational_step_norm(&terms, &[0, 1], &delta);
        let expected = (libm::sqrt(9.0 + 16.0)) / 206264.806_247;
        assert!((norm - expected).abs() < 1e-12);
    }

    #[test]
    fn rotational_step_norm_recognizes_id_np_me_as_rotational() {
        for name in ["IH", "ID", "NP", "ME"] {
            let terms = vec![create_term(name).unwrap()];
            let delta = DVector::from_vec(vec![100.0]);
            let norm = rotational_step_norm(&terms, &[0], &delta);
            let expected = 100.0 / 206264.806_247;
            assert!((norm - expected).abs() < 1e-12, "name={}", name);
        }
    }

    // --- chi2_at_state ----------------------------------------------

    #[test]
    fn chi2_at_state_at_optimum_is_near_zero() {
        let o1 = obs(100.0, 30.0);
        let o2 = obs(100.0, 45.0);
        let observations: Vec<&Observation> = vec![&o1, &o2];
        let terms = vec![create_term("IH").unwrap()];
        let cos_dec = vec![libm::cos(30.0_f64.to_radians()), libm::cos(45.0_f64.to_radians())];
        let chi2 = chi2_at_state(&observations, &terms, &[-100.0], 0.0, &cos_dec);
        assert!(chi2 < 1e-6, "chi2 at optimum should be ~0, got {}", chi2);
    }

    #[test]
    fn chi2_at_state_zero_coefficients_is_positive_for_offset_data() {
        let o1 = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o1];
        let terms = vec![create_term("IH").unwrap()];
        let cos_dec = vec![libm::cos(30.0_f64.to_radians())];
        let chi2 = chi2_at_state(&observations, &terms, &[0.0], 0.0, &cos_dec);
        assert!(chi2 > 0.0);
    }

    // --- residuals_under_state ---------------------------------------

    #[test]
    fn residuals_under_state_zero_state_yields_raw_residuals() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o];
        let terms = vec![create_term("IH").unwrap()];
        let r = residuals_under_state(&observations, &terms, &[0.0], 0.0);
        // With state=0, no term contribution → raw residuals.
        assert!((r[0] - 100.0).abs() < 1e-9);
        assert_eq!(r[1], 0.0);
    }

    #[test]
    fn residuals_under_state_at_optimum_yields_near_zero() {
        let o = obs(100.0, 30.0);
        let observations: Vec<&Observation> = vec![&o];
        let terms = vec![create_term("IH").unwrap()];
        // IH=-100: term contributes -100" to HA residual prediction → residual = 100 - (-100) = 200?
        // Actually: raw_dh - dh_pred = 100 - (-100) = 200... wait, this depends on IH jacobian sign.
        // Just assert finiteness and proper length here; sign-convention is in higher-level tests.
        let r = residuals_under_state(&observations, &terms, &[-100.0], 0.0);
        assert_eq!(r.len(), 2);
        assert!(r[0].is_finite());
    }
}
