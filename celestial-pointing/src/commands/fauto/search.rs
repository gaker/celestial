use super::bic::{compute_bic, fit_is_well_conditioned, try_fit};
use super::candidates::{build_candidate_pool, filter_remaining, physical_candidates};
use super::moves::{commit_move, score_best_move};
use super::prune::{backward_eliminate, prune_terms, MIN_SIGNIFICANCE};
use super::BASE_TERMS;
use crate::error::{Error, Result};
use crate::observation::Observation;
use crate::solver::FitResult;

pub(super) const OUTLIER_SIGMA: f64 = 2.5;
pub(super) const MAX_OUTLIER_PASSES: usize = 5;

pub(super) struct TraceEntry {
    pub(super) label: String,
    pub(super) delta_bic: f64,
    pub(super) rms: f64,
}

pub(super) fn append_trace(report: &mut String, trace: &[TraceEntry]) {
    for e in trace {
        report.push_str(&format!(
            "+ {} (dBIC={:.1}, RMS={:.2}\")\n",
            e.label, e.delta_bic, e.rms,
        ));
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn run_search_with_outlier_passes(
    prepared: &[Observation],
    masked: &mut std::collections::HashSet<usize>,
    active: &mut Vec<String>,
    parallel_groups: &mut Vec<(String, String)>,
    report: &mut String,
    max_terms: usize,
    bic_threshold: f64,
    latitude: f64,
    fit_tol: f64,
) -> Result<(FitResult, Vec<usize>)> {
    let mut newly_masked = Vec::new();
    let mut last_fit = None;
    let pre_existing: std::collections::HashSet<usize> = masked.iter().copied().collect();
    for pass in 0..MAX_OUTLIER_PASSES {
        let observations = active_observations(prepared, masked);
        if observations.len() < BASE_TERMS.len() {
            return Err(Error::Fit("masking removed too many observations".into()));
        }
        let fit = run_single_pass(
            &observations,
            active,
            parallel_groups,
            report,
            max_terms,
            bic_threshold,
            latitude,
            fit_tol,
            pass,
        )?;
        let outliers = find_outliers(prepared, masked, &fit, OUTLIER_SIGMA);
        if outliers.is_empty() {
            last_fit = Some(fit);
            break;
        }
        report_and_mask_outliers(
            prepared,
            masked,
            &pre_existing,
            &mut newly_masked,
            &outliers,
            &fit,
            report,
        );
        last_fit = Some(fit);
    }
    Ok((
        last_fit.ok_or_else(|| Error::Fit("FAUTO produced no fit".into()))?,
        newly_masked,
    ))
}

fn report_and_mask_outliers(
    prepared: &[Observation],
    masked: &mut std::collections::HashSet<usize>,
    pre_existing: &std::collections::HashSet<usize>,
    newly_masked: &mut Vec<usize>,
    outliers: &[usize],
    fit: &FitResult,
    report: &mut String,
) {
    let scale = if fit.popn_sd > 0.0 {
        fit.popn_sd
    } else {
        fit.sky_rms
    };
    for &idx in outliers {
        let sigma = fit
            .diagnostics
            .get(map_to_active_index(prepared, masked, idx))
            .map(|d| d.residual_sky / scale)
            .unwrap_or(OUTLIER_SIGMA);
        masked.insert(idx);
        if !pre_existing.contains(&idx) {
            newly_masked.push(idx);
        }
        report.push_str(&format!(
            "x obs #{} out by {:.1}\u{03c3}, MASKed\n",
            idx + 1,
            sigma,
        ));
    }
}

pub(super) fn active_observations<'a>(
    prepared: &'a [Observation],
    masked: &std::collections::HashSet<usize>,
) -> Vec<&'a Observation> {
    prepared
        .iter()
        .enumerate()
        .filter(|(i, _)| !masked.contains(i))
        .map(|(_, o)| o)
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn run_single_pass(
    observations: &[&Observation],
    active: &mut Vec<String>,
    parallel_groups: &mut Vec<(String, String)>,
    report: &mut String,
    max_terms: usize,
    bic_threshold: f64,
    latitude: f64,
    fit_tol: f64,
    pass: usize,
) -> Result<FitResult> {
    let base_fit = try_fit(observations, active, latitude, fit_tol)?;
    let n_obs = observations.len();
    let mut current_bic = compute_bic(n_obs, active.len(), base_fit.sky_rms);
    write_pass_header(report, pass, active, &base_fit, n_obs, current_bic);
    let physical_trace = physical_prepass(
        observations,
        active,
        &mut current_bic,
        bic_threshold,
        latitude,
        fit_tol,
    )?;
    let trace = search_loop(
        observations,
        active,
        parallel_groups,
        &mut current_bic,
        bic_threshold,
        max_terms,
        latitude,
        fit_tol,
    )?;
    append_trace(report, &physical_trace);
    append_trace(report, &trace);
    let backed_out =
        backward_eliminate(observations, active, parallel_groups, latitude, fit_tol)?;
    for name in &backed_out {
        report.push_str(&format!("- {} (backward-eliminated)\n", name));
    }
    let pruned = prune_terms(observations, active, latitude, fit_tol)?;
    for name in &pruned {
        report.push_str(&format!(
            "- {} (pruned, significance < {:.1})\n",
            name, MIN_SIGNIFICANCE
        ));
    }
    try_fit(observations, active, latitude, fit_tol)
}

fn write_pass_header(
    report: &mut String,
    pass: usize,
    active: &[String],
    base_fit: &FitResult,
    n_obs: usize,
    current_bic: f64,
) {
    if pass == 0 {
        report.push_str(&format!(
            "Base: {} (BIC={:.1}, RMS={:.2}\")\n",
            active.join(" "),
            current_bic,
            base_fit.sky_rms,
        ));
    } else {
        report.push_str(&format!(
            "\nPass {} (N={}, RMS={:.2}\")\n",
            pass + 1,
            n_obs,
            base_fit.sky_rms,
        ));
    }
}

fn physical_prepass(
    observations: &[&Observation],
    active: &mut Vec<String>,
    current_bic: &mut f64,
    threshold: f64,
    latitude: f64,
    fit_tol: f64,
) -> Result<Vec<TraceEntry>> {
    let mut trace = Vec::new();
    let n_obs = observations.len();
    for &candidate in physical_candidates() {
        let mut trial = active.clone();
        trial.push(candidate.to_string());
        let fit = match try_fit(observations, &trial, latitude, fit_tol) {
            Ok(f) => f,
            Err(_) => continue,
        };
        if !fit_is_well_conditioned(&fit) {
            continue;
        }
        let new_bic = compute_bic(n_obs, trial.len(), fit.sky_rms);
        let delta = new_bic - *current_bic;
        if delta < threshold {
            trace.push(TraceEntry {
                label: candidate.to_string(),
                delta_bic: delta,
                rms: fit.sky_rms,
            });
            active.push(candidate.to_string());
            *current_bic = new_bic;
        }
    }
    Ok(trace)
}

#[allow(clippy::too_many_arguments)]
fn search_loop(
    observations: &[&Observation],
    active: &mut Vec<String>,
    parallel_groups: &mut Vec<(String, String)>,
    current_bic: &mut f64,
    threshold: f64,
    max_terms: usize,
    latitude: f64,
    fit_tol: f64,
) -> Result<Vec<TraceEntry>> {
    let mut trace = Vec::new();
    let candidates = build_candidate_pool();
    loop {
        if active.len() >= max_terms {
            break;
        }
        let remaining = filter_remaining(&candidates, active);
        if remaining.singles.is_empty()
            && remaining.pairs.is_empty()
            && remaining.triples.is_empty()
        {
            break;
        }
        let best = score_best_move(observations, active, &remaining, latitude, fit_tol)?;
        let Some(move_) = best else { break };
        let new_bic = compute_bic(
            observations.len(),
            active.len() + move_.added_count(),
            move_.rms,
        );
        let delta = new_bic - *current_bic;
        if delta >= threshold {
            break;
        }
        commit_move(active, parallel_groups, &move_);
        *current_bic = new_bic;
        trace.push(TraceEntry {
            label: move_.label(),
            delta_bic: delta,
            rms: move_.rms,
        });
    }
    Ok(trace)
}

pub(super) fn find_outliers(
    prepared: &[Observation],
    masked: &std::collections::HashSet<usize>,
    fit: &FitResult,
    threshold: f64,
) -> Vec<usize> {
    let scale = if fit.popn_sd > 0.0 {
        fit.popn_sd
    } else {
        fit.sky_rms
    };
    if scale <= 0.0 {
        return Vec::new();
    }
    let active_indices: Vec<usize> = prepared
        .iter()
        .enumerate()
        .filter(|(i, _)| !masked.contains(i))
        .map(|(i, _)| i)
        .collect();
    let mut out = Vec::new();
    for (active_pos, &abs_idx) in active_indices.iter().enumerate() {
        let Some(diag) = fit.diagnostics.get(active_pos) else {
            continue;
        };
        let sigma_dev = diag.residual_sky / scale;
        if sigma_dev > threshold {
            out.push(abs_idx);
        }
    }
    out
}

pub(super) fn map_to_active_index(
    prepared: &[Observation],
    masked: &std::collections::HashSet<usize>,
    abs_idx: usize,
) -> usize {
    prepared
        .iter()
        .enumerate()
        .filter(|(i, _)| !masked.contains(i))
        .position(|(i, _)| i == abs_idx)
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::ObsDiagnostic;
    use crate::test_support::{diag, obs, FitResultBuilder};

    fn fit(
        sky_rms: f64,
        popn_sd: f64,
        diagnostics: Vec<ObsDiagnostic>,
    ) -> FitResult {
        FitResultBuilder::new()
            .sky_rms(sky_rms)
            .popn_sd(popn_sd)
            .diagnostics(diagnostics)
            .build()
    }

    fn masked_set(indices: &[usize]) -> std::collections::HashSet<usize> {
        indices.iter().copied().collect()
    }

    // active_observations should yield references in original order, skipping masked.
    #[test]
    fn active_observations_skips_masked_indices() {
        let prepared = vec![obs(), obs(), obs(), obs()];
        let masked = masked_set(&[1, 3]);
        let active = active_observations(&prepared, &masked);
        assert_eq!(active.len(), 2);
        // Pointer identity confirms we got refs to prepared[0] and prepared[2].
        assert!(std::ptr::eq(active[0], &prepared[0]));
        assert!(std::ptr::eq(active[1], &prepared[2]));
    }

    #[test]
    fn active_observations_empty_mask_returns_everything() {
        let prepared = vec![obs(), obs()];
        let masked = std::collections::HashSet::new();
        assert_eq!(active_observations(&prepared, &masked).len(), 2);
    }

    #[test]
    fn active_observations_all_masked_returns_empty() {
        let prepared = vec![obs(), obs()];
        let masked = masked_set(&[0, 1]);
        assert!(active_observations(&prepared, &masked).is_empty());
    }

    // map_to_active_index converts an absolute prepared-index into its
    // position within the unmasked subsequence.
    #[test]
    fn map_to_active_index_with_no_masking_is_identity() {
        let prepared = vec![obs(); 5];
        let masked = std::collections::HashSet::new();
        assert_eq!(map_to_active_index(&prepared, &masked, 0), 0);
        assert_eq!(map_to_active_index(&prepared, &masked, 3), 3);
    }

    #[test]
    fn map_to_active_index_skips_over_masked_entries() {
        let prepared = vec![obs(); 5];
        let masked = masked_set(&[1, 3]);
        // Active sequence is [0, 2, 4].
        assert_eq!(map_to_active_index(&prepared, &masked, 0), 0);
        assert_eq!(map_to_active_index(&prepared, &masked, 2), 1);
        assert_eq!(map_to_active_index(&prepared, &masked, 4), 2);
    }

    // The function defaults to 0 if the index isn't in the active set.
    // Documenting this so anyone hardening the contract sees the test fail.
    #[test]
    fn map_to_active_index_masked_index_falls_back_to_zero() {
        let prepared = vec![obs(); 3];
        let masked = masked_set(&[1]);
        assert_eq!(map_to_active_index(&prepared, &masked, 1), 0);
    }

    // find_outliers: scale = popn_sd when > 0, else sky_rms. Threshold applied.
    #[test]
    fn find_outliers_returns_empty_when_scale_is_zero() {
        let prepared = vec![obs(), obs()];
        let masked = std::collections::HashSet::new();
        let f = fit(0.0, 0.0, vec![diag(100.0), diag(100.0)]);
        assert!(find_outliers(&prepared, &masked, &f, 2.5).is_empty());
    }

    #[test]
    fn find_outliers_uses_popn_sd_when_available() {
        // popn_sd = 1.0, threshold = 2.5. Residuals: 1.0 (not outlier), 3.0 (outlier).
        // If sky_rms were used as scale, 3.0/5.0 = 0.6 would not trigger.
        let prepared = vec![obs(), obs()];
        let masked = std::collections::HashSet::new();
        let f = fit(5.0, 1.0, vec![diag(1.0), diag(3.0)]);
        let out = find_outliers(&prepared, &masked, &f, 2.5);
        assert_eq!(out, vec![1]);
    }

    #[test]
    fn find_outliers_falls_back_to_sky_rms_when_popn_sd_zero() {
        // popn_sd = 0 forces sky_rms = 2.0 scale. Threshold 2.5 → cutoff 5.0.
        let prepared = vec![obs(), obs()];
        let masked = std::collections::HashSet::new();
        let f = fit(2.0, 0.0, vec![diag(1.0), diag(10.0)]);
        let out = find_outliers(&prepared, &masked, &f, 2.5);
        assert_eq!(out, vec![1]);
    }

    #[test]
    fn find_outliers_returns_absolute_indices_not_active_positions() {
        // prepared[0] masked → active sequence is [prepared[1], prepared[2]].
        // diag[0] corresponds to prepared[1], diag[1] to prepared[2].
        // diag[1] is the outlier → reported as abs_idx 2.
        let prepared = vec![obs(), obs(), obs()];
        let masked = masked_set(&[0]);
        let f = fit(0.0, 1.0, vec![diag(0.5), diag(10.0)]);
        let out = find_outliers(&prepared, &masked, &f, 2.5);
        assert_eq!(out, vec![2]);
    }

    #[test]
    fn find_outliers_missing_diagnostic_skips_silently() {
        let prepared = vec![obs(), obs()];
        let masked = std::collections::HashSet::new();
        let f = fit(0.0, 1.0, vec![diag(10.0)]); // diagnostics shorter than active
        let out = find_outliers(&prepared, &masked, &f, 2.5);
        assert_eq!(out, vec![0]);
    }

    #[test]
    fn append_trace_formats_each_entry() {
        let mut report = String::new();
        let trace = vec![
            TraceEntry {
                label: "TF".into(),
                delta_bic: -10.0,
                rms: 1.23,
            },
            TraceEntry {
                label: "HDSH&HDCH".into(),
                delta_bic: -5.5,
                rms: 0.99,
            },
        ];
        append_trace(&mut report, &trace);
        assert!(report.contains("+ TF (dBIC=-10.0, RMS=1.23\")"));
        assert!(report.contains("+ HDSH&HDCH (dBIC=-5.5, RMS=0.99\")"));
    }

    #[test]
    fn append_trace_empty_is_a_noop() {
        let mut report = String::from("preserved");
        append_trace(&mut report, &[]);
        assert_eq!(report, "preserved");
    }

    // Constants we depend on — pinning so any change is intentional.
    #[test]
    fn outlier_sigma_is_two_point_five() {
        assert_eq!(OUTLIER_SIGMA, 2.5);
    }

    #[test]
    fn max_outlier_passes_is_five() {
        assert_eq!(MAX_OUTLIER_PASSES, 5);
    }
}
