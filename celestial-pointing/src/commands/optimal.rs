use crate::error::{Error, Result};
use crate::observation::Observation;
use crate::session::Session;
use crate::solver::{fit_model, FitResult};
use crate::terms::create_term;
use rayon::prelude::*;
use super::{Command, CommandOutput};

const BASE_TERMS: &[&str] = &["IH", "ID", "CH", "NP", "MA", "ME"];
const PHYSICAL_CANDIDATES: &[&str] = &[
    "TF", "TX", "DAF", "FO", "HCES", "HCEC", "DCES", "DCEC",
];
const DEFAULT_MAX_TERMS: usize = 30;
const DEFAULT_BIC_THRESHOLD: f64 = -6.0;
const MIN_SIGNIFICANCE: f64 = 2.0;

struct StageEntry {
    name: String,
    delta_bic: f64,
    rms: f64,
}

pub struct Optimal;

impl Command for Optimal {
    fn name(&self) -> &str { "OPTIMAL" }
    fn description(&self) -> &str { "Auto-build optimal model using BIC" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        let (max_terms, bic_threshold) = parse_args(args)?;
        let observations = active_observations(&session.observations);
        let n_obs = observations.len();
        if n_obs < BASE_TERMS.len() {
            return Err(Error::Fit("insufficient observations for OPTIMAL".into()));
        }
        let latitude = session.latitude();

        let mut report = String::from("OPTIMAL model search...\n");
        let mut active: Vec<String> = BASE_TERMS.iter().map(|s| s.to_string()).collect();

        let base_fit = try_fit(&observations, &active, latitude)?;
        let mut current_bic = compute_bic(n_obs, active.len(), base_fit.sky_rms);
        append_base_report(&mut report, &active, current_bic, base_fit.sky_rms);

        let mut stage_log = Vec::new();
        current_bic = run_physical_stage(
            &observations, &mut active, current_bic,
            bic_threshold, latitude, &mut stage_log,
        )?;
        let _final_stage_bic = run_harmonic_stage(
            &observations, &mut active, current_bic,
            bic_threshold, max_terms, latitude, &mut stage_log,
        )?;

        for entry in &stage_log {
            report.push_str(&format!(
                "+ {} (dBIC={:.1}, RMS={:.2}\")\n", entry.name, entry.delta_bic, entry.rms,
            ));
        }

        let pruned = prune_terms(&observations, &mut active, latitude)?;
        for name in &pruned {
            report.push_str(&format!("- {} (pruned, significance < {:.1})\n", name, MIN_SIGNIFICANCE));
        }

        let final_fit = try_fit(&observations, &active, latitude)?;
        let _final_bic = compute_bic(n_obs, active.len(), final_fit.sky_rms);
        append_final_report(&mut report, &active, &final_fit);
        load_into_session(session, &active, &final_fit)?;

        Ok(CommandOutput::Text(report))
    }
}

fn parse_args(args: &[&str]) -> Result<(usize, f64)> {
    let max_terms = match args.first() {
        Some(s) => s.parse::<usize>()
            .map_err(|e| Error::Parse(format!("invalid max_terms: {}", e)))?,
        None => DEFAULT_MAX_TERMS,
    };
    let bic_threshold = match args.get(1) {
        Some(s) => s.parse::<f64>()
            .map_err(|e| Error::Parse(format!("invalid bic_threshold: {}", e)))?,
        None => DEFAULT_BIC_THRESHOLD,
    };
    Ok((max_terms, bic_threshold))
}

fn active_observations(observations: &[Observation]) -> Vec<&Observation> {
    observations.iter().filter(|o| !o.masked).collect()
}

fn compute_bic(n_obs: usize, n_terms: usize, sky_rms: f64) -> f64 {
    let n = n_obs as f64;
    let k = n_terms as f64;
    let weighted_rss = sky_rms * sky_rms * n;
    n * (weighted_rss / n).ln() + k * n.ln()
}

fn try_fit(
    observations: &[&Observation],
    term_names: &[String],
    latitude: f64,
) -> Result<FitResult> {
    let terms: Vec<_> = term_names.iter()
        .map(|n| create_term(n))
        .collect::<Result<Vec<_>>>()?;
    let fixed = vec![false; terms.len()];
    let coeffs = vec![0.0; terms.len()];
    fit_model(observations, &terms, &fixed, &coeffs, latitude)
}

fn run_physical_stage(
    observations: &[&Observation],
    active: &mut Vec<String>,
    mut current_bic: f64,
    threshold: f64,
    latitude: f64,
    log: &mut Vec<StageEntry>,
) -> Result<f64> {
    let n_obs = observations.len();
    for &candidate in PHYSICAL_CANDIDATES {
        if active.len() >= DEFAULT_MAX_TERMS {
            break;
        }
        let mut trial = active.clone();
        trial.push(candidate.to_string());
        let fit = match try_fit(observations, &trial, latitude) {
            Ok(f) => f,
            Err(_) => continue,
        };
        let trial_bic = compute_bic(n_obs, trial.len(), fit.sky_rms);
        let delta = trial_bic - current_bic;
        if delta < threshold {
            log.push(StageEntry {
                name: candidate.to_string(),
                delta_bic: delta,
                rms: fit.sky_rms,
            });
            active.push(candidate.to_string());
            current_bic = trial_bic;
        }
    }
    Ok(current_bic)
}

fn generate_harmonic_candidates() -> Vec<String> {
    let results = ["H", "D", "X"];
    let funcs = ["S", "C"];
    let coords = ["H", "D"];
    let mut candidates = Vec::with_capacity(96);
    for r in &results {
        for f in &funcs {
            for c in &coords {
                for n in 1..=8u8 {
                    let suffix = if n == 1 { String::new() } else { n.to_string() };
                    candidates.push(format!("H{}{}{}{}", r, f, c, suffix));
                }
            }
        }
    }
    candidates
}

fn run_harmonic_stage(
    observations: &[&Observation],
    active: &mut Vec<String>,
    mut current_bic: f64,
    threshold: f64,
    max_terms: usize,
    latitude: f64,
    log: &mut Vec<StageEntry>,
) -> Result<f64> {
    let all_candidates = generate_harmonic_candidates();
    let n_obs = observations.len();

    loop {
        if active.len() >= max_terms {
            break;
        }
        let candidates: Vec<&String> = all_candidates.iter()
            .filter(|c| !active.contains(c))
            .collect();
        if candidates.is_empty() {
            break;
        }

        let best = find_best_harmonic(observations, active, &candidates, n_obs, latitude);
        match best {
            Some((name, bic, rms)) => {
                let delta = bic - current_bic;
                if delta < threshold {
                    log.push(StageEntry {
                        name: name.clone(),
                        delta_bic: delta,
                        rms,
                    });
                    active.push(name);
                    current_bic = bic;
                } else {
                    break;
                }
            }
            None => break,
        }
    }
    Ok(current_bic)
}

fn find_best_harmonic(
    observations: &[&Observation],
    active: &[String],
    candidates: &[&String],
    n_obs: usize,
    latitude: f64,
) -> Option<(String, f64, f64)> {
    candidates.par_iter()
        .filter_map(|candidate| {
            let mut trial = active.to_vec();
            trial.push((*candidate).clone());
            let fit = try_fit(observations, &trial, latitude).ok()?;
            let bic = compute_bic(n_obs, trial.len(), fit.sky_rms);
            Some(((*candidate).clone(), bic, fit.sky_rms))
        })
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
}

fn prune_terms(
    observations: &[&Observation],
    active: &mut Vec<String>,
    latitude: f64,
) -> Result<Vec<String>> {
    let fit = try_fit(observations, active, latitude)?;
    let mut pruned = Vec::new();
    let base_set: Vec<String> = BASE_TERMS.iter().map(|s| s.to_string()).collect();

    let to_remove: Vec<String> = active.iter().enumerate()
        .filter(|(i, name)| {
            if base_set.contains(name) {
                return false;
            }
            let sigma = fit.sigma[*i];
            if sigma == 0.0 {
                return false;
            }
            (fit.coefficients[*i] / sigma).abs() < MIN_SIGNIFICANCE
        })
        .map(|(_, name)| name.clone())
        .collect();

    for name in &to_remove {
        active.retain(|n| n != name);
        pruned.push(name.clone());
    }
    Ok(pruned)
}

fn load_into_session(
    session: &mut Session,
    active: &[String],
    result: &FitResult,
) -> Result<()> {
    session.model.remove_all();
    for name in active {
        session.model.add_term(name)?;
    }
    session.model.set_coefficients(&result.coefficients)?;
    session.last_fit = Some(result.clone());
    Ok(())
}

fn append_base_report(report: &mut String, terms: &[String], bic: f64, rms: f64) {
    report.push_str(&format!(
        "Base: {} (BIC={:.1}, RMS={:.2}\")\n",
        terms.join(" "), bic, rms,
    ));
}

fn append_final_report(report: &mut String, terms: &[String], fit: &FitResult) {
    report.push_str(&format!(
        "\nFinal model: {} terms, RMS={:.2}\"\n", terms.len(), fit.sky_rms,
    ));
    report.push_str("Terms: ");
    report.push_str(&terms.join(" "));
    report.push('\n');
}
