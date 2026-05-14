use super::{Command, CommandOutput, FitDisplay};
use crate::error::{Error, Result};
use crate::session::Session;
use crate::solver;

pub struct Fit;

impl Command for Fit {
    fn name(&self) -> &str {
        "FIT"
    }
    
    fn description(&self) -> &str {
        "Fit model to observations"
    }

    fn execute(&self, session: &mut Session, _args: &[&str]) -> Result<CommandOutput> {
        if session.model.term_count() == 0 {
            return raw_rms(session);
        }
        let result = session.fit()?.clone();
        let mut notes: Vec<String> = Vec::new();
        if let Some(n) = ill_conditioning_note(&result.rank_info) {
            notes.push(n);
        }
        if let Some(n) = iteration_note(&result.iter_report) {
            notes.push(n);
        }
        notes.extend(outlier_candidate_notes(&result.diagnostics, result.sky_rms));
        let note = if notes.is_empty() {
            None
        } else {
            Some(notes.join("\n"))
        };
        let fixed = session.model.fixed_flags().to_vec();
        let parallel = session.model.parallel_flags().to_vec();
        Ok(CommandOutput::FitDisplay(FitDisplay {
            term_names: result.term_names,
            coefficients: result.coefficients,
            change: result.change,
            sigma: result.sigma,
            fixed,
            parallel,
            sky_rms: result.sky_rms,
            popn_sd: result.popn_sd,
            note,
        }))
    }
}

pub(crate) const BADNESS_K: f64 = 2.5;

fn ill_conditioning_note(ri: &solver::RankInfo) -> Option<String> {
    if !ri.is_deficient() {
        return None;
    }
    let suppressed = ri.n_free - ri.rank;
    Some(format!(
        "FIT: ill-conditioned at FITTOL={:.0e}; {} of {} direction{} suppressed (cond = {:.2e})",
        ri.tol,
        suppressed,
        ri.n_free,
        if suppressed == 1 { "" } else { "s" },
        ri.condition_number(),
    ))
}

fn iteration_note(report: &solver::IterReport) -> Option<String> {
    if !report.converged {
        return Some(format!(
            "FIT: did not converge after {} iteration{} (final step {:.2}\" on rotational terms)",
            report.iterations,
            if report.iterations == 1 { "" } else { "s" },
            report.final_step_norm * 206264.806_247,
        ));
    }
    if report.iterations >= 2 {
        return Some(format!(
            "FIT: converged in {} iterations",
            report.iterations,
        ));
    }
    None
}

fn outlier_candidate_notes(
    diagnostics: &[solver::ObsDiagnostic],
    sky_rms: f64,
) -> Vec<String> {
    if sky_rms <= 0.0 {
        return Vec::new();
    }
    let cutoff = BADNESS_K * sky_rms;
    diagnostics
        .iter()
        .enumerate()
        .filter(|(_, d)| d.residual_sky > cutoff)
        .map(|(i, d)| {
            format!(
                "Observation #{} is a possible outlier candidate (badness = {:.3}).",
                i + 1,
                d.residual_sky / cutoff,
            )
        })
        .collect()
}

fn raw_rms(session: &Session) -> Result<CommandOutput> {
    let prepared = session.prepared_observations();
    let active: Vec<&_> = prepared.iter().filter(|o| !o.masked).collect();
    if active.is_empty() {
        return Err(Error::Fit("no observations loaded".into()));
    }
    let residuals = solver::build_residuals(&active);
    let rms = solver::compute_sky_rms(&residuals, &active);
    Ok(CommandOutput::Text(format!("Raw sky RMS = {:.2}\"", rms)))
}
