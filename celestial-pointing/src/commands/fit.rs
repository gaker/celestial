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
        let result = session.fit()?;
        Ok(CommandOutput::FitDisplay(FitDisplay {
            term_names: result.term_names.clone(),
            coefficients: result.coefficients.clone(),
            sigma: result.sigma.clone(),
            sky_rms: result.sky_rms,
        }))
    }
}

fn raw_rms(session: &Session) -> Result<CommandOutput> {
    let active: Vec<&_> = session.observations.iter().filter(|o| !o.masked).collect();
    if active.is_empty() {
        return Err(Error::Fit("no observations loaded".into()));
    }
    let residuals = solver::build_residuals(&active);
    let rms = solver::compute_sky_rms(&residuals, &active);
    Ok(CommandOutput::Text(format!("Raw sky RMS = {:.2}\"", rms)))
}
