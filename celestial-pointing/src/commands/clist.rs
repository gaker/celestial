use super::{Command, CommandOutput, FitDisplay};
use crate::error::Result;
use crate::session::Session;

pub struct Clist;

impl Command for Clist {
    fn name(&self) -> &str {
        "CLIST"
    }
    fn description(&self) -> &str {
        "List current coefficients"
    }

    fn execute(&self, session: &mut Session, _args: &[&str]) -> Result<CommandOutput> {
        let names = session
            .model
            .term_names()
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        if names.is_empty() {
            return Ok(CommandOutput::Text("No terms in model".to_string()));
        }
        let coeffs = session.model.coefficients().to_vec();
        let sigma = session
            .last_fit
            .as_ref()
            .map(|f| f.sigma.clone())
            .unwrap_or_else(|| vec![0.0; names.len()]);
        let sky_rms = session.last_fit.as_ref().map(|f| f.sky_rms).unwrap_or(0.0);

        Ok(CommandOutput::FitDisplay(FitDisplay {
            term_names: names,
            coefficients: coeffs,
            sigma,
            sky_rms,
        }))
    }
}
