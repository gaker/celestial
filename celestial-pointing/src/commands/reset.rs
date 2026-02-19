use super::{Command, CommandOutput};
use crate::error::Result;
use crate::session::Session;

pub struct Reset;

impl Command for Reset {
    fn name(&self) -> &str {
        "RESET"
    }
    fn description(&self) -> &str {
        "Zero all coefficients"
    }

    fn execute(&self, session: &mut Session, _args: &[&str]) -> Result<CommandOutput> {
        session.model.zero_coefficients();
        session.last_fit = None;
        let count = session.model.term_count();
        Ok(CommandOutput::Text(format!(
            "Reset {} coefficients to zero",
            count
        )))
    }
}
