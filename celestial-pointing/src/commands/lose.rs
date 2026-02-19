use super::{Command, CommandOutput};
use crate::error::Result;
use crate::session::Session;

pub struct Lose;

impl Command for Lose {
    fn name(&self) -> &str {
        "LOSE"
    }
    fn description(&self) -> &str {
        "Remove term(s) from model"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(crate::error::Error::Parse(
                "LOSE requires term name(s) or ALL".into(),
            ));
        }
        if args.len() == 1 && args[0].eq_ignore_ascii_case("ALL") {
            session.model.remove_all();
            return Ok(CommandOutput::Text("All terms removed".into()));
        }
        for name in args {
            session.model.remove_term(&name.to_uppercase());
        }
        Ok(CommandOutput::Text(format!(
            "Removed: {}",
            args.join(", ").to_uppercase()
        )))
    }
}
