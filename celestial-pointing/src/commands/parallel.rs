use crate::error::{Error, Result};
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Parallel;
pub struct Chain;

impl Command for Parallel {
    fn name(&self) -> &str { "PARALLEL" }
    fn description(&self) -> &str { "Apply terms in parallel" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("PARALLEL requires term names or ALL".into()));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            session.model.set_all_parallel();
            return Ok(CommandOutput::Text(format!(
                "All {} terms set to parallel", session.model.term_count()
            )));
        }
        let mut set = Vec::new();
        for name in args {
            let upper = name.to_uppercase();
            if session.model.set_parallel(&upper) {
                set.push(upper);
            } else {
                return Err(Error::Parse(format!("term {} not in model", name)));
            }
        }
        Ok(CommandOutput::Text(format!("Parallel: {}", set.join(" "))))
    }
}

impl Command for Chain {
    fn name(&self) -> &str { "CHAIN" }
    fn description(&self) -> &str { "Apply terms sequentially (chained)" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("CHAIN requires term names or ALL".into()));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            session.model.set_all_chained();
            return Ok(CommandOutput::Text(format!(
                "All {} terms set to chained", session.model.term_count()
            )));
        }
        let mut set = Vec::new();
        for name in args {
            let upper = name.to_uppercase();
            if session.model.set_chained(&upper) {
                set.push(upper);
            } else {
                return Err(Error::Parse(format!("term {} not in model", name)));
            }
        }
        Ok(CommandOutput::Text(format!("Chained: {}", set.join(" "))))
    }
}
