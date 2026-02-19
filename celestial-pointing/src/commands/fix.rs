use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;

pub struct Fix;
pub struct Unfix;

impl Command for Fix {
    fn name(&self) -> &str {
        "FIX"
    }
    fn description(&self) -> &str {
        "Fix terms at current values during fit"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("FIX requires term names or ALL".into()));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            session.model.fix_all();
            return Ok(CommandOutput::Text(format!(
                "Fixed all {} terms",
                session.model.term_count()
            )));
        }
        let mut fixed = Vec::new();
        for name in args {
            let upper = name.to_uppercase();
            if session.model.fix_term(&upper) {
                fixed.push(upper);
            } else {
                return Err(Error::Parse(format!("term {} not in model", name)));
            }
        }
        Ok(CommandOutput::Text(format!("Fixed: {}", fixed.join(" "))))
    }
}

impl Command for Unfix {
    fn name(&self) -> &str {
        "UNFIX"
    }
    fn description(&self) -> &str {
        "Allow terms to be fitted"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("UNFIX requires term names or ALL".into()));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            session.model.unfix_all();
            return Ok(CommandOutput::Text(format!(
                "Unfixed all {} terms",
                session.model.term_count()
            )));
        }
        let mut unfixed = Vec::new();
        for name in args {
            let upper = name.to_uppercase();
            if session.model.unfix_term(&upper) {
                unfixed.push(upper);
            } else {
                return Err(Error::Parse(format!("term {} not in model", name)));
            }
        }
        Ok(CommandOutput::Text(format!(
            "Unfixed: {}",
            unfixed.join(" ")
        )))
    }
}
