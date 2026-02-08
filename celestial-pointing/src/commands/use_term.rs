use crate::error::Result;
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Use;

impl Command for Use {
    fn name(&self) -> &str { "USE" }
    fn description(&self) -> &str { "Add term(s) to model" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(crate::error::Error::Parse("USE requires term name(s)".into()));
        }
        let mut added = Vec::new();
        for name in args {
            session.model.add_term(name)?;
            added.push(name.to_uppercase());
        }
        Ok(CommandOutput::Text(format!("Added: {}", added.join(", "))))
    }
}
