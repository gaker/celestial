use crate::error::Result;
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Outmod;

impl Command for Outmod {
    fn name(&self) -> &str { "OUTMOD" }
    fn description(&self) -> &str { "Save model to file" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(crate::error::Error::Parse("OUTMOD requires a filename".into()));
        }
        let mut output = String::new();
        for (name, &coeff) in session.model.term_names().iter()
            .zip(session.model.coefficients().iter())
        {
            output += &format!("{} {:.6}\n", name, coeff);
        }
        output += "END\n";
        std::fs::write(args[0], &output).map_err(crate::error::Error::Io)?;
        Ok(CommandOutput::Text(format!("Model saved to {}", args[0])))
    }
}
