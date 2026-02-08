use crate::error::{Error, Result};
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Inmod;

impl Command for Inmod {
    fn name(&self) -> &str { "INMOD" }
    fn description(&self) -> &str { "Load model from file" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("INMOD requires a filename".into()));
        }
        let content = std::fs::read_to_string(args[0])
            .map_err(Error::Io)?;

        session.model.remove_all();
        session.last_fit = None;

        let mut term_coeffs = Vec::new();
        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.eq_ignore_ascii_case("END") {
                break;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() < 2 {
                return Err(Error::Parse(format!("invalid model line: {}", trimmed)));
            }
            let name = parts[0];
            let coeff: f64 = parts[1].parse()
                .map_err(|e| Error::Parse(format!("invalid coefficient: {}", e)))?;
            session.model.add_term(name)?;
            term_coeffs.push(coeff);
        }

        if !term_coeffs.is_empty() {
            session.model.set_coefficients(&term_coeffs)?;
        }

        Ok(CommandOutput::Text(format!(
            "Loaded {} terms from {}", term_coeffs.len(), args[0]
        )))
    }
}
