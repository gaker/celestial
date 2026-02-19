use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;

pub struct Fauto;

impl Command for Fauto {
    fn name(&self) -> &str {
        "FAUTO"
    }
    fn description(&self) -> &str {
        "Auto-add harmonics up to Nth order"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("FAUTO requires a harmonic order".into()));
        }
        let order: usize = args[0]
            .parse()
            .map_err(|e| Error::Parse(format!("invalid order: {}", e)))?;
        if order == 0 {
            return Err(Error::Parse("harmonic order must be >= 1".into()));
        }

        let mut added = Vec::new();
        for n in 1..=order {
            let suffix = if n == 1 { String::new() } else { n.to_string() };
            let names = [format!("HDSH{}", suffix), format!("HDCH{}", suffix)];
            for name in &names {
                if !session.model.term_names().contains(&name.as_str()) {
                    session.model.add_term(name)?;
                    added.push(name.clone());
                }
            }
        }

        if added.is_empty() {
            Ok(CommandOutput::Text(format!(
                "All harmonics up to order {} already in model",
                order
            )))
        } else {
            Ok(CommandOutput::Text(format!(
                "Added {} harmonics: {}",
                added.len(),
                added.join(" ")
            )))
        }
    }
}
