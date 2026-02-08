use crate::error::{Error, Result};
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Mask;
pub struct Unmask;

impl Command for Mask {
    fn name(&self) -> &str { "MASK" }
    fn description(&self) -> &str { "Mask observations (exclude from fit)" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("MASK requires observation numbers".into()));
        }
        let indices = parse_obs_indices(args, session.observations.len())?;
        let mut count = 0;
        for idx in &indices {
            if !session.observations[*idx].masked {
                session.observations[*idx].masked = true;
                count += 1;
            }
        }
        Ok(CommandOutput::Text(format!("Masked {} observations", count)))
    }
}

impl Command for Unmask {
    fn name(&self) -> &str { "UNMASK" }
    fn description(&self) -> &str { "Unmask observations (include in fit)" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("UNMASK requires observation numbers or ALL".into()));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            let count = session.observations.iter().filter(|o| o.masked).count();
            for obs in &mut session.observations {
                obs.masked = false;
            }
            return Ok(CommandOutput::Text(format!("Unmasked {} observations", count)));
        }
        let indices = parse_obs_indices(args, session.observations.len())?;
        let mut count = 0;
        for idx in &indices {
            if session.observations[*idx].masked {
                session.observations[*idx].masked = false;
                count += 1;
            }
        }
        Ok(CommandOutput::Text(format!("Unmasked {} observations", count)))
    }
}

fn parse_obs_indices(args: &[&str], total: usize) -> Result<Vec<usize>> {
    let mut indices = Vec::new();
    for arg in args {
        if arg.contains('-') {
            let parts: Vec<&str> = arg.splitn(2, '-').collect();
            let start: usize = parts[0].parse()
                .map_err(|e| Error::Parse(format!("invalid range start: {}", e)))?;
            let end: usize = parts[1].parse()
                .map_err(|e| Error::Parse(format!("invalid range end: {}", e)))?;
            if start < 1 || end < 1 || start > total || end > total {
                return Err(Error::Parse(format!(
                    "range {}-{} out of bounds (1-{})", start, end, total
                )));
            }
            for i in start..=end {
                indices.push(i - 1);
            }
        } else {
            let num: usize = arg.parse()
                .map_err(|e| Error::Parse(format!("invalid observation number: {}", e)))?;
            if num < 1 || num > total {
                return Err(Error::Parse(format!(
                    "observation {} out of bounds (1-{})", num, total
                )));
            }
            indices.push(num - 1);
        }
    }
    Ok(indices)
}
