use crate::error::{Error, Result};
use crate::session::{AdjustDirection, Session};
use super::{Command, CommandOutput};

pub struct Adjust;

impl Command for Adjust {
    fn name(&self) -> &str { "ADJUST" }
    fn description(&self) -> &str { "Set model correction direction" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            let current = match session.adjust_direction {
                AdjustDirection::TelescopeToStar => "T (telescope to star)",
                AdjustDirection::StarToTelescope => "S (star to telescope)",
            };
            return Ok(CommandOutput::Text(format!("Current direction: {}", current)));
        }
        match args[0].to_uppercase().as_str() {
            "T" => {
                session.adjust_direction = AdjustDirection::TelescopeToStar;
                Ok(CommandOutput::Text("Direction: telescope to star".to_string()))
            }
            "S" => {
                session.adjust_direction = AdjustDirection::StarToTelescope;
                Ok(CommandOutput::Text("Direction: star to telescope".to_string()))
            }
            _ => Err(Error::Parse(format!(
                "ADJUST requires T or S, got {}", args[0]
            ))),
        }
    }
}
