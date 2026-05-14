use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::{AdjustDirection, Session};

pub struct Adjust;

impl Command for Adjust {
    fn name(&self) -> &str {
        "ADJUST"
    }
    
    fn description(&self) -> &str {
        "Set model correction direction"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            let current = match session.adjust_direction {
                AdjustDirection::TelescopeToStar => "T (telescope to star)",
                AdjustDirection::StarToTelescope => "S (star to telescope)",
            };
            return Ok(CommandOutput::Text(format!(
                "Current direction: {}",
                current
            )));
        }
        match args[0].to_uppercase().as_str() {
            "T" => {
                session.adjust_direction = AdjustDirection::TelescopeToStar;
                Ok(CommandOutput::Text(
                    "Direction: telescope to star".to_string(),
                ))
            }
            "S" => {
                session.adjust_direction = AdjustDirection::StarToTelescope;
                Ok(CommandOutput::Text(
                    "Direction: star to telescope".to_string(),
                ))
            }
            _ => Err(Error::Parse(format!(
                "ADJUST requires T or S, got {}",
                args[0]
            ))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unwrap_text(output: CommandOutput) -> String {
        match output {
            CommandOutput::Text(s) => s,
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn metadata() {
        assert_eq!(Adjust.name(), "ADJUST");
        assert_eq!(Adjust.description(), "Set model correction direction");
    }

    #[test]
    fn no_args_reports_default_direction() {
        let mut session = Session::new();
        assert_eq!(
            session.adjust_direction,
            AdjustDirection::TelescopeToStar,
            "default direction should be telescope-to-star",
        );
        let text = unwrap_text(Adjust.execute(&mut session, &[]).unwrap());
        assert!(text.contains("Current direction:"));
        assert!(text.contains("T (telescope to star)"));
        assert_eq!(session.adjust_direction, AdjustDirection::TelescopeToStar);
    }

    #[test]
    fn no_args_reports_star_to_telescope_when_set() {
        let mut session = Session::new();
        session.adjust_direction = AdjustDirection::StarToTelescope;
        let text = unwrap_text(Adjust.execute(&mut session, &[]).unwrap());
        assert!(text.contains("S (star to telescope)"));
        assert_eq!(session.adjust_direction, AdjustDirection::StarToTelescope);
    }

    #[test]
    fn sets_direction_for_each_valid_token() {
        let cases = [
            ("T", AdjustDirection::TelescopeToStar, "telescope to star"),
            ("t", AdjustDirection::TelescopeToStar, "telescope to star"),
            ("S", AdjustDirection::StarToTelescope, "star to telescope"),
            ("s", AdjustDirection::StarToTelescope, "star to telescope"),
        ];
        for (arg, expected_dir, expected_phrase) in cases {
            let mut session = Session::new();
            let text = unwrap_text(Adjust.execute(&mut session, &[arg]).unwrap());
            assert_eq!(
                session.adjust_direction, expected_dir,
                "arg {:?} should set direction to {:?}",
                arg, expected_dir,
            );
            assert!(
                text.contains(expected_phrase),
                "arg {:?} produced unexpected text: {}",
                arg,
                text,
            );
        }
    }

    #[test]
    fn switching_direction_overwrites_previous_value() {
        let mut session = Session::new();
        Adjust.execute(&mut session, &["S"]).unwrap();
        assert_eq!(session.adjust_direction, AdjustDirection::StarToTelescope);
        Adjust.execute(&mut session, &["T"]).unwrap();
        assert_eq!(session.adjust_direction, AdjustDirection::TelescopeToStar);
    }

    #[test]
    fn invalid_token_errors_and_leaves_direction_unchanged() {
        let mut session = Session::new();
        session.adjust_direction = AdjustDirection::StarToTelescope;
        for arg in ["X", "TS", "", "telescope", "1"] {
            let err = Adjust.execute(&mut session, &[arg]).unwrap_err();
            match err {
                Error::Parse(msg) => {
                    assert!(msg.contains("ADJUST requires T or S"));
                    assert!(
                        msg.contains(arg),
                        "error message {:?} should echo offending arg {:?}",
                        msg,
                        arg,
                    );
                }
                other => panic!("expected Error::Parse, got {:?}", other),
            }
            assert_eq!(
                session.adjust_direction,
                AdjustDirection::StarToTelescope,
                "failed parse must not mutate session state",
            );
        }
    }

    #[test]
    fn only_first_arg_is_inspected() {
        let mut session = Session::new();
        Adjust
            .execute(&mut session, &["S", "garbage", "T"])
            .unwrap();
        assert_eq!(session.adjust_direction, AdjustDirection::StarToTelescope);
    }
}
