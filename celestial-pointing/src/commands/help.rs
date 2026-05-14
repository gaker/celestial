use super::{Command, CommandOutput};
use crate::error::Result;
use crate::session::Session;

pub struct Help;

impl Command for Help {
    fn name(&self) -> &str {
        "HELP"
    }
    fn description(&self) -> &str {
        "Show available commands"
    }

    fn execute(&self, _session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if let Some(cmd) = args.first() {
            Ok(CommandOutput::Text(command_help(cmd)))
        } else {
            Ok(CommandOutput::Text(general_help()))
        }
    }
}

fn command_help(cmd: &str) -> String {
    match cmd.to_uppercase().as_str() {
        "APPLY" => "APPLY <ra> <dec>\n  Compute commanded encoder position for a target\n  Args: h m s d m s  OR  decimal_hours decimal_degrees".into(),
        "INDAT" => "INDAT <file>\n  Load observations from file".into(),
        "INMOD" => "INMOD <file>\n  Load model from file".into(),
        "OUTMOD" => "OUTMOD <file>\n  Save model to file".into(),
        "USE" => "USE <term> [term...]\n  Add terms to model\n  Example: USE IH ID CH NP MA ME".into(),
        "LOSE" => "LOSE <term> [term...] | LOSE ALL\n  Remove terms from model".into(),
        "FIT" => "FIT\n  Fit model to observations".into(),
        "CLIST" => "CLIST\n  List coefficients with uncertainties".into(),
        "RESET" => "RESET\n  Zero all coefficients".into(),
        "SLIST" => "SLIST\n  List observations with residuals".into(),
        "MASK" => "MASK <obs> [obs...] | MASK <n>-<m>\n  Exclude observations from fit".into(),
        "UNMASK" => "UNMASK <obs> [obs...] | UNMASK ALL\n  Include masked observations".into(),
        "MVET" => "MVET <sigma> [R]\n  Find weak terms (R to remove)".into(),
        "OUTL" => "OUTL <sigma> [M]\n  Find outliers (M to mask)".into(),
        "FIX" => "FIX <term> [term...] | FIX ALL\n  Fix terms at current values during fit".into(),
        "UNFIX" => "UNFIX <term> [term...] | UNFIX ALL\n  Allow fixed terms to be fitted".into(),
        "PARAL" => "PARAL <term> [term...] | PARAL ALL\n  Apply terms in parallel (default)".into(),
        "CHAIN" => "CHAIN <term> [term...] | CHAIN ALL\n  Apply terms sequentially (rigorous)".into(),
        "ADJUST" => "ADJUST T|S\n  T = telescope to star (default)\n  S = star to telescope".into(),
        "FAUTO" => "FAUTO [max_terms] [bic_threshold]\n  Automatic modeling: pair-aware BIC search over harmonics and physical terms\n  Defaults: max 30 terms, threshold -6.0".into(),
        "LST" => "LST [h m s | decimal_hours | CLEAR]\n  Show/set local sidereal time".into(),
        "CORRECT" => "CORRECT <ra> <dec>\n  Compute actual sky position from encoder reading\n  Args: h m s d m s  OR  decimal_hours decimal_degrees".into(),
        "PREDICT" => "PREDICT <ra> <dec>\n  Show per-term correction breakdown\n  Args: h m s d m s  OR  decimal_hours decimal_degrees".into(),
        "GSCAT" => "GSCAT [file.svg]\n  Scatter plot of residuals (dX vs dDec)\n  No args = terminal, with file = SVG output".into(),
        "GDIST" => "GDIST [file.svg] [D]\n  Histogram of residual distribution\n  No args = terminal (both dX and dDec)\n  D = declination residuals (default = dX)".into(),
        "GMAP" => "GMAP [file.svg] [scale]\n  Sky map with residual vectors\n  No args = terminal, scale = arrow scale factor (default 10)".into(),
        "GHA" => "GHA [file.svg]\n  Residuals vs hour angle\n  No args = terminal, with file = two SVGs (_dx, _dd)".into(),
        "GDEC" => "GDEC [file.svg]\n  Residuals vs declination\n  No args = terminal, with file = two SVGs (_dx, _dd)".into(),
        "GHYST" => "GHYST [file.svg]\n  Hysteresis plot (residuals by sequence and pier side)\n  No args = terminal, with file = two SVGs (_east, _west)".into(),
        "SHOW" => "SHOW\n  Display session state".into(),
        "HELP" => "HELP [command]\n  Show help for a command".into(),
        "QUIT" => "QUIT\n  Exit the program".into(),
        _ => format!("Unknown command: {}", cmd),
    }
}

fn general_help() -> String {
    "\
Commands:
  APPLY <ra> <dec>   Compute commanded position for target
  INDAT <file>       Load observations
  INMOD <file>       Load model
  OUTMOD <file>      Save model

  USE <terms>        Add terms to model
  LOSE <terms>       Remove terms (or ALL)
  FIT                Fit model
  CLIST              List coefficients
  RESET              Zero all coefficients

  SLIST              List observations
  MASK <obs>         Exclude observations
  UNMASK <obs>       Include observations
  MVET <sigma>       Find/remove weak terms
  OUTL <sigma>       Find/mask outliers

  FIX <terms>        Fix terms during fit
  UNFIX <terms>      Unfix terms
  PARAL <terms>      Apply terms in parallel
  CHAIN <terms>      Apply terms sequentially
  ADJUST T|S         Set model direction

  FAUTO              Automatic modeling (BIC + pairs + triples)
  LST [time|CLEAR]   Set/show local sidereal time

  CORRECT <ra> <dec> Actual sky position from encoders
  PREDICT <ra> <dec> Per-term correction breakdown

  GSCAT [file]       Scatter plot of residuals
  GDIST [file]       Histogram of residuals
  GMAP [file]        Sky map with residual vectors
  GHA [file]         Residuals vs hour angle
  GDEC [file]        Residuals vs declination
  GHYST [file]       Hysteresis plot

  SHOW               Display session state
  HELP [cmd]         Show help
  QUIT               Exit

Type HELP <command> for details."
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_help(args: &[&str]) -> String {
        let mut session = Session::new();
        match Help.execute(&mut session, args).unwrap() {
            CommandOutput::Text(s) => s,
            other => panic!("expected Text, got {:?}", other),
        }
    }

    #[test]
    fn metadata() {
        assert_eq!(Help.name(), "HELP");
        assert_eq!(Help.description(), "Show available commands");
    }

    #[test]
    fn no_args_returns_general_help() {
        let body = run_help(&[]);
        assert!(body.starts_with("Commands:"));
        assert!(body.contains("Type HELP <command> for details."));
        assert!(body.contains("APPLY"));
        assert!(body.contains("FAUTO"));
    }

    // Every command listed in dispatch() should have a matching help entry.
    // If a new command is added without a help entry, this test catches it.
    #[test]
    fn every_known_command_has_a_help_entry() {
        let commands = [
            "APPLY", "INDAT", "INMOD", "OUTMOD", "USE", "LOSE", "FIT", "CLIST", "RESET",
            "SLIST", "MASK", "UNMASK", "MVET", "OUTL", "FIX", "UNFIX", "PARAL", "CHAIN",
            "ADJUST", "FAUTO", "LST", "CORRECT", "PREDICT", "GSCAT", "GDIST", "GMAP",
            "GHA", "GDEC", "GHYST", "SHOW", "HELP", "QUIT",
        ];
        for cmd in commands {
            let body = run_help(&[cmd]);
            assert!(
                !body.starts_with("Unknown command:"),
                "command {} missing help entry",
                cmd,
            );
            assert!(
                body.contains(cmd),
                "help body for {} should mention the command name: {}",
                cmd,
                body,
            );
        }
    }

    #[test]
    fn command_lookup_is_case_insensitive() {
        for token in ["apply", "Apply", "APPLY", "aPpLy"] {
            let body = run_help(&[token]);
            assert!(body.starts_with("APPLY"), "token {:?} → {}", token, body);
        }
    }

    #[test]
    fn unknown_command_reports_with_original_casing() {
        let body = run_help(&["bogus"]);
        assert_eq!(body, "Unknown command: bogus");
    }

    #[test]
    fn extra_args_are_ignored() {
        let body = run_help(&["FIT", "MORE", "STUFF"]);
        assert!(body.starts_with("FIT"));
        assert!(!body.contains("MORE"));
    }

    #[test]
    fn session_is_not_mutated() {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        let term_count_before = session.model.term_count();
        Help.execute(&mut session, &[]).unwrap();
        Help.execute(&mut session, &["FIT"]).unwrap();
        Help.execute(&mut session, &["bogus"]).unwrap();
        assert_eq!(session.model.term_count(), term_count_before);
    }
}
