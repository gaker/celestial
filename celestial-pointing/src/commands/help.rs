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
        "PARALLEL" => "PARALLEL <term> [term...] | PARALLEL ALL\n  Apply terms in parallel (default)".into(),
        "CHAIN" => "CHAIN <term> [term...] | CHAIN ALL\n  Apply terms sequentially (rigorous)".into(),
        "ADJUST" => "ADJUST T|S\n  T = telescope to star (default)\n  S = star to telescope".into(),
        "FAUTO" => "FAUTO <order> [H|D]\n  Add harmonics up to Nth order\n  H = HA only, D = Dec only".into(),
        "OPTIMAL" => "OPTIMAL [max_terms] [bic_threshold]\n  Auto-build optimal model using BIC selection\n  Defaults: max 30 terms, threshold -6.0".into(),
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
  PARALLEL <terms>   Apply terms in parallel
  CHAIN <terms>      Apply terms sequentially
  ADJUST T|S         Set model direction

  FAUTO <n>          Add harmonics to nth order
  OPTIMAL            Auto-build optimal model
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
