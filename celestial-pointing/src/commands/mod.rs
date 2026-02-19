pub mod adjust;
pub mod apply;
pub mod clist;
pub mod correct;
pub mod fauto;
pub mod fit;
pub mod fix;
pub mod gdec;
pub mod gdist;
pub mod gha;
pub mod ghyst;
pub mod gmap;
pub mod gscat;
pub mod help;
pub mod indat;
pub mod inmod;
pub mod lose;
pub mod lst;
pub mod mask;
pub mod mvet;
pub mod optimal;
pub mod outl;
pub mod outmod;
pub mod parallel;
pub mod predict;
pub mod reset;
pub mod show;
pub mod slist;
pub mod use_term;

use crate::error::Result;
use crate::session::Session;

pub enum CommandOutput {
    Text(String),
    Table {
        headers: Vec<String>,
        rows: Vec<Vec<String>>,
    },
    FitDisplay(FitDisplay),
    None,
}

pub struct FitDisplay {
    pub term_names: Vec<String>,
    pub coefficients: Vec<f64>,
    pub sigma: Vec<f64>,
    pub sky_rms: f64,
}

pub trait Command {
    fn name(&self) -> &str;
    fn description(&self) -> &str;
    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput>;
}

pub fn dispatch(session: &mut Session, input: &str) -> Result<CommandOutput> {
    let parts: Vec<&str> = input.split_whitespace().collect();
    if parts.is_empty() {
        return Ok(CommandOutput::None);
    }
    let cmd_name = parts[0].to_uppercase();
    let args = &parts[1..];
    match cmd_name.as_str() {
        "ADJUST" => adjust::Adjust.execute(session, args),
        "APPLY" => apply::Apply.execute(session, args),
        "CHAIN" => parallel::Chain.execute(session, args),
        "CLIST" => clist::Clist.execute(session, args),
        "CORRECT" => correct::Correct.execute(session, args),
        "FAUTO" => fauto::Fauto.execute(session, args),
        "FIT" => fit::Fit.execute(session, args),
        "FIX" => fix::Fix.execute(session, args),
        "GDEC" => gdec::Gdec.execute(session, args),
        "GDIST" => gdist::Gdist.execute(session, args),
        "GHA" => gha::Gha.execute(session, args),
        "GHYST" => ghyst::Ghyst.execute(session, args),
        "GMAP" => gmap::Gmap.execute(session, args),
        "GSCAT" => gscat::Gscat.execute(session, args),
        "HELP" => help::Help.execute(session, args),
        "INDAT" => indat::Indat.execute(session, args),
        "INMOD" => inmod::Inmod.execute(session, args),
        "LOSE" => lose::Lose.execute(session, args),
        "LST" => lst::Lst.execute(session, args),
        "MASK" => mask::Mask.execute(session, args),
        "MVET" => mvet::Mvet.execute(session, args),
        "OPTIMAL" => optimal::Optimal.execute(session, args),
        "OUTL" => outl::Outl.execute(session, args),
        "OUTMOD" => outmod::Outmod.execute(session, args),
        "PARALLEL" => parallel::Parallel.execute(session, args),
        "PREDICT" => predict::Predict.execute(session, args),
        "QUIT" => Ok(CommandOutput::Text("Use Ctrl-D to exit".to_string())),
        "RESET" => reset::Reset.execute(session, args),
        "SHOW" => show::Show.execute(session, args),
        "SLIST" => slist::Slist.execute(session, args),
        "UNFIX" => fix::Unfix.execute(session, args),
        "UNMASK" => mask::Unmask.execute(session, args),
        "USE" => use_term::Use.execute(session, args),
        _ => Err(crate::error::Error::Parse(format!(
            "unknown command: {}",
            parts[0]
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::session::Session;

    #[test]
    fn dispatch_use_adds_terms() {
        let mut session = Session::new();
        let result = dispatch(&mut session, "USE IH ID").unwrap();
        assert_eq!(session.model.term_count(), 2);
        assert_eq!(session.model.term_names(), vec!["IH", "ID"]);
        match result {
            CommandOutput::Text(s) => assert!(s.contains("IH")),
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn dispatch_lose_removes_term() {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        session.model.add_term("ID").unwrap();
        dispatch(&mut session, "LOSE IH").unwrap();
        assert_eq!(session.model.term_count(), 1);
        assert_eq!(session.model.term_names(), vec!["ID"]);
    }

    #[test]
    fn dispatch_unknown_command_errors() {
        let mut session = Session::new();
        let result = dispatch(&mut session, "ZZZNOTACMD");
        assert!(result.is_err());
    }

    #[test]
    fn dispatch_fit_no_observations_errors() {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        let result = dispatch(&mut session, "FIT");
        assert!(result.is_err());
    }

    #[test]
    fn dispatch_empty_input_returns_none() {
        let mut session = Session::new();
        let result = dispatch(&mut session, "").unwrap();
        matches!(result, CommandOutput::None);
    }

    #[test]
    fn dispatch_case_insensitive() {
        let mut session = Session::new();
        let result = dispatch(&mut session, "use IH");
        assert!(result.is_ok());
        assert_eq!(session.model.term_count(), 1);
    }
}
