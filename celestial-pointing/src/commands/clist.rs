use super::{Command, CommandOutput, FitDisplay};
use crate::error::Result;
use crate::session::Session;

pub struct Clist;

impl Command for Clist {
    fn name(&self) -> &str {
        "CLIST"
    }
    
    fn description(&self) -> &str {
        "List current coefficients"
    }

    fn execute(&self, session: &mut Session, _args: &[&str]) -> Result<CommandOutput> {
        let names = session
            .model
            .term_names()
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        if names.is_empty() {
            return Ok(CommandOutput::Text("No terms in model".to_string()));
        }
        let coeffs = session.model.coefficients().to_vec();
        let sigma = session
            .last_fit
            .as_ref()
            .map(|f| f.sigma.clone())
            .unwrap_or_else(|| vec![0.0; names.len()]);
        let sky_rms = session.last_fit.as_ref().map(|f| f.sky_rms).unwrap_or(0.0);
        let popn_sd = session.last_fit.as_ref().map(|f| f.popn_sd).unwrap_or(0.0);
        let fixed = session.model.fixed_flags().to_vec();
        let parallel = session.model.parallel_flags().to_vec();
        let change = vec![0.0; names.len()];

        Ok(CommandOutput::FitDisplay(FitDisplay {
            term_names: names,
            coefficients: coeffs,
            change,
            sigma,
            fixed,
            parallel,
            sky_rms,
            popn_sd,
            note: None,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::FitResult;
    use crate::test_support::FitResultBuilder;

    fn unwrap_fit(output: CommandOutput) -> FitDisplay {
        match output {
            CommandOutput::FitDisplay(fd) => fd,
            other => panic!("expected FitDisplay, got {:?}", other),
        }
    }

    fn fake_fit(sigma: Vec<f64>, sky_rms: f64, popn_sd: f64) -> FitResult {
        FitResultBuilder::new()
            .sigma(sigma)
            .sky_rms(sky_rms)
            .popn_sd(popn_sd)
            .build()
    }

    #[test]
    fn metadata() {
        assert_eq!(Clist.name(), "CLIST");
        assert_eq!(Clist.description(), "List current coefficients");
    }

    #[test]
    fn empty_model_returns_text_message() {
        let mut session = Session::new();
        match Clist.execute(&mut session, &[]).unwrap() {
            CommandOutput::Text(s) => assert_eq!(s, "No terms in model"),
            other => panic!("expected Text, got {:?}", other),
        }
    }

    #[test]
    fn args_are_ignored() {
        let mut session = Session::new();
        let with_args = Clist.execute(&mut session, &["junk", "more"]).unwrap();
        match with_args {
            CommandOutput::Text(s) => assert_eq!(s, "No terms in model"),
            other => panic!("expected Text, got {:?}", other),
        }
    }

    #[test]
    fn no_fit_produces_zero_sigma_and_zero_rms() {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        session.model.add_term("ID").unwrap();
        session.model.set_coefficients(&[10.0, -20.0]).unwrap();

        let fd = unwrap_fit(Clist.execute(&mut session, &[]).unwrap());

        assert_eq!(fd.term_names, vec!["IH", "ID"]);
        assert_eq!(fd.coefficients, vec![10.0, -20.0]);
        assert_eq!(fd.sigma, vec![0.0, 0.0]);
        assert_eq!(fd.change, vec![0.0, 0.0]);
        assert_eq!(fd.fixed, vec![false, false]);
        assert_eq!(fd.parallel, vec![false, false]);
        assert_eq!(fd.sky_rms, 0.0);
        assert_eq!(fd.popn_sd, 0.0);
        assert!(fd.note.is_none());
    }

    #[test]
    fn last_fit_supplies_sigma_and_summary_stats() {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        session.model.add_term("ID").unwrap();
        session.model.set_coefficients(&[1.5, 2.5]).unwrap();
        session.last_fit = Some(fake_fit(vec![0.25, 0.75], 1.23, 4.56));

        let fd = unwrap_fit(Clist.execute(&mut session, &[]).unwrap());

        assert_eq!(fd.sigma, vec![0.25, 0.75]);
        assert_eq!(fd.sky_rms, 1.23);
        assert_eq!(fd.popn_sd, 4.56);
        assert_eq!(fd.coefficients, vec![1.5, 2.5]);
    }

    #[test]
    fn fixed_and_parallel_flags_are_surfaced() {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        session.model.add_term("ID").unwrap();
        session.model.add_term("CH").unwrap();
        assert!(session.model.fix_term("ID"));
        assert!(session.model.set_parallel("CH"));

        let fd = unwrap_fit(Clist.execute(&mut session, &[]).unwrap());

        assert_eq!(fd.fixed, vec![false, true, false]);
        assert_eq!(fd.parallel, vec![false, false, true]);
    }

    #[test]
    fn change_vector_always_zero_for_clist() {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        session.last_fit = Some(fake_fit(vec![0.1], 0.0, 0.0));

        let fd = unwrap_fit(Clist.execute(&mut session, &[]).unwrap());

        assert_eq!(fd.change, vec![0.0]);
    }
}
