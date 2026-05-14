use super::{Command, CommandOutput};
use crate::error::Result;
use crate::session::Session;

pub struct Reset;

impl Command for Reset {
    fn name(&self) -> &str {
        "RESET"
    }
    fn description(&self) -> &str {
        "Zero all coefficients"
    }

    fn execute(&self, session: &mut Session, _args: &[&str]) -> Result<CommandOutput> {
        session.model.zero_coefficients();
        session.last_fit = None;
        let count = session.model.term_count();
        Ok(CommandOutput::Text(format!(
            "Reset {} coefficients to zero",
            count
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::FitResult;
    use crate::test_support::FitResultBuilder;

    fn text(out: CommandOutput) -> String {
        match out {
            CommandOutput::Text(s) => s,
            other => panic!("expected Text, got {:?}", other),
        }
    }

    fn fake_fit() -> FitResult {
        FitResultBuilder::new()
            .coefficients(vec![1.0])
            .sky_rms(1.0)
            .popn_sd(1.0)
            .term_names(["IH"])
            .build()
    }

    #[test]
    fn metadata() {
        assert_eq!(Reset.name(), "RESET");
        assert_eq!(Reset.description(), "Zero all coefficients");
    }

    #[test]
    fn empty_model_reports_zero() {
        let mut s = Session::new();
        let body = text(Reset.execute(&mut s, &[]).unwrap());
        assert_eq!(body, "Reset 0 coefficients to zero");
    }

    #[test]
    fn zeros_coefficients_and_keeps_terms() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.model.add_term("ID").unwrap();
        s.model.set_coefficients(&[10.0, -20.0]).unwrap();

        let body = text(Reset.execute(&mut s, &[]).unwrap());

        assert_eq!(body, "Reset 2 coefficients to zero");
        assert_eq!(s.model.term_names(), vec!["IH", "ID"]);
        assert_eq!(s.model.coefficients(), &[0.0, 0.0]);
    }

    #[test]
    fn clears_last_fit() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.last_fit = Some(fake_fit());
        Reset.execute(&mut s, &[]).unwrap();
        assert!(s.last_fit.is_none());
    }

    #[test]
    fn args_are_ignored() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.model.set_coefficients(&[7.0]).unwrap();
        text(Reset.execute(&mut s, &["junk", "ALL", "ignored"]).unwrap());
        assert_eq!(s.model.coefficients(), &[0.0]);
    }

    #[test]
    fn idempotent() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.model.set_coefficients(&[5.0]).unwrap();
        Reset.execute(&mut s, &[]).unwrap();
        let body = text(Reset.execute(&mut s, &[]).unwrap());
        assert_eq!(body, "Reset 1 coefficients to zero");
        assert_eq!(s.model.coefficients(), &[0.0]);
    }
}
