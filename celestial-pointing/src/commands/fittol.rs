use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;

pub struct Fittol;

impl Command for Fittol {
    fn name(&self) -> &str {
        "FITTOL"
    }
    
    fn description(&self) -> &str {
        "Set or report the SVD ill-conditioning tolerance"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Ok(CommandOutput::Text(format!(
                "FITTOL = {:.3e}",
                session.fit_tol
            )));
        }
        let v: f64 = args[0]
            .parse()
            .map_err(|e| Error::Parse(format!("invalid FITTOL value: {}", e)))?;
        if !v.is_finite() || v < 0.0 {
            return Err(Error::Parse(
                "FITTOL must be a non-negative finite number".into(),
            ));
        }
        if v > 0.01 {
            return Err(Error::Parse(
                "FITTOL above 0.01 risks excluding respectable terms".into(),
            ));
        }
        session.fit_tol = v;
        let msg = if v == 0.0 {
            "FITTOL = 0 (ill-conditioning protection disabled)".to_string()
        } else {
            format!("FITTOL = {:.3e}", v)
        };
        Ok(CommandOutput::Text(msg))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::session::DEFAULT_FIT_TOL;

    #[test]
    fn report_default() {
        let mut session = Session::new();
        let result = Fittol.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => assert!(s.contains("1.000e-3")),
            _ => panic!("expected text"),
        }
        assert_eq!(session.fit_tol, DEFAULT_FIT_TOL);
    }

    #[test]
    fn set_zero_disables() {
        let mut session = Session::new();
        Fittol.execute(&mut session, &["0"]).unwrap();
        assert_eq!(session.fit_tol, 0.0);
    }

    #[test]
    fn set_valid_value() {
        let mut session = Session::new();
        Fittol.execute(&mut session, &["0.005"]).unwrap();
        assert_eq!(session.fit_tol, 0.005);
    }

    #[test]
    fn reject_negative() {
        let mut session = Session::new();
        assert!(Fittol.execute(&mut session, &["-0.001"]).is_err());
    }

    #[test]
    fn reject_above_limit() {
        let mut session = Session::new();
        assert!(Fittol.execute(&mut session, &["0.1"]).is_err());
    }

    #[test]
    fn reject_nan() {
        let mut session = Session::new();
        assert!(Fittol.execute(&mut session, &["nan"]).is_err());
    }

    #[test]
    fn reject_garbage() {
        let mut session = Session::new();
        assert!(Fittol.execute(&mut session, &["abc"]).is_err());
    }
}
