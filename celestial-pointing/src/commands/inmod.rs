use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;

pub struct Inmod;

impl Command for Inmod {
    fn name(&self) -> &str {
        "INMOD"
    }
    fn description(&self) -> &str {
        "Load model from file"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("INMOD requires a filename".into()));
        }
        let content = std::fs::read_to_string(args[0]).map_err(Error::Io)?;

        session.model.remove_all();
        session.last_fit = None;

        let mut term_coeffs = Vec::new();
        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.eq_ignore_ascii_case("END") {
                break;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() < 2 {
                return Err(Error::Parse(format!("invalid model line: {}", trimmed)));
            }
            let name = parts[0];
            let coeff: f64 = parts[1]
                .parse()
                .map_err(|e| Error::Parse(format!("invalid coefficient: {}", e)))?;
            session.model.add_term(name)?;
            term_coeffs.push(coeff);
        }

        if !term_coeffs.is_empty() {
            session.model.set_coefficients(&term_coeffs)?;
        }

        Ok(CommandOutput::Text(format!(
            "Loaded {} terms from {}",
            term_coeffs.len(),
            args[0]
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::FitResult;
    use crate::test_support::FitResultBuilder;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static SEQ: AtomicU64 = AtomicU64::new(0);

    fn write_tmp(contents: &str, tag: &str) -> PathBuf {
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let mut path = std::env::temp_dir();
        path.push(format!("celpoint-inmod-{}-{}-{}.mod", tag, std::process::id(), n));
        std::fs::write(&path, contents).expect("write temp model");
        path
    }

    fn cleanup(p: &PathBuf) {
        let _ = std::fs::remove_file(p);
    }

    fn fake_fit() -> FitResult {
        FitResultBuilder::new()
            .coefficients(vec![1.0])
            .sky_rms(1.0)
            .popn_sd(1.0)
            .term_names(["IH"])
            .build()
    }

    fn run(session: &mut Session, contents: &str, tag: &str) -> Result<String> {
        let p = write_tmp(contents, tag);
        let result = Inmod
            .execute(session, &[p.to_str().unwrap()])
            .map(|out| match out {
                CommandOutput::Text(s) => s,
                other => panic!("expected Text, got {:?}", other),
            });
        cleanup(&p);
        result
    }

    #[test]
    fn metadata() {
        assert_eq!(Inmod.name(), "INMOD");
        assert_eq!(Inmod.description(), "Load model from file");
    }

    #[test]
    fn no_args_errors() {
        let mut session = Session::new();
        let err = Inmod.execute(&mut session, &[]).unwrap_err();
        match err {
            Error::Parse(msg) => assert!(msg.contains("INMOD requires a filename")),
            other => panic!("expected Parse, got {:?}", other),
        }
    }

    #[test]
    fn missing_file_returns_io_error() {
        let mut session = Session::new();
        let err = Inmod
            .execute(&mut session, &["/no/such/file/at/all.mod"])
            .unwrap_err();
        assert!(matches!(err, Error::Io(_)));
    }

    #[test]
    fn loads_terms_and_coefficients() {
        let mut session = Session::new();
        let body = run(&mut session, "IH 10.5\nID -20.0\n", "basic").unwrap();
        assert_eq!(session.model.term_names(), vec!["IH", "ID"]);
        assert_eq!(session.model.coefficients(), &[10.5, -20.0]);
        assert!(body.starts_with("Loaded 2 terms from"));
    }

    #[test]
    fn empty_file_loads_zero_terms() {
        let mut session = Session::new();
        let body = run(&mut session, "", "empty").unwrap();
        assert_eq!(session.model.term_count(), 0);
        assert!(body.starts_with("Loaded 0 terms from"));
    }

    #[test]
    fn whitespace_only_lines_are_skipped() {
        let mut session = Session::new();
        run(&mut session, "\n   \nIH 1.0\n\n", "blank").unwrap();
        assert_eq!(session.model.term_names(), vec!["IH"]);
        assert_eq!(session.model.coefficients(), &[1.0]);
    }

    #[test]
    fn end_marker_terminates_parsing() {
        let mut session = Session::new();
        run(
            &mut session,
            "IH 1.0\nID 2.0\nEND\nCH 999.0\nignored garbage\n",
            "end",
        )
        .unwrap();
        assert_eq!(session.model.term_names(), vec!["IH", "ID"]);
    }

    #[test]
    fn end_marker_is_case_insensitive() {
        let mut session = Session::new();
        run(&mut session, "IH 1.0\nend\nCH 2.0\n", "end_lower").unwrap();
        assert_eq!(session.model.term_count(), 1);
    }

    #[test]
    fn previous_model_and_fit_are_cleared_before_load() {
        let mut session = Session::new();
        session.model.add_term("AW").unwrap();
        session.last_fit = Some(fake_fit());
        run(&mut session, "IH 5.0\n", "reset").unwrap();
        assert_eq!(session.model.term_names(), vec!["IH"]);
        assert!(session.last_fit.is_none());
    }

    #[test]
    fn line_with_only_one_token_errors() {
        let mut session = Session::new();
        let err = run(&mut session, "IH\n", "short").unwrap_err();
        match err {
            Error::Parse(msg) => assert!(msg.contains("invalid model line")),
            other => panic!("expected Parse, got {:?}", other),
        }
    }

    #[test]
    fn unparseable_coefficient_errors() {
        let mut session = Session::new();
        let err = run(&mut session, "IH not_a_number\n", "bad_coef").unwrap_err();
        match err {
            Error::Parse(msg) => assert!(msg.contains("invalid coefficient")),
            other => panic!("expected Parse, got {:?}", other),
        }
    }

    #[test]
    fn unknown_term_name_propagates_error() {
        let mut session = Session::new();
        let err = run(&mut session, "NOSUCHTERM 1.0\n", "bad_term").unwrap_err();
        // add_term failure surfaces as Error::UnknownTerm, not Parse.
        assert!(matches!(err, Error::UnknownTerm(_)));
    }

    // Even when parsing fails partway through, the model has already been
    // cleared. This is intentional — INMOD is a "replace everything" op.
    #[test]
    fn parse_failure_after_clearing_leaves_model_in_intermediate_state() {
        let mut session = Session::new();
        session.model.add_term("AW").unwrap();
        let err = run(&mut session, "IH 1.0\nbroken\n", "midfail").unwrap_err();
        assert!(matches!(err, Error::Parse(_)));
        assert!(!session.model.term_names().contains(&"AW"));
    }

    #[test]
    fn extra_columns_after_coefficient_are_ignored() {
        let mut session = Session::new();
        run(&mut session, "IH 3.0 extra junk fields\n", "extra").unwrap();
        assert_eq!(session.model.coefficients(), &[3.0]);
    }
}
