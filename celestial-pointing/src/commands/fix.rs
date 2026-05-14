use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;

pub struct Fix;
pub struct Unfix;

impl Command for Fix {
    fn name(&self) -> &str {
        "FIX"
    }
    fn description(&self) -> &str {
        "Fix terms at current values during fit"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("FIX requires term names or ALL".into()));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            session.model.fix_all();
            return Ok(CommandOutput::Text(format!(
                "Fixed all {} terms",
                session.model.term_count()
            )));
        }
        let mut fixed = Vec::new();
        for name in args {
            let upper = name.to_uppercase();
            if session.model.fix_term(&upper) {
                fixed.push(upper);
            } else {
                return Err(Error::Parse(format!("term {} not in model", name)));
            }
        }
        Ok(CommandOutput::Text(format!("Fixed: {}", fixed.join(" "))))
    }
}

impl Command for Unfix {
    fn name(&self) -> &str {
        "UNFIX"
    }
    fn description(&self) -> &str {
        "Allow terms to be fitted"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("UNFIX requires term names or ALL".into()));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            session.model.unfix_all();
            return Ok(CommandOutput::Text(format!(
                "Unfixed all {} terms",
                session.model.term_count()
            )));
        }
        let mut unfixed = Vec::new();
        for name in args {
            let upper = name.to_uppercase();
            if session.model.unfix_term(&upper) {
                unfixed.push(upper);
            } else {
                return Err(Error::Parse(format!("term {} not in model", name)));
            }
        }
        Ok(CommandOutput::Text(format!(
            "Unfixed: {}",
            unfixed.join(" ")
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn session_with_terms(names: &[&str]) -> Session {
        let mut s = Session::new();
        for name in names {
            s.model.add_term(name).unwrap();
        }
        s
    }

    fn text_body(out: CommandOutput) -> String {
        match out {
            CommandOutput::Text(s) => s,
            other => panic!("expected Text, got {:?}", other),
        }
    }

    fn parse_message(err: &Error) -> &str {
        match err {
            Error::Parse(msg) => msg.as_str(),
            other => panic!("expected Error::Parse, got {:?}", other),
        }
    }

    #[test]
    fn fix_metadata() {
        assert_eq!(Fix.name(), "FIX");
        assert_eq!(Fix.description(), "Fix terms at current values during fit");
    }

    #[test]
    fn unfix_metadata() {
        assert_eq!(Unfix.name(), "UNFIX");
        assert_eq!(Unfix.description(), "Allow terms to be fitted");
    }

    #[test]
    fn fix_no_args_errors() {
        let mut s = session_with_terms(&["IH"]);
        let err = Fix.execute(&mut s, &[]).unwrap_err();
        assert!(parse_message(&err).contains("FIX requires term names or ALL"));
        assert_eq!(s.model.fixed_flags(), &[false]);
    }

    #[test]
    fn unfix_no_args_errors() {
        let mut s = session_with_terms(&["IH"]);
        s.model.fix_all();
        let err = Unfix.execute(&mut s, &[]).unwrap_err();
        assert!(parse_message(&err).contains("UNFIX requires term names or ALL"));
        assert_eq!(s.model.fixed_flags(), &[true]);
    }

    #[test]
    fn fix_all_marks_every_term_regardless_of_case() {
        for token in ["ALL", "all", "All", "aLL"] {
            let mut s = session_with_terms(&["IH", "ID", "CH"]);
            let body = text_body(Fix.execute(&mut s, &[token]).unwrap());
            assert_eq!(body, "Fixed all 3 terms", "token {:?}", token);
            assert_eq!(s.model.fixed_flags(), &[true, true, true]);
        }
    }

    #[test]
    fn unfix_all_clears_every_term_regardless_of_case() {
        for token in ["ALL", "all", "All"] {
            let mut s = session_with_terms(&["IH", "ID"]);
            s.model.fix_all();
            let body = text_body(Unfix.execute(&mut s, &[token]).unwrap());
            assert_eq!(body, "Unfixed all 2 terms", "token {:?}", token);
            assert_eq!(s.model.fixed_flags(), &[false, false]);
        }
    }

    #[test]
    fn fix_all_on_empty_model_reports_zero() {
        let mut s = Session::new();
        let body = text_body(Fix.execute(&mut s, &["ALL"]).unwrap());
        assert_eq!(body, "Fixed all 0 terms");
    }

    #[test]
    fn unfix_all_on_empty_model_reports_zero() {
        let mut s = Session::new();
        let body = text_body(Unfix.execute(&mut s, &["ALL"]).unwrap());
        assert_eq!(body, "Unfixed all 0 terms");
    }

    #[test]
    fn fix_named_terms_uppercases_input() {
        let mut s = session_with_terms(&["IH", "ID", "CH"]);
        let body = text_body(Fix.execute(&mut s, &["ih", "Id"]).unwrap());
        assert_eq!(body, "Fixed: IH ID");
        assert_eq!(s.model.fixed_flags(), &[true, true, false]);
    }

    #[test]
    fn unfix_named_terms_uppercases_input() {
        let mut s = session_with_terms(&["IH", "ID", "CH"]);
        s.model.fix_all();
        let body = text_body(Unfix.execute(&mut s, &["ih", "ch"]).unwrap());
        assert_eq!(body, "Unfixed: IH CH");
        assert_eq!(s.model.fixed_flags(), &[false, true, false]);
    }

    #[test]
    fn fix_unknown_term_errors_with_original_casing_in_message() {
        let mut s = session_with_terms(&["IH"]);
        let err = Fix.execute(&mut s, &["bogus"]).unwrap_err();
        let msg = parse_message(&err);
        assert!(msg.contains("term bogus not in model"), "got {}", msg);
    }

    #[test]
    fn unfix_unknown_term_errors_with_original_casing_in_message() {
        let mut s = session_with_terms(&["IH"]);
        let err = Unfix.execute(&mut s, &["bogus"]).unwrap_err();
        let msg = parse_message(&err);
        assert!(msg.contains("term bogus not in model"), "got {}", msg);
    }

    // Partial-failure behavior: earlier terms in the arg list are mutated
    // before the unknown term is hit. This is documented behavior — if it
    // ever changes to all-or-nothing, this test should be updated deliberately.
    #[test]
    fn fix_stops_at_first_unknown_term_leaving_prior_changes_in_place() {
        let mut s = session_with_terms(&["IH", "ID"]);
        let err = Fix.execute(&mut s, &["IH", "NOPE", "ID"]).unwrap_err();
        assert!(parse_message(&err).contains("NOPE"));
        assert_eq!(s.model.fixed_flags(), &[true, false]);
    }

    #[test]
    fn unfix_stops_at_first_unknown_term_leaving_prior_changes_in_place() {
        let mut s = session_with_terms(&["IH", "ID"]);
        s.model.fix_all();
        let err = Unfix.execute(&mut s, &["IH", "NOPE", "ID"]).unwrap_err();
        assert!(parse_message(&err).contains("NOPE"));
        assert_eq!(s.model.fixed_flags(), &[false, true]);
    }

    // "ALL" is only treated as the keyword when it's the *first* arg —
    // everywhere else it's just a term name, which won't exist in the model.
    #[test]
    fn all_only_special_in_first_position() {
        let mut s = session_with_terms(&["IH"]);
        let err = Fix.execute(&mut s, &["IH", "ALL"]).unwrap_err();
        assert!(parse_message(&err).contains("ALL"));
        assert_eq!(s.model.fixed_flags(), &[true]);
    }

    #[test]
    fn fix_then_unfix_round_trip() {
        let mut s = session_with_terms(&["IH", "ID"]);
        Fix.execute(&mut s, &["IH", "ID"]).unwrap();
        assert_eq!(s.model.fixed_flags(), &[true, true]);
        Unfix.execute(&mut s, &["IH", "ID"]).unwrap();
        assert_eq!(s.model.fixed_flags(), &[false, false]);
    }
}
