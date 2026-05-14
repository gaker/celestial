use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;

pub struct Parallel;
pub struct Chain;

impl Command for Parallel {
    fn name(&self) -> &str {
        "PARAL"
    }

    fn description(&self) -> &str {
        "Apply terms in parallel"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("PARAL requires term names or ALL".into()));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            session.model.set_all_parallel();
            return Ok(CommandOutput::Text(format!(
                "All {} terms set to parallel",
                session.model.term_count()
            )));
        }
        let mut set = Vec::new();
        for name in args {
            let upper = name.to_uppercase();
            if session.model.set_parallel(&upper) {
                set.push(upper);
            } else {
                return Err(Error::Parse(format!("term {} not in model", name)));
            }
        }
        Ok(CommandOutput::Text(format!("Parallel: {}", set.join(" "))))
    }
}

impl Command for Chain {
    fn name(&self) -> &str {
        "CHAIN"
    }
    
    fn description(&self) -> &str {
        "Apply terms sequentially (chained)"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("CHAIN requires term names or ALL".into()));
        }
        if args[0].eq_ignore_ascii_case("ALL") {
            session.model.set_all_chained();
            return Ok(CommandOutput::Text(format!(
                "All {} terms set to chained",
                session.model.term_count()
            )));
        }
        let mut set = Vec::new();
        for name in args {
            let upper = name.to_uppercase();
            if session.model.set_chained(&upper) {
                set.push(upper);
            } else {
                return Err(Error::Parse(format!("term {} not in model", name)));
            }
        }
        Ok(CommandOutput::Text(format!("Chained: {}", set.join(" "))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn session_with(terms: &[&str]) -> Session {
        let mut s = Session::new();
        for t in terms {
            s.model.add_term(t).unwrap();
        }
        s
    }

    fn text(out: CommandOutput) -> String {
        match out {
            CommandOutput::Text(s) => s,
            other => panic!("expected Text, got {:?}", other),
        }
    }

    fn parse_msg(err: &Error) -> &str {
        match err {
            Error::Parse(m) => m,
            other => panic!("expected Parse, got {:?}", other),
        }
    }

    #[test]
    fn metadata() {
        assert_eq!(Parallel.name(), "PARAL");
        assert_eq!(Parallel.description(), "Apply terms in parallel");
        assert_eq!(Chain.name(), "CHAIN");
        assert_eq!(Chain.description(), "Apply terms sequentially (chained)");
    }

    #[test]
    fn paral_no_args_errors() {
        let mut s = session_with(&["IH"]);
        let err = Parallel.execute(&mut s, &[]).unwrap_err();
        assert!(parse_msg(&err).contains("PARAL requires"));
        // add_term default is chained (parallel = false). Confirm no mutation.
        assert_eq!(s.model.parallel_flags(), &[false]);
    }

    #[test]
    fn chain_no_args_errors() {
        let mut s = session_with(&["IH"]);
        let err = Chain.execute(&mut s, &[]).unwrap_err();
        assert!(parse_msg(&err).contains("CHAIN requires"));
    }

    #[test]
    fn paral_all_sets_every_term_parallel_case_insensitive() {
        for token in ["ALL", "all", "All"] {
            let mut s = session_with(&["IH", "ID", "CH"]);
            // Put them in chained state first so the assertion proves PARAL ALL flipped them.
            s.model.set_all_chained();
            assert_eq!(s.model.parallel_flags(), &[false, false, false]);
            let body = text(Parallel.execute(&mut s, &[token]).unwrap());
            assert_eq!(body, "All 3 terms set to parallel", "token {:?}", token);
            assert_eq!(s.model.parallel_flags(), &[true, true, true]);
        }
    }

    #[test]
    fn chain_all_sets_every_term_chained_case_insensitive() {
        for token in ["ALL", "all"] {
            let mut s = session_with(&["IH", "ID"]);
            let body = text(Chain.execute(&mut s, &[token]).unwrap());
            assert_eq!(body, "All 2 terms set to chained", "token {:?}", token);
            assert_eq!(s.model.parallel_flags(), &[false, false]);
        }
    }

    #[test]
    fn paral_all_on_empty_model_reports_zero() {
        let mut s = Session::new();
        let body = text(Parallel.execute(&mut s, &["ALL"]).unwrap());
        assert_eq!(body, "All 0 terms set to parallel");
    }

    #[test]
    fn chain_all_on_empty_model_reports_zero() {
        let mut s = Session::new();
        let body = text(Chain.execute(&mut s, &["ALL"]).unwrap());
        assert_eq!(body, "All 0 terms set to chained");
    }

    #[test]
    fn paral_named_terms_uppercases_and_flips_only_those() {
        let mut s = session_with(&["IH", "ID", "CH"]);
        s.model.set_all_chained();
        let body = text(Parallel.execute(&mut s, &["ih", "ch"]).unwrap());
        assert_eq!(body, "Parallel: IH CH");
        assert_eq!(s.model.parallel_flags(), &[true, false, true]);
    }

    #[test]
    fn chain_named_terms_uppercases_and_flips_only_those() {
        let mut s = session_with(&["IH", "ID", "CH"]);
        // Set them all parallel first so CHAIN actually flips them.
        s.model.set_all_parallel();
        let body = text(Chain.execute(&mut s, &["ih", "ID"]).unwrap());
        assert_eq!(body, "Chained: IH ID");
        assert_eq!(s.model.parallel_flags(), &[false, false, true]);
    }

    #[test]
    fn paral_unknown_term_errors_with_original_casing() {
        let mut s = session_with(&["IH"]);
        let err = Parallel.execute(&mut s, &["bogus"]).unwrap_err();
        assert!(parse_msg(&err).contains("term bogus not in model"));
    }

    #[test]
    fn chain_unknown_term_errors_with_original_casing() {
        let mut s = session_with(&["IH"]);
        let err = Chain.execute(&mut s, &["bogus"]).unwrap_err();
        assert!(parse_msg(&err).contains("term bogus not in model"));
    }

    // Like FIX/UNFIX, mid-list failure leaves earlier changes applied.
    #[test]
    fn paral_partial_failure_leaves_prior_changes() {
        let mut s = session_with(&["IH", "ID"]);
        s.model.set_all_chained();
        let _ = Parallel.execute(&mut s, &["IH", "NOPE", "ID"]).unwrap_err();
        assert_eq!(s.model.parallel_flags(), &[true, false]);
    }

    #[test]
    fn chain_partial_failure_leaves_prior_changes() {
        let mut s = session_with(&["IH", "ID"]);
        s.model.set_all_parallel();
        let _ = Chain.execute(&mut s, &["IH", "NOPE", "ID"]).unwrap_err();
        assert_eq!(s.model.parallel_flags(), &[false, true]);
    }

    // "ALL" only special as first arg; in any other position it's a term lookup.
    #[test]
    fn paral_all_only_special_in_first_position() {
        let mut s = session_with(&["IH"]);
        s.model.set_all_chained();
        let err = Parallel.execute(&mut s, &["IH", "ALL"]).unwrap_err();
        assert!(parse_msg(&err).contains("ALL"));
        assert_eq!(s.model.parallel_flags(), &[true]);
    }

    #[test]
    fn chain_all_only_special_in_first_position() {
        let mut s = session_with(&["IH"]);
        let err = Chain.execute(&mut s, &["IH", "ALL"]).unwrap_err();
        assert!(parse_msg(&err).contains("ALL"));
        assert_eq!(s.model.parallel_flags(), &[false]);
    }

    #[test]
    fn paral_then_chain_round_trip() {
        let mut s = session_with(&["IH", "ID"]);
        Chain.execute(&mut s, &["ALL"]).unwrap();
        assert_eq!(s.model.parallel_flags(), &[false, false]);
        Parallel.execute(&mut s, &["IH", "ID"]).unwrap();
        assert_eq!(s.model.parallel_flags(), &[true, true]);
    }
}
