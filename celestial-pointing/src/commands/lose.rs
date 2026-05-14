use super::{Command, CommandOutput};
use crate::error::Result;
use crate::session::Session;

pub struct Lose;

impl Command for Lose {
    fn name(&self) -> &str {
        "LOSE"
    }
    fn description(&self) -> &str {
        "Remove term(s) from model"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() || (args.len() == 1 && args[0].eq_ignore_ascii_case("ALL")) {
            session.model.remove_all();
            return Ok(CommandOutput::Text("All terms removed".into()));
        }
        for name in args {
            session.model.remove_term(&name.to_uppercase());
        }
        Ok(CommandOutput::Text(format!(
            "Removed: {}",
            args.join(", ").to_uppercase()
        )))
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

    #[test]
    fn metadata() {
        assert_eq!(Lose.name(), "LOSE");
        assert_eq!(Lose.description(), "Remove term(s) from model");
    }

    #[test]
    fn no_args_removes_everything() {
        let mut s = session_with(&["IH", "ID", "CH"]);
        let body = text(Lose.execute(&mut s, &[]).unwrap());
        assert_eq!(body, "All terms removed");
        assert_eq!(s.model.term_count(), 0);
    }

    #[test]
    fn all_keyword_removes_everything_case_insensitive() {
        for token in ["ALL", "all", "All", "aLL"] {
            let mut s = session_with(&["IH", "ID"]);
            let body = text(Lose.execute(&mut s, &[token]).unwrap());
            assert_eq!(body, "All terms removed", "token {:?}", token);
            assert_eq!(s.model.term_count(), 0);
        }
    }

    #[test]
    fn all_keyword_on_empty_model_succeeds() {
        let mut s = Session::new();
        let body = text(Lose.execute(&mut s, &["ALL"]).unwrap());
        assert_eq!(body, "All terms removed");
    }

    #[test]
    fn removes_named_term() {
        let mut s = session_with(&["IH", "ID", "CH"]);
        let body = text(Lose.execute(&mut s, &["ID"]).unwrap());
        assert_eq!(body, "Removed: ID");
        assert_eq!(s.model.term_names(), vec!["IH", "CH"]);
    }

    #[test]
    fn removes_multiple_terms_in_one_call() {
        let mut s = session_with(&["IH", "ID", "CH", "NP"]);
        let body = text(Lose.execute(&mut s, &["IH", "CH"]).unwrap());
        assert_eq!(body, "Removed: IH, CH");
        assert_eq!(s.model.term_names(), vec!["ID", "NP"]);
    }

    #[test]
    fn lowercase_input_is_uppercased() {
        let mut s = session_with(&["IH", "ID"]);
        let body = text(Lose.execute(&mut s, &["ih"]).unwrap());
        assert_eq!(body, "Removed: IH");
        assert_eq!(s.model.term_names(), vec!["ID"]);
    }

    // Removing a term that isn't in the model is a silent no-op — model.remove_term
    // returns nothing, and LOSE still reports "Removed: <name>". This is intentional
    // (matches TPOINT's lenient behavior); pinned here so it can't change silently.
    #[test]
    fn removing_unknown_term_is_silent_noop() {
        let mut s = session_with(&["IH"]);
        let body = text(Lose.execute(&mut s, &["BOGUS"]).unwrap());
        assert_eq!(body, "Removed: BOGUS");
        assert_eq!(s.model.term_names(), vec!["IH"]);
    }

    // The `args.len() == 1 && ALL` guard means "ALL" mixed with other args goes
    // through the per-term removal path — and since no term is named "ALL", it's
    // a no-op for that arg. Documenting the quirk.
    #[test]
    fn all_alongside_other_args_is_not_special() {
        let mut s = session_with(&["IH", "ID"]);
        let body = text(Lose.execute(&mut s, &["IH", "ALL"]).unwrap());
        assert_eq!(body, "Removed: IH, ALL");
        assert_eq!(s.model.term_names(), vec!["ID"]);
    }
}
