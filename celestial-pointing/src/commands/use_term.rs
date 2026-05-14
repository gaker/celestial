use super::{Command, CommandOutput};
use crate::error::Result;
use crate::session::Session;

pub struct Use;

impl Command for Use {
    fn name(&self) -> &str {
        "USE"
    }
    
    fn description(&self) -> &str {
        "Add term(s) to model"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(crate::error::Error::Parse(
                "USE requires term name(s)".into(),
            ));
        }
        let mut added = Vec::new();
        for name in args {
            session.model.add_term(name)?;
            added.push(name.to_uppercase());
        }
        Ok(CommandOutput::Text(format!("Added: {}", added.join(", "))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Error;

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
        assert_eq!(Use.name(), "USE");
        assert_eq!(Use.description(), "Add term(s) to model");
    }

    #[test]
    fn no_args_errors() {
        let mut s = Session::new();
        let err = Use.execute(&mut s, &[]).unwrap_err();
        assert!(parse_msg(&err).contains("USE requires term name(s)"));
        assert_eq!(s.model.term_count(), 0);
    }

    #[test]
    fn adds_single_term() {
        let mut s = Session::new();
        let body = text(Use.execute(&mut s, &["IH"]).unwrap());
        assert_eq!(body, "Added: IH");
        assert_eq!(s.model.term_names(), vec!["IH"]);
    }

    #[test]
    fn adds_multiple_terms_in_order() {
        let mut s = Session::new();
        let body = text(Use.execute(&mut s, &["IH", "ID", "CH"]).unwrap());
        assert_eq!(body, "Added: IH, ID, CH");
        assert_eq!(s.model.term_names(), vec!["IH", "ID", "CH"]);
    }

    #[test]
    fn lowercase_args_are_uppercased_in_output_and_storage() {
        let mut s = Session::new();
        let body = text(Use.execute(&mut s, &["ih", "Id"]).unwrap());
        assert_eq!(body, "Added: IH, ID");
        assert_eq!(s.model.term_names(), vec!["IH", "ID"]);
    }

    #[test]
    fn unknown_term_errors_with_unknown_term_variant() {
        let mut s = Session::new();
        let err = Use.execute(&mut s, &["NOSUCH"]).unwrap_err();
        assert!(
            matches!(err, Error::UnknownTerm(_)),
            "expected UnknownTerm, got {:?}",
            err,
        );
    }

    // Partial-failure: earlier valid terms are added before the bad one is hit.
    // Documenting so any future tightening to all-or-nothing is intentional.
    #[test]
    fn unknown_term_mid_list_leaves_prior_additions_in_place() {
        let mut s = Session::new();
        let _ = Use.execute(&mut s, &["IH", "BOGUS", "ID"]).unwrap_err();
        assert_eq!(s.model.term_names(), vec!["IH"]);
    }

    // Adding a term that's already in the model is a silent no-op. The success
    // message still lists the duplicate (because USE reports what was *asked
    // for*, not what was newly inserted), but the model is unchanged.
    #[test]
    fn duplicate_term_is_silent_noop() {
        let mut s = Session::new();
        Use.execute(&mut s, &["IH", "ID"]).unwrap();
        let body = text(Use.execute(&mut s, &["IH"]).unwrap());
        assert_eq!(body, "Added: IH");
        assert_eq!(s.model.term_names(), vec!["IH", "ID"]);
        assert_eq!(s.model.term_count(), 2);
    }

    #[test]
    fn duplicate_within_single_call_is_silent_noop() {
        let mut s = Session::new();
        let body = text(Use.execute(&mut s, &["IH", "IH", "ID"]).unwrap());
        assert_eq!(body, "Added: IH, IH, ID");
        assert_eq!(s.model.term_names(), vec!["IH", "ID"]);
    }

    #[test]
    fn duplicate_check_is_case_insensitive() {
        let mut s = Session::new();
        Use.execute(&mut s, &["IH"]).unwrap();
        Use.execute(&mut s, &["ih"]).unwrap();
        Use.execute(&mut s, &["Ih"]).unwrap();
        assert_eq!(s.model.term_count(), 1);
    }

    #[test]
    fn duplicate_preserves_existing_coefficient() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.model.set_coefficients(&[42.0]).unwrap();
        Use.execute(&mut s, &["IH"]).unwrap();
        assert_eq!(s.model.coefficients(), &[42.0]);
    }

    #[test]
    fn parametric_term_works() {
        let mut s = Session::new();
        let body = text(Use.execute(&mut s, &["TX3"]).unwrap());
        assert_eq!(body, "Added: TX3");
        assert_eq!(s.model.term_names(), vec!["TX3"]);
    }
}
