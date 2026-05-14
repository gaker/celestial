use super::{Command, CommandOutput};
use crate::error::Result;
use crate::session::Session;

pub struct Outmod;

impl Command for Outmod {
    fn name(&self) -> &str {
        "OUTMOD"
    }
    fn description(&self) -> &str {
        "Save model to file"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(crate::error::Error::Parse(
                "OUTMOD requires a filename".into(),
            ));
        }
        let mut output = String::new();
        for (name, &coeff) in session
            .model
            .term_names()
            .iter()
            .zip(session.model.coefficients().iter())
        {
            output += &format!("{} {:.6}\n", name, coeff);
        }
        output += "END\n";
        std::fs::write(args[0], &output).map_err(crate::error::Error::Io)?;
        Ok(CommandOutput::Text(format!("Model saved to {}", args[0])))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::inmod::Inmod;
    use crate::error::Error;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static SEQ: AtomicU64 = AtomicU64::new(0);

    fn tmp_path(tag: &str) -> PathBuf {
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let mut p = std::env::temp_dir();
        p.push(format!("celpoint-outmod-{}-{}-{}.mod", tag, std::process::id(), n));
        p
    }

    fn cleanup(p: &PathBuf) {
        let _ = std::fs::remove_file(p);
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
        assert_eq!(Outmod.name(), "OUTMOD");
        assert_eq!(Outmod.description(), "Save model to file");
    }

    #[test]
    fn no_args_errors() {
        let mut s = Session::new();
        let err = Outmod.execute(&mut s, &[]).unwrap_err();
        assert!(parse_msg(&err).contains("OUTMOD requires a filename"));
    }

    #[test]
    fn write_to_unwritable_path_returns_io_error() {
        let mut s = Session::new();
        let err = Outmod
            .execute(&mut s, &["/no/such/directory/path/file.mod"])
            .unwrap_err();
        assert!(matches!(err, Error::Io(_)), "got {:?}", err);
    }

    #[test]
    fn empty_model_writes_only_end_marker() {
        let mut s = Session::new();
        let p = tmp_path("empty");
        let body = text(Outmod.execute(&mut s, &[p.to_str().unwrap()]).unwrap());
        let contents = std::fs::read_to_string(&p).unwrap();
        cleanup(&p);
        assert_eq!(contents, "END\n");
        assert!(body.starts_with("Model saved to "));
        assert!(body.contains(p.to_str().unwrap()));
    }

    #[test]
    fn writes_terms_with_six_decimal_coefficients() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.model.add_term("ID").unwrap();
        s.model.set_coefficients(&[10.5, -20.25]).unwrap();

        let p = tmp_path("basic");
        Outmod.execute(&mut s, &[p.to_str().unwrap()]).unwrap();
        let contents = std::fs::read_to_string(&p).unwrap();
        cleanup(&p);
        assert_eq!(contents, "IH 10.500000\nID -20.250000\nEND\n");
    }

    #[test]
    fn coefficient_formatting_pads_to_six_decimals() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.model.set_coefficients(&[1.0]).unwrap();
        let p = tmp_path("pad");
        Outmod.execute(&mut s, &[p.to_str().unwrap()]).unwrap();
        let contents = std::fs::read_to_string(&p).unwrap();
        cleanup(&p);
        assert!(contents.contains("IH 1.000000\n"), "got {:?}", contents);
    }

    #[test]
    fn zero_coefficient_is_emitted_verbatim() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        // Default coefficient is 0.0; do not call set_coefficients.
        let p = tmp_path("zero");
        Outmod.execute(&mut s, &[p.to_str().unwrap()]).unwrap();
        let contents = std::fs::read_to_string(&p).unwrap();
        cleanup(&p);
        assert!(contents.contains("IH 0.000000\n"));
    }

    #[test]
    fn overwrites_existing_file() {
        let p = tmp_path("overwrite");
        std::fs::write(&p, "STALE CONTENT THAT SHOULD BE REPLACED").unwrap();

        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.model.set_coefficients(&[1.0]).unwrap();
        Outmod.execute(&mut s, &[p.to_str().unwrap()]).unwrap();

        let contents = std::fs::read_to_string(&p).unwrap();
        cleanup(&p);
        assert!(!contents.contains("STALE"));
        assert_eq!(contents, "IH 1.000000\nEND\n");
    }

    // Round-trip: a file written by OUTMOD must load cleanly via INMOD and
    // restore identical terms + coefficients. This pins the file format
    // contract between the two commands.
    #[test]
    fn round_trip_with_inmod_restores_model() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.model.add_term("ID").unwrap();
        s.model.add_term("CH").unwrap();
        s.model.set_coefficients(&[12.5, -0.125, 3.0]).unwrap();

        let p = tmp_path("roundtrip");
        Outmod.execute(&mut s, &[p.to_str().unwrap()]).unwrap();

        let mut s2 = Session::new();
        Inmod.execute(&mut s2, &[p.to_str().unwrap()]).unwrap();
        cleanup(&p);

        assert_eq!(s2.model.term_names(), vec!["IH", "ID", "CH"]);
        assert_eq!(s2.model.coefficients(), &[12.5, -0.125, 3.0]);
    }

    #[test]
    fn extra_args_after_filename_are_ignored() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        let p = tmp_path("extra");
        Outmod
            .execute(&mut s, &[p.to_str().unwrap(), "ignored", "args"])
            .unwrap();
        let contents = std::fs::read_to_string(&p).unwrap();
        cleanup(&p);
        assert_eq!(contents, "IH 0.000000\nEND\n");
    }

    #[test]
    fn success_message_echoes_filename() {
        let mut s = Session::new();
        let p = tmp_path("msg");
        let body = text(Outmod.execute(&mut s, &[p.to_str().unwrap()]).unwrap());
        cleanup(&p);
        assert_eq!(body, format!("Model saved to {}", p.to_str().unwrap()));
    }
}
