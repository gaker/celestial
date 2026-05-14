use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;
use celestial_core::Angle;

pub struct Lst;

impl Command for Lst {
    fn name(&self) -> &str {
        "LST"
    }

    fn description(&self) -> &str {
        "Set or show local sidereal time"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return show_lst(session);
        }
        if args[0].eq_ignore_ascii_case("CLEAR") {
            session.lst_override = None;
            return Ok(CommandOutput::Text("LST override cleared".to_string()));
        }
        let angle = parse_lst_args(args)?;
        session.lst_override = Some(angle);
        Ok(CommandOutput::Text(format_lst(angle)))
    }
}

fn show_lst(session: &Session) -> Result<CommandOutput> {
    match session.current_lst() {
        Ok(lst) => Ok(CommandOutput::Text(format_lst(lst))),
        Err(_) => Ok(CommandOutput::Text("No LST set".to_string())),
    }
}

fn format_lst(lst: Angle) -> String {
    let h = lst.hours();
    let hh = libm::floor(h) as u32;
    let mm = libm::floor((h - hh as f64) * 60.0) as u32;
    let ss = (h - hh as f64) * 3600.0 - mm as f64 * 60.0;
    format!("LST = {:02}h {:02}m {:06.3}s", hh, mm, ss)
}

fn parse_lst_args(args: &[&str]) -> Result<Angle> {
    match args.len() {
        1 => parse_decimal_hours(args[0]),
        3 => parse_hms(args[0], args[1], args[2]),
        _ => Err(Error::Parse(
            "LST expects decimal hours (e.g. 14.5) or h m s (e.g. 14 30 00)".to_string(),
        )),
    }
}

fn parse_decimal_hours(s: &str) -> Result<Angle> {
    let hours: f64 = s
        .parse()
        .map_err(|_| Error::Parse(format!("invalid LST value: {}", s)))?;
    validate_hours(hours)?;
    Ok(Angle::from_hours(hours))
}

fn parse_hms(h: &str, m: &str, s: &str) -> Result<Angle> {
    let hh: f64 = h
        .parse()
        .map_err(|_| Error::Parse(format!("invalid hours: {}", h)))?;
    let mm: f64 = m
        .parse()
        .map_err(|_| Error::Parse(format!("invalid minutes: {}", m)))?;
    let ss: f64 = s
        .parse()
        .map_err(|_| Error::Parse(format!("invalid seconds: {}", s)))?;
    let hours = hh + mm / 60.0 + ss / 3600.0;
    validate_hours(hours)?;
    Ok(Angle::from_hours(hours))
}

fn validate_hours(hours: f64) -> Result<()> {
    if !(0.0..24.0).contains(&hours) {
        return Err(Error::Parse(format!(
            "LST must be in range [0, 24), got {}",
            hours
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(session: &mut Session, args: &[&str]) -> Result<String> {
        Lst.execute(session, args).map(|out| match out {
            CommandOutput::Text(s) => s,
            other => panic!("expected Text, got {:?}", other),
        })
    }

    fn parse_msg(err: &Error) -> &str {
        match err {
            Error::Parse(m) => m,
            other => panic!("expected Parse, got {:?}", other),
        }
    }

    #[test]
    fn metadata() {
        assert_eq!(Lst.name(), "LST");
        assert_eq!(Lst.description(), "Set or show local sidereal time");
    }

    #[test]
    fn no_args_with_no_override_reports_no_lst() {
        let mut s = Session::new();
        assert_eq!(run(&mut s, &[]).unwrap(), "No LST set");
    }

    #[test]
    fn no_args_with_override_shows_current_value() {
        let mut s = Session::new();
        s.lst_override = Some(Angle::from_hours(14.5));
        let body = run(&mut s, &[]).unwrap();
        assert_eq!(body, "LST = 14h 30m 00.000s");
    }

    #[test]
    fn clear_resets_override_case_insensitive() {
        for token in ["CLEAR", "clear", "Clear"] {
            let mut s = Session::new();
            s.lst_override = Some(Angle::from_hours(10.0));
            let body = run(&mut s, &[token]).unwrap();
            assert_eq!(body, "LST override cleared");
            assert!(s.lst_override.is_none());
        }
    }

    #[test]
    fn clear_with_no_prior_override_still_succeeds() {
        let mut s = Session::new();
        let body = run(&mut s, &["CLEAR"]).unwrap();
        assert_eq!(body, "LST override cleared");
        assert!(s.lst_override.is_none());
    }

    #[test]
    fn decimal_hours_sets_override() {
        let mut s = Session::new();
        let body = run(&mut s, &["14.5"]).unwrap();
        assert_eq!(body, "LST = 14h 30m 00.000s");
        let lst = s.lst_override.unwrap();
        assert!((lst.hours() - 14.5).abs() < 1e-12);
    }

    #[test]
    fn hms_sets_override() {
        let mut s = Session::new();
        let body = run(&mut s, &["14", "30", "15.5"]).unwrap();
        assert_eq!(body, "LST = 14h 30m 15.500s");
    }

    #[test]
    fn zero_is_a_valid_lst() {
        let mut s = Session::new();
        let body = run(&mut s, &["0"]).unwrap();
        assert_eq!(body, "LST = 00h 00m 00.000s");
    }

    #[test]
    fn exactly_24_is_rejected() {
        let mut s = Session::new();
        let err = run(&mut s, &["24"]).unwrap_err();
        assert!(parse_msg(&err).contains("LST must be in range"));
    }

    #[test]
    fn negative_is_rejected() {
        let mut s = Session::new();
        let err = run(&mut s, &["-0.5"]).unwrap_err();
        assert!(parse_msg(&err).contains("LST must be in range"));
    }

    #[test]
    fn hms_above_24_rejected() {
        let mut s = Session::new();
        let err = run(&mut s, &["23", "59", "61"]).unwrap_err();
        assert!(parse_msg(&err).contains("LST must be in range"));
    }

    #[test]
    fn unparseable_decimal_errors() {
        let mut s = Session::new();
        let err = run(&mut s, &["banana"]).unwrap_err();
        assert!(parse_msg(&err).contains("invalid LST value"));
    }

    #[test]
    fn unparseable_hms_pieces_error_with_specific_field() {
        let mut s = Session::new();
        let err = run(&mut s, &["bad", "30", "0"]).unwrap_err();
        assert!(parse_msg(&err).contains("invalid hours"));

        let err = run(&mut s, &["14", "bad", "0"]).unwrap_err();
        assert!(parse_msg(&err).contains("invalid minutes"));

        let err = run(&mut s, &["14", "30", "bad"]).unwrap_err();
        assert!(parse_msg(&err).contains("invalid seconds"));
    }

    #[test]
    fn wrong_arg_count_errors() {
        let mut s = Session::new();
        for args in [vec!["14", "30"], vec!["14", "30", "0", "extra"]] {
            let err = run(&mut s, &args).unwrap_err();
            assert!(
                parse_msg(&err).contains("decimal hours"),
                "args {:?}: {}",
                args,
                parse_msg(&err),
            );
        }
    }

    #[test]
    fn failed_parse_does_not_overwrite_existing_override() {
        let mut s = Session::new();
        s.lst_override = Some(Angle::from_hours(12.0));
        let _ = run(&mut s, &["banana"]).unwrap_err();
        assert!((s.lst_override.unwrap().hours() - 12.0).abs() < 1e-12);
    }

    #[test]
    fn format_lst_rounds_seconds_consistently() {
        let s1 = format_lst(Angle::from_hours(0.0));
        assert_eq!(s1, "LST = 00h 00m 00.000s");
        let s2 = format_lst(Angle::from_hours(23.999999));
        assert!(s2.starts_with("LST = 23h 59m"), "got {}", s2);
    }
}
