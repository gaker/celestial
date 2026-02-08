use celestial_core::Angle;
use crate::error::{Error, Result};
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Lst;

impl Command for Lst {
    fn name(&self) -> &str { "LST" }
    fn description(&self) -> &str { "Set or show local sidereal time" }

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
    let hh = h.floor() as u32;
    let mm = ((h - hh as f64) * 60.0).floor() as u32;
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
    let hours: f64 = s.parse().map_err(|_| {
        Error::Parse(format!("invalid LST value: {}", s))
    })?;
    validate_hours(hours)?;
    Ok(Angle::from_hours(hours))
}

fn parse_hms(h: &str, m: &str, s: &str) -> Result<Angle> {
    let hh: f64 = h.parse().map_err(|_| Error::Parse(format!("invalid hours: {}", h)))?;
    let mm: f64 = m.parse().map_err(|_| Error::Parse(format!("invalid minutes: {}", m)))?;
    let ss: f64 = s.parse().map_err(|_| Error::Parse(format!("invalid seconds: {}", s)))?;
    let hours = hh + mm / 60.0 + ss / 3600.0;
    validate_hours(hours)?;
    Ok(Angle::from_hours(hours))
}

fn validate_hours(hours: f64) -> Result<()> {
    if !(0.0..24.0).contains(&hours) {
        return Err(Error::Parse(
            format!("LST must be in range [0, 24), got {}", hours),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn show_when_no_lst_set() {
        let mut session = Session::new();
        let result = Lst.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => assert_eq!(s, "No LST set"),
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn set_decimal_hours() {
        let mut session = Session::new();
        Lst.execute(&mut session, &["14.5"]).unwrap();
        let lst = session.current_lst().unwrap();
        assert_eq!(lst.hours(), 14.5);
    }

    #[test]
    fn set_hms() {
        let mut session = Session::new();
        Lst.execute(&mut session, &["14", "30", "00"]).unwrap();
        let lst = session.current_lst().unwrap();
        assert_eq!(lst.hours(), 14.5);
    }

    #[test]
    fn show_after_set() {
        let mut session = Session::new();
        Lst.execute(&mut session, &["14", "30", "00"]).unwrap();
        let result = Lst.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => assert!(s.starts_with("LST = 14h 30m")),
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn clear_override() {
        let mut session = Session::new();
        Lst.execute(&mut session, &["14.5"]).unwrap();
        Lst.execute(&mut session, &["CLEAR"]).unwrap();
        assert!(session.lst_override.is_none());
        assert!(session.current_lst().is_err());
    }

    #[test]
    fn clear_case_insensitive() {
        let mut session = Session::new();
        Lst.execute(&mut session, &["14.5"]).unwrap();
        Lst.execute(&mut session, &["clear"]).unwrap();
        assert!(session.lst_override.is_none());
    }

    #[test]
    fn reject_out_of_range() {
        let mut session = Session::new();
        assert!(Lst.execute(&mut session, &["25.0"]).is_err());
        assert!(Lst.execute(&mut session, &["-1.0"]).is_err());
    }

    #[test]
    fn reject_invalid_input() {
        let mut session = Session::new();
        assert!(Lst.execute(&mut session, &["abc"]).is_err());
    }

    #[test]
    fn reject_wrong_arg_count() {
        let mut session = Session::new();
        assert!(Lst.execute(&mut session, &["14", "30"]).is_err());
    }
}
