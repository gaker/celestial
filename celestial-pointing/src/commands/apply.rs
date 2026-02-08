use celestial_core::Angle;
use crate::error::Result;
use crate::observation::PierSide;
use crate::parser::parse_coordinates;
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Apply;

impl Command for Apply {
    fn name(&self) -> &str { "APPLY" }
    fn description(&self) -> &str { "Compute commanded position for target" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        let (ra, dec) = parse_coordinates(args)?;
        let lst = session.current_lst()?;
        let lat = Angle::from_radians(session.latitude());
        let ha = lst - ra;
        let pier = pier_from_ha(ha);
        let (cmd_ra, cmd_dec) = session.model.target_to_command(ra, dec, lst, lat, pier);
        let delta_ra = (cmd_ra - ra).wrapped();
        let delta_dec = (cmd_dec - dec).wrapped();
        Ok(CommandOutput::Text(format_result(ra, dec, cmd_ra, cmd_dec, delta_ra, delta_dec)))
    }
}

fn pier_from_ha(ha: Angle) -> PierSide {
    if ha.radians() >= 0.0 { PierSide::East } else { PierSide::West }
}

fn format_result(
    ra: Angle,
    dec: Angle,
    cmd_ra: Angle,
    cmd_dec: Angle,
    dra: Angle,
    ddec: Angle,
) -> String {
    format!(
        "Target:   {}  {}\nCommand:  {}  {}\n  \u{0394}RA:  {:+.2}s\n  \u{0394}Dec: {:+.1}\"",
        format_ra(ra),
        format_dec(dec),
        format_ra(cmd_ra),
        format_dec(cmd_dec),
        dra.arcseconds() / 15.0,
        ddec.arcseconds(),
    )
}

fn format_ra(a: Angle) -> String {
    let total_h = a.normalized().hours();
    let h = total_h.floor() as u32;
    let rem = (total_h - h as f64) * 60.0;
    let m = rem.floor() as u32;
    let s = (rem - m as f64) * 60.0;
    format!("{:02}h {:02}m {:05.2}s", h, m, s)
}

fn format_dec(a: Angle) -> String {
    let deg = a.degrees();
    let sign = if deg < 0.0 { '-' } else { '+' };
    let abs = deg.abs();
    let d = abs.floor() as u32;
    let rem = (abs - d as f64) * 60.0;
    let m = rem.floor() as u32;
    let s = (rem - m as f64) * 60.0;
    format!("{}{:02}\u{00b0} {:02}' {:04.1}\"", sign, d, m, s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::session::Session;

    #[test]
    fn empty_model_returns_target_equals_command() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        let args = vec!["12", "30", "00", "+45", "00", "00"];
        let result = Apply.execute(&mut session, &args).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(s.contains("Target:"));
                assert!(s.contains("Command:"));
                assert!(s.contains("\u{0394}RA:  +0.00s"));
                assert!(s.contains("\u{0394}Dec: +0.0\""));
            }
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn pier_east_when_ha_positive() {
        let ha = Angle::from_hours(2.0);
        assert_eq!(pier_from_ha(ha), PierSide::East);
    }

    #[test]
    fn pier_west_when_ha_negative() {
        let ha = Angle::from_hours(-2.0);
        assert_eq!(pier_from_ha(ha), PierSide::West);
    }

    #[test]
    fn pier_east_when_ha_zero() {
        let ha = Angle::from_hours(0.0);
        assert_eq!(pier_from_ha(ha), PierSide::East);
    }

    #[test]
    fn apply_requires_lst() {
        let mut session = Session::new();
        let args = vec!["12.5", "45.0"];
        let result = Apply.execute(&mut session, &args);
        assert!(result.is_err());
    }

    #[test]
    fn apply_with_model_produces_nonzero_deltas() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        session.model.add_term("IH").unwrap();
        session.model.set_coefficients(&[30.0]).unwrap();
        let args = vec!["12.5", "45.0"];
        let result = Apply.execute(&mut session, &args).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(!s.contains("\u{0394}RA:  +0.00s"));
            }
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn apply_decimal_args() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        let args = vec!["12.5", "45.0"];
        let result = Apply.execute(&mut session, &args);
        assert!(result.is_ok());
    }

    #[test]
    fn format_ra_zero() {
        let s = format_ra(Angle::from_hours(0.0));
        assert_eq!(s, "00h 00m 00.00s");
    }

    #[test]
    fn format_ra_12h() {
        let s = format_ra(Angle::from_hours(12.5));
        assert_eq!(s, "12h 30m 00.00s");
    }

    #[test]
    fn format_dec_positive() {
        let s = format_dec(Angle::from_degrees(45.5));
        assert_eq!(s, "+45\u{00b0} 30' 00.0\"");
    }

    #[test]
    fn format_dec_negative() {
        let s = format_dec(Angle::from_degrees(-30.25));
        assert_eq!(s, "-30\u{00b0} 15' 00.0\"");
    }
}
