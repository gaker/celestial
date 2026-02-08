use celestial_core::Angle;
use crate::error::Result;
use crate::observation::PierSide;
use crate::parser::parse_coordinates;
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Correct;

impl Command for Correct {
    fn name(&self) -> &str { "CORRECT" }
    fn description(&self) -> &str { "Compute actual sky position from encoder reading" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        let (enc_ra, enc_dec) = parse_coordinates(args)?;
        let lst = session.current_lst()?;
        let lat = Angle::from_radians(session.latitude());
        let ha = lst - enc_ra;
        let pier = pier_from_ha(ha);
        let (true_ra, true_dec) = session.model.command_to_target(enc_ra, enc_dec, lst, lat, pier);
        let delta_ra = (true_ra - enc_ra).wrapped();
        let delta_dec = (true_dec - enc_dec).wrapped();
        Ok(CommandOutput::Text(format_result(enc_ra, enc_dec, true_ra, true_dec, delta_ra, delta_dec)))
    }
}

fn pier_from_ha(ha: Angle) -> PierSide {
    if ha.radians() >= 0.0 { PierSide::East } else { PierSide::West }
}

fn format_result(
    enc_ra: Angle,
    enc_dec: Angle,
    true_ra: Angle,
    true_dec: Angle,
    dra: Angle,
    ddec: Angle,
) -> String {
    format!(
        "Encoder:  {}  {}\nActual:   {}  {}\n  \u{0394}RA:  {:+.2}s\n  \u{0394}Dec: {:+.1}\"",
        format_ra(enc_ra),
        format_dec(enc_dec),
        format_ra(true_ra),
        format_dec(true_dec),
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
    fn empty_model_encoder_equals_actual() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        let args = vec!["12", "30", "00", "+45", "00", "00"];
        let result = Correct.execute(&mut session, &args).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(s.contains("Encoder:"));
                assert!(s.contains("Actual:"));
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
    fn correct_requires_lst() {
        let mut session = Session::new();
        let args = vec!["12.5", "45.0"];
        let result = Correct.execute(&mut session, &args);
        assert!(result.is_err());
    }

    #[test]
    fn correct_with_model_produces_nonzero_deltas() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        session.model.add_term("IH").unwrap();
        session.model.set_coefficients(&[30.0]).unwrap();
        let args = vec!["12.5", "45.0"];
        let result = Correct.execute(&mut session, &args).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(!s.contains("\u{0394}RA:  +0.00s"));
            }
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn correct_decimal_args() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        let args = vec!["12.5", "45.0"];
        let result = Correct.execute(&mut session, &args);
        assert!(result.is_ok());
    }
}
