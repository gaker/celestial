use super::{Command, CommandOutput};
use crate::error::Result;
use crate::observation::PierSide;
use crate::parser::parse_coordinates;
use crate::session::Session;
use celestial_core::Angle;

pub struct Predict;

impl Command for Predict {
    fn name(&self) -> &str {
        "PREDICT"
    }
    fn description(&self) -> &str {
        "Show correction breakdown by term"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        let (ra, dec) = parse_coordinates(args)?;
        let lst = session.current_lst()?;
        let lat = session.latitude();
        let ha = lst - ra;
        let pier = pier_from_ha(ha);
        let breakdown =
            session
                .model
                .predict_breakdown(ha.radians(), dec.radians(), lat, pier.sign());
        let (cmd_ra, cmd_dec) =
            session
                .model
                .target_to_command(ra, dec, lst, Angle::from_radians(lat), pier);
        Ok(CommandOutput::Text(format_predict(
            ra, dec, ha, &breakdown, cmd_ra, cmd_dec,
        )))
    }
}

fn pier_from_ha(ha: Angle) -> PierSide {
    if ha.radians() >= 0.0 {
        PierSide::East
    } else {
        PierSide::West
    }
}

fn format_predict(
    ra: Angle,
    dec: Angle,
    ha: Angle,
    breakdown: &[(String, f64, f64)],
    cmd_ra: Angle,
    cmd_dec: Angle,
) -> String {
    let mut lines = Vec::new();
    lines.push(format!("Target: {}  {}", format_ra(ra), format_dec(dec)));
    lines.push(format!("HA: {}  Dec: {}", format_ha(ha), format_dec(dec)));
    lines.push(String::new());
    lines.push(format!(
        "{:<12} {:>10} {:>10}",
        "Term", "\u{0394}HA (\")", "\u{0394}Dec (\")"
    ));
    lines.push("\u{2500}".repeat(34));
    let (total_dh, total_dd) = append_breakdown(&mut lines, breakdown);
    lines.push("\u{2500}".repeat(34));
    lines.push(format!(
        "{:<12} {:>10.2} {:>10.2}",
        "Total", total_dh, total_dd
    ));
    lines.push(String::new());
    lines.push(format!(
        "Command: {}  {}",
        format_ra(cmd_ra),
        format_dec(cmd_dec)
    ));
    lines.join("\n")
}

fn append_breakdown(lines: &mut Vec<String>, breakdown: &[(String, f64, f64)]) -> (f64, f64) {
    let mut total_dh = 0.0;
    let mut total_dd = 0.0;
    for (name, dh, dd) in breakdown {
        lines.push(format!("{:<12} {:>10.2} {:>10.2}", name, dh, dd));
        total_dh += dh;
        total_dd += dd;
    }
    (total_dh, total_dd)
}

fn format_ra(angle: Angle) -> String {
    let h = angle.hours().abs();
    let hh = libm::floor(h) as u32;
    let remainder = (h - hh as f64) * 60.0;
    let mm = libm::floor(remainder) as u32;
    let ss = (remainder - mm as f64) * 60.0;
    format!("{:02}h {:02}m {:05.2}s", hh, mm, ss)
}

fn format_dec(angle: Angle) -> String {
    let deg = angle.degrees();
    let sign = if deg < 0.0 { "-" } else { "+" };
    let total = deg.abs();
    let dd = libm::floor(total) as u32;
    let remainder = (total - dd as f64) * 60.0;
    let mm = libm::floor(remainder) as u32;
    let ss = (remainder - mm as f64) * 60.0;
    format!("{}{:02}\u{00b0} {:02}' {:04.1}\"", sign, dd, mm, ss)
}

fn format_ha(angle: Angle) -> String {
    let h = angle.hours();
    let sign = if h < 0.0 { "-" } else { "+" };
    let total = h.abs();
    let hh = libm::floor(total) as u32;
    let remainder = (total - hh as f64) * 60.0;
    let mm = libm::floor(remainder) as u32;
    let ss = (remainder - mm as f64) * 60.0;
    format!("{}{:02}h {:02}m {:05.2}s", sign, hh, mm, ss)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::session::Session;

    #[test]
    fn empty_model_zero_correction() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        let result = Predict.execute(&mut session, &["12.5", "45.0"]).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(s.contains("Total"));
                assert!(s.contains("0.00"));
                assert!(s.contains("Command:"));
            }
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn single_ih_term() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        session.model.add_term("IH").unwrap();
        session.model.set_coefficients(&[10.0]).unwrap();
        let result = Predict.execute(&mut session, &["12.5", "45.0"]).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(s.contains("IH"));
                assert!(s.contains("-10.00"));
            }
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn total_matches_sum() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        session.model.add_term("IH").unwrap();
        session.model.add_term("ID").unwrap();
        session.model.set_coefficients(&[10.0, 20.0]).unwrap();
        let ha = Angle::from_hours(14.0) - Angle::from_hours(12.5);
        let breakdown = session.model.predict_breakdown(
            ha.radians(),
            Angle::from_degrees(45.0).radians(),
            0.0,
            PierSide::East.sign(),
        );
        let sum_dh: f64 = breakdown.iter().map(|(_, dh, _)| dh).sum();
        let sum_dd: f64 = breakdown.iter().map(|(_, _, dd)| dd).sum();
        let (total_dh, total_dd) = session.model.apply_equatorial(
            ha.radians(),
            Angle::from_degrees(45.0).radians(),
            0.0,
            PierSide::East.sign(),
        );
        assert_eq!(sum_dh, total_dh);
        assert_eq!(sum_dd, total_dd);
    }

    #[test]
    fn requires_lst() {
        let mut session = Session::new();
        let result = Predict.execute(&mut session, &["12.5", "45.0"]);
        assert!(result.is_err());
    }

    #[test]
    fn requires_coordinates() {
        let mut session = Session::new();
        session.lst_override = Some(Angle::from_hours(14.0));
        let result = Predict.execute(&mut session, &[]);
        assert!(result.is_err());
    }

    #[test]
    fn pier_from_ha_positive_is_east() {
        let ha = Angle::from_hours(2.0);
        assert_eq!(pier_from_ha(ha), PierSide::East);
    }

    #[test]
    fn pier_from_ha_negative_is_west() {
        let ha = Angle::from_hours(-2.0);
        assert_eq!(pier_from_ha(ha), PierSide::West);
    }

    #[test]
    fn format_ra_basic() {
        let ra = Angle::from_hours(12.5);
        assert_eq!(format_ra(ra), "12h 30m 00.00s");
    }

    #[test]
    fn format_dec_positive() {
        let dec = Angle::from_degrees(45.0);
        assert_eq!(format_dec(dec), "+45\u{00b0} 00' 00.0\"");
    }

    #[test]
    fn format_dec_negative() {
        let dec = Angle::from_degrees(-30.5);
        assert_eq!(format_dec(dec), "-30\u{00b0} 30' 00.0\"");
    }

    #[test]
    fn format_ha_positive() {
        let ha = Angle::from_hours(2.25);
        assert_eq!(format_ha(ha), "+02h 15m 00.00s");
    }

    #[test]
    fn format_ha_negative() {
        let ha = Angle::from_hours(-3.0);
        assert_eq!(format_ha(ha), "-03h 00m 00.00s");
    }
}
