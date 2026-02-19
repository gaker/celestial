use super::{Command, CommandOutput};
use crate::error::Result;
use crate::model::PointingModel;
use crate::observation::{Observation, PierSide};
use crate::session::Session;
use celestial_core::constants::{DEG_TO_RAD, RAD_TO_DEG};
use celestial_core::Angle;

pub struct Slist;

impl Command for Slist {
    fn name(&self) -> &str {
        "SLIST"
    }
    fn description(&self) -> &str {
        "List observations with residuals"
    }

    fn execute(&self, session: &mut Session, _args: &[&str]) -> Result<CommandOutput> {
        let lat = session.latitude();
        let header = format!(
            "{:>5} {:>15} {:>15} {:>7} {:>7} {:>8} {:>8} {:>8} {:>8} {:>8}",
            "", "*HA", "*Dec", "*Az", "*ZD", "dX", "dD", "dS", "dZ", "dR"
        );
        let mut output = header + "\n\n";
        for (i, obs) in session.observations.iter().enumerate() {
            let row = format_row(i + 1, obs, &session.model, lat);
            output += &row;
            output += "\n";
        }
        Ok(CommandOutput::Text(output))
    }
}

fn format_row(num: usize, obs: &Observation, model: &PointingModel, lat: f64) -> String {
    let h = obs.commanded_ha.radians();
    let dec = obs.catalog_dec.radians();
    let pier = obs.pier_side.sign();
    let (model_dh, model_dd) = model.apply_equatorial(h, dec, lat, pier);
    let (raw_dh, raw_dd) = compute_raw_residuals(obs);
    let dh = raw_dh - model_dh;
    let dd = raw_dd - model_dd;
    let dx = dh * libm::cos(dec);
    let (az, zd) = compute_az_zd(h, dec, lat);
    let dr = libm::sqrt(dx * dx + dd * dd);
    let pier_char = pier_indicator(obs.pier_side);
    let mask_char = if obs.masked { "*" } else { "" };

    format!(
        "{:>4}{}{} {:>15} {:>15} {:>7.1} {:>7.1} {:>8.1} {:>8.1} {:>8.1} {:>8.1} {:>8.1}",
        num,
        pier_char,
        mask_char,
        format_hms(obs.commanded_ha),
        format_dms(obs.catalog_dec),
        az * RAD_TO_DEG,
        zd * RAD_TO_DEG,
        dx,
        dd,
        dx,
        zd * RAD_TO_DEG,
        dr
    )
}

fn compute_raw_residuals(obs: &Observation) -> (f64, f64) {
    let raw_dh = (obs.actual_ha - obs.commanded_ha).arcseconds();
    let raw_dd = (obs.observed_dec - obs.catalog_dec).arcseconds();
    (raw_dh, raw_dd)
}

fn compute_az_zd(h: f64, dec: f64, lat: f64) -> (f64, f64) {
    let sin_alt = libm::sin(lat) * libm::sin(dec) + libm::cos(lat) * libm::cos(dec) * libm::cos(h);
    let alt = libm::asin(sin_alt);
    let zd = 90.0 * DEG_TO_RAD - alt;
    let cos_alt = libm::cos(alt);
    let (sin_az, cos_az) = if cos_alt.abs() < 1e-10 {
        (0.0, 1.0)
    } else {
        let sa = -(libm::cos(dec) * libm::sin(h)) / cos_alt;
        let ca = (libm::sin(dec) - libm::sin(lat) * sin_alt) / (libm::cos(lat) * cos_alt);
        (sa, ca)
    };
    let az = libm::atan2(sin_az, cos_az);
    (az, zd)
}

fn pier_indicator(pier_side: PierSide) -> &'static str {
    match pier_side {
        PierSide::West => "b",
        PierSide::East => " ",
        PierSide::Unknown => "?",
    }
}

fn format_hms(angle: Angle) -> String {
    let total_sec = angle.hours().abs() * 3600.0;
    let h = (total_sec / 3600.0) as i32;
    let m = ((total_sec - h as f64 * 3600.0) / 60.0) as i32;
    let s = total_sec - h as f64 * 3600.0 - m as f64 * 60.0;
    let sign = if angle.radians() < 0.0 { "-" } else { "+" };
    format!("{}{:02} {:02} {:05.2}", sign, h, m, s)
}

fn format_dms(angle: Angle) -> String {
    let total_arcsec = angle.degrees().abs() * 3600.0;
    let d = (total_arcsec / 3600.0) as i32;
    let m = ((total_arcsec - d as f64 * 3600.0) / 60.0) as i32;
    let s = total_arcsec - d as f64 * 3600.0 - m as f64 * 60.0;
    let sign = if angle.radians() < 0.0 { "-" } else { "+" };
    format!("{}{:02} {:02} {:04.1}", sign, d, m, s)
}
