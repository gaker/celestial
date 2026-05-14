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
        let prepared = session.prepared_observations();
        let header = format!(
            "{:>5} {:>15} {:>15} {:>7} {:>7} {:>8} {:>8} {:>8} {:>8} {:>8}",
            "", "*HA", "*Dec", "*Az", "*ZD", "dX", "dD", "dS", "dZ", "dR"
        );
        let mut output = header + "\n\n";
        let mut sums = RmsSums::default();
        let mut count = 0_usize;
        for (i, obs) in prepared.iter().enumerate() {
            let row = compute_row(obs, &session.model, lat);
            output += &format_row(i + 1, obs, &row);
            output += "\n";
            if !obs.masked {
                sums.accumulate(&row);
                count += 1;
            }
        }
        if count > 0 {
            output += &format_footer(&sums, count);
        }
        Ok(CommandOutput::Text(output))
    }
}

struct Row {
    az: f64,
    zd: f64,
    dx: f64,
    dd: f64,
    ds: f64,
    dz: f64,
    dr: f64,
}

#[derive(Default)]
struct RmsSums {
    dx: f64,
    dd: f64,
    ds: f64,
    dz: f64,
    dr: f64,
}

impl RmsSums {
    fn accumulate(&mut self, r: &Row) {
        self.dx += r.dx * r.dx;
        self.dd += r.dd * r.dd;
        self.ds += r.ds * r.ds;
        self.dz += r.dz * r.dz;
        self.dr += r.dr * r.dr;
    }
}

fn compute_row(obs: &Observation, model: &PointingModel, lat: f64) -> Row {
    let h = obs.commanded_ha.radians();
    let dec = obs.catalog_dec.radians();
    let pier = obs.pier_side.sign();
    let (model_dh, model_dd) = model.apply_equatorial(h, dec, lat, pier);
    let (raw_dh, raw_dd) = compute_raw_residuals(obs);
    let dh = raw_dh - model_dh;
    let dd = raw_dd - model_dd;
    let dx = dh * libm::cos(dec);
    let (az, zd) = compute_az_zd(h, dec, lat);
    let (ds, dz) = rotate_to_horizon(dx, dd, h, dec, lat);
    let dr = libm::sqrt(dx * dx + dd * dd);
    Row { az, zd, dx, dd, ds, dz, dr }
}

fn format_row(num: usize, obs: &Observation, r: &Row) -> String {
    let pier_char = pier_indicator(obs.pier_side);
    let mask_char = if obs.masked { "*" } else { "" };
    format!(
        "{:>4}{}{} {:>15} {:>15} {:>7.1} {:>7.1} {:>8.1} {:>8.1} {:>8.1} {:>8.1} {:>8.1}",
        num,
        pier_char,
        mask_char,
        format_hms(obs.commanded_ha),
        format_dms(obs.catalog_dec),
        r.az * RAD_TO_DEG,
        r.zd * RAD_TO_DEG,
        r.dx,
        r.dd,
        r.ds,
        r.dz,
        r.dr,
    )
}

fn format_footer(sums: &RmsSums, n: usize) -> String {
    let n_f = n as f64;
    let rms = |sq: f64| libm::sqrt(sq / n_f);
    format!(
        "\n{:>54}{:>9}{:>9}{:>9}{:>9}{:>9}\n{:>54}{:>9.2}{:>9.2}{:>9.2}{:>9.2}{:>9.2}\n",
        "",
        "dX", "dD", "dS", "dZ", "dR",
        "RMS",
        rms(sums.dx), rms(sums.dd), rms(sums.ds), rms(sums.dz), rms(sums.dr),
    )
}

fn rotate_to_horizon(dx: f64, dd: f64, h: f64, dec: f64, lat: f64) -> (f64, f64) {
    let sin_alt =
        libm::sin(lat) * libm::sin(dec) + libm::cos(lat) * libm::cos(dec) * libm::cos(h);
    let cos_alt = libm::sqrt((1.0 - sin_alt * sin_alt).max(0.0));
    if cos_alt < 1e-10 || libm::cos(dec).abs() < 1e-10 {
        return (0.0, 0.0);
    }
    let sin_q = libm::sin(h) * libm::cos(lat) / cos_alt;
    let cos_q = (libm::sin(lat) - libm::sin(dec) * sin_alt) / (libm::cos(dec) * cos_alt);
    let ds = -dx * cos_q - dd * sin_q;
    let dz = dx * sin_q - dd * cos_q;
    (ds, dz)
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

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::HALF_PI;

    fn obs(
        commanded_ha_arcsec: f64,
        actual_ha_arcsec: f64,
        cat_dec_deg: f64,
        obs_dec_deg: f64,
        pier: PierSide,
    ) -> Observation {
        crate::test_support::ObsBuilder::new()
            .commanded_ha_arcsec(commanded_ha_arcsec)
            .actual_ha_arcsec(actual_ha_arcsec)
            .catalog_dec_deg(cat_dec_deg)
            .observed_dec_deg(obs_dec_deg)
            .pier(pier)
            .build()
    }

    fn text(out: CommandOutput) -> String {
        match out {
            CommandOutput::Text(s) => s,
            other => panic!("expected Text, got {:?}", other),
        }
    }

    #[test]
    fn metadata() {
        assert_eq!(Slist.name(), "SLIST");
        assert_eq!(Slist.description(), "List observations with residuals");
    }

    #[test]
    fn empty_session_emits_header_only() {
        let mut s = Session::new();
        let body = text(Slist.execute(&mut s, &[]).unwrap());
        assert!(body.contains("*HA"));
        assert!(body.contains("*Dec"));
        assert!(body.contains("dX"));
        assert!(body.contains("dR"));
        // No data rows and no RMS footer when count == 0.
        assert!(!body.contains("RMS"));
    }

    #[test]
    fn args_are_ignored() {
        let mut s = Session::new();
        let with = text(Slist.execute(&mut s, &["junk"]).unwrap());
        let without = text(Slist.execute(&mut s, &[]).unwrap());
        assert_eq!(with, without);
    }

    #[test]
    fn single_obs_renders_row_and_rms_footer() {
        let mut s = Session::new();
        // 100" HA offset, 30 deg dec, east pier
        s.observations = vec![obs(0.0, 100.0, 30.0, 30.0, PierSide::East)];
        let body = text(Slist.execute(&mut s, &[]).unwrap());

        // Row 1 should be present.
        assert!(body.contains("   1 "));
        // dX should be ~100 * cos(30°) = ~86.6
        assert!(body.contains("86.6"), "expected dX≈86.6, got:\n{}", body);
        // Footer present
        assert!(body.contains("RMS"));
    }

    #[test]
    fn masked_observation_renders_with_asterisk_and_skips_rms() {
        let mut s = Session::new();
        let mut o = obs(0.0, 100.0, 30.0, 30.0, PierSide::East);
        o.masked = true;
        s.observations = vec![o];
        let body = text(Slist.execute(&mut s, &[]).unwrap());

        // Row marker has an asterisk after the number.
        let row1 = body
            .lines()
            .find(|l| l.trim_start().starts_with("1"))
            .unwrap_or("");
        assert!(row1.contains('*'), "expected masked marker, got: {}", row1);
        // No RMS footer because no unmasked obs.
        assert!(!body.contains("RMS"), "should not emit footer:\n{}", body);
    }

    #[test]
    fn west_pier_renders_b_indicator() {
        assert_eq!(pier_indicator(PierSide::West), "b");
        assert_eq!(pier_indicator(PierSide::East), " ");
        assert_eq!(pier_indicator(PierSide::Unknown), "?");
    }

    #[test]
    fn unknown_pier_appears_in_row_marker() {
        let mut s = Session::new();
        s.observations = vec![obs(0.0, 50.0, 30.0, 30.0, PierSide::Unknown)];
        let body = text(Slist.execute(&mut s, &[]).unwrap());
        let row1 = body
            .lines()
            .find(|l| l.trim_start().starts_with("1"))
            .unwrap();
        assert!(row1.contains('?'), "expected ? marker, got: {}", row1);
    }

    #[test]
    fn rms_footer_uses_unmasked_count() {
        let mut s = Session::new();
        let unmasked = obs(0.0, 100.0, 30.0, 30.0, PierSide::East);
        let mut masked = obs(0.0, 9999.0, 30.0, 30.0, PierSide::East);
        masked.masked = true;
        s.observations = vec![unmasked, masked];

        let body = text(Slist.execute(&mut s, &[]).unwrap());
        // RMS computed from a single point with dX≈86.6 should be ≈86.60,
        // not influenced by the masked 9999" obs.
        assert!(body.contains("86.60"), "RMS line missing or wrong:\n{}", body);
    }

    #[test]
    fn rotate_to_horizon_short_circuits_at_horizon() {
        // dec = lat → at HA=0 the star is at zenith → cos_alt ~ 0 → returns (0,0)
        let (ds, dz) = rotate_to_horizon(1.0, 1.0, 0.0, 0.5, 0.5);
        assert_eq!((ds, dz), (0.0, 0.0));
    }

    #[test]
    fn rotate_to_horizon_short_circuits_at_celestial_pole() {
        // dec at the pole → cos(dec) ≈ 0 → returns (0,0)
        let (ds, dz) = rotate_to_horizon(1.0, 1.0, 0.5, HALF_PI - 1e-15, 0.5);
        assert_eq!((ds, dz), (0.0, 0.0));
    }

    #[test]
    fn rotate_to_horizon_nonsingular_case_produces_finite_values() {
        let (ds, dz) = rotate_to_horizon(10.0, 5.0, 0.5, 0.6, 0.7);
        assert!(ds.is_finite());
        assert!(dz.is_finite());
    }

    #[test]
    fn compute_az_zd_pole_singularity_safe() {
        // alt ≈ 90° (zenith) → cos_alt ≈ 0 → fallback (sin_az, cos_az) = (0,1)
        let (az, zd) = compute_az_zd(0.0, 0.5, 0.5);
        assert!(az.is_finite());
        assert!(zd.is_finite());
        assert!(zd >= 0.0);
    }

    #[test]
    fn format_hms_positive_and_negative() {
        assert_eq!(format_hms(Angle::from_hours(0.0)), "+00 00 00.00");
        let neg = format_hms(Angle::from_hours(-1.5));
        assert!(neg.starts_with("-01"), "got {}", neg);
    }

    #[test]
    fn format_dms_positive_and_negative() {
        assert_eq!(format_dms(Angle::from_degrees(0.0)), "+00 00 00.0");
        let neg = format_dms(Angle::from_degrees(-30.5));
        assert!(neg.starts_with("-30"), "got {}", neg);
    }

    #[test]
    fn rms_sums_accumulate_sums_squares() {
        let mut sums = RmsSums::default();
        sums.accumulate(&Row {
            az: 0.0,
            zd: 0.0,
            dx: 3.0,
            dd: 4.0,
            ds: 0.0,
            dz: 0.0,
            dr: 5.0,
        });
        assert_eq!(sums.dx, 9.0);
        assert_eq!(sums.dd, 16.0);
        assert_eq!(sums.dr, 25.0);
    }
}
