//! Writer for TPOINT-style INDAT observation files.
//!
//! Inverse of [`crate::parser::parse_indat`]. Given an [`IndatFile`], emit a
//! whitespace-delimited file the parser can read back into an equivalent
//! `IndatFile`. The writer makes no judgment about reference frames or
//! corrections — what you pass in is what gets written. See the parser docs
//! (and `book/src/pointing/indat-format.md`) for the input-coordinate contract.
//!
//! The output is stable: same input always yields byte-for-byte identical
//! output (no current-time stamps, no nondeterministic ordering).

use std::io::{self, Write};

use celestial_core::Angle;
use celestial_time::scales::conversions::utc_tai::julian_to_calendar;

use crate::observation::{
    IndatFile, IndatOption, MountType, Observation, PierSide, SiteParams,
};

/// Writes `file` to `w` in TPOINT INDAT format. The result round-trips through
/// [`crate::parser::parse_indat`] back to a structurally-identical `IndatFile`
/// (within floating-point representation of DMS/HMS rounding).
pub fn write_indat(file: &IndatFile, w: &mut impl Write) -> io::Result<()> {
    write_header_block(w, file)?;
    write_site_line(w, &file.site, file.date.jd1(), file.date.jd2())?;
    for obs in &file.observations {
        write_observation_line(w, obs)?;
    }
    Ok(())
}

fn write_header_block(w: &mut impl Write, file: &IndatFile) -> io::Result<()> {
    for line in &file.header_lines {
        // Preserve the original prefix if present, otherwise mark it as a
        // comment so any TPOINT-compatible parser treats it as such.
        if line.starts_with('!') {
            writeln!(w, "{line}")?;
        } else {
            writeln!(w, "!{line}")?;
        }
    }
    for opt in &file.options {
        write_option(w, opt, file.mount_type)?;
    }
    // The Equatorial / Altaz mount type is encoded as a `:EQUAT` or `:ALTAZ`
    // option line. If no IndatOption echoing the mount type was present, emit
    // one so the parser recovers the same MountType.
    if !options_set_mount_type(&file.options) {
        write_mount_type_option(w, file.mount_type)?;
    }
    Ok(())
}

fn options_set_mount_type(options: &[IndatOption]) -> bool {
    options
        .iter()
        .any(|o| matches!(o, IndatOption::Equatorial | IndatOption::Altaz))
}

fn write_option(w: &mut impl Write, opt: &IndatOption, _mount: MountType) -> io::Result<()> {
    match opt {
        IndatOption::NoDA => writeln!(w, ":NODA"),
        IndatOption::AllSky => writeln!(w, ":ALLSKY"),
        IndatOption::Equinox => writeln!(w, ":EQUINOX"),
        IndatOption::Equatorial => writeln!(w, ":EQUAT"),
        IndatOption::Altaz => writeln!(w, ":ALTAZ"),
        IndatOption::RotatorTelescope => writeln!(w, ":ROTTEL"),
        IndatOption::RotatorNasmythLeft => writeln!(w, ":ROTNL"),
        IndatOption::RotatorNasmythRight => writeln!(w, ":ROTNR"),
        IndatOption::RotatorCoudeLeft => writeln!(w, ":ROTCL"),
        IndatOption::RotatorCoudeRight => writeln!(w, ":ROTCR"),
        IndatOption::Gimbal { z, y, x } => writeln!(
            w,
            ":GIMBAL {:.6} {:.6} {:.6}",
            z.degrees(),
            y.degrees(),
            x.degrees(),
        ),
    }
}

fn write_mount_type_option(w: &mut impl Write, mount: MountType) -> io::Result<()> {
    match mount {
        MountType::GermanEquatorial | MountType::ForkEquatorial => writeln!(w, ":EQUAT"),
        MountType::Altazimuth => writeln!(w, ":ALTAZ"),
    }
}

fn write_site_line(
    w: &mut impl Write,
    site: &SiteParams,
    jd1: f64,
    jd2: f64,
) -> io::Result<()> {
    let (year, month, day) = jd_to_ymd(jd1, jd2);
    let (lat_d, lat_m, lat_s) = angle_to_dms_components(site.latitude);
    writeln!(
        w,
        "{:+03} {:02} {:02} {} {} {} {:.2} {:.2} {:.2} {:.2} {:.4} {:.4}",
        lat_d,
        lat_m,
        lat_s.round() as i32,
        year,
        month,
        day,
        site.temperature,
        site.pressure,
        site.elevation,
        site.humidity,
        site.wavelength,
        site.lapse_rate,
    )
}

fn write_observation_line(w: &mut impl Write, obs: &Observation) -> io::Result<()> {
    let (tel_ra, tel_dec_raw) = encode_pier_side(obs.observed_ra, obs.observed_dec, obs.pier_side);
    writeln!(
        w,
        "{} {} {} {} {}",
        format_ra_hms(obs.catalog_ra),
        format_dec_dms(obs.catalog_dec.degrees()),
        format_ra_hms(tel_ra),
        format_dec_dms_raw(tel_dec_raw),
        format_lst_hm(obs.lst),
    )
}

/// Inverse of [`crate::observation::decode_pier_side`] for RA, plus the raw
/// declination magnitude with sign convention the parser expects.
fn encode_pier_side(observed_ra: Angle, observed_dec: Angle, pier: PierSide) -> (Angle, f64) {
    let dec_deg = observed_dec.degrees();
    match pier {
        PierSide::West => {
            let tel_ra = (observed_ra - Angle::from_hours(12.0)).normalized();
            let sign = if dec_deg >= 0.0 { 1.0 } else { -1.0 };
            let raw_dec = sign * (180.0 - dec_deg.abs());
            (tel_ra, raw_dec)
        }
        PierSide::East | PierSide::Unknown => (observed_ra, dec_deg),
    }
}

fn jd_to_ymd(jd1: f64, jd2: f64) -> (i32, i32, i32) {
    match julian_to_calendar(jd1, jd2) {
        Ok((y, m, d, _frac)) => (y, m, d),
        Err(_) => (0, 1, 1),
    }
}

fn angle_to_dms_components(a: Angle) -> (i32, u32, f64) {
    let deg = a.degrees();
    let sign = if deg < 0.0 { -1 } else { 1 };
    let mag = deg.abs();
    let d = mag.floor() as i32;
    let rem = (mag - d as f64) * 60.0;
    let m = rem.floor() as u32;
    let s = (rem - m as f64) * 60.0;
    (sign * d, m, s)
}

fn format_ra_hms(ra: Angle) -> String {
    let hours = ra.hours().rem_euclid(24.0);
    let h = hours.floor() as u32;
    let rem = (hours - h as f64) * 60.0;
    let m = rem.floor() as u32;
    let s = (rem - m as f64) * 60.0;
    format!("{:02} {:02} {:07.4}", h, m, s)
}

fn format_dec_dms(dec_deg: f64) -> String {
    let sign = if dec_deg < 0.0 { '-' } else { '+' };
    let mag = dec_deg.abs();
    let d = mag.floor() as u32;
    let rem = (mag - d as f64) * 60.0;
    let m = rem.floor() as u32;
    let s = (rem - m as f64) * 60.0;
    format!("{}{:02} {:02} {:06.3}", sign, d, m, s)
}

/// Like format_dec_dms but allows `|d|` to exceed 90, used for the
/// pier-side-encoded telescope declination column.
fn format_dec_dms_raw(raw_deg: f64) -> String {
    let sign = if raw_deg < 0.0 { '-' } else { '+' };
    let mag = raw_deg.abs();
    let d = mag.floor() as u32;
    let rem = (mag - d as f64) * 60.0;
    let m = rem.floor() as u32;
    let s = (rem - m as f64) * 60.0;
    format!("{}{:03} {:02} {:06.3}", sign, d, m, s)
}

fn format_lst_hm(lst: Angle) -> String {
    let hours = lst.hours().rem_euclid(24.0);
    let h = hours.floor() as u32;
    let m = (hours - h as f64) * 60.0;
    format!("{:02} {:06.3}", h, m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_indat;
    use celestial_time::JulianDate;

    fn make_site() -> SiteParams {
        SiteParams {
            latitude: Angle::from_degrees(39.0 + 26.0 / 3600.0),
            longitude: Angle::from_degrees(0.0),
            temperature: 22.0,
            pressure: 987.0,
            elevation: 231.65,
            humidity: 0.5,
            wavelength: 0.55,
            lapse_rate: 0.0065,
        }
    }

    fn make_obs(
        cat_ra_h: f64,
        cat_dec_deg: f64,
        obs_ra_h: f64,
        obs_dec_deg: f64,
        lst_h: f64,
        pier: PierSide,
    ) -> Observation {
        let lst = Angle::from_hours(lst_h);
        let catalog_ra = Angle::from_hours(cat_ra_h);
        let observed_ra = Angle::from_hours(obs_ra_h);
        Observation {
            catalog_ra,
            catalog_dec: Angle::from_degrees(cat_dec_deg),
            observed_ra,
            observed_dec: Angle::from_degrees(obs_dec_deg),
            lst,
            commanded_ha: (lst - catalog_ra).wrapped(),
            actual_ha: (lst - observed_ra).wrapped(),
            pier_side: pier,
            masked: false,
        }
    }

    fn make_file(options: Vec<IndatOption>, observations: Vec<Observation>) -> IndatFile {
        IndatFile {
            site: make_site(),
            options,
            observations,
            mount_type: MountType::GermanEquatorial,
            header_lines: vec!["!source: writer test".into()],
            date: JulianDate::from_calendar(2024, 7, 14, 0, 0, 0.0),
        }
    }

    fn write_to_string(file: &IndatFile) -> String {
        let mut buf = Vec::new();
        write_indat(file, &mut buf).expect("write");
        String::from_utf8(buf).expect("utf8")
    }

    #[test]
    fn header_comment_lines_get_bang_prefix_if_missing() {
        let mut file = make_file(vec![IndatOption::Equatorial], Vec::new());
        file.header_lines = vec!["plain text".into(), "!already prefixed".into()];
        let out = write_to_string(&file);
        assert!(out.contains("!plain text\n"));
        assert!(out.contains("!already prefixed\n"));
    }

    #[test]
    fn mount_type_inferred_when_no_option_present() {
        let file = make_file(Vec::new(), Vec::new());
        let out = write_to_string(&file);
        assert!(out.contains(":EQUAT\n"));
    }

    #[test]
    fn altaz_mount_emits_altaz_option() {
        let mut file = make_file(Vec::new(), Vec::new());
        file.mount_type = MountType::Altazimuth;
        let out = write_to_string(&file);
        assert!(out.contains(":ALTAZ\n"));
    }

    #[test]
    fn east_pier_passes_through() {
        let obs = make_obs(14.0, 30.0, 14.001, 30.005, 18.0, PierSide::East);
        let file = make_file(vec![IndatOption::Equatorial], vec![obs]);
        let parsed = parse_indat(&write_to_string(&file)).expect("parse");
        assert_eq!(parsed.observations.len(), 1);
        assert_eq!(parsed.observations[0].pier_side, PierSide::East);
    }

    #[test]
    fn west_pier_encodes_and_decodes() {
        let obs = make_obs(5.0, 30.0, 5.0, 30.0, 8.0, PierSide::West);
        let file = make_file(vec![IndatOption::Equatorial], vec![obs]);
        let parsed = parse_indat(&write_to_string(&file)).expect("parse");
        assert_eq!(parsed.observations[0].pier_side, PierSide::West);
        let dec = parsed.observations[0].observed_dec.degrees();
        assert!((dec - 30.0).abs() < 1e-6, "dec round-trip: got {dec}");
    }

    #[test]
    fn observation_round_trip_preserves_angles() {
        let obs = make_obs(
            14.0 + 36.0 / 60.0 + 55.4451 / 3600.0,
            72.0 + 3.0 / 60.0 + 24.120 / 3600.0,
            14.0 + 36.0 / 60.0 + 55.0 / 3600.0,
            72.0 + 3.0 / 60.0 + 24.0 / 3600.0,
            10.0 + 22.219 / 60.0,
            PierSide::East,
        );
        let file = make_file(vec![IndatOption::Equatorial], vec![obs]);
        let parsed = parse_indat(&write_to_string(&file)).expect("parse");
        let p = &parsed.observations[0];
        assert!((p.catalog_ra.hours() - (14.0 + 36.0 / 60.0 + 55.4451 / 3600.0)).abs() < 1e-7);
        assert!((p.catalog_dec.degrees() - (72.0 + 3.0 / 60.0 + 24.120 / 3600.0)).abs() < 1e-6);
    }

    #[test]
    fn site_line_round_trips() {
        let file = make_file(vec![IndatOption::Equatorial], Vec::new());
        let parsed = parse_indat(&write_to_string(&file)).expect("parse");
        assert_eq!(parsed.site.temperature, file.site.temperature);
        assert_eq!(parsed.site.pressure, file.site.pressure);
        assert_eq!(parsed.site.elevation, file.site.elevation);
        assert_eq!(parsed.site.humidity, file.site.humidity);
        assert_eq!(parsed.site.wavelength, file.site.wavelength);
        assert_eq!(parsed.site.lapse_rate, file.site.lapse_rate);
    }

    #[test]
    fn options_are_emitted_in_order() {
        let file = make_file(
            vec![IndatOption::NoDA, IndatOption::Equatorial],
            Vec::new(),
        );
        let out = write_to_string(&file);
        let noda_pos = out.find(":NODA").unwrap();
        let equat_pos = out.find(":EQUAT").unwrap();
        assert!(noda_pos < equat_pos);
    }

    #[test]
    fn noda_option_survives_round_trip() {
        let file = make_file(
            vec![IndatOption::NoDA, IndatOption::Equatorial],
            Vec::new(),
        );
        let parsed = parse_indat(&write_to_string(&file)).expect("parse");
        assert!(parsed.options.contains(&IndatOption::NoDA));
    }

    #[test]
    fn multiple_observations_round_trip() {
        let observations = vec![
            make_obs(14.5, 30.0, 14.5, 30.0, 18.0, PierSide::East),
            make_obs(2.0, -15.0, 2.0, -15.0, 6.0, PierSide::West),
            make_obs(22.0, 60.0, 22.0, 60.0, 1.5, PierSide::East),
        ];
        let file = make_file(vec![IndatOption::Equatorial], observations);
        let parsed = parse_indat(&write_to_string(&file)).expect("parse");
        assert_eq!(parsed.observations.len(), 3);
        assert_eq!(parsed.observations[0].pier_side, PierSide::East);
        assert_eq!(parsed.observations[1].pier_side, PierSide::West);
        assert_eq!(parsed.observations[2].pier_side, PierSide::East);
    }

    #[test]
    fn output_is_deterministic() {
        let obs = make_obs(14.0, 30.0, 14.001, 30.005, 18.0, PierSide::East);
        let file = make_file(vec![IndatOption::Equatorial], vec![obs]);
        let a = write_to_string(&file);
        let b = write_to_string(&file);
        assert_eq!(a, b);
    }
}
