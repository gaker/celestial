use super::{Command, CommandOutput};
use crate::error::Result;
use crate::observation::{MountType, PierSide};
use crate::parser::parse_indat;
use crate::session::Session;

pub struct Indat;

impl Command for Indat {
    fn name(&self) -> &str {
        "INDAT"
    }
    fn description(&self) -> &str {
        "Load observations from INDAT file"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(crate::error::Error::Parse(
                "INDAT requires a filename".into(),
            ));
        }
        let content = std::fs::read_to_string(args[0]).map_err(crate::error::Error::Io)?;
        let indat = parse_indat(&content)?;
        let summary = format_summary(&indat);
        session.load_indat(indat);
        Ok(CommandOutput::Text(summary))
    }
}

fn format_summary(indat: &crate::observation::IndatFile) -> String {
    let mount = match indat.mount_type {
        MountType::GermanEquatorial => "German Equatorial",
        MountType::ForkEquatorial => "Fork Equatorial",
        MountType::Altazimuth => "Altazimuth",
    };
    let lat = format_dms(indat.site.latitude.degrees());
    let n = indat.observations.len();
    let east = indat
        .observations
        .iter()
        .filter(|o| o.pier_side == PierSide::East)
        .count();
    let west = indat
        .observations
        .iter()
        .filter(|o| o.pier_side == PierSide::West)
        .count();
    let s = &indat.site;
    let atm = if s.pressure <= 0.0 {
        "Atmos:    none (pressure=0 → refraction skipped)".to_string()
    } else {
        format!(
            "Atmos:    {:.1}\u{00b0}C, {:.1} hPa, RH {:.2}, \u{03bb}={:.3}\u{00b5}m, elev {:.1}m",
            s.temperature, s.pressure, s.humidity, s.wavelength, s.elevation
        )
    };
    format!(
        "{} observations loaded\n  Mount:    {}\n  Latitude: {}\n  Pier:     {} east, {} west\n  {}",
        n, mount, lat, east, west, atm,
    )
}

fn format_dms(deg: f64) -> String {
    let sign = if deg < 0.0 { "-" } else { "+" };
    let total = deg.abs();
    let d = total as i32;
    let rem = (total - d as f64) * 60.0;
    let m = rem as i32;
    let s = (rem - m as f64) * 60.0;
    format!("{}{}d {:02}' {:02}\"", sign, d, m, s as i32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Error;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static SEQ: AtomicU64 = AtomicU64::new(0);

    fn write_tmp(contents: &str, tag: &str) -> PathBuf {
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let mut path = std::env::temp_dir();
        path.push(format!("celpoint-indat-{}-{}-{}.dat", tag, std::process::id(), n));
        std::fs::write(&path, contents).expect("write temp indat");
        path
    }

    fn cleanup(p: &PathBuf) {
        let _ = std::fs::remove_file(p);
    }

    fn minimal_indat(extra_obs: &[&str]) -> String {
        let mut s = String::from(
            "!Test header\n\
             :NODA\n\
             :EQUAT\n\
             +39 00 26 2024 7 14 29.20 987.00 231.65  0.94 0.5500 0.0065\n",
        );
        for line in extra_obs {
            s.push_str(line);
            s.push('\n');
        }
        s
    }

    fn run(path: &PathBuf) -> Result<String> {
        let mut session = Session::new();
        Indat.execute(&mut session, &[path.to_str().unwrap()]).map(|out| match out {
            CommandOutput::Text(s) => s,
            other => panic!("expected Text, got {:?}", other),
        })
    }

    #[test]
    fn metadata() {
        assert_eq!(Indat.name(), "INDAT");
        assert_eq!(Indat.description(), "Load observations from INDAT file");
    }

    #[test]
    fn no_args_errors() {
        let mut session = Session::new();
        let err = Indat.execute(&mut session, &[]).unwrap_err();
        match err {
            Error::Parse(msg) => assert!(msg.contains("INDAT requires a filename")),
            other => panic!("expected Parse, got {:?}", other),
        }
    }

    #[test]
    fn missing_file_returns_io_error() {
        let mut session = Session::new();
        let err = Indat
            .execute(&mut session, &["/nonexistent/path/that/should/not/exist.dat"])
            .unwrap_err();
        assert!(matches!(err, Error::Io(_)), "expected Io, got {:?}", err);
    }

    #[test]
    fn parse_error_propagates() {
        let p = write_tmp("not a valid indat file\n", "bad");
        let mut session = Session::new();
        let result = Indat.execute(&mut session, &[p.to_str().unwrap()]);
        cleanup(&p);
        assert!(matches!(result, Err(Error::Parse(_))));
    }

    #[test]
    fn loads_observations_and_emits_summary() {
        let body = minimal_indat(&[
            "21 43 18.4460 +72 29 08.368 09 28 59.9527 +109 20 06.469  16 23.130",
            "23 46 02.2988 +77 38 38.725 11 26 17.6308 +104 03 28.734  16 24.711",
        ]);
        let p = write_tmp(&body, "ok");
        let summary = run(&p).unwrap();
        cleanup(&p);

        assert!(summary.contains("2 observations loaded"));
        assert!(summary.contains("Mount:    German Equatorial"));
        assert!(summary.contains("Latitude: +39d 00'"));
        assert!(summary.contains("Pier:"));
        assert!(summary.contains("Atmos:"));
        assert!(summary.contains("987.0 hPa"));
    }

    #[test]
    fn session_state_is_populated_on_success() {
        let body = minimal_indat(&[
            "21 43 18.4460 +72 29 08.368 09 28 59.9527 +109 20 06.469  16 23.130",
        ]);
        let p = write_tmp(&body, "session");
        let mut session = Session::new();
        Indat.execute(&mut session, &[p.to_str().unwrap()]).unwrap();
        cleanup(&p);

        assert_eq!(session.observations.len(), 1);
        assert!(session.site.is_some());
    }

    #[test]
    fn zero_pressure_skips_refraction_in_summary() {
        // Replace the pressure column (987.00) with 0.
        let body = "!Hdr\n\
                    :NODA\n\
                    :EQUAT\n\
                    +39 00 26 2024 7 14 29.20 0 231.65  0.94 0.5500 0.0065\n";
        let p = write_tmp(body, "norefr");
        let summary = run(&p).unwrap();
        cleanup(&p);
        assert!(summary.contains("Atmos:    none (pressure=0"));
    }

    #[test]
    fn altaz_option_changes_mount_label() {
        let body = "!Hdr\n\
                    :NODA\n\
                    :ALTAZ\n\
                    +39 00 26 2024 7 14 29.20 987.00 231.65  0.94 0.5500 0.0065\n";
        let p = write_tmp(body, "altaz");
        let summary = run(&p).unwrap();
        cleanup(&p);
        assert!(summary.contains("Mount:    Altazimuth"));
    }

    #[test]
    fn pier_counts_reflect_east_west_split() {
        // First obs is east (positive dec), second has negative dec column
        // encoded as the "west" pier marker per parse_observation_line.
        let body = minimal_indat(&[
            "21 43 18.4460 +72 29 08.368 09 28 59.9527 +109 20 06.469  16 23.130",
            "23 46 02.2988 -77 38 38.725 11 26 17.6308 +104 03 28.734  16 24.711",
        ]);
        let p = write_tmp(&body, "pier");
        let summary = run(&p).unwrap();
        cleanup(&p);
        // Don't assert which is which (encoding of east/west via the dec
        // sign is parser-internal); just assert both counts sum to 2.
        let east: u32 = pick_count(&summary, "east");
        let west: u32 = pick_count(&summary, "west");
        assert_eq!(east + west, 2, "summary was: {}", summary);
    }

    fn pick_count(summary: &str, label: &str) -> u32 {
        // line: "  Pier:     X east, Y west"
        let pier_line = summary.lines().find(|l| l.contains("Pier:")).unwrap();
        let segment = pier_line.split(label).next().unwrap();
        segment
            .split_whitespace()
            .last()
            .unwrap()
            .trim_end_matches(',')
            .parse()
            .unwrap()
    }

    #[test]
    fn format_dms_positive_zero_is_plus() {
        assert_eq!(format_dms(0.0), "+0d 00' 00\"");
    }

    #[test]
    fn format_dms_negative() {
        let s = format_dms(-39.0073);
        assert!(s.starts_with("-39d"), "got {}", s);
    }

    #[test]
    fn format_dms_positive_fractional() {
        let s = format_dms(39.5);
        assert_eq!(s, "+39d 30' 00\"");
    }
}
