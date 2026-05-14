use super::{Command, CommandOutput};
use crate::error::Result;
use crate::observation::MountType;
use crate::session::Session;

pub struct Show;

impl Command for Show {
    fn name(&self) -> &str {
        "SHOW"
    }
    fn description(&self) -> &str {
        "Display session state"
    }

    fn execute(&self, session: &mut Session, _args: &[&str]) -> Result<CommandOutput> {
        let mount = match session.mount_type {
            MountType::GermanEquatorial => "German Equatorial",
            MountType::ForkEquatorial => "Fork Equatorial",
            MountType::Altazimuth => "Altazimuth",
        };

        let lat_str = session
            .site
            .as_ref()
            .map(|s| format_dms_lat(s.latitude.degrees()))
            .unwrap_or_else(|| "not set".to_string());

        let masked = session.masked_observation_count();
        let total = session.observation_count();
        let obs_str = if masked > 0 {
            format!("{} ({} masked)", total, masked)
        } else {
            format!("{}", total)
        };

        let rms_str = session
            .last_fit
            .as_ref()
            .map(|f| format!("{:.2}\"", f.sky_rms))
            .unwrap_or_else(|| "no fit yet".to_string());

        let output = format!(
            "Mount type: {}\nSite latitude: {}\nObservations: {}\nModel terms: {}\nLast fit RMS: {}",
            mount, lat_str, obs_str, session.model.term_count(), rms_str,
        );

        Ok(CommandOutput::Text(output))
    }
}

fn format_dms_lat(deg: f64) -> String {
    let sign = if deg < 0.0 { "-" } else { "+" };
    let total = deg.abs();
    let d = total as i32;
    let rem = (total - d as f64) * 60.0;
    let m = rem as i32;
    let s = (rem - m as f64) * 60.0;
    format!("{}{} {:02} {:02}", sign, d, m, s as i32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observation::SiteParams;
    use crate::solver::FitResult;
    use crate::test_support::{obs, FitResultBuilder};
    use celestial_core::Angle;

    fn site(lat_deg: f64) -> SiteParams {
        SiteParams {
            latitude: Angle::from_degrees(lat_deg),
            longitude: Angle::from_degrees(0.0),
            temperature: 10.0,
            pressure: 0.0,
            elevation: 0.0,
            humidity: 0.0,
            wavelength: 0.55,
            lapse_rate: 0.0065,
        }
    }

    fn run(s: &mut Session) -> String {
        match Show.execute(s, &[]).unwrap() {
            CommandOutput::Text(t) => t,
            other => panic!("expected Text, got {:?}", other),
        }
    }

    fn fake_fit(sky_rms: f64) -> FitResult {
        FitResultBuilder::new().sky_rms(sky_rms).build()
    }

    #[test]
    fn metadata() {
        assert_eq!(Show.name(), "SHOW");
        assert_eq!(Show.description(), "Display session state");
    }

    #[test]
    fn empty_session_renders_default_state() {
        let mut s = Session::new();
        let body = run(&mut s);
        assert!(body.contains("Mount type: German Equatorial"));
        assert!(body.contains("Site latitude: not set"));
        assert!(body.contains("Observations: 0"));
        assert!(body.contains("Model terms: 0"));
        assert!(body.contains("Last fit RMS: no fit yet"));
    }

    #[test]
    fn each_mount_type_renders_its_label() {
        for (mt, label) in [
            (MountType::GermanEquatorial, "German Equatorial"),
            (MountType::ForkEquatorial, "Fork Equatorial"),
            (MountType::Altazimuth, "Altazimuth"),
        ] {
            let mut s = Session::new();
            s.mount_type = mt;
            let body = run(&mut s);
            assert!(
                body.contains(&format!("Mount type: {}", label)),
                "missing label {}: {}",
                label,
                body,
            );
        }
    }

    #[test]
    fn latitude_renders_when_site_set() {
        let mut s = Session::new();
        s.site = Some(site(39.5));
        let body = run(&mut s);
        assert!(body.contains("Site latitude: +39 30 00"), "got {}", body);
    }

    #[test]
    fn negative_latitude_carries_minus_sign() {
        let mut s = Session::new();
        s.site = Some(site(-33.25));
        let body = run(&mut s);
        assert!(body.contains("Site latitude: -33 15 00"), "got {}", body);
    }

    #[test]
    fn observation_counts_show_masked_breakdown_only_when_nonzero() {
        let mut s = Session::new();
        s.observations = vec![obs(), obs(), obs()];
        let body = run(&mut s);
        assert!(body.contains("Observations: 3"));
        assert!(!body.contains("masked"));

        s.observations[1].masked = true;
        let body = run(&mut s);
        assert!(body.contains("Observations: 3 (1 masked)"), "got {}", body);
    }

    #[test]
    fn model_term_count_is_rendered() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.model.add_term("ID").unwrap();
        let body = run(&mut s);
        assert!(body.contains("Model terms: 2"));
    }

    #[test]
    fn last_fit_renders_sky_rms_with_two_decimals() {
        let mut s = Session::new();
        s.last_fit = Some(fake_fit(1.23456));
        let body = run(&mut s);
        assert!(body.contains("Last fit RMS: 1.23\""), "got {}", body);
    }

    #[test]
    fn args_are_ignored() {
        let mut s = Session::new();
        let with_args = match Show.execute(&mut s, &["whatever"]).unwrap() {
            CommandOutput::Text(t) => t,
            _ => panic!("Text"),
        };
        let no_args = run(&mut s);
        assert_eq!(with_args, no_args);
    }

    #[test]
    fn format_dms_lat_zero_is_positive() {
        assert_eq!(format_dms_lat(0.0), "+0 00 00");
    }

    #[test]
    fn format_dms_lat_fractional_components() {
        assert_eq!(format_dms_lat(39.5083), "+39 30 29");
    }
}
