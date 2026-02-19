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
