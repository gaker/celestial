use crate::error::Result;
use crate::observation::{MountType, PierSide};
use crate::parser::parse_indat;
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Indat;

impl Command for Indat {
    fn name(&self) -> &str { "INDAT" }
    fn description(&self) -> &str { "Load observations from INDAT file" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(crate::error::Error::Parse("INDAT requires a filename".into()));
        }
        let content = std::fs::read_to_string(args[0])
            .map_err(crate::error::Error::Io)?;
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
    let east = indat.observations.iter()
        .filter(|o| o.pier_side == PierSide::East).count();
    let west = indat.observations.iter()
        .filter(|o| o.pier_side == PierSide::West).count();
    format!(
        "{} observations loaded\n  Mount:    {}\n  Latitude: {}\n  Pier:     {} east, {} west",
        n, mount, lat, east, west,
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
