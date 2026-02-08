use crate::error::{Error, Result};
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Outl;

impl Command for Outl {
    fn name(&self) -> &str { "OUTL" }
    fn description(&self) -> &str { "Identify outlier observations" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return Err(Error::Parse("OUTL requires a sigma threshold".into()));
        }
        let threshold: f64 = args[0].parse()
            .map_err(|e| Error::Parse(format!("invalid threshold: {}", e)))?;
        let do_mask = args.get(1)
            .is_some_and(|a| a.eq_ignore_ascii_case("M"));

        let fit = session.last_fit.as_ref()
            .ok_or_else(|| Error::Fit("no fit results available (run FIT first)".into()))?;
        let rms = fit.sky_rms;
        let cutoff = threshold * rms;

        let lat = session.latitude();
        let mut outliers: Vec<(usize, f64)> = Vec::new();

        for (i, obs) in session.observations.iter().enumerate() {
            if obs.masked {
                continue;
            }
            let h = obs.commanded_ha.radians();
            let dec = obs.catalog_dec.radians();
            let pier = obs.pier_side.sign();
            let (model_dh, model_dd) = session.model.apply_equatorial(h, dec, lat, pier);
            let raw_dh = (obs.actual_ha - obs.commanded_ha).arcseconds();
            let raw_dd = (obs.observed_dec - obs.catalog_dec).arcseconds();
            let dh = raw_dh - model_dh;
            let dd = raw_dd - model_dd;
            let dx = dh * dec.cos();
            let dr = (dx * dx + dd * dd).sqrt();
            if dr > cutoff {
                outliers.push((i, dr));
            }
        }

        if outliers.is_empty() {
            return Ok(CommandOutput::Text(format!(
                "No outliers (threshold {:.1} * {:.2}\" = {:.2}\")",
                threshold, rms, cutoff
            )));
        }

        let mut output = format!(
            "Outliers (residual > {:.1} * {:.2}\" = {:.2}\"):\n",
            threshold, rms, cutoff
        );
        for &(idx, dr) in &outliers {
            output += &format!("  obs {:>4}: {:.1}\"\n", idx + 1, dr);
        }

        if do_mask {
            for &(idx, _) in &outliers {
                session.observations[idx].masked = true;
            }
            output += &format!("\nMasked {} observations", outliers.len());
        } else {
            output += &format!("\nUse OUTL {:.1} M to mask", threshold);
        }

        Ok(CommandOutput::Text(output))
    }
}
