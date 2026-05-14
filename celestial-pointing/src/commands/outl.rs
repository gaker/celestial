use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;

pub struct Outl;

const BADNESS_K: f64 = 2.5;

impl Command for Outl {
    fn name(&self) -> &str {
        "OUTL"
    }
    
    fn description(&self) -> &str {
        "Identify outlier observations"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        if args.is_empty() {
            return least_typical(session);
        }
        let threshold: f64 = args[0]
            .parse()
            .map_err(|e| Error::Parse(format!("invalid threshold: {}", e)))?;
        let do_mask = args.get(1).is_some_and(|a| a.eq_ignore_ascii_case("M"));

        let fit = session
            .last_fit
            .as_ref()
            .ok_or_else(|| Error::Fit("no fit results available (run FIT first)".into()))?;
        let rms = fit.sky_rms;
        let cutoff = threshold * rms;

        let lat = session.latitude();
        let prepared = session.prepared_observations();
        let mut outliers: Vec<(usize, f64)> = Vec::new();

        for (i, obs) in prepared.iter().enumerate() {
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
            let dx = dh * libm::cos(dec);
            let dr = libm::sqrt(dx * dx + dd * dd);
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

fn least_typical(session: &Session) -> Result<CommandOutput> {
    let fit = session
        .last_fit
        .as_ref()
        .ok_or_else(|| Error::Fit("no fit results available (run FIT first)".into()))?;
    if fit.diagnostics.is_empty() || fit.sky_rms <= 0.0 {
        return Ok(CommandOutput::Text("No diagnostics available".into()));
    }
    let cutoff = BADNESS_K * fit.sky_rms;
    let worst = fit
        .diagnostics
        .iter()
        .enumerate()
        .max_by(|a, b| {
            a.1.residual_sky
                .partial_cmp(&b.1.residual_sky)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, d)| (i, d.residual_sky / cutoff))
        .unwrap();
    Ok(CommandOutput::Text(format!(
        "Least typical observation is #{}, badness {:.3}.",
        worst.0 + 1,
        worst.1,
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observation::Observation;
    use crate::solver::{FitResult, ObsDiagnostic};
    use crate::test_support::{FitResultBuilder, ObsBuilder};

    fn obs(actual_offset_arcsec: f64) -> Observation {
        ObsBuilder::new()
            .actual_ha_arcsec(actual_offset_arcsec)
            .build()
    }

    fn fit_with(sky_rms: f64, diagnostics: Vec<ObsDiagnostic>) -> FitResult {
        FitResultBuilder::new()
            .sky_rms(sky_rms)
            .popn_sd(1.0)
            .diagnostics(diagnostics)
            .build()
    }

    fn text(out: CommandOutput) -> String {
        match out {
            CommandOutput::Text(s) => s,
            other => panic!("expected Text, got {:?}", other),
        }
    }

    fn parse_msg(err: &Error) -> &str {
        match err {
            Error::Parse(m) => m,
            other => panic!("expected Parse, got {:?}", other),
        }
    }

    fn fit_err(err: &Error) -> &str {
        match err {
            Error::Fit(m) => m,
            other => panic!("expected Fit, got {:?}", other),
        }
    }

    fn diag(residual: f64) -> ObsDiagnostic {
        ObsDiagnostic {
            residual_sky: residual,
            leverage: 0.0,
            studentized: 0.0,
        }
    }

    #[test]
    fn metadata() {
        assert_eq!(Outl.name(), "OUTL");
        assert_eq!(Outl.description(), "Identify outlier observations");
    }

    #[test]
    fn no_args_with_no_fit_errors() {
        let mut s = Session::new();
        let err = Outl.execute(&mut s, &[]).unwrap_err();
        assert!(fit_err(&err).contains("no fit results"));
    }

    #[test]
    fn no_args_with_empty_diagnostics_reports_unavailable() {
        let mut s = Session::new();
        s.last_fit = Some(fit_with(5.0, vec![]));
        let body = text(Outl.execute(&mut s, &[]).unwrap());
        assert_eq!(body, "No diagnostics available");
    }

    #[test]
    fn no_args_with_zero_sky_rms_reports_unavailable() {
        let mut s = Session::new();
        s.last_fit = Some(fit_with(0.0, vec![diag(1.0)]));
        let body = text(Outl.execute(&mut s, &[]).unwrap());
        assert_eq!(body, "No diagnostics available");
    }

    #[test]
    fn no_args_reports_least_typical_by_residual() {
        let mut s = Session::new();
        s.last_fit = Some(fit_with(2.0, vec![diag(1.0), diag(8.0), diag(3.0)]));
        let body = text(Outl.execute(&mut s, &[]).unwrap());
        // worst is index 1 (residual 8.0); badness = 8 / (2.5 * 2.0) = 1.6
        assert!(body.contains("#2"), "got {}", body);
        assert!(body.contains("badness 1.600"), "got {}", body);
    }

    #[test]
    fn unparseable_threshold_errors() {
        let mut s = Session::new();
        let err = Outl.execute(&mut s, &["banana"]).unwrap_err();
        assert!(parse_msg(&err).contains("invalid threshold"));
    }

    #[test]
    fn threshold_path_without_fit_errors() {
        let mut s = Session::new();
        let err = Outl.execute(&mut s, &["3.0"]).unwrap_err();
        assert!(fit_err(&err).contains("no fit results"));
    }

    #[test]
    fn no_outliers_reports_clean() {
        let mut s = Session::new();
        s.observations = vec![obs(10.0)]; // raw residual 10"
        s.last_fit = Some(fit_with(5.0, vec![]));
        // threshold 3.0 * rms 5.0 = 15" cutoff; residual 10" < 15".
        let body = text(Outl.execute(&mut s, &["3.0"]).unwrap());
        assert!(body.contains("No outliers"));
        assert!(body.contains("15.00"));
        assert!(s.observations.iter().all(|o| !o.masked));
    }

    #[test]
    fn outliers_listed_without_m_flag_does_not_mask() {
        let mut s = Session::new();
        s.observations = vec![obs(10.0), obs(100.0)];
        s.last_fit = Some(fit_with(5.0, vec![]));
        let body = text(Outl.execute(&mut s, &["3.0"]).unwrap());
        assert!(body.contains("Outliers"));
        assert!(body.contains("obs    2"));
        assert!(!body.contains("obs    1:"));
        assert!(body.contains("Use OUTL 3.0 M to mask"));
        assert!(s.observations.iter().all(|o| !o.masked));
    }

    #[test]
    fn m_flag_masks_outliers_case_insensitive() {
        for token in ["M", "m"] {
            let mut s = Session::new();
            s.observations = vec![obs(10.0), obs(100.0)];
            s.last_fit = Some(fit_with(5.0, vec![]));
            let body = text(Outl.execute(&mut s, &["3.0", token]).unwrap());
            assert!(body.contains("Masked 1 observations"), "token {:?}: {}", token, body);
            assert!(!s.observations[0].masked);
            assert!(s.observations[1].masked);
        }
    }

    // Non-M second arg is ignored (boolean is `eq_ignore_ascii_case("M")`).
    #[test]
    fn non_m_second_arg_does_not_trigger_masking() {
        let mut s = Session::new();
        s.observations = vec![obs(100.0)];
        s.last_fit = Some(fit_with(5.0, vec![]));
        let body = text(Outl.execute(&mut s, &["3.0", "mask"]).unwrap());
        assert!(body.contains("Use OUTL"));
        assert!(!s.observations[0].masked);
    }

    // Already-masked observations are skipped in the threshold scan.
    #[test]
    fn already_masked_observations_are_skipped() {
        let mut s = Session::new();
        s.observations = vec![obs(100.0), obs(100.0)];
        s.observations[0].masked = true;
        s.last_fit = Some(fit_with(5.0, vec![]));
        let body = text(Outl.execute(&mut s, &["3.0"]).unwrap());
        assert!(body.contains("obs    2"));
        assert!(!body.contains("obs    1:"));
    }
}
