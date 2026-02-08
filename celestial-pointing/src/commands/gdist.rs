use std::path::Path;

use crate::error::Result;
use crate::plot::residuals::{compute_residuals, require_fit};
use crate::session::Session;
use super::{Command, CommandOutput};

pub struct Gdist;

impl Command for Gdist {
    fn name(&self) -> &str { "GDIST" }
    fn description(&self) -> &str { "Histogram of residual distribution" }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        require_fit(session)?;
        let residuals = compute_residuals(session);
        if residuals.is_empty() {
            return Ok(CommandOutput::Text("No active observations".to_string()));
        }
        let dx_vals: Vec<f64> = residuals.iter().map(|r| r.dx).collect();
        let dd_vals: Vec<f64> = residuals.iter().map(|r| r.dd).collect();

        match args.first() {
            Some(path) => svg_output(path, args, &dx_vals, &dd_vals),
            None => terminal_output(&dx_vals, &dd_vals),
        }
    }
}

fn terminal_output(dx: &[f64], dd: &[f64]) -> Result<CommandOutput> {
    let dx_hist = crate::plot::terminal::histogram_terminal(dx, "dX Distribution", "dX");
    let dd_hist = crate::plot::terminal::histogram_terminal(dd, "dDec Distribution", "dDec");
    Ok(CommandOutput::Text(format!("{dx_hist}\n{dd_hist}")))
}

fn svg_output(
    path: &str,
    args: &[&str],
    dx: &[f64],
    dd: &[f64],
) -> Result<CommandOutput> {
    let is_dec = args.get(1).is_some_and(|a| a.eq_ignore_ascii_case("D"));
    let (values, title, label) = if is_dec {
        (dd, "dDec Distribution", "dDec (arcsec)")
    } else {
        (dx, "dX Distribution", "dX (arcsec)")
    };
    write_svg(values, Path::new(path), title, label)
}

fn write_svg(
    values: &[f64],
    path: &Path,
    title: &str,
    x_label: &str,
) -> Result<CommandOutput> {
    crate::plot::svg::histogram_svg(values, path, title, x_label)
        .map_err(|e| crate::error::Error::Io(std::io::Error::other(e.to_string())))?;
    Ok(CommandOutput::Text(format!("Written to {}", path.display())))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::session::Session;

    #[test]
    fn no_fit_returns_error() {
        let mut session = Session::new();
        let result = Gdist.execute(&mut session, &[]);
        let err = result.err().expect("expected error");
        assert!(err.to_string().contains("no fit"));
    }

    #[test]
    fn empty_observations_returns_message() {
        let mut session = Session::new();
        session.last_fit = Some(crate::solver::FitResult {
            coefficients: vec![1.0],
            sigma: vec![0.1],
            sky_rms: 5.0,
            term_names: vec!["IH".to_string()],
        });
        let result = Gdist.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => assert_eq!(s, "No active observations"),
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn terminal_output_contains_both_distributions() {
        let mut session = build_session_with_obs();
        let result = Gdist.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(s.contains("dX Distribution"), "missing dX Distribution");
                assert!(s.contains("dDec Distribution"), "missing dDec Distribution");
            }
            _ => panic!("expected Text output"),
        }
    }

    fn build_session_with_obs() -> Session {
        use celestial_core::Angle;
        use crate::observation::{Observation, PierSide};

        let mut session = Session::new();
        session.last_fit = Some(crate::solver::FitResult {
            coefficients: vec![],
            sigma: vec![],
            sky_rms: 5.0,
            term_names: vec![],
        });
        for i in 0..10 {
            let offset = (i as f64) * 10.0;
            session.observations.push(Observation {
                catalog_ra: Angle::from_hours(0.0),
                catalog_dec: Angle::from_degrees(45.0),
                observed_ra: Angle::from_hours(0.0),
                observed_dec: Angle::from_degrees(45.0 + offset / 3600.0),
                lst: Angle::from_hours(0.0),
                commanded_ha: Angle::from_arcseconds(0.0),
                actual_ha: Angle::from_arcseconds(offset),
                pier_side: PierSide::East,
                masked: false,
            });
        }
        session
    }
}
