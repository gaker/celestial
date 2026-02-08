use std::path::Path;

use crate::error::Result;
use crate::plot::residuals::{compute_residuals, require_fit};
use crate::session::Session;

use super::{Command, CommandOutput};

pub struct Gscat;

impl Command for Gscat {
    fn name(&self) -> &str {
        "GSCAT"
    }

    fn description(&self) -> &str {
        "Scatter plot of residuals (dX vs dDec)"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        require_fit(session)?;
        let residuals = compute_residuals(session);
        if residuals.is_empty() {
            return Ok(CommandOutput::Text("No active observations".to_string()));
        }
        let points: Vec<(f64, f64)> = residuals.iter().map(|r| (r.dx, r.dd)).collect();
        if let Some(path) = args.first() {
            write_svg(&points, Path::new(path))
        } else {
            terminal_output(&points)
        }
    }
}

fn write_svg(points: &[(f64, f64)], path: &Path) -> Result<CommandOutput> {
    crate::plot::svg::scatter_svg(
        points,
        path,
        "Residual Scatter",
        "dX (arcsec)",
        "dDec (arcsec)",
    )
    .map_err(|e| crate::error::Error::Io(std::io::Error::other(e.to_string())))?;
    Ok(CommandOutput::Text(format!("Wrote {}", path.display())))
}

fn terminal_output(points: &[(f64, f64)]) -> Result<CommandOutput> {
    let text = crate::plot::terminal::scatter_terminal(
        points,
        "Residual Scatter (dX vs dDec)",
        "dX (arcsec)",
        "dDec (arcsec)",
    );
    Ok(CommandOutput::Text(text))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observation::{Observation, PierSide};
    use crate::solver::FitResult;
    use celestial_core::Angle;

    fn make_obs(
        cmd_ha_arcsec: f64,
        act_ha_arcsec: f64,
        cat_dec_deg: f64,
        obs_dec_deg: f64,
    ) -> Observation {
        Observation {
            catalog_ra: Angle::from_hours(0.0),
            catalog_dec: Angle::from_degrees(cat_dec_deg),
            observed_ra: Angle::from_hours(0.0),
            observed_dec: Angle::from_degrees(obs_dec_deg),
            lst: Angle::from_hours(0.0),
            commanded_ha: Angle::from_arcseconds(cmd_ha_arcsec),
            actual_ha: Angle::from_arcseconds(act_ha_arcsec),
            pier_side: PierSide::East,
            masked: false,
        }
    }

    fn session_with_fit() -> Session {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        session.model.set_coefficients(&[0.0]).unwrap();
        session.last_fit = Some(FitResult {
            coefficients: vec![0.0],
            sigma: vec![0.1],
            sky_rms: 1.0,
            term_names: vec!["IH".to_string()],
        });
        session
    }

    #[test]
    fn no_fit_returns_error() {
        let mut session = Session::new();
        let result = Gscat.execute(&mut session, &[]);
        assert!(result.is_err());
    }

    #[test]
    fn empty_observations_returns_message() {
        let mut session = session_with_fit();
        let result = Gscat.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => assert!(s.contains("No active observations")),
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn terminal_output_contains_title() {
        let mut session = session_with_fit();
        session.observations.push(make_obs(0.0, 100.0, 45.0, 45.01));
        session
            .observations
            .push(make_obs(0.0, -50.0, 30.0, 30.005));
        let result = Gscat.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(s.contains("Residual Scatter"));
                assert!(s.contains("dX (arcsec)"));
                assert!(s.contains("dDec (arcsec)"));
            }
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn svg_writes_to_temp_file() {
        let mut session = session_with_fit();
        session.observations.push(make_obs(0.0, 100.0, 45.0, 45.01));
        session
            .observations
            .push(make_obs(0.0, -50.0, 30.0, 30.005));
        let dir = std::env::temp_dir();
        let path = dir.join("gscat_test.svg");
        let path_str = path.to_str().unwrap();
        let result = Gscat.execute(&mut session, &[path_str]).unwrap();
        match &result {
            CommandOutput::Text(s) => assert!(s.contains("Wrote")),
            _ => panic!("expected Text output"),
        }
        assert!(path.exists());
        let contents = std::fs::read_to_string(&path).unwrap();
        assert!(contents.contains("<svg"));
        std::fs::remove_file(&path).ok();
    }
}
