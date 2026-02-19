use std::path::Path;

use crate::error::Result;
use crate::plot::residuals::{compute_residuals, require_fit};
use crate::session::Session;

use super::{Command, CommandOutput};

pub struct Gdec;

impl Command for Gdec {
    fn name(&self) -> &str {
        "GDEC"
    }
    fn description(&self) -> &str {
        "Residuals vs declination"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        require_fit(session)?;
        let residuals = compute_residuals(session);
        if residuals.is_empty() {
            return Ok(CommandOutput::Text("No active observations".to_string()));
        }
        let dx_vs_dec: Vec<(f64, f64)> = residuals.iter().map(|r| (r.dec_deg, r.dx)).collect();
        let dd_vs_dec: Vec<(f64, f64)> = residuals.iter().map(|r| (r.dec_deg, r.dd)).collect();
        if let Some(path) = args.first() {
            write_svg(&dx_vs_dec, &dd_vs_dec, Path::new(path))
        } else {
            terminal_output(&dx_vs_dec, &dd_vs_dec)
        }
    }
}

fn terminal_output(dx_vs_dec: &[(f64, f64)], dd_vs_dec: &[(f64, f64)]) -> Result<CommandOutput> {
    let dx_plot = crate::plot::terminal::xy_plot_terminal(
        dx_vs_dec,
        "dX vs Declination",
        "Dec (deg)",
        "dX (arcsec)",
    );
    let dd_plot = crate::plot::terminal::xy_plot_terminal(
        dd_vs_dec,
        "dDec vs Declination",
        "Dec (deg)",
        "dDec (arcsec)",
    );
    Ok(CommandOutput::Text(format!("{dx_plot}\n{dd_plot}")))
}

fn write_svg(
    dx_vs_dec: &[(f64, f64)],
    dd_vs_dec: &[(f64, f64)],
    path: &Path,
) -> Result<CommandOutput> {
    let stem = path
        .file_stem()
        .unwrap_or_default()
        .to_str()
        .unwrap_or("plot");
    let ext = path
        .extension()
        .unwrap_or_default()
        .to_str()
        .unwrap_or("svg");
    let parent = path.parent().unwrap_or(Path::new("."));
    let dx_path = parent.join(format!("{stem}_dx.{ext}"));
    let dd_path = parent.join(format!("{stem}_dd.{ext}"));
    crate::plot::svg::scatter_svg(
        dx_vs_dec,
        &dx_path,
        "dX vs Declination",
        "Dec (deg)",
        "dX (arcsec)",
    )
    .map_err(svg_err)?;
    crate::plot::svg::scatter_svg(
        dd_vs_dec,
        &dd_path,
        "dDec vs Declination",
        "Dec (deg)",
        "dDec (arcsec)",
    )
    .map_err(svg_err)?;
    Ok(CommandOutput::Text(format!(
        "Written to {} and {}",
        dx_path.display(),
        dd_path.display()
    )))
}

fn svg_err(e: Box<dyn std::error::Error>) -> crate::error::Error {
    crate::error::Error::Io(std::io::Error::other(e.to_string()))
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
        let result = Gdec.execute(&mut session, &[]);
        assert!(result.is_err());
    }

    #[test]
    fn empty_observations_returns_message() {
        let mut session = session_with_fit();
        let result = Gdec.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => assert_eq!(s, "No active observations"),
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn terminal_shows_both_dx_and_ddec() {
        let mut session = session_with_fit();
        session.observations.push(make_obs(0.0, 100.0, 45.0, 45.01));
        session
            .observations
            .push(make_obs(0.0, -50.0, 30.0, 30.005));
        let result = Gdec.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(s.contains("dX vs Declination"), "missing dX vs Declination");
                assert!(
                    s.contains("dDec vs Declination"),
                    "missing dDec vs Declination"
                );
                assert!(s.contains("Dec (deg)"), "missing Dec (deg) label");
            }
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn svg_writes_two_files() {
        let mut session = session_with_fit();
        session.observations.push(make_obs(0.0, 100.0, 45.0, 45.01));
        session
            .observations
            .push(make_obs(0.0, -50.0, 30.0, 30.005));
        let dir = std::env::temp_dir();
        let path = dir.join("gdec_test.svg");
        let path_str = path.to_str().unwrap();
        let result = Gdec.execute(&mut session, &[path_str]).unwrap();
        let dx_path = dir.join("gdec_test_dx.svg");
        let dd_path = dir.join("gdec_test_dd.svg");
        match &result {
            CommandOutput::Text(s) => {
                assert!(s.contains("Written to"), "missing Written to");
                assert!(s.contains("_dx.svg"), "missing _dx.svg");
                assert!(s.contains("_dd.svg"), "missing _dd.svg");
            }
            _ => panic!("expected Text output"),
        }
        assert!(dx_path.exists());
        assert!(dd_path.exists());
        let dx_contents = std::fs::read_to_string(&dx_path).unwrap();
        let dd_contents = std::fs::read_to_string(&dd_path).unwrap();
        assert!(dx_contents.contains("<svg"));
        assert!(dd_contents.contains("<svg"));
        std::fs::remove_file(&dx_path).ok();
        std::fs::remove_file(&dd_path).ok();
    }
}
