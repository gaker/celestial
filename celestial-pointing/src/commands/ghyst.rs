use std::path::Path;

use crate::error::Result;
use crate::plot::residuals::{compute_residuals, require_fit};
use crate::session::Session;

use super::{Command, CommandOutput};

pub struct Ghyst;

impl Command for Ghyst {
    fn name(&self) -> &str {
        "GHYST"
    }

    fn description(&self) -> &str {
        "Hysteresis plot (residuals by sequence and pier side)"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        require_fit(session)?;
        let residuals = compute_residuals(session);
        if residuals.is_empty() {
            return Ok(CommandOutput::Text("No active observations".to_string()));
        }
        let (east, west) = split_by_pier(&residuals);
        let all: Vec<(f64, f64)> = residuals.iter().map(|r| (r.index as f64, r.dr)).collect();
        if let Some(path) = args.first() {
            write_svg(&east, &west, Path::new(path))
        } else {
            terminal_output(&all, east.len(), west.len())
        }
    }
}

type PointVec = Vec<(f64, f64)>;

fn split_by_pier(residuals: &[crate::plot::residuals::ObsResidual]) -> (PointVec, PointVec) {
    let east = residuals
        .iter()
        .filter(|r| r.pier_east)
        .map(|r| (r.index as f64, r.dr))
        .collect();
    let west = residuals
        .iter()
        .filter(|r| !r.pier_east)
        .map(|r| (r.index as f64, r.dr))
        .collect();
    (east, west)
}

fn terminal_output(
    all: &[(f64, f64)],
    n_east: usize,
    n_west: usize,
) -> Result<CommandOutput> {
    let plot = crate::plot::terminal::xy_plot_terminal(
        all,
        "Residual vs Observation Sequence",
        "Obs #",
        "dR (arcsec)",
    );
    let summary = format!("  East: {} obs  West: {} obs", n_east, n_west);
    Ok(CommandOutput::Text(format!("{plot}\n{summary}")))
}

fn write_svg(
    east: &[(f64, f64)],
    west: &[(f64, f64)],
    path: &Path,
) -> Result<CommandOutput> {
    let stem = path.file_stem().unwrap_or_default().to_str().unwrap_or("plot");
    let ext = path.extension().unwrap_or_default().to_str().unwrap_or("svg");
    let parent = path.parent().unwrap_or(Path::new("."));
    let east_path = parent.join(format!("{stem}_east.{ext}"));
    let west_path = parent.join(format!("{stem}_west.{ext}"));
    if !east.is_empty() {
        crate::plot::svg::scatter_svg(
            east,
            &east_path,
            "Hysteresis - East",
            "Obs #",
            "dR (arcsec)",
        )
        .map_err(svg_err)?;
    }
    if !west.is_empty() {
        crate::plot::svg::scatter_svg(
            west,
            &west_path,
            "Hysteresis - West",
            "Obs #",
            "dR (arcsec)",
        )
        .map_err(svg_err)?;
    }
    Ok(CommandOutput::Text(format!(
        "Written to {} and {}",
        east_path.display(),
        west_path.display()
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
        pier: PierSide,
    ) -> Observation {
        Observation {
            catalog_ra: Angle::from_hours(0.0),
            catalog_dec: Angle::from_degrees(cat_dec_deg),
            observed_ra: Angle::from_hours(0.0),
            observed_dec: Angle::from_degrees(obs_dec_deg),
            lst: Angle::from_hours(0.0),
            commanded_ha: Angle::from_arcseconds(cmd_ha_arcsec),
            actual_ha: Angle::from_arcseconds(act_ha_arcsec),
            pier_side: pier,
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
        let result = Ghyst.execute(&mut session, &[]);
        assert!(result.is_err());
    }

    #[test]
    fn empty_observations_returns_message() {
        let mut session = session_with_fit();
        let result = Ghyst.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => assert_eq!(s, "No active observations"),
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn pier_side_splitting() {
        let mut session = session_with_fit();
        session
            .observations
            .push(make_obs(0.0, 100.0, 45.0, 45.01, PierSide::East));
        session
            .observations
            .push(make_obs(0.0, -50.0, 30.0, 30.005, PierSide::West));
        session
            .observations
            .push(make_obs(0.0, 200.0, 60.0, 60.02, PierSide::East));
        let residuals = compute_residuals(&session);
        let (east, west) = split_by_pier(&residuals);
        assert_eq!(east.len(), 2);
        assert_eq!(west.len(), 1);
    }

    #[test]
    fn terminal_output_contains_summary() {
        let mut session = session_with_fit();
        session
            .observations
            .push(make_obs(0.0, 100.0, 45.0, 45.01, PierSide::East));
        session
            .observations
            .push(make_obs(0.0, -50.0, 30.0, 30.005, PierSide::West));
        let result = Ghyst.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(s.contains("Residual vs Observation Sequence"));
                assert!(s.contains("East: 1 obs"));
                assert!(s.contains("West: 1 obs"));
            }
            _ => panic!("expected Text output"),
        }
    }

    #[test]
    fn svg_writes_both_files() {
        let mut session = session_with_fit();
        session
            .observations
            .push(make_obs(0.0, 100.0, 45.0, 45.01, PierSide::East));
        session
            .observations
            .push(make_obs(0.0, 200.0, 50.0, 50.02, PierSide::East));
        session
            .observations
            .push(make_obs(0.0, -50.0, 30.0, 30.005, PierSide::West));
        session
            .observations
            .push(make_obs(0.0, -80.0, 35.0, 35.008, PierSide::West));
        let dir = std::env::temp_dir();
        let path = dir.join("ghyst_test.svg");
        let path_str = path.to_str().unwrap();
        let result = Ghyst.execute(&mut session, &[path_str]).unwrap();
        let east_path = dir.join("ghyst_test_east.svg");
        let west_path = dir.join("ghyst_test_west.svg");
        match &result {
            CommandOutput::Text(s) => {
                assert!(s.contains("Written to"));
                assert!(s.contains("east"));
                assert!(s.contains("west"));
            }
            _ => panic!("expected Text output"),
        }
        assert!(east_path.exists());
        assert!(west_path.exists());
        let east_svg = std::fs::read_to_string(&east_path).unwrap();
        let west_svg = std::fs::read_to_string(&west_path).unwrap();
        assert!(east_svg.contains("<svg"));
        assert!(west_svg.contains("<svg"));
        std::fs::remove_file(&east_path).ok();
        std::fs::remove_file(&west_path).ok();
    }
}
