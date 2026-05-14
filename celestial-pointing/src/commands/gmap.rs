use std::path::Path;

use crate::error::Result;
use crate::plot::residuals::{compute_residuals, require_fit};
use crate::session::Session;

use super::{Command, CommandOutput};

pub struct Gmap;

impl Command for Gmap {
    fn name(&self) -> &str {
        "GMAP"
    }

    fn description(&self) -> &str {
        "Sky map with residual vectors"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        require_fit(session)?;
        let residuals = compute_residuals(session);
        if residuals.is_empty() {
            return Ok(CommandOutput::Text("No active observations".to_string()));
        }
        let positions: Vec<(f64, f64)> = residuals.iter().map(|r| (r.ha_deg, r.dec_deg)).collect();
        let vectors: Vec<(f64, f64)> = residuals
            .iter()
            .map(|r| (r.dx / 3600.0, r.dd / 3600.0))
            .collect();

        if let Some(path) = args.first() {
            let scale = parse_scale(args);
            write_svg(&positions, &vectors, Path::new(path), scale)
        } else {
            terminal_output(&positions)
        }
    }
}

fn parse_scale(args: &[&str]) -> f64 {
    args.get(1)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(10.0)
}

fn terminal_output(positions: &[(f64, f64)]) -> Result<CommandOutput> {
    let text =
        crate::plot::terminal::scatter_terminal(positions, "Sky Map", "HA (deg)", "Dec (deg)");
    Ok(CommandOutput::Text(text))
}

fn write_svg(
    positions: &[(f64, f64)],
    vectors: &[(f64, f64)],
    path: &Path,
    scale: f64,
) -> Result<CommandOutput> {
    crate::plot::svg::vector_map_svg(
        positions,
        vectors,
        path,
        "Sky Map - Residual Vectors",
        "HA (deg)",
        "Dec (deg)",
        scale,
    )
    .map_err(|e| crate::error::Error::Io(std::io::Error::other(e.to_string())))?;
    Ok(CommandOutput::Text(format!(
        "Written to {}",
        path.display()
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observation::Observation;
    use crate::test_support::{FitResultBuilder, ObsBuilder};

    fn make_obs(
        cmd_ha_arcsec: f64,
        act_ha_arcsec: f64,
        cat_dec_deg: f64,
        obs_dec_deg: f64,
    ) -> Observation {
        ObsBuilder::new()
            .commanded_ha_arcsec(cmd_ha_arcsec)
            .actual_ha_arcsec(act_ha_arcsec)
            .catalog_dec_deg(cat_dec_deg)
            .observed_dec_deg(obs_dec_deg)
            .build()
    }

    fn session_with_fit() -> Session {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        session.model.set_coefficients(&[0.0]).unwrap();
        session.last_fit = Some(
            FitResultBuilder::new()
                .coefficients(vec![0.0])
                .sky_rms(1.0)
                .popn_sd(1.0)
                .term_names(["IH"])
                .build(),
        );
        session
    }

    #[test]
    fn no_fit_returns_error() {
        let mut session = Session::new();
        let result = Gmap.execute(&mut session, &[]);
        assert!(result.is_err());
    }

    #[test]
    fn empty_observations_returns_message() {
        let mut session = session_with_fit();
        let result = Gmap.execute(&mut session, &[]).unwrap();
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
        let result = Gmap.execute(&mut session, &[]).unwrap();
        match result {
            CommandOutput::Text(s) => {
                assert!(s.contains("Sky Map"));
                assert!(s.contains("HA (deg)"));
                assert!(s.contains("Dec (deg)"));
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
        let path = dir.join("gmap_test.svg");
        let path_str = path.to_str().unwrap();
        let result = Gmap.execute(&mut session, &[path_str]).unwrap();
        match &result {
            CommandOutput::Text(s) => assert!(s.contains("Written to")),
            _ => panic!("expected Text output"),
        }
        assert!(path.exists());
        let contents = std::fs::read_to_string(&path).unwrap();
        assert!(contents.contains("<svg"));
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn scale_parsed_from_args() {
        assert_eq!(parse_scale(&["out.svg"]), 10.0);
        assert_eq!(parse_scale(&["out.svg", "5.0"]), 5.0);
        assert_eq!(parse_scale(&["out.svg", "notanumber"]), 10.0);
        assert_eq!(parse_scale(&[]), 10.0);
    }

    #[test]
    fn svg_with_custom_scale() {
        let mut session = session_with_fit();
        session.observations.push(make_obs(0.0, 100.0, 45.0, 45.01));
        let dir = std::env::temp_dir();
        let path = dir.join("gmap_scale_test.svg");
        let path_str = path.to_str().unwrap();
        let result = Gmap.execute(&mut session, &[path_str, "20.0"]).unwrap();
        match &result {
            CommandOutput::Text(s) => assert!(s.contains("Written to")),
            _ => panic!("expected Text output"),
        }
        assert!(path.exists());
        std::fs::remove_file(&path).ok();
    }
}
