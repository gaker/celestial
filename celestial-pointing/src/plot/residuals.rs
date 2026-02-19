use crate::observation::PierSide;
use crate::session::Session;

pub struct ObsResidual {
    pub ha_deg: f64,
    pub dec_deg: f64,
    pub dh: f64,
    pub dd: f64,
    pub dx: f64,
    pub dr: f64,
    pub index: usize,
    pub pier_east: bool,
}

pub fn compute_residuals(session: &Session) -> Vec<ObsResidual> {
    let lat = session.latitude();
    session
        .observations
        .iter()
        .enumerate()
        .filter(|(_, obs)| !obs.masked)
        .map(|(i, obs)| build_residual(i, obs, &session.model, lat))
        .collect()
}

fn build_residual(
    index: usize,
    obs: &crate::observation::Observation,
    model: &crate::model::PointingModel,
    lat: f64,
) -> ObsResidual {
    let h = obs.commanded_ha.radians();
    let dec = obs.catalog_dec.radians();
    let pier = obs.pier_side.sign();
    let (model_dh, model_dd) = model.apply_equatorial(h, dec, lat, pier);
    let raw_dh = (obs.actual_ha - obs.commanded_ha).arcseconds();
    let raw_dd = (obs.observed_dec - obs.catalog_dec).arcseconds();
    let dh = raw_dh - model_dh;
    let dd = raw_dd - model_dd;
    let dx = dh * libm::cos(dec);
    let dr = libm::sqrt(dx * dx + dd * dd);
    ObsResidual {
        ha_deg: obs.commanded_ha.degrees(),
        dec_deg: obs.catalog_dec.degrees(),
        dh,
        dd,
        dx,
        dr,
        index,
        pier_east: obs.pier_side == PierSide::East,
    }
}

pub fn require_fit(session: &Session) -> crate::error::Result<()> {
    if session.last_fit.is_none() {
        return Err(crate::error::Error::Fit(
            "no fit results - run FIT first".into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observation::{Observation, PierSide};
    use celestial_core::Angle;

    fn make_obs(
        cmd_ha_arcsec: f64,
        act_ha_arcsec: f64,
        cat_dec_deg: f64,
        obs_dec_deg: f64,
        pier: PierSide,
        masked: bool,
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
            masked,
        }
    }

    #[test]
    fn empty_session_returns_empty() {
        let session = Session::new();
        let residuals = compute_residuals(&session);
        assert!(residuals.is_empty());
    }

    #[test]
    fn masked_observations_excluded() {
        let mut session = Session::new();
        session
            .observations
            .push(make_obs(0.0, 100.0, 45.0, 45.0, PierSide::East, false));
        session
            .observations
            .push(make_obs(0.0, 200.0, 30.0, 30.0, PierSide::East, true));
        let residuals = compute_residuals(&session);
        assert_eq!(residuals.len(), 1);
        assert_eq!(residuals[0].index, 0);
    }

    #[test]
    fn residual_no_model_equals_raw() {
        let mut session = Session::new();
        session
            .observations
            .push(make_obs(0.0, 3600.0, 0.0, 2.0, PierSide::East, false));
        let residuals = compute_residuals(&session);
        assert_eq!(residuals.len(), 1);
        let r = &residuals[0];
        assert_eq!(r.dh, 3600.0);
        assert_eq!(r.dd, 7200.0);
        assert_eq!(r.dx, 3600.0 * libm::cos(0.0_f64));
        assert_eq!(r.dr, libm::sqrt(3600.0_f64.powi(2) + 7200.0_f64.powi(2)));
    }

    #[test]
    fn residual_with_ih_model() {
        let mut session = Session::new();
        session
            .observations
            .push(make_obs(0.0, 100.0, 45.0, 45.0, PierSide::East, false));
        session.model.add_term("IH").unwrap();
        session.model.set_coefficients(&[-100.0]).unwrap();
        let residuals = compute_residuals(&session);
        let r = &residuals[0];
        let dec_rad = 45.0_f64.to_radians();
        let model_dh = 100.0;
        let expected_dh = 100.0 - model_dh;
        assert_eq!(r.dh, expected_dh);
        assert_eq!(r.dx, expected_dh * libm::cos(dec_rad));
    }

    #[test]
    fn pier_side_recorded() {
        let mut session = Session::new();
        session
            .observations
            .push(make_obs(0.0, 0.0, 0.0, 0.0, PierSide::East, false));
        session
            .observations
            .push(make_obs(0.0, 0.0, 0.0, 0.0, PierSide::West, false));
        let residuals = compute_residuals(&session);
        assert!(residuals[0].pier_east);
        assert!(!residuals[1].pier_east);
    }

    #[test]
    fn index_tracks_original_position() {
        let mut session = Session::new();
        session
            .observations
            .push(make_obs(0.0, 0.0, 0.0, 0.0, PierSide::East, true));
        session
            .observations
            .push(make_obs(0.0, 0.0, 0.0, 0.0, PierSide::East, false));
        session
            .observations
            .push(make_obs(0.0, 0.0, 0.0, 0.0, PierSide::East, false));
        let residuals = compute_residuals(&session);
        assert_eq!(residuals.len(), 2);
        assert_eq!(residuals[0].index, 1);
        assert_eq!(residuals[1].index, 2);
    }

    #[test]
    fn require_fit_no_fit() {
        let session = Session::new();
        assert!(require_fit(&session).is_err());
    }

    #[test]
    fn require_fit_with_fit() {
        let mut session = Session::new();
        session.last_fit = Some(crate::solver::FitResult {
            coefficients: vec![1.0],
            sigma: vec![0.1],
            sky_rms: 5.0,
            term_names: vec!["IH".to_string()],
        });
        assert!(require_fit(&session).is_ok());
    }

    #[test]
    fn ha_and_dec_degrees_populated() {
        let mut session = Session::new();
        let cmd_ha_arcsec = 3600.0 * 15.0;
        session.observations.push(make_obs(
            cmd_ha_arcsec,
            cmd_ha_arcsec,
            45.0,
            45.0,
            PierSide::East,
            false,
        ));
        let residuals = compute_residuals(&session);
        let r = &residuals[0];
        assert_eq!(r.ha_deg, Angle::from_arcseconds(cmd_ha_arcsec).degrees());
        assert_eq!(r.dec_deg, 45.0);
    }

    #[test]
    fn dr_is_sqrt_dx2_dd2() {
        let mut session = Session::new();
        session
            .observations
            .push(make_obs(0.0, 3.0, 0.0, 4.0 / 3600.0, PierSide::East, false));
        let residuals = compute_residuals(&session);
        let r = &residuals[0];
        let expected_dr = libm::sqrt(r.dx * r.dx + r.dd * r.dd);
        assert_eq!(r.dr, expected_dr);
    }
}
