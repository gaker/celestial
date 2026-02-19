use crate::error::{Error, Result};
use crate::model::PointingModel;
use crate::observation::{IndatFile, MountType, Observation, SiteParams};
use crate::solver::{self, FitResult};
use celestial_core::Angle;
use celestial_time::JulianDate;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AdjustDirection {
    #[default]
    TelescopeToStar,
    StarToTelescope,
}

pub struct Session {
    pub observations: Vec<Observation>,
    pub model: PointingModel,
    pub site: Option<SiteParams>,
    pub mount_type: MountType,
    pub last_fit: Option<FitResult>,
    pub header_lines: Vec<String>,
    pub date: Option<JulianDate>,
    pub adjust_direction: AdjustDirection,
    pub lst_override: Option<Angle>,
}

impl Default for Session {
    fn default() -> Self {
        Self {
            observations: Vec::new(),
            model: PointingModel::new(),
            site: None,
            mount_type: MountType::GermanEquatorial,
            last_fit: None,
            header_lines: Vec::new(),
            date: None,
            adjust_direction: AdjustDirection::default(),
            lst_override: None,
        }
    }
}

impl Session {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn load_indat(&mut self, indat: IndatFile) {
        self.observations = indat.observations;
        self.site = Some(indat.site);
        self.mount_type = indat.mount_type;
        self.header_lines = indat.header_lines;
        self.date = Some(indat.date);
        self.last_fit = None;
    }

    pub fn fit(&mut self) -> Result<&FitResult> {
        let lat = self.latitude();
        let active: Vec<&Observation> = self.observations.iter().filter(|o| !o.masked).collect();
        let fixed = self.model.fixed_flags();
        let coefficients = self.model.coefficients();
        let result = solver::fit_model(&active, self.model.terms(), fixed, coefficients, lat)?;
        self.model.set_coefficients(&result.coefficients)?;
        self.last_fit = Some(result);
        Ok(self.last_fit.as_ref().unwrap())
    }

    pub fn active_observation_count(&self) -> usize {
        self.observations.iter().filter(|o| !o.masked).count()
    }

    pub fn masked_observation_count(&self) -> usize {
        self.observations.iter().filter(|o| o.masked).count()
    }

    pub fn observation_count(&self) -> usize {
        self.observations.len()
    }

    pub fn current_lst(&self) -> Result<Angle> {
        if let Some(lst) = self.lst_override {
            return Ok(lst);
        }
        Err(Error::NoLst)
    }

    pub fn latitude(&self) -> f64 {
        self.site.as_ref().map_or(0.0, |s| s.latitude.radians())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observation::PierSide;
    use crate::parser::parse_indat;
    use celestial_core::Angle;

    #[test]
    fn new_session_defaults() {
        let session = Session::new();
        assert_eq!(session.observation_count(), 0);
        assert_eq!(session.mount_type, MountType::GermanEquatorial);
        assert!(session.site.is_none());
        assert!(session.last_fit.is_none());
        assert!(session.header_lines.is_empty());
        assert!(session.date.is_none());
        assert_eq!(session.model.term_count(), 0);
    }

    #[test]
    fn latitude_no_site_returns_zero() {
        let session = Session::new();
        assert_eq!(session.latitude(), 0.0);
    }

    #[test]
    fn load_indat_populates_session() {
        let content = "\
ASCOM Mount
:NODA
:EQUAT
+39 00 26 2024 7 14 29.20 987.00 231.65  0.94 0.5500 0.0065
21 43 18.4460 +72 29 08.368 09 28 59.9527 +109 20 06.469  16 23.130
23 46 02.2988 +77 38 38.725 11 26 17.6308 +104 03 28.734  16 24.711";

        let indat = parse_indat(content).unwrap();
        let mut session = Session::new();
        session.load_indat(indat);

        assert_eq!(session.observation_count(), 2);
        assert_eq!(session.mount_type, MountType::GermanEquatorial);
        assert!(session.site.is_some());
        assert!(session.date.is_some());
        assert_eq!(session.header_lines.len(), 1);
        assert!(session.last_fit.is_none());
    }

    #[test]
    fn load_indat_sets_latitude() {
        let content = "\
!Test
:EQUAT
+39 00 26 2024 7 14 29.20 987.00 231.65 0.94 0.5500 0.0065
21 43 18.4460 +72 29 08.368 09 28 59.9527 +109 20 06.469  16 23.130";

        let indat = parse_indat(content).unwrap();
        let mut session = Session::new();
        session.load_indat(indat);

        let expected = Angle::from_degrees(39.0 + 26.0 / 3600.0).radians();
        assert_eq!(session.latitude(), expected);
    }

    #[test]
    fn load_indat_clears_previous_fit() {
        let content = "\
!Test
:EQUAT
+39 00 26 2024 7 14 29.20 987.00 231.65 0.94 0.5500 0.0065
21 43 18.4460 +72 29 08.368 09 28 59.9527 +109 20 06.469  16 23.130";

        let indat1 = parse_indat(content).unwrap();
        let indat2 = parse_indat(content).unwrap();
        let mut session = Session::new();
        session.load_indat(indat1);
        session.last_fit = Some(FitResult {
            coefficients: vec![1.0],
            sigma: vec![0.1],
            sky_rms: 5.0,
            term_names: vec!["IH".to_string()],
        });
        session.load_indat(indat2);
        assert!(session.last_fit.is_none());
    }

    fn make_obs(cmd_ha_arcsec: f64, act_ha_arcsec: f64, dec_deg: f64) -> Observation {
        Observation {
            catalog_ra: Angle::from_hours(0.0),
            catalog_dec: Angle::from_degrees(dec_deg),
            observed_ra: Angle::from_hours(0.0),
            observed_dec: Angle::from_degrees(dec_deg),
            lst: Angle::from_hours(0.0),
            commanded_ha: Angle::from_arcseconds(cmd_ha_arcsec),
            actual_ha: Angle::from_arcseconds(act_ha_arcsec),
            pier_side: PierSide::East,
            masked: false,
        }
    }

    #[test]
    fn fit_updates_model_and_stores_result() {
        let mut session = Session::new();
        session.observations = vec![
            make_obs(0.0, 100.0, 30.0),
            make_obs(0.0, 100.0, 45.0),
            make_obs(0.0, 100.0, 60.0),
        ];
        session.model.add_term("IH").unwrap();
        let result = session.fit().unwrap();
        assert_eq!(result.term_names, vec!["IH"]);
        assert!((result.coefficients[0] - (-100.0)).abs() < 1e-6);
        assert!(session.last_fit.is_some());
        assert_eq!(session.model.coefficients().len(), 1);
    }

    #[test]
    fn fit_no_terms_returns_error() {
        let mut session = Session::new();
        session.observations = vec![make_obs(0.0, 100.0, 30.0)];
        let result = session.fit();
        assert!(result.is_err());
    }

    #[test]
    fn fit_no_observations_returns_error() {
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        let result = session.fit();
        assert!(result.is_err());
    }
}
