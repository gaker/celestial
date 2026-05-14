use crate::error::{Error, Result};
use crate::model::PointingModel;
use crate::observation::{IndatFile, IndatOption, MountType, Observation, SiteParams};
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
    pub options: Vec<IndatOption>,
    pub fit_tol: f64,
}

pub const DEFAULT_FIT_TOL: f64 = 1.0e-3;

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
            options: Vec::new(),
            fit_tol: DEFAULT_FIT_TOL,
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
        self.options = indat.options;
        self.last_fit = None;
    }

    pub fn prepared_observations(&self) -> Vec<Observation> {
        self.prepare_observations()
    }

    pub fn fit(&mut self) -> Result<&FitResult> {
        let lat = self.latitude();
        let prepared = self.prepare_observations();
        let active: Vec<&Observation> = prepared.iter().filter(|o| !o.masked).collect();
        let fixed = self.model.fixed_flags();
        let parallel = self.model.parallel_flags();
        let coefficients = self.model.coefficients();
        let options = solver::IterOptions::default();
        let inputs = solver::FitInputs {
            observations: &active,
            terms: self.model.terms(),
            fixed,
            parallel,
            coefficients,
            latitude: lat,
            fit_tol: self.fit_tol,
        };
        let result = solver::fit_model_iterated(&inputs, &options)?;
        self.model.set_coefficients(&result.coefficients)?;
        self.last_fit = Some(result);
        Ok(self.last_fit.as_ref().unwrap())
    }

    fn prepare_observations(&self) -> Vec<Observation> {
        let Some(site) = self.site.as_ref() else {
            return self.observations.clone();
        };
        let apply_diurnal = !self.options.contains(&IndatOption::NoDA);
        let apply_refraction = site.pressure > 0.0;
        if !apply_diurnal && !apply_refraction {
            return self.observations.clone();
        }
        let lat = site.latitude;
        let lon = site.longitude;
        let height_m = site.elevation;
        let location = match celestial_core::Location::new(
            lat.radians(),
            lon.radians(),
            height_m,
        ) {
            Ok(l) => l,
            Err(_) => return self.observations.clone(),
        };
        self.observations
            .iter()
            .map(|obs| {
                let mut cat_ra = obs.catalog_ra;
                let mut cat_dec = obs.catalog_dec;
                if apply_diurnal {
                    let (ra, dec) = crate::diurnal::apply_diurnal(
                        cat_ra, cat_dec, obs.lst, lat, lon, height_m,
                    );
                    cat_ra = ra;
                    cat_dec = dec;
                }
                if apply_refraction {
                    let (ra, dec) = crate::prepare::apply_refraction(
                        cat_ra, cat_dec, obs.lst, &location, site,
                    );
                    cat_ra = ra;
                    cat_dec = dec;
                }
                let commanded_ha = (obs.lst - cat_ra).wrapped();
                Observation {
                    catalog_ra: cat_ra,
                    catalog_dec: cat_dec,
                    commanded_ha,
                    ..obs.clone()
                }
            })
            .collect()
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
    use crate::test_support::{obs_with_ha_offset as make_obs, FitResultBuilder, ObsBuilder};
    use celestial_time::JulianDate;

    fn site(lat_deg: f64, pressure: f64) -> SiteParams {
        SiteParams {
            latitude: Angle::from_degrees(lat_deg),
            longitude: Angle::from_degrees(0.0),
            temperature: 10.0,
            pressure,
            elevation: 0.0,
            humidity: 0.5,
            wavelength: 0.55,
            lapse_rate: 0.0065,
        }
    }

    fn indat(observations: Vec<Observation>, options: Vec<IndatOption>) -> IndatFile {
        IndatFile {
            site: site(39.0, 987.0),
            options,
            observations,
            mount_type: MountType::Altazimuth,
            header_lines: vec!["!hdr".into()],
            date: JulianDate::from_calendar(2024, 7, 14, 0, 0, 0.0),
        }
    }

    fn fake_fit(coefs: Vec<f64>) -> FitResult {
        FitResultBuilder::new()
            .coefficients(coefs)
            .sky_rms(1.0)
            .popn_sd(1.0)
            .build()
    }

    // --- constructor + defaults ---------------------------------------

    #[test]
    fn new_matches_default() {
        let a = Session::new();
        let b = Session::default();
        assert_eq!(a.observation_count(), b.observation_count());
        assert_eq!(a.mount_type, b.mount_type);
        assert_eq!(a.fit_tol, b.fit_tol);
    }

    #[test]
    fn default_session_is_empty_with_canonical_defaults() {
        let s = Session::new();
        assert_eq!(s.observation_count(), 0);
        assert_eq!(s.active_observation_count(), 0);
        assert_eq!(s.masked_observation_count(), 0);
        assert_eq!(s.model.term_count(), 0);
        assert!(s.site.is_none());
        assert!(s.last_fit.is_none());
        assert!(s.header_lines.is_empty());
        assert!(s.date.is_none());
        assert!(s.lst_override.is_none());
        assert!(s.options.is_empty());
        assert_eq!(s.mount_type, MountType::GermanEquatorial);
        assert_eq!(s.adjust_direction, AdjustDirection::TelescopeToStar);
        assert_eq!(s.fit_tol, DEFAULT_FIT_TOL);
    }

    #[test]
    fn default_fit_tol_is_one_milliradian_scale() {
        assert_eq!(DEFAULT_FIT_TOL, 1.0e-3);
    }

    #[test]
    fn adjust_direction_default_is_telescope_to_star() {
        assert_eq!(AdjustDirection::default(), AdjustDirection::TelescopeToStar);
    }

    // --- latitude / current_lst ---------------------------------------

    #[test]
    fn latitude_returns_zero_when_no_site() {
        assert_eq!(Session::new().latitude(), 0.0);
    }

    #[test]
    fn latitude_returns_site_value_in_radians() {
        let mut s = Session::new();
        s.site = Some(site(45.0, 0.0));
        let expected = Angle::from_degrees(45.0).radians();
        assert!((s.latitude() - expected).abs() < 1e-15);
    }

    #[test]
    fn current_lst_returns_no_lst_when_unset() {
        let s = Session::new();
        match s.current_lst() {
            Err(Error::NoLst) => {}
            other => panic!("expected NoLst, got {:?}", other.map(|_| "Ok")),
        }
    }

    #[test]
    fn current_lst_returns_override_when_set() {
        let mut s = Session::new();
        s.lst_override = Some(Angle::from_hours(14.5));
        let got = s.current_lst().unwrap();
        assert!((got.hours() - 14.5).abs() < 1e-12);
    }

    // --- observation counters -----------------------------------------

    #[test]
    fn observation_counts_partition_active_and_masked() {
        let mut s = Session::new();
        s.observations = vec![
            make_obs(0.0, 30.0),
            make_obs(0.0, 30.0),
            make_obs(0.0, 30.0),
        ];
        s.observations[1].masked = true;
        assert_eq!(s.observation_count(), 3);
        assert_eq!(s.active_observation_count(), 2);
        assert_eq!(s.masked_observation_count(), 1);
        assert_eq!(
            s.active_observation_count() + s.masked_observation_count(),
            s.observation_count(),
        );
    }

    #[test]
    fn observation_counts_on_empty_session_are_zero() {
        let s = Session::new();
        assert_eq!(s.observation_count(), 0);
        assert_eq!(s.active_observation_count(), 0);
        assert_eq!(s.masked_observation_count(), 0);
    }

    // --- load_indat ---------------------------------------------------

    #[test]
    fn load_indat_replaces_observations_site_mount_options() {
        let mut s = Session::new();
        // Pre-existing state that should be replaced.
        s.observations = vec![make_obs(0.0, 0.0)];
        s.mount_type = MountType::ForkEquatorial;

        let file = indat(
            vec![make_obs(0.0, 30.0), make_obs(0.0, 45.0)],
            vec![IndatOption::NoDA],
        );
        s.load_indat(file);

        assert_eq!(s.observations.len(), 2);
        assert!(s.site.is_some());
        assert_eq!(s.mount_type, MountType::Altazimuth);
        assert_eq!(s.header_lines, vec!["!hdr"]);
        assert!(s.date.is_some());
        assert_eq!(s.options, vec![IndatOption::NoDA]);
    }

    #[test]
    fn load_indat_clears_previous_last_fit() {
        let mut s = Session::new();
        s.last_fit = Some(fake_fit(vec![1.0]));
        s.load_indat(indat(vec![], vec![]));
        assert!(s.last_fit.is_none());
    }

    // load_indat replaces observations/site/mount/options but does NOT touch
    // the model, lst_override, adjust_direction, or fit_tol. Pinning this
    // explicitly: a user loading new data shouldn't lose their model setup.
    #[test]
    fn load_indat_preserves_model_lst_override_and_user_settings() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        s.lst_override = Some(Angle::from_hours(12.0));
        s.adjust_direction = AdjustDirection::StarToTelescope;
        s.fit_tol = 1.0e-5;

        s.load_indat(indat(vec![], vec![]));

        assert_eq!(s.model.term_names(), vec!["IH"]);
        assert!(s.lst_override.is_some());
        assert_eq!(s.adjust_direction, AdjustDirection::StarToTelescope);
        assert_eq!(s.fit_tol, 1.0e-5);
    }

    // --- prepared_observations ----------------------------------------

    #[test]
    fn prepared_observations_with_no_site_is_unchanged_clone() {
        let mut s = Session::new();
        s.observations = vec![ObsBuilder::new()
            .commanded_ha_arcsec(100.0)
            .actual_ha_arcsec(200.0)
            .catalog_dec_deg(30.0)
            .observed_dec_deg(30.0)
            .build()];
        let prepared = s.prepared_observations();
        assert_eq!(prepared.len(), 1);
        assert_eq!(
            prepared[0].catalog_dec.degrees(),
            s.observations[0].catalog_dec.degrees(),
        );
    }

    // With NoDA *and* pressure=0, both flags fire false → early return clone.
    #[test]
    fn prepared_observations_noda_and_zero_pressure_is_unchanged_clone() {
        let mut s = Session::new();
        s.site = Some(site(39.0, 0.0));
        s.options = vec![IndatOption::NoDA];
        let original = ObsBuilder::new()
            .commanded_ha_arcsec(100.0)
            .actual_ha_arcsec(200.0)
            .catalog_dec_deg(30.0)
            .observed_dec_deg(30.0)
            .build();
        s.observations = vec![original.clone()];
        let prepared = s.prepared_observations();
        assert_eq!(prepared.len(), 1);
        assert_eq!(
            prepared[0].catalog_dec.degrees(),
            original.catalog_dec.degrees(),
        );
        assert_eq!(
            prepared[0].catalog_ra.hours(),
            original.catalog_ra.hours(),
        );
    }

    // Pressure > 0 triggers refraction → catalog coords get nudged.
    #[test]
    fn prepared_observations_with_pressure_applies_refraction() {
        let mut s = Session::new();
        s.site = Some(site(39.0, 987.0));
        s.options = vec![IndatOption::NoDA]; // disable diurnal so only refraction fires
        let target_dec = 45.0_f64;
        s.observations = vec![ObsBuilder::new()
            .lst_hours(6.0)
            .catalog_ra_hours(6.0) // HA = 0, near meridian
            .catalog_dec_deg(target_dec)
            .observed_dec_deg(target_dec)
            .build()];
        let prepared = s.prepared_observations();
        let diff_arcsec =
            (prepared[0].catalog_dec.degrees() - target_dec) * 3600.0;
        assert!(
            diff_arcsec.abs() > 0.1,
            "refraction should shift dec by more than 0.1\"; got {}\"",
            diff_arcsec,
        );
    }

    // Diurnal aberration on by default (no :NODA) at the same site.
    #[test]
    fn prepared_observations_default_options_applies_diurnal() {
        let mut s = Session::new();
        s.site = Some(site(39.0, 0.0)); // pressure 0 → no refraction
        let target_ra = 6.0;
        s.observations = vec![ObsBuilder::new()
            .lst_hours(6.0)
            .catalog_ra_hours(target_ra)
            .build()];
        let prepared = s.prepared_observations();
        let diff_seconds =
            (prepared[0].catalog_ra.hours() - target_ra).abs() * 3600.0;
        assert!(
            diff_seconds < 1.0,
            "diurnal aberration is sub-arcsec; got {} sec of HA",
            diff_seconds,
        );
        // It still ran (commanded_ha was recomputed) — verify by checking
        // commanded_ha matches (lst - cat_ra).wrapped().
        let expected_ha =
            (prepared[0].lst - prepared[0].catalog_ra).wrapped();
        assert!(
            (prepared[0].commanded_ha.arcseconds() - expected_ha.arcseconds()).abs() < 1e-6,
        );
    }

    // commanded_ha gets recomputed in the transform path, overriding whatever
    // was on the original Observation. The non-transform clone path leaves
    // commanded_ha alone — pinning the asymmetry.
    #[test]
    fn prepared_observations_recomputes_commanded_ha_in_transform_path() {
        let mut s = Session::new();
        s.site = Some(site(39.0, 987.0));
        // commanded_ha=999" is bogus; prepared path should recompute.
        s.observations = vec![ObsBuilder::new()
            .lst_hours(6.0)
            .catalog_ra_hours(6.0)
            .catalog_dec_deg(45.0)
            .observed_dec_deg(45.0)
            .commanded_ha_arcsec(999.0)
            .build()];
        let prepared = s.prepared_observations();
        assert!(
            (prepared[0].commanded_ha.arcseconds() - 999.0).abs() > 1.0,
            "commanded_ha should be recomputed, not echoed from input",
        );
    }

    #[test]
    fn prepared_observations_preserves_input_order_and_count() {
        let mut s = Session::new();
        s.site = Some(site(39.0, 987.0));
        s.observations = (0..5)
            .map(|i| make_obs(0.0, 10.0 * i as f64))
            .collect();
        let prepared = s.prepared_observations();
        assert_eq!(prepared.len(), 5);
        for i in 0..5 {
            assert_eq!(
                prepared[i].observed_dec.degrees(),
                s.observations[i].observed_dec.degrees(),
            );
        }
    }

    // --- fit ----------------------------------------------------------

    #[test]
    fn fit_with_no_observations_errors() {
        let mut s = Session::new();
        s.model.add_term("IH").unwrap();
        assert!(s.fit().is_err());
    }

    #[test]
    fn fit_with_no_terms_errors() {
        let mut s = Session::new();
        s.observations = vec![make_obs(100.0, 30.0)];
        assert!(s.fit().is_err());
    }

    #[test]
    fn fit_recovers_ih_offset_and_stores_result() {
        // Three observations all with actual_ha = +100" and commanded_ha = 0.
        // An IH-only model should recover coefficient ≈ -100" because the
        // residual convention is (commanded - actual).
        let mut s = Session::new();
        s.observations = vec![
            make_obs(100.0, 30.0),
            make_obs(100.0, 45.0),
            make_obs(100.0, 60.0),
        ];
        s.model.add_term("IH").unwrap();

        let coefficients;
        let term_names;
        {
            let fit = s.fit().unwrap();
            coefficients = fit.coefficients.clone();
            term_names = fit.term_names.clone();
        }
        assert_eq!(term_names, vec!["IH"]);
        assert!(
            (coefficients[0] - (-100.0)).abs() < 1e-6,
            "IH coefficient should recover -100, got {}",
            coefficients[0],
        );
        assert!(s.last_fit.is_some());
        assert_eq!(s.model.coefficients(), &coefficients[..]);
    }

    #[test]
    fn fit_skips_masked_observations() {
        let mut s = Session::new();
        // Two normal obs at 100", one masked outlier at 9999".
        s.observations = vec![
            make_obs(100.0, 30.0),
            make_obs(100.0, 60.0),
            ObsBuilder::new()
                .catalog_dec_deg(45.0)
                .observed_dec_deg(45.0)
                .actual_ha_arcsec(9999.0)
                .masked(true)
                .build(),
        ];
        s.model.add_term("IH").unwrap();
        let result = s.fit().unwrap();
        // If the masked obs had been included, the LS fit would have been
        // pulled far away from -100".
        assert!(
            (result.coefficients[0] - (-100.0)).abs() < 1.0,
            "masked obs should not influence fit; got {}",
            result.coefficients[0],
        );
    }

    #[test]
    fn fit_overwrites_previous_last_fit() {
        let mut s = Session::new();
        s.last_fit = Some(fake_fit(vec![999.0]));
        s.observations = vec![
            make_obs(100.0, 30.0),
            make_obs(100.0, 45.0),
        ];
        s.model.add_term("IH").unwrap();
        s.fit().unwrap();
        assert!(s.last_fit.is_some());
        assert!(
            (s.last_fit.as_ref().unwrap().coefficients[0] - (-100.0)).abs() < 1e-6,
        );
    }

    #[test]
    fn fit_writes_coefficients_back_to_model() {
        let mut s = Session::new();
        s.observations = vec![
            make_obs(50.0, 30.0),
            make_obs(50.0, 45.0),
        ];
        s.model.add_term("IH").unwrap();
        s.fit().unwrap();
        // Model coefficients should match fit result, not the initial 0.0.
        assert!((s.model.coefficients()[0] - (-50.0)).abs() < 1e-6);
    }

    #[test]
    fn fit_with_all_observations_masked_errors() {
        let mut s = Session::new();
        s.observations = vec![ObsBuilder::new()
            .catalog_dec_deg(30.0)
            .observed_dec_deg(30.0)
            .actual_ha_arcsec(100.0)
            .masked(true)
            .build()];
        s.model.add_term("IH").unwrap();
        assert!(s.fit().is_err());
    }
}
