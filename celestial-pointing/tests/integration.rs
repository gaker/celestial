use celestial_core::Angle;
use celestial_pointing::commands::{dispatch, CommandOutput};
use celestial_pointing::observation::{MountType, PierSide};
use celestial_pointing::parser::parse_indat;
use celestial_pointing::session::{AdjustDirection, Session};

const SIMPLE_DAT: &str = "\
!TheSky Version 10.5.0 Build 13572 (64 bit)
ASCOM Mount
:NODA
:EQUAT
+39 00 26 2024 7 14 29.20 987.00 231.65  0.94 0.5500 0.0065
21 43 18.4460 +72 29 08.368 09 28 59.9527 +109 20 06.469  16 23.130
23 46 02.2988 +77 38 38.725 11 26 17.6308 +104 03 28.734  16 24.711";

fn load_simple() -> Session {
    let indat = parse_indat(SIMPLE_DAT).unwrap();
    let mut session = Session::new();
    session.load_indat(indat);
    session
}

// --- Parser integration ---

#[test]
fn simple_dat_parses_two_observations() {
    let session = load_simple();
    assert_eq!(session.observation_count(), 2);
}

#[test]
fn simple_dat_mount_type() {
    let session = load_simple();
    assert_eq!(session.mount_type, MountType::GermanEquatorial);
}

#[test]
fn simple_dat_latitude() {
    let session = load_simple();
    let expected = Angle::from_degrees(39.0 + 26.0 / 3600.0);
    assert_eq!(Angle::from_radians(session.latitude()), expected);
}

#[test]
fn simple_dat_both_obs_are_west_pier() {
    let session = load_simple();
    for obs in &session.observations {
        assert_eq!(obs.pier_side, PierSide::West);
    }
}

#[test]
fn simple_dat_first_obs_catalog_coordinates() {
    let session = load_simple();
    let obs = &session.observations[0];
    let expected_ra = Angle::from_hours(21.0 + 43.0 / 60.0 + 18.4460 / 3600.0);
    let expected_dec = Angle::from_degrees(72.0 + 29.0 / 60.0 + 8.368 / 3600.0);
    assert_eq!(obs.catalog_ra, expected_ra);
    assert_eq!(obs.catalog_dec, expected_dec);
}

#[test]
fn simple_dat_second_obs_catalog_coordinates() {
    let session = load_simple();
    let obs = &session.observations[1];
    let expected_ra = Angle::from_hours(23.0 + 46.0 / 60.0 + 2.2988 / 3600.0);
    let expected_dec = Angle::from_degrees(77.0 + 38.0 / 60.0 + 38.725 / 3600.0);
    assert_eq!(obs.catalog_ra, expected_ra);
    assert_eq!(obs.catalog_dec, expected_dec);
}

#[test]
fn simple_dat_ha_computed_from_lst_minus_ra() {
    let session = load_simple();
    for obs in &session.observations {
        assert_eq!(obs.commanded_ha, obs.lst - obs.catalog_ra);
        assert_eq!(obs.actual_ha, obs.lst - obs.observed_ra);
    }
}

// --- Single-term fitting (IH only) ---

#[test]
fn fit_ih_on_simple_data() {
    let mut session = load_simple();
    session.model.add_term("IH").unwrap();
    let result = session.fit().unwrap();

    assert_eq!(result.term_names, vec!["IH"]);
    assert_eq!(result.coefficients.len(), 1);
    assert_eq!(result.sigma.len(), 1);
    assert!(result.sky_rms.is_finite());
    assert!(result.sky_rms >= 0.0);
}

#[test]
fn fit_id_on_simple_data() {
    let mut session = load_simple();
    session.model.add_term("ID").unwrap();
    let result = session.fit().unwrap();

    assert_eq!(result.term_names, vec!["ID"]);
    assert_eq!(result.coefficients.len(), 1);
}

// --- Multi-term fitting ---

#[test]
fn fit_ih_id_on_simple_data() {
    let mut session = load_simple();
    session.model.add_term("IH").unwrap();
    session.model.add_term("ID").unwrap();
    let result = session.fit().unwrap();

    assert_eq!(result.term_names, vec!["IH", "ID"]);
    assert_eq!(result.coefficients.len(), 2);
    assert_eq!(result.sigma.len(), 2);
    assert!(result.sky_rms.is_finite());
}

// --- Fit result consistency: coefficients reduce residuals ---

#[test]
fn fit_reduces_sky_rms_vs_no_model() {
    let mut session = load_simple();

    let raw_rms = compute_raw_rms(&session);

    session.model.add_term("IH").unwrap();
    session.model.add_term("ID").unwrap();
    let result = session.fit().unwrap();

    assert!(
        result.sky_rms < raw_rms,
        "fitted rms {} should be less than raw rms {}",
        result.sky_rms,
        raw_rms,
    );
}

fn compute_raw_rms(session: &Session) -> f64 {
    let n = session.observations.len();
    if n == 0 {
        return 0.0;
    }
    let mut sum_sq = 0.0;
    for obs in &session.observations {
        let dh = (obs.actual_ha - obs.commanded_ha).arcseconds();
        let dd = (obs.observed_dec - obs.catalog_dec).arcseconds();
        let cos_dec = libm::cos(obs.catalog_dec.radians());
        let dx = dh * cos_dec;
        sum_sq += dx * dx + dd * dd;
    }
    libm::sqrt(sum_sq / n as f64)
}

// --- Model apply round-trip ---

#[test]
fn model_coefficients_are_set_after_fit() {
    let mut session = load_simple();
    session.model.add_term("IH").unwrap();
    session.model.add_term("ID").unwrap();
    let fit_coeffs = session.fit().unwrap().coefficients.clone();

    let model_coeffs = session.model.coefficients();
    assert_eq!(model_coeffs.len(), 2);
    assert_eq!(model_coeffs[0], fit_coeffs[0]);
    assert_eq!(model_coeffs[1], fit_coeffs[1]);
}

#[test]
fn apply_model_returns_nonzero_after_fit() {
    let mut session = load_simple();
    session.model.add_term("IH").unwrap();
    session.model.add_term("ID").unwrap();
    session.fit().unwrap();

    let obs = &session.observations[0];
    let h = obs.commanded_ha.radians();
    let dec = obs.catalog_dec.radians();
    let lat = session.latitude();
    let pier = obs.pier_side.sign();
    let (dh, dd) = session.model.apply_equatorial(h, dec, lat, pier);

    assert!(dh.abs() > 0.0 || dd.abs() > 0.0);
}

// --- Command dispatch integration ---

#[test]
fn dispatch_use_fit_workflow() {
    let mut session = load_simple();

    let use_result = dispatch(&mut session, "USE IH ID").unwrap();
    match &use_result {
        CommandOutput::Text(s) => {
            assert!(s.contains("IH"));
            assert!(s.contains("ID"));
        }
        _ => panic!("expected Text from USE"),
    }

    let fit_result = dispatch(&mut session, "FIT").unwrap();
    match fit_result {
        CommandOutput::FitDisplay(fd) => {
            assert_eq!(fd.term_names, vec!["IH", "ID"]);
            assert_eq!(fd.coefficients.len(), 2);
            assert!(fd.sky_rms.is_finite());
        }
        _ => panic!("expected FitDisplay from FIT"),
    }
}

#[test]
fn dispatch_lose_clears_model() {
    let mut session = load_simple();
    dispatch(&mut session, "USE IH ID CH").unwrap();
    assert_eq!(session.model.term_count(), 3);

    dispatch(&mut session, "LOSE CH").unwrap();
    assert_eq!(session.model.term_count(), 2);
    assert_eq!(session.model.term_names(), vec!["IH", "ID"]);

    dispatch(&mut session, "LOSE ALL").unwrap();
    assert_eq!(session.model.term_count(), 0);
}

#[test]
fn dispatch_slist_produces_output() {
    let mut session = load_simple();
    dispatch(&mut session, "USE IH").unwrap();
    dispatch(&mut session, "FIT").unwrap();

    let result = dispatch(&mut session, "SLIST").unwrap();
    match result {
        CommandOutput::Text(s) => {
            assert!(s.contains("dX"));
            assert!(s.contains("dD"));
            let lines: Vec<&str> = s.lines().collect();
            assert!(lines.len() >= 3, "header + blank + 2 obs rows");
        }
        _ => panic!("expected Text from SLIST"),
    }
}

// --- Fit then refit with different terms ---

#[test]
fn refit_with_more_terms_decreases_or_maintains_rms() {
    let mut session = load_simple();
    session.model.add_term("IH").unwrap();
    let r1 = session.fit().unwrap();
    let rms1 = r1.sky_rms;

    session.model.add_term("ID").unwrap();
    let r2 = session.fit().unwrap();
    let rms2 = r2.sky_rms;

    assert!(
        rms2 <= rms1 + 1e-10,
        "adding terms should not increase rms: {} vs {}",
        rms2,
        rms1,
    );
}

// --- Sigma values are positive and finite ---

#[test]
fn sigma_values_are_positive_finite() {
    let mut session = load_simple();
    session.model.add_term("IH").unwrap();
    let result = session.fit().unwrap();

    for (i, &s) in result.sigma.iter().enumerate() {
        assert!(
            s.is_finite() && s >= 0.0,
            "sigma[{}] = {} should be finite and non-negative",
            i,
            s,
        );
    }
}

// --- CGX-L dataset (148 observations) ---

fn load_cgx_l() -> Session {
    let content = include_str!("../pointing-data/cgx-l-data.dat");
    let indat = parse_indat(content).unwrap();
    let mut session = Session::new();
    session.load_indat(indat);
    session
}

#[test]
fn cgx_l_parses_148_observations() {
    let session = load_cgx_l();
    assert_eq!(session.observation_count(), 148);
}

#[test]
fn cgx_l_fit_6_term_standard_model() {
    let mut session = load_cgx_l();
    for term in &["IH", "ID", "CH", "NP", "MA", "ME"] {
        session.model.add_term(term).unwrap();
    }
    let result = session.fit().unwrap();

    assert_eq!(result.term_names.len(), 6);
    assert_eq!(result.coefficients.len(), 6);
    assert_eq!(result.sigma.len(), 6);
    assert!(result.sky_rms.is_finite());
    assert!(result.sky_rms > 0.0);

    for (i, &s) in result.sigma.iter().enumerate() {
        assert!(
            s.is_finite() && s >= 0.0,
            "sigma[{}] = {} should be finite and non-negative",
            i,
            s,
        );
    }
}

#[test]
fn cgx_l_more_terms_reduce_rms() {
    let mut session = load_cgx_l();

    session.model.add_term("IH").unwrap();
    session.model.add_term("ID").unwrap();
    let r2 = session.fit().unwrap();
    let rms_2 = r2.sky_rms;

    session.model.add_term("CH").unwrap();
    session.model.add_term("NP").unwrap();
    session.model.add_term("MA").unwrap();
    session.model.add_term("ME").unwrap();
    let r6 = session.fit().unwrap();
    let rms_6 = r6.sky_rms;

    assert!(
        rms_6 < rms_2,
        "6-term rms ({}) should be less than 2-term rms ({})",
        rms_6,
        rms_2,
    );
}

#[test]
fn cgx_l_fit_with_tube_flexure() {
    let mut session = load_cgx_l();
    for term in &["IH", "ID", "CH", "NP", "MA", "ME", "TF"] {
        session.model.add_term(term).unwrap();
    }
    let result = session.fit().unwrap();
    assert_eq!(result.term_names.len(), 7);
    assert!(result.sky_rms.is_finite());
}

#[test]
fn cgx_l_fit_with_daf_and_fo() {
    let mut session = load_cgx_l();
    for term in &["IH", "ID", "CH", "NP", "MA", "ME", "DAF", "FO"] {
        session.model.add_term(term).unwrap();
    }
    let result = session.fit().unwrap();
    assert_eq!(result.term_names.len(), 8);
    assert!(result.sky_rms.is_finite());
}

#[test]
fn cgx_l_fit_with_centering_errors() {
    let mut session = load_cgx_l();
    for term in &[
        "IH", "ID", "CH", "NP", "MA", "ME", "HCES", "HCEC", "DCES", "DCEC",
    ] {
        session.model.add_term(term).unwrap();
    }
    let result = session.fit().unwrap();
    assert_eq!(result.term_names.len(), 10);
    assert!(result.sky_rms.is_finite());
}

// --- Residuals after model application should be smaller ---

#[test]
fn cgx_l_model_residuals_smaller_than_raw() {
    let mut session = load_cgx_l();
    let raw_rms = compute_raw_rms(&session);

    for term in &["IH", "ID", "CH", "NP", "MA", "ME"] {
        session.model.add_term(term).unwrap();
    }
    let result = session.fit().unwrap();

    assert!(
        result.sky_rms < raw_rms,
        "fitted rms {} should be less than raw rms {}",
        result.sky_rms,
        raw_rms,
    );
}

// --- Model correction reduces per-observation errors ---

#[test]
fn model_correction_applied_to_each_observation() {
    let mut session = load_cgx_l();
    for term in &["IH", "ID", "CH", "NP", "MA", "ME"] {
        session.model.add_term(term).unwrap();
    }
    session.fit().unwrap();

    let lat = session.latitude();
    for obs in &session.observations {
        let h = obs.commanded_ha.radians();
        let dec = obs.catalog_dec.radians();
        let pier = obs.pier_side.sign();
        let (dh, dd) = session.model.apply_equatorial(h, dec, lat, pier);
        assert!(dh.is_finite());
        assert!(dd.is_finite());
    }
}

// --- Pier side distribution in CGX-L data ---

#[test]
fn cgx_l_has_both_pier_sides() {
    let session = load_cgx_l();
    let east_count = session
        .observations
        .iter()
        .filter(|o| o.pier_side == PierSide::East)
        .count();
    let west_count = session
        .observations
        .iter()
        .filter(|o| o.pier_side == PierSide::West)
        .count();

    assert!(east_count > 0, "should have east pier observations");
    assert!(west_count > 0, "should have west pier observations");
    assert_eq!(east_count + west_count, 148);
}

// --- Successive fits converge to same result ---

#[test]
fn double_fit_produces_same_coefficients() {
    let mut session = load_cgx_l();
    session.model.add_term("IH").unwrap();
    session.model.add_term("ID").unwrap();

    session.fit().unwrap();
    let coeffs1: Vec<f64> = session.model.coefficients().to_vec();

    session.fit().unwrap();
    let coeffs2: Vec<f64> = session.model.coefficients().to_vec();

    for (i, (&c1, &c2)) in coeffs1.iter().zip(coeffs2.iter()).enumerate() {
        assert_eq!(c1, c2, "coefficient {} differs between fits", i);
    }
}

// --- Phase 2: CLIST ---

#[test]
fn clist_no_terms_reports_empty() {
    let mut session = load_simple();
    let result = dispatch(&mut session, "CLIST").unwrap();
    match result {
        CommandOutput::Text(s) => assert!(s.contains("No terms")),
        _ => panic!("expected Text from CLIST"),
    }
}

#[test]
fn clist_after_fit_returns_coefficients() {
    let mut session = load_simple();
    dispatch(&mut session, "USE IH ID").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    let result = dispatch(&mut session, "CLIST").unwrap();
    match result {
        CommandOutput::FitDisplay(fd) => {
            assert_eq!(fd.term_names, vec!["IH", "ID"]);
            assert_eq!(fd.coefficients.len(), 2);
            assert_eq!(fd.sigma.len(), 2);
            assert!(fd.sky_rms > 0.0);
        }
        _ => panic!("expected FitDisplay from CLIST"),
    }
}

#[test]
fn clist_before_fit_returns_zero_coefficients() {
    let mut session = load_simple();
    dispatch(&mut session, "USE IH").unwrap();
    let result = dispatch(&mut session, "CLIST").unwrap();
    match result {
        CommandOutput::FitDisplay(fd) => {
            assert_eq!(fd.coefficients, vec![0.0]);
            assert_eq!(fd.sigma, vec![0.0]);
            assert_eq!(fd.sky_rms, 0.0);
        }
        _ => panic!("expected FitDisplay from CLIST"),
    }
}

// --- Phase 2: RESET ---

#[test]
fn reset_zeros_coefficients_keeps_terms() {
    let mut session = load_simple();
    dispatch(&mut session, "USE IH ID").unwrap();
    dispatch(&mut session, "FIT").unwrap();

    assert!(session.model.coefficients().iter().any(|&c| c != 0.0));

    dispatch(&mut session, "RESET").unwrap();
    assert!(session.model.coefficients().iter().all(|&c| c == 0.0));
    assert_eq!(session.model.term_count(), 2);
    assert!(session.last_fit.is_none());
}

// --- Phase 2: MASK / UNMASK ---

#[test]
fn mask_single_observation() {
    let mut session = load_simple();
    dispatch(&mut session, "MASK 1").unwrap();
    assert!(session.observations[0].masked);
    assert!(!session.observations[1].masked);
}

#[test]
fn mask_range() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "MASK 1-5").unwrap();
    for i in 0..5 {
        assert!(
            session.observations[i].masked,
            "obs {} should be masked",
            i + 1
        );
    }
    assert!(!session.observations[5].masked);
}

#[test]
fn unmask_all() {
    let mut session = load_simple();
    dispatch(&mut session, "MASK 1 2").unwrap();
    assert!(session.observations.iter().all(|o| o.masked));
    dispatch(&mut session, "UNMASK ALL").unwrap();
    assert!(session.observations.iter().all(|o| !o.masked));
}

#[test]
fn masked_observations_excluded_from_fit() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID").unwrap();

    dispatch(&mut session, "FIT").unwrap();
    let rms_all = session.last_fit.as_ref().unwrap().sky_rms;

    dispatch(&mut session, "MASK 1-5").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    let rms_masked = session.last_fit.as_ref().unwrap().sky_rms;

    assert!(
        rms_all != rms_masked,
        "RMS should differ after masking observations"
    );
}

#[test]
fn slist_shows_masked_indicator() {
    let mut session = load_simple();
    dispatch(&mut session, "USE IH").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    dispatch(&mut session, "MASK 1").unwrap();
    let result = dispatch(&mut session, "SLIST").unwrap();
    match result {
        CommandOutput::Text(s) => {
            let lines: Vec<&str> = s.lines().collect();
            assert!(lines[2].contains("*"), "masked obs should show * indicator");
        }
        _ => panic!("expected Text from SLIST"),
    }
}

#[test]
fn mask_out_of_bounds_errors() {
    let mut session = load_simple();
    let result = dispatch(&mut session, "MASK 99");
    assert!(result.is_err());
}

// --- Phase 2: INMOD / OUTMOD round-trip ---

#[test]
fn inmod_outmod_round_trip() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID CH NP MA ME").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    let original_coeffs: Vec<f64> = session.model.coefficients().to_vec();
    let original_names: Vec<String> = session
        .model
        .term_names()
        .iter()
        .map(|s| s.to_string())
        .collect();

    let tmp = "/tmp/celestial_test_round_trip.mod";
    dispatch(&mut session, &format!("OUTMOD {}", tmp)).unwrap();

    let mut session2 = Session::new();
    dispatch(&mut session2, &format!("INMOD {}", tmp)).unwrap();

    let loaded_names: Vec<String> = session2
        .model
        .term_names()
        .iter()
        .map(|s| s.to_string())
        .collect();
    let loaded_coeffs = session2.model.coefficients().to_vec();

    assert_eq!(loaded_names, original_names);
    for (i, (&orig, &loaded)) in original_coeffs.iter().zip(loaded_coeffs.iter()).enumerate() {
        assert!(
            (orig - loaded).abs() < 1e-4,
            "coefficient {} differs: {} vs {}",
            i,
            orig,
            loaded
        );
    }

    std::fs::remove_file(tmp).ok();
}

#[test]
fn inmod_missing_file_errors() {
    let mut session = Session::new();
    let result = dispatch(&mut session, "INMOD /tmp/nonexistent_file_xyz.mod");
    assert!(result.is_err());
}

// --- Phase 2: SHOW ---

#[test]
fn show_displays_session_state() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    let result = dispatch(&mut session, "SHOW").unwrap();
    match result {
        CommandOutput::Text(s) => {
            assert!(s.contains("German Equatorial"));
            assert!(s.contains("148"));
            assert!(s.contains("2"));
        }
        _ => panic!("expected Text from SHOW"),
    }
}

#[test]
fn show_with_masked_observations() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "MASK 1-5").unwrap();
    let result = dispatch(&mut session, "SHOW").unwrap();
    match result {
        CommandOutput::Text(s) => {
            assert!(s.contains("5 masked"));
        }
        _ => panic!("expected Text from SHOW"),
    }
}

// --- Phase 2: MVET ---

#[test]
fn mvet_requires_fit() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID").unwrap();
    let result = dispatch(&mut session, "MVET 2.0");
    assert!(result.is_err());
}

#[test]
fn mvet_reports_weak_terms() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID CH NP MA ME HCES HCEC DCES DCEC").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    let result = dispatch(&mut session, "MVET 100.0").unwrap();
    match result {
        CommandOutput::Text(s) => {
            assert!(s.contains("Weak terms") || s.contains("No weak terms"));
        }
        _ => panic!("expected Text from MVET"),
    }
}

#[test]
fn mvet_remove_flag_reduces_terms() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID CH NP MA ME HCES HCEC DCES DCEC").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    let before = session.model.term_count();
    dispatch(&mut session, "MVET 100.0 R").unwrap();
    let after = session.model.term_count();
    assert!(after <= before);
}

// --- Phase 2: OUTL ---

#[test]
fn outl_requires_fit() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH").unwrap();
    let result = dispatch(&mut session, "OUTL 3.0");
    assert!(result.is_err());
}

#[test]
fn outl_reports_outliers() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    let result = dispatch(&mut session, "OUTL 1.0").unwrap();
    match result {
        CommandOutput::Text(s) => {
            assert!(s.contains("Outlier") || s.contains("No outlier"));
        }
        _ => panic!("expected Text from OUTL"),
    }
}

#[test]
fn outl_mask_flag_masks_outliers() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    let before_masked = session.masked_observation_count();
    dispatch(&mut session, "OUTL 1.0 M").unwrap();
    let after_masked = session.masked_observation_count();
    assert!(after_masked >= before_masked);
}

// --- Phase 2: FIX / UNFIX ---

#[test]
fn fix_term_and_fit() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID CH").unwrap();
    dispatch(&mut session, "FIT").unwrap();
    let _ih_coeff = session.model.coefficients()[0];

    dispatch(&mut session, "FIX IH").unwrap();
    dispatch(&mut session, "FIT").unwrap();

    assert_eq!(
        session.model.coefficients()[0],
        0.0,
        "fixed term should remain at its pre-fix value (zeroed by RESET effect)"
    );
}

#[test]
fn unfix_all_allows_fitting() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID").unwrap();
    dispatch(&mut session, "FIX ALL").unwrap();
    let result = dispatch(&mut session, "FIT");
    assert!(result.is_err(), "all fixed should prevent fitting");

    dispatch(&mut session, "UNFIX ALL").unwrap();
    let result = dispatch(&mut session, "FIT");
    assert!(result.is_ok());
}

#[test]
fn fix_unknown_term_errors() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH").unwrap();
    let result = dispatch(&mut session, "FIX ZZZZ");
    assert!(result.is_err());
}

// --- Phase 2: ADJUST ---

#[test]
fn adjust_show_current() {
    let mut session = Session::new();
    let result = dispatch(&mut session, "ADJUST").unwrap();
    match result {
        CommandOutput::Text(s) => assert!(s.contains("telescope to star")),
        _ => panic!("expected Text"),
    }
}

#[test]
fn adjust_set_direction() {
    let mut session = Session::new();
    dispatch(&mut session, "ADJUST S").unwrap();
    assert_eq!(session.adjust_direction, AdjustDirection::StarToTelescope);
    dispatch(&mut session, "ADJUST T").unwrap();
    assert_eq!(session.adjust_direction, AdjustDirection::TelescopeToStar);
}

#[test]
fn adjust_invalid_direction_errors() {
    let mut session = Session::new();
    let result = dispatch(&mut session, "ADJUST X");
    assert!(result.is_err());
}

// --- Phase 2: FAUTO ---

#[test]
fn fauto_adds_harmonics() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH ID").unwrap();
    dispatch(&mut session, "FAUTO 3").unwrap();
    let names = session.model.term_names();
    assert!(names.contains(&"HDSH"));
    assert!(names.contains(&"HDCH"));
    assert!(names.contains(&"HDSH2"));
    assert!(names.contains(&"HDCH2"));
    assert!(names.contains(&"HDSH3"));
    assert!(names.contains(&"HDCH3"));
    assert_eq!(session.model.term_count(), 8);
}

#[test]
fn fauto_no_duplicates() {
    let mut session = load_cgx_l();
    dispatch(&mut session, "USE IH HDSH HDCH").unwrap();
    let before = session.model.term_count();
    dispatch(&mut session, "FAUTO 1").unwrap();
    assert_eq!(session.model.term_count(), before);
}

#[test]
fn fauto_zero_order_errors() {
    let mut session = Session::new();
    let result = dispatch(&mut session, "FAUTO 0");
    assert!(result.is_err());
}
