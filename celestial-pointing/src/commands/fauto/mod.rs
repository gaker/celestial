mod bic;
mod candidates;
mod moves;
mod prune;
mod search;

use super::{Command, CommandOutput};
use crate::error::{Error, Result};
use crate::session::Session;
use crate::solver::FitResult;
use search::run_search_with_outlier_passes;

pub(super) const BASE_TERMS: &[&str] = &["IH", "ID", "CH", "NP", "MA", "ME"];
const DEFAULT_MAX_TERMS: usize = 30;
const DEFAULT_BIC_THRESHOLD: f64 = -6.0;

pub struct Fauto;

impl Command for Fauto {
    fn name(&self) -> &str {
        "FAUTO"
    }
    fn description(&self) -> &str {
        "Automatic modeling: pair-aware BIC search over harmonics + physical terms"
    }

    fn execute(&self, session: &mut Session, args: &[&str]) -> Result<CommandOutput> {
        let (max_terms, bic_threshold) = parse_args(args)?;
        // Per the TPOINT manual: any MASKed observations are reactivated before
        // modeling. Clear masks on the session so the prepared list is clean,
        // and the outlier pass below builds the new mask set from scratch.
        for o in session.observations.iter_mut() {
            o.masked = false;
        }
        let active: Vec<String> = seed_terms_from_session(session);
        let prepared = session.prepared_observations();
        if prepared.len() < active.len() {
            return Err(Error::Fit("insufficient observations for FAUTO".into()));
        }
        let latitude = session.latitude();
        let fit_tol = session.fit_tol;

        let mut masked: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut active = active;
        let mut parallel_groups: Vec<(String, String)> = Vec::new();
        let mut report = String::from("FAUTO model search...\n");

        let (final_fit, newly_masked) = run_search_with_outlier_passes(
            &prepared,
            &mut masked,
            &mut active,
            &mut parallel_groups,
            &mut report,
            max_terms,
            bic_threshold,
            latitude,
            fit_tol,
        )?;

        for idx in &newly_masked {
            session.observations[*idx].masked = true;
        }
        write_summary(&mut report, &active, &final_fit);
        load_into_session(session, &active, &parallel_groups, &final_fit)?;
        Ok(CommandOutput::Text(report))
    }
}

/// Pick the starting term set. Matches TPOINT's "accept any pre-existing model
/// as a starting point" behavior: extend whatever the session already has, only
/// falling back to the canonical six-term equatorial base when the session is
/// empty. The base terms are also unioned in so a partial model can't omit
/// them; FAUTO treats them as a structural floor.
fn seed_terms_from_session(session: &Session) -> Vec<String> {
    let mut out: Vec<String> = BASE_TERMS.iter().map(|s| s.to_string()).collect();
    for name in session.model.term_names() {
        if !out.iter().any(|t| t == name) {
            out.push(name.to_string());
        }
    }
    out
}

fn parse_args(args: &[&str]) -> Result<(usize, f64)> {
    let max_terms = match args.first() {
        Some(s) => s
            .parse::<usize>()
            .map_err(|e| Error::Parse(format!("invalid max_terms: {}", e)))?,
        None => DEFAULT_MAX_TERMS,
    };
    let bic_threshold = match args.get(1) {
        Some(s) => s
            .parse::<f64>()
            .map_err(|e| Error::Parse(format!("invalid bic_threshold: {}", e)))?,
        None => DEFAULT_BIC_THRESHOLD,
    };
    Ok((max_terms, bic_threshold))
}

fn write_summary(report: &mut String, active: &[String], final_fit: &FitResult) {
    report.push_str(&format!(
        "\nFinal model: {} terms, RMS={:.2}\"\n",
        active.len(),
        final_fit.sky_rms,
    ));
    report.push_str("Terms: ");
    report.push_str(&active.join(" "));
    report.push('\n');
}

fn load_into_session(
    session: &mut Session,
    active: &[String],
    parallel_groups: &[(String, String)],
    result: &FitResult,
) -> Result<()> {
    session.model.remove_all();
    for name in active {
        session.model.add_term(name)?;
    }
    for (_, partner) in parallel_groups {
        if active.contains(partner) {
            session.model.set_parallel(partner);
        }
    }
    session.model.set_coefficients(&result.coefficients)?;
    session.last_fit = Some(result.clone());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::FitResultBuilder;

    fn parse_err(err: &Error) -> &str {
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

    fn fake_fit(coefficients: Vec<f64>) -> FitResult {
        FitResultBuilder::new()
            .coefficients(coefficients)
            .sky_rms(1.0)
            .popn_sd(1.0)
            .build()
    }

    #[test]
    fn metadata() {
        assert_eq!(Fauto.name(), "FAUTO");
        assert!(Fauto.description().contains("BIC"));
    }

    #[test]
    fn parse_args_no_args_returns_defaults() {
        let (max_terms, threshold) = parse_args(&[]).unwrap();
        assert_eq!(max_terms, DEFAULT_MAX_TERMS);
        assert_eq!(threshold, DEFAULT_BIC_THRESHOLD);
    }

    #[test]
    fn parse_args_one_arg_sets_max_terms_only() {
        let (max_terms, threshold) = parse_args(&["12"]).unwrap();
        assert_eq!(max_terms, 12);
        assert_eq!(threshold, DEFAULT_BIC_THRESHOLD);
    }

    #[test]
    fn parse_args_two_args_sets_both() {
        let (max_terms, threshold) = parse_args(&["15", "-4.5"]).unwrap();
        assert_eq!(max_terms, 15);
        assert_eq!(threshold, -4.5);
    }

    #[test]
    fn parse_args_extra_args_are_ignored() {
        let (max_terms, threshold) = parse_args(&["10", "-3.0", "junk", "more"]).unwrap();
        assert_eq!(max_terms, 10);
        assert_eq!(threshold, -3.0);
    }

    #[test]
    fn parse_args_unparseable_max_terms_errors() {
        let err = parse_args(&["banana"]).unwrap_err();
        assert!(parse_err(&err).contains("invalid max_terms"));
    }

    #[test]
    fn parse_args_negative_max_terms_errors() {
        // usize::parse rejects negative integers.
        let err = parse_args(&["-1"]).unwrap_err();
        assert!(parse_err(&err).contains("invalid max_terms"));
    }

    #[test]
    fn parse_args_unparseable_bic_threshold_errors() {
        let err = parse_args(&["10", "banana"]).unwrap_err();
        assert!(parse_err(&err).contains("invalid bic_threshold"));
    }

    #[test]
    fn execute_with_no_observations_errors_with_fit() {
        let mut session = Session::new();
        let err = Fauto.execute(&mut session, &[]).unwrap_err();
        assert!(fit_err(&err).contains("insufficient observations"));
    }

    // The "insufficient" check counts unmasked observations against BASE_TERMS.len()
    // (6). Fewer than 6 unmasked observations should bail before any fitting.
    #[test]
    fn execute_with_too_few_unmasked_observations_errors() {
        use crate::test_support::obs;

        let mut session = Session::new();
        session.observations = (0..5).map(|_| obs()).collect();
        let err = Fauto.execute(&mut session, &[]).unwrap_err();
        assert!(fit_err(&err).contains("insufficient observations"));
    }

    #[test]
    fn parse_args_argument_errors_propagate_through_execute() {
        let mut session = Session::new();
        let err = Fauto.execute(&mut session, &["banana"]).unwrap_err();
        assert!(parse_err(&err).contains("invalid max_terms"));
    }

    #[test]
    fn load_into_session_clears_existing_model() {
        let mut session = Session::new();
        session.model.add_term("AW").unwrap();
        let result = fake_fit(vec![1.0, 2.0]);
        load_into_session(
            &mut session,
            &["IH".to_string(), "ID".to_string()],
            &[],
            &result,
        )
        .unwrap();
        assert_eq!(session.model.term_names(), vec!["IH", "ID"]);
        assert!(!session.model.term_names().contains(&"AW"));
    }

    #[test]
    fn load_into_session_records_parallel_partner_only_if_active() {
        let mut session = Session::new();
        let active = vec!["IH".to_string(), "HDSH".to_string(), "HDCH".to_string()];
        let groups = vec![("HDSH".to_string(), "HDCH".to_string())];
        let result = fake_fit(vec![1.0, 2.0, 3.0]);
        load_into_session(&mut session, &active, &groups, &result).unwrap();
        // The partner side of each pair should be flagged parallel.
        let idx_hdch = session
            .model
            .term_names()
            .iter()
            .position(|n| *n == "HDCH")
            .unwrap();
        assert!(session.model.parallel_flags()[idx_hdch]);
    }

    #[test]
    fn load_into_session_skips_parallel_if_partner_not_active() {
        let mut session = Session::new();
        let active = vec!["IH".to_string()];
        // Partner "GONE" isn't in active; this would no-op rather than error.
        let groups = vec![("IH".to_string(), "GONE".to_string())];
        let result = fake_fit(vec![1.0]);
        load_into_session(&mut session, &active, &groups, &result).unwrap();
        assert!(!session.model.parallel_flags().iter().any(|&p| p));
    }

    #[test]
    fn load_into_session_stores_fit_in_last_fit() {
        let mut session = Session::new();
        let result = fake_fit(vec![7.5]);
        load_into_session(&mut session, &["IH".to_string()], &[], &result).unwrap();
        let stored = session.last_fit.as_ref().unwrap();
        assert_eq!(stored.coefficients, vec![7.5]);
    }

    #[test]
    fn write_summary_includes_term_count_rms_and_term_list() {
        let mut report = String::new();
        let active = vec!["IH".to_string(), "ID".to_string(), "TF".to_string()];
        let result = fake_fit(vec![1.0, 2.0, 3.0]);
        write_summary(&mut report, &active, &result);
        assert!(report.contains("Final model: 3 terms"));
        assert!(report.contains("RMS=1.00\""));
        assert!(report.contains("Terms: IH ID TF"));
    }

    #[test]
    fn base_terms_includes_six_canonical_equatorial_terms() {
        // Constant pin — changing this changes the FAUTO baseline model.
        assert_eq!(BASE_TERMS, &["IH", "ID", "CH", "NP", "MA", "ME"]);
    }

    #[test]
    fn seed_terms_returns_base_when_session_model_is_empty() {
        let session = Session::new();
        let seed = seed_terms_from_session(&session);
        assert_eq!(seed, BASE_TERMS);
    }

    #[test]
    fn seed_terms_extends_session_model_with_base_terms_floor() {
        // TPOINT: "accept any pre-existing model as a starting point". Session
        // already has IH + HDCH; FAUTO should keep HDCH and ensure the rest of
        // the six-term base is present too.
        let mut session = Session::new();
        session.model.add_term("IH").unwrap();
        session.model.add_term("HDCH").unwrap();
        let seed = seed_terms_from_session(&session);
        for base in BASE_TERMS {
            assert!(seed.iter().any(|t| t == base), "missing base {}: {:?}", base, seed);
        }
        assert!(seed.iter().any(|t| t == "HDCH"));
    }

    #[test]
    fn fauto_reactivates_masked_observations_before_bailing() {
        // Per the TPOINT manual, FAUTO must reactivate masked observations
        // *before* modeling. Even when FAUTO bails (here: too few obs), the
        // input masks should already have been cleared.
        use crate::test_support::ObsBuilder;
        let mut session = Session::new();
        session.observations = (0..3)
            .map(|_| ObsBuilder::new().masked(true).build())
            .collect();
        let _ = Fauto.execute(&mut session, &[]); // insufficient obs, errors
        assert!(
            session.observations.iter().all(|o| !o.masked),
            "FAUTO did not reactivate masks before entry"
        );
    }

    #[test]
    fn default_max_terms_is_thirty() {
        assert_eq!(DEFAULT_MAX_TERMS, 30);
    }

    #[test]
    fn default_bic_threshold_is_negative_six() {
        assert_eq!(DEFAULT_BIC_THRESHOLD, -6.0);
    }
}
