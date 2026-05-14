use super::FitInputs;
use crate::error::{Error, Result};

pub(super) fn validate_fit_inputs(inputs: &FitInputs<'_>) -> Result<()> {
    let free_count = inputs.fixed.iter().filter(|&&f| !f).count();
    if free_count == 0 && !inputs.terms.is_empty() {
        return Err(Error::Fit("all terms are fixed".into()));
    }
    if inputs.terms.is_empty() {
        return Err(Error::Fit("no terms to fit".into()));
    }
    if inputs.observations.len() < free_count {
        return Err(Error::Fit("insufficient observations".into()));
    }
    if inputs.parallel.len() != inputs.terms.len() {
        return Err(Error::Fit("parallel flag length mismatch".into()));
    }
    Ok(())
}

pub(super) fn collect_indices(fixed: &[bool], match_fixed: bool) -> Vec<usize> {
    fixed
        .iter()
        .enumerate()
        .filter(|(_, &f)| f == match_fixed)
        .map(|(i, _)| i)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observation::Observation;
    use crate::terms::create_term;
    use crate::test_support::obs;

    fn make_inputs<'a>(
        observations: &'a [&'a Observation],
        terms: &'a [Box<dyn crate::terms::Term>],
        fixed: &'a [bool],
        parallel: &'a [bool],
        coefficients: &'a [f64],
    ) -> FitInputs<'a> {
        FitInputs {
            observations,
            terms,
            fixed,
            parallel,
            coefficients,
            latitude: 0.0,
            fit_tol: 1.0e-3,
        }
    }

    fn fit_msg(err: &Error) -> &str {
        match err {
            Error::Fit(m) => m,
            other => panic!("expected Fit, got {:?}", other),
        }
    }

    // --- validate_fit_inputs ------------------------------------------

    #[test]
    fn validate_empty_terms_errors_with_no_terms_to_fit() {
        let observations: Vec<&Observation> = Vec::new();
        let terms: Vec<Box<dyn crate::terms::Term>> = Vec::new();
        let inputs = make_inputs(&observations, &terms, &[], &[], &[]);
        let err = validate_fit_inputs(&inputs).unwrap_err();
        assert!(fit_msg(&err).contains("no terms to fit"));
    }

    #[test]
    fn validate_all_fixed_with_nonempty_terms_errors() {
        let o = obs();
        let observations: Vec<&Observation> = vec![&o];
        let terms = vec![create_term("IH").unwrap()];
        let fixed = vec![true];
        let parallel = vec![true];
        let coefficients = vec![0.0];
        let inputs = make_inputs(&observations, &terms, &fixed, &parallel, &coefficients);
        let err = validate_fit_inputs(&inputs).unwrap_err();
        assert!(fit_msg(&err).contains("all terms are fixed"));
    }

    #[test]
    fn validate_insufficient_observations_errors() {
        // 2 free terms but only 1 observation.
        let o = obs();
        let observations: Vec<&Observation> = vec![&o];
        let terms = vec![create_term("IH").unwrap(), create_term("ID").unwrap()];
        let fixed = vec![false, false];
        let parallel = vec![true, true];
        let coefficients = vec![0.0, 0.0];
        let inputs = make_inputs(&observations, &terms, &fixed, &parallel, &coefficients);
        let err = validate_fit_inputs(&inputs).unwrap_err();
        assert!(fit_msg(&err).contains("insufficient observations"));
    }

    #[test]
    fn validate_parallel_length_mismatch_errors() {
        let o = obs();
        let observations: Vec<&Observation> = vec![&o, &o];
        let terms = vec![create_term("IH").unwrap()];
        let fixed = vec![false];
        let parallel = vec![true, true]; // mismatched length
        let coefficients = vec![0.0];
        let inputs = make_inputs(&observations, &terms, &fixed, &parallel, &coefficients);
        let err = validate_fit_inputs(&inputs).unwrap_err();
        assert!(fit_msg(&err).contains("parallel flag length mismatch"));
    }

    #[test]
    fn validate_well_formed_inputs_succeed() {
        let o = obs();
        let observations: Vec<&Observation> = vec![&o, &o, &o];
        let terms = vec![create_term("IH").unwrap()];
        let fixed = vec![false];
        let parallel = vec![true];
        let coefficients = vec![0.0];
        let inputs = make_inputs(&observations, &terms, &fixed, &parallel, &coefficients);
        assert!(validate_fit_inputs(&inputs).is_ok());
    }

    // Edge case: zero observations + zero terms — currently errors with "no
    // terms to fit" because the empty-terms check fires before observation
    // counting. Pinning the precedence.
    #[test]
    fn validate_empty_everything_errors_on_terms_check_first() {
        let observations: Vec<&Observation> = Vec::new();
        let terms: Vec<Box<dyn crate::terms::Term>> = Vec::new();
        let inputs = make_inputs(&observations, &terms, &[], &[], &[]);
        let err = validate_fit_inputs(&inputs).unwrap_err();
        assert!(fit_msg(&err).contains("no terms to fit"));
    }

    // --- collect_indices ----------------------------------------------

    #[test]
    fn collect_indices_returns_matching_positions() {
        let fixed = vec![true, false, true, false, false];
        assert_eq!(collect_indices(&fixed, true), vec![0, 2]);
        assert_eq!(collect_indices(&fixed, false), vec![1, 3, 4]);
    }

    #[test]
    fn collect_indices_partitions_completely() {
        let fixed = vec![true, false, true, false, false];
        let free = collect_indices(&fixed, false);
        let fxd = collect_indices(&fixed, true);
        let mut combined: Vec<usize> = free.iter().chain(fxd.iter()).copied().collect();
        combined.sort();
        assert_eq!(combined, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn collect_indices_empty_input_returns_empty() {
        assert!(collect_indices(&[], true).is_empty());
        assert!(collect_indices(&[], false).is_empty());
    }

    #[test]
    fn collect_indices_all_match_returns_all() {
        let fixed = vec![false, false, false];
        assert_eq!(collect_indices(&fixed, false), vec![0, 1, 2]);
        assert!(collect_indices(&fixed, true).is_empty());
    }
}
