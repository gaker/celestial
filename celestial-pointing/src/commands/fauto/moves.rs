use super::bic::{compute_bic, fit_is_well_conditioned, try_fit};
use super::candidates::CandidatePool;
use crate::error::Result;
use crate::observation::Observation;
use rayon::prelude::*;

pub(super) struct Move {
    pub(super) names: Vec<String>,
    pub(super) pair_idx: Option<(usize, usize)>,
    pub(super) rms: f64,
}

impl Move {
    pub(super) fn added_count(&self) -> usize {
        self.names.len()
    }

    pub(super) fn label(&self) -> String {
        let mut parts: Vec<String> = self.names.to_vec();
        if let Some((i, j)) = self.pair_idx {
            if i < parts.len() && j < parts.len() {
                let merged = format!("{}&{}", parts[i], parts[j]);
                let (lo, hi) = if i < j { (i, j) } else { (j, i) };
                parts.remove(hi);
                parts[lo] = merged;
            }
        }
        parts.join(" ")
    }
}

pub(super) fn commit_move(
    active: &mut Vec<String>,
    parallel_groups: &mut Vec<(String, String)>,
    m: &Move,
) {
    let start = active.len();
    for n in &m.names {
        active.push(n.clone());
    }
    if let Some((i, j)) = m.pair_idx {
        let a = active[start + i].clone();
        let b = active[start + j].clone();
        parallel_groups.push((a, b));
    }
}

pub(super) fn cmp_by_bic(
    a: &Move,
    b: &Move,
    base_len: usize,
    n_obs: usize,
) -> std::cmp::Ordering {
    let bic_a = compute_bic(n_obs, base_len + a.added_count(), a.rms);
    let bic_b = compute_bic(n_obs, base_len + b.added_count(), b.rms);
    bic_a.partial_cmp(&bic_b).unwrap_or(std::cmp::Ordering::Equal)
}

pub(super) fn trial_combo(
    observations: &[&Observation],
    active: &[String],
    names: &[String],
    pair_idx: Option<(usize, usize)>,
    latitude: f64,
    fit_tol: f64,
) -> Option<Move> {
    let mut trial = active.to_vec();
    for n in names {
        trial.push(n.clone());
    }
    let fit = try_fit(observations, &trial, latitude, fit_tol).ok()?;
    if !fit_is_well_conditioned(&fit) {
        return None;
    }
    Some(Move {
        names: names.to_vec(),
        pair_idx,
        rms: fit.sky_rms,
    })
}

pub(super) fn score_best_move(
    observations: &[&Observation],
    active: &[String],
    candidates: &CandidatePool,
    latitude: f64,
    fit_tol: f64,
) -> Result<Option<Move>> {
    let n_obs = observations.len();
    let base_len = active.len();
    let singles_best = candidates
        .singles
        .par_iter()
        .filter_map(|name| {
            trial_combo(
                observations,
                active,
                std::slice::from_ref(name),
                None,
                latitude,
                fit_tol,
            )
        })
        .min_by(|a, b| cmp_by_bic(a, b, base_len, n_obs));
    let pairs_best = candidates
        .pairs
        .par_iter()
        .filter_map(|(a, b)| {
            trial_combo(
                observations,
                active,
                &[a.clone(), b.clone()],
                Some((0, 1)),
                latitude,
                fit_tol,
            )
        })
        .min_by(|a, b| cmp_by_bic(a, b, base_len, n_obs));
    let triples_best = candidates
        .triples
        .par_iter()
        .filter_map(|(a, b, c)| {
            trial_combo(
                observations,
                active,
                &[a.clone(), b.clone(), c.clone()],
                Some((0, 1)),
                latitude,
                fit_tol,
            )
        })
        .min_by(|a, b| cmp_by_bic(a, b, base_len, n_obs));
    Ok([singles_best, pairs_best, triples_best]
        .into_iter()
        .flatten()
        .min_by(|a, b| cmp_by_bic(a, b, base_len, n_obs)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mv(names: Vec<&str>, pair_idx: Option<(usize, usize)>, rms: f64) -> Move {
        Move {
            names: names.into_iter().map(String::from).collect(),
            pair_idx,
            rms,
        }
    }

    #[test]
    fn added_count_matches_names_len() {
        assert_eq!(mv(vec!["A"], None, 1.0).added_count(), 1);
        assert_eq!(mv(vec!["A", "B"], Some((0, 1)), 1.0).added_count(), 2);
        assert_eq!(mv(vec!["A", "B", "C"], Some((0, 1)), 1.0).added_count(), 3);
    }

    #[test]
    fn label_single_term_has_no_ampersand() {
        let m = mv(vec!["HDSH"], None, 1.0);
        assert_eq!(m.label(), "HDSH");
    }

    #[test]
    fn label_pair_joins_with_ampersand() {
        let m = mv(vec!["HDSH", "HDCH"], Some((0, 1)), 1.0);
        assert_eq!(m.label(), "HDSH&HDCH");
    }

    #[test]
    fn label_triple_merges_first_pair_keeping_third_separate() {
        // names = [A, B, C], pair_idx = (0,1) → "A&B C"
        let m = mv(vec!["HDSH", "HDCH", "TF"], Some((0, 1)), 1.0);
        assert_eq!(m.label(), "HDSH&HDCH TF");
    }

    #[test]
    fn label_reversed_pair_indices_normalize_order_of_removal() {
        // pair_idx (1,0) → still merges names[0] & names[1], hi=1 is removed first.
        let m = mv(vec!["A", "B"], Some((1, 0)), 1.0);
        // After: parts[0] = "A&B" (because i=1, j=0 → parts[i]&parts[j] = "B&A"?)
        // The code computes `format!("{}&{}", parts[i], parts[j])` then drops `hi`.
        // For (1,0): merged = "B&A", lo=0, hi=1 → remove(1), parts[0]="B&A".
        assert_eq!(m.label(), "B&A");
    }

    #[test]
    fn label_out_of_bounds_pair_idx_is_a_noop() {
        let m = mv(vec!["A"], Some((0, 5)), 1.0);
        assert_eq!(m.label(), "A");
    }

    #[test]
    fn commit_move_appends_single_term_without_parallel_group() {
        let mut active = vec!["IH".to_string()];
        let mut groups: Vec<(String, String)> = Vec::new();
        commit_move(&mut active, &mut groups, &mv(vec!["TF"], None, 1.0));
        assert_eq!(active, vec!["IH", "TF"]);
        assert!(groups.is_empty());
    }

    #[test]
    fn commit_move_appends_pair_and_records_parallel_group() {
        let mut active = vec!["IH".to_string()];
        let mut groups: Vec<(String, String)> = Vec::new();
        commit_move(
            &mut active,
            &mut groups,
            &mv(vec!["HDSH", "HDCH"], Some((0, 1)), 1.0),
        );
        assert_eq!(active, vec!["IH", "HDSH", "HDCH"]);
        assert_eq!(groups, vec![("HDSH".to_string(), "HDCH".to_string())]);
    }

    #[test]
    fn commit_move_pair_idx_references_appended_positions() {
        // pair_idx refers to indices within the newly-appended block, not the
        // full `active` slice. With 2 pre-existing terms and a 3-name move
        // pairing indices (1,2), groups should pick names[1] and names[2].
        let mut active = vec!["IH".to_string(), "ID".to_string()];
        let mut groups: Vec<(String, String)> = Vec::new();
        commit_move(
            &mut active,
            &mut groups,
            &mv(vec!["TF", "HDSH", "HDCH"], Some((1, 2)), 1.0),
        );
        assert_eq!(active, vec!["IH", "ID", "TF", "HDSH", "HDCH"]);
        assert_eq!(groups, vec![("HDSH".to_string(), "HDCH".to_string())]);
    }

    // Lower BIC = better. Smaller RMS at same term count → lower BIC → "Less" ordering.
    #[test]
    fn cmp_by_bic_prefers_smaller_rms_at_same_added_count() {
        let a = mv(vec!["X"], None, 0.5);
        let b = mv(vec!["Y"], None, 2.0);
        assert_eq!(
            cmp_by_bic(&a, &b, 6, 100),
            std::cmp::Ordering::Less,
        );
    }

    // BIC penalizes term count. A pair with the same RMS as a single should
    // sort higher (worse) because it adds two terms.
    #[test]
    fn cmp_by_bic_penalizes_more_added_terms_at_equal_rms() {
        let single = mv(vec!["X"], None, 1.0);
        let pair = mv(vec!["X", "Y"], Some((0, 1)), 1.0);
        assert_eq!(
            cmp_by_bic(&single, &pair, 6, 100),
            std::cmp::Ordering::Less,
        );
    }

    #[test]
    fn cmp_by_bic_equal_inputs_compare_equal() {
        let a = mv(vec!["X"], None, 1.0);
        let b = mv(vec!["Y"], None, 1.0);
        assert_eq!(cmp_by_bic(&a, &b, 6, 100), std::cmp::Ordering::Equal);
    }

    // NaN RMS sorts as Equal (we use unwrap_or(Equal)). This is what
    // protects the rayon min_by chain from panicking.
    #[test]
    fn cmp_by_bic_nan_is_treated_as_equal() {
        let nan_move = mv(vec!["X"], None, f64::NAN);
        let ok_move = mv(vec!["Y"], None, 1.0);
        assert_eq!(
            cmp_by_bic(&nan_move, &ok_move, 6, 100),
            std::cmp::Ordering::Equal,
        );
    }
}
