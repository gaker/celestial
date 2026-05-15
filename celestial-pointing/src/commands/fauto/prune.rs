use super::bic::{bic_of, compute_bic, fit_is_well_conditioned, try_fit};
use super::BASE_TERMS;
use crate::error::Result;
use crate::observation::Observation;
use rayon::prelude::*;

pub(super) const MIN_SIGNIFICANCE: f64 = 2.0;

pub(super) fn backward_eliminate(
    observations: &[&Observation],
    active: &mut Vec<String>,
    parallel_groups: &mut Vec<(String, String)>,
    latitude: f64,
    fit_tol: f64,
) -> Result<Vec<String>> {
    let mut removed = Vec::new();
    let mut current_bic = bic_of(observations, active, latitude, fit_tol)?;
    loop {
        let groups = elimination_groups(active, parallel_groups);
        let drop = find_best_elimination(
            observations,
            active,
            &groups,
            current_bic,
            latitude,
            fit_tol,
        );
        let Some((names, bic)) = drop else { break };
        for name in &names {
            active.retain(|n| n != name);
            removed.push(name.clone());
        }
        parallel_groups.retain(|(a, b)| !names.contains(a) && !names.contains(b));
        current_bic = bic;
    }
    Ok(removed)
}

pub(super) fn elimination_groups(
    active: &[String],
    parallel_groups: &[(String, String)],
) -> Vec<Vec<String>> {
    let base_set: std::collections::HashSet<&str> = BASE_TERMS.iter().copied().collect();
    let pair_lookup: std::collections::HashMap<&str, &str> = parallel_groups
        .iter()
        .flat_map(|(a, b)| [(a.as_str(), b.as_str()), (b.as_str(), a.as_str())])
        .collect();
    let mut groups = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for name in active {
        if base_set.contains(name.as_str()) || seen.contains(name.as_str()) {
            continue;
        }
        match pair_lookup.get(name.as_str()) {
            Some(&partner) if active.iter().any(|n| n == partner) => {
                seen.insert(name.as_str());
                seen.insert(partner);
                groups.push(vec![name.clone(), partner.to_string()]);
            }
            _ => {
                seen.insert(name.as_str());
                groups.push(vec![name.clone()]);
            }
        }
    }
    groups
}

fn find_best_elimination(
    observations: &[&Observation],
    active: &[String],
    groups: &[Vec<String>],
    current_bic: f64,
    latitude: f64,
    fit_tol: f64,
) -> Option<(Vec<String>, f64)> {
    let n_obs = observations.len();
    groups
        .par_iter()
        .filter_map(|group| {
            let trial: Vec<String> = active
                .iter()
                .filter(|n| !group.contains(n))
                .cloned()
                .collect();
            let fit = try_fit(observations, &trial, latitude, fit_tol).ok()?;
            if !fit_is_well_conditioned(&fit) {
                return None;
            }
            Some((group.clone(), compute_bic(n_obs, trial.len(), fit.sky_rms)))
        })
        .filter(|(_, bic)| *bic <= current_bic)
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
}

pub(super) fn prune_terms(
    observations: &[&Observation],
    active: &mut Vec<String>,
    latitude: f64,
    fit_tol: f64,
) -> Result<Vec<String>> {
    let fit = try_fit(observations, active, latitude, fit_tol)?;
    let base_set: std::collections::HashSet<&str> = BASE_TERMS.iter().copied().collect();
    let to_remove: Vec<String> = active
        .iter()
        .enumerate()
        .filter(|(i, name)| {
            if base_set.contains(name.as_str()) {
                return false;
            }
            let sigma = fit.sigma[*i];
            sigma > 0.0 && (fit.coefficients[*i] / sigma).abs() < MIN_SIGNIFICANCE
        })
        .map(|(_, n)| n.clone())
        .collect();
    let mut pruned = Vec::new();
    for name in &to_remove {
        active.retain(|n| n != name);
        pruned.push(name.clone());
    }
    Ok(pruned)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(strs: &[&str]) -> Vec<String> {
        strs.iter().map(|x| x.to_string()).collect()
    }

    fn pair(a: &str, b: &str) -> (String, String) {
        (a.to_string(), b.to_string())
    }

    // Base terms (IH, ID, CH, NP, MA, ME) are never proposed for elimination.
    #[test]
    fn elimination_groups_excludes_base_terms() {
        let active = s(&["IH", "ID", "CH", "NP", "MA", "ME"]);
        let groups = elimination_groups(&active, &[]);
        assert!(groups.is_empty(), "got {:?}", groups);
    }

    #[test]
    fn elimination_groups_emits_singleton_for_unpaired_extra_term() {
        let active = s(&["IH", "ID", "CH", "NP", "MA", "ME", "TF"]);
        let groups = elimination_groups(&active, &[]);
        assert_eq!(groups, vec![s(&["TF"])]);
    }

    #[test]
    fn elimination_groups_pairs_partners_together() {
        let active = s(&["IH", "ID", "CH", "NP", "MA", "ME", "HDSH", "HDCH"]);
        let groups = elimination_groups(&active, &[pair("HDSH", "HDCH")]);
        // Should emit one group of two, not two groups of one.
        assert_eq!(groups.len(), 1);
        let g = &groups[0];
        assert!(g.contains(&"HDSH".to_string()));
        assert!(g.contains(&"HDCH".to_string()));
    }

    // If a pair was registered but its partner was already removed,
    // the remaining one is treated as a singleton group.
    #[test]
    fn elimination_groups_treats_orphan_pair_member_as_singleton() {
        let active = s(&["IH", "ID", "CH", "NP", "MA", "ME", "HDSH"]);
        let groups = elimination_groups(&active, &[pair("HDSH", "HDCH")]);
        assert_eq!(groups, vec![s(&["HDSH"])]);
    }

    #[test]
    fn elimination_groups_handles_mixed_pairs_and_singletons() {
        let active = s(&[
            "IH", "ID", "CH", "NP", "MA", "ME", "TF", "HDSH", "HDCH", "DAF",
        ]);
        let groups = elimination_groups(&active, &[pair("HDSH", "HDCH")]);
        // Order follows active-list order: TF, then HDSH/HDCH pair, then DAF.
        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0], s(&["TF"]));
        assert!(groups[1].contains(&"HDSH".to_string()));
        assert!(groups[1].contains(&"HDCH".to_string()));
        assert_eq!(groups[2], s(&["DAF"]));
    }

    #[test]
    fn elimination_groups_deduplicates_when_both_pair_members_seen() {
        // Even with multiple pair declarations, each name should only appear once.
        let active = s(&["IH", "ID", "CH", "NP", "MA", "ME", "HDSH", "HDCH"]);
        let groups = elimination_groups(
            &active,
            &[pair("HDSH", "HDCH"), pair("HDCH", "HDSH")],
        );
        assert_eq!(groups.len(), 1);
    }

    #[test]
    fn min_significance_threshold_is_two_sigma() {
        // Pinning the constant so any change is intentional.
        assert_eq!(MIN_SIGNIFICANCE, 2.0);
    }
}
