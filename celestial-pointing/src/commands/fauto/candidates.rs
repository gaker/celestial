const PHYSICAL_CANDIDATES: &[&str] = &[
    "TF", "TX", "DAF", "FO", "HCES", "HCEC", "DCES", "DCEC",
];
const SINGLE_HARMONIC_MAX_FREQ: u8 = 8;
const CROSS_HARMONIC_MAX_HA_FREQ: u8 = 4;
const CROSS_HARMONIC_MAX_DEC_FREQ: u8 = 8;

pub(super) fn physical_candidates() -> &'static [&'static str] {
    PHYSICAL_CANDIDATES
}

#[derive(Clone)]
pub(super) struct CandidatePool {
    pub singles: Vec<String>,
    pub pairs: Vec<(String, String)>,
    pub triples: Vec<(String, String, String)>,
}

pub(super) fn build_candidate_pool() -> CandidatePool {
    let mut singles = Vec::new();
    for &p in PHYSICAL_CANDIDATES {
        singles.push(p.to_string());
    }
    let single_h = single_component_harmonics();
    let cross_h = cross_harmonics();
    let (s_only, pair_singles, pair_list) = split_singles_and_pairs(&single_h);
    let (cs_only, cpair_singles, cpair_list) = split_singles_and_pairs(&cross_h);
    let all_harmonic_singles: Vec<String> = s_only
        .iter()
        .chain(cs_only.iter())
        .chain(pair_singles.iter())
        .chain(cpair_singles.iter())
        .cloned()
        .collect();
    singles.extend(s_only);
    singles.extend(cs_only);
    singles.extend(pair_singles);
    singles.extend(cpair_singles);
    let mut pairs = pair_list;
    pairs.extend(cpair_list);
    let triples = build_triples(&pairs, &all_harmonic_singles);
    CandidatePool {
        singles,
        pairs,
        triples,
    }
}

pub(super) fn filter_remaining(
    candidates: &CandidatePool,
    active: &[String],
) -> CandidatePool {
    let active_set: std::collections::HashSet<&str> =
        active.iter().map(|s| s.as_str()).collect();
    CandidatePool {
        singles: candidates
            .singles
            .iter()
            .filter(|c| !active_set.contains(c.as_str()))
            .cloned()
            .collect(),
        pairs: candidates
            .pairs
            .iter()
            .filter(|(a, b)| {
                !active_set.contains(a.as_str()) && !active_set.contains(b.as_str())
            })
            .cloned()
            .collect(),
        triples: candidates
            .triples
            .iter()
            .filter(|(a, b, c)| {
                !active_set.contains(a.as_str())
                    && !active_set.contains(b.as_str())
                    && !active_set.contains(c.as_str())
            })
            .cloned()
            .collect(),
    }
}

fn single_component_harmonics() -> Vec<String> {
    let results = ["H", "D", "X", "P"];
    let funcs = ["S", "C"];
    let coords = ["H", "D"];
    let mut out = Vec::new();
    for r in &results {
        for c in &coords {
            for n in 1..=SINGLE_HARMONIC_MAX_FREQ {
                for f in &funcs {
                    out.push(format!("H{}{}{}{}", r, f, c, freq_suffix(n)));
                }
            }
        }
    }
    out
}

fn cross_harmonics() -> Vec<String> {
    let results = ["H", "D", "X"];
    let funcs = ["S", "C"];
    let mut out = Vec::new();
    for r in &results {
        for nh in 1..=CROSS_HARMONIC_MAX_HA_FREQ {
            for nd in 1..=CROSS_HARMONIC_MAX_DEC_FREQ {
                for fd in &funcs {
                    for fh in &funcs {
                        out.push(format!(
                            "H{}{}H{}{}D{}",
                            r,
                            fh,
                            freq_suffix(nh),
                            fd,
                            freq_suffix(nd),
                        ));
                    }
                }
            }
        }
    }
    out
}

fn split_singles_and_pairs(
    names: &[String],
) -> (Vec<String>, Vec<String>, Vec<(String, String)>) {
    let mut singletons = Vec::new();
    let mut pair_singles = Vec::new();
    let mut pairs = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for name in names {
        if seen.contains(name) {
            continue;
        }
        match pair_partner(name) {
            Some(partner)
                if names.iter().any(|n| n == &partner) && !seen.contains(&partner) =>
            {
                seen.insert(name.clone());
                seen.insert(partner.clone());
                pair_singles.push(name.clone());
                pair_singles.push(partner.clone());
                pairs.push((name.clone(), partner));
            }
            _ => {
                seen.insert(name.clone());
                singletons.push(name.clone());
            }
        }
    }
    (singletons, pair_singles, pairs)
}

pub(super) fn pair_partner(name: &str) -> Option<String> {
    let bytes = name.as_bytes();
    if bytes.len() < 4 || bytes[0] != b'H' {
        return None;
    }
    let func_idx = 2;
    let flipped = match bytes[func_idx] {
        b'S' => b'C',
        b'C' => b'S',
        _ => return None,
    };
    let mut out = bytes.to_vec();
    out[func_idx] = flipped;
    String::from_utf8(out).ok()
}

fn freq_suffix(n: u8) -> String {
    if n == 1 {
        String::new()
    } else {
        n.to_string()
    }
}

fn build_triples(
    pairs: &[(String, String)],
    harmonic_singles: &[String],
) -> Vec<(String, String, String)> {
    let mut out = Vec::new();
    for (a, b) in pairs {
        let Some(pair_sig) = harmonic_signature(a) else {
            continue;
        };
        for third in harmonic_singles {
            if third == a || third == b {
                continue;
            }
            let Some(t_sig) = harmonic_signature(third) else {
                continue;
            };
            if !triple_companions(&pair_sig, &t_sig) {
                continue;
            }
            out.push((a.clone(), b.clone(), third.clone()));
        }
    }
    out
}

struct HarmonicSig {
    result: char,
    freqs: Vec<u8>,
}

fn harmonic_signature(name: &str) -> Option<HarmonicSig> {
    let bytes = name.as_bytes();
    if bytes.is_empty() || bytes[0] != b'H' {
        return None;
    }
    let result = bytes.get(1).copied()? as char;
    let mut freqs = Vec::new();
    let mut i = 2;
    while i < bytes.len() {
        if i + 1 >= bytes.len() {
            return None;
        }
        i += 2;
        let freq = if i < bytes.len() && bytes[i].is_ascii_digit() {
            let d = bytes[i] - b'0';
            i += 1;
            d
        } else {
            1
        };
        freqs.push(freq);
    }
    Some(HarmonicSig { result, freqs })
}

fn triple_companions(pair: &HarmonicSig, third: &HarmonicSig) -> bool {
    if pair.result != third.result {
        return false;
    }
    third.freqs.iter().any(|f| pair.freqs.contains(f))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pair_partner_hdsh_is_hdch() {
        assert_eq!(pair_partner("HDSH"), Some("HDCH".to_string()));
        assert_eq!(pair_partner("HDCH"), Some("HDSH".to_string()));
    }

    #[test]
    fn pair_partner_cross_harmonic() {
        assert_eq!(pair_partner("HDSH2CD3"), Some("HDCH2CD3".to_string()));
    }

    #[test]
    fn pair_partner_non_harmonic_is_none() {
        assert_eq!(pair_partner("IH"), None);
        assert_eq!(pair_partner("TF"), None);
    }

    #[test]
    fn candidate_pool_contains_known_pairs() {
        let pool = build_candidate_pool();
        let has_hdsh_hdch = pool
            .pairs
            .iter()
            .any(|(a, b)| (a == "HDSH" && b == "HDCH") || (a == "HDCH" && b == "HDSH"));
        assert!(has_hdsh_hdch, "HDSH/HDCH pair should be enumerated");
    }

    #[test]
    fn candidate_pool_has_physical_singles() {
        let pool = build_candidate_pool();
        for name in PHYSICAL_CANDIDATES {
            assert!(
                pool.singles.iter().any(|s| s == name),
                "missing physical candidate {}",
                name,
            );
        }
    }

    #[test]
    fn triples_have_matching_result_code() {
        let pool = build_candidate_pool();
        for (a, _b, c) in &pool.triples {
            let sa = harmonic_signature(a).unwrap();
            let sc = harmonic_signature(c).unwrap();
            assert_eq!(sa.result, sc.result, "{} {} should share result code", a, c);
        }
    }
}
