use crate::parser::{parse_file, planet_name, Variable, Vsop2013File};
use std::collections::HashMap;
use std::path::Path;

pub struct VariableStats {
    pub variable: Variable,
    pub total_terms: usize,
    pub terms_above_threshold: usize,
    pub max_amplitude: f64,
    pub min_amplitude: f64,
    pub time_power_distribution: HashMap<u8, usize>,
}

pub struct PlanetAnalysis {
    pub planet: u8,
    pub planet_name: String,
    pub total_terms: usize,
    pub variable_stats: Vec<VariableStats>,
}

impl PlanetAnalysis {
    pub fn terms_above_threshold(&self) -> usize {
        self.variable_stats
            .iter()
            .map(|v| v.terms_above_threshold)
            .sum()
    }
}

pub fn analyze_file(path: &Path, threshold: f64) -> Result<PlanetAnalysis, String> {
    let vsop = parse_file(path).map_err(|e| format!("Parse error: {}", e))?;
    Ok(analyze_vsop(&vsop, threshold))
}

pub fn analyze_vsop(vsop: &Vsop2013File, threshold: f64) -> PlanetAnalysis {
    let mut variable_stats = Vec::new();

    for var in [
        Variable::A,
        Variable::Lambda,
        Variable::K,
        Variable::H,
        Variable::Q,
        Variable::P,
    ] {
        let blocks: Vec<_> = vsop.blocks_for_variable(var).collect();
        if blocks.is_empty() {
            continue;
        }

        let mut total_terms = 0usize;
        let mut terms_above = 0usize;
        let mut max_amp = 0.0f64;
        let mut min_amp = f64::MAX;
        let mut time_dist: HashMap<u8, usize> = HashMap::new();

        for block in &blocks {
            let count = block.terms.len();
            total_terms += count;
            *time_dist.entry(block.header.time_power).or_insert(0) += count;

            for term in &block.terms {
                let amp = term.amplitude();
                if amp > threshold {
                    terms_above += 1;
                }
                if amp > max_amp {
                    max_amp = amp;
                }
                if amp < min_amp {
                    min_amp = amp;
                }
            }
        }

        variable_stats.push(VariableStats {
            variable: var,
            total_terms,
            terms_above_threshold: terms_above,
            max_amplitude: max_amp,
            min_amplitude: if min_amp == f64::MAX { 0.0 } else { min_amp },
            time_power_distribution: time_dist,
        });
    }

    PlanetAnalysis {
        planet: vsop.planet,
        planet_name: planet_name(vsop.planet).to_string(),
        total_terms: vsop.total_terms(),
        variable_stats,
    }
}

pub fn print_analysis(analysis: &PlanetAnalysis, threshold: f64) {
    println!("\n{}", "=".repeat(70));
    println!(
        "Planet {}: {} ({} total terms)",
        analysis.planet, analysis.planet_name, analysis.total_terms
    );
    println!("{}", "=".repeat(70));

    for stats in &analysis.variable_stats {
        println!(
            "\n  Variable: {} ({})",
            stats.variable,
            stats.variable.name()
        );
        println!("    Total terms: {}", stats.total_terms);
        println!(
            "    Terms above threshold ({:.0e}): {} ({:.1}%)",
            threshold,
            stats.terms_above_threshold,
            (stats.terms_above_threshold as f64 / stats.total_terms as f64) * 100.0
        );
        println!(
            "    Amplitude range: {:.3e} to {:.3e}",
            stats.min_amplitude, stats.max_amplitude
        );

        println!("    Time power distribution:");
        let mut powers: Vec<_> = stats.time_power_distribution.iter().collect();
        powers.sort_by_key(|(k, _)| *k);
        for (power, count) in powers {
            println!("      T^{}: {} terms", power, count);
        }
    }

    println!(
        "\n  Summary: {} terms above threshold ({:.1}%)",
        analysis.terms_above_threshold(),
        (analysis.terms_above_threshold() as f64 / analysis.total_terms as f64) * 100.0
    );
}

pub fn print_summary_table(analyses: &[PlanetAnalysis], threshold: f64) {
    println!("\n{}", "=".repeat(80));
    println!("VSOP2013 Summary (threshold: {:.0e})", threshold);
    println!("{}", "=".repeat(80));
    println!(
        "{:<25} {:>12} {:>12} {:>12} {:>10}",
        "Planet", "Total Terms", "Above Thresh", "Reduction", "%"
    );
    println!("{}", "-".repeat(80));

    let mut grand_total = 0usize;
    let mut grand_above = 0usize;

    for analysis in analyses {
        let above = analysis.terms_above_threshold();
        let reduction = analysis.total_terms - above;
        let pct = (above as f64 / analysis.total_terms as f64) * 100.0;

        grand_total += analysis.total_terms;
        grand_above += above;

        println!(
            "{:<25} {:>12} {:>12} {:>12} {:>9.1}%",
            analysis.planet_name, analysis.total_terms, above, reduction, pct
        );
    }

    println!("{}", "-".repeat(80));
    println!(
        "{:<25} {:>12} {:>12} {:>12} {:>9.1}%",
        "TOTAL",
        grand_total,
        grand_above,
        grand_total - grand_above,
        (grand_above as f64 / grand_total as f64) * 100.0
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{Vsop2013Block, Vsop2013Header, Vsop2013Term};

    fn make_term(s: f64, c: f64) -> Vsop2013Term {
        Vsop2013Term {
            multipliers: [0; 17],
            s_coeff: s,
            c_coeff: c,
        }
    }

    fn make_block(var: Variable, time_power: u8, terms: Vec<Vsop2013Term>) -> Vsop2013Block {
        Vsop2013Block {
            header: Vsop2013Header {
                planet: 3,
                variable: var,
                time_power,
                term_count: terms.len() as u32,
            },
            terms,
        }
    }

    fn make_test_vsop() -> Vsop2013File {
        Vsop2013File {
            planet: 3,
            blocks: vec![
                make_block(
                    Variable::A,
                    0,
                    vec![
                        make_term(3.0, 4.0),   // amplitude = 5
                        make_term(0.6, 0.8),   // amplitude = 1
                        make_term(0.03, 0.04), // amplitude = 0.05
                    ],
                ),
                make_block(
                    Variable::A,
                    1,
                    vec![
                        make_term(6.0, 8.0), // amplitude = 10
                    ],
                ),
                make_block(
                    Variable::Lambda,
                    0,
                    vec![
                        make_term(30.0, 40.0), // amplitude = 50
                        make_term(0.3, 0.4),   // amplitude = 0.5
                    ],
                ),
                make_block(
                    Variable::K,
                    0,
                    vec![
                        make_term(0.12, 0.16), // amplitude = 0.2
                    ],
                ),
            ],
        }
    }

    #[test]
    fn test_variable_stats_creation() {
        let vsop = make_test_vsop();
        let analysis = analyze_vsop(&vsop, 0.1);

        assert_eq!(analysis.planet, 3);
        assert_eq!(analysis.planet_name, "Earth-Moon Barycenter");
        assert_eq!(analysis.total_terms, 7);
        assert_eq!(analysis.variable_stats.len(), 3); // A, Lambda, K
    }

    #[test]
    fn test_analyze_vsop_a_variable() {
        let vsop = make_test_vsop();
        let analysis = analyze_vsop(&vsop, 0.5);

        let a_stats = analysis
            .variable_stats
            .iter()
            .find(|s| s.variable == Variable::A)
            .unwrap();

        assert_eq!(a_stats.total_terms, 4);
        assert_eq!(a_stats.terms_above_threshold, 3); // 5, 1, 10 are above 0.5
        assert!((a_stats.max_amplitude - 10.0).abs() < 1e-10);
        assert!((a_stats.min_amplitude - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_analyze_vsop_time_power_distribution() {
        let vsop = make_test_vsop();
        let analysis = analyze_vsop(&vsop, 0.0);

        let a_stats = analysis
            .variable_stats
            .iter()
            .find(|s| s.variable == Variable::A)
            .unwrap();

        assert_eq!(a_stats.time_power_distribution.get(&0), Some(&3));
        assert_eq!(a_stats.time_power_distribution.get(&1), Some(&1));
    }

    #[test]
    fn test_terms_above_threshold_aggregation() {
        let vsop = make_test_vsop();
        let analysis = analyze_vsop(&vsop, 0.3);

        let total_above = analysis.terms_above_threshold();
        // A: 5, 1, 10 above 0.3 = 3
        // Lambda: 50, 0.5 above 0.3 = 2
        // K: 0.2 not above 0.3 = 0
        assert_eq!(total_above, 5);
    }

    #[test]
    fn test_analyze_vsop_empty_variable() {
        let vsop = Vsop2013File {
            planet: 1,
            blocks: vec![make_block(Variable::A, 0, vec![make_term(1.0, 0.0)])],
        };
        let analysis = analyze_vsop(&vsop, 0.0);

        // Only A should be in stats (Lambda, K, etc. have no blocks)
        assert_eq!(analysis.variable_stats.len(), 1);
        assert_eq!(analysis.variable_stats[0].variable, Variable::A);
    }

    #[test]
    fn test_analyze_vsop_min_amplitude_with_empty_blocks() {
        let vsop = Vsop2013File {
            planet: 2,
            blocks: vec![],
        };
        let analysis = analyze_vsop(&vsop, 0.0);

        assert_eq!(analysis.variable_stats.len(), 0);
        assert_eq!(analysis.total_terms, 0);
    }

    #[test]
    fn test_print_analysis_runs_without_panic() {
        let vsop = make_test_vsop();
        let analysis = analyze_vsop(&vsop, 1e-5);
        print_analysis(&analysis, 1e-5);
    }

    #[test]
    fn test_print_analysis_with_high_threshold() {
        let vsop = make_test_vsop();
        let analysis = analyze_vsop(&vsop, 1000.0);
        print_analysis(&analysis, 1000.0);
    }

    #[test]
    fn test_print_summary_table_single_planet() {
        let vsop = make_test_vsop();
        let analysis = analyze_vsop(&vsop, 1.0);
        print_summary_table(&[analysis], 1.0);
    }

    #[test]
    fn test_print_summary_table_multiple_planets() {
        let vsop1 = Vsop2013File {
            planet: 1,
            blocks: vec![make_block(Variable::A, 0, vec![make_term(10.0, 0.0)])],
        };
        let vsop2 = Vsop2013File {
            planet: 2,
            blocks: vec![make_block(Variable::A, 0, vec![make_term(5.0, 0.0)])],
        };

        let analysis1 = analyze_vsop(&vsop1, 1.0);
        let analysis2 = analyze_vsop(&vsop2, 1.0);

        print_summary_table(&[analysis1, analysis2], 1.0);
    }

    #[test]
    fn test_planet_analysis_struct() {
        let analysis = PlanetAnalysis {
            planet: 5,
            planet_name: "Jupiter".to_string(),
            total_terms: 100,
            variable_stats: vec![
                VariableStats {
                    variable: Variable::A,
                    total_terms: 50,
                    terms_above_threshold: 30,
                    max_amplitude: 1.0,
                    min_amplitude: 0.001,
                    time_power_distribution: HashMap::new(),
                },
                VariableStats {
                    variable: Variable::Lambda,
                    total_terms: 50,
                    terms_above_threshold: 20,
                    max_amplitude: 0.5,
                    min_amplitude: 0.0001,
                    time_power_distribution: HashMap::new(),
                },
            ],
        };

        assert_eq!(analysis.terms_above_threshold(), 50);
    }

    #[test]
    fn test_variable_stats_struct() {
        let mut dist = HashMap::new();
        dist.insert(0u8, 100usize);
        dist.insert(1u8, 50usize);

        let stats = VariableStats {
            variable: Variable::K,
            total_terms: 150,
            terms_above_threshold: 75,
            max_amplitude: 0.01,
            min_amplitude: 1e-10,
            time_power_distribution: dist,
        };

        assert_eq!(stats.variable, Variable::K);
        assert_eq!(stats.total_terms, 150);
        assert_eq!(stats.time_power_distribution.len(), 2);
    }

    #[test]
    #[ignore = "requires local VSOP2013 data files"]
    fn test_analyze_earth() {
        let path = Path::new("references/ephemeris/vsop2013/VSOP2013p3.dat");
        let analysis = analyze_file(path, 1e-10).unwrap();
        assert_eq!(analysis.planet, 3);
        assert!(!analysis.variable_stats.is_empty());
        assert!(analysis.total_terms > 100_000);
    }

    #[test]
    fn test_analyze_file_with_mock_data() {
        use tempfile::TempDir;

        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("test_vsop.dat");

        let content = " VSOP2013  5  1  0  2    JUPITER VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
    2   1  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.2000000000000000 +00  0.0000000000000000 +00
";
        std::fs::write(&file_path, content).unwrap();

        let analysis = analyze_file(&file_path, 0.0).unwrap();
        assert_eq!(analysis.planet, 5);
        assert_eq!(analysis.planet_name, "Jupiter");
        assert_eq!(analysis.total_terms, 2);
        assert!(!analysis.variable_stats.is_empty());
    }

    #[test]
    fn test_analyze_file_parse_error() {
        use tempfile::TempDir;

        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("bad_vsop.dat");

        // Create a file with an invalid VSOP2013 header that will fail parsing
        // The header claims 100 terms but we provide none - this triggers MissingTerms error
        std::fs::write(&file_path, " VSOP2013  3  1  0  100    BAD\n").unwrap();

        let result = analyze_file(&file_path, 0.0);
        match result {
            Err(msg) => assert!(
                msg.contains("Parse error"),
                "Expected Parse error, got: {}",
                msg
            ),
            Ok(_) => panic!("Expected an error"),
        }
    }

    #[test]
    fn test_analyze_file_not_found() {
        let result = analyze_file(Path::new("/nonexistent/path/file.dat"), 0.0);
        match result {
            Err(msg) => assert!(
                msg.contains("Parse error"),
                "Expected Parse error, got: {}",
                msg
            ),
            Ok(_) => panic!("Expected an error"),
        }
    }
}
