use super::*;
use crate::parser::{Vsop2013Block, Vsop2013Header, Vsop2013Term};
use tempfile::TempDir;

fn block(variable: Variable, time_power: u8, terms: &[(f64, f64)]) -> Vsop2013Block {
    Vsop2013Block {
        header: Vsop2013Header {
            planet: 3,
            variable,
            time_power,
            term_count: terms.len() as u32,
        },
        terms: terms
            .iter()
            .map(|&(s_coeff, c_coeff)| Vsop2013Term {
                multipliers: [0; 17],
                s_coeff,
                c_coeff,
            })
            .collect(),
    }
}

// Amplitudes: A 5, 1, 0.05 at T^0 and 10 at T^1; Lambda 50 and 0.5; K 0.2.
fn sample() -> Vsop2013File {
    Vsop2013File {
        planet: 3,
        blocks: vec![
            block(Variable::A, 0, &[(3.0, 4.0), (0.6, 0.8), (0.03, 0.04)]),
            block(Variable::A, 1, &[(6.0, 8.0)]),
            block(Variable::Lambda, 0, &[(30.0, 40.0), (0.3, 0.4)]),
            block(Variable::K, 0, &[(0.12, 0.16)]),
        ],
    }
}

#[test]
fn planet_totals() {
    let analysis = analyze_vsop(&sample(), 0.3);
    assert_eq!(analysis.planet, 3);
    assert_eq!(analysis.planet_name, "Earth-Moon Barycenter");
    assert_eq!(analysis.total_terms, 7);
    let variables: Vec<Variable> = analysis.variable_stats.iter().map(|s| s.variable).collect();
    assert_eq!(variables, [Variable::A, Variable::Lambda, Variable::K]);
    assert_eq!(analysis.terms_above_threshold(), 5);
}

#[test]
fn variable_statistics() {
    let analysis = analyze_vsop(&sample(), 0.5);
    let a = &analysis.variable_stats[0];
    assert_eq!((a.total_terms, a.terms_above_threshold), (4, 3));
    assert_eq!(a.max_amplitude, 10.0);
    assert_eq!(a.min_amplitude, libm::sqrt(0.03 * 0.03 + 0.04 * 0.04));
    assert_eq!(a.time_power_distribution, BTreeMap::from([(0, 3), (1, 1)]));
}

#[test]
fn empty_file_has_no_statistics() {
    let vsop = Vsop2013File {
        planet: 2,
        blocks: vec![],
    };
    let analysis = analyze_vsop(&vsop, 0.0);
    assert!(analysis.variable_stats.is_empty());
    assert_eq!(analysis.total_terms, 0);
    assert!(format_analysis(&analysis).ends_with("Summary: 0 terms above threshold (0.0%)\n"));
}

#[test]
fn empty_blocks_report_zero_amplitudes() {
    let vsop = Vsop2013File {
        planet: 1,
        blocks: vec![block(Variable::Q, 2, &[])],
    };
    let stats = &analyze_vsop(&vsop, 0.0).variable_stats[0];
    assert_eq!((stats.min_amplitude, stats.max_amplitude), (0.0, 0.0));
    assert_eq!(stats.time_power_distribution, BTreeMap::from([(2, 0)]));
}

#[test]
fn report() {
    let vsop = Vsop2013File {
        planet: 5,
        blocks: vec![
            block(Variable::A, 0, &[(3.0, 4.0), (0.6, 0.8)]),
            block(Variable::A, 2, &[(6.0, 8.0)]),
        ],
    };
    let rule = "=".repeat(70);
    let expected = format!(
        "\n{rule}\nPlanet 5: Jupiter (3 total terms)\n{rule}\n\
         \n  Variable: A (A (semi-major axis))\n    Total terms: 3\n    \
         Terms above threshold (2e0): 2 (66.7%)\n    \
         Amplitude range: 1.000e0 to 1.000e1\n    Time power distribution:\n      \
         T^0: 2 terms\n      T^2: 1 terms\n\
         \n  Summary: 2 terms above threshold (66.7%)\n"
    );
    assert_eq!(format_analysis(&analyze_vsop(&vsop, 2.0)), expected);
}

#[test]
fn summary_table() {
    let one = |planet, s| Vsop2013File {
        planet,
        blocks: vec![block(Variable::A, 0, &[(s, 0.0), (0.5, 0.0)])],
    };
    let analyses = [
        analyze_vsop(&one(1, 10.0), 1.0),
        analyze_vsop(&one(2, 5.0), 1e-7),
    ];
    let table = format_summary_table(&analyses);
    let lines: Vec<&str> = table.lines().collect();
    assert_eq!(lines[0], "");
    assert_eq!(lines[2], "VSOP2013 Summary");
    assert_eq!(
        lines[4],
        "Planet                    Threshold  Total Terms Above Thresh    Reduction      %"
    );
    assert_eq!(
        lines[6],
        "Mercury                         1e0            2            1            1  50.0%"
    );
    assert_eq!(
        lines[7],
        "Venus                          1e-7            2            2            0 100.0%"
    );
    assert_eq!(
        lines[9],
        "TOTAL                                          4            3            1  75.0%"
    );
    assert_eq!(lines.len(), 10);
}

#[test]
fn analyzes_a_file() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("test_vsop.dat");
    let content = " VSOP2013  5  1  0  2    JUPITER VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
    2   1  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.2000000000000000 +00  0.0000000000000000 +00
";
    std::fs::write(&path, content).unwrap();
    let analysis = analyze_file(&path, 0.0).unwrap();
    assert_eq!((analysis.planet, analysis.total_terms), (5, 2));
    assert_eq!(analysis.planet_name, "Jupiter");
    assert_eq!(analysis.variable_stats.len(), 1);
}

#[test]
fn file_errors_are_parse_errors() {
    let dir = TempDir::new().unwrap();
    let short = dir.path().join("short.dat");
    std::fs::write(&short, " VSOP2013  3  1  0  100    BAD\n").unwrap();
    assert_eq!(
        analyze_file(&short, 0.0).err(),
        Some("Parse error: Expected 100 terms, found 0".to_string())
    );
    let missing = analyze_file(Path::new("/nonexistent/path/file.dat"), 0.0)
        .err()
        .unwrap();
    assert!(missing.starts_with("Parse error: IO error: "), "{missing}");
}

#[test]
#[ignore = "requires local VSOP2013 data files"]
fn shipped_emb_file() {
    let path = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../vsop2013/VSOP2013p3.dat"
    ));
    let analysis = analyze_file(path, 1e-10).unwrap();
    assert_eq!(analysis.planet, 3);
    assert_eq!(analysis.variable_stats.len(), 6);
    assert_eq!(analysis.total_terms, 294426);
}
