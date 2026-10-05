use super::source::*;
use super::*;
use crate::parser::{MainSeries, PertSeries};
use std::path::Path;
use tempfile::TempDir;
use Coordinate::{Distance, Latitude, Longitude};

fn main_term(delaunay: [i32; 4], a0: f64) -> MainTerm {
    MainTerm {
        delaunay,
        coeffs: [a0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    }
}

fn pert_term(sin_coeff: f64, cos_coeff: f64, multipliers: [i32; 16]) -> PertTerm {
    PertTerm {
        sin_coeff,
        cos_coeff,
        multipliers,
    }
}

fn block(time_power: u8, terms: Vec<PertTerm>) -> PertBlock {
    PertBlock { time_power, terms }
}

const COUNTING: [i32; 16] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];

// At a cut of 1: an empty main series, an empty block, and terms either side
// of the cut.
fn sample() -> ElpData {
    let main = |coordinate, terms| MainSeries { coordinate, terms };
    let pert = |coordinate, blocks| PertSeries { coordinate, blocks };
    ElpData {
        main: [
            main(Longitude, vec![main_term([1, -2, 3, -4], 100.0)]),
            main(Latitude, vec![]),
            main(Distance, vec![main_term([0, 0, 0, 0], 0.5)]),
        ],
        pert: [
            pert(
                Longitude,
                vec![
                    block(0, vec![pert_term(0.0, 2.0, COUNTING)]),
                    block(1, vec![pert_term(0.3, 0.4, [0; 16])]),
                ],
            ),
            pert(Latitude, vec![block(3, vec![])]),
            pert(
                Distance,
                vec![block(2, vec![pert_term(-2.0, 0.0, [0; 16])])],
            ),
        ],
    }
}

fn empty() -> ElpData {
    let mut elp = sample();
    for series in &mut elp.main {
        series.terms.clear();
    }
    for series in &mut elp.pert {
        series.blocks.clear();
    }
    elp
}

fn config(threshold: f64, output_dir: &Path) -> GenerateConfig {
    GenerateConfig {
        threshold,
        output_dir: output_dir.to_path_buf(),
    }
}

fn generate(elp: &ElpData, threshold: f64) -> (String, String) {
    let dir = TempDir::new().unwrap();
    let report = generate_moon_module(elp, &config(threshold, dir.path())).unwrap();
    let source = std::fs::read_to_string(dir.path().join("moon.rs")).unwrap();
    let report = report.replace(&dir.path().display().to_string(), "<dir>");
    (report, source)
}

#[test]
fn coordinate_names() {
    assert_eq!(coord_name(Longitude), "LONGITUDE");
    assert_eq!(coord_name(Latitude), "LATITUDE");
    assert_eq!(coord_name(Distance), "DISTANCE");
}

#[test]
fn cuts_keep_file_order_and_the_threshold_itself() {
    let series = MainSeries {
        coordinate: Latitude,
        terms: [10.0, 100.0, 50.0, -20.0]
            .map(|a0| main_term([0; 4], a0))
            .to_vec(),
    };
    let kept: Vec<f64> = (filter_main_terms(&series, 20.0).iter())
        .map(|t| t.coeffs[0])
        .collect();
    assert_eq!(kept, [100.0, 50.0, -20.0]);
    let terms = [(3.0, 4.0), (30.0, 40.0), (6.0, 8.0)].map(|(s, c)| pert_term(s, c, [0; 16]));
    let block = block(1, terms.to_vec());
    let kept: Vec<f64> = (filter_pert_terms(&block, 10.0).iter())
        .map(|t| t.amplitude())
        .collect();
    assert_eq!(kept, [50.0, 10.0]);
}

#[test]
fn floats_are_shortest_round_trip() {
    let cases = [
        (0.0, "0.0"),
        (1.5, "1.5e0"),
        (-411.60287, "-4.1160287e2"),
        (0.4, "4e-1"),
        (1e-300, "1e-300"),
    ];
    for (value, text) in cases {
        assert_eq!(format_float(value), text);
        assert_eq!(text.parse::<f64>().unwrap(), value);
    }
    assert_eq!(format_floats(&[0.0, -0.18]), "0.0, -1.8e-1");
    assert_eq!(format_ints(&[0, -18, 16]), "0, -18, 16");
}

#[test]
fn literals_round_trip() {
    let mut elp = sample();
    let coeffs = [-411.60287, 168.48, -18433.81, -121.62, 0.4, -0.18, 0.0];
    elp.main[0].terms = vec![MainTerm {
        delaunay: [0, 2, 0, 0],
        coeffs,
    }];
    let term = pert_term(-0.1274921554086e2, 0.6368794709728e1, [0; 16]);
    elp.pert[0].blocks[0].terms = vec![term.clone()];
    let (_, source) = generate(&elp, 0.0);
    assert!(source
        .contains("coeffs: [-4.1160287e2, 1.6848e2, -1.843381e4, -1.2162e2, 4e-1, -1.8e-1, 0.0]"));
    let line = (source.lines())
        .find(|l| l.contains("PertTerm { amplitude: 1.4"))
        .unwrap();
    let field = |name: &str| -> f64 {
        let rest = &line[line.find(name).unwrap() + name.len()..];
        rest[..rest.find(',').unwrap()].parse().unwrap()
    };
    assert_eq!(
        (field("amplitude: "), field("phase: ")),
        (term.amplitude(), term.phase())
    );
}

#[test]
fn header_records_the_command_and_no_date() {
    let lines: Vec<String> = (header(5.0, 12, 13).lines()).map(String::from).collect();
    assert_eq!(
        lines,
        [
            "//! ELP/MPP02 coefficients for the Moon",
            "//!",
            "//! Generated from ELP_MAIN.S1-S3 and ELP_PERT.S1-S3 by",
            "//! `elpmpp02-gen generate --input <dir> --output <dir> --threshold 5e0`",
            "//! Terms retained: 12 of 13 (92.3%)",
            "//!",
            "//! Reference: Chapront & Francou (2003)",
            "//! \"The lunar theory ELP revisited. Introduction of new planetary perturbations\"",
            "//! Astronomy & Astrophysics, 404, 735-742",
            "",
        ]
    );
    assert!(header(1e-3, 0, 0).contains("Terms retained: 0 of 0 (0.0%)"));
}

const SAMPLE_BODY: &str = "\
/// Main problem terms for LONGITUDE (1 of 1 terms)
pub(crate) const MAIN_LONGITUDE: &[MainTerm] = &[
    MainTerm { delaunay: [1, -2, 3, -4], coeffs: [1e2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] },
];

/// Main problem terms for LATITUDE (0 of 0 terms)
pub(crate) const MAIN_LATITUDE: &[MainTerm] = &[
];

/// Main problem terms for DISTANCE (0 of 1 terms)
pub(crate) const MAIN_DISTANCE: &[MainTerm] = &[
];

/// Perturbation terms for LONGITUDE T^0 (1 of 1 terms)
const PERT_LONGITUDE_T0: &[PertTerm] = &[
    PertTerm { amplitude: 2e0, phase: 1.5707963267948966e0, multipliers: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] },
];

/// Perturbation terms for LONGITUDE T^1 (0 of 1 terms)
const PERT_LONGITUDE_T1: &[PertTerm] = &[
];

/// All perturbation blocks for LONGITUDE
pub(crate) const PERT_LONGITUDE: &[PertBlock] = &[
    PertBlock { power: 0, terms: PERT_LONGITUDE_T0 },
    PertBlock { power: 1, terms: PERT_LONGITUDE_T1 },
];

/// Perturbation terms for LATITUDE T^3 (0 of 0 terms)
const PERT_LATITUDE_T3: &[PertTerm] = &[
];

/// All perturbation blocks for LATITUDE
pub(crate) const PERT_LATITUDE: &[PertBlock] = &[
    PertBlock { power: 3, terms: PERT_LATITUDE_T3 },
];

/// Perturbation terms for DISTANCE T^2 (1 of 1 terms)
const PERT_DISTANCE_T2: &[PertTerm] = &[
    PertTerm { amplitude: 2e0, phase: 3.141592653589793e0, multipliers: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] },
];

/// All perturbation blocks for DISTANCE
pub(crate) const PERT_DISTANCE: &[PertBlock] = &[
    PertBlock { power: 2, terms: PERT_DISTANCE_T2 },
];

";

#[test]
fn moon_module() {
    let (report, source) = generate(&sample(), 1.0);
    assert_eq!(
        report,
        "Generated <dir>/moon.rs (3 of 5 terms, 60.0%, threshold 1e0)"
    );
    assert_eq!(source, header(1.0, 3, 5) + LINT + TYPES + SAMPLE_BODY);
    assert_eq!(moon_source(&sample(), 1.0), (source, 3, 5));
}

#[test]
fn empty_data_reports_zero_percent() {
    let (report, _) = generate(&empty(), 1e-3);
    assert_eq!(
        report,
        "Generated <dir>/moon.rs (0 of 0 terms, 0.0%, threshold 1e-3)"
    );
}

#[test]
fn output_dirs_are_created_and_files_replaced() {
    let dir = TempDir::new().unwrap();
    let output_dir = dir.path().join("a/b");
    std::fs::create_dir_all(&output_dir).unwrap();
    std::fs::write(output_dir.join("moon.rs"), "old content").unwrap();
    generate_moon_module(&sample(), &config(1.0, &output_dir)).unwrap();
    generate_moon_module(&sample(), &config(1.0, &dir.path().join("c/d"))).unwrap();
    for path in [output_dir.join("moon.rs"), dir.path().join("c/d/moon.rs")] {
        let source = std::fs::read_to_string(path).unwrap();
        assert_eq!(source, moon_source(&sample(), 1.0).0);
    }
}

#[test]
fn write_errors_name_the_step() {
    let dir = TempDir::new().unwrap();
    let file = dir.path().join("file");
    std::fs::write(&file, "").unwrap();
    let error = generate_moon_module(&sample(), &config(1.0, &file)).unwrap_err();
    assert!(
        error.starts_with("Failed to create output directory: "),
        "{error}"
    );
    std::fs::create_dir(dir.path().join("moon.rs")).unwrap();
    let error = generate_moon_module(&sample(), &config(1.0, dir.path())).unwrap_err();
    let prefix = format!("Failed to write {}: ", dir.path().join("moon.rs").display());
    assert!(error.starts_with(&prefix), "{error}");
}

#[test]
fn analysis() {
    let expected = format!(
        "\nELP/MPP02 Analysis (threshold: 1e0):\n{}\n\
         \nMain Problem Series:\n\
         \x20 LONGITUDE: 1 -> 1 terms (100.0%), amp range: 1.00e2 to 1.00e2\n\
         \x20 LATITUDE: 0 -> 0 terms (0.0%), amp range: 0.00e0 to 0.00e0\n\
         \x20 DISTANCE: 1 -> 0 terms (0.0%), amp range: 5.00e-1 to 5.00e-1\n\
         \nPerturbation Series:\n\
         \x20 LONGITUDE:\n\
         \x20   T^0: 1 -> 1 terms (100.0%), max amp: 2.00e0\n\
         \x20   T^1: 1 -> 0 terms (0.0%), max amp: 5.00e-1\n\
         \x20 LATITUDE:\n\
         \x20 DISTANCE:\n\
         \x20   T^2: 1 -> 1 terms (100.0%), max amp: 2.00e0\n\
         \nTotals:\n  Main problem: 2 terms\n  Perturbations: 3 terms\n  Total: 5 terms\n",
        "-".repeat(70)
    );
    assert_eq!(format_analysis(&sample(), 1.0), expected);
}
