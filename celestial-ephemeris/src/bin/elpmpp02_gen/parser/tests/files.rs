use super::*;
use std::path::PathBuf;
use tempfile::TempDir;

fn write(dir: &TempDir, name: &str, lines: &[&str]) -> PathBuf {
    let path = dir.path().join(name);
    let content: String = lines.iter().map(|line| format!("{line}\n")).collect();
    std::fs::write(&path, content).unwrap();
    path
}

fn main_file(dir: &TempDir, name: &str, terms: usize) -> PathBuf {
    let header = format!("MAIN PROBLEM. LONGITUDE. {terms}");
    let mut lines = vec![header.as_str()];
    lines.extend([MAIN_LINE].repeat(terms));
    write(dir, name, &lines)
}

fn pert_file(dir: &TempDir, name: &str, blocks: &[usize]) -> PathBuf {
    let mut lines = Vec::new();
    for (power, &terms) in blocks.iter().enumerate() {
        lines.push(format!("PERTURBATIONS. LONGITUDE. {terms} {power}"));
        lines.extend(std::iter::repeat_n(PERT_LINE.to_string(), terms));
    }
    let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
    write(dir, name, &lines)
}

#[test]
fn main_files() {
    let dir = TempDir::new().unwrap();
    let second = with_field(MAIN_LINE, 0, "  1");
    let header = "MAIN PROBLEM. LATITUDE. 2";
    let path = write(&dir, "main", &[header, MAIN_LINE, "", "   ", &second]);
    let series = parse_main_file(&path, Latitude).unwrap();
    assert_eq!(series.coordinate, Latitude);
    let delaunay: Vec<[i32; 4]> = series.terms.iter().map(|t| t.delaunay).collect();
    assert_eq!(delaunay, [[0, 2, 0, 0], [1, 2, 0, 0]]);
}

#[test]
fn main_file_errors() {
    let dir = TempDir::new().unwrap();
    let cases: [(&[&str], &str); 4] = [
        (&[], "Invalid format: Empty file"),
        (
            &["MAIN 1"],
            "Invalid header: Main header too short: 'MAIN 1'",
        ),
        (
            &["MAIN PROBLEM. 5", MAIN_LINE],
            "Invalid format: Expected 5 terms, found 1",
        ),
        (
            &["MAIN PROBLEM. 1", MAIN_LINE, MAIN_LINE],
            "Invalid format: Expected 1 terms, found 2",
        ),
    ];
    for (lines, message) in cases {
        let path = write(&dir, "main", lines);
        let error = parse_main_file(&path, Longitude).unwrap_err();
        assert_eq!(error.to_string(), message);
    }
    let missing = parse_main_file(&dir.path().join("missing"), Longitude);
    assert!(matches!(missing, Err(ParseError::IoError(_))));
}

#[test]
fn pert_files() {
    let dir = TempDir::new().unwrap();
    let lines = [
        "PERTURBATIONS. LONGITUDE.  2  0",
        PERT_LINE,
        PERT_LINE,
        "",
        "LATITUDE 0 1",
        "DISTANCE 1 3",
        PERT_LINE,
    ];
    let series = parse_pert_file(&write(&dir, "pert", &lines), Distance).unwrap();
    assert_eq!(series.coordinate, Distance);
    let blocks: Vec<(u8, usize)> = (series.blocks.iter())
        .map(|b| (b.time_power, b.terms.len()))
        .collect();
    assert_eq!(blocks, [(0, 2), (1, 0), (3, 1)]);
}

#[test]
fn pert_file_errors() {
    let dir = TempDir::new().unwrap();
    let header = "PERTURBATIONS. LONGITUDE.         2         0";
    let outside = format!("Invalid format: Line outside any block: '{PERT_LINE}'");
    let cases: [(&[&str], &str); 4] = [
        (
            &[header, PERT_LINE],
            "Invalid format: Block T^0 expected 2 terms, found 1",
        ),
        (&[header, PERT_LINE, PERT_LINE, PERT_LINE], &outside),
        (&[PERT_LINE], &outside),
        (
            &["PERTURBATIONS 1"],
            "Invalid header: Pert header too short: 'PERTURBATIONS 1'",
        ),
    ];
    for (lines, message) in cases {
        let path = write(&dir, "pert", lines);
        let error = parse_pert_file(&path, Longitude).unwrap_err();
        assert_eq!(error.to_string(), message);
    }
    let missing = parse_pert_file(&dir.path().join("missing"), Longitude);
    assert!(matches!(missing, Err(ParseError::IoError(_))));
}

fn elp_files(dir: &TempDir) -> ElpFilePaths {
    ElpFilePaths {
        main_longitude: main_file(dir, "ELP_MAIN.S1", 2),
        main_latitude: main_file(dir, "ELP_MAIN.S2", 3),
        main_distance: main_file(dir, "ELP_MAIN.S3", 1),
        pert_longitude: pert_file(dir, "ELP_PERT.S1", &[2]),
        pert_latitude: pert_file(dir, "ELP_PERT.S2", &[1, 1]),
        pert_distance: pert_file(dir, "ELP_PERT.S3", &[3]),
    }
}

#[test]
fn all_six_files() {
    let dir = TempDir::new().unwrap();
    let data = parse_files(&elp_files(&dir)).unwrap();
    let main: Vec<(Coordinate, usize)> = (data.main.iter())
        .map(|s| (s.coordinate, s.terms.len()))
        .collect();
    assert_eq!(main, [(Longitude, 2), (Latitude, 3), (Distance, 1)]);
    let pert: Vec<(Coordinate, usize)> = (data.pert.iter())
        .map(|s| (s.coordinate, s.blocks.len()))
        .collect();
    assert_eq!(pert, [(Longitude, 1), (Latitude, 2), (Distance, 1)]);
    assert_eq!((data.total_main_terms(), data.total_pert_terms()), (6, 7));
}

#[test]
fn any_missing_file_is_an_error() {
    let dir = TempDir::new().unwrap();
    let good = elp_files(&dir);
    let missing = dir.path().join("nonexistent");
    let broken = [
        ElpFilePaths {
            main_longitude: missing.clone(),
            ..good.clone()
        },
        ElpFilePaths {
            pert_latitude: missing,
            ..good
        },
    ];
    for paths in broken {
        assert!(matches!(parse_files(&paths), Err(ParseError::IoError(_))));
    }
}

#[test]
#[ignore = "requires local ELP/MPP02 data files"]
fn shipped_files() {
    let dir = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../references/ephemeris/elpmpp02"
    ));
    let data = parse_files(&crate::download::find_elp_files(dir).unwrap()).unwrap();
    let main: Vec<usize> = data.main.iter().map(|s| s.terms.len()).collect();
    assert_eq!(main, [1023, 918, 704]);
    let pert: Vec<Vec<usize>> = (data.pert.iter())
        .map(|s| s.blocks.iter().map(|b| b.terms.len()).collect())
        .collect();
    assert_eq!(
        pert,
        [
            [11314, 1199, 219, 2],
            [6462, 516, 52, 0],
            [12115, 1165, 210, 2]
        ]
    );
}
