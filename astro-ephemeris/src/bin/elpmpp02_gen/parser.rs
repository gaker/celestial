use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::download::ElpFilePaths;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Coordinate {
    Longitude,
    Latitude,
    Distance,
}

impl Coordinate {
    #[allow(dead_code)]
    pub fn index(&self) -> usize {
        match self {
            Coordinate::Longitude => 0,
            Coordinate::Latitude => 1,
            Coordinate::Distance => 2,
        }
    }
}

impl fmt::Display for Coordinate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Coordinate::Longitude => write!(f, "Longitude"),
            Coordinate::Latitude => write!(f, "Latitude"),
            Coordinate::Distance => write!(f, "Distance"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct MainTerm {
    pub delaunay: [i32; 4],
    pub coeffs: [f64; 7],
}

impl MainTerm {
    pub fn amplitude(&self) -> f64 {
        self.coeffs[0].abs()
    }
}

#[derive(Debug, Clone)]
pub struct PertTerm {
    pub sin_coeff: f64,
    pub cos_coeff: f64,
    pub multipliers: [i32; 16],
}

impl PertTerm {
    pub fn amplitude(&self) -> f64 {
        (self.sin_coeff.powi(2) + self.cos_coeff.powi(2)).sqrt()
    }

    pub fn phase(&self) -> f64 {
        self.cos_coeff.atan2(self.sin_coeff)
    }
}

#[derive(Debug, Clone)]
pub struct MainSeries {
    pub coordinate: Coordinate,
    pub terms: Vec<MainTerm>,
}

#[derive(Debug, Clone)]
pub struct PertBlock {
    pub time_power: u8,
    pub terms: Vec<PertTerm>,
}

#[derive(Debug, Clone)]
pub struct PertSeries {
    pub coordinate: Coordinate,
    pub blocks: Vec<PertBlock>,
}

#[derive(Debug, Clone)]
pub struct ElpData {
    pub main: [MainSeries; 3],
    pub pert: [PertSeries; 3],
}

impl ElpData {
    pub fn total_main_terms(&self) -> usize {
        self.main.iter().map(|s| s.terms.len()).sum()
    }

    pub fn total_pert_terms(&self) -> usize {
        self.pert
            .iter()
            .flat_map(|s| s.blocks.iter())
            .map(|b| b.terms.len())
            .sum()
    }

    pub fn total_terms(&self) -> usize {
        self.total_main_terms() + self.total_pert_terms()
    }
}

#[derive(Debug)]
pub enum ParseError {
    IoError(std::io::Error),
    InvalidHeader(String),
    InvalidMainTerm(String),
    InvalidPertTerm(String),
    InvalidFormat(String),
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::IoError(e) => write!(f, "IO error: {}", e),
            ParseError::InvalidHeader(s) => write!(f, "Invalid header: {}", s),
            ParseError::InvalidMainTerm(s) => write!(f, "Invalid main term: {}", s),
            ParseError::InvalidPertTerm(s) => write!(f, "Invalid pert term: {}", s),
            ParseError::InvalidFormat(s) => write!(f, "Invalid format: {}", s),
        }
    }
}

impl std::error::Error for ParseError {}

impl From<std::io::Error> for ParseError {
    fn from(e: std::io::Error) -> Self {
        ParseError::IoError(e)
    }
}

fn parse_main_header(line: &str) -> Result<(String, usize), ParseError> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 3 {
        return Err(ParseError::InvalidHeader(format!(
            "Main header too short: '{}'",
            line
        )));
    }
    let term_count: usize = parts
        .last()
        .ok_or_else(|| ParseError::InvalidHeader("No term count".to_string()))?
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid term count in: '{}'", line)))?;

    Ok((line.to_string(), term_count))
}

fn parse_main_term(line: &str) -> Result<MainTerm, ParseError> {
    if line.len() < 80 {
        return Err(ParseError::InvalidMainTerm(format!(
            "Line too short ({}): '{}'",
            line.len(),
            line
        )));
    }

    let mut delaunay = [0i32; 4];
    for (i, d) in delaunay.iter_mut().enumerate() {
        let start = i * 3;
        let end = start + 3;
        let s = &line[start..end];
        *d = s.trim().parse().map_err(|_| {
            ParseError::InvalidMainTerm(format!("Invalid delaunay[{}]: '{}'", i, s))
        })?;
    }

    let mut coeffs = [0.0f64; 7];
    let coeff_start = 14;
    coeffs[0] = line[coeff_start..coeff_start + 13]
        .trim()
        .parse()
        .map_err(|e| {
            ParseError::InvalidMainTerm(format!(
                "Invalid coeff[0]: '{}' - {}",
                &line[coeff_start..coeff_start + 13],
                e
            ))
        })?;

    for (i, c) in coeffs.iter_mut().enumerate().skip(1) {
        let start = coeff_start + 13 + (i - 1) * 12;
        let end = start + 12;
        if end <= line.len() {
            let s = &line[start..end];
            *c = s.trim().parse().unwrap_or(0.0);
        }
    }

    Ok(MainTerm { delaunay, coeffs })
}

fn parse_main_file(path: &Path, coord: Coordinate) -> Result<MainSeries, ParseError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let first_line = lines
        .next()
        .ok_or_else(|| ParseError::InvalidFormat("Empty file".to_string()))??;
    let (_, term_count) = parse_main_header(&first_line)?;

    let mut terms = Vec::with_capacity(term_count);
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let term = parse_main_term(&line)?;
        terms.push(term);
    }

    if terms.len() != term_count {
        return Err(ParseError::InvalidFormat(format!(
            "Expected {} terms, found {}",
            term_count,
            terms.len()
        )));
    }

    Ok(MainSeries {
        coordinate: coord,
        terms,
    })
}

fn parse_pert_header(line: &str) -> Result<(usize, u8), ParseError> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 3 {
        return Err(ParseError::InvalidHeader(format!(
            "Pert header too short: '{}'",
            line
        )));
    }

    let term_count: usize = parts[parts.len() - 2]
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid term count in: '{}'", line)))?;

    let time_power: u8 = parts[parts.len() - 1]
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid time power in: '{}'", line)))?;

    Ok((term_count, time_power))
}

fn parse_fortran_double(s: &str) -> Result<f64, ParseError> {
    let s = s.trim();
    let s = s.replace('D', "E").replace('d', "e");
    s.parse()
        .map_err(|_| ParseError::InvalidPertTerm(format!("Invalid double: '{}'", s)))
}

fn parse_pert_term(line: &str) -> Result<PertTerm, ParseError> {
    if line.len() < 90 {
        return Err(ParseError::InvalidPertTerm(format!(
            "Line too short ({}): '{}'",
            line.len(),
            line
        )));
    }

    let sin_coeff = parse_fortran_double(&line[5..25])?;
    let cos_coeff = parse_fortran_double(&line[25..45])?;

    let mut multipliers = [0i32; 16];
    for (i, m) in multipliers.iter_mut().enumerate() {
        let start = 45 + i * 3;
        let end = start + 3;
        if end <= line.len() {
            let s = &line[start..end];
            *m = s.trim().parse().unwrap_or(0);
        }
    }

    Ok(PertTerm {
        sin_coeff,
        cos_coeff,
        multipliers,
    })
}

fn parse_pert_file(path: &Path, coord: Coordinate) -> Result<PertSeries, ParseError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines().peekable();

    let mut blocks = Vec::new();

    while let Some(line_result) = lines.next() {
        let line = line_result?;

        if line.contains("PERTURBATIONS")
            || line.contains("LONGITUDE")
            || line.contains("LATITUDE")
            || line.contains("DISTANCE")
        {
            let (term_count, time_power) = parse_pert_header(&line)?;

            let mut terms = Vec::with_capacity(term_count);
            for _ in 0..term_count {
                if let Some(term_line_result) = lines.next() {
                    let term_line = term_line_result?;
                    let term = parse_pert_term(&term_line)?;
                    terms.push(term);
                }
            }

            blocks.push(PertBlock { time_power, terms });
        }
    }

    Ok(PertSeries {
        coordinate: coord,
        blocks,
    })
}

pub fn parse_files(paths: &ElpFilePaths) -> Result<ElpData, ParseError> {
    println!("  Parsing ELP_MAIN.S1 (Longitude)...");
    let main_lon = parse_main_file(&paths.main_longitude, Coordinate::Longitude)?;
    println!("    {} terms", main_lon.terms.len());

    println!("  Parsing ELP_MAIN.S2 (Latitude)...");
    let main_lat = parse_main_file(&paths.main_latitude, Coordinate::Latitude)?;
    println!("    {} terms", main_lat.terms.len());

    println!("  Parsing ELP_MAIN.S3 (Distance)...");
    let main_dist = parse_main_file(&paths.main_distance, Coordinate::Distance)?;
    println!("    {} terms", main_dist.terms.len());

    println!("  Parsing ELP_PERT.S1 (Longitude perturbations)...");
    let pert_lon = parse_pert_file(&paths.pert_longitude, Coordinate::Longitude)?;
    let pert_lon_count: usize = pert_lon.blocks.iter().map(|b| b.terms.len()).sum();
    println!(
        "    {} blocks, {} terms",
        pert_lon.blocks.len(),
        pert_lon_count
    );

    println!("  Parsing ELP_PERT.S2 (Latitude perturbations)...");
    let pert_lat = parse_pert_file(&paths.pert_latitude, Coordinate::Latitude)?;
    let pert_lat_count: usize = pert_lat.blocks.iter().map(|b| b.terms.len()).sum();
    println!(
        "    {} blocks, {} terms",
        pert_lat.blocks.len(),
        pert_lat_count
    );

    println!("  Parsing ELP_PERT.S3 (Distance perturbations)...");
    let pert_dist = parse_pert_file(&paths.pert_distance, Coordinate::Distance)?;
    let pert_dist_count: usize = pert_dist.blocks.iter().map(|b| b.terms.len()).sum();
    println!(
        "    {} blocks, {} terms",
        pert_dist.blocks.len(),
        pert_dist_count
    );

    Ok(ElpData {
        main: [main_lon, main_lat, main_dist],
        pert: [pert_lon, pert_lat, pert_dist],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_coordinate_index() {
        assert_eq!(Coordinate::Longitude.index(), 0);
        assert_eq!(Coordinate::Latitude.index(), 1);
        assert_eq!(Coordinate::Distance.index(), 2);
    }

    #[test]
    fn test_coordinate_display() {
        assert_eq!(format!("{}", Coordinate::Longitude), "Longitude");
        assert_eq!(format!("{}", Coordinate::Latitude), "Latitude");
        assert_eq!(format!("{}", Coordinate::Distance), "Distance");
    }

    #[test]
    fn test_main_term_amplitude() {
        let term = MainTerm {
            delaunay: [0, 0, 1, 0],
            coeffs: [22639.55, 0.0, 0.0, 412529.61, 0.0, 0.0, 0.0],
        };
        assert!((term.amplitude() - 22639.55).abs() < 0.01);
    }

    #[test]
    fn test_main_term_amplitude_negative() {
        let term = MainTerm {
            delaunay: [0, 0, 0, 0],
            coeffs: [-12345.67, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        };
        assert!((term.amplitude() - 12345.67).abs() < 0.01);
    }

    #[test]
    fn test_pert_term_amplitude() {
        let term = PertTerm {
            sin_coeff: -12.749,
            cos_coeff: 6.369,
            multipliers: [0; 16],
        };
        let amp = term.amplitude();
        assert!(amp > 14.0 && amp < 15.0);
    }

    #[test]
    fn test_pert_term_phase() {
        let term = PertTerm {
            sin_coeff: 1.0,
            cos_coeff: 0.0,
            multipliers: [0; 16],
        };
        assert!((term.phase() - 0.0).abs() < 1e-10);

        let term2 = PertTerm {
            sin_coeff: 0.0,
            cos_coeff: 1.0,
            multipliers: [0; 16],
        };
        assert!((term2.phase() - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_elp_data_total_main_terms() {
        let data = ElpData {
            main: [
                MainSeries {
                    coordinate: Coordinate::Longitude,
                    terms: vec![
                        MainTerm {
                            delaunay: [0; 4],
                            coeffs: [0.0; 7],
                        },
                        MainTerm {
                            delaunay: [0; 4],
                            coeffs: [0.0; 7],
                        },
                    ],
                },
                MainSeries {
                    coordinate: Coordinate::Latitude,
                    terms: vec![MainTerm {
                        delaunay: [0; 4],
                        coeffs: [0.0; 7],
                    }],
                },
                MainSeries {
                    coordinate: Coordinate::Distance,
                    terms: vec![
                        MainTerm {
                            delaunay: [0; 4],
                            coeffs: [0.0; 7],
                        },
                        MainTerm {
                            delaunay: [0; 4],
                            coeffs: [0.0; 7],
                        },
                        MainTerm {
                            delaunay: [0; 4],
                            coeffs: [0.0; 7],
                        },
                    ],
                },
            ],
            pert: [
                PertSeries {
                    coordinate: Coordinate::Longitude,
                    blocks: vec![],
                },
                PertSeries {
                    coordinate: Coordinate::Latitude,
                    blocks: vec![],
                },
                PertSeries {
                    coordinate: Coordinate::Distance,
                    blocks: vec![],
                },
            ],
        };
        assert_eq!(data.total_main_terms(), 6);
    }

    #[test]
    fn test_elp_data_total_pert_terms() {
        let data = ElpData {
            main: [
                MainSeries {
                    coordinate: Coordinate::Longitude,
                    terms: vec![],
                },
                MainSeries {
                    coordinate: Coordinate::Latitude,
                    terms: vec![],
                },
                MainSeries {
                    coordinate: Coordinate::Distance,
                    terms: vec![],
                },
            ],
            pert: [
                PertSeries {
                    coordinate: Coordinate::Longitude,
                    blocks: vec![PertBlock {
                        time_power: 0,
                        terms: vec![
                            PertTerm {
                                sin_coeff: 0.0,
                                cos_coeff: 0.0,
                                multipliers: [0; 16],
                            },
                            PertTerm {
                                sin_coeff: 0.0,
                                cos_coeff: 0.0,
                                multipliers: [0; 16],
                            },
                        ],
                    }],
                },
                PertSeries {
                    coordinate: Coordinate::Latitude,
                    blocks: vec![PertBlock {
                        time_power: 1,
                        terms: vec![PertTerm {
                            sin_coeff: 0.0,
                            cos_coeff: 0.0,
                            multipliers: [0; 16],
                        }],
                    }],
                },
                PertSeries {
                    coordinate: Coordinate::Distance,
                    blocks: vec![],
                },
            ],
        };
        assert_eq!(data.total_pert_terms(), 3);
    }

    #[test]
    fn test_elp_data_total_terms() {
        let data = ElpData {
            main: [
                MainSeries {
                    coordinate: Coordinate::Longitude,
                    terms: vec![MainTerm {
                        delaunay: [0; 4],
                        coeffs: [0.0; 7],
                    }],
                },
                MainSeries {
                    coordinate: Coordinate::Latitude,
                    terms: vec![],
                },
                MainSeries {
                    coordinate: Coordinate::Distance,
                    terms: vec![],
                },
            ],
            pert: [
                PertSeries {
                    coordinate: Coordinate::Longitude,
                    blocks: vec![PertBlock {
                        time_power: 0,
                        terms: vec![
                            PertTerm {
                                sin_coeff: 0.0,
                                cos_coeff: 0.0,
                                multipliers: [0; 16],
                            },
                            PertTerm {
                                sin_coeff: 0.0,
                                cos_coeff: 0.0,
                                multipliers: [0; 16],
                            },
                        ],
                    }],
                },
                PertSeries {
                    coordinate: Coordinate::Latitude,
                    blocks: vec![],
                },
                PertSeries {
                    coordinate: Coordinate::Distance,
                    blocks: vec![],
                },
            ],
        };
        assert_eq!(data.total_terms(), 3);
    }

    #[test]
    fn test_parse_error_display_io_error() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
        let err = ParseError::IoError(io_err);
        let msg = format!("{}", err);
        assert!(msg.contains("IO error"));
        assert!(msg.contains("file not found"));
    }

    #[test]
    fn test_parse_error_display_invalid_header() {
        let err = ParseError::InvalidHeader("bad header".to_string());
        assert_eq!(format!("{}", err), "Invalid header: bad header");
    }

    #[test]
    fn test_parse_error_display_invalid_main_term() {
        let err = ParseError::InvalidMainTerm("bad term".to_string());
        assert_eq!(format!("{}", err), "Invalid main term: bad term");
    }

    #[test]
    fn test_parse_error_display_invalid_pert_term() {
        let err = ParseError::InvalidPertTerm("bad pert".to_string());
        assert_eq!(format!("{}", err), "Invalid pert term: bad pert");
    }

    #[test]
    fn test_parse_error_display_invalid_format() {
        let err = ParseError::InvalidFormat("bad format".to_string());
        assert_eq!(format!("{}", err), "Invalid format: bad format");
    }

    #[test]
    fn test_parse_error_from_io_error() {
        let io_err = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "access denied");
        let err: ParseError = io_err.into();
        match err {
            ParseError::IoError(e) => assert_eq!(e.kind(), std::io::ErrorKind::PermissionDenied),
            _ => panic!("Expected IoError variant"),
        }
    }

    #[test]
    fn test_parse_fortran_double() {
        let val = parse_fortran_double("-0.1274921554086D+02").unwrap();
        assert!((val - (-12.74921554086)).abs() < 1e-10);

        let val2 = parse_fortran_double("0.6368794709728D+01").unwrap();
        assert!((val2 - 6.368794709728).abs() < 1e-10);
    }

    #[test]
    fn test_parse_fortran_double_lowercase() {
        let val = parse_fortran_double("1.5d+00").unwrap();
        assert!((val - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_parse_fortran_double_standard_notation() {
        let val = parse_fortran_double("1.234e+05").unwrap();
        assert!((val - 123400.0).abs() < 1e-10);
    }

    #[test]
    fn test_parse_fortran_double_invalid() {
        let result = parse_fortran_double("not_a_number");
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidPertTerm(msg)) => assert!(msg.contains("Invalid double")),
            _ => panic!("Expected InvalidPertTerm error"),
        }
    }

    #[test]
    fn test_parse_main_header_valid() {
        let line = "ELP_MAIN.S1  Longitude  1023";
        let result = parse_main_header(line).unwrap();
        assert_eq!(result.0, line);
        assert_eq!(result.1, 1023);
    }

    #[test]
    fn test_parse_main_header_too_short() {
        let result = parse_main_header("AB");
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidHeader(msg)) => assert!(msg.contains("too short")),
            _ => panic!("Expected InvalidHeader error"),
        }
    }

    #[test]
    fn test_parse_main_header_invalid_term_count() {
        let result = parse_main_header("ELP MAIN abc");
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidHeader(msg)) => assert!(msg.contains("Invalid term count")),
            _ => panic!("Expected InvalidHeader error"),
        }
    }

    #[test]
    fn test_parse_main_term_valid() {
        // Format: 4 delaunay args (3 chars each), skip some, then 7 coeffs
        // Minimum 80 chars
        let line = "  0  0  1  0      22639.55000   0.00000   0.00000 412529.61   0.00000   0.00000   0.00000";
        let result = parse_main_term(line).unwrap();
        assert_eq!(result.delaunay, [0, 0, 1, 0]);
        assert!((result.coeffs[0] - 22639.55).abs() < 0.01);
    }

    #[test]
    fn test_parse_main_term_too_short() {
        let result = parse_main_term("short");
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidMainTerm(msg)) => assert!(msg.contains("too short")),
            _ => panic!("Expected InvalidMainTerm error"),
        }
    }

    #[test]
    fn test_parse_main_term_invalid_delaunay() {
        let line = "abc  0  1  0      22639.55000   0.00000   0.00000 412529.61   0.00000   0.00000   0.00000";
        let result = parse_main_term(line);
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidMainTerm(msg)) => assert!(msg.contains("Invalid delaunay")),
            _ => panic!("Expected InvalidMainTerm error"),
        }
    }

    #[test]
    fn test_parse_main_term_invalid_coeff() {
        // Make a line that's long enough but has invalid first coefficient
        let line = "  0  0  1  0      not_a_number   0.00000   0.00000 412529.61   0.00000   0.00000   0.00000";
        let result = parse_main_term(line);
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidMainTerm(msg)) => assert!(msg.contains("Invalid coeff[0]")),
            _ => panic!("Expected InvalidMainTerm error"),
        }
    }

    #[test]
    fn test_parse_pert_header_valid() {
        let line = "PERTURBATIONS LONGITUDE 234 0";
        let result = parse_pert_header(line).unwrap();
        assert_eq!(result.0, 234);
        assert_eq!(result.1, 0);
    }

    #[test]
    fn test_parse_pert_header_with_time_power() {
        let line = "PERTURBATIONS LONGITUDE 567 2";
        let result = parse_pert_header(line).unwrap();
        assert_eq!(result.0, 567);
        assert_eq!(result.1, 2);
    }

    #[test]
    fn test_parse_pert_header_too_short() {
        let result = parse_pert_header("AB");
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidHeader(msg)) => assert!(msg.contains("too short")),
            _ => panic!("Expected InvalidHeader error"),
        }
    }

    #[test]
    fn test_parse_pert_header_invalid_term_count() {
        let result = parse_pert_header("PERT LON abc 0");
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidHeader(msg)) => assert!(msg.contains("Invalid term count")),
            _ => panic!("Expected InvalidHeader error"),
        }
    }

    #[test]
    fn test_parse_pert_header_invalid_time_power() {
        let result = parse_pert_header("PERT LON 100 xyz");
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidHeader(msg)) => assert!(msg.contains("Invalid time power")),
            _ => panic!("Expected InvalidHeader error"),
        }
    }

    #[test]
    fn test_parse_pert_term_valid() {
        // Format: 5 skip, 20 sin, 20 cos, 16 multipliers (3 each) = minimum ~90 chars
        let line = "    1-0.1274921554086D+02 0.6368794709728D+01  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0";
        let result = parse_pert_term(line).unwrap();
        assert!((result.sin_coeff - (-12.74921554086)).abs() < 1e-10);
        assert!((result.cos_coeff - 6.368794709728).abs() < 1e-10);
    }

    #[test]
    fn test_parse_pert_term_too_short() {
        let result = parse_pert_term("short");
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidPertTerm(msg)) => assert!(msg.contains("too short")),
            _ => panic!("Expected InvalidPertTerm error"),
        }
    }

    #[test]
    fn test_parse_pert_term_with_multipliers() {
        let line = "    1 0.1000000000000D+01 0.2000000000000D+01  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16";
        let result = parse_pert_term(line).unwrap();
        assert_eq!(result.multipliers[0], 1);
        assert_eq!(result.multipliers[1], 2);
        assert_eq!(result.multipliers[15], 16);
    }

    fn create_mock_main_file(dir: &Path, name: &str, term_count: usize) -> std::path::PathBuf {
        let path = dir.join(name);
        let mut file = File::create(&path).unwrap();

        writeln!(file, "ELP_MAIN.S1 Longitude {}", term_count).unwrap();
        for i in 0..term_count {
            let coeff = 1000.0 - (i as f64);
            writeln!(
                file,
                "  {}  0  1  0      {:12.5}   0.00000   0.00000 412529.61   0.00000   0.00000   0.00000",
                i % 10, coeff
            ).unwrap();
        }
        path
    }

    fn create_mock_pert_file(dir: &Path, name: &str, blocks: &[(usize, u8)]) -> std::path::PathBuf {
        let path = dir.join(name);
        let mut file = File::create(&path).unwrap();

        for (term_count, time_power) in blocks {
            writeln!(
                file,
                "PERTURBATIONS LONGITUDE {} {}",
                term_count, time_power
            )
            .unwrap();
            for _ in 0..*term_count {
                writeln!(
                    file,
                    "    1-0.1274921554086D+02 0.6368794709728D+01  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0"
                ).unwrap();
            }
        }
        path
    }

    #[test]
    fn test_parse_main_file_valid() {
        let temp_dir = TempDir::new().unwrap();
        let path = create_mock_main_file(temp_dir.path(), "test_main.dat", 3);

        let result = parse_main_file(&path, Coordinate::Longitude).unwrap();
        assert_eq!(result.coordinate, Coordinate::Longitude);
        assert_eq!(result.terms.len(), 3);
    }

    #[test]
    fn test_parse_main_file_empty() {
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("empty.dat");
        File::create(&path).unwrap();

        let result = parse_main_file(&path, Coordinate::Longitude);
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidFormat(msg)) => assert!(msg.contains("Empty file")),
            _ => panic!("Expected InvalidFormat error"),
        }
    }

    #[test]
    fn test_parse_main_file_term_count_mismatch() {
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("mismatch.dat");
        let mut file = File::create(&path).unwrap();
        writeln!(file, "ELP_MAIN.S1 Longitude 5").unwrap();
        writeln!(
            file,
            "  0  0  1  0      22639.55000   0.00000   0.00000 412529.61   0.00000   0.00000   0.00000"
        ).unwrap();
        // Only wrote 1 term but header says 5

        let result = parse_main_file(&path, Coordinate::Longitude);
        assert!(result.is_err());
        match result {
            Err(ParseError::InvalidFormat(msg)) => {
                assert!(msg.contains("Expected 5 terms"));
                assert!(msg.contains("found 1"));
            }
            _ => panic!("Expected InvalidFormat error"),
        }
    }

    #[test]
    fn test_parse_main_file_skips_empty_lines() {
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("with_empty.dat");
        let mut file = File::create(&path).unwrap();
        writeln!(file, "ELP_MAIN.S1 Longitude 2").unwrap();
        writeln!(
            file,
            "  0  0  1  0      22639.55000   0.00000   0.00000 412529.61   0.00000   0.00000   0.00000"
        ).unwrap();
        writeln!(file).unwrap(); // empty line
        writeln!(file, "   ").unwrap(); // whitespace-only line
        writeln!(
            file,
            "  1  0  1  0      12345.00000   0.00000   0.00000 412529.61   0.00000   0.00000   0.00000"
        ).unwrap();

        let result = parse_main_file(&path, Coordinate::Longitude).unwrap();
        assert_eq!(result.terms.len(), 2);
    }

    #[test]
    fn test_parse_main_file_not_found() {
        let result = parse_main_file(Path::new("/nonexistent/file.dat"), Coordinate::Longitude);
        assert!(result.is_err());
        match result {
            Err(ParseError::IoError(_)) => {}
            _ => panic!("Expected IoError"),
        }
    }

    #[test]
    fn test_parse_pert_file_valid() {
        let temp_dir = TempDir::new().unwrap();
        let path = create_mock_pert_file(temp_dir.path(), "test_pert.dat", &[(3, 0), (2, 1)]);

        let result = parse_pert_file(&path, Coordinate::Longitude).unwrap();
        assert_eq!(result.coordinate, Coordinate::Longitude);
        assert_eq!(result.blocks.len(), 2);
        assert_eq!(result.blocks[0].terms.len(), 3);
        assert_eq!(result.blocks[0].time_power, 0);
        assert_eq!(result.blocks[1].terms.len(), 2);
        assert_eq!(result.blocks[1].time_power, 1);
    }

    #[test]
    fn test_parse_pert_file_multiple_keywords() {
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("multi_keyword.dat");
        let mut file = File::create(&path).unwrap();

        writeln!(file, "LONGITUDE 2 0").unwrap();
        writeln!(
            file,
            "    1-0.1000000000000D+01 0.2000000000000D+01  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0"
        ).unwrap();
        writeln!(
            file,
            "    1-0.3000000000000D+01 0.4000000000000D+01  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0"
        ).unwrap();
        writeln!(file, "LATITUDE 1 1").unwrap();
        writeln!(
            file,
            "    1-0.5000000000000D+01 0.6000000000000D+01  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0"
        ).unwrap();

        let result = parse_pert_file(&path, Coordinate::Longitude).unwrap();
        assert_eq!(result.blocks.len(), 2);
    }

    #[test]
    fn test_parse_pert_file_not_found() {
        let result = parse_pert_file(Path::new("/nonexistent/file.dat"), Coordinate::Longitude);
        assert!(result.is_err());
        match result {
            Err(ParseError::IoError(_)) => {}
            _ => panic!("Expected IoError"),
        }
    }

    #[test]
    fn test_parse_pert_file_with_distance_keyword() {
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join("distance.dat");
        let mut file = File::create(&path).unwrap();

        writeln!(file, "DISTANCE 1 0").unwrap();
        writeln!(
            file,
            "    1-0.1000000000000D+01 0.2000000000000D+01  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0"
        ).unwrap();

        let result = parse_pert_file(&path, Coordinate::Distance).unwrap();
        assert_eq!(result.blocks.len(), 1);
        assert_eq!(result.coordinate, Coordinate::Distance);
    }

    fn create_full_mock_elp_files(dir: &Path) -> ElpFilePaths {
        ElpFilePaths {
            main_longitude: create_mock_main_file(dir, "ELP_MAIN.S1", 2),
            main_latitude: create_mock_main_file(dir, "ELP_MAIN.S2", 3),
            main_distance: create_mock_main_file(dir, "ELP_MAIN.S3", 1),
            pert_longitude: create_mock_pert_file(dir, "ELP_PERT.S1", &[(2, 0)]),
            pert_latitude: create_mock_pert_file(dir, "ELP_PERT.S2", &[(1, 0), (1, 1)]),
            pert_distance: create_mock_pert_file(dir, "ELP_PERT.S3", &[(3, 0)]),
        }
    }

    #[test]
    fn test_parse_files_valid() {
        let temp_dir = TempDir::new().unwrap();
        let paths = create_full_mock_elp_files(temp_dir.path());

        let result = parse_files(&paths).unwrap();

        assert_eq!(result.main[0].coordinate, Coordinate::Longitude);
        assert_eq!(result.main[0].terms.len(), 2);
        assert_eq!(result.main[1].coordinate, Coordinate::Latitude);
        assert_eq!(result.main[1].terms.len(), 3);
        assert_eq!(result.main[2].coordinate, Coordinate::Distance);
        assert_eq!(result.main[2].terms.len(), 1);

        assert_eq!(result.pert[0].coordinate, Coordinate::Longitude);
        assert_eq!(result.pert[0].blocks.len(), 1);
        assert_eq!(result.pert[1].coordinate, Coordinate::Latitude);
        assert_eq!(result.pert[1].blocks.len(), 2);
        assert_eq!(result.pert[2].coordinate, Coordinate::Distance);
        assert_eq!(result.pert[2].blocks.len(), 1);

        assert_eq!(result.total_main_terms(), 6);
        // 2 (lon) + 1+1 (lat) + 3 (dist) = 7
        assert_eq!(result.total_pert_terms(), 7);
        assert_eq!(result.total_terms(), 13);
    }

    #[test]
    fn test_parse_files_missing_main_file() {
        let temp_dir = TempDir::new().unwrap();
        let paths = ElpFilePaths {
            main_longitude: temp_dir.path().join("nonexistent.dat"),
            main_latitude: create_mock_main_file(temp_dir.path(), "ELP_MAIN.S2", 1),
            main_distance: create_mock_main_file(temp_dir.path(), "ELP_MAIN.S3", 1),
            pert_longitude: create_mock_pert_file(temp_dir.path(), "ELP_PERT.S1", &[(1, 0)]),
            pert_latitude: create_mock_pert_file(temp_dir.path(), "ELP_PERT.S2", &[(1, 0)]),
            pert_distance: create_mock_pert_file(temp_dir.path(), "ELP_PERT.S3", &[(1, 0)]),
        };

        let result = parse_files(&paths);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_files_missing_pert_file() {
        let temp_dir = TempDir::new().unwrap();
        let paths = ElpFilePaths {
            main_longitude: create_mock_main_file(temp_dir.path(), "ELP_MAIN.S1", 1),
            main_latitude: create_mock_main_file(temp_dir.path(), "ELP_MAIN.S2", 1),
            main_distance: create_mock_main_file(temp_dir.path(), "ELP_MAIN.S3", 1),
            pert_longitude: create_mock_pert_file(temp_dir.path(), "ELP_PERT.S1", &[(1, 0)]),
            pert_latitude: temp_dir.path().join("nonexistent.dat"),
            pert_distance: create_mock_pert_file(temp_dir.path(), "ELP_PERT.S3", &[(1, 0)]),
        };

        let result = parse_files(&paths);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_main_term_short_line_coeffs() {
        // Line is >= 80 chars but may be shorter for later coefficients
        let line =
            "  0  0  1  0      22639.55000                                                       ";
        let result = parse_main_term(line).unwrap();
        assert!((result.coeffs[0] - 22639.55).abs() < 0.01);
        // Later coefficients should default to 0.0
        assert!((result.coeffs[1] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_parse_pert_term_short_multipliers() {
        // Line has valid sin/cos but maybe not all 16 multipliers
        let line = "    1 0.1000000000000D+01 0.2000000000000D+01  1  2  3  4                                      ";
        let result = parse_pert_term(line).unwrap();
        assert_eq!(result.multipliers[0], 1);
        assert_eq!(result.multipliers[1], 2);
        assert_eq!(result.multipliers[2], 3);
        assert_eq!(result.multipliers[3], 4);
        // Remaining should be 0
        assert_eq!(result.multipliers[4], 0);
    }

    #[test]
    fn test_parse_error_is_std_error() {
        fn assert_std_error<T: std::error::Error>() {}
        assert_std_error::<ParseError>();
    }
}
