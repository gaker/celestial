use celestial_core::constants::TWOPI;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::download::ElpFilePaths;
use records::{
    is_pert_header, parse_main_header, parse_main_term, parse_pert_header, parse_pert_term,
};

mod records;
#[cfg(test)]
mod tests;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Coordinate {
    Longitude,
    Latitude,
    Distance,
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
pub(crate) struct MainTerm {
    pub(crate) delaunay: [i32; 4],
    pub(crate) coeffs: [f64; 7],
}

impl MainTerm {
    pub(crate) fn amplitude(&self) -> f64 {
        self.coeffs[0].abs()
    }
}

#[derive(Debug, Clone)]
pub(crate) struct PertTerm {
    pub(crate) sin_coeff: f64,
    pub(crate) cos_coeff: f64,
    pub(crate) multipliers: [i32; 16],
}

impl PertTerm {
    pub(crate) fn amplitude(&self) -> f64 {
        libm::sqrt(self.cos_coeff * self.cos_coeff + self.sin_coeff * self.sin_coeff)
    }

    // From 0 to 2π as the authors' code stores it, so the phase sums match theirs.
    pub(crate) fn phase(&self) -> f64 {
        let phase = libm::atan2(self.cos_coeff, self.sin_coeff);
        if phase < 0.0 {
            phase + TWOPI
        } else {
            phase
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct MainSeries {
    pub(crate) coordinate: Coordinate,
    pub(crate) terms: Vec<MainTerm>,
}

#[derive(Debug, Clone)]
pub(crate) struct PertBlock {
    pub(crate) time_power: u8,
    pub(crate) terms: Vec<PertTerm>,
}

#[derive(Debug, Clone)]
pub(crate) struct PertSeries {
    pub(crate) coordinate: Coordinate,
    pub(crate) blocks: Vec<PertBlock>,
}

#[derive(Debug, Clone)]
pub(crate) struct ElpData {
    pub(crate) main: [MainSeries; 3],
    pub(crate) pert: [PertSeries; 3],
}

impl ElpData {
    pub(crate) fn total_main_terms(&self) -> usize {
        self.main.iter().map(|s| s.terms.len()).sum()
    }

    pub(crate) fn total_pert_terms(&self) -> usize {
        self.pert
            .iter()
            .flat_map(|s| s.blocks.iter())
            .map(|b| b.terms.len())
            .sum()
    }

    pub(crate) fn total_terms(&self) -> usize {
        self.total_main_terms() + self.total_pert_terms()
    }
}

#[derive(Debug)]
pub(crate) enum ParseError {
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

fn parse_main_file(path: &Path, coordinate: Coordinate) -> Result<MainSeries, ParseError> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    let header = lines
        .next()
        .ok_or_else(|| ParseError::InvalidFormat("Empty file".to_string()))??;
    let term_count = parse_main_header(&header)?;
    let terms = lines
        .filter(|line| !line.as_ref().is_ok_and(|l| l.trim().is_empty()))
        .map(|line| parse_main_term(&line?))
        .collect::<Result<Vec<_>, _>>()?;
    if terms.len() != term_count {
        return Err(ParseError::InvalidFormat(format!(
            "Expected {} terms, found {}",
            term_count,
            terms.len()
        )));
    }
    Ok(MainSeries { coordinate, terms })
}

fn parse_pert_file(path: &Path, coordinate: Coordinate) -> Result<PertSeries, ParseError> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    let mut blocks = Vec::new();
    while let Some(line) = lines.next() {
        let line = line?;
        if is_pert_header(&line) {
            let (term_count, time_power) = parse_pert_header(&line)?;
            let terms = parse_pert_block(&mut lines, term_count, time_power)?;
            blocks.push(PertBlock { time_power, terms });
        } else if !line.trim().is_empty() {
            return Err(ParseError::InvalidFormat(format!(
                "Line outside any block: '{}'",
                line
            )));
        }
    }
    Ok(PertSeries { coordinate, blocks })
}

fn parse_pert_block(
    lines: &mut impl Iterator<Item = std::io::Result<String>>,
    term_count: usize,
    time_power: u8,
) -> Result<Vec<PertTerm>, ParseError> {
    let mut terms = Vec::with_capacity(term_count);
    while terms.len() < term_count {
        let Some(line) = lines.next() else {
            return Err(ParseError::InvalidFormat(format!(
                "Block T^{} expected {} terms, found {}",
                time_power,
                term_count,
                terms.len()
            )));
        };
        terms.push(parse_pert_term(&line?)?);
    }
    Ok(terms)
}

fn parse_main(path: &Path, file: u8, coordinate: Coordinate) -> Result<MainSeries, ParseError> {
    println!("  Parsing ELP_MAIN.S{} ({})...", file, coordinate);
    let series = parse_main_file(path, coordinate)?;
    println!("    {} terms", series.terms.len());
    Ok(series)
}

fn parse_pert(path: &Path, file: u8, coordinate: Coordinate) -> Result<PertSeries, ParseError> {
    println!(
        "  Parsing ELP_PERT.S{} ({} perturbations)...",
        file, coordinate
    );
    let series = parse_pert_file(path, coordinate)?;
    let terms: usize = series.blocks.iter().map(|b| b.terms.len()).sum();
    println!("    {} blocks, {} terms", series.blocks.len(), terms);
    Ok(series)
}

pub(crate) fn parse_files(paths: &ElpFilePaths) -> Result<ElpData, ParseError> {
    use Coordinate::{Distance, Latitude, Longitude};
    Ok(ElpData {
        main: [
            parse_main(&paths.main_longitude, 1, Longitude)?,
            parse_main(&paths.main_latitude, 2, Latitude)?,
            parse_main(&paths.main_distance, 3, Distance)?,
        ],
        pert: [
            parse_pert(&paths.pert_longitude, 1, Longitude)?,
            parse_pert(&paths.pert_latitude, 2, Latitude)?,
            parse_pert(&paths.pert_distance, 3, Distance)?,
        ],
    })
}
