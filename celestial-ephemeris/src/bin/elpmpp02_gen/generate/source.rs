use super::{coord_name, filter_main_terms, filter_pert_terms};
use crate::parser::{ElpData, MainSeries, PertBlock, PertSeries};
use crate::report::percent;

pub(super) const LINT: &str = "\
// Coefficients such as 3.14, and phases of exactly pi or pi/2 from terms with
// only a sine or only a cosine part, are data rather than stand-ins for the
// constants.
#![allow(clippy::approx_constant)]

";

pub(super) const TYPES: &str = "\
/// Main problem term with Delaunay argument multipliers
#[derive(Debug, Clone, Copy)]
pub(crate) struct MainTerm {
    /// Multipliers for D, F, l, l' (Delaunay arguments)
    pub(crate) delaunay: [i8; 4],
    /// Amplitude A, then B1 through B6, its derivatives with respect to the
    /// fitted constants (B6 is unused)
    pub(crate) coeffs: [f64; 7],
}

/// Perturbation term with full argument multipliers
#[derive(Debug, Clone, Copy)]
pub(crate) struct PertTerm {
    /// Amplitude (sqrt(S^2 + C^2))
    pub(crate) amplitude: f64,
    /// Phase angle atan2(C, S), with 0 ≤ phase < 2π
    pub(crate) phase: f64,
    /// Multipliers: [D, F, l, l', Me, Ve, Te, Ma, Ju, Sa, Ur, Ne, zeta, ?, ?, ?]
    pub(crate) multipliers: [i8; 16],
}

/// Perturbation block for a specific time power
#[derive(Debug, Clone, Copy)]
pub(crate) struct PertBlock {
    /// Power of T (0, 1, 2, 3)
    pub(crate) power: u8,
    /// Terms for this power
    pub(crate) terms: &'static [PertTerm],
}

";

// Shortest text that parses back to the same f64, so the table holds exactly
// the values the series files define.
pub(super) fn format_float(f: f64) -> String {
    if f == 0.0 {
        "0.0".to_string()
    } else {
        format!("{:e}", f)
    }
}

pub(super) fn format_floats(values: &[f64]) -> String {
    values
        .iter()
        .map(|&v| format_float(v))
        .collect::<Vec<_>>()
        .join(", ")
}

pub(super) fn format_ints(values: &[i32]) -> String {
    values
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

// Records the command that reproduces the file, and nothing that changes
// between runs, so regenerating unchanged data gives an identical file.
pub(super) fn header(threshold: f64, retained: usize, total: usize) -> String {
    format!(
        "//! ELP/MPP02 coefficients for the Moon\n\
         //!\n\
         //! Generated from ELP_MAIN.S1-S3 and ELP_PERT.S1-S3 by\n\
         //! `elpmpp02-gen generate --input <dir> --output <dir> --threshold {:e}`\n\
         //! Terms retained: {} of {} ({:.1}%)\n\
         //!\n\
         //! Reference: Chapront & Francou (2003)\n\
         //! \"The lunar theory ELP revisited. Introduction of new planetary perturbations\"\n\
         //! Astronomy & Astrophysics, 404, 735-742\n\n",
        threshold,
        retained,
        total,
        percent(retained, total)
    )
}

fn main_source(series: &MainSeries, threshold: f64, out: &mut String) -> (usize, usize) {
    let filtered = filter_main_terms(series, threshold);
    let name = coord_name(series.coordinate);
    out.push_str(&format!(
        "/// Main problem terms for {} ({} of {} terms)\npub(crate) const MAIN_{}: &[MainTerm] = &[\n",
        name,
        filtered.len(),
        series.terms.len(),
        name
    ));
    for term in &filtered {
        out.push_str(&format!(
            "    MainTerm {{ delaunay: [{}], coeffs: [{}] }},\n",
            format_ints(&term.delaunay),
            format_floats(&term.coeffs)
        ));
    }
    out.push_str("];\n\n");
    (filtered.len(), series.terms.len())
}

fn pert_block_source(name: &str, block: &PertBlock, threshold: f64, out: &mut String) -> usize {
    let filtered = filter_pert_terms(block, threshold);
    let (power, kept, total) = (block.time_power, filtered.len(), block.terms.len());
    out.push_str(&format!(
        "/// Perturbation terms for {name} T^{power} ({kept} of {total} terms)\n\
         const PERT_{name}_T{power}: &[PertTerm] = &[\n"
    ));
    for term in &filtered {
        out.push_str(&format!(
            "    PertTerm {{ amplitude: {}, phase: {}, multipliers: [{}] }},\n",
            format_float(term.amplitude()),
            format_float(term.phase()),
            format_ints(&term.multipliers)
        ));
    }
    out.push_str("];\n\n");
    filtered.len()
}

fn pert_source(series: &PertSeries, threshold: f64, out: &mut String) -> (usize, usize) {
    let name = coord_name(series.coordinate);
    let mut retained = 0;
    for block in &series.blocks {
        retained += pert_block_source(name, block, threshold, out);
    }
    out.push_str(&format!(
        "/// All perturbation blocks for {}\npub(crate) const PERT_{}: &[PertBlock] = &[\n",
        name, name
    ));
    for block in &series.blocks {
        out.push_str(&format!(
            "    PertBlock {{ power: {}, terms: PERT_{}_T{} }},\n",
            block.time_power, name, block.time_power
        ));
    }
    out.push_str("];\n\n");
    let total = series.blocks.iter().map(|b| b.terms.len()).sum();
    (retained, total)
}

pub(super) fn moon_source(elp: &ElpData, threshold: f64) -> (String, usize, usize) {
    let mut body = String::from(LINT) + TYPES;
    let (mut retained, mut total) = (0, 0);
    for series in &elp.main {
        let (r, t) = main_source(series, threshold, &mut body);
        (retained, total) = (retained + r, total + t);
    }
    for series in &elp.pert {
        let (r, t) = pert_source(series, threshold, &mut body);
        (retained, total) = (retained + r, total + t);
    }
    (header(threshold, retained, total) + &body, retained, total)
}
