use std::collections::BTreeMap;

use super::planet_module_name;
use crate::parser::{planet_name, Variable, Vsop2013File, Vsop2013Term};
use crate::report::percent;

// No term of the full theory uses more than six of the 17 arguments, so the
// tables store the nonzero multipliers in six slots, zero-padded.
pub(super) const ARGUMENTS_PER_TERM: usize = 6;

pub(super) struct FilteredTerm {
    pub(super) mult: [i32; ARGUMENTS_PER_TERM],
    pub(super) index: [u8; ARGUMENTS_PER_TERM],
    pub(super) s_coeff: f64,
    pub(super) c_coeff: f64,
    pub(super) amplitude: f64,
}

impl TryFrom<&Vsop2013Term> for FilteredTerm {
    type Error = String;

    fn try_from(term: &Vsop2013Term) -> Result<Self, String> {
        let mut mult = [0; ARGUMENTS_PER_TERM];
        let mut index = [0; ARGUMENTS_PER_TERM];
        let nonzero = (0u8..).zip(term.multipliers).filter(|&(_, m)| m != 0);
        for (slot, (i, m)) in nonzero.enumerate() {
            if slot == ARGUMENTS_PER_TERM {
                return Err(format!(
                    "a term uses more than {} arguments: {:?}",
                    ARGUMENTS_PER_TERM, term.multipliers
                ));
            }
            (mult[slot], index[slot]) = (m, i);
        }
        Ok(FilteredTerm {
            mult,
            index,
            s_coeff: term.s_coeff,
            c_coeff: term.c_coeff,
            amplitude: term.amplitude(),
        })
    }
}

pub(super) struct TimeBlockData {
    pub(super) power: u8,
    pub(super) terms: Vec<FilteredTerm>,
}

impl TimeBlockData {
    fn largest_first(power: u8, mut terms: Vec<FilteredTerm>) -> Self {
        terms.sort_by(|a, b| b.amplitude.total_cmp(&a.amplitude));
        TimeBlockData { power, terms }
    }
}

pub(super) struct VariableData {
    pub(super) variable: Variable,
    pub(super) blocks: Vec<TimeBlockData>,
    pub(super) total_terms: usize,
    pub(super) retained_terms: usize,
}

pub(super) fn filter_and_group_terms(
    vsop: &Vsop2013File,
    threshold: f64,
) -> Result<Vec<VariableData>, String> {
    let mut data = Vec::new();
    for variable in Variable::ALL {
        let variable_data = variable_data(vsop, variable, threshold)?;
        if !variable_data.blocks.is_empty() {
            data.push(variable_data);
        }
    }
    Ok(data)
}

fn variable_data(
    vsop: &Vsop2013File,
    variable: Variable,
    threshold: f64,
) -> Result<VariableData, String> {
    let mut by_power: BTreeMap<u8, Vec<FilteredTerm>> = BTreeMap::new();
    let mut total_terms = 0;
    for block in vsop.blocks_for_variable(variable) {
        total_terms += block.terms.len();
        let terms = by_power.entry(block.header.time_power).or_default();
        for term in block.terms.iter().filter(|t| t.amplitude() > threshold) {
            terms.push(FilteredTerm::try_from(term)?);
        }
    }
    by_power.retain(|_, terms| !terms.is_empty());
    Ok(VariableData {
        variable,
        retained_terms: by_power.values().map(Vec::len).sum(),
        blocks: by_power
            .into_iter()
            .map(|(power, terms)| TimeBlockData::largest_first(power, terms))
            .collect(),
        total_terms,
    })
}

pub(super) fn format_list<T: std::fmt::Display>(values: &[T]) -> String {
    let parts: Vec<String> = values.iter().map(T::to_string).collect();
    format!("[{}]", parts.join(","))
}

pub(super) fn format_float(f: f64) -> String {
    if f == 0.0 {
        "0.0".to_string()
    } else {
        format!("{:e}", f)
    }
}

pub(super) fn variable_const_name(var: Variable) -> &'static str {
    match var {
        Variable::A => "A",
        Variable::Lambda => "LAMBDA",
        Variable::K => "K",
        Variable::H => "H",
        Variable::Q => "Q",
        Variable::P => "P",
    }
}

pub(super) fn variable_description(var: Variable) -> &'static str {
    match var {
        Variable::A => "Semi-major axis (A)",
        Variable::Lambda => "Mean longitude (Lambda)",
        Variable::K => "e*cos(perihelion) (K)",
        Variable::H => "e*sin(perihelion) (H)",
        Variable::Q => "sin(i/2)*cos(node) (Q)",
        Variable::P => "sin(i/2)*sin(node) (P)",
    }
}

// Records the command that reproduces the file, and nothing that changes
// between runs, so regenerating unchanged data gives an identical file.
fn header(planet: u8, threshold: f64, retained: usize, total: usize) -> String {
    format!(
        "//! VSOP2013 coefficients for {}\n\
         //!\n\
         //! Generated from VSOP2013p{}.dat by\n\
         //! `vsop2013-gen generate --input <dir> --output <dir> --planet {} --threshold {:e}`\n\
         //! Terms retained: {} of {} ({:.1}%)\n\n",
        planet_name(planet),
        planet,
        planet,
        threshold,
        retained,
        total,
        percent(retained, total)
    )
}

pub(super) fn generate_planet_source(
    planet: u8,
    vsop: &Vsop2013File,
    threshold: f64,
) -> Result<(String, usize, usize), String> {
    let variable_data = filter_and_group_terms(vsop, threshold)?;
    let total: usize = variable_data.iter().map(|v| v.total_terms).sum();
    let retained: usize = variable_data.iter().map(|v| v.retained_terms).sum();
    let mut out = header(planet, threshold, retained, total);
    out.push_str("use super::{Term, TimeBlock};\n\n");
    for data in &variable_data {
        variable_source(data, &mut out);
    }
    Ok((out, retained, total))
}

fn variable_source(data: &VariableData, out: &mut String) {
    out.push_str(&format!(
        "/// {} coefficients\npub(crate) const {}: &[TimeBlock] = &[\n",
        variable_description(data.variable),
        variable_const_name(data.variable)
    ));
    for block in &data.blocks {
        time_block_source(block, out);
    }
    out.push_str("];\n\n");
}

fn time_block_source(block: &TimeBlockData, out: &mut String) {
    out.push_str(&format!(
        "    // T^{} terms\n    TimeBlock {{\n        power: {},\n        terms: &[\n",
        block.power, block.power
    ));
    for term in &block.terms {
        out.push_str(&format!(
            "            Term {{ s: {}, c: {}, mult: {}, index: {} }},\n",
            format_float(term.s_coeff),
            format_float(term.c_coeff),
            format_list(&term.mult),
            format_list(&term.index)
        ));
    }
    out.push_str("        ],\n    },\n");
}

const MOD_HEADER: &str = "\
//! VSOP2013 planetary coefficients
//!
//! Truncated coefficient tables for analytical planetary ephemeris.
//! These coefficients are used to compute heliocentric positions of planets.

";

fn mod_types() -> String {
    format!(
        "
/// A single Fourier term in the VSOP2013 series
#[derive(Debug, Clone, Copy)]
pub(crate) struct Term {{
    /// Sine coefficient
    pub(crate) s: f64,
    /// Cosine coefficient
    pub(crate) c: f64,
    // The nonzero multipliers, zero-padded, and which of the 17 fundamental
    // arguments each one multiplies.
    pub(crate) mult: [i32; {n}],
    pub(crate) index: [u8; {n}],
}}

/// Terms grouped by power of T (time)
#[derive(Debug, Clone, Copy)]
pub(crate) struct TimeBlock {{
    /// Power of T (0, 1, 2, ...)
    pub(crate) power: u8,
    /// Terms for this power, sorted by amplitude descending
    pub(crate) terms: &'static [Term],
}}

",
        n = ARGUMENTS_PER_TERM
    )
}

pub(super) fn generate_mod_rs(planets: &[u8]) -> String {
    let mut out = String::from(MOD_HEADER);
    for &p in planets {
        out.push_str(&format!("pub(crate) mod {};\n", planet_module_name(p)));
    }
    out + &mod_types()
}
