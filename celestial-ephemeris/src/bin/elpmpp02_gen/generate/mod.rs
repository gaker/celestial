mod source;
#[cfg(test)]
mod tests;

use crate::parser::{Coordinate, ElpData, MainSeries, MainTerm, PertBlock, PertSeries, PertTerm};
use crate::report::percent;
use source::moon_source;
use std::fs;

// The cut src/lunar_coefficients/moon.rs was generated with. It drops at
// most 19 m over 1950-2050, under the full series' own miss of DE432s with
// either fit (60 m and 21 m).
pub(crate) const DEFAULT_THRESHOLD: f64 = 1e-4;

pub(crate) struct GenerateConfig {
    pub(crate) threshold: f64,
    pub(crate) output_dir: std::path::PathBuf,
}

// Terms stay in file order, the order the authors' code sums them in.
fn filter_main_terms(series: &MainSeries, threshold: f64) -> Vec<&MainTerm> {
    series
        .terms
        .iter()
        .filter(|t| t.amplitude() >= threshold)
        .collect()
}

fn filter_pert_terms(block: &PertBlock, threshold: f64) -> Vec<&PertTerm> {
    block
        .terms
        .iter()
        .filter(|t| t.amplitude() >= threshold)
        .collect()
}

fn coord_name(coord: Coordinate) -> &'static str {
    match coord {
        Coordinate::Longitude => "LONGITUDE",
        Coordinate::Latitude => "LATITUDE",
        Coordinate::Distance => "DISTANCE",
    }
}

// Returns the line reporting what was written, for the caller to print.
pub(crate) fn generate_moon_module(
    elp: &ElpData,
    config: &GenerateConfig,
) -> Result<String, String> {
    let (source, retained, total) = moon_source(elp, config.threshold);
    fs::create_dir_all(&config.output_dir)
        .map_err(|e| format!("Failed to create output directory: {}", e))?;
    let output_path = config.output_dir.join("moon.rs");
    fs::write(&output_path, source)
        .map_err(|e| format!("Failed to write {}: {}", output_path.display(), e))?;
    Ok(format!(
        "Generated {} ({} of {} terms, {:.1}%, threshold {:e})",
        output_path.display(),
        retained,
        total,
        percent(retained, total),
        config.threshold
    ))
}

pub(crate) fn format_analysis(elp: &ElpData, threshold: f64) -> String {
    let mut out = format!(
        "\nELP/MPP02 Analysis (threshold: {:e}):\n{:-<70}\n\nMain Problem Series:\n",
        threshold, ""
    );
    for series in &elp.main {
        out += &main_summary(series, threshold);
    }
    out += "\nPerturbation Series:\n";
    for series in &elp.pert {
        out += &pert_summary(series, threshold);
    }
    out + &totals(elp)
}

fn totals(elp: &ElpData) -> String {
    format!(
        "\nTotals:\n  Main problem: {} terms\n  Perturbations: {} terms\n  Total: {} terms\n",
        elp.total_main_terms(),
        elp.total_pert_terms(),
        elp.total_terms()
    )
}

fn main_summary(series: &MainSeries, threshold: f64) -> String {
    let amplitudes = || series.terms.iter().map(MainTerm::amplitude);
    let (total, kept) = (
        series.terms.len(),
        filter_main_terms(series, threshold).len(),
    );
    format!(
        "  {}: {} -> {} terms ({:.1}%), amp range: {:.2e} to {:.2e}\n",
        coord_name(series.coordinate),
        total,
        kept,
        percent(kept, total),
        amplitudes().reduce(f64::min).unwrap_or(0.0),
        amplitudes().reduce(f64::max).unwrap_or(0.0)
    )
}

fn pert_summary(series: &PertSeries, threshold: f64) -> String {
    let mut out = format!("  {}:\n", coord_name(series.coordinate));
    for block in series.blocks.iter().filter(|b| !b.terms.is_empty()) {
        out += &block_summary(block, threshold);
    }
    out
}

fn block_summary(block: &PertBlock, threshold: f64) -> String {
    let (total, kept) = (block.terms.len(), filter_pert_terms(block, threshold).len());
    format!(
        "    T^{}: {} -> {} terms ({:.1}%), max amp: {:.2e}\n",
        block.time_power,
        total,
        kept,
        percent(kept, total),
        block
            .terms
            .iter()
            .map(PertTerm::amplitude)
            .fold(0.0, f64::max)
    )
}
