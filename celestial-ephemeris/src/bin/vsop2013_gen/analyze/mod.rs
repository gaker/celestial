#[cfg(test)]
mod tests;

use crate::parser::{parse_file, planet_name, Variable, Vsop2013Block, Vsop2013File};
use crate::report::percent;
use std::collections::BTreeMap;
use std::path::Path;

pub(crate) struct VariableStats {
    pub(crate) variable: Variable,
    pub(crate) total_terms: usize,
    pub(crate) terms_above_threshold: usize,
    pub(crate) max_amplitude: f64,
    pub(crate) min_amplitude: f64,
    pub(crate) time_power_distribution: BTreeMap<u8, usize>,
}

pub(crate) struct PlanetAnalysis {
    pub(crate) planet: u8,
    pub(crate) planet_name: String,
    pub(crate) threshold: f64,
    pub(crate) total_terms: usize,
    pub(crate) variable_stats: Vec<VariableStats>,
}

impl PlanetAnalysis {
    pub(crate) fn terms_above_threshold(&self) -> usize {
        self.variable_stats
            .iter()
            .map(|v| v.terms_above_threshold)
            .sum()
    }
}

pub(crate) fn analyze_file(path: &Path, threshold: f64) -> Result<PlanetAnalysis, String> {
    let vsop = parse_file(path).map_err(|e| format!("Parse error: {}", e))?;
    Ok(analyze_vsop(&vsop, threshold))
}

pub(crate) fn analyze_vsop(vsop: &Vsop2013File, threshold: f64) -> PlanetAnalysis {
    PlanetAnalysis {
        planet: vsop.planet,
        planet_name: planet_name(vsop.planet).to_string(),
        threshold,
        total_terms: vsop.total_terms(),
        variable_stats: Variable::ALL
            .into_iter()
            .filter_map(|variable| variable_stats(vsop, variable, threshold))
            .collect(),
    }
}

fn variable_stats(
    vsop: &Vsop2013File,
    variable: Variable,
    threshold: f64,
) -> Option<VariableStats> {
    let blocks: Vec<&Vsop2013Block> = vsop.blocks_for_variable(variable).collect();
    if blocks.is_empty() {
        return None;
    }
    let terms = blocks.iter().flat_map(|b| &b.terms);
    let amplitudes: Vec<f64> = terms.map(|t| t.amplitude()).collect();
    Some(VariableStats {
        variable,
        total_terms: amplitudes.len(),
        terms_above_threshold: amplitudes.iter().filter(|&&a| a > threshold).count(),
        max_amplitude: amplitudes.iter().copied().reduce(f64::max).unwrap_or(0.0),
        min_amplitude: amplitudes.iter().copied().reduce(f64::min).unwrap_or(0.0),
        time_power_distribution: time_power_distribution(&blocks),
    })
}

fn time_power_distribution(blocks: &[&Vsop2013Block]) -> BTreeMap<u8, usize> {
    let mut distribution = BTreeMap::new();
    for block in blocks {
        *distribution.entry(block.header.time_power).or_insert(0) += block.terms.len();
    }
    distribution
}

pub(crate) fn format_analysis(analysis: &PlanetAnalysis) -> String {
    let rule = "=".repeat(70);
    let mut out = format!(
        "\n{rule}\nPlanet {}: {} ({} total terms)\n{rule}\n",
        analysis.planet, analysis.planet_name, analysis.total_terms
    );
    for stats in &analysis.variable_stats {
        out.push_str(&format_variable(stats, analysis.threshold));
    }
    let above = analysis.terms_above_threshold();
    out.push_str(&format!(
        "\n  Summary: {} terms above threshold ({:.1}%)\n",
        above,
        percent(above, analysis.total_terms)
    ));
    out
}

fn format_variable(stats: &VariableStats, threshold: f64) -> String {
    let above = stats.terms_above_threshold;
    let mut out = format!(
        "\n  Variable: {} ({})\n    Total terms: {}\n    \
         Terms above threshold ({:e}): {} ({:.1}%)\n    \
         Amplitude range: {:.3e} to {:.3e}\n    Time power distribution:\n",
        stats.variable,
        stats.variable.name(),
        stats.total_terms,
        threshold,
        above,
        percent(above, stats.total_terms),
        stats.min_amplitude,
        stats.max_amplitude
    );
    for (power, count) in &stats.time_power_distribution {
        out.push_str(&format!("      T^{}: {} terms\n", power, count));
    }
    out
}

pub(crate) fn format_summary_table(analyses: &[PlanetAnalysis]) -> String {
    let (equals, dashes) = ("=".repeat(80), "-".repeat(80));
    let mut out = format!(
        "\n{equals}\nVSOP2013 Summary\n{equals}\n{:<25} {:>9} {:>12} {:>12} {:>12} {:>6}\n{dashes}\n",
        "Planet", "Threshold", "Total Terms", "Above Thresh", "Reduction", "%"
    );
    for analysis in analyses {
        let threshold = format!("{:e}", analysis.threshold);
        let above = analysis.terms_above_threshold();
        out += &summary_row(
            &analysis.planet_name,
            &threshold,
            analysis.total_terms,
            above,
        );
    }
    let total = analyses.iter().map(|a| a.total_terms).sum();
    let above = analyses.iter().map(|a| a.terms_above_threshold()).sum();
    out + &dashes + "\n" + &summary_row("TOTAL", "", total, above)
}

fn summary_row(name: &str, threshold: &str, total: usize, above: usize) -> String {
    format!(
        "{:<25} {:>9} {:>12} {:>12} {:>12} {:>5.1}%\n",
        name,
        threshold,
        total,
        above,
        total - above,
        percent(above, total)
    )
}
