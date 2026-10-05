mod source;
#[cfg(test)]
mod tests;

use crate::parser::Vsop2013File;
use crate::report::percent;
use source::{generate_mod_rs, generate_planet_source};
use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

pub(crate) struct GenerateConfig {
    pub(crate) threshold: Option<f64>,
    pub(crate) output_dir: std::path::PathBuf,
}

// The shipped tables. Each cut moves the body's geocentric direction by at
// most 10 mas over 1950-2050, or by no more than the full series already
// misses DE432s where that is worse (Jupiter outward). 10 mas is about the
// aberration error from using heliocentric rather than barycentric velocity.
// Earth comes from the EMB series, and every geocentric position carries
// its error, so EMB keeps more terms than the rule alone would ask.
pub(crate) fn default_threshold(planet: u8) -> f64 {
    match planet {
        1 | 5 => 1e-9,
        2..=4 => 1e-10,
        6 => 3e-10,
        7 => 3e-8,
        8 => 1e-8,
        _ => 1e-7,
    }
}

const MODULE_NAMES: [&str; 9] = [
    "mercury", "venus", "emb", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto",
];

fn planet_module_name(planet: u8) -> &'static str {
    match planet {
        1..=9 => MODULE_NAMES[usize::from(planet) - 1],
        _ => "unknown",
    }
}

fn module_name_to_planet(name: &str) -> Option<u8> {
    let index = MODULE_NAMES.iter().position(|&m| m == name)?;
    Some(index as u8 + 1)
}

fn discover_existing_planets(output_dir: &Path) -> Result<BTreeSet<u8>, String> {
    let read_error = |e: std::io::Error| format!("Failed to read {}: {}", output_dir.display(), e);
    let mut planets = BTreeSet::new();
    for entry in fs::read_dir(output_dir).map_err(read_error)? {
        let path = entry.map_err(read_error)?.path();
        if path.extension().is_some_and(|ext| ext == "rs") {
            let stem = path.file_stem().and_then(|s| s.to_str());
            planets.extend(stem.and_then(module_name_to_planet));
        }
    }
    Ok(planets)
}

// Returns the line reporting what was written, for the caller to print.
fn generate_planet_module(
    planet: u8,
    vsop: &Vsop2013File,
    config: &GenerateConfig,
) -> Result<String, String> {
    let threshold = config.threshold.unwrap_or(default_threshold(planet));
    let (source, retained, total) = generate_planet_source(planet, vsop, threshold)?;
    let filename = format!("{}.rs", planet_module_name(planet));
    let filepath = config.output_dir.join(&filename);
    fs::write(&filepath, &source)
        .map_err(|e| format!("Failed to write {}: {}", filepath.display(), e))?;
    let percentage = percent(retained, total);
    Ok(format!(
        "  Generated {filename} ({retained} terms of {total} = {percentage:.1}%, threshold {threshold:e})"
    ))
}

fn write_mod_rs(output_dir: &Path, planets: BTreeSet<u8>) -> Result<(), String> {
    let planets: Vec<u8> = planets.into_iter().collect();
    fs::write(output_dir.join("mod.rs"), generate_mod_rs(&planets))
        .map_err(|e| format!("Failed to write mod.rs: {}", e))?;
    let names: Vec<&str> = planets.iter().map(|&p| planet_module_name(p)).collect();
    println!("  Generated mod.rs (planets: {:?})", names);
    Ok(())
}

pub(crate) fn generate_all(
    vsop_files: &[(u8, Vsop2013File)],
    config: &GenerateConfig,
) -> Result<(), String> {
    fs::create_dir_all(&config.output_dir)
        .map_err(|e| format!("Failed to create output dir: {}", e))?;
    let mut planets = discover_existing_planets(&config.output_dir)?;
    planets.extend(vsop_files.iter().map(|(p, _)| *p));
    write_mod_rs(&config.output_dir, planets)?;
    for (planet, vsop) in vsop_files {
        println!("{}", generate_planet_module(*planet, vsop, config)?);
    }
    Ok(())
}
