use crate::parser::{planet_name, Variable, Vsop2013File, Vsop2013Term};
use chrono::Utc;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

pub struct GenerateConfig {
    pub threshold: f64,
    pub output_dir: std::path::PathBuf,
}

struct FilteredTerm {
    multipliers: [i16; 17],
    s_coeff: f64,
    c_coeff: f64,
    amplitude: f64,
}

impl From<&Vsop2013Term> for FilteredTerm {
    fn from(term: &Vsop2013Term) -> Self {
        let mut mults = [0i16; 17];
        for (i, &m) in term.multipliers.iter().enumerate() {
            mults[i] = m as i16;
        }
        FilteredTerm {
            multipliers: mults,
            s_coeff: term.s_coeff,
            c_coeff: term.c_coeff,
            amplitude: term.amplitude(),
        }
    }
}

struct TimeBlockData {
    power: u8,
    terms: Vec<FilteredTerm>,
}

struct VariableData {
    variable: Variable,
    blocks: Vec<TimeBlockData>,
    total_terms: usize,
    retained_terms: usize,
}

fn filter_and_group_terms(vsop: &Vsop2013File, threshold: f64) -> Vec<VariableData> {
    let mut result = Vec::new();

    for var in [
        Variable::A,
        Variable::Lambda,
        Variable::K,
        Variable::H,
        Variable::Q,
        Variable::P,
    ] {
        let mut blocks_by_power: BTreeMap<u8, Vec<FilteredTerm>> = BTreeMap::new();
        let mut total = 0usize;
        let mut retained = 0usize;

        for block in vsop.blocks_for_variable(var) {
            total += block.terms.len();

            let filtered: Vec<FilteredTerm> = block
                .terms
                .iter()
                .filter(|t| t.amplitude() > threshold)
                .map(FilteredTerm::from)
                .collect();

            retained += filtered.len();

            if !filtered.is_empty() {
                blocks_by_power
                    .entry(block.header.time_power)
                    .or_default()
                    .extend(filtered);
            }
        }

        let mut blocks: Vec<TimeBlockData> = blocks_by_power
            .into_iter()
            .map(|(power, mut terms)| {
                terms.sort_by(|a, b| {
                    b.amplitude
                        .partial_cmp(&a.amplitude)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                TimeBlockData { power, terms }
            })
            .collect();

        blocks.sort_by_key(|b| b.power);

        if !blocks.is_empty() {
            result.push(VariableData {
                variable: var,
                blocks,
                total_terms: total,
                retained_terms: retained,
            });
        }
    }

    result
}

fn format_multipliers(mults: &[i16; 17]) -> String {
    let parts: Vec<String> = mults.iter().map(|m| format!("{}", m)).collect();
    format!("[{}]", parts.join(","))
}

fn format_float(f: f64) -> String {
    if f == 0.0 {
        "0.0".to_string()
    } else {
        format!("{:.15e}", f)
    }
}

fn variable_const_name(var: Variable) -> &'static str {
    match var {
        Variable::A => "A",
        Variable::Lambda => "LAMBDA",
        Variable::K => "K",
        Variable::H => "H",
        Variable::Q => "Q",
        Variable::P => "P",
    }
}

fn variable_description(var: Variable) -> &'static str {
    match var {
        Variable::A => "Semi-major axis (A)",
        Variable::Lambda => "Mean longitude (Lambda)",
        Variable::K => "e*cos(perihelion) (K)",
        Variable::H => "e*sin(perihelion) (H)",
        Variable::Q => "sin(i/2)*cos(node) (Q)",
        Variable::P => "sin(i/2)*sin(node) (P)",
    }
}

fn generate_planet_source(
    planet: u8,
    vsop: &Vsop2013File,
    threshold: f64,
) -> (String, usize, usize) {
    let variable_data = filter_and_group_terms(vsop, threshold);

    let total_original: usize = variable_data.iter().map(|v| v.total_terms).sum();
    let total_retained: usize = variable_data.iter().map(|v| v.retained_terms).sum();

    let percentage = if total_original > 0 {
        (total_retained as f64 / total_original as f64) * 100.0
    } else {
        0.0
    };

    let date = Utc::now().format("%Y-%m-%d").to_string();

    let mut out = String::new();

    out.push_str("#![allow(clippy::excessive_precision)]\n");
    out.push('\n');
    out.push_str(&format!(
        "//! VSOP2013 coefficients for {}\n",
        planet_name(planet)
    ));
    out.push_str("//!\n");
    out.push_str(&format!("//! Generated from VSOP2013p{}.dat\n", planet));
    out.push_str(&format!("//! Threshold: {:.0e}\n", threshold));
    out.push_str(&format!(
        "//! Terms retained: {} of {} ({:.1}%)\n",
        total_retained, total_original, percentage
    ));
    out.push_str(&format!("//! Generated: {}\n", date));
    out.push('\n');
    out.push_str("use super::{Term, TimeBlock};\n");
    out.push('\n');

    for var_data in &variable_data {
        out.push_str(&format!(
            "/// {} coefficients\n",
            variable_description(var_data.variable)
        ));
        out.push_str(&format!(
            "pub const {}: &[TimeBlock] = &[\n",
            variable_const_name(var_data.variable)
        ));

        for block in &var_data.blocks {
            out.push_str(&format!("    // T^{} terms\n", block.power));
            out.push_str("    TimeBlock {\n");
            out.push_str(&format!("        power: {},\n", block.power));
            out.push_str("        terms: &[\n");

            for term in &block.terms {
                out.push_str(&format!(
                    "            Term {{ s: {}, c: {}, mult: {} }},\n",
                    format_float(term.s_coeff),
                    format_float(term.c_coeff),
                    format_multipliers(&term.multipliers)
                ));
            }

            out.push_str("        ],\n");
            out.push_str("    },\n");
        }

        out.push_str("];\n");
        out.push('\n');
    }

    (out, total_retained, total_original)
}

fn planet_module_name(planet: u8) -> &'static str {
    match planet {
        1 => "mercury",
        2 => "venus",
        3 => "emb",
        4 => "mars",
        5 => "jupiter",
        6 => "saturn",
        7 => "uranus",
        8 => "neptune",
        9 => "pluto",
        _ => "unknown",
    }
}

fn module_name_to_planet(name: &str) -> Option<u8> {
    match name {
        "mercury" => Some(1),
        "venus" => Some(2),
        "emb" => Some(3),
        "mars" => Some(4),
        "jupiter" => Some(5),
        "saturn" => Some(6),
        "uranus" => Some(7),
        "neptune" => Some(8),
        "pluto" => Some(9),
        _ => None,
    }
}

fn discover_existing_planets(output_dir: &Path) -> BTreeSet<u8> {
    let mut planets = BTreeSet::new();
    if let Ok(entries) = fs::read_dir(output_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.extension().is_some_and(|ext| ext == "rs") {
                if let Some(stem) = path.file_stem().and_then(|s| s.to_str()) {
                    if stem != "mod" {
                        if let Some(planet_num) = module_name_to_planet(stem) {
                            planets.insert(planet_num);
                        }
                    }
                }
            }
        }
    }
    planets
}

pub fn generate_mod_rs(planets: &[u8]) -> String {
    let mut out = String::new();

    out.push_str("//! VSOP2013 planetary coefficients\n");
    out.push_str("//!\n");
    out.push_str("//! Truncated coefficient tables for analytical planetary ephemeris.\n");
    out.push_str("//! These coefficients are used to compute heliocentric positions of planets.\n");
    out.push('\n');

    for &p in planets {
        out.push_str(&format!("pub mod {};\n", planet_module_name(p)));
    }

    out.push('\n');
    out.push_str("/// A single Fourier term in the VSOP2013 series\n");
    out.push_str("#[derive(Debug, Clone, Copy)]\n");
    out.push_str("pub struct Term {\n");
    out.push_str("    /// Sine coefficient\n");
    out.push_str("    pub s: f64,\n");
    out.push_str("    /// Cosine coefficient\n");
    out.push_str("    pub c: f64,\n");
    out.push_str("    /// Argument multipliers for 17 fundamental arguments\n");
    out.push_str("    pub mult: [i16; 17],\n");
    out.push_str("}\n");
    out.push('\n');
    out.push_str("/// Terms grouped by power of T (time)\n");
    out.push_str("#[derive(Debug, Clone, Copy)]\n");
    out.push_str("pub struct TimeBlock {\n");
    out.push_str("    /// Power of T (0, 1, 2, ...)\n");
    out.push_str("    pub power: u8,\n");
    out.push_str("    /// Terms for this power, sorted by amplitude descending\n");
    out.push_str("    pub terms: &'static [Term],\n");
    out.push_str("}\n");
    out.push('\n');

    out
}

pub fn generate_planet_module(
    planet: u8,
    vsop: &Vsop2013File,
    config: &GenerateConfig,
    output_dir: &Path,
) -> Result<(), String> {
    let (source, retained, total) = generate_planet_source(planet, vsop, config.threshold);

    let filename = format!("{}.rs", planet_module_name(planet));
    let filepath = output_dir.join(&filename);

    fs::write(&filepath, &source)
        .map_err(|e| format!("Failed to write {}: {}", filepath.display(), e))?;

    let percentage = if total > 0 {
        (retained as f64 / total as f64) * 100.0
    } else {
        0.0
    };

    println!(
        "  Generated {} ({} terms of {} = {:.1}%)",
        filename, retained, total, percentage
    );

    Ok(())
}

pub fn generate_all(
    vsop_files: &[(u8, Vsop2013File)],
    config: &GenerateConfig,
) -> Result<(), String> {
    fs::create_dir_all(&config.output_dir)
        .map_err(|e| format!("Failed to create output dir: {}", e))?;

    let mut all_planets: BTreeSet<u8> = discover_existing_planets(&config.output_dir);
    for (p, _) in vsop_files {
        all_planets.insert(*p);
    }

    let planets: Vec<u8> = all_planets.into_iter().collect();

    let mod_source = generate_mod_rs(&planets);
    let mod_path = config.output_dir.join("mod.rs");
    fs::write(&mod_path, &mod_source).map_err(|e| format!("Failed to write mod.rs: {}", e))?;
    println!(
        "  Generated mod.rs (planets: {:?})",
        planets
            .iter()
            .map(|p| planet_module_name(*p))
            .collect::<Vec<_>>()
    );

    for (planet, vsop) in vsop_files {
        generate_planet_module(*planet, vsop, config, &config.output_dir)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{Vsop2013Block, Vsop2013Header};
    use tempfile::TempDir;

    fn make_term(s: f64, c: f64) -> Vsop2013Term {
        Vsop2013Term {
            multipliers: [0; 17],
            s_coeff: s,
            c_coeff: c,
        }
    }

    fn make_term_with_mults(s: f64, c: f64, mults: [i32; 17]) -> Vsop2013Term {
        Vsop2013Term {
            multipliers: mults,
            s_coeff: s,
            c_coeff: c,
        }
    }

    fn make_block(var: Variable, time_power: u8, terms: Vec<Vsop2013Term>) -> Vsop2013Block {
        Vsop2013Block {
            header: Vsop2013Header {
                planet: 3,
                variable: var,
                time_power,
                term_count: terms.len() as u32,
            },
            terms,
        }
    }

    fn make_test_vsop() -> Vsop2013File {
        Vsop2013File {
            planet: 3,
            blocks: vec![
                make_block(
                    Variable::A,
                    0,
                    vec![
                        make_term(3.0, 4.0),   // amplitude = 5
                        make_term(0.6, 0.8),   // amplitude = 1
                        make_term(0.03, 0.04), // amplitude = 0.05
                    ],
                ),
                make_block(
                    Variable::A,
                    1,
                    vec![
                        make_term(6.0, 8.0), // amplitude = 10
                    ],
                ),
                make_block(
                    Variable::Lambda,
                    0,
                    vec![
                        make_term(30.0, 40.0), // amplitude = 50
                    ],
                ),
            ],
        }
    }

    #[test]
    fn test_format_float() {
        assert_eq!(format_float(0.0), "0.0");
        assert!(format_float(1.5e-10).contains("1.5"));
        assert!(format_float(-123.456789).contains("-1.23456789"));
    }

    #[test]
    fn test_format_float_negative_zero() {
        // -0.0 should be treated as 0.0
        assert_eq!(format_float(-0.0), "0.0");
    }

    #[test]
    fn test_format_multipliers() {
        let mults = [0i16, 1, -2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let result = format_multipliers(&mults);
        assert_eq!(result, "[0,1,-2,3,0,0,0,0,0,0,0,0,0,0,0,0,0]");
    }

    #[test]
    fn test_variable_const_name() {
        assert_eq!(variable_const_name(Variable::A), "A");
        assert_eq!(variable_const_name(Variable::Lambda), "LAMBDA");
        assert_eq!(variable_const_name(Variable::K), "K");
        assert_eq!(variable_const_name(Variable::H), "H");
        assert_eq!(variable_const_name(Variable::Q), "Q");
        assert_eq!(variable_const_name(Variable::P), "P");
    }

    #[test]
    fn test_variable_description() {
        assert!(variable_description(Variable::A).contains("Semi-major"));
        assert!(variable_description(Variable::Lambda).contains("Mean longitude"));
        assert!(variable_description(Variable::K).contains("cos"));
        assert!(variable_description(Variable::H).contains("sin"));
        assert!(variable_description(Variable::Q).contains("node"));
        assert!(variable_description(Variable::P).contains("node"));
    }

    #[test]
    fn test_planet_module_name() {
        assert_eq!(planet_module_name(1), "mercury");
        assert_eq!(planet_module_name(2), "venus");
        assert_eq!(planet_module_name(3), "emb");
        assert_eq!(planet_module_name(4), "mars");
        assert_eq!(planet_module_name(5), "jupiter");
        assert_eq!(planet_module_name(6), "saturn");
        assert_eq!(planet_module_name(7), "uranus");
        assert_eq!(planet_module_name(8), "neptune");
        assert_eq!(planet_module_name(9), "pluto");
        assert_eq!(planet_module_name(99), "unknown");
    }

    #[test]
    fn test_module_name_to_planet() {
        assert_eq!(module_name_to_planet("mercury"), Some(1));
        assert_eq!(module_name_to_planet("venus"), Some(2));
        assert_eq!(module_name_to_planet("emb"), Some(3));
        assert_eq!(module_name_to_planet("mars"), Some(4));
        assert_eq!(module_name_to_planet("jupiter"), Some(5));
        assert_eq!(module_name_to_planet("saturn"), Some(6));
        assert_eq!(module_name_to_planet("uranus"), Some(7));
        assert_eq!(module_name_to_planet("neptune"), Some(8));
        assert_eq!(module_name_to_planet("pluto"), Some(9));
        assert_eq!(module_name_to_planet("unknown"), None);
    }

    #[test]
    fn test_generate_mod_rs() {
        let planets = vec![3, 9];
        let source = generate_mod_rs(&planets);
        assert!(source.contains("pub mod emb;"));
        assert!(source.contains("pub mod pluto;"));
        assert!(source.contains("pub struct Term"));
        assert!(source.contains("pub struct TimeBlock"));
        assert!(source.contains("VSOP2013 planetary coefficients"));
    }

    #[test]
    fn test_generate_mod_rs_all_planets() {
        let planets: Vec<u8> = (1..=9).collect();
        let source = generate_mod_rs(&planets);
        assert!(source.contains("pub mod mercury;"));
        assert!(source.contains("pub mod venus;"));
        assert!(source.contains("pub mod emb;"));
        assert!(source.contains("pub mod mars;"));
        assert!(source.contains("pub mod jupiter;"));
        assert!(source.contains("pub mod saturn;"));
        assert!(source.contains("pub mod uranus;"));
        assert!(source.contains("pub mod neptune;"));
        assert!(source.contains("pub mod pluto;"));
    }

    #[test]
    fn test_filtered_term_from_vsop_term() {
        let vsop_term = make_term_with_mults(
            3.0,
            4.0,
            [1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        );
        let filtered: FilteredTerm = (&vsop_term).into();

        assert!((filtered.s_coeff - 3.0).abs() < 1e-10);
        assert!((filtered.c_coeff - 4.0).abs() < 1e-10);
        assert!((filtered.amplitude - 5.0).abs() < 1e-10);
        assert_eq!(filtered.multipliers[0], 1);
        assert_eq!(filtered.multipliers[1], 2);
        assert_eq!(filtered.multipliers[2], 3);
    }

    #[test]
    fn test_filter_and_group_terms() {
        let vsop = make_test_vsop();
        let var_data = filter_and_group_terms(&vsop, 0.5);

        // Should have 2 variables: A and Lambda (K and others are empty)
        assert_eq!(var_data.len(), 2);

        let a_data = var_data.iter().find(|v| v.variable == Variable::A).unwrap();
        assert_eq!(a_data.total_terms, 4);
        assert_eq!(a_data.retained_terms, 3); // 5, 1, 10 above threshold 0.5

        // Check blocks are sorted by power
        assert_eq!(a_data.blocks.len(), 2);
        assert_eq!(a_data.blocks[0].power, 0);
        assert_eq!(a_data.blocks[1].power, 1);

        // Check terms are sorted by amplitude descending
        let t0_terms = &a_data.blocks[0].terms;
        assert!(t0_terms[0].amplitude > t0_terms[1].amplitude);
    }

    #[test]
    fn test_filter_and_group_terms_high_threshold() {
        let vsop = make_test_vsop();
        let var_data = filter_and_group_terms(&vsop, 100.0);

        // No terms above threshold, so no variable data
        assert!(var_data.is_empty());
    }

    #[test]
    fn test_generate_planet_source() {
        let vsop = make_test_vsop();
        let (source, retained, total) = generate_planet_source(3, &vsop, 0.5);

        assert!(source.contains("VSOP2013 coefficients for"));
        assert!(source.contains("Earth-Moon Barycenter"));
        assert!(source.contains("Threshold:"));
        assert!(source.contains("pub const A:"));
        assert!(source.contains("pub const LAMBDA:"));
        assert!(source.contains("TimeBlock {"));
        assert!(source.contains("Term {"));

        assert_eq!(total, 5); // 4 in A + 1 in Lambda
        assert!(retained <= total);
    }

    #[test]
    fn test_generate_planet_source_contains_multipliers() {
        let mut mults = [0i32; 17];
        mults[0] = 1;
        mults[1] = -2;
        let vsop = Vsop2013File {
            planet: 1,
            blocks: vec![make_block(
                Variable::A,
                0,
                vec![make_term_with_mults(10.0, 0.0, mults)],
            )],
        };

        let (source, _, _) = generate_planet_source(1, &vsop, 0.0);
        assert!(source.contains("1,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"));
    }

    #[test]
    fn test_generate_planet_source_empty_vsop() {
        let vsop = Vsop2013File {
            planet: 5,
            blocks: vec![],
        };

        let (source, retained, total) = generate_planet_source(5, &vsop, 0.0);
        assert!(source.contains("Jupiter"));
        assert_eq!(retained, 0);
        assert_eq!(total, 0);
    }

    #[test]
    fn test_discover_existing_planets_empty_dir() {
        let temp_dir = TempDir::new().unwrap();
        let planets = discover_existing_planets(temp_dir.path());
        assert!(planets.is_empty());
    }

    #[test]
    fn test_discover_existing_planets_with_files() {
        let temp_dir = TempDir::new().unwrap();

        // Create some planet files
        fs::write(temp_dir.path().join("mercury.rs"), "").unwrap();
        fs::write(temp_dir.path().join("venus.rs"), "").unwrap();
        fs::write(temp_dir.path().join("mod.rs"), "").unwrap(); // Should be ignored

        let planets = discover_existing_planets(temp_dir.path());
        assert!(planets.contains(&1)); // mercury
        assert!(planets.contains(&2)); // venus
        assert!(!planets.contains(&3)); // emb not present
    }

    #[test]
    fn test_discover_existing_planets_ignores_non_rs_files() {
        let temp_dir = TempDir::new().unwrap();

        fs::write(temp_dir.path().join("mercury.rs"), "").unwrap();
        fs::write(temp_dir.path().join("venus.txt"), "").unwrap(); // Should be ignored
        fs::write(temp_dir.path().join("mars"), "").unwrap(); // Should be ignored

        let planets = discover_existing_planets(temp_dir.path());
        assert!(planets.contains(&1));
        assert!(!planets.contains(&2)); // venus.txt ignored
        assert!(!planets.contains(&4)); // mars (no extension) ignored
    }

    #[test]
    fn test_generate_planet_module_creates_file() {
        let temp_dir = TempDir::new().unwrap();
        let vsop = make_test_vsop();
        let config = GenerateConfig {
            threshold: 0.0,
            output_dir: temp_dir.path().to_path_buf(),
        };

        generate_planet_module(3, &vsop, &config, temp_dir.path()).unwrap();

        let file_path = temp_dir.path().join("emb.rs");
        assert!(file_path.exists());

        let content = fs::read_to_string(file_path).unwrap();
        assert!(content.contains("VSOP2013"));
    }

    #[test]
    fn test_generate_all_creates_mod_and_planet_files() {
        let temp_dir = TempDir::new().unwrap();
        let vsop1 = Vsop2013File {
            planet: 1,
            blocks: vec![make_block(Variable::A, 0, vec![make_term(1.0, 0.0)])],
        };
        let vsop3 = Vsop2013File {
            planet: 3,
            blocks: vec![make_block(Variable::A, 0, vec![make_term(1.0, 0.0)])],
        };

        let config = GenerateConfig {
            threshold: 0.0,
            output_dir: temp_dir.path().to_path_buf(),
        };

        generate_all(&[(1, vsop1), (3, vsop3)], &config).unwrap();

        assert!(temp_dir.path().join("mod.rs").exists());
        assert!(temp_dir.path().join("mercury.rs").exists());
        assert!(temp_dir.path().join("emb.rs").exists());

        let mod_content = fs::read_to_string(temp_dir.path().join("mod.rs")).unwrap();
        assert!(mod_content.contains("pub mod mercury;"));
        assert!(mod_content.contains("pub mod emb;"));
    }

    #[test]
    fn test_generate_all_preserves_existing_planets() {
        let temp_dir = TempDir::new().unwrap();

        // Create an existing planet file
        fs::write(temp_dir.path().join("venus.rs"), "// existing").unwrap();

        let vsop = Vsop2013File {
            planet: 3,
            blocks: vec![make_block(Variable::A, 0, vec![make_term(1.0, 0.0)])],
        };

        let config = GenerateConfig {
            threshold: 0.0,
            output_dir: temp_dir.path().to_path_buf(),
        };

        generate_all(&[(3, vsop)], &config).unwrap();

        let mod_content = fs::read_to_string(temp_dir.path().join("mod.rs")).unwrap();
        // Should include both venus (existing) and emb (newly generated)
        assert!(mod_content.contains("pub mod venus;"));
        assert!(mod_content.contains("pub mod emb;"));
    }

    #[test]
    fn test_generate_config_struct() {
        let config = GenerateConfig {
            threshold: 1e-10,
            output_dir: std::path::PathBuf::from("/tmp/test"),
        };
        assert_eq!(config.threshold, 1e-10);
    }

    #[test]
    fn test_generate_planet_source_time_power_comments() {
        let vsop = Vsop2013File {
            planet: 1,
            blocks: vec![
                make_block(Variable::A, 0, vec![make_term(10.0, 0.0)]),
                make_block(Variable::A, 2, vec![make_term(5.0, 0.0)]),
            ],
        };

        let (source, _, _) = generate_planet_source(1, &vsop, 0.0);
        assert!(source.contains("// T^0 terms"));
        assert!(source.contains("// T^2 terms"));
    }

    #[test]
    fn test_generate_planet_module_empty_vsop() {
        // This test exercises the else branch at line 321 (percentage = 0.0 when total == 0)
        let temp_dir = TempDir::new().unwrap();
        let vsop = Vsop2013File {
            planet: 7,
            blocks: vec![],
        };
        let config = GenerateConfig {
            threshold: 0.0,
            output_dir: temp_dir.path().to_path_buf(),
        };

        let result = generate_planet_module(7, &vsop, &config, temp_dir.path());
        assert!(result.is_ok());

        let file_path = temp_dir.path().join("uranus.rs");
        assert!(file_path.exists());

        let content = fs::read_to_string(file_path).unwrap();
        assert!(content.contains("Uranus"));
        assert!(content.contains("0 of 0")); // Terms retained: 0 of 0
    }
}
