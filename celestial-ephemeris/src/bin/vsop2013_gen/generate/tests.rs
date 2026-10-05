use super::source::*;
use super::*;
use crate::parser::{Variable, Vsop2013Block, Vsop2013Header, Vsop2013Term};
use std::path::Path;
use tempfile::TempDir;

fn term(s_coeff: f64, c_coeff: f64, multipliers: [i32; 17]) -> Vsop2013Term {
    Vsop2013Term {
        multipliers,
        s_coeff,
        c_coeff,
    }
}

fn block(variable: Variable, time_power: u8, terms: Vec<Vsop2013Term>) -> Vsop2013Block {
    Vsop2013Block {
        header: Vsop2013Header {
            planet: 3,
            variable,
            time_power,
            term_count: terms.len() as u32,
        },
        terms,
    }
}

fn single_term(planet: u8, multipliers: [i32; 17]) -> Vsop2013File {
    Vsop2013File {
        planet,
        blocks: vec![block(Variable::A, 0, vec![term(1.0, 0.0, multipliers)])],
    }
}

// Amplitudes 5, 1 and 0.05 at T^0 and 10 at T^1 for A; 50 for Lambda.
fn sample() -> Vsop2013File {
    let t = |s, c| term(s, c, [0; 17]);
    Vsop2013File {
        planet: 3,
        blocks: vec![
            block(
                Variable::A,
                0,
                vec![t(3.0, 4.0), t(0.6, 0.8), t(0.03, 0.04)],
            ),
            block(Variable::A, 1, vec![t(6.0, 8.0)]),
            block(Variable::Lambda, 0, vec![t(30.0, 40.0)]),
        ],
    }
}

fn config(output_dir: &Path) -> GenerateConfig {
    GenerateConfig {
        threshold: Some(0.0),
        output_dir: output_dir.to_path_buf(),
    }
}

#[test]
fn floats_are_shortest_round_trip() {
    assert_eq!(format_float(0.0), "0.0");
    assert_eq!(format_float(-0.0), "0.0");
    assert_eq!(format_float(1.5e-10), "1.5e-10");
    assert_eq!(format_float(-123.456789), "-1.23456789e2");
    for x in [0.1 + 0.2, -3.1072908423321636e-5, f64::MIN_POSITIVE] {
        assert_eq!(format_float(x).parse::<f64>().unwrap(), x);
    }
}

#[test]
fn lists_format() {
    assert_eq!(format_list(&[0, 1, -2, 72490]), "[0,1,-2,72490]");
    assert_eq!(format_list::<u8>(&[]), "[]");
}

#[test]
fn variable_names() {
    let expected = [
        ("A", "Semi-major axis (A)"),
        ("LAMBDA", "Mean longitude (Lambda)"),
        ("K", "e*cos(perihelion) (K)"),
        ("H", "e*sin(perihelion) (H)"),
        ("Q", "sin(i/2)*cos(node) (Q)"),
        ("P", "sin(i/2)*sin(node) (P)"),
    ];
    for (variable, (name, description)) in Variable::ALL.into_iter().zip(expected) {
        assert_eq!(variable_const_name(variable), name);
        assert_eq!(variable_description(variable), description);
    }
}

#[test]
fn module_names() {
    let names = [
        "mercury", "venus", "emb", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto",
    ];
    for (i, name) in names.into_iter().enumerate() {
        let planet = i as u8 + 1;
        assert_eq!(planet_module_name(planet), name);
        assert_eq!(module_name_to_planet(name), Some(planet));
    }
    assert_eq!(planet_module_name(0), "unknown");
    assert_eq!(planet_module_name(99), "unknown");
    assert_eq!(module_name_to_planet("unknown"), None);
    assert_eq!(module_name_to_planet("mod"), None);
}

#[test]
fn mod_rs_lists_the_planets_and_types() {
    let source = generate_mod_rs(&[3, 9]);
    assert!(source.starts_with("//! VSOP2013 planetary coefficients\n"));
    assert!(source.contains("\npub(crate) mod emb;\npub(crate) mod pluto;\n\n"));
    assert!(source.contains("pub(crate) struct Term {"));
    assert!(source.contains("    pub(crate) mult: [i32; 6],\n    pub(crate) index: [u8; 6],\n"));
    assert!(source.contains("pub(crate) struct TimeBlock {"));
    let all = generate_mod_rs(&(1..=9).collect::<Vec<u8>>());
    for planet in 1..=9 {
        let line = format!("pub(crate) mod {};\n", planet_module_name(planet));
        assert!(all.contains(&line), "{line}");
    }
}

#[test]
fn filtered_terms_carry_the_amplitude() {
    let mults = [1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3];
    let filtered = FilteredTerm::try_from(&term(3.0, 4.0, mults)).unwrap();
    assert_eq!((filtered.s_coeff, filtered.c_coeff), (3.0, 4.0));
    assert_eq!(filtered.amplitude, 5.0);
    assert_eq!(filtered.mult, [1, 2, -3, 0, 0, 0]);
    assert_eq!(filtered.index, [0, 3, 16, 0, 0, 0]);
}

#[test]
fn terms_fill_every_slot_and_no_more() {
    let mut mults = [0; 17];
    for i in [0, 2, 4, 6, 8, 10] {
        mults[i] = i as i32 + 1;
    }
    let filtered = FilteredTerm::try_from(&term(1.0, 0.0, mults)).unwrap();
    assert_eq!(filtered.mult, [1, 3, 5, 7, 9, 11]);
    assert_eq!(filtered.index, [0, 2, 4, 6, 8, 10]);
    mults[16] = 1;
    let error = FilteredTerm::try_from(&term(1.0, 0.0, mults))
        .err()
        .unwrap();
    assert_eq!(
        error,
        "a term uses more than 6 arguments: [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 0, 0, 0, 0, 1]"
    );
    let vsop = single_term(3, mults);
    assert_eq!(generate_planet_source(3, &vsop, 0.0).err().unwrap(), error);
    let dir = TempDir::new().unwrap();
    assert_eq!(generate_all(&[(3, vsop)], &config(dir.path())), Err(error));
}

#[test]
fn terms_are_cut_and_grouped_by_power() {
    let data = filter_and_group_terms(&sample(), 0.5).unwrap();
    let variables: Vec<Variable> = data.iter().map(|v| v.variable).collect();
    assert_eq!(variables, [Variable::A, Variable::Lambda]);
    let a = &data[0];
    assert_eq!((a.total_terms, a.retained_terms), (4, 3));
    let powers: Vec<u8> = a.blocks.iter().map(|b| b.power).collect();
    assert_eq!(powers, [0, 1]);
    let amplitudes: Vec<f64> = a.blocks[0].terms.iter().map(|t| t.amplitude).collect();
    assert_eq!(amplitudes, [5.0, 1.0]);
    assert!(filter_and_group_terms(&sample(), 100.0).unwrap().is_empty());
}

#[test]
fn planet_source() {
    let (source, retained, total) = generate_planet_source(3, &sample(), 0.5).unwrap();
    assert_eq!((retained, total), (4, 5));
    for expected in [
        "pub(crate) const A: &[TimeBlock] = &[",
        "pub(crate) const LAMBDA: &[TimeBlock] = &[",
        "    // T^0 terms\n    TimeBlock {\n        power: 0,\n        terms: &[\n",
        "    // T^1 terms\n",
        "            Term { s: 3e0, c: 4e0, mult: [0,0,0,0,0,0], index: [0,0,0,0,0,0] },\n",
    ] {
        assert!(source.contains(expected), "{expected}");
    }
}

#[test]
fn header_records_the_command_and_no_date() {
    let (source, _, _) = generate_planet_source(3, &sample(), 0.5).unwrap();
    let header: Vec<&str> = source
        .lines()
        .take_while(|l| !l.starts_with("use "))
        .collect();
    assert_eq!(
        header,
        [
            "//! VSOP2013 coefficients for Earth-Moon Barycenter",
            "//!",
            "//! Generated from VSOP2013p3.dat by",
            "//! `vsop2013-gen generate --input <dir> --output <dir> --planet 3 --threshold 5e-1`",
            "//! Terms retained: 4 of 5 (80.0%)",
            "",
        ]
    );
}

#[test]
fn default_thresholds_are_the_shipped_ones() {
    for planet in 1..=9 {
        let path = format!(
            "{}/src/planetary_coefficients/{}.rs",
            env!("CARGO_MANIFEST_DIR"),
            planet_module_name(planet)
        );
        let shipped = fs::read_to_string(path).unwrap();
        let command = format!(
            "--planet {} --threshold {:e}`",
            planet,
            default_threshold(planet)
        );
        assert!(shipped.contains(&command), "{}", command);
    }
}

#[test]
fn multipliers_beyond_i16_are_kept() {
    let mut mults = [0i32; 17];
    mults[13] = 72490;
    let (source, _, _) = generate_planet_source(9, &single_term(9, mults), 1e-7).unwrap();
    assert!(source.contains("mult: [72490,0,0,0,0,0], index: [13,0,0,0,0,0]"));
}

#[test]
fn empty_series() {
    let vsop = Vsop2013File {
        planet: 5,
        blocks: vec![],
    };
    let (source, retained, total) = generate_planet_source(5, &vsop, 0.0).unwrap();
    assert_eq!((retained, total), (0, 0));
    assert!(source.contains("//! VSOP2013 coefficients for Jupiter\n"));
    assert!(source.contains("//! Terms retained: 0 of 0 (0.0%)\n"));
}

#[test]
fn existing_planet_files_are_discovered() {
    let dir = TempDir::new().unwrap();
    assert_eq!(discover_existing_planets(dir.path()), Ok(BTreeSet::new()));
    for name in "mercury.rs venus.rs mod.rs mars.txt saturn other.rs".split(' ') {
        fs::write(dir.path().join(name), "").unwrap();
    }
    assert_eq!(
        discover_existing_planets(dir.path()),
        Ok(BTreeSet::from([1, 2]))
    );
}

#[test]
fn missing_output_dir_is_an_error() {
    let dir = TempDir::new().unwrap();
    let missing = dir.path().join("missing");
    let error = discover_existing_planets(&missing).unwrap_err();
    assert!(
        error.starts_with(&format!("Failed to read {}: ", missing.display())),
        "{error}"
    );
}

#[test]
fn planet_module_is_written() {
    let dir = TempDir::new().unwrap();
    let report = generate_planet_module(3, &sample(), &config(dir.path())).unwrap();
    assert_eq!(
        report,
        "  Generated emb.rs (5 terms of 5 = 100.0%, threshold 0e0)"
    );
    let content = fs::read_to_string(dir.path().join("emb.rs")).unwrap();
    assert_eq!(
        content,
        generate_planet_source(3, &sample(), 0.0).unwrap().0
    );
}

#[test]
fn generate_all_writes_planets_and_keeps_existing_ones() {
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("venus.rs"), "// existing").unwrap();
    let files = [(1, single_term(1, [0; 17])), (3, single_term(3, [0; 17]))];
    generate_all(&files, &config(dir.path())).unwrap();
    assert!(dir.path().join("mercury.rs").exists());
    assert!(dir.path().join("emb.rs").exists());
    let mod_rs = fs::read_to_string(dir.path().join("mod.rs")).unwrap();
    assert_eq!(mod_rs, generate_mod_rs(&[1, 2, 3]));
}

#[test]
fn write_errors_name_the_step() {
    let dir = TempDir::new().unwrap();
    let files = [(3, single_term(3, [0; 17]))];
    let error = |output_dir: &Path| generate_all(&files, &config(output_dir)).unwrap_err();
    let file = dir.path().join("file");
    fs::write(&file, "").unwrap();
    assert!(error(&file).starts_with("Failed to create output dir: "));
    fs::create_dir(dir.path().join("mod.rs")).unwrap();
    assert!(error(dir.path()).starts_with("Failed to write mod.rs: "));
    fs::remove_dir(dir.path().join("mod.rs")).unwrap();
    let planet = dir.path().join("emb.rs");
    fs::create_dir(&planet).unwrap();
    let prefix = format!("Failed to write {}: ", planet.display());
    assert!(error(dir.path()).starts_with(&prefix));
}

// Write and search but no read: the planets already there can't be listed,
// and a mod.rs written anyway would drop them.
#[cfg(unix)]
#[test]
fn unreadable_output_dir_is_an_error() {
    use std::os::unix::fs::PermissionsExt;
    let dir = TempDir::new().unwrap();
    fs::write(dir.path().join("venus.rs"), "// existing").unwrap();
    fs::set_permissions(dir.path(), fs::Permissions::from_mode(0o300)).unwrap();
    let result = generate_all(&[(3, single_term(3, [0; 17]))], &config(dir.path()));
    fs::set_permissions(dir.path(), fs::Permissions::from_mode(0o700)).unwrap();
    let error = result.unwrap_err();
    assert!(error.starts_with("Failed to read "), "{error}");
    assert!(!dir.path().join("mod.rs").exists());
}
