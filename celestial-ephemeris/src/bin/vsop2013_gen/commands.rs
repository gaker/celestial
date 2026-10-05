use std::path::{Path, PathBuf};

use crate::analyze::{analyze_file, format_analysis, format_summary_table};
use crate::download::{
    default_client, download_all, download_planet, find_planet_files, planet_filename,
};
use crate::generate::{default_threshold, generate_all, GenerateConfig};
use crate::parser::{parse_file, planet_name, Vsop2013File};

fn parse_planet_arg(arg: &Option<String>) -> Result<Option<u8>, String> {
    match arg {
        None => Ok(None),
        Some(s) if s == "all" => Ok(None),
        Some(s) => {
            let p: u8 = s.parse().map_err(|_| format!("Invalid planet: {}", s))?;
            if !(1..=9).contains(&p) {
                return Err(format!("Planet must be 1-9, got {}", p));
            }
            Ok(Some(p))
        }
    }
}

fn input_files(input: &Path, planet: Option<u8>) -> Result<Vec<(u8, PathBuf)>, String> {
    let Some(p) = planet else {
        let found = find_planet_files(input);
        if found.is_empty() {
            return Err(format!("No VSOP2013 files found in {}", input.display()));
        }
        return Ok(found);
    };
    let path = input.join(planet_filename(p));
    if !path.exists() {
        return Err(format!("File not found: {}", path.display()));
    }
    Ok(vec![(p, path)])
}

pub(crate) fn cmd_download(output: PathBuf, planet_arg: Option<String>) -> Result<(), String> {
    std::fs::create_dir_all(&output).map_err(|e| format!("Failed to create output dir: {}", e))?;
    let planet = parse_planet_arg(&planet_arg)?;
    let client = default_client()?;
    match planet {
        None => download_all(&client, &output),
        Some(p) => {
            println!("Downloading planet {} ({})...", p, planet_name(p));
            download_planet(&client, p, &output)
        }
    }
}

pub(crate) fn cmd_analyze(
    input: PathBuf,
    planet_arg: Option<String>,
    threshold: Option<f64>,
) -> Result<(), String> {
    let files = input_files(&input, parse_planet_arg(&planet_arg)?)?;
    let mut analyses = Vec::new();
    for (p, path) in &files {
        println!("Parsing planet {} ({})...", p, planet_name(*p));
        let analysis = analyze_file(path, threshold.unwrap_or_else(|| default_threshold(*p)))?;
        print!("{}", format_analysis(&analysis));
        analyses.push(analysis);
    }
    if analyses.len() > 1 {
        print!("{}", format_summary_table(&analyses));
    }
    Ok(())
}

fn parse_all(files: &[(u8, PathBuf)]) -> Result<Vec<(u8, Vsop2013File)>, String> {
    println!("Parsing VSOP2013 files...");
    let mut vsop_files = Vec::new();
    for (p, path) in files {
        println!("  Parsing planet {} ({})...", p, planet_name(*p));
        let vsop = parse_file(path).map_err(|e| format!("Parse error: {}", e))?;
        vsop_files.push((*p, vsop));
    }
    Ok(vsop_files)
}

pub(crate) fn cmd_generate(
    input: PathBuf,
    output: PathBuf,
    threshold: Option<f64>,
    planet_arg: Option<String>,
) -> Result<(), String> {
    let files = input_files(&input, parse_planet_arg(&planet_arg)?)?;
    let vsop_files = parse_all(&files)?;
    let config = GenerateConfig {
        threshold,
        output_dir: output,
    };
    println!("\nGenerating Rust code...");
    generate_all(&vsop_files, &config)?;
    println!("\nGeneration complete!");
    println!("Output directory: {}", config.output_dir.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    fn write_vsop_file(dir: &Path, planet: u8) {
        let content = format!(" VSOP2013  {}  1  0  2    PLANET VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
    2   1  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.2000000000000000 +00  0.0000000000000000 +00
", planet);
        fs::write(dir.join(planet_filename(planet)), content).unwrap();
    }

    #[test]
    fn planet_args() {
        assert_eq!(parse_planet_arg(&None), Ok(None));
        assert_eq!(parse_planet_arg(&Some("all".to_string())), Ok(None));
        for p in 1..=9 {
            assert_eq!(parse_planet_arg(&Some(p.to_string())), Ok(Some(p)));
        }
        let cases = [
            ("0", "Planet must be 1-9, got 0"),
            ("10", "Planet must be 1-9, got 10"),
            ("mars", "Invalid planet: mars"),
            ("300", "Invalid planet: 300"),
        ];
        for (arg, message) in cases {
            assert_eq!(
                parse_planet_arg(&Some(arg.to_string())),
                Err(message.to_string())
            );
        }
    }

    #[test]
    fn download_creates_the_output_dir_before_checking_the_planet() {
        let dir = TempDir::new().unwrap();
        let output = dir.path().join("new_subdir");
        let result = cmd_download(output.clone(), Some("invalid".to_string()));
        assert_eq!(result, Err("Invalid planet: invalid".to_string()));
        assert!(output.exists());
    }

    #[test]
    fn missing_or_bad_inputs() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().to_path_buf();
        let output = dir.path().join("output");
        let cases = [
            (
                None,
                format!("No VSOP2013 files found in {}", input.display()),
            ),
            (
                Some("3"),
                format!("File not found: {}", input.join("VSOP2013p3.dat").display()),
            ),
            (Some("bad"), "Invalid planet: bad".to_string()),
            (Some("15"), "Planet must be 1-9, got 15".to_string()),
        ];
        for (planet, message) in cases {
            let planet = planet.map(str::to_string);
            let analyze = cmd_analyze(input.clone(), planet.clone(), Some(1e-10));
            let generate = cmd_generate(input.clone(), output.clone(), Some(1e-10), planet);
            assert_eq!(analyze, Err(message.clone()));
            assert_eq!(generate, Err(message));
        }
        assert!(!output.exists());
    }

    #[test]
    fn analyzes_one_or_all_planets() {
        let dir = TempDir::new().unwrap();
        for planet in [1, 3, 5] {
            write_vsop_file(dir.path(), planet);
        }
        let input = dir.path().to_path_buf();
        assert_eq!(
            cmd_analyze(input.clone(), Some("3".to_string()), Some(1e-10)),
            Ok(())
        );
        assert_eq!(cmd_analyze(input, None, None), Ok(()));
    }

    #[test]
    fn generates_one_planet() {
        let dir = TempDir::new().unwrap();
        write_vsop_file(dir.path(), 5);
        let output = dir.path().join("output");
        let result = cmd_generate(
            dir.path().to_path_buf(),
            output.clone(),
            Some(1e-10),
            Some("5".to_string()),
        );
        assert_eq!(result, Ok(()));
        assert!(output.join("mod.rs").exists());
        let jupiter = fs::read_to_string(output.join("jupiter.rs")).unwrap();
        assert!(jupiter.contains("--planet 5 --threshold 1e-10`"));
    }

    #[test]
    fn generates_all_planets_at_their_default_cuts() {
        let dir = TempDir::new().unwrap();
        for planet in [1, 3] {
            write_vsop_file(dir.path(), planet);
        }
        let output = dir.path().join("output");
        assert_eq!(
            cmd_generate(dir.path().to_path_buf(), output.clone(), None, None),
            Ok(())
        );
        let mercury = fs::read_to_string(output.join("mercury.rs")).unwrap();
        assert!(mercury.contains("--planet 1 --threshold 1e-9`"));
        let emb = fs::read_to_string(output.join("emb.rs")).unwrap();
        assert!(emb.contains("--planet 3 --threshold 1e-10`"));
    }
}
