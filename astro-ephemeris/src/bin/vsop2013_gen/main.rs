#[cfg(feature = "cli")]
use clap::{Parser, Subcommand};

mod analyze;
mod download;
mod generate;
mod parser;

use analyze::{analyze_file, print_analysis, print_summary_table};
use download::{default_client, download_all, download_planet, find_planet_files};
use generate::{generate_all, GenerateConfig};
use parser::{parse_file, planet_name, Vsop2013File};
use std::path::PathBuf;

#[cfg(feature = "cli")]
#[derive(Parser)]
#[command(name = "vsop2013-gen")]
#[command(about = "VSOP2013 ephemeris data processor and Rust code generator")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[cfg(feature = "cli")]
#[derive(Subcommand)]
enum Commands {
    Download {
        #[arg(short, long, default_value = "./vsop2013")]
        output: PathBuf,
        #[arg(short, long, help = "Planet number (1-9) or 'all'")]
        planet: Option<String>,
    },
    Analyze {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long, help = "Planet number (1-9) or 'all'")]
        planet: Option<String>,
        #[arg(short, long, default_value = "1e-10")]
        threshold: f64,
    },
    Generate {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long)]
        threshold: f64,
        #[arg(short, long, help = "Planet number (1-9) or 'all'")]
        planet: Option<String>,
    },
}

#[cfg(feature = "cli")]
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

#[cfg(feature = "cli")]
fn cmd_download(output: PathBuf, planet_arg: Option<String>) -> Result<(), String> {
    std::fs::create_dir_all(&output).map_err(|e| format!("Failed to create output dir: {}", e))?;

    let planet = parse_planet_arg(&planet_arg)?;
    let client = default_client()?;

    match planet {
        None => download_all(&client, &output),
        Some(p) => {
            println!("Downloading planet {} ({})...", p, planet_name(p));
            download_planet(p, &output)
        }
    }
}

#[cfg(feature = "cli")]
fn cmd_analyze(input: PathBuf, planet_arg: Option<String>, threshold: f64) -> Result<(), String> {
    let planet = parse_planet_arg(&planet_arg)?;

    let files = match planet {
        Some(p) => {
            let path = input.join(download::planet_filename(p));
            if !path.exists() {
                return Err(format!("File not found: {}", path.display()));
            }
            vec![(p, path)]
        }
        None => {
            let found = find_planet_files(&input);
            if found.is_empty() {
                return Err(format!("No VSOP2013 files found in {}", input.display()));
            }
            found
        }
    };

    let mut analyses = Vec::new();

    for (p, path) in &files {
        println!("Parsing planet {} ({})...", p, planet_name(*p));
        let analysis = analyze_file(path, threshold)?;
        print_analysis(&analysis, threshold);
        analyses.push(analysis);
    }

    if analyses.len() > 1 {
        print_summary_table(&analyses, threshold);
    }

    Ok(())
}

#[cfg(feature = "cli")]
fn cmd_generate(
    input: PathBuf,
    output: PathBuf,
    threshold: f64,
    planet_arg: Option<String>,
) -> Result<(), String> {
    let planet = parse_planet_arg(&planet_arg)?;

    let files = match planet {
        Some(p) => {
            let path = input.join(download::planet_filename(p));
            if !path.exists() {
                return Err(format!("File not found: {}", path.display()));
            }
            vec![(p, path)]
        }
        None => {
            let found = find_planet_files(&input);
            if found.is_empty() {
                return Err(format!("No VSOP2013 files found in {}", input.display()));
            }
            found
        }
    };

    println!("Parsing VSOP2013 files...");
    let mut vsop_files: Vec<(u8, Vsop2013File)> = Vec::new();
    for (p, path) in &files {
        println!("  Parsing planet {} ({})...", p, planet_name(*p));
        let vsop = parse_file(path).map_err(|e| format!("Parse error: {}", e))?;
        vsop_files.push((*p, vsop));
    }

    let config = GenerateConfig {
        threshold,
        output_dir: output.clone(),
    };

    println!("\nGenerating Rust code (threshold: {:.0e})...", threshold);
    generate_all(&vsop_files, &config)?;

    println!("\nGeneration complete!");
    println!("Output directory: {}", output.display());
    Ok(())
}

#[cfg(feature = "cli")]
fn main() {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Download { output, planet } => cmd_download(output, planet),
        Commands::Analyze {
            input,
            planet,
            threshold,
        } => cmd_analyze(input, planet, threshold),
        Commands::Generate {
            input,
            output,
            threshold,
            planet,
        } => cmd_generate(input, output, threshold, planet),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

#[cfg(not(feature = "cli"))]
fn main() {
    eprintln!("vsop2013-gen requires the 'cli' feature.");
    eprintln!("Run with: cargo run --features cli --bin vsop2013-gen -- <args>");
    std::process::exit(1);
}

#[cfg(all(test, feature = "cli"))]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn test_parse_planet_arg_none() {
        let result = parse_planet_arg(&None);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), None);
    }

    #[test]
    fn test_parse_planet_arg_all() {
        let result = parse_planet_arg(&Some("all".to_string()));
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), None);
    }

    #[test]
    fn test_parse_planet_arg_valid_planets() {
        for p in 1..=9 {
            let result = parse_planet_arg(&Some(p.to_string()));
            assert!(result.is_ok(), "Failed for planet {}", p);
            assert_eq!(result.unwrap(), Some(p));
        }
    }

    #[test]
    fn test_parse_planet_arg_invalid_planet_zero() {
        let result = parse_planet_arg(&Some("0".to_string()));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Planet must be 1-9"));
    }

    #[test]
    fn test_parse_planet_arg_invalid_planet_ten() {
        let result = parse_planet_arg(&Some("10".to_string()));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Planet must be 1-9"));
    }

    #[test]
    fn test_parse_planet_arg_invalid_string() {
        let result = parse_planet_arg(&Some("mars".to_string()));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Invalid planet"));
    }

    #[test]
    fn test_cmd_download_creates_output_dir() {
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().join("new_subdir");

        assert!(!output_path.exists());

        // This will fail because we can't actually download in tests,
        // but it should at least create the directory first
        let _ = cmd_download(output_path.clone(), Some("invalid".to_string()));

        // Directory should have been created before the invalid planet error
        assert!(output_path.exists());
    }

    #[test]
    fn test_cmd_download_invalid_planet() {
        let temp_dir = TempDir::new().unwrap();
        let result = cmd_download(temp_dir.path().to_path_buf(), Some("invalid".to_string()));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Invalid planet"));
    }

    #[test]
    fn test_cmd_analyze_no_files() {
        let temp_dir = TempDir::new().unwrap();
        let result = cmd_analyze(temp_dir.path().to_path_buf(), None, 1e-10);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("No VSOP2013 files found"));
    }

    #[test]
    fn test_cmd_analyze_specific_planet_not_found() {
        let temp_dir = TempDir::new().unwrap();
        let result = cmd_analyze(temp_dir.path().to_path_buf(), Some("3".to_string()), 1e-10);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("File not found"));
    }

    #[test]
    fn test_cmd_analyze_invalid_planet() {
        let temp_dir = TempDir::new().unwrap();
        let result = cmd_analyze(
            temp_dir.path().to_path_buf(),
            Some("invalid".to_string()),
            1e-10,
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Invalid planet"));
    }

    #[test]
    fn test_cmd_generate_no_files() {
        let temp_dir = TempDir::new().unwrap();
        let input = temp_dir.path().join("input");
        let output = temp_dir.path().join("output");
        fs::create_dir_all(&input).unwrap();

        let result = cmd_generate(input, output, 1e-10, None);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("No VSOP2013 files found"));
    }

    #[test]
    fn test_cmd_generate_specific_planet_not_found() {
        let temp_dir = TempDir::new().unwrap();
        let input = temp_dir.path().join("input");
        let output = temp_dir.path().join("output");
        fs::create_dir_all(&input).unwrap();

        let result = cmd_generate(input, output, 1e-10, Some("5".to_string()));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("File not found"));
    }

    #[test]
    fn test_cmd_generate_invalid_planet() {
        let temp_dir = TempDir::new().unwrap();
        let input = temp_dir.path().join("input");
        let output = temp_dir.path().join("output");
        fs::create_dir_all(&input).unwrap();

        let result = cmd_generate(input, output, 1e-10, Some("bad".to_string()));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Invalid planet"));
    }

    #[test]
    fn test_cmd_generate_planet_out_of_range() {
        let temp_dir = TempDir::new().unwrap();
        let input = temp_dir.path().join("input");
        let output = temp_dir.path().join("output");
        fs::create_dir_all(&input).unwrap();

        let result = cmd_generate(input, output, 1e-10, Some("15".to_string()));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Planet must be 1-9"));
    }

    fn create_mock_vsop_file(path: &std::path::Path, planet: u8) {
        let content = format!(" VSOP2013  {}  1  0  2    PLANET VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
    2   1  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.2000000000000000 +00  0.0000000000000000 +00
", planet);
        std::fs::write(path, content).unwrap();
    }

    #[test]
    fn test_cmd_analyze_with_mock_file() {
        let temp_dir = TempDir::new().unwrap();
        let input = temp_dir.path().to_path_buf();

        // Create a mock VSOP file for planet 3
        let file_path = input.join("VSOP2013p3.dat");
        create_mock_vsop_file(&file_path, 3);

        // Analyze specific planet
        let result = cmd_analyze(input.clone(), Some("3".to_string()), 1e-10);
        assert!(result.is_ok(), "cmd_analyze failed: {:?}", result.err());
    }

    #[test]
    fn test_cmd_analyze_all_planets_mock() {
        let temp_dir = TempDir::new().unwrap();
        let input = temp_dir.path().to_path_buf();

        // Create mock files for multiple planets
        for planet in [1, 3, 5] {
            let file_path = input.join(format!("VSOP2013p{}.dat", planet));
            create_mock_vsop_file(&file_path, planet);
        }

        // Analyze all available planets (None triggers the find_planet_files path)
        let result = cmd_analyze(input, None, 1e-10);
        assert!(result.is_ok(), "cmd_analyze failed: {:?}", result.err());
    }

    #[test]
    fn test_cmd_generate_with_mock_file() {
        let temp_dir = TempDir::new().unwrap();
        let input = temp_dir.path().join("input");
        let output = temp_dir.path().join("output");
        fs::create_dir_all(&input).unwrap();

        // Create a mock VSOP file
        let file_path = input.join("VSOP2013p5.dat");
        create_mock_vsop_file(&file_path, 5);

        // Generate for specific planet
        let result = cmd_generate(input.clone(), output.clone(), 1e-10, Some("5".to_string()));
        assert!(result.is_ok(), "cmd_generate failed: {:?}", result.err());

        // Check output files
        assert!(output.join("mod.rs").exists());
        assert!(output.join("jupiter.rs").exists());
    }

    #[test]
    fn test_cmd_generate_all_planets_mock() {
        let temp_dir = TempDir::new().unwrap();
        let input = temp_dir.path().join("input");
        let output = temp_dir.path().join("output");
        fs::create_dir_all(&input).unwrap();

        // Create mock files for multiple planets
        for planet in [1, 2] {
            let file_path = input.join(format!("VSOP2013p{}.dat", planet));
            create_mock_vsop_file(&file_path, planet);
        }

        // Generate for all available planets
        let result = cmd_generate(input, output.clone(), 1e-10, None);
        assert!(result.is_ok(), "cmd_generate failed: {:?}", result.err());

        // Check output files
        assert!(output.join("mod.rs").exists());
        assert!(output.join("mercury.rs").exists());
        assert!(output.join("venus.rs").exists());
    }
}
