#[cfg(feature = "cli")]
use clap::{Parser, Subcommand};

mod download;
mod generate;
mod parser;

use download::{default_client, download_all, find_elp_files};
use generate::{generate_moon_module, print_analysis, GenerateConfig};
use parser::parse_files;
use std::path::PathBuf;

#[cfg(feature = "cli")]
#[derive(Parser)]
#[command(name = "elpmpp02-gen")]
#[command(about = "ELP/MPP02 lunar ephemeris data processor and Rust code generator")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[cfg(feature = "cli")]
#[derive(Subcommand)]
enum Commands {
    /// Download ELP/MPP02 data files from IMCCE
    Download {
        /// Output directory for downloaded files
        #[arg(short, long, default_value = "./elpmpp02")]
        output: PathBuf,
    },
    /// Analyze ELP/MPP02 data files
    Analyze {
        /// Directory containing ELP/MPP02 data files
        #[arg(short, long)]
        input: PathBuf,
        /// Amplitude threshold for filtering terms
        #[arg(short, long, default_value = "1e-5")]
        threshold: f64,
    },
    /// Generate Rust code from ELP/MPP02 data
    Generate {
        /// Directory containing ELP/MPP02 data files
        #[arg(short, long)]
        input: PathBuf,
        /// Output directory for generated Rust code
        #[arg(short, long)]
        output: PathBuf,
        /// Amplitude threshold for filtering terms
        #[arg(short, long, default_value = "1e-5")]
        threshold: f64,
    },
}

#[cfg(feature = "cli")]
fn cmd_download(output: PathBuf) -> Result<(), String> {
    std::fs::create_dir_all(&output).map_err(|e| format!("Failed to create output dir: {}", e))?;
    let client = default_client()?;
    download_all(&client, &output)
}

#[cfg(feature = "cli")]
fn cmd_analyze(input: PathBuf, threshold: f64) -> Result<(), String> {
    let paths = find_elp_files(&input)
        .ok_or_else(|| format!("ELP/MPP02 files not found in {}", input.display()))?;

    println!("Parsing ELP/MPP02 files from {}...", input.display());
    let elp = parse_files(&paths).map_err(|e| format!("Parse error: {}", e))?;

    print_analysis(&elp, threshold);
    Ok(())
}

#[cfg(feature = "cli")]
fn cmd_generate(input: PathBuf, output: PathBuf, threshold: f64) -> Result<(), String> {
    let paths = find_elp_files(&input)
        .ok_or_else(|| format!("ELP/MPP02 files not found in {}", input.display()))?;

    println!("Parsing ELP/MPP02 files from {}...", input.display());
    let elp = parse_files(&paths).map_err(|e| format!("Parse error: {}", e))?;

    let config = GenerateConfig {
        threshold,
        output_dir: output.clone(),
    };

    println!("\nGenerating Rust code (threshold: {:.0e})...", threshold);
    generate_moon_module(&elp, &config)?;

    println!("\nGeneration complete!");
    println!("Output directory: {}", output.display());
    Ok(())
}

#[cfg(feature = "cli")]
fn main() {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Download { output } => cmd_download(output),
        Commands::Analyze { input, threshold } => cmd_analyze(input, threshold),
        Commands::Generate {
            input,
            output,
            threshold,
        } => cmd_generate(input, output, threshold),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

#[cfg(not(feature = "cli"))]
fn main() {
    eprintln!("elpmpp02-gen requires the 'cli' feature.");
    eprintln!("Run with: cargo run --features cli --bin elpmpp02-gen -- <args>");
    std::process::exit(1);
}

#[cfg(all(test, feature = "cli"))]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_cmd_analyze_files_not_found() {
        let temp_dir = TempDir::new().unwrap();
        let result = cmd_analyze(temp_dir.path().to_path_buf(), 1e-5);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    #[test]
    fn test_cmd_generate_files_not_found() {
        let temp_dir = TempDir::new().unwrap();
        let input = temp_dir.path().to_path_buf();
        let output = temp_dir.path().join("output");

        let result = cmd_generate(input, output, 1e-5);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }
}
