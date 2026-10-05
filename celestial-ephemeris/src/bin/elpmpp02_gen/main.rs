use clap::{Parser, Subcommand};

mod commands;
mod download;
mod generate;
#[path = "../shared/http.rs"]
mod http;
mod parser;
#[path = "../shared/report.rs"]
mod report;

use commands::{cmd_analyze, cmd_download, cmd_generate};
use generate::DEFAULT_THRESHOLD;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "elpmpp02-gen")]
#[command(about = "ELP/MPP02 lunar ephemeris data processor and Rust code generator")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

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
        #[arg(short, long, default_value_t = DEFAULT_THRESHOLD)]
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
        #[arg(short, long, default_value_t = DEFAULT_THRESHOLD)]
        threshold: f64,
    },
}

impl Commands {
    fn run(self) -> Result<(), String> {
        match self {
            Commands::Download { output } => cmd_download(output),
            Commands::Analyze { input, threshold } => cmd_analyze(input, threshold),
            Commands::Generate {
                input,
                output,
                threshold,
            } => cmd_generate(input, output, threshold),
        }
    }
}

fn main() {
    if let Err(e) = Cli::parse().command.run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;
    use tempfile::TempDir;

    fn default_threshold(args: &[&str]) -> f64 {
        match Cli::try_parse_from(args).unwrap().command {
            Commands::Analyze { threshold, .. } | Commands::Generate { threshold, .. } => threshold,
            Commands::Download { .. } => panic!("no threshold"),
        }
    }

    #[test]
    fn default_threshold_is_the_shipped_one() {
        let analyze = default_threshold(&["elpmpp02-gen", "analyze", "-i", "in"]);
        let generate = default_threshold(&["elpmpp02-gen", "generate", "-i", "in", "-o", "out"]);
        let shipped = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/lunar_coefficients/moon.rs"
        ))
        .unwrap();
        assert_eq!(analyze, generate);
        assert!(shipped.contains(&format!("--threshold {:e}`", generate)));
    }

    fn run(args: &[&str]) -> Result<(), String> {
        let args = std::iter::once("elpmpp02-gen").chain(args.iter().copied());
        let cli = Cli::try_parse_from(args).map_err(|e| e.to_string())?;
        cli.command.run()
    }

    #[test]
    fn commands_reach_their_handlers() {
        Cli::command().debug_assert();
        let dir = TempDir::new().unwrap();
        let path = dir.path().display().to_string();
        let file = dir.path().join("file");
        std::fs::write(&file, "").unwrap();
        let download = run(&["download", "-o", file.to_str().unwrap()]).unwrap_err();
        assert!(
            download.starts_with("Failed to create output dir: "),
            "{download}"
        );
        let missing = Err(format!("ELP/MPP02 files not found in {path}"));
        assert_eq!(run(&["analyze", "-i", &path, "-t", "1e-5"]), missing);
        assert_eq!(run(&["generate", "-i", &path, "-o", &path]), missing);
    }
}
