use clap::{Parser, Subcommand};

mod analyze;
mod commands;
mod download;
mod generate;
#[path = "../shared/http.rs"]
mod http;
mod parser;
#[path = "../shared/report.rs"]
mod report;

use commands::{cmd_analyze, cmd_download, cmd_generate};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "vsop2013-gen")]
#[command(about = "VSOP2013 ephemeris data processor and Rust code generator")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

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
        #[arg(short, long, help = "Amplitude cut; defaults to each shipped table's")]
        threshold: Option<f64>,
    },
    Generate {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long, help = "Amplitude cut; defaults to each shipped table's")]
        threshold: Option<f64>,
        #[arg(short, long, help = "Planet number (1-9) or 'all'")]
        planet: Option<String>,
    },
}

impl Commands {
    fn run(self) -> Result<(), String> {
        match self {
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

    fn run(args: &[&str]) -> Result<(), String> {
        let args = std::iter::once("vsop2013-gen").chain(args.iter().copied());
        let cli = Cli::try_parse_from(args).map_err(|e| e.to_string())?;
        cli.command.run()
    }

    #[test]
    fn commands_reach_their_handlers() {
        Cli::command().debug_assert();
        let dir = TempDir::new().unwrap();
        let path = dir.path().display().to_string();
        let download = run(&["download", "-o", &path, "-p", "10"]);
        assert_eq!(download, Err("Planet must be 1-9, got 10".to_string()));
        let analyze = run(&["analyze", "-i", &path, "-p", "all", "-t", "1e-7"]);
        assert_eq!(analyze, Err(format!("No VSOP2013 files found in {path}")));
        let generate = run(&["generate", "-i", &path, "-o", &path, "-p", "3"]);
        let missing = dir.path().join("VSOP2013p3.dat");
        assert_eq!(
            generate,
            Err(format!("File not found: {}", missing.display()))
        );
    }
}
