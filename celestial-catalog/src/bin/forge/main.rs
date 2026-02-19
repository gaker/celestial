//! Forge: astronomical catalog data pipeline CLI
//!
//! Ingests raw catalog data (Gaia DR3, Hipparcos) and produces
//! a unified HEALPix-indexed binary catalog.

mod build_index;
mod cli;
mod download_gaia;
mod gaia;
mod ingest_gaia;
mod ingest_hipparcos;
mod merge;

use clap::Parser;
use cli::{Cli, Commands};

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    if cli.verbose {
        eprintln!("Verbose mode enabled");
    }

    match &cli.command {
        Commands::DownloadGaia(args) => download_gaia::run(args, &cli),
        Commands::IngestGaia(args) => ingest_gaia::run(args, &cli),
        Commands::IngestHipparcos(args) => ingest_hipparcos::run(args, &cli),
        Commands::Merge(args) => merge::run(args, &cli),
        Commands::BuildIndex(args) => build_index::run(args, &cli),
    }
}
