//! query-catalog: header information and cone searches on a HEALPix-indexed star catalog.

mod cli;
mod commands;
mod output;

use clap::Parser;
use cli::{Cli, Commands};

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Info => commands::info(&cli.catalog),
        Commands::Search(args) => commands::search(&cli.catalog, args),
    }
}
