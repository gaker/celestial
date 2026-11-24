//! CLI argument definitions for forge

use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "forge")]
#[command(about = "Astronomical catalog data pipeline")]
#[command(version)]
pub(crate) struct Cli {
    /// Enable verbose output
    #[arg(short, long, global = true)]
    pub(crate) verbose: bool,

    #[command(subcommand)]
    pub(crate) command: Commands,
}

#[derive(Subcommand)]
pub(crate) enum Commands {
    /// Download Gaia DR3 source files from ESA CDN
    DownloadGaia(DownloadGaiaArgs),

    /// Ingest Gaia DR3 catalog from gzipped CSV files
    IngestGaia(IngestGaiaArgs),

    /// Ingest Hipparcos catalog with epoch propagation to J2016.0
    IngestHipparcos(IngestHipparcosArgs),

    /// Merge ingested catalogs with cross-match deduplication
    Merge(MergeArgs),

    /// Build HEALPix-indexed binary catalog
    BuildIndex(BuildIndexArgs),
}

#[derive(Parser)]
pub(crate) struct DownloadGaiaArgs {
    /// Output directory for downloaded .csv.gz files
    #[arg(long)]
    pub(crate) output: PathBuf,

    /// Maximum concurrent downloads
    #[arg(long, default_value = "4")]
    pub(crate) concurrency: usize,

    /// Download only the first N files (for testing)
    #[arg(long)]
    pub(crate) limit: Option<usize>,

    /// Retry failed downloads up to N times
    #[arg(long, default_value = "3")]
    pub(crate) retries: u32,
}

#[derive(Parser)]
pub(crate) struct IngestGaiaArgs {
    /// Directory containing gzipped Gaia CSV files
    #[arg(long)]
    pub(crate) path: PathBuf,

    /// Magnitude limit (keep stars brighter than this)
    #[arg(long, default_value = "15.0")]
    pub(crate) mag_limit: f32,

    /// Output working directory for intermediate files
    #[arg(long)]
    pub(crate) output: PathBuf,

    /// Skip final concatenation (for incremental ingestion)
    #[arg(long)]
    pub(crate) no_concat: bool,

    /// Number of threads for parallel processing (0 = all cores)
    #[arg(short, long, default_value = "0")]
    pub(crate) threads: usize,
}

#[derive(Parser)]
pub(crate) struct IngestHipparcosArgs {
    /// Working directory for source data (hip2.dat, crossmatch CSV).
    /// Files are downloaded automatically if not present.
    #[arg(long)]
    pub(crate) workdir: PathBuf,

    /// Magnitude limit (keep stars brighter than this)
    #[arg(long, default_value = "7.0")]
    pub(crate) mag_limit: f32,

    /// Output working directory for ingested binary
    #[arg(long)]
    pub(crate) output: PathBuf,
}

#[derive(Parser)]
pub(crate) struct MergeArgs {
    /// Working directory containing ingested catalogs
    #[arg(long)]
    pub(crate) workdir: PathBuf,
}

#[derive(Parser)]
pub(crate) struct BuildIndexArgs {
    /// Working directory containing merged catalog
    #[arg(long)]
    pub(crate) workdir: PathBuf,

    /// HEALPix order (nside = 2^order), at most 10
    #[arg(long, default_value = "8", value_parser = parse_healpix_order)]
    pub(crate) healpix_order: u32,

    /// Output binary catalog file
    #[arg(long)]
    pub(crate) output: PathBuf,

    /// Maximum stars per HEALPix cell (brightest kept, rest discarded)
    #[arg(long)]
    pub(crate) max_per_cell: Option<u32>,
}

fn parse_healpix_order(text: &str) -> Result<u32, String> {
    let order = text.parse::<u32>().map_err(|e| e.to_string())?;
    crate::build_index::check_build_order(order)?;
    Ok(order)
}
