//! CLI argument definitions for forge

use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "forge")]
#[command(about = "Astronomical catalog data pipeline")]
#[command(version)]
pub struct Cli {
    /// Enable verbose output
    #[arg(short, long, global = true)]
    pub verbose: bool,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
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
pub struct DownloadGaiaArgs {
    /// Output directory for downloaded .csv.gz files
    #[arg(long)]
    pub output: PathBuf,

    /// Maximum concurrent downloads
    #[arg(long, default_value = "4")]
    pub concurrency: usize,

    /// Download only the first N files (for testing)
    #[arg(long)]
    pub limit: Option<usize>,

    /// Retry failed downloads up to N times
    #[arg(long, default_value = "3")]
    pub retries: u32,
}

#[derive(Parser)]
pub struct IngestGaiaArgs {
    /// Directory containing gzipped Gaia CSV files
    #[arg(long)]
    pub path: PathBuf,

    /// Magnitude limit (keep stars brighter than this)
    #[arg(long, default_value = "15.0")]
    pub mag_limit: f32,

    /// Output working directory for intermediate files
    #[arg(long)]
    pub output: PathBuf,

    /// Skip final concatenation (for incremental ingestion)
    #[arg(long)]
    pub no_concat: bool,

    /// Number of threads for parallel processing (0 = all cores)
    #[arg(short, long, default_value = "0")]
    pub threads: usize,
}

#[derive(Parser)]
pub struct IngestHipparcosArgs {
    /// Working directory for source data (hip2.dat, crossmatch CSV).
    /// Files are downloaded automatically if not present.
    #[arg(long)]
    pub workdir: PathBuf,

    /// Magnitude limit (keep stars brighter than this)
    #[arg(long, default_value = "7.0")]
    pub mag_limit: f32,

    /// Output working directory for ingested binary
    #[arg(long)]
    pub output: PathBuf,
}

#[derive(Parser)]
pub struct MergeArgs {
    /// Working directory containing ingested catalogs
    #[arg(long)]
    pub workdir: PathBuf,
}

#[derive(Parser)]
pub struct BuildIndexArgs {
    /// Working directory containing merged catalog
    #[arg(long)]
    pub workdir: PathBuf,

    /// HEALPix order (nside = 2^order)
    #[arg(long, default_value = "8")]
    pub healpix_order: u32,

    /// Output binary catalog file
    #[arg(long)]
    pub output: PathBuf,

    /// Maximum stars per HEALPix cell (brightest kept, rest discarded)
    #[arg(long)]
    pub max_per_cell: Option<u32>,
}
