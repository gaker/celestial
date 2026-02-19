use celestial_catalog::query::catalog::FLAG_SOURCE_HIPPARCOS;
use celestial_catalog::query::{cone_search, Catalog, ConeSearchParams, ConeSearchResult};
use celestial_core::angle::{AngleUnits, DmsFmt, HmsFmt};
use celestial_core::Angle;
use celestial_time::JulianDate;
use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;
use std::time::Instant;

#[derive(Clone, ValueEnum)]
enum OutputFormat {
    Table,
    Json,
    Csv,
}

#[derive(Parser)]
#[command(name = "query-healpix")]
#[command(about = "Query HEALPix-indexed star catalogs")]
struct Cli {
    /// Path to the catalog file
    #[arg(long)]
    catalog: PathBuf,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Print catalog information
    Info,
    /// Perform a cone search
    Search {
        /// Right ascension (degrees, or HMS e.g. 18h36m56s, 18:36:56)
        ra: String,
        /// Declination (degrees, or DMS e.g. +38d47m01s, -5:22:30)
        dec: String,
        /// Search radius in degrees
        #[arg(long, default_value = "1.0")]
        radius: f64,
        /// Maximum magnitude filter
        #[arg(long)]
        mag_max: Option<f64>,
        /// Maximum number of results
        #[arg(long)]
        limit: Option<usize>,
        /// Observation epoch as Julian Date (conflicts with --date)
        #[arg(long, conflicts_with = "date")]
        epoch: Option<f64>,
        /// Observation date as ISO 8601 YYYY-MM-DD (conflicts with --epoch)
        #[arg(long, conflicts_with = "epoch")]
        date: Option<String>,
        /// Print query timing
        #[arg(long)]
        timing: bool,
        /// Output decimal degrees instead of HMS/DMS
        #[arg(long)]
        raw: bool,
        /// Output format
        #[arg(long, value_enum, default_value = "table")]
        format: OutputFormat,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Info => {
            let catalog = Catalog::open(&cli.catalog)?;
            let size_mb = catalog.file_size() as f64 / 1_048_576.0;
            println!("{}", catalog.header());
            println!(
                "File size: {} bytes ({:.2} MB)",
                catalog.file_size(),
                size_mb
            );
        }
        Commands::Search {
            ra,
            dec,
            radius,
            mag_max,
            limit,
            epoch,
            date,
            timing,
            raw,
            format,
        } => {
            let catalog = Catalog::open(&cli.catalog)?;

            let ra_deg = parse_ra(&ra)?;
            let dec_deg = parse_dec(&dec)?;

            let epoch = if let Some(date_str) = date {
                Some(parse_date_to_jd(&date_str)?)
            } else {
                epoch.map(JulianDate::from_f64)
            };

            let params = ConeSearchParams {
                ra_deg,
                dec_deg,
                radius_deg: radius,
                max_mag: mag_max,
                max_results: limit,
                epoch,
            };

            let start = if timing { Some(Instant::now()) } else { None };

            let results = cone_search(&catalog, &params);

            if let Some(start_time) = start {
                let elapsed = start_time.elapsed();
                eprintln!(
                    "Query completed in {:.2} ms",
                    elapsed.as_secs_f64() * 1000.0
                );
            }

            match format {
                OutputFormat::Table => print_table(&results, raw),
                OutputFormat::Json => print_json(&results),
                OutputFormat::Csv => print_csv(&results),
            }
        }
    }

    Ok(())
}

fn print_table(results: &[ConeSearchResult], raw: bool) {
    let hms = HmsFmt { frac_digits: 4 };
    let dms = DmsFmt { frac_digits: 4 };

    for (i, result) in results.iter().enumerate() {
        let source = source_name(result.star.flags);

        if raw {
            println!(
                "{:4}: {:>20} RA={:.6}° Dec={:+.6}° Mag={:5.2} Dist={:.4}° Source={}",
                i + 1,
                result.star.source_id,
                result.ra_deg,
                result.dec_deg,
                result.star.mag,
                result.distance_deg,
                source
            );
        } else {
            let ra_str = hms.fmt(Angle::from_degrees(result.ra_deg));
            let dec_str = dms.fmt(Angle::from_degrees(result.dec_deg));
            println!(
                "{:4}: {:>20} RA={} Dec={} Mag={:5.2} Dist={:.4}° Source={}",
                i + 1,
                result.star.source_id,
                ra_str,
                dec_str,
                result.star.mag,
                result.distance_deg,
                source
            );
        }
    }

    if results.is_empty() {
        println!("No stars found matching the search criteria.");
    } else {
        println!("\nTotal results: {}", results.len());
    }
}

#[derive(serde::Serialize)]
struct JsonStar {
    source_id: i64,
    ra_deg: f64,
    dec_deg: f64,
    mag: f32,
    distance_deg: f64,
    source: &'static str,
}

fn print_json(results: &[ConeSearchResult]) {
    let stars: Vec<JsonStar> = results
        .iter()
        .map(|r| JsonStar {
            source_id: r.star.source_id,
            ra_deg: r.ra_deg,
            dec_deg: r.dec_deg,
            mag: r.star.mag,
            distance_deg: r.distance_deg,
            source: source_name(r.star.flags),
        })
        .collect();

    println!("{}", serde_json::to_string_pretty(&stars).unwrap());
}

fn print_csv(results: &[ConeSearchResult]) {
    println!("source_id,ra_deg,dec_deg,mag,distance_deg,source");
    for r in results {
        println!(
            "{},{},{},{},{},{}",
            r.star.source_id,
            r.ra_deg,
            r.dec_deg,
            r.star.mag,
            r.distance_deg,
            source_name(r.star.flags)
        );
    }
}

fn source_name(flags: u16) -> &'static str {
    if (flags & FLAG_SOURCE_HIPPARCOS) != 0 {
        "HIPPARCOS"
    } else {
        "GAIA"
    }
}

fn parse_ra(s: &str) -> anyhow::Result<f64> {
    s.hms()
        .or_else(|_| s.deg())
        .map(|a| a.degrees())
        .map_err(|e| anyhow::anyhow!("Cannot parse RA '{}': {}", s, e))
}

fn parse_dec(s: &str) -> anyhow::Result<f64> {
    s.dms()
        .or_else(|_| s.deg())
        .map(|a| a.degrees())
        .map_err(|e| anyhow::anyhow!("Cannot parse Dec '{}': {}", s, e))
}

fn parse_date_to_jd(date_str: &str) -> anyhow::Result<JulianDate> {
    let parts: Vec<&str> = date_str.split('-').collect();
    if parts.len() != 3 {
        anyhow::bail!("Invalid date format, expected YYYY-MM-DD");
    }
    let year: i32 = parts[0].parse()?;
    let month: u8 = parts[1].parse()?;
    let day: u8 = parts[2].parse()?;

    Ok(JulianDate::from_calendar(year, month, day, 0, 0, 0.0))
}
