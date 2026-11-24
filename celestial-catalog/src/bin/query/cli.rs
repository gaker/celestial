//! CLI argument definitions for query-catalog

use celestial_core::angle::AngleUnits;
use celestial_time::julian::JulianDate;
use celestial_time::scales::conversions::{ToTAI, ToTT};
use celestial_time::scales::utc::UTC;
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Clone, Copy, ValueEnum)]
pub(crate) enum OutputFormat {
    Table,
    Json,
    Csv,
}

#[derive(Parser)]
#[command(name = "query-catalog")]
#[command(about = "Query HEALPix-indexed star catalogs")]
pub(crate) struct Cli {
    /// Path to the catalog file
    #[arg(long)]
    pub(crate) catalog: PathBuf,

    #[command(subcommand)]
    pub(crate) command: Commands,
}

#[derive(Subcommand)]
pub(crate) enum Commands {
    /// Print catalog information
    Info,
    /// Perform a cone search
    Search(SearchArgs),
}

#[derive(Args)]
pub(crate) struct SearchArgs {
    /// Right ascension (degrees, or HMS e.g. 18h36m56s, 18:36:56)
    pub(crate) ra: String,
    /// Declination (degrees, or DMS e.g. +38d47m01s, -5:22:30)
    pub(crate) dec: String,
    /// Search radius in degrees
    #[arg(long, default_value = "1.0")]
    pub(crate) radius: f64,
    /// Maximum magnitude filter
    #[arg(long)]
    pub(crate) mag_max: Option<f64>,
    /// Maximum number of results
    #[arg(long)]
    pub(crate) limit: Option<usize>,
    /// Observation epoch as a TT Julian Date (conflicts with --date)
    #[arg(long, conflicts_with = "date")]
    pub(crate) epoch: Option<f64>,
    /// Observation date as ISO 8601 YYYY-MM-DD, read as 00:00:00 UTC (conflicts with --epoch)
    #[arg(long, conflicts_with = "epoch", value_parser = parse_date)]
    pub(crate) date: Option<JulianDate>,
    /// Print query timing
    #[arg(long)]
    pub(crate) timing: bool,
    /// Output decimal degrees instead of HMS/DMS
    #[arg(long)]
    pub(crate) raw: bool,
    /// Output format
    #[arg(long, value_enum, default_value = "table")]
    pub(crate) format: OutputFormat,
}

pub(crate) fn parse_ra(s: &str) -> anyhow::Result<f64> {
    s.hms()
        .or_else(|_| s.deg())
        .map(|a| a.degrees())
        .map_err(|e| anyhow::anyhow!("Cannot parse RA '{}': {}", s, e))
}

pub(crate) fn parse_dec(s: &str) -> anyhow::Result<f64> {
    s.dms()
        .or_else(|_| s.deg())
        .map(|a| a.degrees())
        .map_err(|e| anyhow::anyhow!("Cannot parse Dec '{}': {}", s, e))
}

// Returns the TT Julian Date of 00:00:00 UTC on the date, because cone_search reads its epoch
// as TT.
fn parse_date(text: &str) -> Result<JulianDate, String> {
    // UTC::from_str trims whitespace, which a date given as one argument must not have.
    if text.chars().any(char::is_whitespace) {
        return Err("expected a date as YYYY-MM-DD with no spaces".to_string());
    }
    let utc: UTC = format!("{}T00:00:00", text)
        .parse()
        .map_err(|e| format!("expected a date as YYYY-MM-DD: {}", e))?;
    let tt = utc
        .to_tai()
        .and_then(|tai| tai.to_tt())
        .map_err(|e| e.to_string())?;
    Ok(tt.to_julian_date())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    fn search_args(extra: &[&str]) -> Result<Cli, clap::Error> {
        let base = ["query-catalog", "--catalog", "c.bin", "search", "10", "20"];
        Cli::try_parse_from(base.iter().chain(extra))
    }

    #[test]
    fn test_program_name_is_query_catalog() {
        assert_eq!(Cli::command().get_name(), "query-catalog");
    }

    #[test]
    fn test_date_without_zero_padding_is_rejected() {
        let err = parse_date("2024-3-5").unwrap_err();
        assert!(err.contains("expected a date as YYYY-MM-DD"), "{err}");
        for text in [
            " 2024-03-05",
            "2024-03-05 ",
            "\t2024-03-05",
            "2024-03-05T12:00:00",
            "2024-03",
            "2024-02-30",
            "24-03-05",
            "",
        ] {
            assert!(parse_date(text).is_err(), "{text}");
        }
    }

    #[test]
    fn test_date_is_tt_at_midnight_utc() {
        // TT - UTC is 37 s of leap seconds plus 32.184 s, 69.184 s or 0.00080074 days.
        let jd = parse_date("2024-03-05").unwrap();
        assert_eq!((jd.jd1(), jd.jd2()), (2460374.5, 0.0008007407407407408));
    }

    #[test]
    fn test_date_flag_parses_and_conflicts_with_epoch() {
        let Commands::Search(args) = search_args(&["--date", "2024-03-05"]).unwrap().command else {
            panic!("expected the search subcommand");
        };
        let jd = args.date.unwrap();
        let expected = parse_date("2024-03-05").unwrap();
        assert_eq!((jd.jd1(), jd.jd2()), (expected.jd1(), expected.jd2()));
        assert!(search_args(&["--date", "2024-3-5"]).is_err());
        assert!(search_args(&["--date", "2024-03-05", "--epoch", "2460000.5"]).is_err());
    }
}
