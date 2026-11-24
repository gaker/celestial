//! The info and search subcommands.

use crate::cli::{parse_dec, parse_ra, OutputFormat, SearchArgs};
use crate::output::{print_csv, print_json, print_table};
use celestial_catalog::query::catalog::Catalog;
use celestial_catalog::query::cone::{cone_search, ConeSearchParams};
use celestial_time::julian::JulianDate;
use std::path::Path;
use std::time::Instant;

pub(crate) fn info(catalog_path: &Path) -> anyhow::Result<()> {
    let catalog = Catalog::open(catalog_path)?;
    let size_mb = catalog.file_size() as f64 / 1_048_576.0;
    println!("{}", catalog.header());
    println!(
        "File size: {} bytes ({:.2} MB)",
        catalog.file_size(),
        size_mb
    );
    Ok(())
}

pub(crate) fn search(catalog_path: &Path, args: &SearchArgs) -> anyhow::Result<()> {
    let catalog = Catalog::open(catalog_path)?;
    let params = search_params(args)?;
    let start = Instant::now();
    let results = cone_search(&catalog, &params)?;
    if args.timing {
        let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
        eprintln!("Query completed in {:.2} ms", elapsed_ms);
    }
    match args.format {
        OutputFormat::Table => print_table(&results, args.raw),
        OutputFormat::Json => print_json(&results)?,
        OutputFormat::Csv => print_csv(&results),
    }
    Ok(())
}

fn search_params(args: &SearchArgs) -> anyhow::Result<ConeSearchParams> {
    Ok(ConeSearchParams {
        ra_deg: parse_ra(&args.ra)?,
        dec_deg: parse_dec(&args.dec)?,
        radius_deg: args.radius,
        max_mag: args.mag_max,
        max_results: args.limit,
        epoch: args.date.or(args.epoch.map(JulianDate::from_f64)),
    })
}
