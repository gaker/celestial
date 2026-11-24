//! Build HEALPix-indexed binary catalog from merged catalog.
//!
//! Memory-mapped passes for handling arbitrarily large catalogs:
//! 1. Count stars per HEALPix pixel
//! 2. Scatter stars to their final positions in output file
//! 3. Sort each pixel's stars by magnitude in-place
//! 4. With a per-pixel cap, copy the brightest stars of each pixel to the output

mod layout;
mod passes;
mod report;

use crate::cli::{BuildIndexArgs, Cli};
use crate::intermediate;
use anyhow::Context;
use celestial_catalog::query::catalog::CatalogHeader;
use celestial_catalog::query::healpix::nside_and_npix;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

// At order 10 the library's disc query answers every cone up to the whole sky, which is
// 12 * 4^10 pixels, its limit. At order 11 it rejects cones wider than 60 degrees. Order 10
// already has 12.6 million pixels: a 201 MB pixel table, and in memory 50 MB of u32 counts
// plus 100 MB of usize write positions during the scatter pass.
const MAX_BUILD_ORDER: u32 = 10;

pub(crate) fn check_build_order(order: u32) -> Result<(), String> {
    if order > MAX_BUILD_ORDER {
        return Err(format!(
            "HEALPix order {order} is above the build maximum of {MAX_BUILD_ORDER}; above it \
             the pixel table grows past 805 MB and cone searches wider than 60 degrees fail"
        ));
    }
    Ok(())
}

// The HEALPix grid of the catalog being built.
#[derive(Debug, Clone, Copy)]
struct Grid {
    order: u32,
    nside: u32,
    npix: usize,
}

impl Grid {
    fn new(order: u32) -> anyhow::Result<Self> {
        check_build_order(order).map_err(anyhow::Error::msg)?;
        let (nside, npix) = nside_and_npix(order)?;
        Ok(Self {
            order,
            nside: u32::try_from(nside)?,
            npix: usize::try_from(npix)?,
        })
    }

    fn catalog_header(&self, total_stars: u64, epoch: f64, mag_limit: f32) -> CatalogHeader {
        CatalogHeader {
            order: self.order,
            nside: self.nside,
            npix: self.npix as u64,
            total_stars,
            epoch,
            mag_limit,
        }
    }
}

pub(crate) fn run(args: &BuildIndexArgs, cli: &Cli) -> anyhow::Result<()> {
    let grid = Grid::new(args.healpix_order)?;
    create_output_dir(&args.output)?;
    let temps = TempPaths::beside(&args.output)?;
    print_plan(args, &grid, cli);
    let start = Instant::now();
    let input = passes::MergedInput::open(&args.workdir)?;
    let result = build(args, &grid, input, &temps, start);
    if result.is_err() {
        intermediate::remove_temp_file(&temps.sorted);
        intermediate::remove_temp_file(&temps.capped);
    }
    result
}

// Temporary files in the output's directory, so that the final rename does not cross file
// systems.
struct TempPaths {
    sorted: PathBuf,
    capped: PathBuf,
}

impl TempPaths {
    fn beside(output: &Path) -> anyhow::Result<Self> {
        let name = output
            .file_name()
            .with_context(|| format!("Output path {output:?} does not name a file"))?;
        let with_suffix = |suffix: &str| {
            let mut temp = name.to_os_string();
            temp.push(suffix);
            output.with_file_name(temp)
        };
        Ok(Self {
            sorted: with_suffix(".sorted.tmp"),
            capped: with_suffix(".tmp"),
        })
    }
}

// Builds the catalog in a temporary file, checks it with Catalog::open and only then renames
// it over the output, so a failed build leaves any earlier catalog in place.
fn build(
    args: &BuildIndexArgs,
    grid: &Grid,
    input: passes::MergedInput,
    temps: &TempPaths,
    start: Instant,
) -> anyhow::Result<()> {
    let merged = input.header;
    println!("Memory-mapped input file ({} stars)", merged.records);
    let header = grid.catalog_header(merged.records, merged.epoch, merged.mag_limit);
    let counts = passes::write_sorted(&input, grid, &header, &temps.sorted)?;
    drop(input);
    let built = match args.max_per_cell {
        Some(cap) => {
            println!("Pass 4: Compacting to {} max stars per cell...", cap);
            passes::compact_with_cap(&temps.sorted, &temps.capped, &header, &counts, cap)?;
            fs::remove_file(&temps.sorted)?;
            &temps.capped
        }
        None => &temps.sorted,
    };
    report::check_and_report(built, merged.records, args.max_per_cell, start)?;
    fs::rename(built, &args.output)
        .with_context(|| format!("Failed to move {built:?} to {:?}", args.output))?;
    println!("Written {:?}", args.output);
    Ok(())
}

fn create_output_dir(output: &Path) -> anyhow::Result<()> {
    match output.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create output directory {parent:?}")),
        _ => Ok(()),
    }
}

fn print_plan(args: &BuildIndexArgs, grid: &Grid, cli: &Cli) {
    println!("=== Build HEALPix Index ===");
    println!("Working directory: {:?}", args.workdir);
    println!("HEALPix order: {}", grid.order);
    println!("nside: {}", grid.nside);
    println!("npix: {}", grid.npix);
    println!("Output: {:?}", args.output);
    match args.max_per_cell {
        Some(cap) => println!("Max per cell: {}", cap),
        None => println!("Max per cell: unlimited"),
    }
    println!("Verbose: {}", cli.verbose);
    println!();
}

#[cfg(test)]
mod tests;
