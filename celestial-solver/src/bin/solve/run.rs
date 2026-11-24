use anyhow::{ensure, Result};
use std::path::Path;

use celestial_catalog::query::catalog::Catalog;
use celestial_images::formats::format::ImageFormat;
use celestial_images::formats::image::Image;
use celestial_solver::match_field::MatchParams;
use celestial_solver::solve::{SolveParams, SolveResult};

use crate::cli::Cli;
use crate::timing::Timings;
use crate::{debug_images, report};

// The report is the result of the run and goes to stdout, so it can be redirected to a
// file on its own. Status lines, timings and log messages go to stderr.
pub(crate) fn run(cli: &Cli) -> Result<()> {
    check_paths(cli)?;
    let mut t = Timings::new();

    t.open.start();
    let img = Image::open(&cli.image)?;
    let catalog = Catalog::open(&cli.catalog)?;
    t.open.stop();

    t.solve.start();
    let result = solve(cli, &img, &catalog)?;
    t.solve.stop();

    println!("{}", report::report(&img, &result)?);

    t.save.start();
    save(cli, &img, &result)?;
    t.save.stop();

    eprintln!("{t}");
    Ok(())
}

fn check_paths(cli: &Cli) -> Result<()> {
    ensure!(cli.image.exists(), "{} does not exist", cli.image.display());
    ensure!(
        cli.catalog.exists(),
        "{} does not exist",
        cli.catalog.display()
    );
    if let Some(output) = &cli.output {
        validate_output_path(output)?;
    }
    Ok(())
}

fn validate_output_path(path: &Path) -> Result<()> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    match ImageFormat::from_extension(ext) {
        Some(ImageFormat::Fits) | Some(ImageFormat::Xisf) => Ok(()),
        _ => anyhow::bail!(
            "--output must end in .fits/.fit/.fts or .xisf (got {:?})",
            path.display().to_string()
        ),
    }
}

fn solve(cli: &Cli, img: &Image, catalog: &Catalog) -> Result<SolveResult> {
    let mut solver = celestial_solver::solve::solve(img, catalog);
    if let Some(v) = cli.focal_length {
        solver = solver.focal_length_mm(v);
    }
    if let Some(v) = cli.pixel_size {
        solver = solver.pixel_size_um(v);
    }
    if let Some(n) = cli.max_stars {
        solver = solver.params(SolveParams {
            matching: MatchParams {
                max_stars: n,
                ..MatchParams::default()
            },
            ..SolveParams::default()
        });
    }
    Ok(solver.run()?)
}

fn save(cli: &Cli, img: &Image, result: &SolveResult) -> Result<()> {
    if let Some(output) = &cli.output {
        result.save_with(img, output)?;
        eprintln!("wrote {}", output.display());
    }
    if let Some(dir) = &cli.debug {
        debug_images::write(dir, img, result)?;
    }
    Ok(())
}
