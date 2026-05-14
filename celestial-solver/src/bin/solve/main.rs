mod timing;

use std::path::{Path, PathBuf};
use std::fs::create_dir_all;
use anyhow::{ensure, Result};
use clap::Parser;

use celestial_images::formats::Image;
use celestial_images::ImageFormat;
use celestial_catalog::query::Catalog;
use celestial_solver::annotate::Annotation;
use celestial_solver::fit_wcs::SolveDisplay;
use celestial_solver::match_field::MatchParams;
use celestial_solver::observation::{date_obs_from_image, geodetic_from_image};
use celestial_solver::SolveParams;

#[derive(Parser)]
#[command(name = "image-solver")]
#[command(about = "Solve an astronomical image against the Celestial catalog")]
struct Cli {
    #[arg(long, short, help = "Input image (FITS or XISF). Must contain a position hint (CRVAL1/2, RA/DEC, or OBJCTRA/OBJCTDEC) and plate scale info (FOCALLEN + XPIXSZ) unless overridden by flags.")]
    image: PathBuf,

    #[arg(long, short, help = "Catalog binary (celestial-catalog .bin file).")]
    catalog: PathBuf,

    #[arg(long, short, value_name = "PATH", help = "Write the solved image to PATH with WCS (and SIP, if fit) written to the header. Extension determines format: .fits/.fit or .xisf. Omit to just print the solution.")]
    output: Option<PathBuf>,

    #[arg(long, value_name = "MM", help = "Override telescope focal length in mm. Takes priority over the FOCALLEN header.")]
    focal_length: Option<f64>,

    #[arg(long, value_name = "UM", help = "Override pixel size in microns. Takes priority over the XPIXSZ header.")]
    pixel_size: Option<f64>,

    #[arg(long, value_name = "DIR", help = "Write debug PNGs (detections, quad matches, solve overlay) to DIR. Directory is created if missing.")]
    debug: Option<PathBuf>,

    #[arg(long, value_name = "N", help = "Override the brightest-star cap used for quad matching on both image and catalog sides. Default 100. Higher values try more candidate quads at quadratic cost; lower values speed match-up at the risk of missing the field.")]
    max_stars: Option<usize>,

    #[arg(long, help = "Enable debug-level logging from the solver pipeline.")]
    verbose: bool,
}

fn main() -> Result<()> {
    let mut t = timing::Timings::new();
    let cli = Cli::parse();

    let level = if cli.verbose { log::LevelFilter::Debug } else { log::LevelFilter::Warn };
    env_logger::Builder::new()
        .filter_level(level)
        .parse_default_env()
        .format_target(false)
        .format_timestamp(None)
        .init();

    ensure!(cli.image.exists(), "{} does not exist", cli.image.display());
    ensure!(cli.catalog.exists(), "{} does not exist", cli.catalog.display());
    if let Some(output) = &cli.output {
        validate_output_path(output)?;
    }

    t.open.start();
    let img = Image::open(&cli.image)?;
    let catalog = Catalog::open(&cli.catalog)?;
    t.open.stop();

    t.solve.start();
    let mut solver = celestial_solver::solve(&img, &catalog);
    if let Some(v) = cli.focal_length {
        solver = solver.focal_length_mm(v);
    }
    if let Some(v) = cli.pixel_size {
        solver = solver.pixel_size_um(v);
    }
    if let Some(n) = cli.max_stars {
        solver = solver.params(SolveParams {
            matching: MatchParams { max_stars: n, ..MatchParams::default() },
            ..SolveParams::default()
        });
    }
    let result = solver.run()?;
    t.solve.stop();

    let software = format!("Celestial Image Solver v{}", env!("CARGO_PKG_VERSION"));
    let obs_time = date_obs_from_image(&img).map(format_date_for_display);
    let geodetic = geodetic_from_image(&img);

    let mut display = SolveDisplay::new(&result.wcs, result.sip.as_ref())
        .with_software(&software);
    if let Some(time) = &obs_time {
        display = display.with_observation_time(time);
    }
    if let Some(geo) = geodetic {
        display = display.with_geodetic(geo.lon_deg, geo.lat_deg, geo.alt_m);
    }
    eprintln!("{}", display);

    if let Some(output) = &cli.output {
        result.save_with(&img, output)?;
        eprintln!("wrote {}", output.display());
    }

    if let Some(debug_dir) = &cli.debug {
        write_debug_images(debug_dir, &img, &result)?;
    }

    println!("{}", t);
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

fn write_debug_images(
    dir: &PathBuf,
    img: &Image,
    result: &celestial_solver::solve::SolveResult,
) -> Result<()> {
    create_dir_all(dir)?;

    if let Some(mut ann) = Annotation::from_image(img) {
        ann.draw_detections(&result.stars, &result.pairs, 5);
        ann.to_image().save(dir.join("detections.png"))?;
    }

    if let Some(mut ann) = Annotation::from_image(img) {
        ann.draw_solve_overlay(&result.wcs, 20.0);
        ann.to_image().save(dir.join("solve.png"))?;
    }

    eprintln!("debug images written to {}", dir.display());
    Ok(())
}

fn format_date_for_display(date: String) -> String {
    let swapped = date.replace('T', " ");
    swapped.split('.').next().unwrap_or(&swapped).to_string()
}
