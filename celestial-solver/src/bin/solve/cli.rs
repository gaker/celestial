use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "image-solver")]
#[command(about = "Solve an astronomical image against the Celestial catalog")]
pub(crate) struct Cli {
    #[arg(
        long,
        short,
        help = "Input image (FITS or XISF). Must contain a position hint (a celestial WCS, CRVAL1/2, RA/DEC, or OBJCTRA/OBJCTDEC) and plate scale info (FOCALLEN + XPIXSZ) unless --focal-length and --pixel-size override them."
    )]
    pub(crate) image: PathBuf,

    #[arg(long, short, help = "Catalog binary (celestial-catalog .bin file).")]
    pub(crate) catalog: PathBuf,

    #[arg(
        long,
        short,
        value_name = "PATH",
        help = "Write the solved image to PATH with WCS (and SIP, if fit) written to the header. Extension determines format: .fits/.fit or .xisf. Omit to just print the solution."
    )]
    pub(crate) output: Option<PathBuf>,

    #[arg(
        long,
        value_name = "MM",
        help = "Override telescope focal length in mm. Takes priority over the FOCALLEN header."
    )]
    pub(crate) focal_length: Option<f64>,

    #[arg(
        long,
        value_name = "UM",
        help = "Override pixel size in microns. Takes priority over the XPIXSZ header."
    )]
    pub(crate) pixel_size: Option<f64>,

    #[arg(
        long,
        value_name = "DIR",
        help = "Write debug PNGs (detections and solve overlay) to DIR. Directory is created if missing."
    )]
    pub(crate) debug: Option<PathBuf>,

    #[arg(
        long,
        value_name = "N",
        help = "Override the brightest-star cap used for quad matching on both image and catalog sides. Default 100. Higher values try more candidate quads at quadratic cost; lower values speed match-up at the risk of missing the field."
    )]
    pub(crate) max_stars: Option<usize>,

    #[arg(long, help = "Enable debug-level logging from the solver pipeline.")]
    pub(crate) verbose: bool,
}
