# celestial-solver

Astronomical plate solver. Given an astronomical image and a star catalog, recovers the
image's WCS solution — the mapping from pixel coordinates to sky coordinates.

[![Crates.io](https://img.shields.io/crates/v/celestial-solver)](https://crates.io/crates/celestial-solver)
[![Documentation](https://docs.rs/celestial-solver/badge.svg)](https://docs.rs/celestial-solver)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-solver)](https://github.com/gaker/celestial)

Matches star quads against a HEALPix-indexed Gaia/Hipparcos catalog, fits a TAN WCS, refines
with weighted least squares, and optionally fits SIP polynomial distortion. Accepts any pixel
type the FITS/XISF reader produces — no pre-processing required.

## Installation

```toml
[dependencies]
celestial-solver = "0.1"
celestial-images = "0.1"   # opens FITS / XISF
celestial-catalog = "0.1"  # opens the star catalog
```

## Example

For any FITS image from a capture application (NINA, SharpCap, Ekos, SGP, PixInsight, etc.),
the headers carry enough metadata that solving is one line:

```rust,ignore
use celestial_images::formats::image::Image;
use celestial_catalog::query::catalog::Catalog;

let img = Image::open("lights/m31_001.fits")?;
let catalog = Catalog::open("celestial-37M.bin")?;

let result = celestial_solver::solve::solve(&img, &catalog).run()?;

println!("RMS: {:.3} px, {} stars matched", result.wcs.rms_px, result.wcs.n_stars);
println!("Projection origin: {:.4}°, {:.4}°", result.wcs.crval1, result.wcs.crval2);
```

The solver reads the header's celestial WCS at the image centre (or `CRVAL1/2`, `RA/DEC`,
or `OBJCTRA/OBJCTDEC`) for the position hint, `MJD-OBS`/`DATE-OBS` for the epoch, and
`FOCALLEN` + `XPIXSZ` for the plate scale.

### Writing the solution back to the image

```rust,ignore
// Save a new file with WCS (and SIP, if fit) baked into the header.
// Format dispatched on extension: .fits / .fit / .xisf
result.save_with(&img, "m31_001.solved.fits")?;

// Or apply keywords to an in-memory Image you already have.
let mut annotated = img.clone();
result.write_into(&mut annotated);
```

### Overriding headers

The builder exposes one setter per piece of metadata. Missing fields fall back to headers.

```rust,ignore
// Scope has a focal reducer; FOCALLEN in the header reflects the bare OTA
let result = celestial_solver::solve::solve(&img, &catalog)
    .focal_length_mm(500)      // accepts any T: Into<f64>
    .run()?;
```

### Images without headers

For PNG, raw TIFF, or any image where the metadata lives elsewhere:

```rust,ignore
use celestial_solver::metadata::image_metadata::ImageMetadata;
use celestial_coords::frames::icrs::ICRSPosition;
use celestial_time::scales::utc::utc_from_calendar;

let meta = ImageMetadata {
    hint: ICRSPosition::from_degrees(83.633, 22.014)?,    // M1 (Crab)
    scale_arcsec: 1.5,
    epoch: utc_from_calendar(2026, 4, 21, 20, 0, 0.0)?.to_julian_date(),
    focal_mm: None,
    pixel_um: None,
};

let result = celestial_solver::solve::solve(&img, &catalog)
    .with_metadata(meta)
    .run()?;
```

## Pixel types

`solve()` accepts any `Image` regardless of its pixel type. Raw8 from a planetary cam, raw16
from a DSO session, stacked f32 masters — all handled natively, no conversion required.

Saturation is estimated from the data for every pixel type: when two or more pixels share
the image maximum and it is above the background, the limit is 90% of that maximum. Set
`DetectionParams::saturation_adu` to override it. NaN pixels are treated as blank and
skipped.

## CLI

The `image-solver` binary is behind the `cli` feature:

```sh
cargo install celestial-solver --features cli
```

It prints the solution report to stdout, and timings and log messages to stderr. An
excerpt of the report:

```text
image-solver --image m106.fit --catalog celestial-37M.bin --output m106_solved.fits

Projection ............... Gnomonic (TAN-SIP)
Projection origin ........ [1164.500000 880.500000] px -> [RA: 12 18 54.266  Dec: +47 18 38.81]
Resolution ............... 0.575 arcsec/px (linear approximation)
Rotation ................. -85.785 deg (linear approximation)
Focal distance ........... 2725.60 mm
Pixel size ............... 7.60 um
Control points ........... 24
SIP RMS .................. 0.2055 px
```

Pass `--focal-length` / `--pixel-size` to override the headers, `--debug <DIR>` to write
detection and solve-overlay PNGs for diagnostics.

## Pipeline

1. **Detection** — multi-scale wavelet structure map → connected components → thresholded
   barycenter centroiding.
2. **Quad matching** — four-star hash matching against a HEALPix-indexed catalog.
3. **WCS fit** — weighted least squares on matched centroid pairs.
4. **Refinement** — re-match all detected stars against the initial WCS, refit with outlier
   rejection.
5. **SIP** — optional polynomial distortion fit (default order 4).

## License

Licensed under either of Apache 2.0 or MIT.
