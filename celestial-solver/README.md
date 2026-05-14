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
use celestial_images::formats::Image;
use celestial_catalog::query::Catalog;

let img = Image::open("lights/m31_001.fits")?;
let catalog = Catalog::open("celestial-37M.bin")?;

let result = celestial_solver::solve(&img, &catalog).run()?;

println!("RMS: {:.3} px, {} stars matched", result.wcs.rms_px, result.wcs.n_stars);
println!("Center: {:.4}°, {:.4}°", result.wcs.crval1, result.wcs.crval2);
```

The solver reads `CRVAL1/2` (or `OBJCTRA/OBJCTDEC`, or `RA/DEC`) for the position hint,
`DATE-OBS`/`MJD-OBS` for the epoch, and `FOCALLEN` + `XPIXSZ` for the plate scale.

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
let result = celestial_solver::solve(&img, &catalog)
    .focal_length_mm(500)      // accepts any T: Into<f64>
    .run()?;
```

### Images without headers

For PNG, raw TIFF, or any image where the metadata lives elsewhere:

```rust,ignore
use celestial_solver::metadata::ImageMetadata;
use celestial_coords::ICRSPosition;
use celestial_time::utc_from_calendar;

let meta = ImageMetadata {
    hint: ICRSPosition::from_degrees(83.633, 22.014)?,    // M1 (Crab)
    scale_arcsec: 1.5,
    epoch: utc_from_calendar(2026, 4, 21, 20, 0, 0.0).to_julian_date(),
    focal_mm: None,
    pixel_um: None,
};

let result = celestial_solver::solve(&img, &catalog)
    .with_metadata(meta)
    .run()?;
```

## Pixel types

`solve()` accepts any `Image` regardless of its pixel type. Raw8 from a planetary cam, raw16
from a DSO session, stacked f32 masters — all handled natively, no conversion required.

Saturation detection uses the pixel type's native ceiling (`u8::MAX`, `u16::MAX`, …) for
integer types; for floats it falls back to a tied-maximum heuristic.

## CLI

`cargo install celestial-solver` installs the `image-solver` binary:

```text
image-solver --image m31.fits --catalog celestial-37M.bin --output m31_solved.fits

Projection origin .... RA: 00h42m44.3s  Dec: +41°16'09"
Resolution ........... 1.447 arcsec/px
Rotation ............. 0.142° (flipped)
Focal distance ....... 2032 mm
Pixel size ........... 3.76 µm
Control points ....... 84
RMS ................. 0.21 px
```

Pass `--focal-length` / `--pixel-size` to override the headers, `--debug <DIR>` to dump
detection/quad/overlay PNGs for diagnostics.

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
