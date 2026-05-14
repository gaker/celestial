//! Astronomical plate solver. Recovers the WCS (World Coordinate System) of an image
//! by detecting stars, matching them against a star catalog, and fitting a gnomonic
//! projection with optional SIP polynomial distortion.
//!
//! The primary entry point is [`solve()`], which returns a [`Solver`] builder:
//!
//! ```rust,ignore
//! use celestial_images::formats::Image;
//! use celestial_catalog::query::Catalog;
//!
//! let img = Image::open("m31.fits")?;
//! let catalog = Catalog::open("celestial-37M.bin")?;
//!
//! let result = celestial_solver::solve(&img, &catalog).run()?;
//!
//! println!("RMS: {:.3} px, {} stars matched", result.wcs.rms_px, result.wcs.n_stars);
//! result.save_with(&img, "m31_solved.fits")?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! For FITS or XISF images from capture software, the headers carry everything the
//! solver needs: a position hint (`CRVAL1/2` or `OBJCTRA/OBJCTDEC`), an epoch
//! (`DATE-OBS` or `MJD-OBS`), and plate scale (`FOCALLEN` + `XPIXSZ`). The builder
//! reads these automatically.
//!
//! For images without FITS headers (PNG, raw TIFF), pass an [`ImageMetadata`][metadata::ImageMetadata] directly:
//!
//! ```rust,ignore
//! use celestial_solver::metadata::ImageMetadata;
//! use celestial_coords::ICRSPosition;
//! use celestial_time::utc_from_calendar;
//!
//! # let img = celestial_images::formats::Image::open("frame.png")?;
//! # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
//! let meta = ImageMetadata {
//!     hint: ICRSPosition::from_degrees(83.633, 22.014)?,
//!     scale_arcsec: 1.5,
//!     epoch: utc_from_calendar(2026, 4, 21, 20, 0, 0.0).to_julian_date(),
//!     focal_mm: None,
//!     pixel_um: None,
//! };
//! let result = celestial_solver::solve(&img, &catalog)
//!     .with_metadata(meta)
//!     .run()?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Pixel types
//!
//! The solver accepts any pixel type the FITS/XISF reader produces — `u8`, `u16`,
//! `i16`, `i32`, `f32`, `f64`. Raw8 from a planetary camera, raw16 from a DSO session,
//! stacked f32 masters all work without pre-conversion. Saturation detection uses the
//! native ceiling of integer types (`u16::MAX` for 16-bit frames) and falls back to a
//! tied-maximum heuristic for float types.
//!
//! # Pipeline
//!
//! 1. Multi-scale wavelet structure map over [`detect`] produces centroids.
//! 2. [`match_field`] hashes four-star quads against a HEALPix-indexed catalog.
//! 3. [`fit_wcs`] fits a weighted-least-squares TAN WCS from matched pairs.
//! 4. The WCS is refined against all detected stars, with outlier rejection.
//! 5. Optional SIP polynomial distortion is fit over the refined pairs.

pub mod annotate;
pub mod detect;
pub mod fit_wcs;
pub mod match_field;
pub mod metadata;
pub mod observation;
pub mod output;
pub mod solve;

pub use solve::{solve, SolveParams, SolveResult, Solver};
