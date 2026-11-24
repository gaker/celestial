//! Star detection and centroiding.
//!
//! [`find_bright_stars`] runs the full detection pipeline: tiled background estimation,
//! wavelet structure map, connected components, centroiding, SNR filtering, and
//! deduplication. It's generic over [`Pixel`] — no pixel-type conversion is required
//! at the call site.
//!
//! ```rust,ignore
//! use celestial_solver::detect::{find_bright_stars, DetectionParams};
//!
//! let image: Vec<u16> = vec![/* raw16 from an astro camera */];
//! let stars = find_bright_stars(&image, 4656, 3520, &DetectionParams::default())?;
//! ```

mod background;
mod centroid;
mod components;
mod dedup;
mod params;
pub mod pixel;
mod saturation;
pub(crate) mod stats;
mod structure;
#[cfg(test)]
pub(crate) mod test_images;

use crate::error::SolveError;
use background::{estimate_background, BackgroundMap};
use centroid::{centroid_regions, CentroidSettings};
use components::StarRegion;
use dedup::deduplicate;
use params::check_params;
use pixel::Pixel;
use rayon::prelude::*;
use saturation::estimate_saturation;

/// A star detected in the image.
///
/// Position is in pixel coordinates (subpixel precision via centroiding). Flux is the
/// background-subtracted total signal in ADU. SNR is peak-over-noise from the local
/// annulus. `saturated` is `true` if any pixel in the region hit the saturation
/// threshold.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::{find_bright_stars, DetectionParams};
///
/// let image = vec![0u16; 128 * 128];
/// let stars = find_bright_stars(&image, 128, 128, &DetectionParams::default())?;
/// for s in &stars {
///     println!("({:.2}, {:.2}) flux={:.0} snr={:.1}", s.x, s.y, s.flux, s.snr);
/// }
/// ```
#[derive(Debug, Clone)]
pub struct DetectedStar {
    /// Centroid x in pixel coordinates.
    pub x: f64,
    /// Centroid y in pixel coordinates.
    pub y: f64,
    /// Background-subtracted total flux.
    pub flux: f64,
    /// Signal-to-noise ratio at the peak.
    pub snr: f64,
    /// `true` if any pixel in the region exceeded the saturation threshold.
    pub saturated: bool,
    /// Number of saturated pixels in the region.
    pub saturated_count: u32,
    /// Estimated local background level.
    pub background: f64,
}

/// Tuning knobs for [`find_bright_stars`].
///
/// Defaults work for typical deep-sky frames. Lower `min_snr` to catch fainter stars
/// at the cost of false positives; raise `min_structure_size` to reject compact
/// features (hot pixels, cosmic rays).
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::DetectionParams;
///
/// let params = DetectionParams {
///     min_snr: 8.0,
///     min_spacing: 15.0,
///     ..DetectionParams::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct DetectionParams {
    /// Number of wavelet scale layers to accumulate for the structure map.
    pub structure_layers: u32,
    /// Minimum connected-component size (in pixels) to keep.
    pub min_structure_size: usize,
    /// Maximum connected-component size (in pixels) — rejects nebulosity and glows.
    pub max_structure_size: usize,
    /// Minimum fraction of bounding-box area that must be lit inside a component.
    pub min_coverage: f64,
    /// Minimum signal-to-noise ratio to keep a detection.
    pub min_snr: f64,
    /// If `true`, don't reject components with multiple intensity maxima as unresolved
    /// pairs.
    pub allow_clustered: bool,
    /// Border margin in pixels. A region whose bounding box is closer than this to an
    /// edge is dropped, because its PSF is truncated. 0 keeps regions that touch an edge.
    pub border_margin: usize,
    /// Minimum spacing (pixels) between detections after deduplication.
    pub min_spacing: f64,
    /// Explicit saturation threshold in ADU. `None` estimates it from the pixel data for
    /// every pixel type: when two or more pixels share the image maximum and it is above
    /// the background, the threshold is 90% of that maximum; otherwise there is none.
    pub saturation_adu: Option<f64>,
    /// Maximum number of saturated pixels in a region before the whole region is
    /// rejected.
    pub max_saturated_pixels: u32,
    /// Background mesh cell size in pixels. Smaller means finer gradient following;
    /// larger means smoother.
    pub mesh_size: usize,
}

impl Default for DetectionParams {
    fn default() -> Self {
        Self {
            structure_layers: 3,
            min_structure_size: 5,
            max_structure_size: 2000,
            min_coverage: 0.5,
            min_snr: 5.0,
            allow_clustered: false,
            border_margin: 5,
            min_spacing: 20.0,
            saturation_adu: None,
            max_saturated_pixels: 500,
            mesh_size: 64,
        }
    }
}

/// Detects stars in an image. The main entry point to this module.
///
/// Runs background estimation, structure-map construction, connected-components
/// labeling, centroiding, SNR filtering, and deduplication. Output is sorted by
/// descending flux. Non-finite pixels (NaN marks blank pixels) are ignored, and no
/// star is reported with its centroid on one.
///
/// Generic over any type implementing [`Pixel`]: `u8`, `i8`, `u16`, `i16`, `u32`,
/// `i32`, `u64`, `i64`, `f32` and `f64` are all supported.
///
/// # Errors
///
/// [`SolveError::InvalidDimensions`] for a zero width or height,
/// [`SolveError::BufferLength`] when `image.len() != width * height`,
/// [`SolveError::NoFinitePixels`] when every pixel is NaN or infinite,
/// [`SolveError::NonFiniteStatistic`] when the background or its noise overflows, and
/// [`SolveError::InvalidParameter`] when a field of `params` is out of range:
/// `mesh_size` of 0, `structure_layers` above 10,
/// `min_structure_size` above `max_structure_size`, `min_coverage` outside 0 to 1,
/// `min_snr` or `min_spacing` negative or not finite, or `saturation_adu` negative or
/// NaN. `saturation_adu` of `+inf` means no pixel counts as saturated.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::{find_bright_stars, DetectionParams};
///
/// // u16 raw from an astro camera
/// let image: Vec<u16> = vec![0; 4656 * 3520];
/// let stars = find_bright_stars(&image, 4656, 3520, &DetectionParams::default())?;
/// ```
pub fn find_bright_stars<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    params: &DetectionParams,
) -> Result<Vec<DetectedStar>, SolveError> {
    check_params(params)?;
    let t = std::time::Instant::now();
    let bg = estimate_background(image, width, height, params.mesh_size)?;
    log::info!(
        "[detect]   background: {:?} (noise={:.2})",
        t.elapsed(),
        bg.noise
    );
    // A noise of 0 means every finite residual is equal (or at most one pixel is
    // finite), so no pixel stands above the background. Returning here also keeps the
    // SNR from dividing by 0.
    if bg.noise == 0.0 {
        return Ok(Vec::new());
    }

    let regions = find_regions(image, &bg, params);
    let settings = centroid_settings(image, &bg, params);
    let stars = centroid_regions(image, &bg, &regions, &settings, params.max_saturated_pixels);
    Ok(brightest_first(stars, params.min_spacing))
}

fn centroid_settings<T: Pixel>(
    image: &[T],
    bg: &BackgroundMap,
    params: &DetectionParams,
) -> CentroidSettings {
    CentroidSettings {
        saturation_limit: params
            .saturation_adu
            .unwrap_or_else(|| estimate_saturation(image, &bg.map)),
        min_snr: params.min_snr,
        allow_clustered: params.allow_clustered,
        background_noise: bg.noise,
    }
}

// `stars` by descending flux, without those within `min_spacing` of a brighter one.
fn brightest_first(mut stars: Vec<DetectedStar>, min_spacing: f64) -> Vec<DetectedStar> {
    let t = std::time::Instant::now();
    stars.sort_by(|a, b| b.flux.total_cmp(&a.flux));
    deduplicate(&mut stars, min_spacing);
    log::info!(
        "[detect]   sort+dedup: {:?} ({} after dedup)",
        t.elapsed(),
        stars.len()
    );
    stars
}

// The reference pixel (CRPIX) of every WCS fit and the origin of the quad-star tangent
// plane. Detection coordinates are 0-based: the centre of the first pixel is 0.0, so
// the image spans -0.5 to width - 0.5 and this point is half a pixel right of and
// below its geometric centre. In PixInsight's frame, where the first pixel's outer
// corner is 0.0, it is (width + 1, height + 1) / 2.
pub(crate) fn image_centre(width: usize, height: usize) -> (f64, f64) {
    (width as f64 / 2.0, height as f64 / 2.0)
}

fn find_regions<T: Pixel>(
    image: &[T],
    bg: &BackgroundMap,
    params: &DetectionParams,
) -> Vec<StarRegion> {
    let t = std::time::Instant::now();
    let mask = structure::build_structure_map(image, bg, params.structure_layers);
    let mask_count: usize = mask.par_iter().map(|&v| v as usize).sum();
    log::info!(
        "[detect]   structure:  {:?} ({mask_count} px set)",
        t.elapsed()
    );

    let t = std::time::Instant::now();
    let regions = components::find_connected_components(
        &mask,
        bg.width,
        bg.height,
        params.min_structure_size,
        params.max_structure_size,
        params.min_coverage,
        params.border_margin,
    );
    log::info!(
        "[detect]   components: {:?} ({} regions)",
        t.elapsed(),
        regions.len()
    );
    regions
}

#[cfg(test)]
mod tests {
    use super::test_images::{
        add_gaussian, gaussian_image, gaussian_profile, quantised_frame, sample_sigma_about_10,
    };
    use super::*;

    #[test]
    fn image_centre_is_half_the_size_in_detection_pixels() {
        assert_eq!(image_centre(4944, 3284), (2472.0, 1642.0));
        assert_eq!(image_centre(2327, 1761), (1163.5, 880.5));
    }

    #[test]
    fn find_bright_stars_flat_field_returns_empty() {
        let w = 128;
        let h = 128;
        let image = vec![100.0_f32; w * h];
        let stars = find_bright_stars(&image, w, h, &DetectionParams::default()).unwrap();
        assert!(
            stars.is_empty(),
            "flat field should produce no stars, got {}",
            stars.len()
        );
    }

    #[test]
    fn find_bright_stars_detects_single_star() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, (64.0, 64.0), 5000.0, 2.0, 100.0);
        let stars = find_bright_stars(&image, w, h, &DetectionParams::default()).unwrap();

        assert_eq!(stars.len(), 1, "expected 1 star, got {}", stars.len());
        let s = &stars[0];
        // The Gaussian is centred on (64, 64).
        assert_eq!(
            (s.x, s.y, s.flux, s.snr),
            (64.0, 64.0, 89720.26708984375, 6116156.423986896)
        );
    }

    #[test]
    fn find_bright_stars_sorts_by_descending_flux() {
        let w = 256;
        let h = 256;
        let mut image = vec![100.0_f32; w * h];
        add_gaussian(&mut image, w, h, (64.0, 64.0), 1000.0, 2.0);
        add_gaussian(&mut image, w, h, (192.0, 64.0), 5000.0, 2.0);
        add_gaussian(&mut image, w, h, (128.0, 192.0), 3000.0, 2.0);

        let stars = find_bright_stars(&image, w, h, &DetectionParams::default()).unwrap();
        assert_eq!(stars.len(), 3);

        for pair in stars.windows(2) {
            assert!(
                pair[0].flux >= pair[1].flux,
                "stars must be sorted by descending flux: {} then {}",
                pair[0].flux,
                pair[1].flux,
            );
        }
    }

    #[test]
    fn find_bright_stars_respects_min_snr() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, (64.0, 64.0), 5000.0, 2.0, 100.0);

        let strict = DetectionParams {
            min_snr: 1.0e9,
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, w, h, &strict).unwrap();
        assert!(
            stars.is_empty(),
            "impossibly high min_snr should reject all stars"
        );
    }

    #[test]
    fn find_bright_stars_deduplicates_close_stars() {
        let w = 128;
        let h = 128;
        let mut image = vec![100.0_f32; w * h];
        add_gaussian(&mut image, w, h, (64.0, 64.0), 5000.0, 2.0);
        add_gaussian(&mut image, w, h, (67.0, 64.0), 3000.0, 2.0);

        let params = DetectionParams {
            min_spacing: 20.0,
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, w, h, &params).unwrap();
        assert_eq!(
            stars.len(),
            1,
            "stars within min_spacing should be deduplicated"
        );
    }

    #[test]
    fn find_bright_stars_detects_on_u16_image() {
        let w = 128;
        let h = 128;
        let image: Vec<u16> = gaussian_profile(w, h, (64.0, 64.0), 2.0)
            .map(|g| 500u16.saturating_add((30_000.0 * g) as u16))
            .collect();
        let stars = find_bright_stars(&image, w, h, &DetectionParams::default()).unwrap();
        assert_eq!(stars.len(), 1, "u16 image should detect one star");
        // The Gaussian is centred on (64, 64); the centroid's y is 1.4e-14 px from it.
        assert_eq!((stars[0].x, stars[0].y), (64.0, 63.999999999999986));
    }

    #[test]
    fn find_bright_stars_detects_on_u8_image() {
        let w = 128;
        let h = 128;
        let image: Vec<u8> = gaussian_profile(w, h, (64.0, 64.0), 2.0)
            .map(|g| 20u8.saturating_add((200.0 * g) as u8))
            .collect();
        let stars = find_bright_stars(&image, w, h, &DetectionParams::default()).unwrap();
        assert_eq!(stars.len(), 1, "u8 image should detect one star");
    }

    fn one_centred_star(w: usize, h: usize) -> Vec<f32> {
        gaussian_image(w, h, (w as f64 / 2.0, h as f64 / 2.0), 5000.0, 2.0, 100.0)
    }

    #[test]
    fn find_bright_stars_on_one_mesh_cell_64x64() {
        let stars = find_bright_stars(
            &one_centred_star(64, 64),
            64,
            64,
            &DetectionParams::default(),
        )
        .unwrap();
        assert_eq!(stars.len(), 1);
    }

    #[test]
    fn find_bright_stars_on_one_mesh_column_64x200() {
        let stars = find_bright_stars(
            &one_centred_star(64, 200),
            64,
            200,
            &DetectionParams::default(),
        )
        .unwrap();
        assert_eq!(stars.len(), 1);
    }

    #[test]
    fn find_bright_stars_on_one_mesh_row_200x64() {
        let stars = find_bright_stars(
            &one_centred_star(200, 64),
            200,
            64,
            &DetectionParams::default(),
        )
        .unwrap();
        assert_eq!(stars.len(), 1);
    }

    #[test]
    fn find_bright_stars_on_1x1_image() {
        let stars = find_bright_stars(&[100.0_f32], 1, 1, &DetectionParams::default()).unwrap();
        assert_eq!(stars.len(), 0);
    }

    #[test]
    fn find_bright_stars_rejects_0x0_image() {
        let image: Vec<f32> = Vec::new();
        let err = find_bright_stars(&image, 0, 0, &DetectionParams::default()).unwrap_err();
        assert!(matches!(
            err,
            SolveError::InvalidDimensions {
                width: 0,
                height: 0
            }
        ));
    }

    #[test]
    fn find_bright_stars_rejects_zero_mesh_size() {
        let params = DetectionParams {
            mesh_size: 0,
            ..DetectionParams::default()
        };
        let err = find_bright_stars(&one_centred_star(128, 128), 128, 128, &params).unwrap_err();
        assert!(matches!(
            err,
            SolveError::InvalidParameter {
                name: "mesh_size",
                ..
            }
        ));
    }

    #[test]
    fn find_bright_stars_rejects_short_buffer() {
        let image = one_centred_star(128, 128);
        let err =
            find_bright_stars(&image[..1000], 128, 128, &DetectionParams::default()).unwrap_err();
        assert!(matches!(
            err,
            SolveError::BufferLength {
                width: 128,
                height: 128,
                expected: 16384,
                actual: 1000
            }
        ));
    }

    #[test]
    fn find_bright_stars_drops_stars_inside_border_margin() {
        let image = gaussian_image(128, 128, (6.0, 64.0), 5000.0, 2.0, 100.0);
        let stars = find_bright_stars(&image, 128, 128, &DetectionParams::default()).unwrap();
        assert_eq!(stars.len(), 0);
        let margin_1 = DetectionParams {
            border_margin: 1,
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, 128, 128, &margin_1).unwrap();
        assert_eq!(stars.len(), 1);
    }

    #[test]
    fn find_bright_stars_rejects_nan_and_oversized_params() {
        let image = one_centred_star(128, 128);
        let nan_snr = DetectionParams {
            min_snr: f64::NAN,
            ..DetectionParams::default()
        };
        let err = find_bright_stars(&image, 128, 128, &nan_snr).unwrap_err();
        assert!(matches!(
            err,
            SolveError::InvalidParameter {
                name: "min_snr",
                ..
            }
        ));
        let layers = DetectionParams {
            structure_layers: 11,
            ..DetectionParams::default()
        };
        let err = find_bright_stars(&image, 128, 128, &layers).unwrap_err();
        assert!(matches!(
            err,
            SolveError::InvalidParameter {
                name: "structure_layers",
                ..
            }
        ));
    }

    #[test]
    fn find_bright_stars_ignores_a_blank_border() {
        let mut image = one_centred_star(128, 128);
        let clean = find_bright_stars(&image, 128, 128, &DetectionParams::default()).unwrap();
        for row in image.chunks_mut(128) {
            row[..16].fill(f32::NAN);
        }
        let stars = find_bright_stars(&image, 128, 128, &DetectionParams::default()).unwrap();
        assert_eq!(stars.len(), 1);
        assert_eq!((stars[0].x, stars[0].y), (clean[0].x, clean[0].y));
    }

    #[test]
    fn find_bright_stars_never_reports_a_star_at_a_nan_pixel() {
        let mut image = one_centred_star(128, 128);
        image[64 * 128 + 64] = f32::NAN;
        let params = DetectionParams {
            allow_clustered: true,
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, 128, 128, &params).unwrap();
        assert_eq!(stars.len(), 0);
    }

    #[test]
    fn find_bright_stars_respects_saturation_adu_override() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, (64.0, 64.0), 5000.0, 2.0, 100.0);

        let params = DetectionParams {
            saturation_adu: Some(500.0),
            max_saturated_pixels: 0,
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, w, h, &params).unwrap();
        assert!(
            stars.is_empty(),
            "star above saturation_adu with zero tolerance should be rejected, got {} stars",
            stars.len(),
        );

        let params = DetectionParams {
            saturation_adu: Some(50000.0),
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, w, h, &params).unwrap();
        assert_eq!(stars.len(), 1, "star below saturation_adu should survive");
        assert!(!stars[0].saturated);
    }

    #[test]
    fn find_bright_stars_snr_is_finite_on_a_quantised_frame() {
        let image = quantised_frame(&[(40, 40, 22), (90, 80, 16)]);
        let noise = sample_sigma_about_10(&image);
        let stars = find_bright_stars(&image, 128, 128, &DetectionParams::default()).unwrap();
        let snrs: Vec<f64> = stars.iter().map(|s| s.snr).collect();
        assert_eq!(snrs, vec![(22.0 - 10.0) / noise, (16.0 - 10.0) / noise]);
    }

    #[test]
    fn find_bright_stars_rejects_negative_and_meaningless_params() {
        let image = one_centred_star(128, 128);
        let d = DetectionParams::default();
        let bad = [
            DetectionParams {
                min_snr: -1.0,
                ..d.clone()
            },
            DetectionParams {
                min_snr: f64::INFINITY,
                ..d.clone()
            },
            DetectionParams {
                min_coverage: -1.0,
                ..d.clone()
            },
            DetectionParams {
                min_coverage: 1.5,
                ..d.clone()
            },
            DetectionParams {
                min_spacing: -5.0,
                ..d.clone()
            },
            DetectionParams {
                min_spacing: f64::INFINITY,
                ..d.clone()
            },
            DetectionParams {
                saturation_adu: Some(-1.0),
                ..d.clone()
            },
            DetectionParams {
                saturation_adu: Some(f64::NEG_INFINITY),
                ..d.clone()
            },
            DetectionParams {
                min_structure_size: 10,
                max_structure_size: 9,
                ..d.clone()
            },
        ];
        for p in &bad {
            assert!(find_bright_stars(&image, 128, 128, p).is_err(), "{p:?}");
        }
        let no_limit = DetectionParams {
            saturation_adu: Some(f64::INFINITY),
            ..d
        };
        assert_eq!(
            find_bright_stars(&image, 128, 128, &no_limit)
                .unwrap()
                .len(),
            1
        );
    }
}
