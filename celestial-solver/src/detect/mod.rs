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
//! let stars = find_bright_stars(&image, 4656, 3520, &DetectionParams::default());
//! ```

mod background;
mod centroid;
pub mod components;
mod dedup;
mod pixel;
pub mod structure;

use rayon::prelude::*;

pub use background::estimate_background;
pub use dedup::deduplicate;
pub use pixel::Pixel;

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
/// let stars = find_bright_stars(&image, 128, 128, &DetectionParams::default());
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
    /// Border margin in pixels — detections closer than this to the edge are dropped.
    pub border_margin: usize,
    /// Minimum spacing (pixels) between detections after deduplication.
    pub min_spacing: f64,
    /// Explicit saturation threshold in ADU. `None` auto-estimates from pixel data for
    /// float types or uses the native type ceiling for integer types.
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

/// Per-pixel background estimate across the image.
///
/// Produced by [`estimate_background`]. The `map` is row-major, same dimensions as the
/// source image, in f64. `noise` is the sigma-clipped standard deviation of the
/// background residuals — useful as a SNR denominator.
#[derive(Debug, Clone)]
pub struct BackgroundMap {
    /// Background estimate at every pixel, row-major.
    pub map: Vec<f64>,
    /// Sigma-clipped standard deviation of background residuals.
    pub noise: f64,
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
}

/// Detects stars in an image. The main entry point to this module.
///
/// Runs background estimation, structure-map construction, connected-components
/// labeling, centroiding, SNR filtering, and deduplication. Output is sorted by
/// descending flux.
///
/// Generic over any type implementing [`Pixel`] — `u8`, `u16`, `i16`, `i32`, `f32`,
/// `f64` are all supported.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::{find_bright_stars, DetectionParams};
///
/// // u16 raw from an astro camera
/// let image: Vec<u16> = vec![0; 4656 * 3520];
/// let stars = find_bright_stars(&image, 4656, 3520, &DetectionParams::default());
/// ```
pub fn find_bright_stars<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    params: &DetectionParams,
) -> Vec<DetectedStar> {
    let t = std::time::Instant::now();
    let bg = estimate_background(image, width, height, params.mesh_size);
    log::info!("[detect]   background: {:?} (noise={:.2})", t.elapsed(), bg.noise);

    let t = std::time::Instant::now();
    let mask = structure::build_structure_map(image, &bg, params.structure_layers);
    let mask_count: usize = mask.par_iter().map(|&v| v as usize).sum();
    log::info!("[detect]   structure:  {:?} ({mask_count} px set)", t.elapsed());

    let t = std::time::Instant::now();
    let regions = components::find_connected_components(
        &mask,
        width,
        height,
        params.min_structure_size,
        params.max_structure_size,
        params.min_coverage,
    );
    log::info!("[detect]   components: {:?} ({} regions)", t.elapsed(), regions.len());

    let t = std::time::Instant::now();
    let saturation_limit = match params.saturation_adu {
        Some(adu) => adu,
        None => estimate_saturation(image),
    };

    enum CentroidOutcome {
        Star(DetectedStar),
        RejectedSnr,
        RejectedSat,
    }

    let outcomes: Vec<CentroidOutcome> = regions
        .par_iter()
        .map(|region| {
            let star = centroid::centroid_region(
                image,
                width,
                height,
                region,
                saturation_limit,
                params.min_snr,
                params.allow_clustered,
            );
            match star {
                None => CentroidOutcome::RejectedSnr,
                Some(s) if s.saturated_count > params.max_saturated_pixels => {
                    CentroidOutcome::RejectedSat
                }
                Some(s) => CentroidOutcome::Star(s),
            }
        })
        .collect();

    let mut rejected_snr = 0_usize;
    let mut rejected_sat = 0_usize;
    let mut stars: Vec<DetectedStar> = Vec::with_capacity(outcomes.len());
    for o in outcomes {
        match o {
            CentroidOutcome::Star(s) => stars.push(s),
            CentroidOutcome::RejectedSnr => rejected_snr += 1,
            CentroidOutcome::RejectedSat => rejected_sat += 1,
        }
    }
    log::info!("[detect]   centroid:   {:?} ({} passed, {} rej_snr, {} rej_sat)",
        t.elapsed(), stars.len(), rejected_snr, rejected_sat);

    let t = std::time::Instant::now();
    stars.sort_by(|a, b| {
        b.flux.partial_cmp(&a.flux).unwrap_or(std::cmp::Ordering::Equal)
    });

    deduplicate(&mut stars, params.min_spacing);
    log::info!("[detect]   sort+dedup: {:?} ({} after dedup)", t.elapsed(), stars.len());

    stars
}

fn estimate_saturation<T: Pixel>(image: &[T]) -> f64 {
    let mut max_val = 0.0_f64;
    let mut count_at_max: u32 = 0;
    for &v in image {
        let v = v.to_f64();
        if v > max_val {
            max_val = v;
            count_at_max = 1;
        } else if v == max_val {
            count_at_max += 1;
        }
    }
    if count_at_max >= 2 {
        max_val * 0.9
    } else {
        f64::INFINITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gaussian_image(
        w: usize,
        h: usize,
        cx: f64,
        cy: f64,
        peak: f32,
        sigma: f64,
        bg_level: f32,
    ) -> Vec<f32> {
        let mut image = vec![bg_level; w * h];
        for y in 0..h {
            for x in 0..w {
                let dx = x as f64 - cx;
                let dy = y as f64 - cy;
                let r_sq = dx * dx + dy * dy;
                let g = libm::exp(-r_sq / (2.0 * sigma * sigma));
                image[y * w + x] = bg_level + peak * g as f32;
            }
        }
        image
    }

    fn add_gaussian(
        image: &mut [f32],
        w: usize,
        h: usize,
        cx: f64,
        cy: f64,
        peak: f32,
        sigma: f64,
    ) {
        for y in 0..h {
            for x in 0..w {
                let dx = x as f64 - cx;
                let dy = y as f64 - cy;
                let r_sq = dx * dx + dy * dy;
                let g = libm::exp(-r_sq / (2.0 * sigma * sigma));
                image[y * w + x] += peak * g as f32;
            }
        }
    }

    #[test]
    fn estimate_saturation_unique_max_returns_infinity() {
        let image: Vec<f32> = (0..100).map(|i| i as f32).collect();
        assert_eq!(estimate_saturation(&image), f64::INFINITY);
    }

    #[test]
    fn estimate_saturation_tied_max_returns_fraction() {
        let image = vec![10.0_f32, 50.0, 30.0, 50.0, 20.0];
        let limit = estimate_saturation(&image);
        assert!((limit - 45.0).abs() < 1e-5, "expected ~45.0, got {limit}");
    }

    #[test]
    fn estimate_saturation_many_at_max() {
        let mut image = vec![0.0_f32; 1000];
        for v in image.iter_mut().take(50) {
            *v = 65535.0;
        }
        let limit = estimate_saturation(&image);
        assert!((limit - 65535.0 * 0.9).abs() < 1e-1);
    }

    #[test]
    fn estimate_saturation_empty_image() {
        let image: Vec<f32> = vec![];
        assert_eq!(estimate_saturation(&image), f64::INFINITY);
    }

    #[test]
    fn estimate_saturation_single_pixel() {
        let image = vec![1000.0_f32];
        assert_eq!(estimate_saturation(&image), f64::INFINITY);
    }

    #[test]
    fn estimate_saturation_all_zero() {
        let image = vec![0.0_f32; 100];
        let limit = estimate_saturation(&image);
        assert_eq!(limit, 0.0);
    }

    #[test]
    fn find_bright_stars_flat_field_returns_empty() {
        let w = 128;
        let h = 128;
        let image = vec![100.0_f32; w * h];
        let stars = find_bright_stars(&image, w, h, &DetectionParams::default());
        assert!(stars.is_empty(), "flat field should produce no stars, got {}", stars.len());
    }

    #[test]
    fn find_bright_stars_detects_single_star() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, 64.0, 64.0, 5000.0, 2.0, 100.0);
        let stars = find_bright_stars(&image, w, h, &DetectionParams::default());

        assert_eq!(stars.len(), 1, "expected 1 star, got {}", stars.len());
        let s = &stars[0];
        assert!((s.x - 64.0).abs() < 0.5, "cx={}, expected ~64.0", s.x);
        assert!((s.y - 64.0).abs() < 0.5, "cy={}, expected ~64.0", s.y);
        assert!(s.flux > 0.0);
        assert!(s.snr > 5.0);
    }

    #[test]
    fn find_bright_stars_sorts_by_descending_flux() {
        let w = 256;
        let h = 256;
        let mut image = vec![100.0_f32; w * h];
        add_gaussian(&mut image, w, h, 64.0, 64.0, 1000.0, 2.0);
        add_gaussian(&mut image, w, h, 192.0, 64.0, 5000.0, 2.0);
        add_gaussian(&mut image, w, h, 128.0, 192.0, 3000.0, 2.0);

        let stars = find_bright_stars(&image, w, h, &DetectionParams::default());
        assert_eq!(stars.len(), 3);

        for pair in stars.windows(2) {
            assert!(
                pair[0].flux >= pair[1].flux,
                "stars must be sorted by descending flux: {} then {}",
                pair[0].flux, pair[1].flux,
            );
        }
    }

    #[test]
    fn find_bright_stars_respects_min_snr() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, 64.0, 64.0, 5000.0, 2.0, 100.0);

        let strict = DetectionParams {
            min_snr: 1.0e9,
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, w, h, &strict);
        assert!(stars.is_empty(), "impossibly high min_snr should reject all stars");
    }

    #[test]
    fn find_bright_stars_deduplicates_close_stars() {
        let w = 128;
        let h = 128;
        let mut image = vec![100.0_f32; w * h];
        add_gaussian(&mut image, w, h, 64.0, 64.0, 5000.0, 2.0);
        add_gaussian(&mut image, w, h, 67.0, 64.0, 3000.0, 2.0);

        let params = DetectionParams {
            min_spacing: 20.0,
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, w, h, &params);
        assert_eq!(stars.len(), 1, "stars within min_spacing should be deduplicated");
    }

    #[test]
    fn find_bright_stars_detects_on_u16_image() {
        let w = 128;
        let h = 128;
        let mut image = vec![500u16; w * h];
        for y in 0..h {
            for x in 0..w {
                let dx = x as f64 - 64.0;
                let dy = y as f64 - 64.0;
                let r_sq = dx * dx + dy * dy;
                let g = libm::exp(-r_sq / 8.0);
                image[y * w + x] = image[y * w + x].saturating_add((30_000.0 * g) as u16);
            }
        }
        let stars = find_bright_stars(&image, w, h, &DetectionParams::default());
        assert_eq!(stars.len(), 1, "u16 image should detect one star");
        assert!((stars[0].x - 64.0).abs() < 0.5);
        assert!((stars[0].y - 64.0).abs() < 0.5);
    }

    #[test]
    fn find_bright_stars_detects_on_u8_image() {
        let w = 128;
        let h = 128;
        let mut image = vec![20u8; w * h];
        for y in 0..h {
            for x in 0..w {
                let dx = x as f64 - 64.0;
                let dy = y as f64 - 64.0;
                let r_sq = dx * dx + dy * dy;
                let g = libm::exp(-r_sq / 8.0);
                image[y * w + x] = image[y * w + x].saturating_add((200.0 * g) as u8);
            }
        }
        let stars = find_bright_stars(&image, w, h, &DetectionParams::default());
        assert_eq!(stars.len(), 1, "u8 image should detect one star");
    }

    #[test]
    fn find_bright_stars_respects_saturation_adu_override() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, 64.0, 64.0, 5000.0, 2.0, 100.0);

        let params = DetectionParams {
            saturation_adu: Some(500.0),
            max_saturated_pixels: 0,
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, w, h, &params);
        assert!(
            stars.is_empty(),
            "star above saturation_adu with zero tolerance should be rejected, got {} stars",
            stars.len(),
        );

        let params = DetectionParams {
            saturation_adu: Some(50000.0),
            ..DetectionParams::default()
        };
        let stars = find_bright_stars(&image, w, h, &params);
        assert_eq!(stars.len(), 1, "star below saturation_adu should survive");
        assert!(!stars[0].saturated);
    }
}
