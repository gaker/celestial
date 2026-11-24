mod barycenter;
mod local_background;

use super::background::BackgroundMap;
use super::components::StarRegion;
use super::{DetectedStar, Pixel};
use barycenter::thresholded_barycenter;
use local_background::iterative_local_background;
use rayon::prelude::*;

// Per-image settings for `centroid_region`.
pub(super) struct CentroidSettings {
    pub(super) saturation_limit: f64,
    pub(super) min_snr: f64,
    pub(super) allow_clustered: bool,
    // Sigma-clipped noise of the whole image's background residuals. Positive:
    // `find_bright_stars` returns before centroiding when it is 0.
    pub(super) background_noise: f64,
}

pub(super) fn centroid_region<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    region: &StarRegion,
    settings: &CentroidSettings,
) -> Option<DetectedStar> {
    // No finite pixel around the region (it fills the image, or the annulus is blank)
    // leaves no background to measure flux or SNR against.
    let (local_bg, local_noise) =
        iterative_local_background(image, width, height, region, settings.background_noise)?;

    if !settings.allow_clustered
        && is_unresolved_blend(image, width, height, region, local_bg, settings)
    {
        return None;
    }

    // A flat annulus (clipped or quantised data) has a MAD of 0; the floor keeps its SNR
    // on the scale of the rest of the image.
    let noise = local_noise.max(settings.background_noise);
    let snr = (region_peak(image, width, region) - local_bg) / noise;
    if snr < settings.min_snr {
        return None;
    }

    let (x, y) = centroid(image, width, region, local_bg)?;
    let (flux, saturated, saturated_count) =
        measure_region(image, width, region, local_bg, settings.saturation_limit);

    Some(DetectedStar {
        x,
        y,
        flux,
        snr,
        saturated,
        saturated_count,
        background: local_bg,
    })
}

// Centroids every region in parallel; the output keeps the region order.
pub(super) fn centroid_regions<T: Pixel>(
    image: &[T],
    bg: &BackgroundMap,
    regions: &[StarRegion],
    settings: &CentroidSettings,
    max_saturated_pixels: u32,
) -> Vec<DetectedStar> {
    let t = std::time::Instant::now();
    let found: Vec<Option<DetectedStar>> = regions
        .par_iter()
        .map(|region| centroid_region(image, bg.width, bg.height, region, settings))
        .collect();
    let rejected_snr = found.iter().filter(|s| s.is_none()).count();
    let (stars, saturated): (Vec<DetectedStar>, Vec<DetectedStar>) = found
        .into_iter()
        .flatten()
        .partition(|s| s.saturated_count <= max_saturated_pixels);
    log::info!(
        "[detect]   centroid:   {:?} ({} passed, {rejected_snr} rej_snr, {} rej_sat)",
        t.elapsed(),
        stars.len(),
        saturated.len()
    );
    stars
}

// The thresholded barycenter of `region`, when its pixels above the background have
// weight and it falls on a finite pixel.
fn centroid<T: Pixel>(
    image: &[T],
    width: usize,
    region: &StarRegion,
    local_bg: f64,
) -> Option<(f64, f64)> {
    let (cx, cy, weight) = thresholded_barycenter(image, width, region, local_bg);
    if weight <= 0.0 {
        log::warn!("[centroid] REJECT zero weight npix={}", region.pixels.len());
        return None;
    }
    if !pixel_at(image, width, cx, cy).is_finite() {
        return None;
    }
    Some((cx, cy))
}

// A region of more than 20 pixels with two or more significant maxima is two stars
// that the structure map merged. Saturated regions are exempt: a flat clipped core has
// many tied maxima.
fn is_unresolved_blend<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    region: &StarRegion,
    local_bg: f64,
    settings: &CentroidSettings,
) -> bool {
    let has_saturated = region
        .pixels
        .iter()
        .any(|&(px, py)| is_saturated(image[py * width + px].to_f64(), settings.saturation_limit));
    !has_saturated
        && region.pixels.len() > 20
        && count_significant_maxima(image, width, height, region, local_bg) > 1
}

// The pixel that contains the centroid `(cx, cy)`. The centroid is a weighted mean of
// region pixels, so it rounds to a pixel inside the image.
fn pixel_at<T: Pixel>(image: &[T], width: usize, cx: f64, cy: f64) -> f64 {
    let x = libm::round(cx) as usize;
    let y = libm::round(cy) as usize;
    image[y * width + x].to_f64()
}

fn count_significant_maxima<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    region: &StarRegion,
    local_bg: f64,
) -> usize {
    let peak = region_peak(image, width, region);
    let range = peak - local_bg;
    if range <= 0.0 {
        return 0;
    }
    let threshold = local_bg + range * 0.5;
    region
        .pixels
        .iter()
        .filter(|&&(px, py)| {
            let val = image[py * width + px].to_f64();
            if !val.is_finite() || val < threshold {
                return false;
            }
            is_local_max(image, width, height, (px, py), val)
        })
        .count()
}

// Whether no pixel in the 3x3 neighbourhood of `(px, py)`, clipped to the image, is
// brighter than its value `val`.
fn is_local_max<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    (px, py): (usize, usize),
    val: f64,
) -> bool {
    for dy in -1_i32..=1 {
        for dx in -1_i32..=1 {
            if dx == 0 && dy == 0 {
                continue;
            }
            let nx = px as i32 + dx;
            let ny = py as i32 + dy;
            if nx < 0 || ny < 0 {
                continue;
            }
            let (ux, uy) = (nx as usize, ny as usize);
            if ux < width && uy < height && image[uy * width + ux].to_f64() > val {
                return false;
            }
        }
    }
    true
}

fn region_peak<T: Pixel>(image: &[T], width: usize, region: &StarRegion) -> f64 {
    let mut peak = f64::NEG_INFINITY;
    for &(px, py) in &region.pixels {
        let v = image[py * width + px].to_f64();
        if v.is_finite() && v > peak {
            peak = v;
        }
    }
    peak
}

fn measure_region<T: Pixel>(
    image: &[T],
    width: usize,
    region: &StarRegion,
    local_bg: f64,
    saturation_limit: f64,
) -> (f64, bool, u32) {
    let mut flux = 0.0_f64;
    let mut saturated_count = 0_u32;
    for &(px, py) in &region.pixels {
        let raw = image[py * width + px].to_f64();
        if !raw.is_finite() {
            continue;
        }
        let sub = raw - local_bg;
        if sub > 0.0 {
            flux += sub;
        }
        if is_saturated(raw, saturation_limit) {
            saturated_count += 1;
        }
    }
    (flux, saturated_count > 0, saturated_count)
}

fn is_saturated(v: f64, saturation_limit: f64) -> bool {
    v.is_finite() && v >= saturation_limit
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::background::estimate_background;
    use crate::detect::components::find_connected_components;
    use crate::detect::structure::build_structure_map;
    use crate::detect::test_images::gaussian_image;

    fn settings(min_snr: f64, background_noise: f64) -> CentroidSettings {
        CentroidSettings {
            saturation_limit: f64::INFINITY,
            min_snr,
            allow_clustered: false,
            background_noise,
        }
    }

    #[test]
    fn centroid_on_centered_gaussian() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, (64.0, 64.0), 5000.0, 2.0, 100.0);
        let bg = estimate_background(&image, w, h, 32).unwrap();
        let mask = build_structure_map(&image, &bg, 5);
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);

        assert!(!regions.is_empty(), "should find at least one region");

        let star = centroid_region(&image, w, h, &regions[0], &settings(5.0, bg.noise));
        assert!(star.is_some(), "should produce a star");
        let star = star.unwrap();
        // The Gaussian is centred on (64, 64); the centroid is 1.4e-14 px from it.
        assert_eq!((star.x, star.y), (63.999999999999986, 63.999999999999986));
    }

    #[test]
    fn centroid_on_offset_gaussian() {
        let w = 128;
        let h = 128;
        let true_x = 64.3;
        let true_y = 64.7;
        let image = gaussian_image(w, h, (true_x, true_y), 5000.0, 2.5, 100.0);
        let bg = estimate_background(&image, w, h, 32).unwrap();
        let mask = build_structure_map(&image, &bg, 5);
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);

        assert!(!regions.is_empty());
        let star = centroid_region(&image, w, h, &regions[0], &settings(5.0, bg.noise));
        assert!(star.is_some());
        let star = star.unwrap();
        // The centroid is 0.0063 px from the true centre (64.3, 64.7) on each axis.
        assert_eq!((star.x, star.y), (64.30633811413355, 64.69366188586643));
    }

    #[test]
    fn rejects_low_snr() {
        let w = 128;
        let h = 128;
        let mut image = gaussian_image(w, h, (64.0, 64.0), 50.0, 2.0, 100.0);
        let mut rng_state = 12345_u64;
        for px in &mut image {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let noise = ((rng_state >> 33) as f64 / u32::MAX as f64 - 0.5) * 20.0;
            *px += noise as f32;
        }
        let bg = estimate_background(&image, w, h, 32).unwrap();
        let mask = build_structure_map(&image, &bg, 5);
        let regions = find_connected_components(&mask, w, h, 3, 2000, 0.3, 5);

        for region in &regions {
            let star = centroid_region(&image, w, h, region, &settings(50.0, bg.noise));
            assert!(star.is_none(), "low SNR star should be rejected");
        }
    }

    // A 5 x 5 pyramid peaking at 5100 on a background of exactly 100, so the annulus
    // MAD is 0.
    fn flat_annulus_star() -> (Vec<f32>, StarRegion) {
        let (w, h) = (64, 64);
        let mut image = vec![100.0_f32; w * h];
        let mut pixels = Vec::new();
        for y in 30..35_usize {
            for x in 30..35_usize {
                let ring = x.abs_diff(32).max(y.abs_diff(32));
                image[y * w + x] = 5100.0 - 1000.0 * ring as f32;
                pixels.push((x, y));
            }
        }
        let region = StarRegion {
            x_min: 30,
            y_min: 30,
            x_max: 34,
            y_max: 34,
            pixels,
        };
        (image, region)
    }

    #[test]
    fn snr_noise_is_floored_at_the_background_noise() {
        let (image, region) = flat_annulus_star();
        let star = centroid_region(&image, 64, 64, &region, &settings(5.0, 10.0)).unwrap();
        // (5100 - 100) / background noise 10.
        assert_eq!(star.snr, 500.0);
    }

    #[test]
    fn measures_flux() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, (64.0, 64.0), 5000.0, 2.0, 100.0);
        let bg = estimate_background(&image, w, h, 32).unwrap();
        let mask = build_structure_map(&image, &bg, 5);
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);

        assert!(!regions.is_empty());
        let star = centroid_region(&image, w, h, &regions[0], &settings(5.0, bg.noise)).unwrap();
        assert!(star.flux > 0.0, "flux should be positive");
        assert!(!star.saturated, "should not be saturated");
    }

    // A 7 x 3 region whose peak is on the image's last row. The blend check runs on its 21
    // pixels and tests the peak's neighbours, including the row below it.
    #[test]
    fn blend_check_reads_no_row_below_the_image() {
        let (w, h) = (64, 64);
        let mut image = vec![100.0_f32; w * h];
        let mut pixels = Vec::new();
        for y in 61..64_usize {
            for x in 29..36_usize {
                let ring = x.abs_diff(32).max(y.abs_diff(63));
                image[y * w + x] = 5100.0 - 1000.0 * ring as f32;
                pixels.push((x, y));
            }
        }
        let region = StarRegion {
            x_min: 29,
            y_min: 61,
            x_max: 35,
            y_max: 63,
            pixels,
        };
        let star = centroid_region(&image, w, h, &region, &settings(5.0, 10.0)).unwrap();
        assert_eq!(star.snr, 500.0);
    }

    #[test]
    fn region_with_no_annulus_pixels_is_not_a_star() {
        let (w, h) = (3, 3);
        let mut image = vec![1000.0_f64; w * h];
        image[4] = 2000.0;
        let region = StarRegion {
            x_min: 0,
            y_min: 0,
            x_max: 2,
            y_max: 2,
            pixels: (0..h).flat_map(|y| (0..w).map(move |x| (x, y))).collect(),
        };
        let mut s = settings(0.0, 10.0);
        s.allow_clustered = true;
        let star = centroid_region(&image, w, h, &region, &s);
        assert_eq!(star.map(|s| s.background), None);
    }
}
