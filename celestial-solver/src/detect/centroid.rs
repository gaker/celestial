use super::components::StarRegion;
use super::{DetectedStar, Pixel};

pub fn centroid_region<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    region: &StarRegion,
    saturation_limit: f64,
    min_snr: f64,
    allow_clustered: bool,
) -> Option<DetectedStar> {
    let (local_bg, noise) = iterative_local_background(
        image, width, height, region,
    );

    let has_saturated = region.pixels.iter()
        .any(|&(px, py)| image[py * width + px].to_f64() >= saturation_limit);

    if !allow_clustered && !has_saturated && region.pixels.len() > 20 {
        let n_maxima = count_significant_maxima(image, width, region, local_bg);
        if n_maxima > 1 {
            return None;
        }
    }

    let peak = region_peak(image, width, region);
    let snr = (peak - local_bg) / noise.max(1e-20);
    if snr < min_snr {
        return None;
    }

    let (cx, cy, weight) = thresholded_barycenter(
        image, width, region, local_bg,
    );
    if weight <= 0.0 {
        log::warn!("[centroid] REJECT zero weight npix={}", region.pixels.len());
        return None;
    }

    let (flux, saturated, saturated_count) = measure_region(
        image, width, region, local_bg, saturation_limit,
    );

    Some(DetectedStar {
        x: cx,
        y: cy,
        flux,
        snr,
        saturated,
        saturated_count,
        background: local_bg,
    })
}

fn iterative_local_background<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    region: &StarRegion,
) -> (f64, f64) {
    let mut inflate = 4_usize;
    let mut prev_median = f64::NAN;

    for _ in 0..200 {
        let samples = sample_annular_frame(
            image, width, height, region, inflate,
        );
        if samples.is_empty() {
            break;
        }
        let med = median_of(&samples);

        if prev_median.is_finite() {
            let denom = libm::fabs(prev_median).max(1e-20);
            let change = libm::fabs(prev_median - med) / denom;
            if change < 0.01 {
                let mad = mad_of(&samples, med);
                return (med, 1.4826 * mad);
            }
        }

        prev_median = med;
        inflate += 1;
    }

    if prev_median.is_finite() {
        let samples = sample_annular_frame(
            image, width, height, region, inflate.saturating_sub(1),
        );
        let mad = if samples.is_empty() {
            1.0
        } else {
            mad_of(&samples, prev_median)
        };
        return (prev_median, 1.4826 * mad);
    }

    (0.0, 1.0)
}

fn sample_annular_frame<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    region: &StarRegion,
    inflate: usize,
) -> Vec<f64> {
    let x0 = region.x_min.saturating_sub(inflate);
    let y0 = region.y_min.saturating_sub(inflate);
    let x1 = (region.x_max + inflate).min(width - 1);
    let y1 = (region.y_max + inflate).min(height - 1);

    let ix0 = region.x_min.saturating_sub(1);
    let iy0 = region.y_min.saturating_sub(1);
    let ix1 = (region.x_max + 1).min(width - 1);
    let iy1 = (region.y_max + 1).min(height - 1);

    let mut samples = Vec::new();
    for y in y0..=y1 {
        for x in x0..=x1 {
            if x >= ix0 && x <= ix1 && y >= iy0 && y <= iy1 {
                continue;
            }
            let v = image[y * width + x].to_f64();
            if v.is_finite() {
                samples.push(v);
            }
        }
    }
    samples
}

fn count_significant_maxima<T: Pixel>(
    image: &[T],
    width: usize,
    region: &StarRegion,
    local_bg: f64,
) -> usize {
    let peak = region_peak(image, width, region);
    let range = peak - local_bg;
    if range <= 0.0 {
        return 0;
    }
    let threshold = local_bg + range * 0.5;

    let mut count = 0;
    for &(px, py) in &region.pixels {
        let val = image[py * width + px].to_f64();
        if val < threshold {
            continue;
        }
        let mut is_max = true;
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
                if ux >= width {
                    continue;
                }
                let nval = image[uy * width + ux].to_f64();
                if nval > val {
                    is_max = false;
                    break;
                }
            }
            if !is_max {
                break;
            }
        }
        if is_max {
            count += 1;
        }
    }
    count
}

fn region_peak<T: Pixel>(image: &[T], width: usize, region: &StarRegion) -> f64 {
    let mut peak = f64::NEG_INFINITY;
    for &(px, py) in &region.pixels {
        let v = image[py * width + px].to_f64();
        if v > peak {
            peak = v;
        }
    }
    peak
}

fn thresholded_barycenter<T: Pixel>(
    image: &[T],
    width: usize,
    region: &StarRegion,
    local_bg: f64,
) -> (f64, f64, f64) {
    let mut values = Vec::with_capacity(region.pixels.len());
    let mut coords = Vec::with_capacity(region.pixels.len());

    for &(px, py) in &region.pixels {
        let v = image[py * width + px].to_f64() - local_bg;
        if v > 1e-7 {
            values.push(v);
            coords.push((px, py));
        }
    }

    if values.is_empty() {
        return (0.0, 0.0, 0.0);
    }

    let med = median_of(&values);
    let std = stddev_of(&values);
    let lo = med + 1.5 * std;

    let max_val = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max_val <= lo || max_val <= 0.0 {
        return flux_weighted_center(&values, &coords);
    }

    let range = max_val - lo;
    let mut sum_x = 0.0_f64;
    let mut sum_y = 0.0_f64;
    let mut sum_w = 0.0_f64;

    for (i, &v) in values.iter().enumerate() {
        if v < lo {
            continue;
        }
        let normalized = (v - lo) / range;
        let (px, py) = coords[i];
        sum_x += px as f64 * normalized;
        sum_y += py as f64 * normalized;
        sum_w += normalized;
    }

    if sum_w > 0.0 {
        (sum_x / sum_w, sum_y / sum_w, sum_w)
    } else {
        flux_weighted_center(&values, &coords)
    }
}

fn flux_weighted_center(
    values: &[f64],
    coords: &[(usize, usize)],
) -> (f64, f64, f64) {
    let mut sum_x = 0.0_f64;
    let mut sum_y = 0.0_f64;
    let mut sum_w = 0.0_f64;
    for (i, &v) in values.iter().enumerate() {
        if v > 0.0 {
            let (px, py) = coords[i];
            sum_x += px as f64 * v;
            sum_y += py as f64 * v;
            sum_w += v;
        }
    }
    if sum_w > 0.0 {
        (sum_x / sum_w, sum_y / sum_w, sum_w)
    } else {
        (0.0, 0.0, 0.0)
    }
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
        let sub = raw - local_bg;
        if sub > 0.0 {
            flux += sub;
        }
        if raw >= saturation_limit {
            saturated_count += 1;
        }
    }
    (flux, saturated_count > 0, saturated_count)
}

fn median_of(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut sorted: Vec<f64> = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) * 0.5
    }
}

fn mad_of(data: &[f64], median: f64) -> f64 {
    let deviations: Vec<f64> = data.iter().map(|&v| libm::fabs(v - median)).collect();
    median_of(&deviations)
}

fn stddev_of(data: &[f64]) -> f64 {
    let n = data.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let mean = data.iter().sum::<f64>() / n;
    let var = data.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0);
    libm::sqrt(var)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::background::estimate_background;
    use crate::detect::components::find_connected_components;
    use crate::detect::structure::build_structure_map;

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

    #[test]
    fn centroid_on_centered_gaussian() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, 64.0, 64.0, 5000.0, 2.0, 100.0);
        let bg = estimate_background(&image, w, h, 32);
        let mask = build_structure_map(&image, &bg, 5);
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5);

        assert!(!regions.is_empty(), "should find at least one region");

        let star = centroid_region(&image, w, h, &regions[0], f64::INFINITY, 5.0, false);
        assert!(star.is_some(), "should produce a star");
        let star = star.unwrap();
        assert!((star.x - 64.0).abs() < 0.5, "cx={}, expected ~64.0", star.x);
        assert!((star.y - 64.0).abs() < 0.5, "cy={}, expected ~64.0", star.y);
    }

    #[test]
    fn centroid_on_offset_gaussian() {
        let w = 128;
        let h = 128;
        let true_x = 64.3;
        let true_y = 64.7;
        let image = gaussian_image(w, h, true_x, true_y, 5000.0, 2.5, 100.0);
        let bg = estimate_background(&image, w, h, 32);
        let mask = build_structure_map(&image, &bg, 5);
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5);

        assert!(!regions.is_empty());
        let star = centroid_region(&image, w, h, &regions[0], f64::INFINITY, 5.0, false);
        assert!(star.is_some());
        let star = star.unwrap();
        assert!((star.x - true_x).abs() < 0.5, "cx={}, expected ~{true_x}", star.x);
        assert!((star.y - true_y).abs() < 0.5, "cy={}, expected ~{true_y}", star.y);
    }

    #[test]
    fn rejects_low_snr() {
        let w = 128;
        let h = 128;
        let mut image = gaussian_image(w, h, 64.0, 64.0, 50.0, 2.0, 100.0);
        let mut rng_state = 12345_u64;
        for px in &mut image {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let noise = ((rng_state >> 33) as f64 / u32::MAX as f64 - 0.5) * 20.0;
            *px += noise as f32;
        }
        let bg = estimate_background(&image, w, h, 32);
        let mask = build_structure_map(&image, &bg, 5);
        let regions = find_connected_components(&mask, w, h, 3, 2000, 0.3);

        for region in &regions {
            let star = centroid_region(&image, w, h, region, f64::INFINITY, 50.0, false);
            assert!(star.is_none(), "low SNR star should be rejected");
        }
    }

    #[test]
    fn measures_flux() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, 64.0, 64.0, 5000.0, 2.0, 100.0);
        let bg = estimate_background(&image, w, h, 32);
        let mask = build_structure_map(&image, &bg, 5);
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5);

        assert!(!regions.is_empty());
        let star = centroid_region(&image, w, h, &regions[0], f64::INFINITY, 5.0, false).unwrap();
        assert!(star.flux > 0.0, "flux should be positive");
        assert!(!star.saturated, "should not be saturated");
    }
}
