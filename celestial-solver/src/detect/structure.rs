//! Wavelet structure-map construction.

mod gaussian;
mod median;

use super::background::BackgroundMap;
use super::stats::median_in_place;
use super::Pixel;
use gaussian::separable_gaussian;
use median::median_filter_3x3;
use rayon::prelude::*;

/// Builds a binary structure map from an image and its background.
///
/// Applies a multi-scale wavelet-style highpass filter, dilates and thresholds the
/// result, then erodes to clean up noise. Output is a row-major `u8` buffer of 0/1
/// values with the same dimensions as `bg`. Input pixels pass through
/// [`Pixel::to_f64`].
///
/// `structure_layers` controls the wavelet depth. Higher values catch larger, more
/// diffuse structures at the cost of runtime and false positives.
pub(super) fn build_structure_map<T: Pixel>(
    image: &[T],
    bg: &BackgroundMap,
    structure_layers: u32,
) -> Vec<u8> {
    let w = bg.width;
    let h = bg.height;

    let t = std::time::Instant::now();
    let highpass = highpass_filter(image, bg, w, h, structure_layers);
    log::info!("[structure]     highpass: {:?}", t.elapsed());

    let t = std::time::Instant::now();
    let dilated = dilate_3x3(&highpass, w, h);
    log::info!("[structure]     dilate:   {:?}", t.elapsed());

    let t = std::time::Instant::now();
    let binary = adaptive_binarize(&dilated, w, h);
    log::info!("[structure]     binarize: {:?}", t.elapsed());

    let t = std::time::Instant::now();
    let out = erode_3x3(&binary, w, h);
    log::info!("[structure]     erode:    {:?}", t.elapsed());
    out
}

fn highpass_filter<T: Pixel>(
    image: &[T],
    bg: &BackgroundMap,
    w: usize,
    h: usize,
    structure_layers: u32,
) -> Vec<f32> {
    let t = std::time::Instant::now();
    let subtracted = subtract_background(image, bg, w * h);
    log::info!("[structure]       sub_bg:   {:?}", t.elapsed());

    let t = std::time::Instant::now();
    let subtracted = median_filter_3x3(&subtracted, w, h);
    log::info!("[structure]       median3:  {:?}", t.elapsed());

    let filter_size = 1 + (1_u32 << structure_layers);
    let sigma = filter_size as f64 / 6.0;
    let t = std::time::Instant::now();
    let smoothed = separable_gaussian(&subtracted, w, h, sigma);
    log::info!(
        "[structure]       gaussian: {:?} (sigma={sigma:.2})",
        t.elapsed()
    );

    let t = std::time::Instant::now();
    let result = normalised_excess(&subtracted, &smoothed);
    log::info!("[structure]       hp_diff:  {:?}", t.elapsed());
    result
}

// The positive part of each image pixel minus its background, over `n` pixels.
fn subtract_background<T: Pixel>(image: &[T], bg: &BackgroundMap, n: usize) -> Vec<f32> {
    let mut subtracted = vec![0.0_f32; n];
    subtracted
        .par_iter_mut()
        .zip(image.par_iter())
        .zip(bg.map.par_iter())
        .for_each(|((out, &px), &bg_val)| {
            // Non-finite pixels (blank, or beyond f32 range) carry no structure.
            let v = (px.to_f64() - bg_val) as f32;
            *out = if v.is_finite() && v > 0.0 { v } else { 0.0 };
        });
    subtracted
}

// The positive part of `subtracted` minus `smoothed`, scaled so its largest value is 1.
fn normalised_excess(subtracted: &[f32], smoothed: &[f32]) -> Vec<f32> {
    let mut result = vec![0.0_f32; subtracted.len()];
    result
        .par_iter_mut()
        .zip(subtracted.par_iter())
        .zip(smoothed.par_iter())
        .for_each(|((out, &s), &g)| {
            let v = s - g;
            *out = if v > 0.0 { v } else { 0.0 };
        });

    let max_val = result.par_iter().copied().reduce(|| 0.0_f32, f32::max);

    if max_val > 0.0 {
        let inv = 1.0 / max_val;
        result.par_iter_mut().for_each(|v| *v *= inv);
    }
    result
}

fn dilate_3x3(data: &[f32], w: usize, h: usize) -> Vec<f32> {
    let mut out = vec![0.0_f32; w * h];
    out.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let y0 = if y > 0 { y - 1 } else { 0 };
        let y1 = if y + 1 < h { y + 1 } else { h - 1 };
        for x in 0..w {
            let x0 = if x > 0 { x - 1 } else { 0 };
            let x1 = if x + 1 < w { x + 1 } else { w - 1 };
            let mut max_v = data[y * w + x];
            for ny in y0..=y1 {
                for nx in x0..=x1 {
                    let v = data[ny * w + nx];
                    if v > max_v {
                        max_v = v;
                    }
                }
            }
            row[x] = max_v;
        }
    });
    out
}

fn adaptive_binarize(data: &[f32], w: usize, h: usize) -> Vec<u8> {
    let n = w * h;
    let (median, noise) = subsample_median_and_noise(data);
    let threshold = if median < 1e-10 && noise <= 0.0 {
        median
    } else {
        median + 3.0 * noise
    };

    let mut mask = vec![0_u8; n];
    mask.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        let off = y * w;
        for x in 0..w {
            if data[off + x] > threshold {
                row[x] = 1;
            }
        }
    });
    mask
}

const SUBSAMPLE_STRIDE: usize = 21;

fn subsample_median_and_noise(data: &[f32]) -> (f32, f32) {
    let mut sample: Vec<f32> = data
        .iter()
        .step_by(SUBSAMPLE_STRIDE)
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    let Some(median) = median_in_place(&mut sample) else {
        return (0.0, 0.0);
    };
    for v in sample.iter_mut() {
        *v = (*v - median).abs();
    }
    let mad = median_in_place(&mut sample).unwrap_or(0.0);
    (median, 1.4826 * mad)
}

fn erode_3x3(mask: &[u8], w: usize, h: usize) -> Vec<u8> {
    let mut out = vec![0_u8; w * h];
    out.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        for x in 0..w {
            if mask[y * w + x] != 0 && neighbours_all_set(mask, x, y, w, h) {
                row[x] = 1;
            }
        }
    });
    out
}

// Whether every pixel of `mask` in the 3x3 neighbourhood of (x, y), clipped to the
// image, is set.
fn neighbours_all_set(mask: &[u8], x: usize, y: usize, w: usize, h: usize) -> bool {
    let y0 = if y > 0 { y - 1 } else { 0 };
    let y1 = if y + 1 < h { y + 1 } else { h - 1 };
    let x0 = if x > 0 { x - 1 } else { 0 };
    let x1 = if x + 1 < w { x + 1 } else { w - 1 };
    (y0..=y1).all(|ny| (x0..=x1).all(|nx| mask[ny * w + nx] != 0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::background::estimate_background;
    use crate::detect::test_images::{add_gaussian, gaussian_image};

    #[test]
    fn structure_map_detects_star() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, (64.0, 64.0), 5000.0, 2.0, 100.0);
        let bg = estimate_background(&image, w, h, 32).unwrap();
        let mask = build_structure_map(&image, &bg, 5);

        assert_eq!(mask.len(), w * h);
        assert_eq!(mask[64 * w + 64], 1, "star center should be detected");
    }

    #[test]
    fn structure_map_rejects_flat_field() {
        let w = 128;
        let h = 128;
        let image = vec![100.0_f32; w * h];
        let bg = estimate_background(&image, w, h, 32).unwrap();
        let mask = build_structure_map(&image, &bg, 5);

        let count: usize = mask.iter().map(|&v| v as usize).sum();
        assert_eq!(count, 0, "flat field should produce empty mask");
    }

    #[test]
    fn multiple_stars_detected() {
        let w = 256;
        let h = 256;
        let bg_level = 100.0_f32;
        let mut image = vec![bg_level; w * h];
        let sigma = 2.0_f64;

        let positions = [(64.0, 64.0), (192.0, 64.0), (64.0, 192.0), (192.0, 192.0)];
        for &centre in &positions {
            add_gaussian(&mut image, w, h, centre, 5000.0, sigma);
        }

        let bg = estimate_background(&image, w, h, 32).unwrap();
        let mask = build_structure_map(&image, &bg, 5);

        for &(sx, sy) in &positions {
            let ix = sx as usize;
            let iy = sy as usize;
            assert_eq!(
                mask[iy * w + ix],
                1,
                "star at ({sx}, {sy}) should be in structure map"
            );
        }
    }
}
