//! Wavelet structure-map construction.

use super::{BackgroundMap, Pixel};
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
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::{estimate_background, structure::build_structure_map};
///
/// let image = vec![0u16; 1024 * 1024];
/// let bg = estimate_background(&image, 1024, 1024, 64);
/// let mask = build_structure_map(&image, &bg, 3);
/// ```
pub fn build_structure_map<T: Pixel>(
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
    let n = w * h;

    let t = std::time::Instant::now();
    let mut subtracted = vec![0.0_f32; n];
    subtracted
        .par_iter_mut()
        .zip(image.par_iter())
        .zip(bg.map.par_iter())
        .for_each(|((out, &px), &bg_val)| {
            let v = (px.to_f64() - bg_val) as f32;
            *out = if v > 0.0 { v } else { 0.0 };
        });
    log::info!("[structure]       sub_bg:   {:?}", t.elapsed());

    let t = std::time::Instant::now();
    let subtracted = median_filter_3x3(&subtracted, w, h);
    log::info!("[structure]       median3:  {:?}", t.elapsed());

    let filter_size = 1 + (1_u32 << structure_layers);
    let sigma = filter_size as f64 / 6.0;
    let t = std::time::Instant::now();
    let smoothed = separable_gaussian(&subtracted, w, h, sigma);
    log::info!("[structure]       gaussian: {:?} (sigma={sigma:.2})", t.elapsed());

    let t = std::time::Instant::now();
    let mut result = vec![0.0_f32; n];
    result
        .par_iter_mut()
        .zip(subtracted.par_iter())
        .zip(smoothed.par_iter())
        .for_each(|((out, &s), &g)| {
            let v = s - g;
            *out = if v > 0.0 { v } else { 0.0 };
        });

    let max_val = result
        .par_iter()
        .copied()
        .reduce(|| 0.0_f32, f32::max);

    if max_val > 0.0 {
        let inv = 1.0 / max_val;
        result.par_iter_mut().for_each(|v| *v *= inv);
    }
    log::info!("[structure]       hp_diff:  {:?}", t.elapsed());

    result
}

fn separable_gaussian(data: &[f32], w: usize, h: usize, sigma: f64) -> Vec<f32> {
    let kernel_f64 = build_gaussian_kernel(sigma);
    let inv_wt = 1.0 / kernel_f64.iter().sum::<f64>();
    let kernel: Vec<f32> = kernel_f64.iter().map(|&k| (k * inv_wt) as f32).collect();
    let half = kernel.len() / 2;

    let mut temp = vec![0.0_f32; w * h];
    temp.par_chunks_mut(w)
        .enumerate()
        .for_each(|(y, row)| {
            let row_off = y * w;
            for (x, out) in row.iter_mut().enumerate() {
                let mut sum = 0.0_f32;
                for (ki, &kv) in kernel.iter().enumerate() {
                    let sx = x as i64 + ki as i64 - half as i64;
                    let sx = clamp_coord(sx, w);
                    sum += data[row_off + sx] * kv;
                }
                *out = sum;
            }
        });

    let mut out = vec![0.0_f32; w * h];
    out.par_chunks_mut(w)
        .enumerate()
        .for_each(|(y, row)| {
            for x in 0..w {
                let mut sum = 0.0_f32;
                for (ki, &kv) in kernel.iter().enumerate() {
                    let sy = y as i64 + ki as i64 - half as i64;
                    let sy = clamp_coord(sy, h);
                    sum += temp[sy * w + x] * kv;
                }
                row[x] = sum;
            }
        });

    out
}

fn build_gaussian_kernel(sigma: f64) -> Vec<f64> {
    let radius = libm::ceil(3.0 * sigma) as usize;
    let size = 2 * radius + 1;
    let mut kernel = Vec::with_capacity(size);
    let inv_2sig2 = -1.0 / (2.0 * sigma * sigma);
    for i in 0..size {
        let d = i as f64 - radius as f64;
        kernel.push(libm::exp(d * d * inv_2sig2));
    }
    kernel
}

fn clamp_coord(c: i64, max: usize) -> usize {
    if c < 0 {
        0
    } else if c >= max as i64 {
        max - 1
    } else {
        c as usize
    }
}

fn median_filter_3x3(data: &[f32], w: usize, h: usize) -> Vec<f32> {
    let mut out = vec![0.0_f32; w * h];
    if w < 3 || h < 3 {
        return median_filter_3x3_edge_only(data, w, h);
    }
    out.par_chunks_mut(w)
        .enumerate()
        .for_each(|(y, row)| {
            if y == 0 || y == h - 1 {
                fill_edge_row(data, row, y, w, h);
                return;
            }
            let off_t = (y - 1) * w;
            let off_m = y * w;
            let off_b = (y + 1) * w;
            row[0] = edge_median(data, 0, y, w, h);
            for x in 1..w - 1 {
                row[x] = median9_network([
                    data[off_t + x - 1], data[off_t + x], data[off_t + x + 1],
                    data[off_m + x - 1], data[off_m + x], data[off_m + x + 1],
                    data[off_b + x - 1], data[off_b + x], data[off_b + x + 1],
                ]);
            }
            row[w - 1] = edge_median(data, w - 1, y, w, h);
        });
    out
}

#[inline(always)]
fn median9_network(pixels: [f32; 9]) -> f32 {
    let mut buf = pixels;
    let (_, m, _) = buf.select_nth_unstable_by(4, |a, b| {
        a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
    });
    *m
}

fn fill_edge_row(data: &[f32], row: &mut [f32], y: usize, w: usize, h: usize) {
    for (x, cell) in row.iter_mut().enumerate() {
        *cell = edge_median(data, x, y, w, h);
    }
}

fn edge_median(data: &[f32], x: usize, y: usize, w: usize, h: usize) -> f32 {
    let x0 = x.saturating_sub(1);
    let x1 = if x + 1 < w { x + 1 } else { w - 1 };
    let y0 = y.saturating_sub(1);
    let y1 = if y + 1 < h { y + 1 } else { h - 1 };
    let mut buf = [0.0_f32; 9];
    let mut n = 0;
    for ny in y0..=y1 {
        for nx in x0..=x1 {
            buf[n] = data[ny * w + nx];
            n += 1;
        }
    }
    buf[..n].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    buf[n / 2]
}

fn median_filter_3x3_edge_only(data: &[f32], w: usize, h: usize) -> Vec<f32> {
    let mut out = vec![0.0_f32; w * h];
    for y in 0..h {
        for x in 0..w {
            out[y * w + x] = edge_median(data, x, y, w, h);
        }
    }
    out
}

fn dilate_3x3(data: &[f32], w: usize, h: usize) -> Vec<f32> {
    let mut out = vec![0.0_f32; w * h];
    out.par_chunks_mut(w)
        .enumerate()
        .for_each(|(y, row)| {
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
    mask.par_chunks_mut(w)
        .enumerate()
        .for_each(|(y, row)| {
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
    if sample.is_empty() {
        return (0.0, 0.0);
    }
    let median = median_in_place(&mut sample);
    for v in sample.iter_mut() {
        *v = (*v - median).abs();
    }
    let mad = median_in_place(&mut sample);
    (median, 1.4826 * mad)
}

fn erode_3x3(mask: &[u8], w: usize, h: usize) -> Vec<u8> {
    let mut out = vec![0_u8; w * h];
    out.par_chunks_mut(w)
        .enumerate()
        .for_each(|(y, row)| {
            let y0 = if y > 0 { y - 1 } else { 0 };
            let y1 = if y + 1 < h { y + 1 } else { h - 1 };
            for x in 0..w {
                if mask[y * w + x] == 0 {
                    continue;
                }
                let x0 = if x > 0 { x - 1 } else { 0 };
                let x1 = if x + 1 < w { x + 1 } else { w - 1 };
                let mut all_set = true;
                for ny in y0..=y1 {
                    for nx in x0..=x1 {
                        if mask[ny * w + nx] == 0 {
                            all_set = false;
                            break;
                        }
                    }
                    if !all_set {
                        break;
                    }
                }
                if all_set {
                    row[x] = 1;
                }
            }
        });
    out
}

fn median_in_place(buf: &mut [f32]) -> f32 {
    let n = buf.len();
    if n == 0 {
        return 0.0;
    }
    let cmp = |a: &f32, b: &f32| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal);
    let mid = n / 2;
    let (_, m, _) = buf.select_nth_unstable_by(mid, cmp);
    let hi = *m;
    if n % 2 == 1 {
        hi
    } else {
        let lo = buf[..mid]
            .iter()
            .copied()
            .fold(f32::NEG_INFINITY, f32::max);
        (lo + hi) * 0.5
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::background::estimate_background;

    #[test]
    fn median9_network_matches_sort_on_all_permutations() {
        fn permute(v: &mut [f32], k: usize, out: &mut Vec<[f32; 9]>) {
            if k == v.len() {
                let mut arr = [0.0_f32; 9];
                arr.copy_from_slice(v);
                out.push(arr);
                return;
            }
            for i in k..v.len() {
                v.swap(k, i);
                permute(v, k + 1, out);
                v.swap(k, i);
            }
        }
        let mut v = [1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let mut perms = Vec::new();
        permute(&mut v, 0, &mut perms);
        assert_eq!(perms.len(), 362880);
        for p in perms {
            let got = median9_network(p);
            assert_eq!(got, 5.0, "permutation {p:?} produced {got}");
        }
    }

    #[test]
    fn median9_network_matches_sort_on_duplicates() {
        let cases: &[([f32; 9], f32)] = &[
            ([1.0; 9], 1.0),
            ([1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0], 2.0),
            ([0.0, 0.0, 0.0, 0.0, 5.0, 9.0, 9.0, 9.0, 9.0], 5.0),
            ([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0], 1.0),
        ];
        for &(p, want) in cases {
            let got = median9_network(p);
            assert_eq!(got, want, "input {p:?} expected {want}, got {got}");
        }
    }

    #[test]
    fn median_filter_3x3_matches_naive() {
        let w = 17;
        let h = 13;
        let data: Vec<f32> = (0..w * h).map(|i| ((i * 37) % 101) as f32).collect();
        let fast = median_filter_3x3(&data, w, h);
        let naive = median_filter_3x3_edge_only(&data, w, h);
        for i in 0..w * h {
            assert_eq!(fast[i], naive[i], "pixel {i}");
        }
    }

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
    fn structure_map_detects_star() {
        let w = 128;
        let h = 128;
        let image = gaussian_image(w, h, 64.0, 64.0, 5000.0, 2.0, 100.0);
        let bg = estimate_background(&image, w, h, 32);
        let mask = build_structure_map(&image, &bg, 5);

        assert_eq!(mask.len(), w * h);
        assert_eq!(mask[64 * w + 64], 1, "star center should be detected");
    }

    #[test]
    fn structure_map_rejects_flat_field() {
        let w = 128;
        let h = 128;
        let image = vec![100.0_f32; w * h];
        let bg = estimate_background(&image, w, h, 32);
        let mask = build_structure_map(&image, &bg, 5);

        let count: usize = mask.iter().map(|&v| v as usize).sum();
        assert_eq!(count, 0, "flat field should produce empty mask");
    }

    #[test]
    fn gaussian_kernel_is_symmetric() {
        let kernel = build_gaussian_kernel(3.0);
        let n = kernel.len();
        for i in 0..n / 2 {
            assert_eq!(kernel[i], kernel[n - 1 - i]);
        }
    }

    #[test]
    fn multiple_stars_detected() {
        let w = 256;
        let h = 256;
        let bg_level = 100.0_f32;
        let mut image = vec![bg_level; w * h];
        let sigma = 2.0_f64;

        let positions = [(64.0, 64.0), (192.0, 64.0), (64.0, 192.0), (192.0, 192.0)];
        for &(sx, sy) in &positions {
            for y in 0..h {
                for x in 0..w {
                    let dx = x as f64 - sx;
                    let dy = y as f64 - sy;
                    let r_sq = dx * dx + dy * dy;
                    let g = libm::exp(-r_sq / (2.0 * sigma * sigma));
                    image[y * w + x] += 5000.0 * g as f32;
                }
            }
        }

        let bg = estimate_background(&image, w, h, 32);
        let mask = build_structure_map(&image, &bg, 5);

        for &(sx, sy) in &positions {
            let ix = sx as usize;
            let iy = sy as usize;
            assert_eq!(
                mask[iy * w + ix], 1,
                "star at ({sx}, {sy}) should be in structure map"
            );
        }
    }
}
