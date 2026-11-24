use super::stats::median_in_place;
use super::Pixel;
use crate::error::SolveError;
use rayon::prelude::*;

/// Per-pixel background estimate across the image.
///
/// Produced by [`estimate_background`]. The `map` is row-major, same dimensions as the
/// source image, in f64. `noise` is the sigma-clipped standard deviation of the
/// background residuals, used as an SNR denominator.
#[derive(Debug, Clone)]
pub(super) struct BackgroundMap {
    /// Background estimate at every pixel, row-major.
    pub(super) map: Vec<f64>,
    /// Sigma-clipped standard deviation of background residuals.
    pub(super) noise: f64,
    /// Image width in pixels.
    pub(super) width: usize,
    /// Image height in pixels.
    pub(super) height: usize,
}

/// Estimates per-pixel background across an image using tile medians.
///
/// Computes a `mesh_size × mesh_size` tile median grid, bilinearly interpolates it to
/// full resolution, then measures the residual noise with sigma-clipped statistics.
/// Non-finite pixels (NaN marks blank pixels) are left out of every statistic.
///
/// # Errors
///
/// [`SolveError::InvalidDimensions`] for a zero width or height,
/// [`SolveError::BufferLength`] when `image.len() != width * height`,
/// [`SolveError::InvalidParameter`] for a zero `mesh_size`,
/// [`SolveError::NoFinitePixels`] when every pixel is NaN or infinite, and
/// [`SolveError::NonFiniteStatistic`] when the background or its noise overflows.
/// A `mesh_size` larger than both image dimensions gives a single mesh cell.
pub(super) fn estimate_background<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    mesh_size: usize,
) -> Result<BackgroundMap, SolveError> {
    check_geometry(image.len(), width, height, mesh_size)?;
    let grid = MeshGrid::new(width, height, mesh_size);
    let mut mesh = build_mesh(image, width, height, &grid);
    fill_empty_cells(&mut mesh)?;
    let map = interpolate_mesh(&mesh, width, height, &grid);
    if !map.par_iter().all(|v| v.is_finite()) {
        return Err(SolveError::NonFiniteStatistic { name: "background" });
    }
    let noise = estimate_noise(image, &map);
    if !noise.is_finite() {
        return Err(SolveError::NonFiniteStatistic {
            name: "background noise",
        });
    }
    Ok(BackgroundMap {
        map,
        noise,
        width,
        height,
    })
}

fn check_geometry(
    len: usize,
    width: usize,
    height: usize,
    mesh_size: usize,
) -> Result<(), SolveError> {
    let expected = width.checked_mul(height);
    if width == 0 || height == 0 || expected.is_none() {
        return Err(SolveError::InvalidDimensions { width, height });
    }
    if let Some(expected) = expected.filter(|&e| e != len) {
        return Err(SolveError::BufferLength {
            width,
            height,
            expected,
            actual: len,
        });
    }
    if mesh_size == 0 {
        return Err(SolveError::InvalidParameter {
            name: "mesh_size",
            reason: "must be at least 1".to_string(),
        });
    }
    Ok(())
}

struct MeshGrid {
    size: usize,
    cols: usize,
    rows: usize,
}

impl MeshGrid {
    // Cells of `mesh_size` pixels, or one cell when it exceeds both image dimensions.
    fn new(width: usize, height: usize, mesh_size: usize) -> Self {
        let size = mesh_size.min(width.max(height));
        MeshGrid {
            size,
            cols: width.div_ceil(size),
            rows: height.div_ceil(size),
        }
    }
}

// Cell medians, row-major. A cell with no finite pixels is NaN.
fn build_mesh<T: Pixel>(image: &[T], width: usize, height: usize, grid: &MeshGrid) -> Vec<f64> {
    let mut mesh = vec![0.0_f64; grid.rows * grid.cols];
    let size = grid.size;

    mesh.par_chunks_mut(grid.cols)
        .enumerate()
        .for_each(|(row, mesh_row)| {
            let y0 = row * size;
            let y1 = (y0 + size).min(height);
            let mut cell_buf: Vec<f64> = Vec::with_capacity(size.min(width) * size.min(height));

            for (col, cell) in mesh_row.iter_mut().enumerate() {
                let x0 = col * size;
                let x1 = (x0 + size).min(width);

                cell_buf.clear();
                for y in y0..y1 {
                    let row_px = &image[y * width + x0..y * width + x1];
                    cell_buf.extend(row_px.iter().map(|v| v.to_f64()).filter(|v| v.is_finite()));
                }
                *cell = median_in_place(&mut cell_buf).unwrap_or(f64::NAN);
            }
        });
    mesh
}

// A cell with no finite pixels (a blank border or mask) takes the median of the other
// cells. Leaving it at 0 would make the finite pixels next to it look bright.
fn fill_empty_cells(mesh: &mut [f64]) -> Result<(), SolveError> {
    if mesh.iter().all(|v| v.is_finite()) {
        return Ok(());
    }
    let mut finite: Vec<f64> = mesh.iter().copied().filter(|v| v.is_finite()).collect();
    let fill = median_in_place(&mut finite).ok_or(SolveError::NoFinitePixels)?;
    for v in mesh.iter_mut().filter(|v| !v.is_finite()) {
        *v = fill;
    }
    Ok(())
}

fn interpolate_mesh(mesh: &[f64], width: usize, height: usize, grid: &MeshGrid) -> Vec<f64> {
    let half = grid.size as f64 * 0.5;
    let inv_mesh = 1.0 / grid.size as f64;
    let cols = grid.cols;
    let mut map = vec![0.0_f64; width * height];

    map.par_chunks_mut(width).enumerate().for_each(|(y, row)| {
        let fy = (y as f64 - half + 0.5) * inv_mesh;
        let (ry0, ry1, ty) = clamp_grid(fy, grid.rows);

        for (x, cell) in row.iter_mut().enumerate() {
            let fx = (x as f64 - half + 0.5) * inv_mesh;
            let (rx0, rx1, tx) = clamp_grid(fx, cols);

            let v00 = mesh[ry0 * cols + rx0];
            let v10 = mesh[ry0 * cols + rx1];
            let v01 = mesh[ry1 * cols + rx0];
            let v11 = mesh[ry1 * cols + rx1];

            let top = v00 + tx * (v10 - v00);
            let bot = v01 + tx * (v11 - v01);
            *cell = top + ty * (bot - top);
        }
    });
    map
}

// The two neighbouring cell indices along one axis and the interpolation weight of the
// second. An axis with one cell returns that cell twice, so the value is not
// interpolated along that axis.
fn clamp_grid(f: f64, n: usize) -> (usize, usize, f64) {
    if n < 2 {
        return (0, 0, 0.0);
    }
    let max = n - 2;
    if f <= 0.0 {
        (0, 1, 0.0)
    } else if f >= max as f64 {
        (max, max + 1, 1.0)
    } else {
        let i = (f as usize).min(max);
        (i, i + 1, f - i as f64)
    }
}

fn estimate_noise<T: Pixel>(image: &[T], background: &[f64]) -> f64 {
    let residuals: Vec<f64> = image
        .par_iter()
        .zip(background.par_iter())
        .map(|(&px, &bg)| px.to_f64() - bg)
        .filter(|r| r.is_finite())
        .collect();

    sigma_clipped_std(&residuals, 3.0, 5)
}

// Quantised data can hold most pixels at one value. Once clipping removes the pixels one
// step above it, the sigma is 0, so the last positive sigma is returned. That is the
// spread of the last pass, sqrt(f(1-f)) steps for a fraction f of pixels one step up, so
// on a nearly noiseless frame it can be far below the step itself. 0 means the first pass
// found every value equal. A non-finite sigma is returned as is for the caller to reject.
fn sigma_clipped_std(data: &[f64], kappa: f64, iterations: usize) -> f64 {
    let mut mask = vec![true; data.len()];
    let mut last_positive = 0.0;

    for _ in 0..iterations {
        let (mean, std) = masked_mean_std(data, &mask);
        if std == 0.0 {
            return last_positive;
        }
        if !std.is_finite() {
            return std;
        }
        last_positive = std;
        let changed = clip_outside(data, &mut mask, mean - kappa * std, mean + kappa * std);
        if changed == 0 {
            return std;
        }
    }

    let std = masked_mean_std(data, &mask).1;
    if std == 0.0 {
        last_positive
    } else {
        std
    }
}

// Masks out the values outside `[lo, hi]` and returns how many it masked.
fn clip_outside(data: &[f64], mask: &mut [bool], lo: f64, hi: f64) -> u64 {
    data.par_iter()
        .zip(mask.par_iter_mut())
        .map(|(&v, m)| {
            if *m && (v < lo || v > hi) {
                *m = false;
                1_u64
            } else {
                0
            }
        })
        .sum::<u64>()
}

fn masked_mean_std(data: &[f64], mask: &[bool]) -> (f64, f64) {
    let (sum, count) = ordered_masked_sum(data, mask, |v| v);
    if count < 2 {
        return (sum / count.max(1) as f64, 0.0);
    }
    let mean = sum / count as f64;
    let (var_sum, _) = ordered_masked_sum(data, mask, |v| {
        let d = v - mean;
        d * d
    });
    (mean, libm::sqrt(var_sum / (count - 1) as f64))
}

// Length of the fixed chunks that the parallel sums split the data into.
const SUM_CHUNK: usize = 1 << 14;

// Sum of `term(v)` over the masked-in values, and their count. Each chunk is summed
// sequentially and the chunk sums are added in index order, so the floating-point result
// does not depend on how rayon schedules the chunks.
fn ordered_masked_sum(data: &[f64], mask: &[bool], term: impl Fn(f64) -> f64 + Sync) -> (f64, u64) {
    let partials: Vec<(f64, u64)> = data
        .par_chunks(SUM_CHUNK)
        .zip(mask.par_chunks(SUM_CHUNK))
        .map(|(d, m)| {
            d.iter()
                .zip(m)
                .filter(|(_, &on)| on)
                .fold((0.0, 0), |(s, c), (&v, _)| (s + term(v), c + 1))
        })
        .collect();
    partials
        .iter()
        .fold((0.0, 0), |(s, c), &(ps, pc)| (s + ps, c + pc))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::test_images::{gaussian_image, quantised_frame, sample_sigma_about_10};

    #[test]
    fn flat_image_background_is_constant() {
        let width = 128;
        let height = 128;
        let image = vec![100.0_f32; width * height];
        let bg = estimate_background(&image, width, height, 32).unwrap();

        assert!(bg.map.iter().all(|&v| v == 100.0));
        assert_eq!(bg.noise, 0.0);
    }

    #[test]
    fn gradient_background_is_smooth() {
        let width = 256;
        let height = 256;
        let mut image = vec![0.0_f32; width * height];
        for y in 0..height {
            for x in 0..width {
                image[y * width + x] = (x + y) as f32;
            }
        }

        let bg = estimate_background(&image, width, height, 32).unwrap();

        // The image value x + y at (128, 128).
        assert_eq!(bg.map[128 * width + 128], 256.0);
    }

    #[test]
    fn star_does_not_dominate_background() {
        let width = 128;
        let height = 128;
        let mut image = vec![100.0_f32; width * height];
        image[64 * width + 64] = 50000.0;

        let bg = estimate_background(&image, width, height, 32).unwrap();

        assert_eq!(bg.map[64 * width + 64], 100.0);
    }

    #[test]
    fn sigma_clipped_std_rejects_outliers() {
        let mut data = vec![0.0_f64; 1000];
        for (i, v) in data.iter_mut().enumerate() {
            *v = (i as f64 * 0.001) - 0.5;
        }
        data[500] = 1000.0;
        data[501] = -1000.0;

        // The sample standard deviation of the 998 values left once the two outliers
        // are clipped, the same to the last bit.
        assert_eq!(sigma_clipped_std(&data, 3.0, 5), 0.289108975119202);
    }

    // 1M values spread over eight decades, so the float sum depends on summation order.
    fn wide_range_data() -> Vec<f64> {
        let mut s: u64 = 7;
        (0..1_000_000)
            .map(|i| {
                s = s
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                let unit = (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5;
                unit * libm::pow(10.0, (i % 8) as f64)
            })
            .collect()
    }

    fn mean_std_bits_with_threads(data: &[f64], mask: &[bool], threads: usize) -> (u64, u64) {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let (mean, std) = pool.install(|| masked_mean_std(data, mask));
        (mean.to_bits(), std.to_bits())
    }

    #[test]
    fn masked_mean_std_is_independent_of_thread_count() {
        let data = wide_range_data();
        let mask: Vec<bool> = (0..data.len()).map(|i| i % 5 != 0).collect();
        let reference = mean_std_bits_with_threads(&data, &mask, 1);
        for threads in [2, 3, 4, 7, 8, 16] {
            for _ in 0..3 {
                assert_eq!(
                    mean_std_bits_with_threads(&data, &mask, threads),
                    reference,
                    "{threads} threads"
                );
            }
        }
    }

    #[test]
    fn single_mesh_cell_uses_the_cell_value() {
        let image = vec![1.0_f32; 32 * 32];
        let bg = estimate_background(&image, 32, 32, 64).unwrap();
        assert_eq!(bg.map, vec![1.0; 32 * 32]);
        assert_eq!(bg.noise, 0.0);
    }

    // 2 x 2 cells of 2 x 2 pixels. Every cell median is 100, so the residuals are seven
    // -1, seven +1 and one 0, with mean 0 and sample standard deviation exactly 1.
    #[test]
    fn noise_skips_nan_pixels() {
        let image = [
            [99.0_f32, 101.0, 99.0, 101.0],
            [99.0, 101.0, 99.0, 101.0],
            [99.0, 101.0, 99.0, 101.0],
            [100.0, f32::NAN, 99.0, 101.0],
        ]
        .concat();
        let bg = estimate_background(&image, 4, 4, 2).unwrap();
        assert_eq!(bg.noise, 1.0);
    }

    #[test]
    fn all_nan_image_is_rejected() {
        let image = vec![f32::NAN; 128 * 128];
        let err = estimate_background(&image, 128, 128, 64).unwrap_err();
        assert!(matches!(err, SolveError::NoFinitePixels));
    }

    // A NaN cell next to cells of 100 takes their median instead of 0.
    #[test]
    fn blank_cell_takes_the_median_of_the_other_cells() {
        let mut image = vec![100.0_f32; 4 * 4];
        for y in 0..2 {
            image[y * 4] = f32::NAN;
            image[y * 4 + 1] = f32::NAN;
        }
        let bg = estimate_background(&image, 4, 4, 2).unwrap();
        assert_eq!(bg.map, vec![100.0; 4 * 4]);
    }

    #[test]
    fn invalid_geometry_is_rejected() {
        let image = vec![1.0_f32; 16];
        assert!(matches!(
            estimate_background(&image, 0, 16, 2).unwrap_err(),
            SolveError::InvalidDimensions {
                width: 0,
                height: 16
            }
        ));
        assert!(matches!(
            estimate_background(&image, 4, 5, 2).unwrap_err(),
            SolveError::BufferLength {
                width: 4,
                height: 5,
                expected: 20,
                actual: 16
            }
        ));
        assert!(matches!(
            estimate_background(&image, usize::MAX, 2, 2).unwrap_err(),
            SolveError::InvalidDimensions {
                width: usize::MAX,
                height: 2
            }
        ));
        assert!(matches!(
            estimate_background(&image, 4, 4, 0).unwrap_err(),
            SolveError::InvalidParameter {
                name: "mesh_size",
                ..
            }
        ));
    }

    #[test]
    fn huge_mesh_size_gives_one_cell() {
        let image = gaussian_image(256, 256, (128.0, 128.0), 5000.0, 2.0, 100.0);
        let one_cell = estimate_background(&image, 256, 256, 256).unwrap();
        for mesh in [1 << 32, usize::MAX] {
            let bg = estimate_background(&image, 256, 256, mesh).unwrap();
            assert_eq!(bg.map, one_cell.map, "mesh {mesh}");
            assert_eq!(bg.noise, one_cell.noise, "mesh {mesh}");
        }
    }

    // The first clipping pass removes every 11, which would leave sigma 0. The noise is
    // the sigma of that first pass, over all pixels.
    #[test]
    fn noise_of_a_quantised_frame_is_the_last_positive_sigma() {
        let image = quantised_frame(&[]);
        let bg = estimate_background(&image, 128, 128, 64).unwrap();
        assert_eq!(bg.noise, sample_sigma_about_10(&image));
    }

    #[test]
    fn overflowing_statistics_are_an_error() {
        let image: Vec<f64> = gaussian_image(128, 128, (64.0, 64.0), 5000.0, 2.0, 100.0)
            .iter()
            .map(|&v| v as f64 * 1e300)
            .collect();
        assert!(matches!(
            estimate_background(&image, 128, 128, 64).unwrap_err(),
            SolveError::NonFiniteStatistic {
                name: "background noise"
            }
        ));
    }
}
