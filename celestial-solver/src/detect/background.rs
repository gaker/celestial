use super::{BackgroundMap, Pixel};
use rayon::prelude::*;

/// Estimates per-pixel background across an image using tile medians.
///
/// Computes a `mesh_size × mesh_size` tile median grid, bilinearly interpolates it to
/// full resolution, then measures the residual noise with sigma-clipped statistics.
///
/// Called internally by [`crate::detect::find_bright_stars`]. Exposed for callers that
/// want to re-use the background map (e.g. for photometry or flat-field analysis).
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::estimate_background;
///
/// let image = vec![0u16; 1024 * 1024];
/// let bg = estimate_background(&image, 1024, 1024, 64);
/// println!("noise sigma: {:.2}", bg.noise);
/// ```
pub fn estimate_background<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    mesh_size: usize,
) -> BackgroundMap {
    let cols = width.div_ceil(mesh_size);
    let rows = height.div_ceil(mesh_size);
    let mesh = build_mesh(image, width, height, mesh_size, cols, rows);
    let map = interpolate_mesh(&mesh, width, height, mesh_size, cols, rows);
    let noise = estimate_noise(image, &map);
    BackgroundMap { map, noise, width, height }
}

fn build_mesh<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    mesh_size: usize,
    cols: usize,
    rows: usize,
) -> Vec<f64> {
    let mut mesh = vec![0.0_f64; rows * cols];

    mesh.par_chunks_mut(cols).enumerate().for_each(|(row, mesh_row)| {
        let y0 = row * mesh_size;
        let y1 = (y0 + mesh_size).min(height);
        let mut cell_buf: Vec<f64> = Vec::with_capacity(mesh_size * mesh_size);

        for (col, cell) in mesh_row.iter_mut().enumerate() {
            let x0 = col * mesh_size;
            let x1 = (x0 + mesh_size).min(width);

            cell_buf.clear();
            for y in y0..y1 {
                for x in x0..x1 {
                    let v = image[y * width + x].to_f64();
                    if v.is_finite() {
                        cell_buf.push(v);
                    }
                }
            }

            *cell = if cell_buf.is_empty() {
                0.0
            } else {
                median(&mut cell_buf)
            };
        }
    });
    mesh
}

fn interpolate_mesh(
    mesh: &[f64],
    width: usize,
    height: usize,
    mesh_size: usize,
    cols: usize,
    rows: usize,
) -> Vec<f64> {
    let half = mesh_size as f64 * 0.5;
    let inv_mesh = 1.0 / mesh_size as f64;
    let mut map = vec![0.0_f64; width * height];

    map.par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            let fy = (y as f64 - half + 0.5) * inv_mesh;
            let (ry, ty) = clamp_grid(fy, rows);

            for (x, cell) in row.iter_mut().enumerate() {
                let fx = (x as f64 - half + 0.5) * inv_mesh;
                let (rx, tx) = clamp_grid(fx, cols);

                let v00 = mesh[ry * cols + rx];
                let v10 = mesh[ry * cols + rx + 1];
                let v01 = mesh[(ry + 1) * cols + rx];
                let v11 = mesh[(ry + 1) * cols + rx + 1];

                let top = v00 + tx * (v10 - v00);
                let bot = v01 + tx * (v11 - v01);
                *cell = top + ty * (bot - top);
            }
        });
    map
}

fn clamp_grid(f: f64, n: usize) -> (usize, f64) {
    if n < 2 {
        return (0, 0.0);
    }
    let max = n - 2;
    if f <= 0.0 {
        (0, 0.0)
    } else if f >= max as f64 {
        (max, 1.0)
    } else {
        let i = f as usize;
        let i = i.min(max);
        (i, f - i as f64)
    }
}

fn estimate_noise<T: Pixel>(image: &[T], background: &[f64]) -> f64 {
    let residuals: Vec<f64> = image
        .par_iter()
        .zip(background.par_iter())
        .map(|(&px, &bg)| px.to_f64() - bg)
        .collect();

    sigma_clipped_std(&residuals, 3.0, 5)
}

fn sigma_clipped_std(data: &[f64], kappa: f64, iterations: usize) -> f64 {
    let mut mask = vec![true; data.len()];

    for _ in 0..iterations {
        let (mean, std) = masked_mean_std(data, &mask);
        if std <= 0.0 {
            return 0.0;
        }
        let lo = mean - kappa * std;
        let hi = mean + kappa * std;
        let changed = data
            .par_iter()
            .zip(mask.par_iter_mut())
            .map(|(&v, m)| {
                if *m && (v < lo || v > hi) {
                    *m = false;
                    1_u64
                } else {
                    0
                }
            })
            .sum::<u64>();
        if changed == 0 {
            return std;
        }
    }

    masked_mean_std(data, &mask).1
}

fn masked_mean_std(data: &[f64], mask: &[bool]) -> (f64, f64) {
    let (sum, count) = data
        .par_iter()
        .zip(mask.par_iter())
        .filter(|(_, &m)| m)
        .map(|(&v, _)| (v, 1_u64))
        .reduce(|| (0.0, 0), |(s1, c1), (s2, c2)| (s1 + s2, c1 + c2));

    if count < 2 {
        return (sum / count.max(1) as f64, 0.0);
    }
    let mean = sum / count as f64;

    let var_sum: f64 = data
        .par_iter()
        .zip(mask.par_iter())
        .filter(|(_, &m)| m)
        .map(|(&v, _)| {
            let d = v - mean;
            d * d
        })
        .sum();
    (mean, libm::sqrt(var_sum / (count - 1) as f64))
}

fn median(buf: &mut [f64]) -> f64 {
    buf.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = buf.len();
    if n % 2 == 1 {
        buf[n / 2]
    } else {
        (buf[n / 2 - 1] + buf[n / 2]) * 0.5
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_image_background_is_constant() {
        let width = 128;
        let height = 128;
        let image = vec![100.0_f32; width * height];
        let bg = estimate_background(&image, width, height, 32);

        for &v in &bg.map {
            assert!((v - 100.0).abs() < 1e-6, "expected ~100, got {v}");
        }
        assert!(bg.noise < 1e-6, "expected ~0 noise, got {}", bg.noise);
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

        let bg = estimate_background(&image, width, height, 32);

        let center = bg.map[128 * width + 128];
        let expected = 256.0;
        assert!(
            (center - expected).abs() < 20.0,
            "center bg {center} not near expected {expected}"
        );
    }

    #[test]
    fn star_does_not_dominate_background() {
        let width = 128;
        let height = 128;
        let mut image = vec![100.0_f32; width * height];
        image[64 * width + 64] = 50000.0;

        let bg = estimate_background(&image, width, height, 32);

        let star_bg = bg.map[64 * width + 64];
        assert!(
            (star_bg - 100.0).abs() < 10.0,
            "background at star position {star_bg} too far from 100"
        );
    }

    #[test]
    fn sigma_clipped_std_rejects_outliers() {
        let mut data = vec![0.0_f64; 1000];
        for (i, v) in data.iter_mut().enumerate() {
            *v = (i as f64 * 0.001) - 0.5;
        }
        data[500] = 1000.0;
        data[501] = -1000.0;

        let std = sigma_clipped_std(&data, 3.0, 5);
        assert!(std < 1.0, "sigma-clipped std {std} should be small");
    }

    #[test]
    fn median_odd() {
        let mut buf = vec![5.0, 1.0, 3.0];
        assert_eq!(median(&mut buf), 3.0);
    }

    #[test]
    fn median_even() {
        let mut buf = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(median(&mut buf), 2.5);
    }
}
