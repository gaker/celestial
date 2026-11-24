// Separable Gaussian smoothing with clamped edges.

use rayon::prelude::*;

pub(super) fn separable_gaussian(data: &[f32], w: usize, h: usize, sigma: f64) -> Vec<f32> {
    let kernel_f64 = build_gaussian_kernel(sigma);
    let inv_wt = 1.0 / kernel_f64.iter().sum::<f64>();
    let kernel: Vec<f32> = kernel_f64.iter().map(|&k| (k * inv_wt) as f32).collect();
    let rows = convolve_rows(data, w, h, &kernel);
    convolve_columns(&rows, w, h, &kernel)
}

fn convolve_rows(data: &[f32], w: usize, h: usize, kernel: &[f32]) -> Vec<f32> {
    let half = kernel.len() / 2;
    let mut temp = vec![0.0_f32; w * h];
    temp.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
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
    temp
}

fn convolve_columns(temp: &[f32], w: usize, h: usize, kernel: &[f32]) -> Vec<f32> {
    let half = kernel.len() / 2;
    let mut out = vec![0.0_f32; w * h];
    out.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gaussian_kernel_is_symmetric() {
        let kernel = build_gaussian_kernel(3.0);
        let n = kernel.len();
        for i in 0..n / 2 {
            assert_eq!(kernel[i], kernel[n - 1 - i]);
        }
    }
}
