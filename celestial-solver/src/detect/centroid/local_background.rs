// The background median and noise around a region, from an annulus grown until its
// median is stable.

use super::super::components::StarRegion;
use super::super::stats::{mad, median_in_place};
use super::super::Pixel;

// The annulus stops growing once its median moves by at most this fraction of the
// background noise in one step. Each step adds a ring one pixel wide, so consecutive
// annuli share most of their pixels and their medians differ by much less than the
// standard error of either (1.25 / sqrt(n), 0.09 to 0.11 sigma for the 130 to 180
// pixels of the first annulus around a region 6 to 10 px across). The step change
// therefore measures how the background varies with radius, not sampling noise. Growth
// stops once one more ring moves the median by at most 0.1 sigma, a shift that would
// change a star's SNR by at most 0.1 against a default min_snr of 5.
const ANNULUS_CONVERGENCE: f64 = 0.1;

// Background median and robust noise (1.4826 x MAD) of the pixels around `region`,
// growing the annulus one pixel per step until the median is stable.
pub(super) fn iterative_local_background<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    region: &StarRegion,
    background_noise: f64,
) -> Option<(f64, f64)> {
    let tolerance = ANNULUS_CONVERGENCE * background_noise;
    let mut prev: Option<(f64, Vec<f64>)> = None;

    for inflate in (4_usize..).take(200) {
        let mut samples = sample_annular_frame(image, width, height, region, inflate);
        let Some(med) = median_in_place(&mut samples) else {
            break;
        };
        if let Some((prev_med, _)) = &prev {
            if libm::fabs(prev_med - med) <= tolerance {
                return Some((med, 1.4826 * mad(&samples, med)));
            }
        }
        prev = Some((med, samples));
    }

    prev.map(|(med, samples)| (med, 1.4826 * mad(&samples, med)))
}

// The finite pixels of the box around `region` grown by `inflate` pixels, outside the
// box grown by 1 pixel, both clipped to the image.
fn sample_annular_frame<T: Pixel>(
    image: &[T],
    width: usize,
    height: usize,
    region: &StarRegion,
    inflate: usize,
) -> Vec<f64> {
    let (x0, y0, x1, y1) = grown_box(region, inflate, width, height);
    let (ix0, iy0, ix1, iy1) = grown_box(region, 1, width, height);
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

// The bounding box of `region` grown by `by` pixels on each side and clipped to the
// image, as (x0, y0, x1, y1) inclusive.
fn grown_box(
    region: &StarRegion,
    by: usize,
    width: usize,
    height: usize,
) -> (usize, usize, usize, usize) {
    (
        region.x_min.saturating_sub(by),
        region.y_min.saturating_sub(by),
        (region.x_max + by).min(width - 1),
        (region.y_max + by).min(height - 1),
    )
}
