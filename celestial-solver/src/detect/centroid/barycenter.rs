// The centre of a region weighted by its pixels above a threshold.

use super::super::components::StarRegion;
use super::super::stats::{median_in_place, stddev};
use super::super::Pixel;

pub(super) fn thresholded_barycenter<T: Pixel>(
    image: &[T],
    width: usize,
    region: &StarRegion,
    local_bg: f64,
) -> (f64, f64, f64) {
    let (values, coords) = values_above_background(image, width, region, local_bg);
    let Some(med) = median_in_place(&mut values.clone()) else {
        return (0.0, 0.0, 0.0);
    };
    let std = stddev(&values);
    let lo = med + 1.5 * std;

    let max_val = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max_val <= lo || max_val <= 0.0 {
        return flux_weighted_center(&values, &coords);
    }

    let (sum_x, sum_y, sum_w) = moments_above(&values, &coords, lo, max_val - lo);
    if sum_w > 0.0 {
        (sum_x / sum_w, sum_y / sum_w, sum_w)
    } else {
        flux_weighted_center(&values, &coords)
    }
}

// The finite values of `region` more than 1e-7 above `local_bg`, with their pixels.
fn values_above_background<T: Pixel>(
    image: &[T],
    width: usize,
    region: &StarRegion,
    local_bg: f64,
) -> (Vec<f64>, Vec<(usize, usize)>) {
    let mut values = Vec::with_capacity(region.pixels.len());
    let mut coords = Vec::with_capacity(region.pixels.len());
    for &(px, py) in &region.pixels {
        let v = image[py * width + px].to_f64() - local_bg;
        if v.is_finite() && v > 1e-7 {
            values.push(v);
            coords.push((px, py));
        }
    }
    (values, coords)
}

// The sums of x, y and weight over the values at or above `lo`, each weighted by its
// height above `lo` divided by `range`.
fn moments_above(
    values: &[f64],
    coords: &[(usize, usize)],
    lo: f64,
    range: f64,
) -> (f64, f64, f64) {
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
    (sum_x, sum_y, sum_w)
}

fn flux_weighted_center(values: &[f64], coords: &[(usize, usize)]) -> (f64, f64, f64) {
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
