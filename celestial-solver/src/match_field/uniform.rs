//! Spatially uniform choice of quad stars (Lang et al. 2010, sec. 2.3): the brightest
//! stars of each cell of a grid, rather than the brightest stars overall, so that both
//! star lists sample the sky at the same density wherever the image falls.

use std::collections::HashMap;

// The size in pixels of the cells of a grid over a `width` x `height` image that holds
// `n_stars` quad stars in cells of about `k_neighbors` stars each: one cell is about the
// area of the `k_neighbors` nearest neighbours a quad is built from. The count along each
// axis is rounded, at least 1, so every cell has the same size.
pub(super) fn image_cell_px(
    width: usize,
    height: usize,
    n_stars: usize,
    k_neighbors: usize,
) -> (f64, f64) {
    let (w, h) = (width as f64, height as f64);
    let n_cells = (n_stars as f64 / k_neighbors.max(1) as f64).max(1.0);
    let side = libm::sqrt(w * h / n_cells);
    let nx = libm::round(w / side).max(1.0);
    let ny = libm::round(h / side).max(1.0);
    (w / nx, h / ny)
}

// Indices of the `n` points that are brightest within their cells of a grid `cell` wide
// and high, with a cell corner at the origin, in ascending order. `points` must be
// brightest first, and so is the result. A cell whose area holds `share(cell)` of a full
// cell's area of the region sampled gives that share of a full cell's points: its k-th
// brightest point ranks k / share, where a full cell's ranks k.
pub(super) fn brightest_in_cells(
    points: &[(f64, f64)],
    cell: (f64, f64),
    n: usize,
    share: impl Fn((i64, i64)) -> f64,
) -> Vec<usize> {
    let mut taken: HashMap<(i64, i64), usize> = HashMap::new();
    let mut ranked: Vec<(f64, usize)> = points
        .iter()
        .enumerate()
        .map(|(i, &(x, y))| {
            let key = (
                libm::floor(x / cell.0) as i64,
                libm::floor(y / cell.1) as i64,
            );
            let sweep = taken.entry(key).or_insert(0);
            *sweep += 1;
            (*sweep as f64 / share(key), i)
        })
        .collect();
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut chosen: Vec<usize> = ranked.into_iter().take(n).map(|(_, i)| i).collect();
    chosen.sort_unstable();
    chosen
}

// The share of each cell, `side` wide, of a grid with a cell corner at the origin, that
// lies inside the disc of `radius` about the origin: the exact area of the intersection
// over the cell's area.
pub(super) fn disc_cell_share(radius: f64, side: f64) -> impl Fn((i64, i64)) -> f64 {
    move |(i, j)| {
        let (x0, y0) = (i as f64 * side, j as f64 * side);
        let (x1, y1) = (x0 + side, y0 + side);
        let area = signed_disc_area(x1, y1, radius)
            - signed_disc_area(x0, y1, radius)
            - signed_disc_area(x1, y0, radius)
            + signed_disc_area(x0, y0, radius);
        area / (side * side)
    }
}

// The area of the disc of `radius` about the origin inside the rectangle with corners at
// the origin and (x, y), negative when exactly one of x and y is. The disc is symmetric
// about both axes, so the area of any rectangle follows by inclusion and exclusion.
fn signed_disc_area(x: f64, y: f64, radius: f64) -> f64 {
    let sign = libm::copysign(1.0, x) * libm::copysign(1.0, y);
    let (x, y) = (libm::fabs(x).min(radius), libm::fabs(y).min(radius));
    if x * x + y * y <= radius * radius {
        return sign * x * y;
    }
    // The rectangle's far edge at height y is inside the disc up to u0, and the circle
    // bounds the area from u0 to x.
    let u0 = libm::sqrt(radius * radius - y * y);
    sign * (y * u0 + circle_integral(x, radius) - circle_integral(u0, radius))
}

// The integral of sqrt(radius^2 - t^2) over t from 0 to `u`, for 0 <= u <= radius.
fn circle_integral(u: f64, radius: f64) -> f64 {
    let height = libm::sqrt((radius * radius - u * u).max(0.0));
    (u * height + radius * radius * libm::asin(u / radius)) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::PI;

    #[test]
    fn image_cell_px_gives_cells_of_k_neighbours() {
        // 100 stars in cells of 10: 10 cells of 316 px for a 1000 px square, rounded to
        // 3 x 3 cells.
        assert_eq!(
            image_cell_px(1000, 1000, 100, 10),
            (1000.0 / 3.0, 1000.0 / 3.0)
        );
        // A strip keeps one cell across its width.
        assert_eq!(image_cell_px(6000, 200, 100, 10), (6000.0 / 17.0, 200.0));
        // Fewer stars than one neighbourhood: one cell.
        assert_eq!(image_cell_px(600, 400, 5, 10), (600.0, 400.0));
    }

    fn full(_: (i64, i64)) -> f64 {
        1.0
    }

    #[test]
    fn brightest_in_cells_takes_one_point_per_cell_per_sweep() {
        // Three points in cell (0, 0), then one in each of cells (1, 0) and (0, 1).
        let points = [(0.1, 0.1), (0.2, 0.2), (0.3, 0.3), (1.5, 0.5), (0.5, 1.5)];
        assert_eq!(brightest_in_cells(&points, (1.0, 1.0), 3, full), [0, 3, 4]);
        assert_eq!(
            brightest_in_cells(&points, (1.0, 1.0), 4, full),
            [0, 1, 3, 4]
        );
    }

    #[test]
    fn brightest_in_cells_keeps_brightness_order_in_one_cell() {
        let points = [(0.5, 0.5), (0.1, 0.9), (0.9, 0.1)];
        assert_eq!(brightest_in_cells(&points, (1.0, 1.0), 2, full), [0, 1]);
    }

    #[test]
    fn brightest_in_cells_puts_negative_coordinates_in_their_own_cells() {
        let points = [(-0.5, 0.5), (-0.2, 0.2), (0.5, 0.5)];
        assert_eq!(brightest_in_cells(&points, (1.0, 1.0), 2, full), [0, 2]);
    }

    // Four points in cell (0, 0) and four in cell (1, 0), which has half the share,
    // alternating in brightness. The full cell's points rank 1 to 4 and the other's 2, 4,
    // 6 and 8, so six points take four from the full cell.
    #[test]
    fn brightest_in_cells_gives_a_cell_its_share_of_the_points() {
        let points: Vec<(f64, f64)> = (0..8)
            .map(|i| (0.5 + f64::from(i % 2), 0.1 * f64::from(i / 2) + 0.1))
            .collect();
        let half = |(i, _): (i64, i64)| if i == 0 { 1.0 } else { 0.5 };
        assert_eq!(
            brightest_in_cells(&points, (1.0, 1.0), 6, half),
            [0, 1, 2, 3, 4, 6]
        );
    }

    #[test]
    fn disc_cell_share_is_the_area_inside_the_disc() {
        let unit = disc_cell_share(1.0, 1.0);
        assert_eq!(unit((0, 0)), PI / 4.0);
        assert_eq!(unit((-1, -1)), PI / 4.0);
        assert_eq!(unit((1, 0)), 0.0);
        assert_eq!(disc_cell_share(10.0, 1.0)((3, -4)), 1.0);
        let edge = disc_cell_share(2.5, 1.0);
        assert_eq!(edge((2, 0)), edge((-3, -1)));
    }
}
