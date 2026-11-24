//! Uniform grid over the image that returns the reference stars whose Gaussian reaches a
//! pixel, so scoring a detection reads a few references instead of all of them.

use super::affine::Point;

pub(super) struct RefGrid {
    cell: f64,
    cols: usize,
    rows: usize,
    // Reference indices per cell in row-major order: cell `c` holds
    // `members[offsets[c]..offsets[c + 1]]`, ascending.
    offsets: Vec<usize>,
    members: Vec<usize>,
}

impl RefGrid {
    // Each disk is a reference's predicted pixel and the radius beyond which its
    // Gaussian is not evaluated. Cells are sized so that the image has about one cell
    // per reference.
    pub(super) fn new(width: f64, height: f64, disks: &[(Point, f64)]) -> Self {
        let n = disks.len().max(1) as f64;
        let cell = libm::sqrt(width * height / n).max(1.0);
        let cols = (libm::ceil(width / cell) as usize).max(1);
        let rows = (libm::ceil(height / cell) as usize).max(1);
        let mut grid = Self {
            cell,
            cols,
            rows,
            offsets: vec![0; cols * rows + 1],
            members: Vec::new(),
        };
        grid.fill(disks);
        grid
    }

    // Counting sort of (cell, reference) entries, so each cell lists its references in
    // ascending index order.
    fn fill(&mut self, disks: &[(Point, f64)]) {
        for disk in disks {
            for c in self.cells_of(disk) {
                self.offsets[c + 1] += 1;
            }
        }
        for c in 0..self.cols * self.rows {
            self.offsets[c + 1] += self.offsets[c];
        }
        let mut next = self.offsets.clone();
        self.members = vec![0; self.offsets[self.cols * self.rows]];
        for (k, disk) in disks.iter().enumerate() {
            for c in self.cells_of(disk) {
                self.members[next[c]] = k;
                next[c] += 1;
            }
        }
    }

    // Every cell that the disk's bounding box overlaps.
    fn cells_of(&self, &(centre, radius): &(Point, f64)) -> impl Iterator<Item = usize> {
        let (c0, c1) = self.span(centre[0] - radius, centre[0] + radius, self.cols);
        let (r0, r1) = self.span(centre[1] - radius, centre[1] + radius, self.rows);
        let cols = self.cols;
        (r0..r1).flat_map(move |r| (c0..c1).map(move |c| r * cols + c))
    }

    // Half-open range of cell indices covering [lo, hi], clamped to the grid; empty when
    // the interval misses the grid.
    fn span(&self, lo: f64, hi: f64, n: usize) -> (usize, usize) {
        let first = libm::floor(lo / self.cell);
        let last = libm::floor(hi / self.cell);
        if last < 0.0 || first >= n as f64 {
            return (0, 0);
        }
        (first.max(0.0) as usize, (last as usize).min(n - 1) + 1)
    }

    // The references whose disks may contain `x`; empty outside the image.
    pub(super) fn candidates(&self, x: Point) -> &[usize] {
        let (c0, c1) = self.span(x[0], x[0], self.cols);
        let (r0, r1) = self.span(x[1], x[1], self.rows);
        if c1 <= c0 || r1 <= r0 {
            return &[];
        }
        let c = r0 * self.cols + c0;
        &self.members[self.offsets[c]..self.offsets[c + 1]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_reference_is_listed_in_every_cell_its_disk_overlaps() {
        // Four references on a 100 x 100 image make 50 px cells.
        let disks = [
            ([10.0, 10.0], 5.0),
            ([50.0, 50.0], 1.0),
            ([90.0, 10.0], 5.0),
            ([10.0, 90.0], 200.0),
        ];
        let grid = RefGrid::new(100.0, 100.0, &disks);
        assert_eq!(grid.candidates([0.0, 0.0]), [0, 1, 3]);
        assert_eq!(grid.candidates([99.0, 0.0]), [1, 2, 3]);
        assert_eq!(grid.candidates([0.0, 99.0]), [1, 3]);
        assert_eq!(grid.candidates([99.0, 99.0]), [1, 3]);
    }

    #[test]
    fn points_off_the_grid_have_no_candidates() {
        let grid = RefGrid::new(100.0, 100.0, &[([50.0, 50.0], 500.0)]);
        assert_eq!(grid.candidates([-1.0, 50.0]), [0usize; 0]);
        assert_eq!(grid.candidates([50.0, 100.0]), [0usize; 0]);
    }
}
