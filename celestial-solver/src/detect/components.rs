//! Connected-components labeling for the structure map.

use std::collections::BTreeMap;

/// A contiguous set of mask pixels, representing a candidate star region.
///
/// Produced by [`find_connected_components`]. The bounding box is inclusive; `pixels`
/// lists every `(x, y)` coordinate inside the component.
pub(super) struct StarRegion {
    /// Minimum x coordinate in the region.
    pub(super) x_min: usize,
    /// Minimum y coordinate in the region.
    pub(super) y_min: usize,
    /// Maximum x coordinate (inclusive).
    pub(super) x_max: usize,
    /// Maximum y coordinate (inclusive).
    pub(super) y_max: usize,
    /// All `(x, y)` pixels in the region.
    pub(super) pixels: Vec<(usize, usize)>,
}

/// Labels connected components in a binary mask and filters by size / coverage.
///
/// Uses 4-connectivity with union-find. Returns only regions whose pixel count falls
/// in `[min_size, max_size]`, whose pixel count ÷ bounding-box area ≥ `min_coverage`,
/// and whose bounding box has at least `border_margin` pixels between it and every image
/// edge. A margin of 0 keeps regions that touch an edge.
pub(super) fn find_connected_components(
    mask: &[u8],
    width: usize,
    height: usize,
    min_size: usize,
    max_size: usize,
    min_coverage: f64,
    border_margin: usize,
) -> Vec<StarRegion> {
    let labels = label_pixels(mask, width, height);
    let limits = RegionLimits {
        min_size,
        max_size,
        min_coverage,
        border_margin,
    };
    collect_regions(&labels, width, height, &limits)
}

// The root label of every set pixel of `mask`, and 0 for every other pixel.
fn label_pixels(mask: &[u8], width: usize, height: usize) -> Vec<u32> {
    let mut labels = vec![0_u32; width * height];
    let mut forest = Forest::new();
    for y in 0..height {
        for x in 0..width {
            if mask[y * width + x] != 0 {
                labels[y * width + x] = label_pixel(&labels, &mut forest, x, y, width);
            }
        }
    }
    for label in labels.iter_mut() {
        if *label > 0 {
            *label = forest.find(*label);
        }
    }
    labels
}

// Up to two labels of already-labelled neighbours, with different roots when found.
struct Neighbours {
    labels: [u32; 2],
    count: usize,
}

impl Neighbours {
    fn push(&mut self, label: u32) {
        self.labels[self.count] = label;
        self.count += 1;
    }
}

// The label for the set pixel (x, y) from its labelled neighbours to the west, north,
// north-west and north-east: a new label when it has none, else the root of the first.
// Neighbours with different roots are joined.
fn label_pixel(labels: &[u32], forest: &mut Forest, x: usize, y: usize, width: usize) -> u32 {
    let mut found = Neighbours {
        labels: [0; 2],
        count: 0,
    };
    if x > 0 && labels[y * width + x - 1] > 0 {
        found.push(labels[y * width + x - 1]);
    }
    if y > 0 && labels[(y - 1) * width + x] > 0 {
        let l = labels[(y - 1) * width + x];
        if found.count == 0 || forest.find(l) != forest.find(found.labels[0]) {
            found.push(l);
        }
    }
    if y > 0 && x > 0 && labels[(y - 1) * width + x - 1] > 0 {
        forest.add_diagonal(&mut found, labels[(y - 1) * width + x - 1]);
    }
    if y > 0 && x + 1 < width && labels[(y - 1) * width + x + 1] > 0 {
        forest.add_diagonal(&mut found, labels[(y - 1) * width + x + 1]);
    }
    if found.count == 0 {
        return forest.add();
    }
    let root = forest.find(found.labels[0]);
    if found.count == 2 {
        forest.union(found.labels[0], found.labels[1]);
    }
    root
}

// Union-find over labels by rank. Label 0 is the background.
struct Forest {
    parent: Vec<u32>,
    rank: Vec<u8>,
}

impl Forest {
    fn new() -> Self {
        Self {
            parent: vec![0],
            rank: vec![0],
        }
    }

    // A new label that is its own root.
    fn add(&mut self) -> u32 {
        let label = self.parent.len() as u32;
        self.parent.push(label);
        self.rank.push(0);
        label
    }

    fn find(&self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            x = self.parent[x as usize];
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.rank[ra as usize] < self.rank[rb as usize] {
            self.parent[ra as usize] = rb;
        } else if self.rank[ra as usize] > self.rank[rb as usize] {
            self.parent[rb as usize] = ra;
        } else {
            self.parent[rb as usize] = ra;
            self.rank[ra as usize] += 1;
        }
    }

    // Adds the diagonal neighbour `l` to `found` when its root differs from theirs and
    // there is room; otherwise joins it to the first of them.
    fn add_diagonal(&mut self, found: &mut Neighbours, l: u32) {
        let r = self.find(l);
        let [n0, n1] = found.labels;
        let new = match found.count {
            0 => true,
            1 => r != self.find(n0),
            _ => r != self.find(n0) && r != self.find(n1),
        };
        if new && found.count < 2 {
            found.push(l);
        } else if found.count > 0 {
            self.union(l, n0);
        }
    }
}

impl StarRegion {
    fn at(x: usize, y: usize) -> Self {
        StarRegion {
            x_min: x,
            y_min: y,
            x_max: x,
            y_max: y,
            pixels: vec![(x, y)],
        }
    }

    fn include(&mut self, x: usize, y: usize) {
        self.x_min = self.x_min.min(x);
        self.x_max = self.x_max.max(x);
        self.y_min = self.y_min.min(y);
        self.y_max = self.y_max.max(y);
        self.pixels.push((x, y));
    }
}

// Regions keyed by root label. BTreeMap iteration runs in label order, which depends
// only on the mask, so the region order is the same on every call.
fn gather_regions(labels: &[u32], width: usize, height: usize) -> BTreeMap<u32, StarRegion> {
    let mut map: BTreeMap<u32, StarRegion> = BTreeMap::new();
    for y in 0..height {
        for x in 0..width {
            let l = labels[y * width + x];
            if l != 0 {
                map.entry(l)
                    .and_modify(|r| r.include(x, y))
                    .or_insert_with(|| StarRegion::at(x, y));
            }
        }
    }
    map
}

struct RegionLimits {
    min_size: usize,
    max_size: usize,
    min_coverage: f64,
    border_margin: usize,
}

// Why `collect_regions` drops a region, in the order it checks.
#[derive(Clone, Copy)]
enum Rejection {
    Size,
    Border,
    Elongated,
    Coverage,
}

fn collect_regions(
    labels: &[u32],
    width: usize,
    height: usize,
    limits: &RegionLimits,
) -> Vec<StarRegion> {
    let map = gather_regions(labels, width, height);
    let total = map.len();
    let mut rejected = [0_usize; 4];
    let result: Vec<StarRegion> = map
        .into_values()
        .filter(|r| match rejection(r, width, height, limits) {
            None => true,
            Some(why) => {
                rejected[why as usize] += 1;
                false
            }
        })
        .collect();
    let [rej_size, rej_border, rej_bbox, rej_coverage] = rejected;
    log::debug!("[components] total={total} rej_size={rej_size} rej_border={rej_border} rej_bbox={rej_bbox} rej_coverage={rej_coverage} passed={}", result.len());
    result
}

fn rejection(
    r: &StarRegion,
    width: usize,
    height: usize,
    limits: &RegionLimits,
) -> Option<Rejection> {
    let count = r.pixels.len();
    if count < limits.min_size || count > limits.max_size {
        return Some(Rejection::Size);
    }
    if near_border(r, width, height, limits.border_margin) {
        return Some(Rejection::Border);
    }
    let bbox_w = r.x_max - r.x_min + 1;
    let bbox_h = r.y_max - r.y_min + 1;
    let area = bbox_w * bbox_h;
    if area > 0 {
        let aspect = bbox_w.max(bbox_h) as f64 / bbox_w.min(bbox_h).max(1) as f64;
        if aspect > 4.0 && count > 8 {
            return Some(Rejection::Elongated);
        }
        let cov = count as f64 / area as f64;
        if count > 8 && cov < limits.min_coverage {
            return Some(Rejection::Coverage);
        }
    }
    None
}

// A region within `margin` pixels of an edge has a truncated PSF and a biased centroid.
fn near_border(r: &StarRegion, width: usize, height: usize, margin: usize) -> bool {
    r.x_min < margin
        || r.y_min < margin
        || r.x_max.saturating_add(margin) >= width
        || r.y_max.saturating_add(margin) >= height
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mask_with_blob(w: usize, h: usize, cx: usize, cy: usize, r: usize) -> Vec<u8> {
        let mut mask = vec![0_u8; w * h];
        let r_sq = (r * r) as i64;
        for y in 0..h {
            for x in 0..w {
                let dx = x as i64 - cx as i64;
                let dy = y as i64 - cy as i64;
                if dx * dx + dy * dy <= r_sq {
                    mask[y * w + x] = 1;
                }
            }
        }
        mask
    }

    #[test]
    fn single_blob_detected() {
        let w = 64;
        let h = 64;
        let mask = mask_with_blob(w, h, 32, 32, 4);
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);

        assert_eq!(regions.len(), 1);
        assert!(regions[0].pixels.len() >= 5);
    }

    #[test]
    fn rejects_border_blob() {
        let w = 64;
        let h = 64;
        let mask = mask_with_blob(w, h, 0, 0, 4);
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);

        assert!(regions.is_empty(), "border blob should be rejected");
    }

    // Blob of radius 2 at x = 7 has x_min = 5: kept with margin 5, dropped with margin 6.
    #[test]
    fn border_margin_is_measured_from_the_bounding_box() {
        let (w, h) = (64, 64);
        let mask = mask_with_blob(w, h, 7, 32, 2);
        assert_eq!(
            find_connected_components(&mask, w, h, 5, 2000, 0.5, 5).len(),
            1
        );
        assert_eq!(
            find_connected_components(&mask, w, h, 5, 2000, 0.5, 6).len(),
            0
        );
        let mask = mask_with_blob(w, h, 56, 32, 2);
        assert_eq!(
            find_connected_components(&mask, w, h, 5, 2000, 0.5, 5).len(),
            1
        );
        assert_eq!(
            find_connected_components(&mask, w, h, 5, 2000, 0.5, 6).len(),
            0
        );
    }

    #[test]
    fn margin_zero_keeps_a_region_touching_an_edge() {
        let (w, h) = (64, 64);
        let mask = mask_with_blob(w, h, 2, 32, 2);
        assert_eq!(
            find_connected_components(&mask, w, h, 5, 2000, 0.5, 0).len(),
            1
        );
        let mask = mask_with_blob(w, h, 3, 32, 2);
        assert_eq!(
            find_connected_components(&mask, w, h, 5, 2000, 0.5, 0).len(),
            1
        );
    }

    #[test]
    fn rejects_tiny_blob() {
        let w = 64;
        let h = 64;
        let mut mask = vec![0_u8; w * h];
        mask[32 * w + 32] = 1;
        mask[32 * w + 33] = 1;
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);

        assert!(regions.is_empty(), "2-pixel blob should be rejected");
    }

    #[test]
    fn rejects_huge_blob() {
        let w = 128;
        let h = 128;
        let mask = mask_with_blob(w, h, 64, 64, 30);
        let count: usize = mask.iter().map(|&v| v as usize).sum();
        let regions = find_connected_components(&mask, w, h, 5, 100, 0.5, 5);

        assert!(
            regions.is_empty(),
            "blob with {count} pixels should exceed max_size=100"
        );
    }

    #[test]
    fn two_separate_blobs() {
        let w = 128;
        let h = 128;
        let mut mask = mask_with_blob(w, h, 30, 30, 4);
        let mask2 = mask_with_blob(w, h, 90, 90, 4);
        for i in 0..mask.len() {
            mask[i] |= mask2[i];
        }
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);

        assert_eq!(regions.len(), 2, "should find 2 separate blobs");
    }

    #[test]
    fn empty_mask_gives_no_regions() {
        let w = 64;
        let h = 64;
        let mask = vec![0_u8; w * h];
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);
        assert!(regions.is_empty());
    }

    type RegionKey = (usize, usize, usize, usize, Vec<(usize, usize)>);

    fn region_keys(regions: &[StarRegion]) -> Vec<RegionKey> {
        regions
            .iter()
            .map(|r| (r.x_min, r.y_min, r.x_max, r.y_max, r.pixels.clone()))
            .collect()
    }

    // 36 identical blobs on a 6 x 6 grid.
    fn blob_grid_mask() -> (Vec<u8>, usize, usize) {
        let (w, h) = (192, 192);
        let mut mask = vec![0_u8; w * h];
        for gy in 0..6 {
            for gx in 0..6 {
                let blob = mask_with_blob(w, h, 16 + 32 * gx, 16 + 32 * gy, 4);
                for (m, b) in mask.iter_mut().zip(&blob) {
                    *m |= b;
                }
            }
        }
        (mask, w, h)
    }

    #[test]
    fn region_order_is_identical_across_calls() {
        let (mask, w, h) = blob_grid_mask();
        let first = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);
        assert_eq!(first.len(), 36);
        for _ in 0..20 {
            let again = find_connected_components(&mask, w, h, 5, 2000, 0.5, 5);
            assert_eq!(region_keys(&again), region_keys(&first));
        }
    }

    #[test]
    fn rejects_elongated_streak() {
        let w = 128;
        let h = 128;
        let mut mask = vec![0_u8; w * h];
        for x in 20..80 {
            mask[64 * w + x] = 1;
        }
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.785, 5);
        assert!(
            regions.is_empty(),
            "1px-tall streak should fail coverage filter"
        );
    }
}
