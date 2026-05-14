//! Connected-components labeling for the structure map.

/// A contiguous set of mask pixels, representing a candidate star region.
///
/// Produced by [`find_connected_components`]. The bounding box is inclusive; `pixels`
/// lists every `(x, y)` coordinate inside the component.
pub struct StarRegion {
    /// Minimum x coordinate in the region.
    pub x_min: usize,
    /// Minimum y coordinate in the region.
    pub y_min: usize,
    /// Maximum x coordinate (inclusive).
    pub x_max: usize,
    /// Maximum y coordinate (inclusive).
    pub y_max: usize,
    /// All `(x, y)` pixels in the region.
    pub pixels: Vec<(usize, usize)>,
}

/// Labels connected components in a binary mask and filters by size / coverage.
///
/// Uses 4-connectivity with union-find. Returns only regions whose pixel count falls
/// in `[min_size, max_size]` and whose pixel count ÷ bounding-box area ≥ `min_coverage`.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::components::find_connected_components;
///
/// let mask = vec![0u8; 128 * 128];
/// let regions = find_connected_components(&mask, 128, 128, 5, 2000, 0.5);
/// ```
pub fn find_connected_components(
    mask: &[u8],
    width: usize,
    height: usize,
    min_size: usize,
    max_size: usize,
    min_coverage: f64,
) -> Vec<StarRegion> {
    let n = width * height;
    let mut labels = vec![0_u32; n];
    let mut parent = vec![0_u32; 1];
    let mut rank = vec![0_u8; 1];
    let mut next_label = 1_u32;

    for y in 0..height {
        for x in 0..width {
            if mask[y * width + x] == 0 {
                continue;
            }

            let mut neighbors = [0_u32; 2];
            let mut nc = 0;

            if x > 0 && labels[y * width + x - 1] > 0 {
                neighbors[nc] = labels[y * width + x - 1];
                nc += 1;
            }
            if y > 0 && labels[(y - 1) * width + x] > 0 {
                let l = labels[(y - 1) * width + x];
                if nc == 0 || find(&parent, l) != find(&parent, neighbors[0]) {
                    neighbors[nc] = l;
                    nc += 1;
                }
            }
            if y > 0 && x > 0 && labels[(y - 1) * width + x - 1] > 0 {
                let l = labels[(y - 1) * width + x - 1];
                let r = find(&parent, l);
                let new = nc == 0
                    || (nc == 1 && r != find(&parent, neighbors[0]))
                    || (nc == 2
                        && r != find(&parent, neighbors[0])
                        && r != find(&parent, neighbors[1]));
                if new && nc < 2 {
                    neighbors[nc] = l;
                    nc += 1;
                } else if nc > 0 {
                    union(&mut parent, &mut rank, l, neighbors[0]);
                }
            }
            if y > 0 && x + 1 < width && labels[(y - 1) * width + x + 1] > 0 {
                let l = labels[(y - 1) * width + x + 1];
                let r = find(&parent, l);
                let new = nc == 0
                    || (nc == 1 && r != find(&parent, neighbors[0]))
                    || (nc == 2
                        && r != find(&parent, neighbors[0])
                        && r != find(&parent, neighbors[1]));
                if new && nc < 2 {
                    neighbors[nc] = l;
                    nc += 1;
                } else if nc > 0 {
                    union(&mut parent, &mut rank, l, neighbors[0]);
                }
            }

            if nc == 0 {
                parent.push(next_label);
                rank.push(0);
                labels[y * width + x] = next_label;
                next_label += 1;
            } else {
                let root = find(&parent, neighbors[0]);
                labels[y * width + x] = root;
                if nc == 2 {
                    union(&mut parent, &mut rank, neighbors[0], neighbors[1]);
                }
            }
        }
    }

    flatten_labels(&mut labels, &parent);
    collect_regions(
        &labels, width, height, min_size, max_size, min_coverage,
    )
}

fn find(parent: &[u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        x = parent[x as usize];
    }
    x
}

fn union(parent: &mut [u32], rank: &mut [u8], a: u32, b: u32) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra == rb {
        return;
    }
    if rank[ra as usize] < rank[rb as usize] {
        parent[ra as usize] = rb;
    } else if rank[ra as usize] > rank[rb as usize] {
        parent[rb as usize] = ra;
    } else {
        parent[rb as usize] = ra;
        rank[ra as usize] += 1;
    }
}

fn flatten_labels(labels: &mut [u32], parent: &[u32]) {
    for label in labels.iter_mut() {
        if *label > 0 {
            *label = find(parent, *label);
        }
    }
}

fn collect_regions(
    labels: &[u32],
    width: usize,
    height: usize,
    min_size: usize,
    max_size: usize,
    min_coverage: f64,
) -> Vec<StarRegion> {
    let mut map: std::collections::HashMap<u32, StarRegion> = std::collections::HashMap::new();

    for y in 0..height {
        for x in 0..width {
            let l = labels[y * width + x];
            if l == 0 {
                continue;
            }
            map.entry(l)
                .and_modify(|r| {
                    if x < r.x_min { r.x_min = x; }
                    if x > r.x_max { r.x_max = x; }
                    if y < r.y_min { r.y_min = y; }
                    if y > r.y_max { r.y_max = y; }
                    r.pixels.push((x, y));
                })
                .or_insert(StarRegion {
                    x_min: x,
                    y_min: y,
                    x_max: x,
                    y_max: y,
                    pixels: vec![(x, y)],
                });
        }
    }

    let mut rej_size = 0_usize;
    let mut rej_border = 0_usize;
    let mut rej_bbox = 0_usize;
    let mut rej_coverage = 0_usize;
    let total = map.len();

    let result: Vec<StarRegion> = map.into_values()
        .filter(|r| {
            let count = r.pixels.len();
            if count < min_size || count > max_size {
                rej_size += 1;
                return false;
            }
            if touches_border(r, width, height) {
                rej_border += 1;
                return false;
            }
            let bbox_w = r.x_max - r.x_min + 1;
            let bbox_h = r.y_max - r.y_min + 1;
            let area = bbox_w * bbox_h;
            if area > 0 {
                let aspect = bbox_w.max(bbox_h) as f64 / bbox_w.min(bbox_h).max(1) as f64;
                if aspect > 4.0 && count > 8 {
                    rej_bbox += 1;
                    return false;
                }
                let cov = count as f64 / area as f64;
                if count > 8 && cov < min_coverage {
                    rej_coverage += 1;
                    return false;
                }
            }
            true
        })
        .collect();

    log::debug!("[components] total={total} rej_size={rej_size} rej_border={rej_border} rej_bbox={rej_bbox} rej_coverage={rej_coverage} passed={}", result.len());
    result
}

fn touches_border(r: &StarRegion, width: usize, height: usize) -> bool {
    r.x_min == 0 || r.y_min == 0 || r.x_max >= width - 1 || r.y_max >= height - 1
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
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5);

        assert_eq!(regions.len(), 1);
        assert!(regions[0].pixels.len() >= 5);
    }

    #[test]
    fn rejects_border_blob() {
        let w = 64;
        let h = 64;
        let mask = mask_with_blob(w, h, 0, 0, 4);
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5);

        assert!(regions.is_empty(), "border blob should be rejected");
    }

    #[test]
    fn rejects_tiny_blob() {
        let w = 64;
        let h = 64;
        let mut mask = vec![0_u8; w * h];
        mask[32 * w + 32] = 1;
        mask[32 * w + 33] = 1;
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5);

        assert!(regions.is_empty(), "2-pixel blob should be rejected");
    }

    #[test]
    fn rejects_huge_blob() {
        let w = 128;
        let h = 128;
        let mask = mask_with_blob(w, h, 64, 64, 30);
        let count: usize = mask.iter().map(|&v| v as usize).sum();
        let regions = find_connected_components(&mask, w, h, 5, 100, 0.5);

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
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5);

        assert_eq!(regions.len(), 2, "should find 2 separate blobs");
    }

    #[test]
    fn empty_mask_gives_no_regions() {
        let w = 64;
        let h = 64;
        let mask = vec![0_u8; w * h];
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.5);
        assert!(regions.is_empty());
    }

    #[test]
    fn rejects_elongated_streak() {
        let w = 128;
        let h = 128;
        let mut mask = vec![0_u8; w * h];
        for x in 20..80 {
            mask[64 * w + x] = 1;
        }
        let regions = find_connected_components(&mask, w, h, 5, 2000, 0.785);
        assert!(
            regions.is_empty(),
            "1px-tall streak should fail coverage filter"
        );
    }
}
