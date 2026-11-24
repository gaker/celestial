use std::cmp::Ordering;
use std::collections::HashSet;

use celestial_core::constants::{DEG_TO_RAD, RAD_TO_DEG};

use super::healpix::reduce_ra_deg;

#[derive(Debug, Clone)]
pub struct QuadStar {
    pub source_id: i64,
    pub x: f64,
    pub y: f64,
    pub mag: f32,
}

#[derive(Debug, Clone)]
pub struct Quad {
    pub hash: [f64; 4],
    pub star_indices: [usize; 4],
}

pub fn tan_project_star(
    ra_deg: f64,
    dec_deg: f64,
    center_ra_deg: f64,
    center_dec_deg: f64,
) -> Option<(f64, f64)> {
    if ![ra_deg, dec_deg, center_ra_deg, center_dec_deg]
        .iter()
        .all(|v| v.is_finite())
    {
        return None;
    }
    let (sd0, cd0) = libm::sincos(center_dec_deg * DEG_TO_RAD);
    let (sd, cd) = libm::sincos(dec_deg * DEG_TO_RAD);
    let (sdra, cdra) = libm::sincos((ra_deg - center_ra_deg) * DEG_TO_RAD);
    // denom is the cosine of the distance from the tangent point. Its rounding error is
    // about one ulp of 1, so at or below f64::EPSILON it does not separate points in
    // front of the tangent plane from points on or behind it.
    let denom = sd * sd0 + cd * cd0 * cdra;
    if denom <= f64::EPSILON {
        return None;
    }
    let xi = (cd * sdra) / denom;
    let eta = (sd * cd0 - cd * sd0 * cdra) / denom;
    Some((xi, eta))
}

pub fn tan_deproject_star(
    xi_rad: f64,
    eta_rad: f64,
    center_ra_deg: f64,
    center_dec_deg: f64,
) -> (f64, f64) {
    let (sd0, cd0) = libm::sincos(center_dec_deg * DEG_TO_RAD);
    let rho = libm::sqrt(xi_rad * xi_rad + eta_rad * eta_rad);
    if rho < 1e-20 {
        return (reduce_ra_deg(center_ra_deg), center_dec_deg);
    }
    let c = libm::atan(rho);
    let (sc, cc) = libm::sincos(c);
    let dec = libm::asin(cc * sd0 + eta_rad * sc * cd0 / rho) * RAD_TO_DEG;
    let ra =
        center_ra_deg + libm::atan2(xi_rad * sc, rho * cd0 * cc - eta_rad * sd0 * sc) * RAD_TO_DEG;
    (reduce_ra_deg(ra), dec)
}

// Each candidate spine pair of a quad's four slots, with the two slots left for C and D.
const SPINE_PAIRS: [([usize; 2], [usize; 2]); 6] = [
    ([0, 1], [2, 3]),
    ([0, 2], [1, 3]),
    ([0, 3], [1, 2]),
    ([1, 2], [0, 3]),
    ([1, 3], [0, 2]),
    ([2, 3], [0, 1]),
];

// Labels a quad from its own geometry (Lang et al. 2010, AJ 139, 1782, Sec. 2.2): A and B
// are the most widely separated pair. Returns the hash and the star indices in A, B, C, D
// order, or None when no pair gives a finite hash.
fn canonical_quad(stars: &[QuadStar], indices: [usize; 4]) -> Option<([f64; 4], [usize; 4])> {
    let mut best: Option<(f64, [f64; 4], [usize; 4])> = None;
    for ([i, j], [k, l]) in SPINE_PAIRS {
        let (a, b) = (indices[i], indices[j]);
        let d2 = distance_sq(&stars[a], &stars[b]);
        if best.as_ref().is_some_and(|(best_d2, _, _)| d2 < *best_d2) {
            continue;
        }
        let Some((hash, order)) = orient_spine(stars, a, b, indices[k], indices[l]) else {
            continue;
        };
        // An exact tie in separation, such as a rectangle's diagonals, takes the smaller
        // hash, so the result does not depend on the order of `indices`.
        let replace = match &best {
            None => true,
            Some((best_d2, best_hash, _)) => {
                d2 > *best_d2 || cmp_hash(&hash, best_hash) == Ordering::Less
            }
        };
        if replace {
            best = Some((d2, hash, order));
        }
    }
    best.map(|(_, hash, order)| (hash, order))
}

// Chooses which end of the spine is A. The paper places B at (1, 1) and requires
// xC + xD <= 1 in that frame. Here B is at (1, 0), where swapping A and B maps (x, y) to
// (1 - x, -y), so requiring cx + dx <= 1 in this frame breaks the same symmetry with a
// different cut through hash space. Comparing the sums of both orientations makes that
// choice even when rounding puts both sums on one side of 1. An exact tie takes the
// smaller hash.
fn orient_spine(
    stars: &[QuadStar],
    a: usize,
    b: usize,
    c: usize,
    d: usize,
) -> Option<([f64; 4], [usize; 4])> {
    let forward = normalize_quad_ordered(stars, a, b, c, d)?;
    let reverse = normalize_quad_ordered(stars, b, a, c, d)?;
    let sum = |hash: &[f64; 4]| hash[0] + hash[2];
    let keep_forward = match sum(&forward.0).total_cmp(&sum(&reverse.0)) {
        Ordering::Less => true,
        Ordering::Greater => false,
        Ordering::Equal => cmp_hash(&forward.0, &reverse.0) != Ordering::Greater,
    };
    Some(if keep_forward { forward } else { reverse })
}

fn distance_sq(s1: &QuadStar, s2: &QuadStar) -> f64 {
    let dx = s2.x - s1.x;
    let dy = s2.y - s1.y;
    dx * dx + dy * dy
}

fn cmp_hash(h1: &[f64; 4], h2: &[f64; 4]) -> Ordering {
    h1.iter()
        .zip(h2)
        .map(|(v1, v2)| v1.total_cmp(v2))
        .find(|o| o.is_ne())
        .unwrap_or(Ordering::Equal)
}

// Hash of C and D in the frame with A at (0, 0) and B at (1, 0), with C and D ordered so
// that cx <= dx (ties by cy <= dy).
fn normalize_quad_ordered(
    stars: &[QuadStar],
    a: usize,
    b: usize,
    c: usize,
    d: usize,
) -> Option<([f64; 4], [usize; 4])> {
    let (cx, cy) = spine_frame(&stars[a], &stars[b], &stars[c])?;
    let (dx, dy) = spine_frame(&stars[a], &stars[b], &stars[d])?;
    let (hash, order) = if cx < dx || (cx == dx && cy <= dy) {
        ([cx, cy, dx, dy], [a, b, c, d])
    } else {
        ([dx, dy, cx, cy], [a, b, d, c])
    };
    hash.iter().all(|v| v.is_finite()).then_some((hash, order))
}

// Coordinates of `p` in the frame with `sa` at (0, 0) and `sb` at (1, 0), or None when
// `sa` and `sb` coincide.
fn spine_frame(sa: &QuadStar, sb: &QuadStar, p: &QuadStar) -> Option<(f64, f64)> {
    let dx_ab = sb.x - sa.x;
    let dy_ab = sb.y - sa.y;
    let scale = dx_ab * dx_ab + dy_ab * dy_ab;
    if scale < 1e-30 {
        return None;
    }
    let px = p.x - sa.x;
    let py = p.y - sa.y;
    Some((
        (px * dx_ab + py * dy_ab) / scale,
        (py * dx_ab - px * dy_ab) / scale,
    ))
}

fn knn_indices(stars: &[QuadStar], seed: usize, k: usize) -> Vec<usize> {
    let sx = stars[seed].x;
    let sy = stars[seed].y;
    let mut dists: Vec<(usize, f64)> = stars
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != seed)
        .map(|(i, s)| {
            let dx = s.x - sx;
            let dy = s.y - sy;
            (i, dx * dx + dy * dy)
        })
        .collect();
    dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
    dists.truncate(k);
    dists.into_iter().map(|(i, _)| i).collect()
}

pub fn neighbor_quads(stars: &[QuadStar], k_neighbors: usize, min_spine_sq: f64) -> Vec<Quad> {
    let n = stars.len();
    if n < 4 {
        return Vec::new();
    }
    let k = k_neighbors.min(n - 1);
    let mut quads = Vec::new();
    let mut seen = HashSet::new();

    for seed in 0..n {
        for idx in seed_combinations(seed, &knn_indices(stars, seed, k)) {
            if !seen.insert(idx) {
                continue;
            }
            let Some((hash, order)) = canonical_quad(stars, idx) else {
                continue;
            };
            if distance_sq(&stars[order[0]], &stars[order[1]]) >= min_spine_sq {
                quads.push(Quad {
                    hash,
                    star_indices: order,
                });
            }
        }
    }
    quads
}

// Every set of `seed` and three of its neighbours, each sorted by star index.
fn seed_combinations(seed: usize, neighbors: &[usize]) -> Vec<[usize; 4]> {
    let nn = neighbors.len();
    let mut sets = Vec::new();
    for a in 0..nn {
        for b in (a + 1)..nn {
            for c in (b + 1)..nn {
                let mut idx = [seed, neighbors[a], neighbors[b], neighbors[c]];
                idx.sort();
                sets.push(idx);
            }
        }
    }
    sets
}

// Bins one hash value into a 0.01-wide cell.
fn hash_bin(value: f64) -> i32 {
    libm::floor(value * 100.0) as i32
}

// Packs four cells into a key, 16 bits each, offset by 500 so small negative cells stay
// distinct.
fn pack_bins(bins: [i32; 4]) -> u64 {
    bins.iter().enumerate().fold(0, |key, (i, &bin)| {
        key | ((bin.wrapping_add(500) as u64 & 0xFFFF) << (16 * i))
    })
}

pub fn discrete_hash_key(hash: &[f64; 4]) -> u64 {
    pack_bins(hash.map(hash_bin))
}

pub fn discrete_hash_key_neighbors(hash: &[f64; 4]) -> [u64; 81] {
    let base = hash.map(hash_bin);
    let mut keys = [0u64; 81];
    for (n, key) in keys.iter_mut().enumerate() {
        // The base-3 digits of n, last value fastest, give each cell's offset of -1, 0 or +1.
        let mut digits = n;
        let mut bins = base;
        for bin in bins.iter_mut().rev() {
            *bin = bin.wrapping_add((digits % 3) as i32 - 1);
            digits /= 3;
        }
        *key = pack_bins(bins);
    }
    keys
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_stars(coords: &[(f64, f64, f32)]) -> Vec<QuadStar> {
        coords
            .iter()
            .enumerate()
            .map(|(i, &(x, y, mag))| QuadStar {
                source_id: i as i64,
                x,
                y,
                mag,
            })
            .collect()
    }

    fn normalize_quad(
        stars: &[QuadStar],
        a: usize,
        b: usize,
        c: usize,
        d: usize,
    ) -> Option<[f64; 4]> {
        normalize_quad_ordered(stars, a, b, c, d).map(|(hash, _)| hash)
    }

    fn find_spine(stars: &[QuadStar], indices: [usize; 4]) -> [usize; 4] {
        canonical_quad(stars, indices).unwrap().1
    }

    #[test]
    fn test_find_spine_horizontal() {
        let stars = make_stars(&[
            (0.0, 0.0, 5.0),
            (10.0, 0.0, 5.0),
            (3.0, 1.0, 5.0),
            (7.0, -1.0, 5.0),
        ]);
        let [a, b, _, _] = find_spine(&stars, [0, 1, 2, 3]);
        assert_eq!(a, 0);
        assert_eq!(b, 1);
    }

    #[test]
    fn test_find_spine_diagonal() {
        let stars = make_stars(&[
            (1.0, 1.0, 5.0),
            (2.0, 2.0, 5.0),
            (0.0, 0.0, 5.0),
            (10.0, 10.0, 5.0),
        ]);
        let [a, b, _, _] = find_spine(&stars, [0, 1, 2, 3]);
        assert_eq!(a, 2);
        assert_eq!(b, 3);
    }

    #[test]
    fn test_normalize_quad_unit_square() {
        let stars = make_stars(&[
            (0.0, 0.0, 5.0),
            (1.0, 0.0, 5.0),
            (0.3, 0.4, 5.0),
            (0.7, -0.2, 5.0),
        ]);
        let hash = normalize_quad(&stars, 0, 1, 2, 3).unwrap();
        assert_eq!(hash, [0.3, 0.4, 0.7, -0.2]);
    }

    #[test]
    fn test_normalize_quad_scaled() {
        let scale = 5.0;
        let stars = make_stars(&[
            (0.0, 0.0, 5.0),
            (scale, 0.0, 5.0),
            (0.3 * scale, 0.4 * scale, 5.0),
            (0.7 * scale, -0.2 * scale, 5.0),
        ]);
        let hash = normalize_quad(&stars, 0, 1, 2, 3).unwrap();
        assert_eq!(hash, [0.3, 0.4, 0.7, -0.2]);
    }

    #[test]
    fn test_normalize_quad_rotation_invariance() {
        let cx = 0.3;
        let cy = 0.4;
        let dx = 0.7;
        let dy = -0.2;

        let angle = 0.73;
        let (sa, ca) = libm::sincos(angle);

        let rotate = |x: f64, y: f64| (x * ca - y * sa, x * sa + y * ca);

        let (ax, ay) = rotate(0.0, 0.0);
        let (bx, by) = rotate(1.0, 0.0);
        let (cx_r, cy_r) = rotate(cx, cy);
        let (dx_r, dy_r) = rotate(dx, dy);

        let offset_x = 3.5;
        let offset_y = -7.2;
        let scale = 4.0;

        let stars = make_stars(&[
            (ax * scale + offset_x, ay * scale + offset_y, 5.0),
            (bx * scale + offset_x, by * scale + offset_y, 5.0),
            (cx_r * scale + offset_x, cy_r * scale + offset_y, 5.0),
            (dx_r * scale + offset_x, dy_r * scale + offset_y, 5.0),
        ]);

        // Rounding in the rotation, offset and scale moves three values by 1 or 2 ulp from
        // (cx, cy, dx, dy). libm's sincos gives the same bits on every target.
        let hash = normalize_quad(&stars, 0, 1, 2, 3).unwrap();
        assert_eq!(
            hash,
            [
                0.30000000000000004,
                cy,
                0.7000000000000001,
                -0.20000000000000007
            ]
        );
    }

    #[test]
    fn test_normalize_quad_cd_canonical_order() {
        let stars = make_stars(&[
            (0.0, 0.0, 5.0),
            (1.0, 0.0, 5.0),
            (0.8, 0.3, 5.0),
            (0.2, 0.1, 5.0),
        ]);
        let hash_cd = normalize_quad(&stars, 0, 1, 2, 3).unwrap();
        let hash_dc = normalize_quad(&stars, 0, 1, 3, 2).unwrap();
        assert_eq!(hash_cd, hash_dc);
        assert!(hash_cd[0] <= hash_cd[2]);
    }

    #[test]
    fn test_normalize_quad_degenerate_returns_none() {
        let stars = make_stars(&[
            (1.0, 1.0, 5.0),
            (1.0, 1.0, 5.0),
            (2.0, 3.0, 5.0),
            (4.0, 5.0, 5.0),
        ]);
        assert!(normalize_quad(&stars, 0, 1, 2, 3).is_none());
    }

    #[test]
    fn test_tan_project_star_center() {
        assert_eq!(tan_project_star(45.0, 30.0, 45.0, 30.0), Some((0.0, 0.0)));
    }

    #[test]
    fn test_tan_project_star_behind_returns_none() {
        let result = tan_project_star(0.0, 0.0, 180.0, 0.0);
        assert!(result.is_none());
    }

    #[test]
    fn neighbor_quads_generates_from_local_neighborhoods() {
        let stars = make_stars(&[
            (0.0, 0.0, 5.0),
            (1.0, 0.0, 5.0),
            (0.5, 1.0, 5.0),
            (0.5, 0.5, 5.0),
            (10.0, 10.0, 5.0),
        ]);
        let quads = neighbor_quads(&stars, 3, 0.0);
        assert!(!quads.is_empty());
        for q in &quads {
            assert_eq!(q.star_indices.len(), 4);
        }
    }

    #[test]
    fn neighbor_quads_no_duplicates() {
        let stars = make_stars(&[
            (0.0, 0.0, 5.0),
            (1.0, 0.0, 5.0),
            (2.0, 0.0, 5.0),
            (0.0, 1.0, 5.0),
            (1.0, 1.0, 5.0),
            (2.0, 1.0, 5.0),
        ]);
        let quads = neighbor_quads(&stars, 5, 0.0);
        let mut hashes: Vec<[usize; 4]> = quads
            .iter()
            .map(|q| {
                let mut s = q.star_indices;
                s.sort();
                s
            })
            .collect();
        let before = hashes.len();
        hashes.sort();
        hashes.dedup();
        assert_eq!(before, hashes.len());
    }

    #[test]
    fn neighbor_quads_respects_min_spine() {
        let stars = make_stars(&[
            (0.0, 0.0, 5.0),
            (0.1, 0.0, 5.0),
            (0.05, 0.1, 5.0),
            (0.05, 0.05, 5.0),
        ]);
        let quads = neighbor_quads(&stars, 3, 1.0);
        assert!(quads.is_empty());
    }

    #[test]
    fn discrete_hash_key_neighbors_find_close_values() {
        let h1 = [0.300, 0.400, 0.700, -0.200];
        let h2 = [0.304, 0.397, 0.702, -0.203];
        let neighbors = discrete_hash_key_neighbors(&h1);
        let target = discrete_hash_key(&h2);
        assert!(neighbors.contains(&target));
    }

    #[test]
    fn discrete_hash_key_distant_values_no_overlap() {
        let h1 = [0.3, 0.4, 0.7, -0.2];
        let h2 = [0.5, 0.1, 0.3, 0.8];
        let neighbors = discrete_hash_key_neighbors(&h1);
        let target = discrete_hash_key(&h2);
        assert!(!neighbors.contains(&target));
    }

    // Canonical frame of the reference quad: A = (0, 0), B = (1, 0). Each hash value is
    // in the middle of its 0.01 bin, so rounding in the rotation does not change the key.
    const REFERENCE_QUAD: [(f64, f64); 4] =
        [(0.0, 0.0), (1.0, 0.0), (0.215, 0.335), (0.605, -0.245)];

    fn placed(points: &[(f64, f64)], deg: f64, scale: f64, offset: (f64, f64)) -> Vec<QuadStar> {
        let (s, c) = libm::sincos(deg * DEG_TO_RAD);
        let rotated: Vec<(f64, f64, f32)> = points
            .iter()
            .map(|&(x, y)| {
                let (rx, ry) = (x * c - y * s, x * s + y * c);
                (rx * scale + offset.0, ry * scale + offset.1, 5.0)
            })
            .collect();
        make_stars(&rotated)
    }

    fn single_quad_key(stars: &[QuadStar]) -> u64 {
        let quads = neighbor_quads(stars, 3, 0.0);
        assert_eq!(quads.len(), 1);
        discrete_hash_key(&quads[0].hash)
    }

    fn reference_key() -> u64 {
        discrete_hash_key(&[0.215, 0.335, 0.605, -0.245])
    }

    #[test]
    fn quad_hash_key_is_invariant_under_rotation() {
        for step in 0..24 {
            let deg = 15.0 * step as f64;
            let stars = placed(&REFERENCE_QUAD, deg, 1.0, (0.0, 0.0));
            assert_eq!(single_quad_key(&stars), reference_key(), "rotation {deg}");
        }
    }

    #[test]
    fn quad_hash_key_is_invariant_under_translation_and_scale() {
        for (scale, offset) in [
            (3.7, (12.5, -4.25)),
            (0.002, (-1.0, 0.5)),
            (250.0, (1e3, 2e3)),
        ] {
            let stars = placed(&REFERENCE_QUAD, 105.0, scale, offset);
            assert_eq!(single_quad_key(&stars), reference_key(), "scale {scale}");
        }
    }

    #[test]
    fn quad_hash_is_invariant_under_input_order() {
        let stars = placed(&REFERENCE_QUAD, 195.0, 3.7, (12.5, -4.25));
        let reference = neighbor_quads(&stars, 3, 0.0)[0].hash;
        for perm in permutations_of_4() {
            let permuted: Vec<QuadStar> = perm.iter().map(|&i| stars[i].clone()).collect();
            let quads = neighbor_quads(&permuted, 3, 0.0);
            assert_eq!(quads[0].hash, reference, "order {perm:?}");
            assert_eq!(
                discrete_hash_key(&quads[0].hash),
                reference_key(),
                "order {perm:?}"
            );
        }
    }

    fn permutations_of_4() -> Vec<[usize; 4]> {
        let mut out = Vec::new();
        for code in 0..256 {
            let p = [code & 3, (code >> 2) & 3, (code >> 4) & 3, code >> 6];
            let mut sorted = p;
            sorted.sort();
            if sorted == [0, 1, 2, 3] {
                out.push(p);
            }
        }
        assert_eq!(out.len(), 24);
        out
    }

    struct Lcg(u64);

    impl Lcg {
        fn uniform(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    // Counts image quads whose hash neighbourhood contains the catalog quad built from
    // the same four stars with the same A, B, C, D labels.
    fn own_quad_matches(catalog: &[Quad], image: &[Quad]) -> usize {
        let mut index: std::collections::HashMap<u64, Vec<[usize; 4]>> =
            std::collections::HashMap::new();
        for q in catalog {
            index
                .entry(discrete_hash_key(&q.hash))
                .or_default()
                .push(q.star_indices);
        }
        image
            .iter()
            .filter(|q| {
                discrete_hash_key_neighbors(&q.hash)
                    .iter()
                    .any(|k| index.get(k).is_some_and(|v| v.contains(&q.star_indices)))
            })
            .count()
    }

    #[test]
    fn field_quads_match_their_catalog_quads_at_any_rotation() {
        let mut rng = Lcg(5);
        let points: Vec<(f64, f64)> = (0..60).map(|_| (rng.uniform(), rng.uniform())).collect();
        let catalog = neighbor_quads(&placed(&points, 0.0, 1.0, (0.0, 0.0)), 10, 0.0);
        assert_eq!(catalog.len(), 4663);
        let mut found = Vec::new();
        for deg in [0.0, 90.0, 150.0, 180.0] {
            let image = neighbor_quads(&placed(&points, deg, 1.0, (0.0, 0.0)), 10, 0.0);
            found.push((deg, image.len(), own_quad_matches(&catalog, &image)));
        }
        assert_eq!(
            found,
            [
                (0.0, 4663, 4663),
                (90.0, 4663, 4663),
                (150.0, 4663, 4663),
                (180.0, 4663, 4663)
            ]
        );
    }

    #[test]
    fn tan_project_star_rejects_non_finite_input() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(tan_project_star(bad, 30.0, 45.0, 30.0), None);
            assert_eq!(tan_project_star(45.0, bad, 45.0, 30.0), None);
            assert_eq!(tan_project_star(45.0, 30.0, bad, 30.0), None);
            assert_eq!(tan_project_star(45.0, 30.0, 45.0, bad), None);
        }
    }

    #[test]
    fn tan_project_star_rejects_points_on_the_tangent_plane() {
        // cos(90 degrees) rounds to 6.1e-17 rather than 0.
        assert_eq!(tan_project_star(90.0, 0.0, 0.0, 0.0), None);
        assert_eq!(tan_project_star(0.0, 0.0, 0.0, 90.0), None);
    }

    #[test]
    fn tan_deproject_star_returns_ra_in_0_to_360() {
        for (xi, center_ra) in [(-0.01, 0.0), (0.01, 359.9), (-0.01, 0.2)] {
            let (ra, _) = tan_deproject_star(xi, 0.0, center_ra, 10.0);
            assert!(
                (0.0..360.0).contains(&ra),
                "xi {xi} centre {center_ra}: {ra}"
            );
        }
    }
}
