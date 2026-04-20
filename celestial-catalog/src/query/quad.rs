use celestial_core::constants::{DEG_TO_RAD, RAD_TO_DEG};
use celestial_time::JulianDate;

use super::catalog::Catalog;
use super::cone::{cone_search, ConeSearchParams};

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
    let (sd0, cd0) = libm::sincos(center_dec_deg * DEG_TO_RAD);
    let (sd, cd) = libm::sincos(dec_deg * DEG_TO_RAD);
    let (sdra, cdra) = libm::sincos((ra_deg - center_ra_deg) * DEG_TO_RAD);
    let denom = sd * sd0 + cd * cd0 * cdra;
    if denom <= 0.0 {
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
        return (center_ra_deg, center_dec_deg);
    }
    let c = libm::atan(rho);
    let (sc, cc) = libm::sincos(c);
    let dec = libm::asin(cc * sd0 + eta_rad * sc * cd0 / rho) * RAD_TO_DEG;
    let ra = center_ra_deg
        + libm::atan2(xi_rad * sc, rho * cd0 * cc - eta_rad * sd0 * sc) * RAD_TO_DEG;
    (ra, dec)
}

pub fn find_spine(stars: &[QuadStar], indices: [usize; 4]) -> (usize, usize, usize, usize) {
    let mut best_dist_sq = -1.0_f64;
    let mut a_idx = 0;
    let mut b_idx = 1;

    for i in 0..4 {
        for j in (i + 1)..4 {
            let si = &stars[indices[i]];
            let sj = &stars[indices[j]];
            let dx = sj.x - si.x;
            let dy = sj.y - si.y;
            let d2 = dx * dx + dy * dy;
            if d2 > best_dist_sq {
                best_dist_sq = d2;
                a_idx = i;
                b_idx = j;
            }
        }
    }

    let a = indices[a_idx];
    let b = indices[b_idx];

    let (a, b) = orient_spine(&stars[a], a, &stars[b], b);

    let mut others = [0usize; 2];
    let mut k = 0;
    for (i, &val) in indices.iter().enumerate() {
        if i != a_idx && i != b_idx {
            others[k] = val;
            k += 1;
        }
    }

    (a, b, others[0], others[1])
}

pub fn orient_spine(sa: &QuadStar, a: usize, sb: &QuadStar, b: usize) -> (usize, usize) {
    if sa.x < sb.x {
        (a, b)
    } else if sa.x > sb.x {
        (b, a)
    } else if sa.y < sb.y {
        (a, b)
    } else if sa.y > sb.y {
        (b, a)
    } else {
        (a.min(b), a.max(b))
    }
}

pub fn normalize_quad(stars: &[QuadStar], a: usize, b: usize, c: usize, d: usize) -> Option<[f64; 4]> {
    normalize_quad_ordered(stars, a, b, c, d).map(|(hash, _)| hash)
}

pub fn normalize_quad_ordered(
    stars: &[QuadStar],
    a: usize,
    b: usize,
    c: usize,
    d: usize,
) -> Option<([f64; 4], [usize; 4])> {
    let ax = stars[a].x;
    let ay = stars[a].y;
    let dx_ab = stars[b].x - ax;
    let dy_ab = stars[b].y - ay;
    let scale = dx_ab * dx_ab + dy_ab * dy_ab;
    if scale < 1e-30 {
        return None;
    }

    let px_c = stars[c].x - ax;
    let py_c = stars[c].y - ay;
    let cx = (px_c * dx_ab + py_c * dy_ab) / scale;
    let cy = (py_c * dx_ab - px_c * dy_ab) / scale;

    let px_d = stars[d].x - ax;
    let py_d = stars[d].y - ay;
    let dx_val = (px_d * dx_ab + py_d * dy_ab) / scale;
    let dy_val = (py_d * dx_ab - px_d * dy_ab) / scale;

    let (cx, cy, dx_val, dy_val, c_out, d_out) =
        if cx < dx_val || (cx == dx_val && cy <= dy_val) {
            (cx, cy, dx_val, dy_val, c, d)
        } else {
            (dx_val, dy_val, cx, cy, d, c)
        };

    Some(([cx, cy, dx_val, dy_val], [a, b, c_out, d_out]))
}

pub fn combinations_4(n: usize) -> Combinations4 {
    let start = if n >= 4 { [0, 1, 2, 3] } else { [0, 0, 0, 0] };
    Combinations4 {
        n,
        current: start,
        done: n < 4,
    }
}

pub struct Combinations4 {
    n: usize,
    current: [usize; 4],
    done: bool,
}

impl Iterator for Combinations4 {
    type Item = [usize; 4];

    fn next(&mut self) -> Option<[usize; 4]> {
        if self.done {
            return None;
        }
        let result = self.current;

        let mut i = 3;
        loop {
            self.current[i] += 1;
            if self.current[i] <= self.n - 4 + i {
                break;
            }
            if i == 0 {
                self.done = true;
                return Some(result);
            }
            i -= 1;
        }
        for j in (i + 1)..4 {
            self.current[j] = self.current[j - 1] + 1;
        }

        Some(result)
    }
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
    dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    dists.truncate(k);
    dists.into_iter().map(|(i, _)| i).collect()
}

pub fn neighbor_quads(
    stars: &[QuadStar],
    k_neighbors: usize,
    min_spine_sq: f64,
) -> Vec<Quad> {
    let n = stars.len();
    if n < 4 {
        return Vec::new();
    }
    let k = k_neighbors.min(n - 1);
    let mut quads = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for seed in 0..n {
        let neighbors = knn_indices(stars, seed, k);
        if neighbors.len() < 3 {
            continue;
        }
        let nn = neighbors.len();
        for a in 0..nn {
            for b in (a + 1)..nn {
                for c in (b + 1)..nn {
                    let mut idx = [seed, neighbors[a], neighbors[b], neighbors[c]];
                    idx.sort();
                    if !seen.insert(idx) {
                        continue;
                    }
                    let (sa, sb, sc, sd) = find_spine(stars, idx);
                    let dx = stars[sb].x - stars[sa].x;
                    let dy = stars[sb].y - stars[sa].y;
                    if dx * dx + dy * dy < min_spine_sq {
                        continue;
                    }
                    if let Some((hash, ordered)) = normalize_quad_ordered(stars, sa, sb, sc, sd) {
                        quads.push(Quad {
                            hash,
                            star_indices: ordered,
                        });
                    }
                }
            }
        }
    }
    quads
}

pub fn discrete_hash_key(hash: &[f64; 4]) -> u64 {
    let quantize = |v: f64| ((v * 100.0).floor() as i32 + 500) as u64 & 0xFFFF;
    quantize(hash[0]) | (quantize(hash[1]) << 16) | (quantize(hash[2]) << 32) | (quantize(hash[3]) << 48)
}

pub fn discrete_hash_key_neighbors(hash: &[f64; 4]) -> [u64; 81] {
    let base = |v: f64| (v * 100.0).floor() as i32;
    let pack = |a: i32, b: i32, c: i32, d: i32| -> u64 {
        let qa = (a + 500) as u64 & 0xFFFF;
        let qb = (b + 500) as u64 & 0xFFFF;
        let qc = (c + 500) as u64 & 0xFFFF;
        let qd = (d + 500) as u64 & 0xFFFF;
        qa | (qb << 16) | (qc << 32) | (qd << 48)
    };
    let b0 = base(hash[0]);
    let b1 = base(hash[1]);
    let b2 = base(hash[2]);
    let b3 = base(hash[3]);
    let mut keys = [0u64; 81];
    let mut idx = 0;
    for d0 in -1..=1i32 {
        for d1 in -1..=1i32 {
            for d2 in -1..=1i32 {
                for d3 in -1..=1i32 {
                    keys[idx] = pack(b0 + d0, b1 + d1, b2 + d2, b3 + d3);
                    idx += 1;
                }
            }
        }
    }
    keys
}

pub fn get_quad_hashes(
    catalog: &Catalog,
    ra_deg: f64,
    dec_deg: f64,
    radius_deg: f64,
    max_stars: usize,
    epoch: Option<JulianDate>,
    min_spine_ratio: Option<f64>,
) -> Vec<Quad> {
    let params = ConeSearchParams {
        ra_deg,
        dec_deg,
        radius_deg,
        max_mag: None,
        max_results: None,
        epoch,
    };
    let mut results = cone_search(catalog, &params);
    results.sort_by(|a, b| {
        a.star
            .mag
            .partial_cmp(&b.star.mag)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    results.truncate(max_stars);

    let stars: Vec<QuadStar> = results
        .iter()
        .filter_map(|r| {
            tan_project_star(r.ra_deg, r.dec_deg, ra_deg, dec_deg).map(|(x, y)| QuadStar {
                source_id: r.star.source_id,
                x,
                y,
                mag: r.star.mag,
            })
        })
        .collect();

    let field_radius_rad = radius_deg * DEG_TO_RAD;
    let ratio = min_spine_ratio.unwrap_or(0.05);
    let min_spine = field_radius_rad * ratio;
    let min_spine_sq = min_spine * min_spine;

    neighbor_quads(&stars, 10, min_spine_sq)
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

    #[test]
    fn test_combinations_4_count() {
        assert_eq!(combinations_4(4).count(), 1);
        assert_eq!(combinations_4(5).count(), 5);
        assert_eq!(combinations_4(10).count(), 210);
        assert_eq!(combinations_4(50).count(), 230300);
        assert_eq!(combinations_4(3).count(), 0);
        assert_eq!(combinations_4(0).count(), 0);
    }

    #[test]
    fn test_combinations_4_values() {
        let combos: Vec<_> = combinations_4(5).collect();
        assert_eq!(combos[0], [0, 1, 2, 3]);
        assert_eq!(combos[1], [0, 1, 2, 4]);
        assert_eq!(combos[2], [0, 1, 3, 4]);
        assert_eq!(combos[3], [0, 2, 3, 4]);
        assert_eq!(combos[4], [1, 2, 3, 4]);
    }

    #[test]
    fn test_find_spine_horizontal() {
        let stars = make_stars(&[
            (0.0, 0.0, 5.0),
            (10.0, 0.0, 5.0),
            (3.0, 1.0, 5.0),
            (7.0, -1.0, 5.0),
        ]);
        let (a, b, _, _) = find_spine(&stars, [0, 1, 2, 3]);
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
        let (a, b, _, _) = find_spine(&stars, [0, 1, 2, 3]);
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
        assert!((hash[0] - 0.3).abs() < 1e-12);
        assert!((hash[1] - 0.4).abs() < 1e-12);
        assert!((hash[2] - 0.7).abs() < 1e-12);
        assert!((hash[3] - (-0.2)).abs() < 1e-12);
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
        assert!((hash[0] - 0.3).abs() < 1e-12);
        assert!((hash[1] - 0.4).abs() < 1e-12);
        assert!((hash[2] - 0.7).abs() < 1e-12);
        assert!((hash[3] - (-0.2)).abs() < 1e-12);
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

        let hash = normalize_quad(&stars, 0, 1, 2, 3).unwrap();
        assert!((hash[0] - cx).abs() < 1e-10);
        assert!((hash[1] - cy).abs() < 1e-10);
        assert!((hash[2] - dx).abs() < 1e-10);
        assert!((hash[3] - dy).abs() < 1e-10);
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
        let result = tan_project_star(45.0, 30.0, 45.0, 30.0);
        let (x, y) = result.unwrap();
        assert!(x.abs() < 1e-15);
        assert!(y.abs() < 1e-15);
    }

    #[test]
    fn test_tan_project_star_behind_returns_none() {
        let result = tan_project_star(0.0, 0.0, 180.0, 0.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_orient_spine_determinism() {
        let s1 = QuadStar { source_id: 0, x: 5.0, y: 3.0, mag: 5.0 };
        let s2 = QuadStar { source_id: 1, x: 1.0, y: 7.0, mag: 5.0 };
        let (a1, b1) = orient_spine(&s1, 0, &s2, 1);
        let (a2, b2) = orient_spine(&s2, 1, &s1, 0);
        assert_eq!(a1, a2);
        assert_eq!(b1, b2);
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
        let mut hashes: Vec<[usize; 4]> = quads.iter().map(|q| {
            let mut s = q.star_indices;
            s.sort();
            s
        }).collect();
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
}
