// 3x3 median filter, with the neighbourhood clipped to the image at its edges.

use rayon::prelude::*;

pub(super) fn median_filter_3x3(data: &[f32], w: usize, h: usize) -> Vec<f32> {
    if w < 3 || h < 3 {
        return median_filter_3x3_edge_only(data, w, h);
    }
    let mut out = vec![0.0_f32; w * h];
    out.par_chunks_mut(w).enumerate().for_each(|(y, row)| {
        if y == 0 || y == h - 1 {
            fill_edge_row(data, row, y, w, h);
        } else {
            fill_interior_row(data, row, y, w, h);
        }
    });
    out
}

// Row `y`, which has a row above and below it in the image.
fn fill_interior_row(data: &[f32], row: &mut [f32], y: usize, w: usize, h: usize) {
    let off_t = (y - 1) * w;
    let off_m = y * w;
    let off_b = (y + 1) * w;
    row[0] = edge_median(data, 0, y, w, h);
    for x in 1..w - 1 {
        row[x] = median9_network([
            data[off_t + x - 1],
            data[off_t + x],
            data[off_t + x + 1],
            data[off_m + x - 1],
            data[off_m + x],
            data[off_m + x + 1],
            data[off_b + x - 1],
            data[off_b + x],
            data[off_b + x + 1],
        ]);
    }
    row[w - 1] = edge_median(data, w - 1, y, w, h);
}

#[inline(always)]
fn median9_network(pixels: [f32; 9]) -> f32 {
    let mut buf = pixels;
    let (_, m, _) = buf.select_nth_unstable_by(4, |a, b| {
        a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
    });
    *m
}

fn fill_edge_row(data: &[f32], row: &mut [f32], y: usize, w: usize, h: usize) {
    for (x, cell) in row.iter_mut().enumerate() {
        *cell = edge_median(data, x, y, w, h);
    }
}

fn edge_median(data: &[f32], x: usize, y: usize, w: usize, h: usize) -> f32 {
    let x0 = x.saturating_sub(1);
    let x1 = if x + 1 < w { x + 1 } else { w - 1 };
    let y0 = y.saturating_sub(1);
    let y1 = if y + 1 < h { y + 1 } else { h - 1 };
    let mut buf = [0.0_f32; 9];
    let mut n = 0;
    for ny in y0..=y1 {
        for nx in x0..=x1 {
            buf[n] = data[ny * w + nx];
            n += 1;
        }
    }
    buf[..n].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    buf[n / 2]
}

fn median_filter_3x3_edge_only(data: &[f32], w: usize, h: usize) -> Vec<f32> {
    let mut out = vec![0.0_f32; w * h];
    for y in 0..h {
        for x in 0..w {
            out[y * w + x] = edge_median(data, x, y, w, h);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median9_network_matches_sort_on_all_permutations() {
        fn permute(v: &mut [f32], k: usize, out: &mut Vec<[f32; 9]>) {
            if k == v.len() {
                let mut arr = [0.0_f32; 9];
                arr.copy_from_slice(v);
                out.push(arr);
                return;
            }
            for i in k..v.len() {
                v.swap(k, i);
                permute(v, k + 1, out);
                v.swap(k, i);
            }
        }
        let mut v = [1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let mut perms = Vec::new();
        permute(&mut v, 0, &mut perms);
        assert_eq!(perms.len(), 362880);
        for p in perms {
            let got = median9_network(p);
            assert_eq!(got, 5.0, "permutation {p:?} produced {got}");
        }
    }

    #[test]
    fn median9_network_matches_sort_on_duplicates() {
        let cases: &[([f32; 9], f32)] = &[
            ([1.0; 9], 1.0),
            ([1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0], 2.0),
            ([0.0, 0.0, 0.0, 0.0, 5.0, 9.0, 9.0, 9.0, 9.0], 5.0),
            ([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0], 1.0),
        ];
        for &(p, want) in cases {
            let got = median9_network(p);
            assert_eq!(got, want, "input {p:?} expected {want}, got {got}");
        }
    }

    #[test]
    fn median_filter_3x3_matches_naive() {
        let w = 17;
        let h = 13;
        let data: Vec<f32> = (0..w * h).map(|i| ((i * 37) % 101) as f32).collect();
        let fast = median_filter_3x3(&data, w, h);
        let naive = median_filter_3x3_edge_only(&data, w, h);
        for i in 0..w * h {
            assert_eq!(fast[i], naive[i], "pixel {i}");
        }
    }
}
