//! The pairs a hypothesis is fitted to.

use super::affine::{AffineFit, LeftOut, Point};
use super::field::Placement;
use super::project;
use crate::match_field::StarPair;

// The fitted pairs that have a tangent-plane projection, with a lookup from a
// detection's exact pixel to its pair.
pub(super) struct FittedPairs {
    // Sky position in degrees, tangent-plane position and detected pixel of each pair.
    points: Vec<((f64, f64), Point, Point)>,
    // Pixel of each pair by bit pattern, sorted, with the pair's index.
    by_pixel: Vec<((u64, u64), usize)>,
}

impl FittedPairs {
    pub(super) fn new(pairs: &[StarPair], tangent: (f64, f64)) -> Self {
        let points: Vec<((f64, f64), Point, Point)> = pairs
            .iter()
            .filter_map(|p| {
                let sky = (p.ra_deg, p.dec_deg);
                Some((sky, project(sky.0, sky.1, tangent)?, [p.px_x, p.px_y]))
            })
            .collect();
        let mut by_pixel: Vec<((u64, u64), usize)> = points
            .iter()
            .enumerate()
            .map(|(i, &(_, _, px))| (pixel_key(px), i))
            .collect();
        by_pixel.sort_unstable();
        Self { points, by_pixel }
    }

    // The affine map from tangent-plane positions to pixels fitted to the pairs.
    pub(super) fn affine_fit(&self) -> Option<AffineFit> {
        let points: Vec<(Point, Point)> = self.points.iter().map(|&(_, q, p)| (q, p)).collect();
        AffineFit::new(&points)
    }

    // The index of the pair whose detection is at `x`.
    pub(super) fn find(&self, x: Point) -> Option<usize> {
        let key = pixel_key(x);
        let at = self.by_pixel.binary_search_by(|e| e.0.cmp(&key)).ok()?;
        Some(self.by_pixel[at].1)
    }

    // `model` without pair `index`, whose residual is its detected pixel minus
    // `placement`'s pixel for it. `None` when `placement` has no pixel for the pair or
    // the other pairs do not determine a model.
    pub(super) fn left_out(
        &self,
        index: usize,
        model: &AffineFit,
        placement: Placement<'_>,
    ) -> Option<LeftOut> {
        let (sky, q, p) = self.points[index];
        let at = placement.pixel(model, sky, q)?;
        model.left_out(q, [p[0] - at[0], p[1] - at[1]])
    }
}

// Pairs carry the exact detected coordinates, so bit equality identifies a detection.
fn pixel_key(p: Point) -> (u64, u64) {
    (p[0].to_bits(), p[1].to_bits())
}
