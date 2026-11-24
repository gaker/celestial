//! The detections and catalog stars a hypothesis is verified against, and the
//! references a hypothesis places in the image.

use super::affine::{AffineFit, Point};
use super::{project, MAX_REFERENCES, MAX_SCORED_STARS};
use crate::detect::DetectedStar;
use crate::fit_wcs::solution::WcsSolution;

// The detections and catalog stars a hypothesis is verified against.
pub(crate) struct Field<'a> {
    pub(super) stars: &'a [DetectedStar],
    // Indices of the `MAX_SCORED_STARS` brightest detections in descending flux; ties
    // keep input order.
    pub(super) scored: Vec<usize>,
    pub(super) width: f64,
    pub(super) height: f64,
    // (ra, dec) in degrees of the tangent point every position is projected about.
    pub(super) tangent: (f64, f64),
    // Candidate reference stars, brightest first: sky position in degrees and
    // tangent-plane position.
    pub(super) references: Vec<((f64, f64), Point)>,
}

impl<'a> Field<'a> {
    // `references` are (ra, dec) in degrees, brightest first.
    pub(crate) fn new(
        stars: &'a [DetectedStar],
        (width, height): (usize, usize),
        tangent: (f64, f64),
        references: &[(f64, f64)],
    ) -> Self {
        let mut scored: Vec<usize> = (0..stars.len()).collect();
        scored.sort_by(|&a, &b| stars[b].flux.total_cmp(&stars[a].flux));
        scored.truncate(MAX_SCORED_STARS);
        let references = references
            .iter()
            .filter_map(|&sky| Some((sky, project(sky.0, sky.1, tangent)?)))
            .collect();
        Self {
            stars,
            scored,
            width: width as f64,
            height: height as f64,
            tangent,
            references,
        }
    }

    pub(super) fn area(&self) -> f64 {
        self.width * self.height
    }

    pub(super) fn contains(&self, p: Point) -> bool {
        (0.0..self.width).contains(&p[0]) && (0.0..self.height).contains(&p[1])
    }
}

// A reference star placed in the image by a hypothesis.
pub(super) struct Placed {
    pub(super) sky: (f64, f64),
    pub(super) q: Point,
    pub(super) pixel: Point,
    // `1 + h(q)`: the per-axis variance of a detection about `pixel`, in units of
    // `POSITION_SIGMA_PX^2`.
    pub(super) variance: f64,
}

// What maps a sky position to the pixel a hypothesis predicts for it.
#[derive(Clone, Copy)]
pub(super) enum Placement<'w> {
    // The hypothesis's affine fit.
    Affine,
    // A WCS fitted to the hypothesis's pairs. A position it does not project to a
    // finite pixel is not in the image.
    Wcs(&'w WcsSolution),
}

impl Placement<'_> {
    // The pixel of `sky`, whose tangent-plane position is `q`.
    pub(super) fn pixel(self, model: &AffineFit, sky: (f64, f64), q: Point) -> Option<Point> {
        match self {
            Placement::Affine => Some(model.predict(q)),
            Placement::Wcs(wcs) => {
                let (x, y) = wcs.sky_to_pixel(sky.0, sky.1)?;
                (x.is_finite() && y.is_finite()).then_some([x, y])
            }
        }
    }
}

// The `MAX_REFERENCES` brightest field references that `placement` puts in the image,
// with the variance of `model`'s prediction.
pub(super) fn place_references(
    field: &Field<'_>,
    model: &AffineFit,
    placement: Placement<'_>,
) -> Vec<Placed> {
    let placed = field.references.iter().filter_map(|&(sky, q)| {
        let pixel = placement.pixel(model, sky, q)?;
        let variance = 1.0 + model.leverage(q);
        field.contains(pixel).then_some(Placed {
            sky,
            q,
            pixel,
            variance,
        })
    });
    placed.take(MAX_REFERENCES).collect()
}
