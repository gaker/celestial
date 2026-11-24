// Gnomonic projection of matched pairs onto the tangent plane of a CRVAL.

use celestial_catalog::query::quad::tan_project_star;
use celestial_core::constants::RAD_TO_DEG;

use crate::match_field::StarPair;

pub(super) type ProjectedPair = (f64, f64, f64, f64);

// The pairs that have a gnomonic projection, index-aligned with their tangent-plane rows
// and weights. Pairs more than 90 degrees from the tangent point are dropped from all
// three, so a mask over `rows` also indexes `pairs` and `weights`.
pub(super) struct TangentProjection<'a> {
    pub(super) pairs: Vec<&'a StarPair>,
    pub(super) rows: Vec<ProjectedPair>,
    pub(super) weights: Vec<f64>,
}

pub(super) fn project_pairs_to_tangent<'a>(
    pairs: &'a [StarPair],
    weights: &[f64],
    crval_ra: f64,
    crval_dec: f64,
    crpix1: f64,
    crpix2: f64,
) -> TangentProjection<'a> {
    let mut tp = TangentProjection {
        pairs: Vec::with_capacity(pairs.len()),
        rows: Vec::with_capacity(pairs.len()),
        weights: Vec::with_capacity(pairs.len()),
    };
    for (p, &wt) in pairs.iter().zip(weights.iter()) {
        if let Some((xi, eta)) = tan_project_star(p.ra_deg, p.dec_deg, crval_ra, crval_dec) {
            tp.pairs.push(p);
            tp.rows.push((
                p.px_x - crpix1,
                p.px_y - crpix2,
                xi * RAD_TO_DEG,
                eta * RAD_TO_DEG,
            ));
            tp.weights.push(wt);
        }
    }
    tp
}

pub(super) fn extract_kept_pairs<'a>(
    tp: &TangentProjection<'a>,
    mask: &[bool],
) -> (Vec<&'a StarPair>, Vec<f64>) {
    tp.pairs
        .iter()
        .zip(&tp.weights)
        .zip(mask)
        .filter(|(_, &m)| m)
        .map(|((&p, &w), _)| (p, w))
        .unzip()
}
