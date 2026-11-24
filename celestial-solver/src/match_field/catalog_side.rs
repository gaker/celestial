// The catalog side of the quad search: the cone's quad stars, their quads and the hash
// index of those quads.

use std::collections::HashMap;

use celestial_catalog::query::catalog::Catalog;
use celestial_catalog::query::cone::ConeSearchResult;
use celestial_catalog::query::quad::{neighbor_quads, Quad, QuadStar};
use celestial_core::constants::{ARCSEC_TO_RAD, DEG_TO_RAD};
use celestial_time::julian::JulianDate;

use crate::error::SolveError;

use super::hash_index::build_hash_index;
use super::quad_stars::{
    catalog_cone_search, cone_share, reference_count, search_radius, QuadSelection,
};
use super::uniform::image_cell_px;
use super::MatchParams;

// The catalog side of the search: the cone's quad stars with their cone-search results,
// the quads built from them, the hash index of those quads and the reference stars.
pub(super) struct CatalogSide {
    pub(super) quad_stars: Vec<QuadStar>,
    pub(super) results: Vec<ConeSearchResult>,
    pub(super) quads: Vec<Quad>,
    pub(super) index: HashMap<u64, Vec<usize>>,
    pub(super) references: Vec<(f64, f64)>,
    // Quads whose squared spine is shorter than this are not used, on either side.
    pub(super) min_spine_sq: f64,
}

// The catalog side for a cone about `centre` covering an image of `(w, h)` pixels at
// `scale_arcsec`, with as many quad stars as `n_detections` give the image.
pub(super) fn search_catalog(
    catalog: &Catalog,
    centre: (f64, f64),
    (w, h): (usize, usize),
    scale_arcsec: f64,
    epoch: JulianDate,
    n_detections: usize,
    params: &MatchParams,
) -> Result<CatalogSide, SolveError> {
    let radius = search_radius(w, h, scale_arcsec, params.search_radius_deg);
    let min_spine_sq = min_spine_sq(w, h, scale_arcsec);
    let selection = catalog_selection(n_detections, (w, h), scale_arcsec, radius, params);
    let n_references = reference_count(w, h, scale_arcsec, radius);
    let cone = catalog_cone_search(catalog, centre, radius, &selection, n_references, epoch)?;
    let quads = neighbor_quads(&cone.quad_stars, params.k_neighbors, min_spine_sq);
    let index = build_hash_index(&quads);
    log::debug!("match_field: radius={radius:.3}\u{00b0}, scale={scale_arcsec:.3}\"/px");
    let (n_stars, n_quads) = (cone.quad_stars.len(), quads.len());
    log::debug!("match_field: {n_stars} catalog stars, {n_quads} catalog quads");
    Ok(CatalogSide {
        quad_stars: cone.quad_stars,
        results: cone.results,
        quads,
        index,
        references: cone.references,
        min_spine_sq,
    })
}

// Quads shorter than 5% of the image's own search radius are not used, whatever cone is
// searched. The square of that length on the tangent plane, for an image of `w` x `h`
// pixels at `scale_arcsec`.
fn min_spine_sq(w: usize, h: usize, scale_arcsec: f64) -> f64 {
    let min_spine = search_radius(w, h, scale_arcsec, None) * DEG_TO_RAD * 0.05;
    min_spine * min_spine
}

// The catalog stars for quads sample the cone at the density of the image's quad stars,
// at most `params.max_stars` of the `n_detections`, over the image, through cells of the
// area the image's cells have. The hint gives no rotation, so the image can lie anywhere
// in the cone, and the whole cone is sampled.
fn catalog_selection(
    n_detections: usize,
    (w, h): (usize, usize),
    scale_arcsec: f64,
    radius_deg: f64,
    params: &MatchParams,
) -> QuadSelection {
    let n_image = n_detections.min(params.max_stars);
    let (cell_w, cell_h) = image_cell_px(w, h, n_image, params.k_neighbors);
    QuadSelection {
        count: cone_share(n_image, w, h, scale_arcsec, radius_deg),
        cell_rad: libm::sqrt(cell_w * cell_h) * scale_arcsec * ARCSEC_TO_RAD,
    }
}
