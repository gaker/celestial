//! Star-quad hash matching against a catalog.
//!
//! Builds four-star quads from detections and from a cone-searched catalog region,
//! then matches them by geometric hash. The output is a set of individual [`StarPair`]s
//! (image pixel ↔ catalog RA/Dec), aggregated across all matched quads with vote
//! counts.

mod aggregate;
mod hash_index;
mod quad_stars;
mod verify;

use anyhow::Result;
use celestial_catalog::query::{neighbor_quads, Quad};
use celestial_coords::ICRSPosition;
use celestial_core::constants::PI;
use celestial_images::formats::Image;
use celestial_time::JulianDate;
use rayon::prelude::*;

use crate::detect::DetectedStar;

use aggregate::extract_pairs;
use hash_index::{build_hash_index, find_matches};
use quad_stars::{catalog_cone_search, search_radius, stars_to_quad_stars};
use verify::verify_pairs;

/// Tuning for [`match_field`].
///
/// `max_stars` caps the number of brightest detections considered — more costs cubic
/// matching time; fewer risks missing the right solution in sparse fields.
/// `search_radius_deg` overrides the cone search radius (default: auto-computed from
/// image dimensions and plate scale).
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::match_field::MatchParams;
///
/// let params = MatchParams {
///     max_stars: 200,
///     ..MatchParams::default()
/// };
/// ```
#[derive(Debug, Clone)]
pub struct MatchParams {
    /// Maximum brightest detections to use when building image quads.
    pub max_stars: usize,
    /// How many nearest-neighbor stars to consider when building quads per anchor.
    pub k_neighbors: usize,
    /// Catalog cone search radius override, in degrees. `None` auto-computes.
    pub search_radius_deg: Option<f64>,
}

impl Default for MatchParams {
    fn default() -> Self {
        Self {
            max_stars: 100,
            k_neighbors: 10,
            search_radius_deg: None,
        }
    }
}

/// A geometric quad match: one image-side quad paired with a catalog-side quad.
///
/// Produced by [`match_field`]. The two quads have similar four-star topology (edge
/// ratios and interior angles) under the current hashing scheme.
#[derive(Debug, Clone)]
pub struct QuadMatch {
    /// Quad in image pixel space.
    pub image_quad: Quad,
    /// Quad in catalog sky space.
    pub catalog_quad: Quad,
}

/// A single matched star pair: detected pixel ↔ catalog sky position.
///
/// Aggregated across all [`QuadMatch`]es. Multiple quads touching the same pair bump
/// its `votes` count; the solver weights pairs by vote count and SNR.
///
/// # Examples
///
/// ```rust,ignore
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve(&img, &catalog).run()?;
/// let top = &result.pairs[0];
/// println!("pixel ({:.1}, {:.1}) = sky ({:.4}°, {:.4}°)", top.px_x, top.px_y, top.ra_deg, top.dec_deg);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct StarPair {
    /// Detected pixel x.
    pub px_x: f64,
    /// Detected pixel y.
    pub px_y: f64,
    /// Catalog RA in degrees.
    pub ra_deg: f64,
    /// Catalog Dec in degrees.
    pub dec_deg: f64,
    /// Number of quad matches that contributed this pair.
    pub votes: usize,
    /// SNR from the detection.
    pub snr: f64,
}

/// Output of [`match_field`].
///
/// Exposes the raw quad-star sets and quad matches alongside the aggregated
/// [`StarPair`]s. Pairs are what the WCS fitter consumes; the other fields are
/// exposed for debug rendering and analysis.
#[derive(Debug)]
pub struct FieldMatch {
    /// Image-side quad-star set (brightest detections with quad-relevant bookkeeping).
    pub image_stars: Vec<celestial_catalog::query::QuadStar>,
    /// Catalog-side quad-star set (cone-searched stars with quad-relevant bookkeeping).
    pub catalog_stars: Vec<celestial_catalog::query::QuadStar>,
    /// Every quad match found.
    pub matches: Vec<QuadMatch>,
    /// Aggregated pixel↔sky pairs derived from the matched quads.
    pub pairs: Vec<StarPair>,
}

/// Matches detected stars against a catalog via four-star quad hashing.
///
/// Builds quads from the brightest `params.max_stars` detections, cone-searches the
/// catalog around `hint` with a radius derived from the image geometry, builds
/// catalog quads, and matches by geometric hash. Tries both parity orientations and
/// returns the better result.
///
/// # Errors
///
/// Returns an error when the cone search fails or no quad matches are found.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::match_field::{match_field, MatchParams};
///
/// # let stars = vec![];
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
/// # let hint = celestial_coords::ICRSPosition::from_degrees(180.0, 30.0)?;
/// # let epoch = celestial_time::JulianDate::new(2_400_000.5, 0.0);
/// let field = match_field(&stars, &img, &catalog, &hint, 1.5, epoch, &MatchParams::default())?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn match_field(
    stars: &[DetectedStar],
    image: &Image,
    catalog: &celestial_catalog::query::Catalog,
    hint: &ICRSPosition,
    scale_arcsec: f64,
    epoch: JulianDate,
    params: &MatchParams,
) -> Result<FieldMatch> {
    let w = image.width();
    let h = image.height();
    let ra = hint.ra().degrees();
    let dec = hint.dec().degrees();

    let radius = search_radius(w, h, scale_arcsec, params.search_radius_deg);
    let field_rad = radius * PI / 180.0;
    let min_spine = field_rad * 0.05;
    let min_spine_sq = min_spine * min_spine;

    let (catalog_stars, catalog_results) = catalog_cone_search(
        catalog, ra, dec, radius, params.max_stars, epoch,
    );
    let catalog_quads = neighbor_quads(&catalog_stars, 10, min_spine_sq);
    let index = build_hash_index(&catalog_quads);

    log::debug!("match_field: radius={:.3}\u{00b0}, scale={:.3}\"/px", radius, scale_arcsec);
    log::debug!("match_field: {} catalog stars, {} catalog quads", catalog_stars.len(), catalog_quads.len());

    let results: Vec<ParityAttempt> = [false, true]
        .par_iter()
        .map(|&parity| {
            let img_stars = stars_to_quad_stars(stars, w, h, scale_arcsec, params.max_stars, parity);
            let img_quads = neighbor_quads(&img_stars, params.k_neighbors, min_spine_sq);
            let m = find_matches(&img_quads, &catalog_quads, &index);
            let p = extract_pairs(&m, stars, &img_stars, &catalog_stars, &catalog_results);
            let score = verify_pairs(&p, ra, dec);
            let label = if parity { "flipped" } else { "normal" };
            log::debug!("match_field [{label}]: {} img quads, {} raw matches, {} pairs, verified={}",
                img_quads.len(), m.len(), p.len(), score);
            ParityAttempt { img_stars, matches: m, pairs: p, score }
        })
        .collect();

    let mut best_score = 0usize;
    let mut best: Option<ParityAttempt> = None;
    for r in results {
        if r.score > best_score {
            best_score = r.score;
            best = Some(r);
        }
    }

    let ParityAttempt { img_stars: image_stars, matches, pairs, .. } = best.unwrap_or_default();
    log::debug!("match_field: best parity => {} pairs (score {})", pairs.len(), best_score);

    Ok(FieldMatch {
        image_stars,
        catalog_stars,
        matches,
        pairs,
    })
}

#[derive(Default)]
struct ParityAttempt {
    img_stars: Vec<celestial_catalog::query::QuadStar>,
    matches: Vec<QuadMatch>,
    pairs: Vec<StarPair>,
    score: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn match_params_default_values() {
        let p = MatchParams::default();
        assert_eq!(p.max_stars, 100);
        assert_eq!(p.k_neighbors, 10);
        assert!(p.search_radius_deg.is_none());
    }

    #[test]
    fn match_params_clone() {
        let p = MatchParams { max_stars: 250, k_neighbors: 20, search_radius_deg: Some(3.5) };
        let q = p.clone();
        assert_eq!(q.max_stars, 250);
        assert_eq!(q.k_neighbors, 20);
        assert_eq!(q.search_radius_deg, Some(3.5));
    }
}
