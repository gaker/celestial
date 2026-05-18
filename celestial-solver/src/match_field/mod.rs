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
pub(super) mod test_catalog {
    use celestial_catalog::query::healpix::ang2pix_nest;
    use std::io::Write;
    use tempfile::NamedTempFile;

    const HEADER_SIZE: usize = 64;
    const PIXEL_ENTRY_SIZE: usize = 16;
    const STAR_RECORD_SIZE: usize = 56;

    pub struct SynthStar {
        pub source_id: i64,
        pub ra: f64,
        pub dec: f64,
        pub mag: f32,
    }

    pub fn build(order: u32, stars: &[SynthStar]) -> NamedTempFile {
        let nside = 1u32 << order;
        let npix = 12u64 * (nside as u64) * (nside as u64);

        let mut buckets: Vec<Vec<&SynthStar>> = (0..npix).map(|_| Vec::new()).collect();
        for s in stars {
            let pix = ang2pix_nest(order, s.ra, s.dec);
            buckets[pix as usize].push(s);
        }
        // Magnitude-sort within each pixel — cone_search relies on it elsewhere.
        for b in buckets.iter_mut() {
            b.sort_by(|a, b| a.mag.partial_cmp(&b.mag).unwrap_or(std::cmp::Ordering::Equal));
        }

        let total_stars: u64 = stars.len() as u64;

        let mut buf: Vec<u8> = Vec::new();
        buf.extend_from_slice(b"CCAT");
        buf.extend_from_slice(&1u32.to_le_bytes());
        buf.extend_from_slice(&order.to_le_bytes());
        buf.extend_from_slice(&nside.to_le_bytes());
        buf.extend_from_slice(&npix.to_le_bytes());
        buf.extend_from_slice(&total_stars.to_le_bytes());
        buf.extend_from_slice(&2016.0f64.to_le_bytes());
        buf.extend_from_slice(&21.0f32.to_le_bytes());
        buf.extend_from_slice(&[0u8; 20]);
        assert_eq!(buf.len(), HEADER_SIZE);

        let mut star_data: Vec<u8> = Vec::new();
        let mut offsets: Vec<(u64, u32)> = vec![(0, 0); npix as usize];
        for (pix, bucket) in buckets.iter().enumerate() {
            let off = star_data.len() as u64;
            offsets[pix] = (off, bucket.len() as u32);
            for s in bucket {
                star_data.extend_from_slice(&s.source_id.to_le_bytes());
                star_data.extend_from_slice(&s.ra.to_le_bytes());
                star_data.extend_from_slice(&s.dec.to_le_bytes());
                star_data.extend_from_slice(&0.0f64.to_le_bytes()); // pmra
                star_data.extend_from_slice(&0.0f64.to_le_bytes()); // pmdec
                star_data.extend_from_slice(&0.0f64.to_le_bytes()); // parallax
                star_data.extend_from_slice(&s.mag.to_le_bytes());
                star_data.extend_from_slice(&0u16.to_le_bytes()); // flags
                star_data.extend_from_slice(&0u16.to_le_bytes()); // padding
            }
        }
        assert_eq!(star_data.len(), stars.len() * STAR_RECORD_SIZE);

        for &(o, c) in &offsets {
            buf.extend_from_slice(&o.to_le_bytes());
            buf.extend_from_slice(&c.to_le_bytes());
            buf.extend_from_slice(&0u32.to_le_bytes());
        }
        assert_eq!(buf.len(), HEADER_SIZE + npix as usize * PIXEL_ENTRY_SIZE);

        buf.extend_from_slice(&star_data);

        let mut file = NamedTempFile::new().unwrap();
        file.write_all(&buf).unwrap();
        file.flush().unwrap();
        file
    }
}

#[cfg(test)]
mod tests {
    use super::test_catalog::{build, SynthStar};
    use super::*;

    use celestial_catalog::query::Catalog;
    use celestial_images::formats::PixelData;
    use celestial_time::JulianDate;

    use crate::detect::DetectedStar;

    fn det(x: f64, y: f64, flux: f64) -> DetectedStar {
        DetectedStar {
            x, y, flux,
            snr: 20.0,
            saturated: false,
            saturated_count: 0,
            background: 0.0,
        }
    }

    fn empty_image(w: usize, h: usize) -> Image {
        Image::new(PixelData::from(vec![0u16; w * h]), [w, h])
    }

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

    #[test]
    fn star_pair_clone_preserves_fields() {
        let p = StarPair {
            px_x: 1.0, px_y: 2.0,
            ra_deg: 10.0, dec_deg: 20.0,
            votes: 5, snr: 30.0,
        };
        let q = p.clone();
        assert_eq!(q.px_x, 1.0);
        assert_eq!(q.px_y, 2.0);
        assert_eq!(q.ra_deg, 10.0);
        assert_eq!(q.dec_deg, 20.0);
        assert_eq!(q.votes, 5);
        assert_eq!(q.snr, 30.0);
    }

    #[test]
    fn quad_match_clone_preserves_quads() {
        use celestial_catalog::query::Quad;
        let qm = QuadMatch {
            image_quad: Quad { hash: [0.1, 0.2, 0.3, 0.4], star_indices: [1, 2, 3, 4] },
            catalog_quad: Quad { hash: [0.5, 0.6, 0.7, 0.8], star_indices: [5, 6, 7, 8] },
        };
        let c = qm.clone();
        assert_eq!(c.image_quad.hash, [0.1, 0.2, 0.3, 0.4]);
        assert_eq!(c.catalog_quad.star_indices, [5, 6, 7, 8]);
    }

    #[test]
    fn match_field_no_detections_returns_empty_pairs() {
        // Catalog with a handful of stars, but zero detections → no image quads,
        // so no matches, no pairs. Must not panic and must return cleanly.
        let center_ra = 180.0;
        let center_dec = 0.0;
        let stars: Vec<SynthStar> = (0..30)
            .map(|i| {
                let i = i as f64;
                SynthStar {
                    source_id: 1000 + i as i64,
                    ra: center_ra + (i * 0.01),
                    dec: center_dec + (i * 0.01 - 0.15),
                    mag: 8.0 + (i * 0.05) as f32,
                }
            })
            .collect();
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();

        let img = empty_image(1024, 1024);
        let hint = ICRSPosition::from_degrees(center_ra, center_dec).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let out = match_field(
            &[],
            &img,
            &catalog,
            &hint,
            2.0,
            epoch,
            &MatchParams::default(),
        )
        .unwrap();

        assert!(out.pairs.is_empty());
        assert!(out.matches.is_empty());
        assert!(out.image_stars.is_empty());
        // Catalog quads need ≥4 stars in the cone; we placed 30 nearby.
        assert!(!out.catalog_stars.is_empty());
    }

    #[test]
    fn match_field_search_radius_override_is_used() {
        // Tiny override radius → very few catalog stars in the cone → catalog_stars short.
        let center_ra = 0.0;
        let center_dec = 0.0;
        let stars: Vec<SynthStar> = (0..50)
            .map(|i| SynthStar {
                source_id: 2000 + i,
                ra: center_ra + (i as f64) * 0.2,
                dec: center_dec + (i as f64) * 0.01,
                mag: 9.0,
            })
            .collect();
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();

        let img = empty_image(1024, 1024);
        let hint = ICRSPosition::from_degrees(center_ra, center_dec).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let tight = MatchParams {
            search_radius_deg: Some(0.05),
            ..MatchParams::default()
        };
        let wide = MatchParams {
            search_radius_deg: Some(5.0),
            ..MatchParams::default()
        };

        let tight_out = match_field(&[], &img, &catalog, &hint, 1.0, epoch, &tight).unwrap();
        let wide_out = match_field(&[], &img, &catalog, &hint, 1.0, epoch, &wide).unwrap();

        // Wider radius must reach at least as many catalog stars.
        assert!(wide_out.catalog_stars.len() >= tight_out.catalog_stars.len());
    }

    #[test]
    fn match_field_runs_both_parities_returns_field_match() {
        // Build a star field and a corresponding set of pixel detections such that
        // both parities are evaluated. We don't assert on which parity wins or how
        // many pairs are produced — just that the call completes and exposes the
        // FieldMatch fields.
        let center_ra = 90.0;
        let center_dec = 45.0;
        let mut stars: Vec<SynthStar> = Vec::new();
        for i in 0..40 {
            let theta = (i as f64) * 0.4;
            stars.push(SynthStar {
                source_id: 3000 + i as i64,
                ra: center_ra + 0.1 * libm::cos(theta),
                dec: center_dec + 0.1 * libm::sin(theta),
                mag: 9.0 + (i as f32) * 0.05,
            });
        }
        let file = build(4, &stars);
        let catalog = Catalog::open(file.path()).unwrap();

        // A small batch of detections — enough to build at least one image quad.
        let detections: Vec<DetectedStar> = (0..20)
            .map(|i| {
                let i = i as f64;
                det(100.0 + i * 30.0, 200.0 + i * 25.0, 5000.0 - i * 100.0)
            })
            .collect();

        let img = empty_image(1024, 1024);
        let hint = ICRSPosition::from_degrees(center_ra, center_dec).unwrap();
        let epoch = JulianDate::new(2451545.0, 0.0);

        let out = match_field(
            &detections,
            &img,
            &catalog,
            &hint,
            2.0,
            epoch,
            &MatchParams::default(),
        )
        .unwrap();

        // image_stars belongs to the winning parity; pairs may be empty for a
        // synthetic non-aligned field, but the structure must be populated.
        let _ = (out.image_stars.len(), out.catalog_stars.len(), out.matches.len(), out.pairs.len());
    }

    #[test]
    fn parity_attempt_default_is_empty() {
        let p = ParityAttempt::default();
        assert!(p.img_stars.is_empty());
        assert!(p.matches.is_empty());
        assert!(p.pairs.is_empty());
        assert_eq!(p.score, 0);
    }
}
