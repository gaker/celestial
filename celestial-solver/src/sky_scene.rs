//! Synthetic star fields for end-to-end tests: a catalog file, the detections a camera
//! would see, and a rendered image of them.

use std::io::Write;

use celestial_catalog::query::healpix::ang2pix_nest;
use celestial_catalog::query::quad::{tan_deproject_star, tan_project_star};
use celestial_catalog::query::test_support::{catalog_bytes, make_star};
use celestial_coords::frames::icrs::ICRSPosition;
use celestial_core::constants::{DEG_TO_RAD, RAD_TO_DEG, TWOPI};
use celestial_images::formats::image::Image;
use celestial_images::formats::pixel_data::PixelData;
use celestial_time::julian::JulianDate;
use tempfile::NamedTempFile;

use crate::detect::DetectedStar;
use crate::match_field::StarPair;
use crate::metadata::image_metadata::ImageMetadata;

pub(crate) const WIDTH: usize = 1024;
pub(crate) const HEIGHT: usize = 1024;
pub(crate) const SCALE_ARCSEC: f64 = 2.0;
// The true field centre; the hint is 0.32 degrees away from it.
pub(crate) const CENTER: (f64, f64) = (180.3, 0.1);
const HINT: (f64, f64) = (180.0, 0.0);
const ROTATION_RAD: f64 = 0.4;
const MARGIN_PX: f64 = 30.0;

pub(crate) struct SynthStar {
    pub(crate) source_id: i64,
    pub(crate) ra: f64,
    pub(crate) dec: f64,
    pub(crate) mag: f32,
}

// A catalog file of HEALPix order `order` holding `stars`, brightest first within each
// pixel.
pub(crate) fn build(order: u32, stars: &[SynthStar]) -> NamedTempFile {
    let mut by_mag: Vec<&SynthStar> = stars.iter().collect();
    by_mag.sort_by(|a, b| a.mag.total_cmp(&b.mag));
    let pixel_stars: Vec<(u64, Vec<_>)> = by_mag
        .iter()
        .map(|s| {
            let pixel = ang2pix_nest(order, s.ra, s.dec).unwrap_or_else(|e| {
                panic!("synthetic star ({}, {}) has no pixel: {e}", s.ra, s.dec)
            });
            (pixel, vec![make_star(s.source_id, s.ra, s.dec, s.mag, 0)])
        })
        .collect();
    let mut file = NamedTempFile::new().expect("a temporary file can be created");
    file.write_all(&catalog_bytes(order, 2016.0, 21.0, &pixel_stars))
        .expect("the temporary catalog file is writable");
    file
}

pub(crate) struct Scene {
    pub(crate) catalog: NamedTempFile,
    // Image width and height in pixels.
    pub(crate) size: (usize, usize),
    // (pixel x, pixel y, ra, dec) of every catalog star placed in the image.
    pub(crate) truth: Vec<(f64, f64, f64, f64)>,
    // Pixel positions of every star with no catalog counterpart.
    pub(crate) junk: Vec<(f64, f64)>,
}

// Deterministic uniform deviates in [0, 1) from a 64-bit LCG.
pub(crate) fn uniform(seed: u64) -> impl FnMut() -> f64 {
    let mut s = seed;
    move || {
        s = s
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (s >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn random_catalog(next: &mut impl FnMut() -> f64, n: usize) -> Vec<SynthStar> {
    (0..n)
        .map(|i| SynthStar {
            source_id: i as i64 + 1,
            ra: CENTER.0 + (next() - 0.5) * 0.9,
            dec: CENTER.1 + (next() - 0.5) * 0.9,
            mag: 8.0 + i as f32 * 0.02,
        })
        .collect()
}

// Pixel position of a sky position under a TAN WCS at CENTER with SCALE_ARCSEC,
// rotated by ROTATION_RAD; `flip` mirrors the y axis.
fn sky_to_pixel(ra: f64, dec: f64, flip: bool) -> Option<(f64, f64)> {
    sky_to_pixel_in(ra, dec, flip, (WIDTH, HEIGHT))
}

// As `sky_to_pixel`, for an image of `size` pixels.
fn sky_to_pixel_in(ra: f64, dec: f64, flip: bool, size: (usize, usize)) -> Option<(f64, f64)> {
    let (width, height) = (size.0 as f64, size.1 as f64);
    let (xi, eta) = tan_project_star(ra, dec, CENTER.0, CENTER.1)?;
    let per_px = SCALE_ARCSEC / 3600.0;
    let (u, v) = (xi * RAD_TO_DEG / per_px, eta * RAD_TO_DEG / per_px);
    let (sn, cs) = libm::sincos(ROTATION_RAD);
    let y_sign = if flip { -1.0 } else { 1.0 };
    let x = width / 2.0 + cs * u - sn * v;
    let y = height / 2.0 + y_sign * (sn * u + cs * v);
    let inside =
        (MARGIN_PX..width - MARGIN_PX).contains(&x) && (MARGIN_PX..height - MARGIN_PX).contains(&y);
    inside.then_some((x, y))
}

fn random_pixels(next: &mut impl FnMut() -> f64, n: usize) -> Vec<(f64, f64)> {
    let span_x = WIDTH as f64 - 2.0 * MARGIN_PX;
    let span_y = HEIGHT as f64 - 2.0 * MARGIN_PX;
    (0..n)
        .map(|_| (MARGIN_PX + next() * span_x, MARGIN_PX + next() * span_y))
        .collect()
}

// `n_catalog` catalog stars around CENTER; every one that lands in the image is
// visible, and `n_junk` stars with no catalog counterpart are added.
pub(crate) fn star_field(seed: u64, n_catalog: usize, n_junk: usize, flip: bool) -> Scene {
    let mut next = uniform(seed);
    let stars = random_catalog(&mut next, n_catalog);
    let truth = stars
        .iter()
        .filter_map(|s| sky_to_pixel(s.ra, s.dec, flip).map(|(x, y)| (x, y, s.ra, s.dec)))
        .collect();
    let junk = random_pixels(&mut next, n_junk);
    Scene {
        catalog: build(6, &stars),
        size: (WIDTH, HEIGHT),
        truth,
        junk,
    }
}

// `n_catalog` catalog stars spread uniformly over a disc of `radius_deg` about CENTER,
// fainter with each index, and an image of `size` pixels that shows every catalog star
// inside it, with no junk. The image ranks its stars by catalog magnitude plus a uniform
// error of up to 0.5 magnitudes, as a camera band differs from the catalog's.
pub(crate) fn wide_field(
    seed: u64,
    n_catalog: usize,
    radius_deg: f64,
    size: (usize, usize),
) -> Scene {
    let mut next = uniform(seed);
    let radius = radius_deg * DEG_TO_RAD;
    let stars: Vec<SynthStar> = (0..n_catalog)
        .map(|i| {
            let (r, theta) = (radius * libm::sqrt(next()), TWOPI * next());
            let (ra, dec) = tan_deproject_star(
                r * libm::cos(theta),
                r * libm::sin(theta),
                CENTER.0,
                CENTER.1,
            );
            SynthStar {
                source_id: i as i64 + 1,
                ra,
                dec,
                mag: 8.0 + 6.0 * i as f32 / n_catalog as f32,
            }
        })
        .collect();
    let mut ranked: Vec<(f64, (f64, f64, f64, f64))> = stars
        .iter()
        .filter_map(|s| {
            let (x, y) = sky_to_pixel_in(s.ra, s.dec, false, size)?;
            Some((s.mag as f64 + next() - 0.5, (x, y, s.ra, s.dec)))
        })
        .collect();
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
    let truth = ranked.into_iter().map(|(_, t)| t).collect();
    Scene {
        catalog: build(6, &stars),
        size,
        truth,
        junk: Vec::new(),
    }
}

// A catalog and a set of stars with no relation between them.
pub(crate) fn unrelated_field(seed: u64, n_stars: usize) -> Scene {
    let mut next = uniform(seed);
    let stars = random_catalog(&mut next, 300);
    let junk = random_pixels(&mut next, n_stars);
    Scene {
        catalog: build(6, &stars),
        size: (WIDTH, HEIGHT),
        truth: Vec::new(),
        junk,
    }
}

impl Scene {
    // Detections in descending flux, alternating true and junk stars while both last.
    pub(crate) fn detections(&self) -> Vec<DetectedStar> {
        let truth: Vec<(f64, f64)> = self.truth.iter().map(|&(x, y, _, _)| (x, y)).collect();
        let n = truth.len().max(self.junk.len());
        let alternating = (0..n).flat_map(|i| [truth.get(i), self.junk.get(i)]);
        alternating
            .flatten()
            .copied()
            .enumerate()
            .map(|(i, (x, y))| DetectedStar {
                x,
                y,
                flux: 30_000.0 - i as f64 * 50.0,
                snr: 50.0,
                saturated: false,
                saturated_count: 0,
                background: 100.0,
            })
            .collect()
    }

    // Gaussian stars (sigma 1.5 px) on a background of 100 with uniform noise of
    // amplitude 10, in the flux order of `detections`.
    pub(crate) fn image(&self) -> Image {
        let mut px = vec![100.0f32; WIDTH * HEIGHT];
        for (i, s) in self.detections().iter().enumerate() {
            add_star(&mut px, s.x, s.y, 20_000.0 - i as f64 * 100.0);
        }
        let mut next = uniform(99);
        for v in px.iter_mut() {
            *v += ((next() - 0.5) * 10.0) as f32;
        }
        Image::new(PixelData::from(px), [WIDTH, HEIGHT])
    }

    // Each pair must join a detection to the catalog star placed there: the scene star
    // nearest the detection is a catalog star, and it is the catalog star nearest the
    // pair's sky position. Positions are compared by nearest neighbour because the cone
    // search returns epoch-propagated positions that can differ in the last bit.
    pub(crate) fn assert_true_pairs(&self, pairs: &[StarPair]) {
        assert!(!pairs.is_empty());
        for p in pairs {
            let by_pixel = nearest(self.truth.iter().map(|t| (t.0, t.1)), (p.px_x, p.px_y));
            let by_sky = nearest(self.truth.iter().map(|t| (t.2, t.3)), (p.ra_deg, p.dec_deg));
            let junk = nearest(self.junk.iter().copied(), (p.px_x, p.px_y));
            assert!(
                by_pixel.1 < junk.1,
                "pair at ({}, {}) is a junk star",
                p.px_x,
                p.px_y
            );
            assert_eq!(by_pixel.0, by_sky.0, "pair at ({}, {})", p.px_x, p.px_y);
        }
    }

    pub(crate) fn metadata(&self) -> ImageMetadata {
        ImageMetadata {
            hint: hint(),
            scale_arcsec: SCALE_ARCSEC,
            epoch: epoch(),
            focal_mm: None,
            pixel_um: None,
        }
    }
}

pub(crate) fn hint() -> ICRSPosition {
    ICRSPosition::from_degrees(HINT.0, HINT.1).expect("the hint is a valid sky position")
}

pub(crate) fn epoch() -> JulianDate {
    JulianDate::new(2451545.0, 0.0)
}

fn add_star(px: &mut [f32], cx: f64, cy: f64, peak: f64) {
    let (x0, x1) = (cx as usize - 8, cx as usize + 9);
    let (y0, y1) = (cy as usize - 8, cy as usize + 9);
    for y in y0..y1 {
        for x in x0..x1 {
            let dx = x as f64 - cx;
            let dy = y as f64 - cy;
            let g = libm::exp(-(dx * dx + dy * dy) / (2.0 * 1.5 * 1.5));
            px[y * WIDTH + x] += (peak * g) as f32;
        }
    }
}

// Index of the point nearest `to` and its squared distance; (usize::MAX, inf) for none.
fn nearest(points: impl Iterator<Item = (f64, f64)>, to: (f64, f64)) -> (usize, f64) {
    points
        .map(|(x, y)| (x - to.0) * (x - to.0) + (y - to.1) * (y - to.1))
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .unwrap_or((usize::MAX, f64::INFINITY))
}
