//! Catalog file contents for tests, in this crate and in crates that enable the
//! `test-support` feature.

use super::catalog::{pixel_entry_bytes, CatalogHeader, StarRecord};
use super::format::STAR_RECORD_SIZE;

// A star with no proper motion or parallax.
pub fn make_star(source_id: i64, ra: f64, dec: f64, mag: f32, flags: u16) -> StarRecord {
    StarRecord {
        source_id,
        ra,
        dec,
        pmra: 0.0,
        pmdec: 0.0,
        parallax: 0.0,
        mag,
        flags,
        _padding: 0,
    }
}

pub fn header(order: u32, total_stars: u64, epoch: f64, mag_limit: f32) -> CatalogHeader {
    let nside = 1u64 << order;
    CatalogHeader {
        order,
        nside: nside as u32,
        npix: 12 * nside * nside,
        total_stars,
        epoch,
        mag_limit,
    }
}

// The bytes of a catalog file with every pixel's records in pixel order. A pixel listed
// twice keeps the stars of both lists. Records keep the order they are given in.
pub fn catalog_bytes(
    order: u32,
    epoch: f64,
    mag_limit: f32,
    pixel_stars: &[(u64, Vec<StarRecord>)],
) -> Vec<u8> {
    let total_stars: u64 = pixel_stars.iter().map(|(_, s)| s.len() as u64).sum();
    let header = header(order, total_stars, epoch, mag_limit);
    let mut pixels: Vec<Vec<StarRecord>> = vec![Vec::new(); header.npix as usize];
    for (pixel, stars) in pixel_stars {
        pixels[*pixel as usize].extend_from_slice(stars);
    }

    let mut buf = header.to_bytes().to_vec();
    let mut offset = 0u64;
    for stars in &pixels {
        buf.extend_from_slice(&pixel_entry_bytes(offset, stars.len() as u32));
        offset += (stars.len() * STAR_RECORD_SIZE) as u64;
    }
    for star in pixels.iter().flatten() {
        buf.extend_from_slice(&star.to_le_bytes());
    }
    buf
}
