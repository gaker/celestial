//! HEALPix utilities for cone search queries.
//!
//! Provides conversion between sky coordinates and HEALPix pixel indices,
//! as well as disc/cone query support for efficient spatial searches.

use celestial_core::constants::{PI, RAD_TO_DEG, TWOPI};
use celestial_core::{math::vincenty_angular_separation, Angle};
use std::collections::HashSet;

/// Convert (RA, Dec) in degrees to HEALPix nested pixel index.
///
/// Implements the Gorski et al. (2005) algorithm for the nested scheme.
///
/// # Arguments
/// * `order` - HEALPix order (nside = 2^order)
/// * `ra_deg` - Right ascension in degrees
/// * `dec_deg` - Declination in degrees
///
/// # Returns
/// Nested pixel index in range [0, 12*nside^2)
pub fn ang2pix_nest(order: u32, ra_deg: f64, dec_deg: f64) -> u64 {
    let ra = Angle::from_degrees(ra_deg);
    let dec = Angle::from_degrees(dec_deg);
    let phi = ra.radians();
    let z = dec.sin();
    let nside = 1u64 << order;
    let (face, ix, iy) = compute_face_and_position(phi, z, nside);
    let ipix_in_face = xy2pix_nest(ix, iy, order);
    face as u64 * nside * nside + ipix_in_face
}

/// Query all HEALPix pixels that overlap a cone/disc on the sphere.
///
/// Returns a conservative set of pixels - may include some pixels that
/// don't actually overlap the cone, but will never miss pixels that do.
///
/// Uses a grid-based approach: samples points within the search cone and
/// collects all unique pixel indices that those points fall into.
///
/// # Arguments
/// * `nside` - HEALPix nside parameter
/// * `ra_deg` - Cone center right ascension in degrees
/// * `dec_deg` - Cone center declination in degrees
/// * `radius_deg` - Cone radius in degrees
///
/// # Returns
/// Vector of nested pixel indices that overlap the cone
pub(crate) fn query_disc_nest(nside: u64, ra_deg: f64, dec_deg: f64, radius_deg: f64) -> Vec<u64> {
    let order = nside.trailing_zeros();

    // Pixel size in degrees (approximate) - for order 8, ~0.22 degrees
    let pixel_size_deg = 58.6 / nside as f64; // 58.6 = sqrt(4*pi*(180/pi)^2 / 12) / 1
    let step = pixel_size_deg * 0.5; // Half-pixel resolution for safety

    let mut pixels = HashSet::new();

    // Declination range — pad by one pixel to catch pixels straddling the boundary
    let dec_min = (dec_deg - radius_deg - pixel_size_deg).max(-90.0);
    let dec_max = (dec_deg + radius_deg + pixel_size_deg).min(90.0);

    // Step through declination
    let mut dec = dec_min;
    while dec <= dec_max {
        // RA range expands near poles due to convergence
        let cos_dec = libm::cos(dec * PI / 180.0).max(0.01);
        let ra_step = step / cos_dec;

        // For very high declinations, we need full RA coverage
        let ra_range = if libm::fabs(dec) > 89.0 {
            360.0
        } else {
            (radius_deg / cos_dec).min(180.0) * 2.0
        };

        let ra_min = ra_deg - ra_range / 2.0;
        let ra_max = ra_deg + ra_range / 2.0;

        let mut ra = ra_min;
        while ra <= ra_max {
            // Normalize RA to [0, 360)
            let ra_norm = ((ra % 360.0) + 360.0) % 360.0;

            // Check if this point is actually within the search radius
            let dist = angular_separation_deg(ra_deg, dec_deg, ra_norm, dec);
            if dist <= radius_deg + pixel_size_deg {
                // Include margin for pixel extent
                let pixel = ang2pix_nest(order, ra_norm, dec);
                pixels.insert(pixel);
            }

            ra += ra_step;
        }

        dec += step;
    }

    pixels.into_iter().collect()
}

/// Compute angular distance between two points on the sphere using Vincenty formula.
///
/// Accurate at all angular separations.
///
/// # Arguments
/// * `ra1_deg`, `dec1_deg` - First point in degrees
/// * `ra2_deg`, `dec2_deg` - Second point in degrees
///
/// # Returns
/// Angular distance in degrees
pub(crate) fn angular_separation_deg(
    ra1_deg: f64,
    dec1_deg: f64,
    ra2_deg: f64,
    dec2_deg: f64,
) -> f64 {
    let dec1 = Angle::from_degrees(dec1_deg);
    let dec2 = Angle::from_degrees(dec2_deg);
    let delta_lon = Angle::from_degrees(ra2_deg - ra1_deg).radians();

    let (d1_sin, d1_cos) = dec1.sin_cos();
    let (d2_sin, d2_cos) = dec2.sin_cos();

    let sep_rad = vincenty_angular_separation(d1_sin, d1_cos, d2_sin, d2_cos, delta_lon);
    sep_rad * RAD_TO_DEG
}

/// Determine which of the 12 HEALPix base faces contains the point,
/// and compute the (ix, iy) position within that face.
fn compute_face_and_position(phi: f64, z: f64, nside: u64) -> (u32, u64, u64) {
    let z_abs = libm::fabs(z);
    let tt = phi_to_tt(phi);
    if z_abs <= 2.0 / 3.0 {
        compute_equatorial_face(tt, z, nside)
    } else {
        compute_polar_face(tt, z, z_abs, nside)
    }
}

/// Convert phi to tt (0..4 range for the 4 quadrants).
fn phi_to_tt(phi: f64) -> f64 {
    let phi_norm = if phi < 0.0 { phi + TWOPI } else { phi };
    phi_norm * 2.0 / PI
}

/// Compute face and position for equatorial belt (-2/3 <= z <= 2/3).
fn compute_equatorial_face(tt: f64, z: f64, nside: u64) -> (u32, u64, u64) {
    let temp1 = nside as f64 * (0.5 + tt);
    let temp2 = nside as f64 * z * 0.75;
    let jp = (temp1 - temp2) as i64;
    let jm = (temp1 + temp2) as i64;
    let nside_i = nside as i64;
    let ifp = jp / nside_i;
    let ifm = jm / nside_i;
    let face = compute_equatorial_face_number(ifp, ifm);
    let ix = jm - (face as i64 % 4) * nside_i;
    let iy = nside_i - 1 - (jp - (face as i64 / 4) * nside_i);
    (face, ix as u64, iy as u64)
}

fn compute_equatorial_face_number(ifp: i64, ifm: i64) -> u32 {
    match (ifp, ifm) {
        (4, _) => ((ifm + 4) % 4) as u32,
        (_, 4) => ((ifp + 4) % 4 + 4) as u32,
        _ if ifp == ifm => (ifp + 4) as u32,
        _ if ifp < ifm => ifp as u32,
        _ => (ifm + 8) as u32,
    }
}

/// Compute face and position for polar caps (|z| > 2/3).
fn compute_polar_face(tt: f64, z: f64, z_abs: f64, nside: u64) -> (u32, u64, u64) {
    let tp = tt - libm::floor(tt);
    let tmp = nside as f64 * libm::sqrt(3.0 * (1.0 - z_abs));
    let jp = (tp * tmp) as i64;
    let jm = ((1.0 - tp) * tmp) as i64;
    let jp = jp.min(nside as i64 - 1);
    let jm = jm.min(nside as i64 - 1);
    let ntt = libm::floor(tt) as u32;
    let face_offset = if z > 0.0 { 0 } else { 8 };
    let face = (ntt % 4) + face_offset;
    let (ix, iy) = if z > 0.0 {
        (nside as i64 - jm - 1, nside as i64 - jp - 1)
    } else {
        (jp, jm)
    };
    (face, ix as u64, iy as u64)
}

/// Convert (ix, iy) to nested pixel index within a base face using Z-order curve.
fn xy2pix_nest(ix: u64, iy: u64, order: u32) -> u64 {
    let mut result: u64 = 0;
    for i in 0..order {
        let bit_x = (ix >> i) & 1;
        let bit_y = (iy >> i) & 1;
        result |= (bit_x << (2 * i)) | (bit_y << (2 * i + 1));
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xy2pix_nest() {
        assert_eq!(xy2pix_nest(0, 0, 2), 0);
        assert_eq!(xy2pix_nest(1, 0, 2), 1);
        assert_eq!(xy2pix_nest(0, 1, 2), 2);
        assert_eq!(xy2pix_nest(1, 1, 2), 3);
    }

    #[test]
    fn test_ang2pix_nest_poles() {
        let north_pole = ang2pix_nest(0, 0.0, 90.0);
        assert!(north_pole < 12);
        let south_pole = ang2pix_nest(0, 0.0, -90.0);
        assert!(south_pole < 12);
    }

    #[test]
    fn test_ang2pix_nest_equator() {
        let pixel = ang2pix_nest(0, 0.0, 0.0);
        assert!(pixel < 12);
    }

    #[test]
    fn test_ang2pix_nest_order8_bounds() {
        let nside = 1u64 << 8;
        let npix = 12 * nside * nside;
        for ra in [0.0, 90.0, 180.0, 270.0] {
            for dec in [-89.0, -45.0, 0.0, 45.0, 89.0] {
                let pixel = ang2pix_nest(8, ra, dec);
                assert!(
                    pixel < npix,
                    "pixel {} >= npix {} for ({}, {})",
                    pixel,
                    npix,
                    ra,
                    dec
                );
            }
        }
    }

    #[test]
    fn test_angular_separation_deg() {
        // Same point
        assert!((angular_separation_deg(0.0, 0.0, 0.0, 0.0) - 0.0).abs() < 1e-10);

        // 90 degrees apart on equator
        let dist = angular_separation_deg(0.0, 0.0, 90.0, 0.0);
        assert!((dist - 90.0).abs() < 1e-10);

        // Pole to pole
        let dist = angular_separation_deg(0.0, 90.0, 0.0, -90.0);
        assert!((dist - 180.0).abs() < 1e-10);

        // Small separation
        let dist = angular_separation_deg(0.0, 0.0, 0.1, 0.1);
        assert!(dist > 0.14 && dist < 0.15);
    }

    #[test]
    fn test_query_disc_nest_basic() {
        let nside = 16u64;
        let ra = 0.0;
        let dec = 0.0;
        let radius = 10.0;

        let pixels = query_disc_nest(nside, ra, dec, radius);

        // Should return some pixels
        assert!(!pixels.is_empty());

        // All pixels should be in valid range
        let npix = 12 * nside * nside;
        for &pix in &pixels {
            assert!(pix < npix);
        }

        // Center pixel should be included
        let order = nside.trailing_zeros();
        let center_pixel = ang2pix_nest(order, ra, dec);
        assert!(pixels.contains(&center_pixel));
    }

    #[test]
    fn test_query_disc_nest_pole() {
        let nside = 16u64;
        let ra = 0.0;
        let dec = 90.0;
        let radius = 5.0;

        let pixels = query_disc_nest(nside, ra, dec, radius);

        // Should return pixels around north pole
        assert!(!pixels.is_empty());

        // Center pixel should be included
        let order = nside.trailing_zeros();
        let center_pixel = ang2pix_nest(order, ra, dec);
        assert!(pixels.contains(&center_pixel));
    }
}
