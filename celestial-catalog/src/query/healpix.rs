//! HEALPix utilities for cone search queries.
//!
//! Provides conversion between sky coordinates and HEALPix pixel indices,
//! as well as disc/cone query support for efficient spatial searches.

use celestial_core::angle::Angle;
use celestial_core::constants::RAD_TO_DEG;
use celestial_core::math::angular_separation;

use crate::error::{CatalogError, CatalogResult};

// The HEALPix library's 64-bit NEST scheme stops at order 29: ix and iy need 29 bits each.
const MAX_ORDER: u32 = 29;

/// Convert (RA, Dec) in degrees to HEALPix nested pixel index.
///
/// Implements the Gorski et al. (2005) algorithm for the nested scheme.
///
/// # Arguments
/// * `order` - HEALPix order (nside = 2^order), at most 29
/// * `ra_deg` - Right ascension in degrees, any finite value
/// * `dec_deg` - Declination in degrees, in [-90, 90]
///
/// # Returns
/// Nested pixel index in range [0, 12*nside^2)
///
/// # Errors
/// Returns an error if `order` is above 29, either angle is not finite, or
/// `dec_deg` is outside [-90, 90].
pub fn ang2pix_nest(order: u32, ra_deg: f64, dec_deg: f64) -> CatalogResult<u64> {
    let (nside, _) = nside_and_npix(order)?;
    check_position(ra_deg, dec_deg)?;
    // Dividing by 90 keeps tt, the longitude in quarter turns, below 4 for every RA below 360.
    let tt = reduce_ra_deg(ra_deg) / 90.0;
    let (face, ix, iy) = compute_face_and_position(tt, dec_deg, nside);
    Ok(u64::from(face) * nside * nside + xy2pix_nest(ix, iy, order))
}

// Returns (nside, npix) = (2^order, 12 nside^2), or an error for an order above 29.
pub fn nside_and_npix(order: u32) -> CatalogResult<(u64, u64)> {
    if order > MAX_ORDER {
        return Err(CatalogError::invalid_input(format!(
            "HEALPix order {order} is above the maximum of {MAX_ORDER}"
        )));
    }
    let nside = 1u64 << order;
    Ok((nside, 12 * nside * nside))
}

// Reduces a finite RA into [0, 360) degrees. fmod is exact, so RA + 360 k reduces to the same
// value as RA whenever both are representable.
pub fn reduce_ra_deg(ra_deg: f64) -> f64 {
    let ra = libm::fmod(ra_deg, 360.0);
    let ra = if ra < 0.0 { ra + 360.0 } else { ra };
    // ra + 360 rounds to 360 for ra in (-3e-14, 0). The zero test maps -0.0 to +0.0.
    if ra >= 360.0 || ra == 0.0 {
        0.0
    } else {
        ra
    }
}

pub(super) fn check_position(ra_deg: f64, dec_deg: f64) -> CatalogResult<()> {
    if !ra_deg.is_finite() || !dec_deg.is_finite() {
        return Err(CatalogError::invalid_input(format!(
            "HEALPix position ({ra_deg}, {dec_deg}) is not finite"
        )));
    }
    if !(-90.0..=90.0).contains(&dec_deg) {
        return Err(CatalogError::invalid_input(format!(
            "declination {dec_deg} is outside [-90, 90] degrees"
        )));
    }
    Ok(())
}

// Gorski et al. (2005), Sec. 4: the ring of each base pixel's southern corner, counted from
// the north pole in units of nside, and the base pixel's central longitude in units of 45
// degrees.
const FACE_RING: [u64; 12] = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4];
const FACE_LONGITUDE: [i64; 12] = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7];

// Returns the centre of a NEST pixel as (RA, Dec) in degrees, with RA in [0, 360).
pub(crate) fn pix2ang_nest(order: u32, pixel: u64) -> CatalogResult<(f64, f64)> {
    let (nside, npix) = nside_and_npix(order)?;
    if pixel >= npix {
        return Err(CatalogError::invalid_input(format!(
            "HEALPix pixel {pixel} is outside [0, {npix}) at order {order}"
        )));
    }
    let face = (pixel >> (2 * order)) as usize;
    let (ix, iy) = pix2xy_nest(pixel & (nside * nside - 1), order);
    // Rings of pixel centres are numbered 1 to 4 nside - 1 from the north pole.
    let ring = FACE_RING[face] * nside - ix - iy - 1;
    let (ring_quarter, shift, dec_deg) = ring_geometry(ring, nside);
    // The sum is even on every face, so the division is exact. It is below 1 only on face 4
    // west of RA 0.
    let jp = (FACE_LONGITUDE[face] * ring_quarter + ix as i64 - iy as i64 + 1 + shift) / 2;
    let jp = if jp < 1 { jp + 4 * nside as i64 } else { jp };
    let ra_deg = (2 * jp - shift - 1) as f64 * 45.0 / ring_quarter as f64;
    Ok((ra_deg, dec_deg))
}

// Returns the number of pixels per quarter turn in the ring, whether its first pixel centre
// is offset by half a pixel from the quarter-turn start, and the ring's declination.
fn ring_geometry(ring: u64, nside: u64) -> (i64, i64, f64) {
    if ring < nside {
        (ring as i64, 0, 90.0 - polar_colatitude_deg(ring, nside))
    } else if ring > 3 * nside {
        let from_pole = 4 * nside - ring;
        (
            from_pole as i64,
            0,
            polar_colatitude_deg(from_pole, nside) - 90.0,
        )
    } else {
        // Both operands are exact, so z is rounded once.
        let z = (2 * (2 * nside as i64 - ring as i64)) as f64 / (3 * nside) as f64;
        let shift = ((ring - nside) & 1) as i64;
        (nside as i64, shift, libm::asin(z) * RAD_TO_DEG)
    }
}

// In the polar caps 1 - |z| = r^2 / (3 nside^2) for ring r from the pole, so
// sin(colatitude / 2) = r / (nside sqrt 6). This keeps the precision 1 - |z| loses near a pole.
fn polar_colatitude_deg(from_pole: u64, nside: u64) -> f64 {
    let sin_half = from_pole as f64 / (nside as f64 * libm::sqrt(6.0));
    2.0 * libm::asin(sin_half) * RAD_TO_DEG
}

/// Compute angular distance between two points on the sphere.
///
/// Converts both points to unit vectors and takes atan2 of the length of their cross
/// product and their dot product, which keeps full precision near 0 and 180 degrees.
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
    let sep_rad = angular_separation(
        Angle::from_degrees(ra1_deg).radians(),
        Angle::from_degrees(dec1_deg).radians(),
        Angle::from_degrees(ra2_deg).radians(),
        Angle::from_degrees(dec2_deg).radians(),
    );
    sep_rad * RAD_TO_DEG
}

/// Determine which of the 12 HEALPix base faces contains the point,
/// and compute the (ix, iy) position within that face.
fn compute_face_and_position(tt: f64, dec_deg: f64, nside: u64) -> (u32, u64, u64) {
    let z = Angle::from_degrees(dec_deg).sin();
    if libm::fabs(z) <= 2.0 / 3.0 {
        compute_equatorial_face(tt, z, nside)
    } else {
        compute_polar_face(tt, dec_deg, nside)
    }
}

/// Compute face and position for equatorial belt (-2/3 <= z <= 2/3).
fn compute_equatorial_face(tt: f64, z: f64, nside: u64) -> (u32, u64, u64) {
    // Working from the start of the quarter turn keeps 0.5 + frac below 1.5, where it is
    // exact for tt >= 1. 0.5 + tt rounds tt = 4 - 4.4e-16 up to 4.5, which loses the point's
    // face.
    let quarter = libm::floor(tt);
    let frac = tt - quarter;
    let temp1 = nside as f64 * (0.5 + frac);
    let temp2 = nside as f64 * z * 0.75;
    // jp and jm count the ascending and descending edge lines crossed since the quarter turn
    // began. Their exact values are in [0, 2 * nside) because frac < 1 and |0.75 z| <= 1/2.
    // The rounded sums can reach 2 * nside only when the exact value lies within a few ulps
    // below it, where the exact floor is 2 * nside - 1. A sum rounded below 0 casts to 0,
    // which is also the exact floor.
    let last = 2 * nside - 1;
    let jp = ((temp1 - temp2) as u64).min(last);
    let jm = ((temp1 + temp2) as u64).min(last);
    let quarter = quarter as u64;
    let face = compute_equatorial_face_number(quarter + jp / nside, quarter + jm / nside);
    let ix = jm % nside;
    let iy = nside - 1 - jp % nside;
    (face, ix, iy)
}

// ifp and ifm are in 0..=4. Equal values mean an equatorial face, and 4 is face 4 again past
// RA 360. Otherwise a smaller ifp names a northern face (0 to 3) and a smaller ifm a southern
// face (8 to 11).
fn compute_equatorial_face_number(ifp: u64, ifm: u64) -> u32 {
    let face = if ifp == ifm {
        ifp | 4
    } else if ifp < ifm {
        ifp
    } else {
        ifm + 8
    };
    face as u32
}

/// Compute face and position for polar caps (|z| > 2/3).
fn compute_polar_face(tt: f64, dec_deg: f64, nside: u64) -> (u32, u64, u64) {
    let tp = tt - libm::floor(tt);
    // 1 - |sin(dec)| = 2 sin^2(colatitude / 2). Subtracting from 1 instead has an absolute
    // error near 1e-16, which is over 0.5% of 1 - |sin(dec)| within 1e-5 degrees of a pole.
    let half_colat = Angle::from_degrees(90.0 - libm::fabs(dec_deg)).radians() * 0.5;
    let sin_half = libm::sin(half_colat);
    let tmp = nside as f64 * libm::sqrt(3.0 * (2.0 * sin_half * sin_half));
    let jp = (tp * tmp) as i64;
    let jm = ((1.0 - tp) * tmp) as i64;
    let jp = jp.min(nside as i64 - 1);
    let jm = jm.min(nside as i64 - 1);
    let ntt = libm::floor(tt) as u32;
    let north = dec_deg > 0.0;
    let face_offset = if north { 0 } else { 8 };
    let face = (ntt % 4) + face_offset;
    let (ix, iy) = if north {
        (nside as i64 - jm - 1, nside as i64 - jp - 1)
    } else {
        (jp, jm)
    };
    (face, ix as u64, iy as u64)
}

// Inverse of xy2pix_nest: ix from the even bits, iy from the odd bits.
fn pix2xy_nest(pixel_in_face: u64, order: u32) -> (u64, u64) {
    let (mut ix, mut iy) = (0, 0);
    for i in 0..order {
        ix |= ((pixel_in_face >> (2 * i)) & 1) << i;
        iy |= ((pixel_in_face >> (2 * i + 1)) & 1) << i;
    }
    (ix, iy)
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
pub(crate) mod tests {
    use super::*;

    #[test]
    fn test_xy2pix_nest() {
        assert_eq!(xy2pix_nest(0, 0, 2), 0);
        assert_eq!(xy2pix_nest(1, 0, 2), 1);
        assert_eq!(xy2pix_nest(0, 1, 2), 2);
        assert_eq!(xy2pix_nest(1, 1, 2), 3);
    }

    // Reference pixels are from independent implementations of the HEALPix projection
    // (Gorski et al. 2005, Sec. 4 and Fig. 4), one in 50-digit decimal arithmetic. Every
    // point is at least 0.05 pixel widths from a pixel edge, so rounding cannot move it,
    // except in test_ang2pix_nest_corner_at_ra_360, which pins near-ties on purpose.
    fn assert_pixels(order: u32, cases: &[(f64, f64, u64)]) {
        for &(ra, dec, expected) in cases {
            let pixel = ang2pix_nest(order, ra, dec).unwrap();
            assert_eq!(pixel, expected, "order {order}, ra {ra}, dec {dec}");
        }
    }

    #[test]
    fn test_ang2pix_nest_order0_every_base_face() {
        assert_pixels(
            0,
            &[
                (45.0, 60.0, 0),
                (45.0, 30.0, 0),
                (135.0, 60.0, 1),
                (225.0, 60.0, 2),
                (315.0, 60.0, 3),
                (332.72, 27.52, 3),
                (330.0, 30.0, 3),
                (0.0, 0.0, 4),
                (355.0, 1.0, 4),
                (340.29, -2.93, 4),
                (90.0, 0.0, 5),
                (180.0, 0.0, 6),
                (270.0, 0.0, 7),
                (45.0, -60.0, 8),
                (135.0, -60.0, 9),
                (225.0, -60.0, 10),
                (315.0, -60.0, 11),
                (329.52, -30.98, 11),
                (300.0, -20.0, 11),
                (330.0, -30.0, 11),
            ],
        );
    }

    #[test]
    fn test_ang2pix_nest_order1() {
        assert_pixels(
            1,
            &[
                (340.29, -2.93, 18),
                (329.52, -30.98, 45),
                (332.72, 27.52, 13),
                (10.0, 20.0, 19),
                (200.0, -70.0, 40),
                (357.0, 50.0, 13),
                (285.0, 10.0, 29),
            ],
        );
    }

    #[test]
    fn test_ang2pix_nest_order8() {
        assert_pixels(
            8,
            &[
                (340.29, -2.93, 298492),
                (329.52, -30.98, 752158),
                (332.72, 27.52, 215149),
                (83.633, -5.375, 364869),
                (10.0, 89.5, 65531),
                (300.0, -85.0, 721481),
                (359.9, 0.3, 311298),
                (275.0, 40.0, 239843),
            ],
        );
    }

    #[test]
    fn test_ang2pix_nest_order29() {
        assert_pixels(
            29,
            &[
                (340.29, -2.93, 1312781905329110755),
                (329.52, -30.98, 3308029027212997022),
                (332.72, 27.52, 946235892846093492),
                (123.456, 45.678, 464737556353774960),
                (299.1, -55.5, 3224518567304799058),
                (0.5, 0.5, 1369116848293459455),
            ],
        );
    }

    #[test]
    fn test_ang2pix_nest_reduces_any_finite_ra() {
        assert_pixels(
            8,
            &[
                (725.0, 10.0, 312927),
                (5.0, 10.0, 312927),
                (-400.0, -20.0, 776784),
                (320.0, -20.0, 776784),
                (1e6, 10.0, 509297),
                (-1e6, 10.0, 379570),
                (-1e-20, 0.4, 311299),
                (360.0, 0.4, 311299),
                (0.0, 0.4, 311299),
            ],
        );
    }

    #[test]
    fn test_reduce_ra_deg_returns_positive_zero() {
        for ra in [-0.0, 0.0, -360.0, 360.0, 720.0] {
            assert_eq!(reduce_ra_deg(ra).to_bits(), 0.0f64.to_bits(), "RA {ra}");
        }
    }

    #[test]
    fn test_ang2pix_nest_rejects_order_above_29() {
        for order in [30, 31, 32, 64] {
            let err = ang2pix_nest(order, 10.0, 20.0).unwrap_err().to_string();
            assert!(
                err.contains(&format!("HEALPix order {order} is above the maximum of 29")),
                "{err}"
            );
        }
    }

    #[test]
    fn test_ang2pix_nest_rejects_invalid_position() {
        let cases = [
            (f64::NAN, 10.0),
            (10.0, f64::NAN),
            (f64::INFINITY, 10.0),
            (f64::NEG_INFINITY, 10.0),
            (10.0, f64::INFINITY),
            (10.0, 95.0),
            (10.0, -90.000001),
        ];
        for (ra, dec) in cases {
            assert!(
                ang2pix_nest(8, ra, dec).is_err(),
                "({ra}, {dec}) was accepted"
            );
        }
    }

    #[test]
    fn test_nside_and_npix() {
        assert_eq!(nside_and_npix(0).unwrap(), (1, 12));
        assert_eq!(nside_and_npix(8).unwrap(), (256, 786432));
        assert_eq!(nside_and_npix(29).unwrap(), (1 << 29, 3458764513820540928));
        assert!(nside_and_npix(30).is_err());
    }

    // RA 359.99999999999994 and |Dec| = 41.810314895778596 put the point within 1e-15 of the
    // corner where faces 3, 4 and 0 (north) or 11, 4 and 8 (south) meet. The points are
    // 6e-16 of a base pixel width inside faces 3 and 11.
    #[test]
    fn test_ang2pix_nest_corner_at_ra_360() {
        let ra = 359.99999999999994;
        let dec = 41.810314895778596;
        assert_pixels(0, &[(ra, -dec, 11), (ra, dec, 3)]);
        assert_pixels(8, &[(ra, -dec, 742741), (ra, dec, 218453)]);
        assert_pixels(
            29,
            &[
                (ra, -dec, 3266610929719399765),
                (ra, dec, 960767920505705813),
            ],
        );
    }

    #[test]
    fn test_ang2pix_nest_order29_within_1e5_degrees_of_poles() {
        assert_pixels(
            29,
            &[
                (26.077063, 89.999999594822, 288230376151711736),
                (227.02533, 89.999999438964, 864691128455135216),
                (207.757062, 89.999999845104, 864691128455135230),
                (222.843454, 89.999999691508, 864691128455135228),
                (341.23153, -89.999999847078, 3170534137668829185),
                (209.803682, -89.999994640649, 2882303761517119890),
                (150.890056, -89.999999581151, 2594073385365405703),
                (201.69262, -89.999998888248, 2882303761517117575),
            ],
        );
    }

    #[test]
    fn test_ang2pix_nest_whole_turns_keep_boundary_pixels() {
        for order in [0u32, 3, 8] {
            let nside = 1u32 << order;
            let step = (nside / 4).max(1) as usize;
            for dec in [0.0, 10.0, -20.0, 60.0, -75.0] {
                for j in (0..4 * nside).step_by(step) {
                    let ra = f64::from(j) * 90.0 / f64::from(nside);
                    let pixel = ang2pix_nest(order, ra, dec).unwrap();
                    for turns in [-1.0, 1.0, -1000.0] {
                        let shifted = ang2pix_nest(order, ra + 360.0 * turns, dec).unwrap();
                        assert_eq!(
                            shifted, pixel,
                            "order {order}, ra {ra}, dec {dec}, {turns} turns"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_ang2pix_nest_poles() {
        let north_pole = ang2pix_nest(0, 0.0, 90.0).unwrap();
        assert!(north_pole < 12);
        let south_pole = ang2pix_nest(0, 0.0, -90.0).unwrap();
        assert!(south_pole < 12);
    }

    #[test]
    fn test_ang2pix_nest_equator() {
        let pixel = ang2pix_nest(0, 0.0, 0.0).unwrap();
        assert!(pixel < 12);
    }

    #[test]
    fn test_ang2pix_nest_order8_bounds() {
        let nside = 1u64 << 8;
        let npix = 12 * nside * nside;
        for ra in [0.0, 90.0, 180.0, 270.0] {
            for dec in [-89.0, -45.0, 0.0, 45.0, 89.0] {
                let pixel = ang2pix_nest(8, ra, dec).unwrap();
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
        // (ra1, dec1, ra2, dec2, separation): ERFA eraSeps on the inputs times
        // DEG_TO_RAD, with the result times RAD_TO_DEG.
        let cases = [
            (0.0, 0.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 90.0, 0.0, 90.0),
            (0.0, 90.0, 0.0, 0.0, 90.0),
            (0.0, 0.0, 180.0, 0.0, 180.0),
            (0.0, 90.0, 0.0, -90.0, 180.0),
            (0.0, 0.0, 0.1, 0.1, 0.14142132033778318),
        ];
        for (ra1, dec1, ra2, dec2, expected) in cases {
            assert_eq!(
                angular_separation_deg(ra1, dec1, ra2, dec2),
                expected,
                "({ra1}, {dec1}) to ({ra2}, {dec2})"
            );
        }
    }

    #[test]
    fn test_pix2ang_nest_round_trips_every_pixel_to_order_4() {
        for order in 0..=4 {
            let (_, npix) = nside_and_npix(order).unwrap();
            for pixel in 0..npix {
                let (ra, dec) = pix2ang_nest(order, pixel).unwrap();
                assert_eq!(
                    ang2pix_nest(order, ra, dec).unwrap(),
                    pixel,
                    "order {order}"
                );
            }
        }
    }

    // Each base pixel's four corner pixels, then pixels spread by a 64-bit LCG.
    pub(crate) fn sample_pixels(order: u32) -> Vec<u64> {
        let (nside, npix) = nside_and_npix(order).unwrap();
        let face_size = nside * nside;
        let last = nside - 1;
        let mut pixels = Vec::new();
        for face in 0..12 {
            for (ix, iy) in [(0, 0), (last, 0), (0, last), (last, last)] {
                pixels.push(face * face_size + xy2pix_nest(ix, iy, order));
            }
        }
        let mut state = 0x9e3779b97f4a7c15u64;
        for _ in 0..20_000 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            pixels.push(state % npix);
        }
        pixels
    }

    #[test]
    fn test_pix2ang_nest_round_trips_sampled_pixels() {
        for order in [8, 20, 29] {
            for pixel in sample_pixels(order) {
                let (ra, dec) = pix2ang_nest(order, pixel).unwrap();
                assert_eq!(
                    ang2pix_nest(order, ra, dec).unwrap(),
                    pixel,
                    "order {order}"
                );
            }
        }
    }

    // Base pixel centres are at z = 2/3, 0 and -2/3 (Gorski et al. 2005, Sec. 4). 41.81... is
    // asin(2/3) in degrees, correctly rounded from a 50-digit evaluation.
    #[test]
    fn test_pix2ang_nest_base_pixel_centres() {
        let cap = 41.810314895778596;
        let expected = [
            (45.0, cap),
            (135.0, cap),
            (225.0, cap),
            (315.0, cap),
            (0.0, 0.0),
            (90.0, 0.0),
            (180.0, 0.0),
            (270.0, 0.0),
            (45.0, -cap),
            (135.0, -cap),
            (225.0, -cap),
            (315.0, -cap),
        ];
        for (pixel, centre) in expected.into_iter().enumerate() {
            assert_eq!(
                pix2ang_nest(0, pixel as u64).unwrap(),
                centre,
                "pixel {pixel}"
            );
        }
    }

    // The ring next to each pole holds four pixels at RA 45 + 90 k and z = +-(1 - 1 / (3 nside^2)),
    // the north one at ix = iy = nside - 1 of base pixels 0 to 3 and the south one at ix = iy = 0
    // of base pixels 8 to 11. The declinations are asin(11/12) at order 1 and
    // 90 - 2 asin(1 / (nside sqrt 6)) at orders 8 and 29, in degrees, each correctly rounded
    // from a 50-digit evaluation.
    #[test]
    fn test_pix2ang_nest_pixels_next_to_the_poles() {
        let rings = [
            (1, 66.44353569089877),
            (8, 89.81725848475484),
            (29, 89.9999999128621),
        ];
        for (order, dec) in rings {
            let face_size = 1u64 << (2 * order);
            for k in 0..4 {
                let ra = 45.0 + 90.0 * k as f64;
                let north = pix2ang_nest(order, k * face_size + face_size - 1).unwrap();
                let south = pix2ang_nest(order, (8 + k) * face_size).unwrap();
                assert_eq!((north, south), ((ra, dec), (ra, -dec)), "order {order}");
            }
        }
    }

    #[test]
    fn test_pix2ang_nest_rejects_bad_input() {
        assert!(pix2ang_nest(30, 0).is_err());
        assert!(pix2ang_nest(0, 12).is_err());
        assert!(pix2ang_nest(8, 786432).is_err());
        assert!(pix2ang_nest(29, u64::MAX).is_err());
    }
}
