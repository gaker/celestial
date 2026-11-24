//! Disc queries over the HEALPix NEST scheme.

use celestial_core::angle::Angle;

use crate::error::{CatalogError, CatalogResult};

use super::healpix::{angular_separation_deg, check_position, nside_and_npix, pix2ang_nest};

// A whole-sky query is allowed up to order 10, 12 * 4^10 pixels or 100 MB of indices. Larger
// discs fail instead of allocating and scanning: the whole sky at order 29 is 3.5e18 pixels.
const MAX_DISC_PIXELS: u64 = 12 << 20;

// An upper bound on the distance from a pixel's centre to any of its points, in units of
// 45 / nside degrees, the half-diagonal h of the pixel in the HEALPix projection plane
// (Calabretta & Roukema 2007, eqs. 1-3). Every point of the pixel is within |dx| + |dy| <= h
// of the centre in the plane, and the straight path between them stays in the pixel. Along
// it the arc on the sphere has ds^2 = ds_theta^2 + ds_phi^2, with
// |ds_theta| <= 8 / (pi sqrt 5) |dy| everywhere, |ds_phi| <= |dx| in the equatorial belt and
// |ds_phi| <= sqrt(2/3) (|dx| + |dy|) in the polar caps. Over |dx| + |dy| <= h both bounds
// peak at dx = 0, so the distance is at most sqrt(64 / (5 pi^2) + 2/3) h = 1.40130 h. This
// covers the curved pixel edges, which the centre-to-corner distance does not.
const PIXEL_RADIUS_PER_HALF_DIAGONAL: f64 = 1.41;

// Covers rounding in pixel centres, in ang2pix_nest near pixel edges and in the separations,
// each below 1e-13 degrees.
const ROUNDING_MARGIN_DEG: f64 = 1e-9;

/// Returns the NEST pixels at `order` that may overlap a cone, in ascending order.
///
/// Every pixel that holds a point within `radius_deg` of the centre is returned. Pixels near
/// the edge that do not overlap the cone may be returned too.
///
/// # Errors
/// Returns an error if `order` is above 29, the centre is not a valid position, the radius is
/// not in [0, 180] degrees, or the cone covers more than 12 * 4^10 pixels' worth of the sky.
pub(crate) fn query_disc_nest(
    order: u32,
    ra_deg: f64,
    dec_deg: f64,
    radius_deg: f64,
) -> CatalogResult<Vec<u64>> {
    let (_, npix) = nside_and_npix(order)?;
    check_position(ra_deg, dec_deg)?;
    check_disc_size(order, npix, radius_deg)?;
    let cone = (ra_deg, dec_deg, radius_deg);
    // Descending from the base pixels keeps every pixel whose bounding circle meets the cone.
    // A pixel holding a point of the cone is always kept, and so is each of its ancestors.
    let mut pixels = pixels_near_cone(0, 0..12, cone)?;
    for level in 1..=order {
        let children = pixels.iter().flat_map(|&p| 4 * p..4 * p + 4);
        pixels = pixels_near_cone(level, children, cone)?;
    }
    Ok(pixels)
}

// Every point of the cone lies in a returned pixel, so a query returns at least
// npix sin^2(radius / 2) pixels, the cone's share of the sphere.
fn check_disc_size(order: u32, npix: u64, radius_deg: f64) -> CatalogResult<()> {
    if !(0.0..=180.0).contains(&radius_deg) {
        return Err(CatalogError::invalid_input(format!(
            "disc radius {radius_deg} is outside [0, 180] degrees"
        )));
    }
    let sin_half = Angle::from_degrees(radius_deg / 2.0).sin();
    let fewest = npix as f64 * sin_half * sin_half;
    if fewest > MAX_DISC_PIXELS as f64 {
        return Err(CatalogError::invalid_input(format!(
            "a disc of radius {radius_deg} degrees covers at least {fewest:.3e} pixels at \
                 HEALPix order {order}, above the limit of {MAX_DISC_PIXELS}"
        )));
    }
    Ok(())
}

fn pixels_near_cone(
    order: u32,
    candidates: impl Iterator<Item = u64>,
    (ra_deg, dec_deg, radius_deg): (f64, f64, f64),
) -> CatalogResult<Vec<u64>> {
    let half_diagonal_deg = 45.0 / (1u64 << order) as f64;
    let reach =
        radius_deg + PIXEL_RADIUS_PER_HALF_DIAGONAL * half_diagonal_deg + ROUNDING_MARGIN_DEG;
    let mut kept = Vec::new();
    for pixel in candidates {
        let (pixel_ra, pixel_dec) = pix2ang_nest(order, pixel)?;
        if angular_separation_deg(ra_deg, dec_deg, pixel_ra, pixel_dec) <= reach {
            kept.push(pixel);
        }
    }
    Ok(kept)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query::healpix::ang2pix_nest;
    use crate::query::healpix::tests::sample_pixels;

    #[test]
    fn test_query_disc_nest_basic() {
        let nside = 16u64;
        let ra = 0.0;
        let dec = 0.0;
        let radius = 10.0;

        let pixels = query_disc_nest(nside.trailing_zeros(), ra, dec, radius).unwrap();

        // Should return some pixels
        assert!(!pixels.is_empty());

        // All pixels should be in valid range
        let npix = 12 * nside * nside;
        for &pix in &pixels {
            assert!(pix < npix);
        }

        // Center pixel should be included
        let order = nside.trailing_zeros();
        let center_pixel = ang2pix_nest(order, ra, dec).unwrap();
        assert!(pixels.contains(&center_pixel));
    }

    #[test]
    fn test_query_disc_nest_pole() {
        let nside = 16u64;
        let ra = 0.0;
        let dec = 90.0;
        let radius = 5.0;

        let pixels = query_disc_nest(nside.trailing_zeros(), ra, dec, radius).unwrap();

        // Should return pixels around north pole
        assert!(!pixels.is_empty());

        // Center pixel should be included
        let order = nside.trailing_zeros();
        let center_pixel = ang2pix_nest(order, ra, dec).unwrap();
        assert!(pixels.contains(&center_pixel));
    }

    // Centres of pixels `depth` orders finer lie inside the pixel, within 1/2^depth of a pixel
    // width of every point on its edge, so they test the bound right up to the edges.
    fn assert_pixel_radius_bounds(order: u32, pixels: impl Iterator<Item = u64>, depth: u32) {
        let bound = PIXEL_RADIUS_PER_HALF_DIAGONAL * 45.0 / (1u64 << order) as f64;
        for pixel in pixels {
            let (ra, dec) = pix2ang_nest(order, pixel).unwrap();
            let first = pixel << (2 * depth);
            for sub in first..first + (1u64 << (2 * depth)) {
                let (sub_ra, sub_dec) = pix2ang_nest(order + depth, sub).unwrap();
                let distance = angular_separation_deg(ra, dec, sub_ra, sub_dec);
                assert!(
                    distance <= bound,
                    "order {order} pixel {pixel}: {distance} > {bound}"
                );
            }
        }
    }

    #[test]
    fn test_pixel_radius_bounds_every_point_of_the_pixel() {
        for order in 0..=2 {
            let (_, npix) = nside_and_npix(order).unwrap();
            assert_pixel_radius_bounds(order, 0..npix, 6);
        }
        // Pixels next to the poles, on the z = 2/3 rings and in the equatorial belt.
        for order in [12, 20] {
            let pixels = sample_pixels(order).into_iter().take(48 + 200);
            assert_pixel_radius_bounds(order, pixels, 5);
        }
    }

    #[test]
    fn test_query_disc_nest_returns_ascending_pixels() {
        for (order, ra, dec, radius) in [
            (4, 10.0, 20.0, 15.0),
            (6, 0.0, -90.0, 3.0),
            (3, 300.0, 41.8, 40.0),
        ] {
            let pixels = query_disc_nest(order, ra, dec, radius).unwrap();
            assert!(pixels.len() > 4, "order {order}");
            assert!(
                pixels.windows(2).all(|w| w[0] < w[1]),
                "order {order}: {pixels:?}"
            );
        }
    }

    // Each pixel `depth` orders finer whose centre is in the cone lies in a pixel that holds a
    // point of the cone, so its parent must be returned.
    #[test]
    fn test_query_disc_nest_covers_every_fine_pixel_in_the_cone() {
        let cones = [
            (0.0, -90.0, 2.0),
            (200.0, 89.5, 1.0),
            (359.95, 0.0, 0.1),
            (45.0, 41.810314895778596, 3.0),
            (100.0, -41.810314895778596, 10.0),
            (80.0, 20.0, 90.0),
            (10.0, 10.0, 0.0),
        ];
        for (order, depth) in [(0, 5), (3, 4), (6, 2)] {
            let (_, fine_npix) = nside_and_npix(order + depth).unwrap();
            for (ra, dec, radius) in cones {
                let pixels = query_disc_nest(order, ra, dec, radius).unwrap();
                for fine in 0..fine_npix {
                    let (fine_ra, fine_dec) = pix2ang_nest(order + depth, fine).unwrap();
                    if angular_separation_deg(ra, dec, fine_ra, fine_dec) <= radius {
                        let parent = fine >> (2 * depth);
                        assert!(
                            pixels.binary_search(&parent).is_ok(),
                            "order {order}, cone ({ra}, {dec}) radius {radius}: no pixel {parent}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_query_disc_nest_whole_sky_is_every_pixel() {
        let pixels = query_disc_nest(8, 10.0, -60.0, 180.0).unwrap();
        assert_eq!(pixels, (0..786432).collect::<Vec<u64>>());
    }

    #[test]
    fn test_query_disc_nest_rejects_discs_above_the_pixel_limit() {
        for (order, radius) in [(29, 180.0), (29, 1.0), (29, 0.001), (11, 180.0), (12, 90.0)] {
            let err = query_disc_nest(order, 0.0, 0.0, radius)
                .unwrap_err()
                .to_string();
            assert!(
                err.contains("above the limit of 12582912"),
                "order {order}, radius {radius}: {err}"
            );
        }
        let at_corner = query_disc_nest(29, 0.0, 0.0, 0.0).unwrap();
        assert_eq!(at_corner.len(), 4);
        assert_eq!(query_disc_nest(29, 0.0, 90.0, 1e-5).unwrap().len(), 27_144);
    }
}
