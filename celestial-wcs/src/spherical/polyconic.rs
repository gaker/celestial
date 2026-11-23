use crate::common::{check_nonzero_param, native_coord};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::WcsResult;
use crate::newton::bracketed_newton;
use crate::trig::{atan2d, cotd, sincosd, sind, to_degrees, to_radians};

use super::pseudocylindrical::{deproject_sfl, parallel_longitude, project_sfl};

// Paper II §5.5.1: BON diverges at θ1 = 0, which is handled as SFL.
// Eqs. 146 to 150, with eq. 150's y = Y0 − Rθ cos A written as
// θ + 2Rθ sin²(A/2), since Rθ = Y0 − θ. Near the central meridian the printed
// difference cancels.
pub(crate) fn project_bon(
    native: NativeCoord,
    theta_1: f64,
    y0: &WcsResult<f64>,
) -> WcsResult<IntermediateCoord> {
    if theta_1 == 0.0 {
        return project_sfl(native);
    }
    let y0 = y0.clone()?;
    let (phi, theta) = (native.phi_deg(), native.theta_deg());
    let r_theta = y0 - theta;
    let a_phi = to_degrees(phi * sincosd(theta).1 / r_theta);
    let sin_half_a = sind(a_phi / 2.0);
    Ok(IntermediateCoord::new(
        r_theta * sind(a_phi),
        theta + 2.0 * r_theta * sin_half_a * sin_half_a,
    ))
}

// Eqs. 151 to 154. Rθ takes the sign of θ1, and dividing both arguments of
// arg() by it only flips their signs with θ1's.
pub(crate) fn deproject_bon(
    inter: IntermediateCoord,
    theta_1: f64,
    y0: &WcsResult<f64>,
) -> WcsResult<NativeCoord> {
    if theta_1 == 0.0 {
        return deproject_sfl(inter);
    }
    let y0 = y0.clone()?;
    let (x, y_offset) = (inter.x_deg(), y0 - inter.y_deg());
    let sign = libm::copysign(1.0, theta_1);
    let r_theta = sign * libm::sqrt(x * x + y_offset * y_offset);
    let theta = y0 - r_theta;
    let a_phi = atan2d(sign * x, sign * y_offset);
    let phi = parallel_longitude(to_radians(a_phi) * r_theta, sincosd(theta).1)?;
    native_coord(phi, theta)
}

pub(super) fn bonne_y0(theta_1: f64) -> WcsResult<f64> {
    check_nonzero_param(theta_1, "BON projection: theta_1")?;
    Ok(to_degrees(cotd(theta_1)) + theta_1)
}

// Eqs. 155 and 156, with 1 − cos E taken as 2 sin²(E/2): near the equator
// E = φ sin θ is tiny and the difference would round to 0.
pub(crate) fn project_pco(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let (phi, theta) = (native.phi_deg(), native.theta_deg());
    if theta == 0.0 {
        return Ok(IntermediateCoord::new(phi, 0.0));
    }
    let cot_theta = cotd(theta);
    let e = phi * sind(theta);
    let sin_half_e = sind(e / 2.0);
    let x = to_degrees(cot_theta * sind(e));
    let y = theta + to_degrees(cot_theta * 2.0 * sin_half_e * sin_half_e);
    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_pco(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (x, y) = (inter.x_deg(), inter.y_deg());
    if y == 0.0 {
        return native_coord(x, 0.0);
    }
    let theta = libm::copysign(solve_pco_inverse(x, libm::fabs(y))?, y);
    // Eq. 158. Its arguments are cos E and sin E for E = φ sin θ, here scaled
    // by cos θ so that no tangent is needed near the pole.
    let (sin_theta, cos_theta) = sincosd(theta);
    let e = atan2d(
        x * sin_theta,
        to_degrees(cos_theta) - (y - theta) * sin_theta,
    );
    native_coord(e / sin_theta, theta)
}

// Eq. 157 multiplied through by sin θ: k(θ) = (x² + w²) sin θ − (360°/π) w cos θ
// with w = y − θ. The projection is symmetric about the equator, so y > 0
// suffices. Then k rises, with slope (π/180° (x² + w²) + 360°/π) cos θ per
// degree, from −(360°/π) y at the equator to x² sin y ≥ 0 at θ = y, or to
// x² + w² at the pole when y is beyond it, and its one root lies between.
fn solve_pco_inverse(x: f64, y: f64) -> WcsResult<f64> {
    let top = y.min(90.0);
    let k_and_slope = |theta| {
        let (sin_theta, cos_theta) = sincosd(theta);
        let w = y - theta;
        let r2 = x * x + w * w;
        let k = r2 * sin_theta - to_degrees(2.0 * w) * cos_theta;
        Ok((k, (to_radians(r2) + to_degrees(2.0)) * cos_theta))
    };
    bracketed_newton((0.0, top), top, k_and_slope, "PCO inverse")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{assert_pinned, PinnedRow};
    use crate::spherical::Projection;

    // Per-projection native_reference checks are covered by
    // spherical::tests::test_all_projections_map_reference_to_origin_and_back.

    // Eq. 145 puts (0°, 0°) at the origin. There Rθ = Y0 and A = 0, so
    // y = Y0 − Rθ is exactly 0.
    #[test]
    fn test_bon_reference_maps_to_origin() {
        for theta_1 in [30.0, 45.0, 60.0, 75.0, -45.0] {
            let inter = Projection::bon(theta_1)
                .project(NativeCoord::new(0.0, 0.0))
                .unwrap();
            assert_eq!((inter.x_deg(), inter.y_deg()), (0.0, 0.0), "{theta_1}");
        }
    }

    // Near the central meridian Y0 − Rθ cos A cancels. Against a 60-digit
    // evaluation y at (1e-6°, 0°) is 0.55 ulp high and at (−30°, 0°) 0.57 ulp
    // low.
    #[test]
    fn test_bon_y_near_central_meridian() {
        let y = |phi| {
            let native = NativeCoord::new(phi, 0.0);
            Projection::bon(45.0).project(native).unwrap().y_deg()
        };
        assert_eq!(y(1e-6), 4.887787183204918e-15);
        assert_eq!(y(-30.0), 4.367570407885135);
    }

    // Against a 60-digit evaluation x is within 18.45 ulp, at (−180°, 70°)
    // with θ1 = 60°: Rθ = 93.08 − 70 carries the rounding of Y0, and at
    // A = −152.8° x = Rθ sin A is 5.2 times as sensitive, relatively, to A.
    // y is within 1.73 ulp. The deprojected values are within 8.95 ulp, θ at
    // (120°, 10°) with θ1 = 60°, where θ = Y0 − Rθ is 93.08 − 83.08. At
    // (−180°, 20°) and (−180°, 70°) the projected point is past the seam, and φ
    // comes back just below 180°.
    #[test]
    fn test_bon_pinned_both_ways() {
        let cases: [(f64, &[PinnedRow]); 3] = [
            (
                45.0,
                &[
                    (
                        (-30.0, 0.0),
                        (-29.5718171176194, 4.367570407885135),
                        (-30.0, 0.0),
                    ),
                    (
                        (120.0, 10.0),
                        (88.43180352003323, 75.86990124722809),
                        (119.99999999999997, 9.999999999999986),
                    ),
                    (
                        (-180.0, 70.0),
                        (-30.49571407629065, 112.92727975221868),
                        (179.99999999999991, 70.0),
                    ),
                    (
                        (-180.0, 20.0),
                        (-72.82304387095974, 140.62852800616764),
                        (179.99999999999994, 20.0),
                    ),
                    ((30.0, 90.0), (0.0, 90.0), (0.0, 90.0)),
                ],
            ),
            (
                60.0,
                &[
                    (
                        (-180.0, 70.0),
                        (-10.538013684345579, 113.61322773922643),
                        (180.0, 70.0),
                    ),
                    (
                        (120.0, 10.0),
                        (82.16728046637347, 80.80047574168647),
                        (119.99999999999999, 10.000000000000014),
                    ),
                ],
            ),
            (
                -45.0,
                &[
                    (
                        (120.0, -10.0),
                        (88.43180352003323, -75.86990124722809),
                        (119.99999999999997, -9.999999999999986),
                    ),
                    (
                        (-60.0, -50.0),
                        (-35.165102950094266, -63.58836304814382),
                        (-60.0, -50.0),
                    ),
                    ((60.0, -90.0), (0.0, -90.0), (0.0, -90.0)),
                ],
            ),
        ];
        for (theta_1, rows) in cases {
            assert_pinned(&Projection::bon(theta_1), rows);
        }
    }

    // On the equator PCO is x = φ, y = 0, and every step is odd in φ, so the
    // central meridian is a mirror bit for bit.
    #[test]
    fn test_pco_equator_and_mirror_are_exact() {
        let proj = Projection::pco();
        let inter = proj.project(NativeCoord::new(45.0, 0.0)).unwrap();
        assert_eq!((inter.x_deg(), inter.y_deg()), (45.0, 0.0));
        let native = proj.deproject(IntermediateCoord::new(30.0, 0.0)).unwrap();
        assert_eq!((native.phi_deg(), native.theta_deg()), (30.0, 0.0));

        let east = proj.project(NativeCoord::new(30.0, 45.0)).unwrap();
        let west = proj.project(NativeCoord::new(-30.0, 45.0)).unwrap();
        assert_eq!((west.x_deg(), west.y_deg()), (-east.x_deg(), east.y_deg()));
    }

    // Against a 60-digit evaluation the projected values are within 1.56 ulp
    // and the deprojected ones within 2.46 ulp. At (120°, 45°) y passes 90°.
    // At (−180°, 30°) the projected point is just past the seam, and φ comes
    // back as 180°.
    #[test]
    fn test_pco_pinned_both_ways() {
        let rows = [
            (
                (-30.0, -45.0),
                (-20.73187096094458, -48.88233666888357),
                (-29.999999999999982, -45.0),
            ),
            (
                (30.0, 20.0),
                (28.04033977032675, 22.51748224146966),
                (30.0, 20.0),
            ),
            (
                (120.0, -20.0),
                (103.3639904954456, -58.689828222688874),
                (120.00000000000003, -19.999999999999993),
            ),
            (
                (120.0, 45.0),
                (57.06473532770469, 97.15551376209666),
                (119.99999999999996, 45.0),
            ),
            (
                (-180.0, 30.0),
                (-99.23920117592256, 129.23920117592255),
                (180.0, 30.000000000000007),
            ),
            ((-60.0, -90.0), (0.0, -90.0), (0.0, -90.0)),
        ];
        assert_pinned(&Projection::pco(), &rows);
    }

    #[test]
    fn test_pco_inverse_reaches_the_poles_and_far_longitudes() {
        // The forwards of (37.5, ±89.9), (150, 60), (120, −70) and (179, 30),
        // which lie near a pole or where φ sin θ passes 90°. Against a 60-digit
        // evaluation, (150, 60) has φ 0.94 ulp high and the rest are correctly
        // rounded. Near the pole an ulp of y moves φ by about 700 ulp, so the
        // rounding of the forward's y leaves φ 1.7e-12 from 37.5.
        let cases = [
            (
                (0.06087612562789088, 89.92066462626875),
                (37.50000000000171, 89.9),
            ),
            (
                (0.06087612562789088, -89.92066462626875),
                (37.50000000000171, -89.9),
            ),
            (
                (25.376207679373938, 114.3004044769716),
                (150.00000000000006, 60.0),
            ),
            ((19.229694214294685, -98.92281427410889), (120.0, -70.0)),
            (
                (99.23542245122769, 128.3731867640306),
                (179.0, 30.000000000000004),
            ),
        ];
        for ((x, y), expected) in cases {
            let native = deproject_pco(IntermediateCoord::new(x, y)).unwrap();
            let got = (native.phi_deg(), native.theta_deg());
            assert_eq!(got, expected, "({x}, {y})");
        }
    }

    // Eqs. 155 and 156 hold for every θ ≠ 0, with x = φ, y = 0 only at θ = 0.
    // Near the equator E = φ sin θ is tiny, and the term (1 − cos E) cot θ
    // is about φ²θ/2, as large as θ itself. Against 60 digits, at φ = 37.5°
    // x is 1.00 ulp below 37.5 and y 0.52 ulp below 1.2141841232875297e-9;
    // at φ = 179° both are correctly rounded. The inverses are of those x and
    // y: at φ = 37.5° correctly rounded, at φ = 179° φ 1.00 ulp low and θ
    // 0.92 ulp high.
    #[test]
    fn test_pco_keeps_latitudes_near_the_equator() {
        let cases = [
            (
                (37.5, 1e-9),
                (37.49999999999999, 1.2141841232875295e-9),
                (37.49999999999999, 1e-9),
            ),
            (
                (37.5, -1e-9),
                (37.49999999999999, -1.2141841232875295e-9),
                (37.49999999999999, -1e-9),
            ),
            (
                (179.0, 1e-12),
                (179.0, 5.880123373692965e-12),
                (178.99999999999997, 1.0000000000000002e-12),
            ),
        ];
        for ((phi, theta), forward, inverse) in cases {
            let inter = project_pco(NativeCoord::new(phi, theta)).unwrap();
            assert_eq!((inter.x_deg(), inter.y_deg()), forward, "({phi}, {theta})");
            let back = deproject_pco(inter).unwrap();
            assert_eq!((back.phi_deg(), back.theta_deg()), inverse, "{forward:?}");
        }
    }

    // Paper II §5.5.1: BON diverges at θ1 = 0, which is handled as SFL.
    #[test]
    fn test_bon_zero_is_sanson_flamsteed() {
        let (bon, sfl) = (Projection::bon(0.0), Projection::sfl());
        let native = NativeCoord::new(37.5, 45.0);
        let inter = sfl.project(native).unwrap();
        let got = bon.project(native).unwrap();
        assert_eq!((got.x_deg(), got.y_deg()), (inter.x_deg(), inter.y_deg()));
        let expected = sfl.deproject(inter).unwrap();
        let got = bon.deproject(inter).unwrap();
        assert_eq!(
            (got.phi_deg(), got.theta_deg()),
            (expected.phi_deg(), expected.theta_deg())
        );
    }
}
