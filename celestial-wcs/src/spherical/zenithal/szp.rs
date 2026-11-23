use crate::common::{pole_native_coord, ZERO_LIMIT};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{sincosd, to_degrees, to_radians};

use super::sin::deproject_slant;

// Paper II eqs. 41 and 42.
pub(crate) fn project_szp(
    native: NativeCoord,
    point: &WcsResult<(f64, f64, f64)>,
) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    if theta == 90.0 {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }

    let point = point.clone()?;
    let (xp, yp, zp) = point;
    let (sin_theta, cos_theta) = sincosd(theta);
    let (sin_phi, cos_phi) = sincosd(native.phi_deg());

    let denom = zp - (1.0 - sin_theta);
    check_projectable(point, (sin_phi, cos_phi), (sin_theta, cos_theta), denom)?;

    let x = to_degrees((zp * cos_theta * sin_phi - xp * (1.0 - sin_theta)) / denom);
    let y = -to_degrees((zp * cos_theta * cos_phi + yp * (1.0 - sin_theta)) / denom);

    Ok(IntermediateCoord::new(x, y))
}

// Of the two points where a ray from P meets the sphere, the inverse
// (eq. 38) returns the one closer to the pole, nearer the plane. Which one
// a point is depends on its side of the limb, where σ cos θ − ρ sin θ = 1
// with ρ and σ from eq. 51 (eq. 48), and on which side of the plane P is.
// Past θ∞ (eq. 53), where the projection diverges, the denominator has
// lost the sign of zp.
fn check_projectable(
    (xp, yp, zp): (f64, f64, f64),
    (sin_phi, cos_phi): (f64, f64),
    (sin_theta, cos_theta): (f64, f64),
    denom: f64,
) -> WcsResult<()> {
    if denom.abs() < ZERO_LIMIT {
        return Err(WcsError::singularity(
            "SZP projection singularity: denominator = 0",
        ));
    }
    let limb = (xp * sin_phi - yp * cos_phi) * cos_theta - (zp - 1.0) * sin_theta;
    if zp * denom < 0.0 || zp * (1.0 - limb) < 0.0 {
        return Err(WcsError::out_of_bounds(
            "SZP projection: point beyond the boundary",
        ));
    }
    Ok(())
}

pub(crate) fn deproject_szp(
    inter: IntermediateCoord,
    point: &WcsResult<(f64, f64, f64)>,
) -> WcsResult<NativeCoord> {
    let x_big = to_radians(inter.x_deg());
    let y_big = to_radians(inter.y_deg());

    if x_big == 0.0 && y_big == 0.0 {
        return Ok(pole_native_coord());
    }

    let (xp, yp, zp) = point.clone()?;
    let x_prime = (x_big - xp) / zp;
    let y_prime = (y_big - yp) / zp;

    deproject_slant(x_big, y_big, x_prime, y_prime, "SZP")
}

// The point of projection P (eqs. 45-47), in a frame with the plane of
// projection at z = 0 and the centre of the sphere at z = 1.
pub(crate) fn projection_point(
    mu: f64,
    phi_c_deg: f64,
    theta_c_deg: f64,
) -> WcsResult<(f64, f64, f64)> {
    let (sin_phi_c, cos_phi_c) = sincosd(phi_c_deg);
    let (sin_theta_c, cos_theta_c) = sincosd(theta_c_deg);

    let xp = -mu * cos_theta_c * sin_phi_c;
    let yp = mu * cos_theta_c * cos_phi_c;
    let zp = mu * sin_theta_c + 1.0;

    if zp.abs() < ZERO_LIMIT {
        return Err(WcsError::singularity("SZP projection singularity: zp = 0"));
    }
    Ok((xp, yp, zp))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // (µ, φc, θc) configurations that take both signs of the inverse's
    // quadratic. Against a 60-digit evaluation the projected values are within
    // 1.62 ulp and the deprojected ones within 1.87 ulp.
    #[test]
    fn test_szp_pinned_both_ways() {
        let cases = [
            (
                (0.0, 0.0, 90.0),
                (
                    (45.0, 60.0),
                    (23.390904037010284, -23.390904037010287),
                    (45.0, 60.00000000000001),
                ),
            ),
            (
                (2.0, 30.0, 60.0),
                (
                    (20.0, 75.0),
                    (5.497763465756207, -14.737590137322195),
                    (19.999999999999996, 75.0),
                ),
            ),
            (
                (0.5, 0.0, 90.0),
                (
                    (30.0, 30.0),
                    (37.214700440970965, -64.4577519522176),
                    (29.999999999999996, 30.000000000000004),
                ),
            ),
            (
                (3.0, 45.0, 45.0),
                (
                    (60.0, 20.0),
                    (82.03800749825356, -57.06714315449885),
                    (59.99999999999999, 20.000000000000004),
                ),
            ),
            (
                (2.0, 30.0, 70.0),
                (
                    (45.0, 50.0),
                    (30.078196475116847, -31.346878198782306),
                    (44.99999999999999, 50.00000000000001),
                ),
            ),
        ];
        for ((mu, phi_c, theta_c), row) in cases {
            assert_pinned(&Projection::szp(mu, phi_c, theta_c), &[row]);
        }
    }

    #[test]
    fn test_szp_singularity_errors() {
        let native_pole = NativeCoord::new(0.0, -90.0);
        let native_normal = NativeCoord::new(0.0, 45.0);
        let inter_normal = IntermediateCoord::new(10.0, 10.0);
        let inter_far = IntermediateCoord::new(500.0, 500.0);

        let err = Projection::szp(1.0, 0.0, -90.0)
            .project(native_normal)
            .unwrap_err();
        assert!(
            matches!(err, WcsError::Singularity { .. }),
            "zp_singularity_project: {:?}",
            err
        );

        let err = Projection::szp(1.0, 0.0, 90.0)
            .project(native_pole)
            .unwrap_err();
        assert!(
            matches!(err, WcsError::Singularity { .. }),
            "denominator_singularity: {:?}",
            err
        );

        let err = Projection::szp(1.0, 0.0, -90.0)
            .deproject(inter_normal)
            .unwrap_err();
        assert!(
            matches!(err, WcsError::Singularity { .. }),
            "zp_singularity_deproject: {:?}",
            err
        );

        let err = Projection::szp(2.0, 45.0, 60.0)
            .deproject(inter_far)
            .unwrap_err();
        assert!(
            matches!(err, WcsError::OutOfBounds { .. }),
            "negative_discriminant: {:?}",
            err
        );
    }

    #[test]
    fn test_szp_rejects_points_beyond_the_boundary() {
        // With θc = 90° these are the AZP cases: the limb at −30° for µ = 2 and
        // 19.47° for µ = −3, the divergence at −30° for µ = 0.5. For µ = 2 and
        // (φc, θc) = (180°, 65°), eq. 48 puts the limb at −5° on φ = 0° and −55°
        // on φ = 180°.
        let cases = [
            (2.0, 0.0, 90.0, 37.5, -25.0, true),
            (2.0, 0.0, 90.0, 37.5, -40.0, false),
            (-3.0, 0.0, 90.0, 37.5, 25.0, true),
            (-3.0, 0.0, 90.0, 37.5, 15.0, false),
            (0.5, 0.0, 90.0, 37.5, -25.0, true),
            (0.5, 0.0, 90.0, 37.5, -40.0, false),
            (2.0, 180.0, 65.0, 0.0, -4.0, true),
            (2.0, 180.0, 65.0, 0.0, -6.0, false),
            (2.0, 180.0, 65.0, 180.0, -54.0, true),
            (2.0, 180.0, 65.0, 180.0, -56.0, false),
            (2.0, 180.0, 65.0, 0.0, -90.0, false),
        ];
        for (mu, phi_c, theta_c, phi, theta, inside) in cases {
            let native = NativeCoord::new(phi, theta);
            let got = Projection::szp(mu, phi_c, theta_c).project(native);
            let label = format!("({mu}, {phi_c}, {theta_c}) at ({phi}, {theta}): {got:?}");
            assert_eq!(got.is_ok(), inside, "{label}");
            assert!(
                inside || matches!(got, Err(WcsError::OutOfBounds { .. })),
                "{label}"
            );
        }
    }
}
