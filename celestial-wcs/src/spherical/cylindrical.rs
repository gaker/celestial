use crate::common::native_coord;
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{asind, atan2d, atand, sincosd, sind, tand, to_degrees, to_radians};

// Paper II eqs. 84 and 85.
pub(crate) fn project_car(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    Ok(IntermediateCoord::new(native.phi_deg(), native.theta_deg()))
}

pub(crate) fn deproject_car(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let phi = cylinder_longitude(inter.x_deg())?;
    let theta = inter.y_deg();
    if theta.abs() > 90.0 {
        return Err(WcsError::out_of_bounds(
            "CAR deprojection: point beyond the poles",
        ));
    }
    Ok(NativeCoord::new(phi, theta))
}

// φ is x itself, so the seam at x = ±180 is exact, and it is returned as
// +180° as everywhere else. Beyond it lies off the end of the cylinder.
fn cylinder_longitude(x: f64) -> WcsResult<f64> {
    if x.abs() > 180.0 {
        return Err(WcsError::out_of_bounds("Point beyond the cylinder's seam"));
    }
    Ok(if x == -180.0 { 180.0 } else { x })
}

// Eqs. 86 and 87. y diverges at the poles and nowhere else.
pub(crate) fn project_mer(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    if theta.abs() == 90.0 {
        return Err(WcsError::singularity(
            "MER projection undefined at theta = +/-90",
        ));
    }
    let y = to_degrees(libm::log(tand((90.0 + theta) / 2.0)));
    Ok(IntermediateCoord::new(native.phi_deg(), y))
}

// Eq. 88.
pub(crate) fn deproject_mer(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let theta = 2.0 * atand(libm::exp(to_radians(inter.y_deg()))) - 90.0;
    Ok(NativeCoord::new(cylinder_longitude(inter.x_deg())?, theta))
}

// Eqs. 79 and 80.
pub(crate) fn project_cea(native: NativeCoord, lambda: f64) -> WcsResult<IntermediateCoord> {
    let y = to_degrees(sind(native.theta_deg()) / lambda);
    Ok(IntermediateCoord::new(native.phi_deg(), y))
}

// Eqs. 81 and 82.
pub(crate) fn deproject_cea(inter: IntermediateCoord, lambda: f64) -> WcsResult<NativeCoord> {
    let sin_theta = lambda * to_radians(inter.y_deg());
    if sin_theta.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "CEA deprojection: |lambda * y| > 1",
        ));
    }
    Ok(NativeCoord::new(
        cylinder_longitude(inter.x_deg())?,
        asind(sin_theta),
    ))
}

// Eqs. 74 and 75. For −1 ≤ μ ≤ 0 the parallel where μ + cos θ = 0 projects to
// infinity.
pub(crate) fn project_cyp(
    native: NativeCoord,
    mu: f64,
    lambda: f64,
) -> WcsResult<IntermediateCoord> {
    let (sin_theta, cos_theta) = sincosd(native.theta_deg());
    let denom = mu + cos_theta;
    if denom == 0.0 {
        return Err(WcsError::singularity(
            "CYP projection singularity: mu + cos(theta) = 0",
        ));
    }
    check_cyp_domain(mu, cos_theta, denom)?;
    let y = to_degrees((mu + lambda) / denom * sin_theta);
    Ok(IntermediateCoord::new(lambda * native.phi_deg(), y))
}

// dy/dθ is a constant times (μ + λ)(μ cos θ + 1)/(μ + cos θ)², so y is
// monotonic on the equator's side of the parallel where μ + cos θ = 0 and, for
// μ < −1, of the limb where μ cos θ + 1 = 0. Beyond either, y takes a value
// that a point nearer the equator also has, and the inverse returns that point.
fn check_cyp_domain(mu: f64, cos_theta: f64, denom: f64) -> WcsResult<()> {
    let beyond_divergence = (mu + 1.0) * denom < 0.0;
    let beyond_limb = (mu * cos_theta + 1.0) * (mu + 1.0) < 0.0;
    if beyond_divergence || beyond_limb {
        return Err(WcsError::out_of_bounds(
            "CYP projection: point beyond the boundary",
        ));
    }
    Ok(())
}

// Eqs. 76 to 78.
pub(crate) fn deproject_cyp(
    inter: IntermediateCoord,
    mu: f64,
    lambda: f64,
) -> WcsResult<NativeCoord> {
    let eta = to_radians(inter.y_deg()) / (mu + lambda);
    if mu == -1.0 {
        return deproject_cyp_on_sphere(inter.x_deg() / lambda, eta);
    }
    let sin_offset = eta * mu / libm::sqrt(eta * eta + 1.0);
    if sin_offset.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "CYP deprojection: point outside valid region",
        ));
    }
    let theta = atan2d(eta, 1.0) + asind(sin_offset);
    native_coord(inter.x_deg() / lambda, theta)
}

// With μ = −1 the two terms of eq. 77 cancel for every η and give the trivial
// root θ = 0. Eq. 75 there reduces to η = −cot(θ/2), which is one-to-one on
// θ ≠ 0 and reaches |θ| ≤ 90° only for |η| ≥ 1.
fn deproject_cyp_on_sphere(phi: f64, eta: f64) -> WcsResult<NativeCoord> {
    if libm::fabs(eta) < 1.0 {
        return Err(WcsError::out_of_bounds(
            "CYP deprojection: point outside valid region",
        ));
    }
    native_coord(phi, -2.0 * atand(1.0 / eta))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // Per-projection native_reference checks are covered by
    // spherical::tests::test_all_projections_map_reference_to_origin_and_back.

    #[test]
    fn test_car_known_values_and_roundtrip() {
        // CAR is a literal pass-through (phi, theta) -> (x, y) in degrees.
        let proj = Projection::car();

        // Spec-anchored: a couple of known points should be byte-exact.
        let native = NativeCoord::new(90.0, 45.0);
        let inter = proj.project(native).unwrap();
        assert_eq!(inter.x_deg(), 90.0);
        assert_eq!(inter.y_deg(), 45.0);

        let native2 = NativeCoord::new(-120.0, -60.0);
        let inter2 = proj.project(native2).unwrap();
        assert_eq!((inter2.x_deg(), inter2.y_deg()), (-120.0, -60.0));

        // CAR is an identity in degrees, so the roundtrip must be byte-exact.
        // φ = −180° comes back as +180°, the same meridian.
        for phi_deg in [-179.5, -90.0, 0.0, 45.0, 90.0, 135.0, 180.0] {
            for theta_deg in [-85.0, -45.0, 0.0, 45.0, 85.0] {
                let original = NativeCoord::new(phi_deg, theta_deg);
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert_eq!(original.phi_deg(), recovered.phi_deg());
                assert_eq!(original.theta_deg(), recovered.theta_deg());
            }
        }
    }

    #[test]
    fn test_mer_known_value_and_pins() {
        let proj = Projection::mer();

        // Eq. 87 at θ = 45: (180°/π) ln tan 67.5° = 50.498986710526211442…
        // This is 1.1 ulp low: the logarithm turns the rounding of tan 67.5°
        // into an absolute error in y.
        let inter = proj.project(NativeCoord::new(0.0, 45.0)).unwrap();
        assert_eq!((inter.x_deg(), inter.y_deg()), (0.0, 50.498986710526204));

        // Against a 60-digit evaluation y is within 1.11 ulp and θ within
        // 0.78 ulp.
        let rows = [
            (
                (-180.0, -45.0),
                (-180.0, -50.498986710526204),
                (180.0, -44.99999999999999),
            ),
            ((-180.0, 45.0), (-180.0, 50.498986710526204), (180.0, 45.0)),
            ((45.0, 80.0), (45.0, 139.58661673332185), (45.0, 80.0)),
            ((135.0, -80.0), (135.0, -139.58661673332185), (135.0, -80.0)),
        ];
        assert_pinned(&proj, &rows);
    }

    #[test]
    fn test_mer_reference_point_is_exact() {
        // Eq. 87 at θ = 0 is ln tan 45° = 0.
        let inter = Projection::mer()
            .project(NativeCoord::new(30.0, 0.0))
            .unwrap();
        assert_eq!((inter.x_deg(), inter.y_deg()), (30.0, 0.0));
        let native = Projection::mer()
            .deproject(IntermediateCoord::new(30.0, 0.0))
            .unwrap();
        assert_eq!((native.phi_deg(), native.theta_deg()), (30.0, 0.0));
    }

    #[test]
    fn test_cea_known_value_and_pins() {
        let proj = Projection::cea();

        // Eq. 80 with λ = 1 at θ = 30: 90°/π = 28.647889756541160438…, which
        // rounds to this value (0.28 ulp away).
        let inter = proj.project(NativeCoord::new(90.0, 30.0)).unwrap();
        assert_eq!((inter.x_deg(), inter.y_deg()), (90.0, 28.64788975654116));

        // Against a 60-digit evaluation y is within 0.29 ulp, so correctly
        // rounded, and θ within 1.03 ulp (at ±85°, where sin⁻¹ steepens).
        let rows = [
            (
                (-180.0, -85.0),
                (-180.0, -57.07775177396626),
                (180.0, -85.0),
            ),
            (
                (-180.0, -45.0),
                (-180.0, -40.51423422706977),
                (180.0, -44.99999999999999),
            ),
            ((135.0, 85.0), (135.0, 57.07775177396626), (135.0, 85.0)),
        ];
        assert_pinned(&proj, &rows);
        let half = (
            (60.0, 45.0),
            (60.0, 81.02846845413954),
            (60.0, 44.99999999999999),
        );
        assert_pinned(&Projection::cea_with_lambda(0.5), &[half]);
    }

    // (µ, λ) = (1, 1) at three points, then three other (µ, λ) at one point.
    // Against a 60-digit evaluation y is within 1.56 ulp and θ within 1.38
    // ulp.
    #[test]
    fn test_cyp_pinned_both_ways() {
        let rows = [
            (
                (-180.0, -60.0),
                (-180.0, -66.15946745061504),
                (180.0, -60.0),
            ),
            (
                (60.0, 45.0),
                (60.0, 47.465377882114446),
                (60.0, 44.999999999999986),
            ),
            ((135.0, 30.0), (135.0, 30.704715700484144), (135.0, 30.0)),
        ];
        assert_pinned(&Projection::cyp(1.0, 1.0), &rows);
        let cases = [
            (
                (0.5, 2.0),
                ((60.0, 45.0), (120.0, 83.90772643006274), (60.0, 45.0)),
            ),
            (
                (2.0, 0.5),
                ((60.0, 45.0), (30.0, 37.41469906971313), (60.0, 45.0)),
            ),
            (
                (5.0, 0.5),
                (
                    (60.0, 45.0),
                    (30.0, 39.04400194218131),
                    (60.0, 44.999999999999986),
                ),
            ),
        ];
        for ((mu, lambda), row) in cases {
            assert_pinned(&Projection::cyp(mu, lambda), &[row]);
        }
    }

    // The seam is one meridian, which every projection returns as +180°.
    #[test]
    fn test_seam_deprojects_to_plus_180() {
        for proj in [
            Projection::car(),
            Projection::mer(),
            Projection::cea(),
            Projection::cyp(1.0, 1.0),
        ] {
            let native = proj.deproject(IntermediateCoord::new(-180.0, 0.0)).unwrap();
            assert_eq!(native, NativeCoord::new(180.0, 0.0), "{proj:?}");
        }
    }

    // Eq. 75 keeps the sign of y that μ + 1 gives it near the equator only up
    // to the parallel where μ + cos θ = 0. For μ < −1 it turns back at the limb
    // cos θ = −1/μ. Points past either fold onto other points.
    #[test]
    fn test_cyp_refuses_points_that_fold_over() {
        let refused = [
            (Projection::cyp(-0.5, 1.0), -85.93, -85.27),
            (Projection::cyp(-0.5, 1.0), 0.0, 90.0),
            (Projection::cyp(-0.5, 0.2), 10.0, 70.0),
            (Projection::cyp(-2.0, 1.0), 10.0, 61.0),
            (Projection::cyp(-2.0, 1.0), 10.0, -80.0),
        ];
        for (proj, phi, theta) in refused {
            let result = proj.project(NativeCoord::new(phi, theta));
            let out = matches!(result, Err(WcsError::OutOfBounds { .. }));
            assert!(out, "{proj:?} ({phi}, {theta}): {result:?}");
        }
        let kept = [
            (Projection::cyp(-0.5, 0.2), 10.0, 20.0),
            (Projection::cyp(-2.0, 1.0), 10.0, 59.0),
            (Projection::cyp(-2.0, 1.0), 10.0, -60.0),
        ];
        for (proj, phi, theta) in kept {
            let result = proj.project(NativeCoord::new(phi, theta));
            assert!(result.is_ok(), "{proj:?} ({phi}, {theta}): {result:?}");
        }
    }

    // With µ = −1 the projection point is on the sphere and eq. 77 gives only
    // θ = 0, so the inverse solves η = −cot(θ/2) instead. Against 60 digits
    // θ is within 0.38 ulp. y at ±40° is 1.74 ulp off: µ + cos θ = cos 40° − 1
    // cancels to −0.234, which magnifies the rounding of cos 40° about 3.3
    // times.
    #[test]
    fn test_cyp_projection_point_on_sphere() {
        let rows = [
            (
                (10.0, 40.0),
                (20.0, -157.41886044224438),
                (10.0, 40.000000000000014),
            ),
            (
                (10.0, -40.0),
                (20.0, 157.41886044224438),
                (10.0, -40.000000000000014),
            ),
            ((10.0, 90.0), (20.0, -57.29577951308232), (10.0, 90.0)),
        ];
        assert_pinned(&Projection::cyp(-1.0, 2.0), &rows);
        let rows = [
            (
                (60.0, 40.0),
                (30.0, 78.70943022112219),
                (60.0, 40.000000000000014),
            ),
            (
                (60.0, -40.0),
                (30.0, -78.70943022112219),
                (60.0, -40.000000000000014),
            ),
            ((60.0, 90.0), (30.0, 28.64788975654116), (60.0, 90.0)),
        ];
        assert_pinned(&Projection::cyp(-1.0, 0.5), &rows);
    }

    // |η| < 1 is cot(θ/2) for some |θ| > 90°, beyond the poles.
    #[test]
    fn test_cyp_on_sphere_refuses_beyond_poles() {
        let proj = Projection::cyp(-1.0, 2.0);
        for y in [10.0, 0.0, -57.0] {
            let result = proj.deproject(IntermediateCoord::new(0.0, y));
            let out = matches!(result, Err(WcsError::OutOfBounds { .. }));
            assert!(out, "y = {y}: {result:?}");
        }
    }

    #[test]
    fn test_cylindrical_error_cases() {
        // MER blows up at the poles (theta = +/- 90).
        let mer = Projection::mer();
        for theta in [90.0, -90.0] {
            let native = NativeCoord::new(0.0, theta);
            assert!(
                mer.project(native).is_err(),
                "MER theta = {} should error",
                theta
            );
        }

        // CEA deproject: |lambda * y| > 1 is out of range.
        let cea = Projection::cea();
        assert!(cea.deproject(IntermediateCoord::new(0.0, 100.0)).is_err());

        // With μ = −0.5 the parallels θ = ±60° (μ + cos θ = 0) project to
        // infinity.
        let cyp = Projection::cyp(-0.5, 1.0);
        for theta in [60.0, -60.0] {
            assert!(cyp.project(NativeCoord::new(0.0, theta)).is_err());
        }

        // μ = 2 puts the edge of the map at |y| = (μ + λ)/√(μ² − 1) rad,
        // 99.2°.
        let edge = Projection::cyp(2.0, 1.0).deproject(IntermediateCoord::new(0.0, 200.0));
        assert!(
            matches!(edge, Err(WcsError::OutOfBounds { .. })),
            "{edge:?}"
        );

        // CYP deproject with lambda = 0 is an invalid parameter.
        let cyp_zero = Projection::cyp(1.0, 0.0);
        assert!(cyp_zero
            .deproject(IntermediateCoord::new(10.0, 10.0))
            .is_err());
    }
}
