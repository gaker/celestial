use celestial_core::constants::RAD_TO_DEG;

use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::newton::NewtonConfig;
use crate::trig::{atan2d, sincosd};

#[inline]
pub(crate) fn pole_native_coord() -> NativeCoord {
    NativeCoord::new(0.0, 90.0)
}

// Paper II eqs. 12 and 13, with Rθ and φ in degrees.
#[inline]
pub(crate) fn radial_to_intermediate(r_theta: f64, phi: f64) -> IntermediateCoord {
    let (sin_phi, cos_phi) = sincosd(phi);
    IntermediateCoord::new(r_theta * sin_phi, -r_theta * cos_phi)
}

// Every deprojection ends here. A point beyond a projection's boundary, off
// the end of a cylinder, past the edge of a bounded map or behind a cone's
// apex, comes back with θ or φ off the native sphere (Paper II §5), and
// wrapping φ would hide that. A point on the boundary can come back just past
// it, so the check allows the tolerance the iterative deprojections converge
// to. It is measured as an arc on the sphere: past the seam that is
// (|φ| − 180°) cos θ, which lets any φ through at a pole, where all are one.
pub(crate) fn native_coord(phi: f64, theta: f64) -> WcsResult<NativeCoord> {
    const BOUNDARY_TOL: f64 = NewtonConfig::DEFAULT_TOL * RAD_TO_DEG;
    if theta.is_nan() {
        return Err(WcsError::invalid_parameter("latitude must be finite"));
    }
    let clamped = theta.clamp(-90.0, 90.0);
    let past_pole = libm::fabs(theta) - 90.0;
    let past_seam = (libm::fabs(phi) - 180.0) * sincosd(clamped).1;
    if past_pole > BOUNDARY_TOL || past_seam > BOUNDARY_TOL {
        return Err(WcsError::out_of_bounds(
            "Point outside the projection boundary",
        ));
    }
    Ok(NativeCoord::new(wrap_longitude(phi)?, clamped))
}

// No equation in Papers I or II has an image for a NaN or infinite value.
pub(crate) fn require_finite(values: &[f64], what: &str) -> WcsResult<()> {
    if values.iter().all(|v| v.is_finite()) {
        return Ok(());
    }
    Err(WcsError::invalid_parameter(format!(
        "{what} must be finite"
    )))
}

// Paper II §2.1 puts every latitude, native or celestial, on [−90°, 90°].
pub(crate) fn require_on_sphere(lon: f64, lat: f64) -> WcsResult<()> {
    require_finite(&[lon, lat], "coordinates")?;
    if libm::fabs(lat) > 90.0 {
        return Err(WcsError::out_of_bounds("latitude lies beyond a pole"));
    }
    Ok(())
}

// Paper II §7.3.4 returns longitudes in the normal range [−180°, 180°], whose
// ends are one meridian; +180° stands for it. The remainder is exact.
pub(crate) fn wrap_longitude(deg: f64) -> WcsResult<f64> {
    if !deg.is_finite() {
        return Err(WcsError::invalid_parameter("longitude must be finite"));
    }
    let wrapped = libm::remainder(deg, 360.0);
    Ok(if wrapped == -180.0 { 180.0 } else { wrapped })
}

// (φ, θ) → (x, y) → (φ, θ): the projection of the first pair and the
// deprojection of the second, each pinned one way.
#[cfg(test)]
pub(crate) type PinnedRow = ((f64, f64), (f64, f64), (f64, f64));

#[cfg(test)]
pub(crate) fn assert_pinned(proj: &crate::spherical::Projection, rows: &[PinnedRow]) {
    for &(native, inter, back) in rows {
        let got = proj.project(NativeCoord::new(native.0, native.1)).unwrap();
        assert_eq!((got.x_deg(), got.y_deg()), inter, "project {native:?}");
        let got = proj
            .deproject(IntermediateCoord::new(inter.0, inter.1))
            .unwrap();
        assert_eq!(
            (got.phi_deg(), got.theta_deg()),
            back,
            "deproject {inter:?}"
        );
    }
}

// A dimensionless quantity of order 1 closer to 0 than this counts as 0 where
// it is a divisor or decides a degenerate case. It is about 4.5 ε, the size of
// the rounding such a quantity carries.
pub(crate) const ZERO_LIMIT: f64 = 1e-15;

// The smallest |θa| (conics) or |θ1| (BON), in degrees, that
// check_nonzero_param accepts. It is not an accuracy bound: Y0 grows as
// cot θa, and y = Y0 − Rθ cos Cφ loses ε Y0, about 0.007° at this limit.
const MIN_PARALLEL_DEG: f64 = 1e-10;

#[inline]
pub(crate) fn check_nonzero_param(value: f64, context: &str) -> WcsResult<()> {
    if value.abs() < MIN_PARALLEL_DEG {
        return Err(WcsError::invalid_parameter(format!(
            "{}: parameter cannot be zero",
            context
        )));
    }
    Ok(())
}

// Paper II eqs. 14 and 15, in degrees. φ has no value at the origin; 0 stands
// for it.
#[inline]
pub(crate) fn intermediate_to_polar(x: f64, y: f64) -> (f64, f64, bool) {
    let r_theta = libm::sqrt(x * x + y * y);
    let is_pole = r_theta == 0.0;
    let phi = if is_pole { 0.0 } else { atan2d(x, -y) };
    (phi, r_theta, is_pole)
}
// A 2×2 determinant is indistinguishable from zero once it is no larger than
// the rounding error of ad − bc, about ε(|ad| + |bc|). An absolute threshold
// would instead reject well-conditioned matrices whose entries are small.
pub(crate) fn is_singular(m: [[f64; 2]; 2], det: f64) -> bool {
    let scale = libm::fabs(m[0][0] * m[1][1]) + libm::fabs(m[0][1] * m[1][0]);
    !det.is_finite() || libm::fabs(det) <= f64::EPSILON * scale
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::SQRT2;

    #[test]
    fn test_pole_native_coord() {
        let pole = pole_native_coord();
        assert_eq!(pole.phi_deg(), 0.0);
        assert_eq!(pole.theta_deg(), 90.0);
    }

    #[test]
    fn test_radial_to_intermediate() {
        // At the origin the projection is (0, 0); the assertion is byte-exact.
        let zero = radial_to_intermediate(0.0, 0.0);
        assert_eq!(zero.x_deg(), 0.0);
        assert_eq!(zero.y_deg(), 0.0);

        let inter = radial_to_intermediate(2.0, 30.0);
        assert_eq!(inter.x_deg(), 1.0);
        assert_eq!(inter.y_deg(), -2.0 * sincosd(30.0).1);
    }

    // φ = ±180° is one meridian, and +180° stands for it.
    #[test]
    fn test_native_coord_wraps_phi_at_the_seam() {
        let phi_at = |phi: f64, theta: f64| native_coord(phi, theta).unwrap().phi_deg();
        assert_eq!(phi_at(180.0, 0.0), 180.0);
        assert_eq!(phi_at(-180.0, 0.0), 180.0);
        assert_eq!(phi_at(270.0, 90.0), -90.0);
    }

    #[test]
    fn test_wrap_longitude_to_normal_range() {
        let cases = [
            (180.0, 180.0),
            (-180.0, 180.0),
            (540.0, 180.0),
            (-179.5, -179.5),
            (359.5, -0.5),
            (-360.0, 0.0),
            (123.456 + 360.0, 123.45600000000002),
        ];
        for (deg, expected) in cases {
            assert_eq!(wrap_longitude(deg).unwrap(), expected, "{deg}");
        }
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(wrap_longitude(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn test_native_coord_rejects_points_off_the_sphere() {
        let past = 1e-9;
        for (phi, theta) in [(180.0 + past, 0.0), (0.0, 90.0 + past), (-185.0, -86.0)] {
            let err = native_coord(phi, theta).unwrap_err();
            let out = matches!(err, WcsError::OutOfBounds { .. });
            assert!(out, "({phi}, {theta}): {err:?}");
        }
    }

    #[test]
    fn test_native_coord_keeps_points_on_the_boundary() {
        let next = |x: f64| f64::from_bits(x.to_bits() + 1);
        let native = native_coord(0.0, -next(90.0)).unwrap();
        assert_eq!(native.theta_deg(), -90.0);
        let native = native_coord(next(180.0), 0.0).unwrap();
        assert_eq!(native.phi_deg(), -179.99999999999997);
    }

    #[test]
    fn test_native_coord_rejects_non_finite_coordinates() {
        assert!(native_coord(f64::NAN, 0.0).is_err());
        assert!(native_coord(f64::INFINITY, 0.0).is_err());
        assert!(native_coord(0.0, f64::NAN).is_err());
        assert!(native_coord(0.0, f64::INFINITY).is_err());
    }

    #[test]
    fn test_check_nonzero_param() {
        assert!(check_nonzero_param(1.0, "test").is_ok());
        assert!(check_nonzero_param(-0.5, "test").is_ok());
        assert!(check_nonzero_param(0.0, "test").is_err());
        assert!(check_nonzero_param(1e-11, "test").is_err());
    }

    #[test]
    fn test_intermediate_to_polar() {
        // Origin is flagged as the pole.
        let (phi, r, is_pole) = intermediate_to_polar(0.0, 0.0);
        assert_eq!(phi, 0.0);
        assert_eq!(r, 0.0);
        assert!(is_pole);

        let (phi, r, is_pole) = intermediate_to_polar(1.0, -1.0);
        assert_eq!((phi, r), (45.0, SQRT2));
        assert!(!is_pole);
    }
}
