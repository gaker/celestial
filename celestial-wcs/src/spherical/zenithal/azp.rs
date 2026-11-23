use celestial_core::constants::RAD_TO_DEG;

use crate::common::{native_coord, pole_native_coord, ZERO_LIMIT};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{asind, atan2d, sincosd, to_degrees};

// Paper II eqs. 20-22. With γ = 0 they reduce exactly to eq. 16.
pub(crate) fn project_azp(
    native: NativeCoord,
    mu: f64,
    (sin_gamma, cos_gamma): (f64, f64),
) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    let (sin_theta, cos_theta) = sincosd(theta);
    let (sin_phi, cos_phi) = sincosd(native.phi_deg());
    let denom = mu + sin_theta + cos_theta * cos_phi * (sin_gamma / cos_gamma);
    check_azp_domain(mu, sin_theta, denom)?;
    if theta == 90.0 {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }
    let r = to_degrees((mu + 1.0) * cos_theta / denom);
    Ok(IntermediateCoord::new(
        r * sin_phi,
        -r * cos_phi / cos_gamma,
    ))
}

// The inverse (eq. 17) returns the planeward of the two points on a ray from
// the point of projection, and when |µ| > 1 the limb at θx = sin⁻¹(−1/µ)
// (eq. 32) divides them. The projection diverges where the denominator of
// eq. 22 passes through zero, and beyond that it has lost the sign µ + 1 it
// has at the pole.
fn check_azp_domain(mu: f64, sin_theta: f64, denom: f64) -> WcsResult<()> {
    if denom.abs() < ZERO_LIMIT {
        return Err(WcsError::singularity("AZP projection diverges"));
    }
    let beyond_limb = mu.abs() > 1.0 && sin_theta < -1.0 / mu;
    if beyond_limb || (mu + 1.0) * denom < 0.0 {
        return Err(WcsError::out_of_bounds(
            "AZP projection: point beyond the boundary",
        ));
    }
    Ok(())
}

// Eqs. 23-28.
pub(crate) fn deproject_azp(
    inter: IntermediateCoord,
    mu: f64,
    (sin_gamma, cos_gamma): (f64, f64),
) -> WcsResult<NativeCoord> {
    let (x, y) = (inter.x_deg(), inter.y_deg());
    if x == 0.0 && y == 0.0 {
        return Ok(pole_native_coord());
    }
    let phi = atan2d(x, -y * cos_gamma);
    let r = libm::sqrt(x * x + (y * cos_gamma).powi(2));
    // The denominator is in degrees.
    let denom = to_degrees(mu + 1.0) + y * sin_gamma;
    if denom.abs() < ZERO_LIMIT * RAD_TO_DEG {
        return Err(outside_azp());
    }
    native_coord(phi, azp_latitude(r / denom, mu)?)
}

// Eq. 24 as θ = ψ − ω, the planeward solution.
fn azp_latitude(rho: f64, mu: f64) -> WcsResult<f64> {
    let s = rho * mu / libm::sqrt(rho * rho + 1.0);
    if s.abs() > 1.0 {
        return Err(outside_azp());
    }
    Ok(atan2d(1.0, rho) - asind(s))
}

fn outside_azp() -> WcsError {
    WcsError::out_of_bounds("Point outside AZP projection boundary")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // (µ, γ) across the µ range and one slant. Against a 60-digit evaluation
    // the projected values are within 2.28 ulp (y for µ = 0) and the
    // deprojected ones within 1.00 ulp.
    #[test]
    fn test_azp_pinned_both_ways() {
        let cases = [
            (
                (0.0, 0.0),
                (
                    (60.0, 45.0),
                    (49.619600587961294, -28.64788975654117),
                    (59.99999999999999, 44.99999999999999),
                ),
            ),
            (
                (2.0, 0.0),
                (
                    (60.0, 45.0),
                    (38.88249584318589, -22.448819441827883),
                    (59.99999999999999, 45.0),
                ),
            ),
            (
                (10.0, 0.0),
                (
                    (60.0, 45.0),
                    (36.0461443504815, -20.811251143998604),
                    (59.99999999999999, 45.000000000000014),
                ),
            ),
            (
                (2.0, 30.0),
                (
                    (30.0, 70.0),
                    (9.4494617415177, -18.8989234830354),
                    (30.000000000000004, 70.0),
                ),
            ),
        ];
        for ((mu, gamma), row) in cases {
            assert_pinned(&Projection::azp(mu, gamma), &[row]);
        }
    }

    // µ + sin θ is 1.1e-12 here: near the divergence but not on it, so the
    // point projects, far out.
    #[test]
    fn test_azp_projects_next_to_its_divergence() {
        let native = NativeCoord::new(0.0, -29.99999999993);
        let inter = Projection::azp(0.5, 0.0).project(native).unwrap();
        assert!(inter.y_deg() < -1e12, "{inter:?}");
    }

    #[test]
    fn test_azp_error_cases() {
        // Project: µ + sin θ = 0 (non-slant denom singularity)
        let native = NativeCoord::new(0.0, -30.0);
        let err = Projection::azp(0.5, 0.0).project(native).unwrap_err();
        assert!(
            matches!(err, WcsError::Singularity { .. }),
            "non-slant: {:?}",
            err
        );

        // Project: slant denom singularity
        let native = NativeCoord::new(0.0, -45.0);
        let err = Projection::azp(0.0, 45.0).project(native).unwrap_err();
        assert!(
            matches!(err, WcsError::Singularity { .. }),
            "slant: {:?}",
            err
        );

        // Deproject: each way past the boundary. With µ = 0 and γ = 30° the
        // denominator of eq. 25 is (180°/π) + y/2, exactly zero at this y.
        let cases = &[
            (10.0, 0.0, 50.0, 50.0),               // no slant, far out
            (0.0, 30.0, 0.0, -114.59155902616465), // denominator zero
            (10.0, 30.0, 80.0, 80.0),              // s out of bounds with slant
        ];
        for &(mu, gamma, x, y) in cases {
            let inter = IntermediateCoord::new(x, y);
            let err = Projection::azp(mu, gamma).deproject(inter).unwrap_err();
            assert!(
                matches!(err, WcsError::OutOfBounds { .. }),
                "(mu={mu}, gamma={gamma}, x={x}, y={y}): got: {err:?}"
            );
        }
    }

    #[test]
    fn test_azp_rejects_points_beyond_the_boundary() {
        // The limb is at θx = sin⁻¹(−1/µ) (eq. 32): −30° for µ = 2, 19.47° for
        // µ = −3. For µ = 0.5 the projection diverges at sin⁻¹(−µ) = −30°;
        // tilted by 25° it diverges at −51.94° on φ = 0° and −1.94° on φ = 180°,
        // and µ + sin θ = 0 is no longer special.
        let cases = [
            (2.0, 0.0, 37.5, -25.0, true),
            (2.0, 0.0, 37.5, -40.0, false),
            (-3.0, 0.0, 37.5, 25.0, true),
            (-3.0, 0.0, 37.5, 15.0, false),
            (0.5, 0.0, 37.5, -25.0, true),
            (0.5, 0.0, 37.5, -40.0, false),
            (2.0, 25.0, 0.0, -90.0, false),
            (0.5, 25.0, 0.0, -30.0, true),
            (0.5, 25.0, 0.0, -50.0, true),
            (0.5, 25.0, 0.0, -55.0, false),
            (0.5, 25.0, 180.0, -1.0, true),
            (0.5, 25.0, 180.0, -3.0, false),
        ];
        for (mu, gamma, phi, theta, inside) in cases {
            let native = NativeCoord::new(phi, theta);
            let got = Projection::azp(mu, gamma).project(native);
            let label = format!("(µ, γ) = ({mu}, {gamma}) at ({phi}, {theta}): {got:?}");
            assert_eq!(got.is_ok(), inside, "{label}");
            assert!(
                inside || matches!(got, Err(WcsError::OutOfBounds { .. })),
                "{label}"
            );
        }
    }
}
