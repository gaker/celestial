use celestial_core::constants::{DEG_TO_RAD, HALF_PI, RAD_TO_DEG};

use crate::common::{
    native_coord_from_radians, pole_native_coord,
    radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

pub(crate) fn project_azp(
    native: NativeCoord,
    mu: f64,
    gamma_deg: f64,
) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let gamma = gamma_deg * DEG_TO_RAD;

    let (sin_theta, cos_theta) = libm::sincos(theta);

    let denom = mu + sin_theta;
    if denom.abs() < 1e-10 {
        return Err(WcsError::singularity(
            "AZP projection singularity: mu + sin(theta) = 0",
        ));
    }

    if gamma_deg.abs() < 1e-10 {
        if theta == HALF_PI {
            return Ok(IntermediateCoord::new(0.0, 0.0));
        }
        let r_theta = (mu + 1.0) * cos_theta / denom;
        Ok(radial_to_intermediate(r_theta, phi))
    } else {
        let (sin_gamma, cos_gamma) = libm::sincos(gamma);
        let tan_gamma = sin_gamma / cos_gamma;

        let denom_full = denom + cos_theta * libm::cos(phi) * tan_gamma;
        if denom_full.abs() < 1e-10 {
            return Err(WcsError::singularity("AZP slant projection singularity"));
        }

        let r = (mu + 1.0) * cos_theta / denom_full;
        let (ps, pc) = libm::sincos(phi);
        let x = r * ps * RAD_TO_DEG;
        let y = -r * pc / cos_gamma * RAD_TO_DEG;
        Ok(IntermediateCoord::new(x, y))
    }
}

pub(crate) fn deproject_azp(
    inter: IntermediateCoord,
    mu: f64,
    gamma_deg: f64,
) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;

    if x == 0.0 && y == 0.0 {
        return Ok(pole_native_coord());
    }

    if gamma_deg.abs() < 1e-10 {
        let r_theta = libm::sqrt(x * x + y * y);
        let phi = libm::atan2(x, -y);
        let rho = r_theta / (mu + 1.0);
        let s = rho * mu / libm::sqrt(rho * rho + 1.0);
        if s.abs() > 1.0 {
            return Err(WcsError::out_of_bounds(
                "Point outside AZP projection boundary",
            ));
        }
        let theta = libm::atan2(1.0_f64, rho) - libm::asin(s);
        Ok(native_coord_from_radians(phi, theta))
    } else {
        let gamma = gamma_deg * DEG_TO_RAD;
        let (sin_gamma, cos_gamma) = libm::sincos(gamma);

        let phi = libm::atan2(x, -y * cos_gamma);

        let r_theta = libm::sqrt(x * x + (y * cos_gamma).powi(2));

        let denom = (mu + 1.0) + y * sin_gamma;
        if denom.abs() < 1e-15 {
            return Err(WcsError::out_of_bounds(
                "Point outside AZP projection boundary",
            ));
        }
        let rho = r_theta / denom;

        let psi = libm::atan2(1.0_f64, rho);
        let s = rho * mu / libm::sqrt(rho * rho + 1.0);
        if s.abs() > 1.0 {
            return Err(WcsError::out_of_bounds(
                "Point outside AZP projection boundary",
            ));
        }
        let omega = libm::asin(s);

        let theta = psi - omega;

        Ok(native_coord_from_radians(phi, theta))
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;
    #[test]
    fn test_azp_roundtrip() {
        // (mu, gamma, phi, theta, ulp). Sweeps both slant (gamma != 0) and
        // non-slant configs across the mu range.
        let cases: &[(f64, f64, f64, f64, u64)] = &[
            (0.0, 0.0, 60.0, 45.0, 3),
            (1.0, 0.0, 60.0, 45.0, 3),
            (2.0, 0.0, 60.0, 45.0, 3),
            (5.0, 0.0, 60.0, 45.0, 3),
            (10.0, 0.0, 60.0, 45.0, 3),
            (2.0, 0.0, 45.0, 60.0, 2),
            (2.0, 30.0, 30.0, 70.0, 5),
        ];
        for (mu, gamma, phi, theta, ulp) in cases {
            let proj = Projection::azp(*mu, *gamma);
            let original = NativeCoord::new(Angle::from_degrees(*phi), Angle::from_degrees(*theta));
            let inter = proj.project(original).unwrap();
            let recovered = proj.deproject(inter).unwrap();
            assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), *ulp,
                "phi (mu={}, gamma={})", mu, gamma);
            assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), *ulp,
                "theta (mu={}, gamma={})", mu, gamma);
        }
    }

    #[test]
    fn test_azp_error_cases() {
        // Project: mu + sin(theta) = 0 (non-slant denom singularity)
        let mu: f64 = 0.5;
        let theta_singular = -mu.asin() * RAD_TO_DEG;
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta_singular));
        let err = project_azp(native, mu, 0.0).unwrap_err();
        assert!(matches!(err, WcsError::Singularity { .. }), "non-slant: {:?}", err);

        // Project: slant denom singularity
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-45.0));
        let err = project_azp(native, 0.0, 45.0).unwrap_err();
        assert!(matches!(err, WcsError::Singularity { .. }), "slant: {:?}", err);

        // Deproject: out-of-bounds / singularity cases
        let cases = &[
            (10.0, 0.0, 50.0, 50.0),               // no slant, far out
            (0.0, 90.0, 0.0, -57.29577951308232),  // denom zero with slant
            (10.0, 30.0, 80.0, 80.0),              // s out of bounds with slant
        ];
        for &(mu, gamma, x, y) in cases {
            let inter = IntermediateCoord::new(x, y);
            let err = deproject_azp(inter, mu, gamma).unwrap_err();
            assert!(
                matches!(err, WcsError::OutOfBounds { .. } | WcsError::Singularity { .. }),
                "(mu={}, gamma={}, x={}, y={}): got: {:?}", mu, gamma, x, y, err,
            );
        }
    }
}
