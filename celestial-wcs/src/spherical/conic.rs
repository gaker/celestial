use celestial_core::constants::{DEG_TO_RAD, HALF_PI};
use celestial_core::Angle;

use crate::common::{
    check_nonzero_param, deproject_conic_polar, native_coord_from_radians, project_conic_xy,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

#[inline]
fn standard_parallels(theta_a: f64, eta: f64) -> (f64, f64) {
    (theta_a - eta, theta_a + eta)
}

#[inline]
fn eta_cot_eta(eta: f64) -> f64 {
    if eta.abs() < 1e-12 {
        1.0
    } else {
        eta * libm::cos(eta) / libm::sin(eta)
    }
}

#[inline]
fn sin_eta_over_eta(eta: f64) -> f64 {
    if eta.abs() < 1e-12 {
        1.0
    } else {
        libm::sin(eta) / eta
    }
}

pub(crate) fn project_cop(
    native: NativeCoord,
    theta_a_deg: f64,
    eta_deg: f64,
) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_a = theta_a_deg * DEG_TO_RAD;
    let eta = eta_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COP projection: theta_a")?;

    let c = libm::sin(theta_a);
    let cos_eta = libm::cos(eta);
    let cot_theta_a = libm::cos(theta_a) / libm::sin(theta_a);

    let delta = theta - theta_a;
    let cos_delta = libm::cos(delta);
    if cos_delta.abs() < 1e-15 {
        return Err(WcsError::singularity(
            "COP projection: theta - theta_a = +/-90",
        ));
    }
    let tan_delta = libm::sin(delta) / cos_delta;

    let r_theta = cos_eta * (cot_theta_a - tan_delta);
    let y0 = cos_eta * cot_theta_a;

    Ok(project_conic_xy(r_theta, y0, c, phi))
}

pub(crate) fn deproject_cop(
    inter: IntermediateCoord,
    theta_a_deg: f64,
    eta_deg: f64,
) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_a = theta_a_deg * DEG_TO_RAD;
    let eta = eta_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COP projection: theta_a")?;

    let cos_eta = libm::cos(eta);
    if cos_eta.abs() < 1e-15 {
        return Err(WcsError::invalid_parameter(
            "COP projection: eta = +/-90",
        ));
    }
    let cot_theta_a = libm::cos(theta_a) / libm::sin(theta_a);
    let y0 = cos_eta * cot_theta_a;
    let c = libm::sin(theta_a);

    let (phi, r_unsigned) = deproject_conic_polar(x, y, y0, theta_a, c);
    let r_signed = theta_a.signum() * r_unsigned;

    let theta = theta_a + libm::atan(cot_theta_a - r_signed / cos_eta);

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_coe(
    native: NativeCoord,
    theta_a_deg: f64,
    eta_deg: f64,
) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_a = theta_a_deg * DEG_TO_RAD;
    let eta = eta_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COE projection: theta_a")?;

    let (theta_1, theta_2) = standard_parallels(theta_a, eta);
    let sin_t1 = libm::sin(theta_1);
    let sin_t2 = libm::sin(theta_2);
    let gamma = sin_t1 + sin_t2;
    if gamma.abs() < 1e-15 {
        return Err(WcsError::invalid_parameter(
            "COE projection: sin(theta_1) + sin(theta_2) = 0",
        ));
    }
    let c = gamma / 2.0;

    let bracket_theta = 1.0 + sin_t1 * sin_t2 - gamma * libm::sin(theta);
    if bracket_theta < 0.0 {
        return Err(WcsError::out_of_bounds(
            "COE projection: point outside valid region",
        ));
    }
    let r_theta = (2.0 / gamma) * libm::sqrt(bracket_theta);

    let bracket_a = 1.0 + sin_t1 * sin_t2 - gamma * libm::sin(theta_a);
    let y0 = (2.0 / gamma) * libm::sqrt(bracket_a.max(0.0));

    Ok(project_conic_xy(r_theta, y0, c, phi))
}

pub(crate) fn deproject_coe(
    inter: IntermediateCoord,
    theta_a_deg: f64,
    eta_deg: f64,
) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_a = theta_a_deg * DEG_TO_RAD;
    let eta = eta_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COE projection: theta_a")?;

    let (theta_1, theta_2) = standard_parallels(theta_a, eta);
    let sin_t1 = libm::sin(theta_1);
    let sin_t2 = libm::sin(theta_2);
    let gamma = sin_t1 + sin_t2;
    if gamma.abs() < 1e-15 {
        return Err(WcsError::invalid_parameter(
            "COE projection: sin(theta_1) + sin(theta_2) = 0",
        ));
    }

    let bracket_a = 1.0 + sin_t1 * sin_t2 - gamma * libm::sin(theta_a);
    let y0 = (2.0 / gamma) * libm::sqrt(bracket_a.max(0.0));
    let c = gamma / 2.0;

    let (phi, r_unsigned) = deproject_conic_polar(x, y, y0, theta_a, c);
    let r_signed = theta_a.signum() * r_unsigned;

    let sin_theta = (1.0 + sin_t1 * sin_t2) / gamma - gamma * (r_signed / 2.0).powi(2);
    if sin_theta.abs() > 1.0 + 1e-12 {
        return Err(WcsError::out_of_bounds(
            "COE deprojection: point outside valid region",
        ));
    }
    let theta = libm::asin(sin_theta.clamp(-1.0, 1.0));

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_cod(
    native: NativeCoord,
    theta_a_deg: f64,
    eta_deg: f64,
) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_a = theta_a_deg * DEG_TO_RAD;
    let eta = eta_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COD projection: theta_a")?;

    let sin_t_a = libm::sin(theta_a);
    let cos_t_a = libm::cos(theta_a);
    let cot_theta_a = cos_t_a / sin_t_a;

    let c = sin_t_a * sin_eta_over_eta(eta);
    let y0 = eta_cot_eta(eta) * cot_theta_a;
    let r_theta = theta_a - theta + y0;

    Ok(project_conic_xy(r_theta, y0, c, phi))
}

pub(crate) fn deproject_cod(
    inter: IntermediateCoord,
    theta_a_deg: f64,
    eta_deg: f64,
) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_a = theta_a_deg * DEG_TO_RAD;
    let eta = eta_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COD projection: theta_a")?;

    let sin_t_a = libm::sin(theta_a);
    let cos_t_a = libm::cos(theta_a);
    let cot_theta_a = cos_t_a / sin_t_a;
    let y0 = eta_cot_eta(eta) * cot_theta_a;
    let c = sin_t_a * sin_eta_over_eta(eta);

    let (phi, r_unsigned) = deproject_conic_polar(x, y, y0, theta_a, c);
    let r_signed = theta_a.signum() * r_unsigned;

    let theta = theta_a + y0 - r_signed;

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_coo(
    native: NativeCoord,
    theta_a_deg: f64,
    eta_deg: f64,
) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();
    let theta_a = theta_a_deg * DEG_TO_RAD;
    let eta = eta_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COO projection: theta_a")?;

    if theta.abs() >= HALF_PI - 1e-10 && theta.signum() != theta_a.signum() {
        return Err(WcsError::singularity(
            "COO projection: singularity at opposite pole",
        ));
    }

    let (c, psi) = coo_c_psi(theta_a, eta)?;
    let tan_half_xi = libm::tan((HALF_PI - theta) / 2.0);
    if tan_half_xi <= 0.0 {
        return Err(WcsError::singularity(
            "COO projection: tan((90-theta)/2) <= 0",
        ));
    }
    let r_theta = psi * tan_half_xi.powf(c);
    let tan_half_xi_a = libm::tan((HALF_PI - theta_a) / 2.0);
    let y0 = psi * tan_half_xi_a.powf(c);

    Ok(project_conic_xy(r_theta, y0, c, phi))
}

pub(crate) fn deproject_coo(
    inter: IntermediateCoord,
    theta_a_deg: f64,
    eta_deg: f64,
) -> WcsResult<NativeCoord> {
    let x = inter.x_deg() * DEG_TO_RAD;
    let y = inter.y_deg() * DEG_TO_RAD;
    let theta_a = theta_a_deg * DEG_TO_RAD;
    let eta = eta_deg * DEG_TO_RAD;

    check_nonzero_param(theta_a, "COO projection: theta_a")?;

    let (c, psi) = coo_c_psi(theta_a, eta)?;
    let tan_half_xi_a = libm::tan((HALF_PI - theta_a) / 2.0);
    let y0 = psi * tan_half_xi_a.powf(c);

    let (phi, r_unsigned) = deproject_conic_polar(x, y, y0, theta_a, c);

    if r_unsigned < 1e-15 {
        return Ok(NativeCoord::new(
            Angle::from_degrees(0.0),
            Angle::from_degrees(90.0 * theta_a.signum()),
        ));
    }

    let r_signed = theta_a.signum() * r_unsigned;
    let tan_half_xi = (r_signed / psi).powf(1.0 / c);
    let theta = HALF_PI - 2.0 * libm::atan(tan_half_xi);

    Ok(native_coord_from_radians(phi, theta))
}

fn coo_c_psi(theta_a: f64, eta: f64) -> WcsResult<(f64, f64)> {
    let (theta_1, theta_2) = standard_parallels(theta_a, eta);
    let cos_t1 = libm::cos(theta_1);
    let cos_t2 = libm::cos(theta_2);
    if cos_t1 <= 0.0 || cos_t2 <= 0.0 {
        return Err(WcsError::invalid_parameter(
            "COO projection: standard parallels reach the poles",
        ));
    }
    let tan_half_xi_1 = libm::tan((HALF_PI - theta_1) / 2.0);
    let tan_half_xi_2 = libm::tan((HALF_PI - theta_2) / 2.0);
    if tan_half_xi_1 <= 0.0 || tan_half_xi_2 <= 0.0 {
        return Err(WcsError::singularity(
            "COO projection: tan((90-theta_i)/2) <= 0",
        ));
    }

    let c = if eta.abs() < 1e-12 {
        libm::sin(theta_a)
    } else {
        libm::log(cos_t2 / cos_t1) / libm::log(tan_half_xi_2 / tan_half_xi_1)
    };
    if c.abs() < 1e-15 {
        return Err(WcsError::invalid_parameter(
            "COO projection: degenerate cone (C = 0)",
        ));
    }

    let psi = cos_t1 / (c * tan_half_xi_1.powf(c));

    Ok((c, psi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;

    fn assert_reaches_origin(proj: &Projection, theta_a: f64) {
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(theta_a));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10, "x not zero for {:?}", proj);
        assert!(inter.y_deg().abs() < 1e-10, "y not zero for {:?}", proj);
    }

    fn assert_roundtrip(proj: &Projection, phi_deg: f64, theta_deg: f64, ulp: u64) {
        let original = NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), ulp,
            "phi mismatch at ({}, {})", phi_deg, theta_deg);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), ulp,
            "theta mismatch at ({}, {})", phi_deg, theta_deg);
    }

    #[test]
    fn test_conic_native_reference_and_origin() {
        let theta_a = 45.0;
        let projs = [
            Projection::cop(theta_a, 0.0),
            Projection::coe(theta_a, 0.0),
            Projection::cod(theta_a, 0.0),
            Projection::coo(theta_a, 0.0),
            Projection::cop(theta_a, 15.0),
            Projection::coe(theta_a, 15.0),
            Projection::cod(theta_a, 15.0),
            Projection::coo(theta_a, 15.0),
        ];
        for proj in &projs {
            let (phi0, theta0) = proj.native_reference();
            assert_eq!(phi0, 0.0);
            assert_eq!(theta0, theta_a);
            assert_reaches_origin(proj, theta_a);
        }
    }

    #[test]
    fn test_cop_roundtrip_one_standard() {
        for theta_a in [30.0, 45.0, 60.0, 75.0] {
            let proj = Projection::cop(theta_a, 0.0);
            for phi_deg in [-120.0, -60.0, 0.0, 60.0, 120.0] {
                for theta_deg in [20.0, 40.0, 60.0, 80.0] {
                    assert_roundtrip(&proj, phi_deg, theta_deg, 15);
                }
            }
        }
    }

    #[test]
    fn test_cop_roundtrip_two_standard() {
        // Paper II Fig. 24: theta_1 = 20, theta_2 = 70  =>  theta_a = 45, eta = 25
        let proj = Projection::cop(45.0, 25.0);
        for phi_deg in [-90.0, -30.0, 0.0, 30.0, 90.0] {
            for theta_deg in [25.0, 45.0, 65.0, 85.0] {
                assert_roundtrip(&proj, phi_deg, theta_deg, 20);
            }
        }
    }

    #[test]
    fn test_coe_roundtrip_one_and_two_standard() {
        for eta in [0.0, 10.0, 25.0] {
            let proj = Projection::coe(45.0, eta);
            for phi_deg in [-120.0, 0.0, 90.0] {
                for theta_deg in [20.0, 50.0, 80.0] {
                    assert_roundtrip(&proj, phi_deg, theta_deg, 20);
                }
            }
        }
    }

    #[test]
    fn test_cod_roundtrip_one_and_two_standard() {
        for eta in [0.0, 5.0, 20.0] {
            let proj = Projection::cod(45.0, eta);
            for phi_deg in [-120.0, 0.0, 90.0] {
                for theta_deg in [20.0, 50.0, 80.0] {
                    assert_roundtrip(&proj, phi_deg, theta_deg, 20);
                }
            }
        }
    }

    #[test]
    fn test_coo_roundtrip_one_and_two_standard() {
        for eta in [0.0, 10.0, 25.0] {
            let proj = Projection::coo(45.0, eta);
            for phi_deg in [-90.0, 0.0, 90.0] {
                for theta_deg in [25.0, 50.0, 80.0] {
                    assert_roundtrip(&proj, phi_deg, theta_deg, 30);
                }
            }
        }
    }

    #[test]
    fn test_conic_error_cases() {
        // theta_a = 0 is invalid for every conic (the cone degenerates to a plane).
        for proj in [
            Projection::cop(0.0, 0.0),
            Projection::coe(0.0, 0.0),
            Projection::cod(0.0, 0.0),
            Projection::coo(0.0, 0.0),
        ] {
            let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
            assert!(proj.project(native).is_err(), "{:?} theta_a=0 should error", proj);
        }

        // COP: theta - theta_a = +/-90 is singular (cos(delta) = 0).  With
        // theta_a = 45, theta = -45 hits the singularity.
        let cop = Projection::cop(45.0, 0.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-45.0));
        assert!(cop.project(native).is_err());

        // COO is singular at the pole opposite theta_a (tan((90-theta)/2) <= 0).
        let coo = Projection::coo(45.0, 0.0);
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-90.0));
        assert!(coo.project(native).is_err());
    }

    #[test]
    fn test_conic_eta_zero_matches_eta_small() {
        // For each conic, eta=0 should agree numerically with eta=1e-8 within
        // a generous tolerance.  This guards against typos in the eta=0 limit
        // branches of c, y0 and the inverse formulae.
        let theta_a = 45.0;
        let pt = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(55.0));
        let pairs = [
            (Projection::cop(theta_a, 0.0), Projection::cop(theta_a, 1e-8)),
            (Projection::coe(theta_a, 0.0), Projection::coe(theta_a, 1e-8)),
            (Projection::cod(theta_a, 0.0), Projection::cod(theta_a, 1e-8)),
            (Projection::coo(theta_a, 0.0), Projection::coo(theta_a, 1e-8)),
        ];
        for (zero, small) in pairs {
            let a = zero.project(pt).unwrap();
            let b = small.project(pt).unwrap();
            assert!((a.x_deg() - b.x_deg()).abs() < 1e-6,
                "{:?}: x diverges at eta -> 0: {} vs {}", zero, a.x_deg(), b.x_deg());
            assert!((a.y_deg() - b.y_deg()).abs() < 1e-6,
                "{:?}: y diverges at eta -> 0: {} vs {}", zero, a.y_deg(), b.y_deg());
        }
    }
}
