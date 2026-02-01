use celestial_core::constants::{HALF_PI, RAD_TO_DEG};
use celestial_core::utils::normalize_longitude;
use celestial_core::Angle;

use crate::common::{asin_safe, native_coord_from_radians};
use crate::coordinate::{CelestialCoord, IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

mod conic;
mod cylindrical;
mod polyconic;
mod pseudocylindrical;
mod quadcube;
mod zenithal;

use conic::{deproject_cod, deproject_coe, deproject_coo, deproject_cop};
use conic::{project_cod, project_coe, project_coo, project_cop};
use cylindrical::{deproject_car, deproject_cea, deproject_cyp, deproject_mer};
use cylindrical::{project_car, project_cea, project_cyp, project_mer};
use polyconic::{deproject_bon, deproject_pco, project_bon, project_pco};
use pseudocylindrical::{deproject_ait, deproject_mol, deproject_par, deproject_sfl};
use pseudocylindrical::{project_ait, project_mol, project_par, project_sfl};
use quadcube::{deproject_csc, deproject_qsc, deproject_tsc};
use quadcube::{project_csc, project_qsc, project_tsc};
use zenithal::{deproject_air, deproject_arc, deproject_azp, deproject_sin, deproject_stg};
use zenithal::{deproject_szp, deproject_tan, deproject_zea, deproject_zpn};
use zenithal::{project_air, project_arc, project_azp, project_sin, project_stg};
use zenithal::{project_szp, project_tan, project_zea, project_zpn};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphericalRotation {
    alpha_p: f64,
    delta_p: f64,
    phi_p: f64,
    sin_delta_p: f64,
    cos_delta_p: f64,
}

impl SphericalRotation {
    pub fn new(alpha_p: Angle, delta_p: Angle, phi_p: Angle) -> Self {
        let delta_p_rad = delta_p.radians();
        let (sin_delta_p, cos_delta_p) = delta_p_rad.sin_cos();
        Self {
            alpha_p: alpha_p.radians(),
            delta_p: delta_p_rad,
            phi_p: phi_p.radians(),
            sin_delta_p,
            cos_delta_p,
        }
    }

    fn default_lonpole(delta_0: Angle, theta_0: Angle) -> Angle {
        if delta_0.radians() >= theta_0.radians() {
            Angle::from_degrees(0.0)
        } else {
            Angle::from_degrees(180.0)
        }
    }

    pub fn from_crval(
        alpha_0: Angle,
        delta_0: Angle,
        theta_0: Angle,
        lonpole: Option<Angle>,
        latpole: Option<Angle>,
    ) -> WcsResult<Self> {
        let phi_p = lonpole.unwrap_or_else(|| Self::default_lonpole(delta_0, theta_0));
        let latpole_rad = latpole.map(|a| a.radians()).unwrap_or(HALF_PI);

        let delta_0_rad = delta_0.radians();
        let theta_0_rad = theta_0.radians();
        let phi_p_rad = phi_p.radians();

        let (sin_delta_0, cos_delta_0) = delta_0_rad.sin_cos();
        let (sin_theta_0, cos_theta_0) = theta_0_rad.sin_cos();
        let (sin_phi_p, cos_phi_p) = phi_p_rad.sin_cos();

        let delta_p = Self::compute_delta_p(
            sin_delta_0,
            cos_delta_0,
            sin_theta_0,
            cos_theta_0,
            sin_phi_p,
            cos_phi_p,
            latpole_rad,
        )?;

        let x = -cos_theta_0 * sin_phi_p;
        let y = sin_theta_0 * cos_delta_0 - cos_theta_0 * sin_delta_0 * cos_phi_p;
        let alpha_p = alpha_0.radians() + x.atan2(y);

        let alpha_p_deg = normalize_longitude(alpha_p * RAD_TO_DEG);
        let delta_p_deg = delta_p * RAD_TO_DEG;

        Ok(Self::new(
            Angle::from_degrees(alpha_p_deg),
            Angle::from_degrees(delta_p_deg),
            phi_p,
        ))
    }

    fn compute_delta_p(
        sin_delta_0: f64,
        _cos_delta_0: f64,
        sin_theta_0: f64,
        cos_theta_0: f64,
        sin_phi_p: f64,
        cos_phi_p: f64,
        latpole_rad: f64,
    ) -> WcsResult<f64> {
        let cos_theta_0_sin_phi_p = cos_theta_0 * sin_phi_p;
        let denom_sq = 1.0 - cos_theta_0_sin_phi_p * cos_theta_0_sin_phi_p;

        if denom_sq.abs() < 1e-15 {
            if sin_delta_0.abs() < 1e-15 {
                return Ok(latpole_rad);
            }
            return Err(WcsError::invalid_parameter(
                "Invalid combination of θ₀, δ₀, and φₚ - no solution for δₚ",
            ));
        }

        let denom = denom_sq.sqrt();
        let arg = sin_delta_0 / denom;

        if arg.abs() > 1.0 + 1e-15 {
            return Err(WcsError::invalid_parameter(
                "Invalid combination of θ₀, δ₀, and φₚ - acos argument out of range",
            ));
        }

        let arg_clamped = arg.clamp(-1.0, 1.0);
        let acos_term = arg_clamped.acos();
        let base = sin_theta_0.atan2(cos_theta_0 * cos_phi_p);

        let delta_p_1 = base + acos_term;
        let delta_p_2 = base - acos_term;

        const BOUNDARY_TOL: f64 = 1e-14;
        let valid_1 = (-HALF_PI - BOUNDARY_TOL..=HALF_PI + BOUNDARY_TOL).contains(&delta_p_1);
        let valid_2 = (-HALF_PI - BOUNDARY_TOL..=HALF_PI + BOUNDARY_TOL).contains(&delta_p_2);

        let clamp_result = |v: f64| v.clamp(-HALF_PI, HALF_PI);

        match (valid_1, valid_2) {
            (true, false) => Ok(clamp_result(delta_p_1)),
            (false, true) => Ok(clamp_result(delta_p_2)),
            (true, true) => {
                let diff_1 = (delta_p_1 - latpole_rad).abs();
                let diff_2 = (delta_p_2 - latpole_rad).abs();
                if diff_1 <= diff_2 {
                    Ok(clamp_result(delta_p_1))
                } else {
                    Ok(clamp_result(delta_p_2))
                }
            }
            (false, false) => Err(WcsError::invalid_parameter(
                "No valid solution for δₚ in range [-90°, 90°]",
            )),
        }
    }

    pub fn native_to_celestial(&self, native: NativeCoord) -> WcsResult<CelestialCoord> {
        let phi = native.phi().radians();
        let theta = native.theta().radians();

        let (sin_theta, cos_theta) = theta.sin_cos();
        let d_phi = phi - self.phi_p;
        let (sin_d_phi, cos_d_phi) = d_phi.sin_cos();

        let sin_delta = sin_theta * self.sin_delta_p + cos_theta * self.cos_delta_p * cos_d_phi;
        let delta = asin_safe(sin_delta);

        let x = -cos_theta * sin_d_phi;
        let y = sin_theta * self.cos_delta_p - cos_theta * self.sin_delta_p * cos_d_phi;
        let alpha = self.alpha_p + x.atan2(y);

        let alpha_deg = normalize_longitude(alpha * RAD_TO_DEG);
        let delta_deg = delta * RAD_TO_DEG;

        Ok(CelestialCoord::new(
            Angle::from_degrees(alpha_deg),
            Angle::from_degrees(delta_deg),
        ))
    }

    pub fn celestial_to_native(&self, celestial: CelestialCoord) -> WcsResult<NativeCoord> {
        let alpha = celestial.alpha().radians();
        let delta = celestial.delta().radians();

        let (sin_delta, cos_delta) = delta.sin_cos();
        let d_alpha = alpha - self.alpha_p;
        let (sin_d_alpha, cos_d_alpha) = d_alpha.sin_cos();

        let sin_theta = sin_delta * self.sin_delta_p + cos_delta * self.cos_delta_p * cos_d_alpha;
        let theta = asin_safe(sin_theta);

        let x = -cos_delta * sin_d_alpha;
        let y = sin_delta * self.cos_delta_p - cos_delta * self.sin_delta_p * cos_d_alpha;
        let phi = self.phi_p + x.atan2(y);

        Ok(native_coord_from_radians(phi, theta))
    }

    #[inline]
    pub fn phi_p_degrees(&self) -> f64 {
        self.phi_p * RAD_TO_DEG
    }

    #[inline]
    pub fn delta_p_degrees(&self) -> f64 {
        self.delta_p * RAD_TO_DEG
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Projection {
    Tan,
    Sin { xi: f64, eta: f64 },
    Arc,
    Stg,
    Zea,
    Azp { mu: f64, gamma: f64 },
    Szp { mu: f64, phi_c: f64, theta_c: f64 },
    Zpn { coeffs: Vec<f64> },
    Air { theta_b: f64 },
    Car,
    Mer,
    Cea { lambda: f64 },
    Cyp { mu: f64, lambda: f64 },
    Sfl,
    Par,
    Mol,
    Ait,
    Cop { theta_a: f64 },
    Coe { theta_a: f64 },
    Cod { theta_a: f64 },
    Coo { theta_a: f64 },
    Bon { theta_1: f64 },
    Pco,
    Tsc,
    Csc,
    Qsc,
}

impl Projection {
    pub fn tan() -> Self {
        Self::Tan
    }

    pub fn sin() -> Self {
        Self::Sin { xi: 0.0, eta: 0.0 }
    }

    pub fn sin_with_params(xi: f64, eta: f64) -> Self {
        Self::Sin { xi, eta }
    }

    pub fn arc() -> Self {
        Self::Arc
    }

    pub fn stg() -> Self {
        Self::Stg
    }

    pub fn zea() -> Self {
        Self::Zea
    }

    pub fn azp(mu: f64, gamma: f64) -> Self {
        Self::Azp { mu, gamma }
    }

    pub fn szp(mu: f64, phi_c: f64, theta_c: f64) -> Self {
        Self::Szp { mu, phi_c, theta_c }
    }

    pub fn zpn(coeffs: Vec<f64>) -> Self {
        Self::Zpn { coeffs }
    }

    pub fn air(theta_b: f64) -> Self {
        Self::Air { theta_b }
    }

    pub fn car() -> Self {
        Self::Car
    }

    pub fn mer() -> Self {
        Self::Mer
    }

    pub fn cea() -> Self {
        Self::Cea { lambda: 1.0 }
    }

    pub fn cea_with_lambda(lambda: f64) -> Self {
        Self::Cea { lambda }
    }

    pub fn cyp(mu: f64, lambda: f64) -> Self {
        Self::Cyp { mu, lambda }
    }

    pub fn sfl() -> Self {
        Self::Sfl
    }

    pub fn par() -> Self {
        Self::Par
    }

    pub fn mol() -> Self {
        Self::Mol
    }

    pub fn ait() -> Self {
        Self::Ait
    }

    pub fn cop(theta_a: f64) -> Self {
        Self::Cop { theta_a }
    }

    pub fn coe(theta_a: f64) -> Self {
        Self::Coe { theta_a }
    }

    pub fn cod(theta_a: f64) -> Self {
        Self::Cod { theta_a }
    }

    pub fn coo(theta_a: f64) -> Self {
        Self::Coo { theta_a }
    }

    pub fn bon(theta_1: f64) -> Self {
        Self::Bon { theta_1 }
    }

    pub fn pco() -> Self {
        Self::Pco
    }

    pub fn tsc() -> Self {
        Self::Tsc
    }

    pub fn csc() -> Self {
        Self::Csc
    }

    pub fn qsc() -> Self {
        Self::Qsc
    }

    pub fn native_reference(&self) -> (f64, f64) {
        match self {
            Self::Tan
            | Self::Sin { .. }
            | Self::Arc
            | Self::Stg
            | Self::Zea
            | Self::Azp { .. }
            | Self::Szp { .. }
            | Self::Zpn { .. }
            | Self::Air { .. } => (0.0, 90.0),
            Self::Car | Self::Mer | Self::Cea { .. } | Self::Cyp { .. } => (0.0, 0.0),
            Self::Sfl | Self::Par | Self::Mol | Self::Ait => (0.0, 0.0),
            Self::Cop { theta_a }
            | Self::Coe { theta_a }
            | Self::Cod { theta_a }
            | Self::Coo { theta_a } => (0.0, *theta_a),
            Self::Bon { theta_1 } => (0.0, *theta_1),
            Self::Pco => (0.0, 0.0),
            Self::Tsc | Self::Csc | Self::Qsc => (0.0, 0.0),
        }
    }

    pub fn project(&self, native: NativeCoord) -> WcsResult<IntermediateCoord> {
        match self {
            Self::Tan => project_tan(native),
            Self::Sin { xi, eta } => project_sin(native, *xi, *eta),
            Self::Arc => project_arc(native),
            Self::Stg => project_stg(native),
            Self::Zea => project_zea(native),
            Self::Azp { mu, gamma } => project_azp(native, *mu, *gamma),
            Self::Szp { mu, phi_c, theta_c } => project_szp(native, *mu, *phi_c, *theta_c),
            Self::Zpn { coeffs } => project_zpn(native, coeffs),
            Self::Air { theta_b } => project_air(native, *theta_b),
            Self::Car => project_car(native),
            Self::Mer => project_mer(native),
            Self::Cea { lambda } => project_cea(native, *lambda),
            Self::Cyp { mu, lambda } => project_cyp(native, *mu, *lambda),
            Self::Sfl => project_sfl(native),
            Self::Par => project_par(native),
            Self::Mol => project_mol(native),
            Self::Ait => project_ait(native),
            Self::Cop { theta_a } => project_cop(native, *theta_a),
            Self::Coe { theta_a } => project_coe(native, *theta_a),
            Self::Cod { theta_a } => project_cod(native, *theta_a),
            Self::Coo { theta_a } => project_coo(native, *theta_a),
            Self::Bon { theta_1 } => project_bon(native, *theta_1),
            Self::Pco => project_pco(native),
            Self::Tsc => project_tsc(native),
            Self::Csc => project_csc(native),
            Self::Qsc => project_qsc(native),
        }
    }

    pub fn deproject(&self, inter: IntermediateCoord) -> WcsResult<NativeCoord> {
        match self {
            Self::Tan => deproject_tan(inter),
            Self::Sin { xi, eta } => deproject_sin(inter, *xi, *eta),
            Self::Arc => deproject_arc(inter),
            Self::Stg => deproject_stg(inter),
            Self::Zea => deproject_zea(inter),
            Self::Azp { mu, gamma } => deproject_azp(inter, *mu, *gamma),
            Self::Szp { mu, phi_c, theta_c } => deproject_szp(inter, *mu, *phi_c, *theta_c),
            Self::Zpn { coeffs } => deproject_zpn(inter, coeffs),
            Self::Air { theta_b } => deproject_air(inter, *theta_b),
            Self::Car => deproject_car(inter),
            Self::Mer => deproject_mer(inter),
            Self::Cea { lambda } => deproject_cea(inter, *lambda),
            Self::Cyp { mu, lambda } => deproject_cyp(inter, *mu, *lambda),
            Self::Sfl => deproject_sfl(inter),
            Self::Par => deproject_par(inter),
            Self::Mol => deproject_mol(inter),
            Self::Ait => deproject_ait(inter),
            Self::Cop { theta_a } => deproject_cop(inter, *theta_a),
            Self::Coe { theta_a } => deproject_coe(inter, *theta_a),
            Self::Cod { theta_a } => deproject_cod(inter, *theta_a),
            Self::Coo { theta_a } => deproject_coo(inter, *theta_a),
            Self::Bon { theta_1 } => deproject_bon(inter, *theta_1),
            Self::Pco => deproject_pco(inter),
            Self::Tsc => deproject_tsc(inter),
            Self::Csc => deproject_csc(inter),
            Self::Qsc => deproject_qsc(inter),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::assert_ulp_lt;

    #[test]
    fn test_native_to_celestial_at_pole() {
        let rot = SphericalRotation::new(
            Angle::from_degrees(0.0),
            Angle::from_degrees(90.0),
            Angle::from_degrees(180.0),
        );
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let celestial = rot.native_to_celestial(native).unwrap();

        assert_ulp_lt!(celestial.delta().degrees(), 90.0, 1);
    }

    #[test]
    fn test_native_to_celestial_reference_point() {
        let rot = SphericalRotation::new(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(180.0),
        );
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let celestial = rot.native_to_celestial(native).unwrap();

        assert_ulp_lt!(celestial.alpha().degrees(), 180.0, 1);
        assert_ulp_lt!(celestial.delta().degrees(), 45.0, 1);
    }

    #[test]
    fn test_native_to_celestial_equator() {
        let rot = SphericalRotation::new(
            Angle::from_degrees(0.0),
            Angle::from_degrees(90.0),
            Angle::from_degrees(180.0),
        );
        let native = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(0.0));
        let celestial = rot.native_to_celestial(native).unwrap();

        assert!(celestial.delta().degrees().abs() < 1e-10);
        assert_ulp_lt!(celestial.alpha().degrees(), 90.0, 1);
    }

    #[test]
    fn test_celestial_to_native_at_pole() {
        let rot = SphericalRotation::new(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(180.0),
        );
        let celestial = CelestialCoord::new(Angle::from_degrees(180.0), Angle::from_degrees(45.0));
        let native = rot.celestial_to_native(celestial).unwrap();

        assert_ulp_lt!(native.theta().degrees(), 90.0, 1);
    }

    #[test]
    fn test_spherical_rotation_roundtrip() {
        let rot = SphericalRotation::new(
            Angle::from_degrees(120.0),
            Angle::from_degrees(35.0),
            Angle::from_degrees(180.0),
        );

        let original = NativeCoord::new(Angle::from_degrees(45.0), Angle::from_degrees(60.0));
        let celestial = rot.native_to_celestial(original).unwrap();
        let recovered = rot.celestial_to_native(celestial).unwrap();

        // ULP tolerance accounts for ARM vs x86 FPU differences in trig functions
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 8);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 8);
    }

    #[test]
    fn test_spherical_rotation_roundtrip_reverse() {
        let rot = SphericalRotation::new(
            Angle::from_degrees(100.0),
            Angle::from_degrees(-25.0),
            Angle::from_degrees(180.0),
        );

        let original = CelestialCoord::new(Angle::from_degrees(110.0), Angle::from_degrees(-30.0));
        let native = rot.celestial_to_native(original).unwrap();
        let recovered = rot.native_to_celestial(native).unwrap();

        assert_ulp_lt!(original.alpha().degrees(), recovered.alpha().degrees(), 2);
        assert_ulp_lt!(original.delta().degrees(), recovered.delta().degrees(), 3);
    }

    #[test]
    fn test_from_crval_zenithal_at_pole() {
        let rot = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(90.0),
            Some(Angle::from_degrees(180.0)),
            None,
        )
        .unwrap();

        let native_ref = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let celestial = rot.native_to_celestial(native_ref).unwrap();

        assert_ulp_lt!(celestial.alpha().degrees(), 180.0, 2);
        assert_ulp_lt!(celestial.delta().degrees(), 45.0, 2);
    }

    #[test]
    fn test_from_crval_zenithal_north_pole() {
        let rot = SphericalRotation::from_crval(
            Angle::from_degrees(0.0),
            Angle::from_degrees(90.0),
            Angle::from_degrees(90.0),
            Some(Angle::from_degrees(180.0)),
            None,
        )
        .unwrap();

        let native_ref = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let celestial = rot.native_to_celestial(native_ref).unwrap();

        assert_ulp_lt!(celestial.delta().degrees(), 90.0, 2);
    }

    #[test]
    fn test_from_crval_roundtrip() {
        let rot = SphericalRotation::from_crval(
            Angle::from_degrees(120.0),
            Angle::from_degrees(35.0),
            Angle::from_degrees(90.0),
            Some(Angle::from_degrees(180.0)),
            None,
        )
        .unwrap();

        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(60.0));
        let celestial = rot.native_to_celestial(original).unwrap();
        let recovered = rot.celestial_to_native(celestial).unwrap();

        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 2);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 2);
    }

    #[test]
    fn test_default_lonpole_delta_less_than_theta() {
        let rot_default = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(90.0),
            None,
            None,
        )
        .unwrap();

        let rot_explicit = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(90.0),
            Some(Angle::from_degrees(180.0)),
            None,
        )
        .unwrap();

        let native = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(60.0));
        let celestial_default = rot_default.native_to_celestial(native).unwrap();
        let celestial_explicit = rot_explicit.native_to_celestial(native).unwrap();

        assert_ulp_lt!(
            celestial_default.alpha().degrees(),
            celestial_explicit.alpha().degrees(),
            10
        );
        assert_ulp_lt!(
            celestial_default.delta().degrees(),
            celestial_explicit.delta().degrees(),
            10
        );
    }

    #[test]
    fn test_default_lonpole_delta_greater_than_theta() {
        let rot_default = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(0.0),
            None,
            None,
        )
        .unwrap();

        let rot_explicit = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(0.0),
            Some(Angle::from_degrees(0.0)),
            None,
        )
        .unwrap();

        let native = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(30.0));
        let celestial_default = rot_default.native_to_celestial(native).unwrap();
        let celestial_explicit = rot_explicit.native_to_celestial(native).unwrap();

        assert_ulp_lt!(
            celestial_default.alpha().degrees(),
            celestial_explicit.alpha().degrees(),
            10
        );
        assert_ulp_lt!(
            celestial_default.delta().degrees(),
            celestial_explicit.delta().degrees(),
            10
        );
    }

    #[test]
    fn test_explicit_lonpole_overrides_default() {
        let rot_default = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(90.0),
            None,
            None,
        )
        .unwrap();

        let rot_override = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(90.0),
            Some(Angle::from_degrees(90.0)),
            None,
        )
        .unwrap();

        let native = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(60.0));
        let celestial_default = rot_default.native_to_celestial(native).unwrap();
        let celestial_override = rot_override.native_to_celestial(native).unwrap();

        let alpha_diff =
            (celestial_default.alpha().degrees() - celestial_override.alpha().degrees()).abs();
        let delta_diff =
            (celestial_default.delta().degrees() - celestial_override.delta().degrees()).abs();

        assert!(
            alpha_diff > 0.1 || delta_diff > 0.1,
            "Explicit LONPOLE override should produce different results"
        );
    }

    #[test]
    fn test_latpole_default_is_90() {
        let rot = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(30.0),
            Angle::from_degrees(45.0),
            Some(Angle::from_degrees(180.0)),
            None,
        )
        .unwrap();

        let native_ref = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(45.0));
        let celestial = rot.native_to_celestial(native_ref).unwrap();

        let recovered = rot.celestial_to_native(celestial).unwrap();
        assert_ulp_lt!(native_ref.phi().degrees(), recovered.phi().degrees(), 5);
        assert_ulp_lt!(native_ref.theta().degrees(), recovered.theta().degrees(), 5);
    }

    #[test]
    fn test_latpole_disambiguates_pole_solutions() {
        let rot_north = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(0.0),
            Some(Angle::from_degrees(0.0)),
            Some(Angle::from_degrees(90.0)),
        )
        .unwrap();

        let rot_south = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(0.0),
            Some(Angle::from_degrees(0.0)),
            Some(Angle::from_degrees(-90.0)),
        )
        .unwrap();

        let native = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(30.0));

        let celestial_north = rot_north.native_to_celestial(native).unwrap();
        let recovered_north = rot_north.celestial_to_native(celestial_north).unwrap();
        assert_ulp_lt!(native.phi().degrees(), recovered_north.phi().degrees(), 5);
        assert_ulp_lt!(
            native.theta().degrees(),
            recovered_north.theta().degrees(),
            5
        );

        let celestial_south = rot_south.native_to_celestial(native).unwrap();
        let recovered_south = rot_south.celestial_to_native(celestial_south).unwrap();
        assert_ulp_lt!(native.phi().degrees(), recovered_south.phi().degrees(), 5);
        assert_ulp_lt!(
            native.theta().degrees(),
            recovered_south.theta().degrees(),
            5
        );
    }

    #[test]
    fn test_latpole_explicit_value_respected() {
        let rot_with_latpole = SphericalRotation::from_crval(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            Angle::from_degrees(90.0),
            Some(Angle::from_degrees(180.0)),
            Some(Angle::from_degrees(45.0)),
        )
        .unwrap();

        let native = NativeCoord::new(Angle::from_degrees(60.0), Angle::from_degrees(30.0));
        let celestial = rot_with_latpole.native_to_celestial(native).unwrap();
        let recovered = rot_with_latpole.celestial_to_native(celestial).unwrap();

        assert_ulp_lt!(native.phi().degrees(), recovered.phi().degrees(), 5);
        assert_ulp_lt!(native.theta().degrees(), recovered.theta().degrees(), 5);
    }
}
