use crate::transforms::CartesianFrame;
use celestial_core::constants::{FRAME_BIAS_PHI_RAD, J2000_OBLIQUITY_RAD};
use celestial_core::Vector3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EclipticCartesian {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl EclipticCartesian {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn from_vector3(v: &Vector3) -> Self {
        Self {
            x: v.x,
            y: v.y,
            z: v.z,
        }
    }
}

impl CartesianFrame for EclipticCartesian {
    fn to_icrs(&self) -> Vector3 {
        let eps = J2000_OBLIQUITY_RAD;
        let phi = FRAME_BIAS_PHI_RAD;
        let (sin_eps, cos_eps) = eps.sin_cos();
        let (sin_phi, cos_phi) = phi.sin_cos();

        let y1 = self.y * cos_eps - self.z * sin_eps;
        let z1 = self.y * sin_eps + self.z * cos_eps;

        Vector3::new(
            self.x * cos_phi + y1 * sin_phi,
            -self.x * sin_phi + y1 * cos_phi,
            z1,
        )
    }

    fn from_icrs(icrs: &Vector3) -> Self {
        let eps = J2000_OBLIQUITY_RAD;
        let phi = FRAME_BIAS_PHI_RAD;
        let (sin_eps, cos_eps) = eps.sin_cos();
        let (sin_phi, cos_phi) = phi.sin_cos();

        let x1 = icrs.x * cos_phi - icrs.y * sin_phi;
        let y1 = icrs.x * sin_phi + icrs.y * cos_phi;

        Self {
            x: x1,
            y: y1 * cos_eps + icrs.z * sin_eps,
            z: -y1 * sin_eps + icrs.z * cos_eps,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip() {
        let ecl = EclipticCartesian::new(-9.8753625435, -27.9588613710, 5.8504463318);
        let icrs = ecl.to_icrs();
        let back = EclipticCartesian::from_icrs(&icrs);

        let tol = 1e-14;
        assert!((ecl.x - back.x).abs() < tol, "X roundtrip error");
        assert!((ecl.y - back.y).abs() < tol, "Y roundtrip error");
        assert!((ecl.z - back.z).abs() < tol, "Z roundtrip error");
    }

    #[test]
    fn test_pluto_vector_signs() {
        let ecl = EclipticCartesian::new(-9.8753625435, -27.9588613710, 5.8504463318);
        let icrs = ecl.to_icrs();

        assert!(icrs.x < 0.0, "X should be negative");
        assert!(icrs.y < 0.0, "Y should be negative");
        assert!(icrs.z < 0.0, "Z should be negative in ICRS");
    }

    #[test]
    fn test_ecliptic_x_axis() {
        let ecl = EclipticCartesian::new(1.0, 0.0, 0.0);
        let icrs = ecl.to_icrs();

        assert!(
            (icrs.x - 1.0).abs() < 1e-10,
            "X-axis X component should be ~1"
        );
        assert!(icrs.y.abs() < 1e-6, "X-axis Y component should be ~0");
        assert!(icrs.z.abs() < 1e-10, "X-axis Z component should be ~0");
    }
}
