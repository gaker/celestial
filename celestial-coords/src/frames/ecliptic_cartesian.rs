use crate::transforms::cartesian::CartesianFrame;
use celestial_core::ecliptic::{icrs_to_vsop2013, vsop2013_to_icrs};
use celestial_core::matrix::Vector3;

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
        vsop2013_to_icrs(&Vector3::new(self.x, self.y, self.z))
    }

    fn from_icrs(icrs: &Vector3) -> Self {
        let v = icrs_to_vsop2013(icrs);
        Self::new(v.x, v.y, v.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::{VSOP2013_OBLIQUITY_RAD, VSOP2013_PHI_RAD};

    #[test]
    fn test_roundtrip() {
        let ecl = EclipticCartesian::new(-9.8753625435, -27.9588613710, 5.8504463318);
        assert_eq!(EclipticCartesian::from_icrs(&ecl.to_icrs()), ecl);
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
    fn test_axes_follow_vsop2013_equation_3() {
        // VSOP2013 eq. (3): ecliptic = R1(eps) R3(phi) ICRS. These are the rows of that matrix.
        let (se, ce) = libm::sincos(VSOP2013_OBLIQUITY_RAD);
        let (sp, cp) = libm::sincos(VSOP2013_PHI_RAD);
        let rows = [
            [cp, sp, 0.0],
            [-sp * ce, cp * ce, se],
            [sp * se, -cp * se, ce],
        ];
        let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        for (i, [x, y, z]) in axes.into_iter().enumerate() {
            let row = rows[i];
            let ecl = EclipticCartesian::new(x, y, z);
            assert_eq!(ecl.to_icrs(), Vector3::new(row[0], row[1], row[2]), "{i}");

            let column = EclipticCartesian::new(rows[0][i], rows[1][i], rows[2][i]);
            let icrs = Vector3::new(x, y, z);
            assert_eq!(EclipticCartesian::from_icrs(&icrs), column, "{i}");
        }
    }
}
