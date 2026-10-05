// The ecliptic frame of VSOP2013: the inertial mean ecliptic and equinox of
// J2000, related to ICRS by eq. 3 of the VSOP2013 paper,
// ecliptic = R1(ε) R3(φ) ICRS.

use crate::constants::{VSOP2013_OBLIQUITY_RAD, VSOP2013_PHI_RAD};
use crate::matrix::Vector3;

pub fn vsop2013_to_icrs(v: &Vector3) -> Vector3 {
    let (sin_eps, cos_eps) = libm::sincos(VSOP2013_OBLIQUITY_RAD);
    let (sin_phi, cos_phi) = libm::sincos(VSOP2013_PHI_RAD);
    let y1 = v.y * cos_eps - v.z * sin_eps;
    let z1 = v.y * sin_eps + v.z * cos_eps;
    Vector3::new(
        v.x * cos_phi - y1 * sin_phi,
        v.x * sin_phi + y1 * cos_phi,
        z1,
    )
}

pub fn icrs_to_vsop2013(v: &Vector3) -> Vector3 {
    let (sin_eps, cos_eps) = libm::sincos(VSOP2013_OBLIQUITY_RAD);
    let (sin_phi, cos_phi) = libm::sincos(VSOP2013_PHI_RAD);
    let x1 = v.x * cos_phi + v.y * sin_phi;
    let y1 = v.y * cos_phi - v.x * sin_phi;
    Vector3::new(
        x1,
        y1 * cos_eps + v.z * sin_eps,
        -y1 * sin_eps + v.z * cos_eps,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // The rows of R1(ε) R3(φ).
    fn rows() -> [[f64; 3]; 3] {
        let (se, ce) = libm::sincos(VSOP2013_OBLIQUITY_RAD);
        let (sp, cp) = libm::sincos(VSOP2013_PHI_RAD);
        [
            [cp, sp, 0.0],
            [-sp * ce, cp * ce, se],
            [sp * se, -cp * se, ce],
        ]
    }

    #[test]
    fn axes_follow_vsop2013_equation_3() {
        let rows = rows();
        let axes = [Vector3::x_axis(), Vector3::y_axis(), Vector3::z_axis()];
        for (i, axis) in axes.iter().enumerate() {
            let row = Vector3::from_array(rows[i]);
            assert_eq!(vsop2013_to_icrs(axis), row, "{i}");
            let column = Vector3::new(rows[0][i], rows[1][i], rows[2][i]);
            assert_eq!(icrs_to_vsop2013(axis), column, "{i}");
        }
    }

    #[test]
    fn round_trip() {
        let v = Vector3::new(-9.8753625435, -27.958861371, 5.8504463318);
        assert_eq!(icrs_to_vsop2013(&vsop2013_to_icrs(&v)), v);
    }
}
