use celestial_core::constants::{DEG_TO_RAD, HALF_PI, RAD_TO_DEG};

use crate::common::native_coord_from_radians;
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};

#[derive(Debug, Clone, Copy, PartialEq)]
struct QuadcubeFace {
    face: u8,
    xi: f64,
    eta: f64,
    zeta: f64,
    phi_c: f64,
    theta_c: f64,
}

fn select_quadcube_face(phi: f64, theta: f64) -> QuadcubeFace {
    let (sin_theta, cos_theta) = theta.sin_cos();
    let (sin_phi, cos_phi) = phi.sin_cos();

    let l = cos_theta * cos_phi;
    let m = cos_theta * sin_phi;
    let n = sin_theta;

    let candidates = [
        (0, m, -l, n, 0.0, HALF_PI),
        (1, m, n, l, 0.0, 0.0),
        (2, -l, n, m, HALF_PI, 0.0),
        (3, -m, n, -l, std::f64::consts::PI, 0.0),
        (4, l, n, -m, -HALF_PI, 0.0),
        (5, m, l, -n, 0.0, -HALF_PI),
    ];

    let mut best_face = 0u8;
    let mut best_zeta = f64::NEG_INFINITY;
    let mut best_xi = 0.0;
    let mut best_eta = 0.0;
    let mut best_phi_c = 0.0;
    let mut best_theta_c = 0.0;

    for (face, xi, eta, zeta, phi_c, theta_c) in candidates {
        if zeta > best_zeta {
            best_face = face;
            best_xi = xi;
            best_eta = eta;
            best_zeta = zeta;
            best_phi_c = phi_c;
            best_theta_c = theta_c;
        }
    }

    QuadcubeFace {
        face: best_face,
        xi: best_xi,
        eta: best_eta,
        zeta: best_zeta,
        phi_c: best_phi_c,
        theta_c: best_theta_c,
    }
}

fn quadcube_face_from_xy(x: f64, y: f64) -> (u8, f64, f64, f64, f64) {
    let x_deg = x;
    let y_deg = y;

    let face;
    let phi_c;
    let theta_c;

    if y_deg > 45.0 {
        face = 0;
        phi_c = 0.0;
        theta_c = 90.0;
    } else if y_deg < -45.0 {
        face = 5;
        phi_c = 0.0;
        theta_c = -90.0;
    } else if (-45.0..45.0).contains(&x_deg) {
        face = 1;
        phi_c = 0.0;
        theta_c = 0.0;
    } else if (45.0..135.0).contains(&x_deg) {
        face = 2;
        phi_c = 90.0;
        theta_c = 0.0;
    } else if !(-135.0..135.0).contains(&x_deg) {
        face = 3;
        phi_c = 180.0;
        theta_c = 0.0;
    } else {
        face = 4;
        phi_c = -90.0;
        theta_c = 0.0;
    }

    (face, phi_c, theta_c, x_deg - phi_c, y_deg - theta_c)
}

fn face_coords_to_direction_cosines(face: u8, xi: f64, eta: f64, zeta: f64) -> (f64, f64, f64) {
    match face {
        0 => (-eta, xi, zeta),
        1 => (zeta, xi, eta),
        2 => (-xi, zeta, eta),
        3 => (-zeta, -xi, eta),
        4 => (xi, -zeta, eta),
        5 => (eta, xi, -zeta),
        _ => (0.0, 0.0, 1.0),
    }
}

pub(crate) fn project_tsc(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let face = select_quadcube_face(phi, theta);

    if face.zeta <= 0.0 {
        return Err(WcsError::singularity(
            "TSC projection: point on back of cube face",
        ));
    }

    let chi = face.xi / face.zeta;
    let psi = face.eta / face.zeta;

    let x = face.phi_c * RAD_TO_DEG + 45.0 * chi;
    let y = face.theta_c * RAD_TO_DEG + 45.0 * psi;

    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_tsc(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (face, _phi_c, _theta_c, x_rel, y_rel) =
        quadcube_face_from_xy(inter.x_deg(), inter.y_deg());

    let chi = x_rel / 45.0;
    let psi = y_rel / 45.0;

    if chi.abs() > 1.0 || psi.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "TSC deprojection: point outside cube face",
        ));
    }

    let zeta = 1.0 / (1.0 + chi * chi + psi * psi).sqrt();
    let xi = chi * zeta;
    let eta = psi * zeta;

    let (l, m, n) = face_coords_to_direction_cosines(face, xi, eta, zeta);

    let theta = n.asin();
    let phi = m.atan2(l);

    Ok(native_coord_from_radians(phi, theta))
}

pub(crate) fn project_csc(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let face = select_quadcube_face(phi, theta);

    if face.zeta <= 0.0 {
        return Err(WcsError::singularity(
            "CSC projection: point on back of cube face",
        ));
    }

    let chi = face.xi / face.zeta;
    let psi = face.eta / face.zeta;

    let x_face = csc_forward_poly(chi, psi);
    let y_face = csc_forward_poly(psi, chi);

    let x = face.phi_c * RAD_TO_DEG + 45.0 * x_face;
    let y = face.theta_c * RAD_TO_DEG + 45.0 * y_face;

    Ok(IntermediateCoord::new(x, y))
}

fn csc_forward_poly(chi: f64, psi: f64) -> f64 {
    const GAMMA_STAR: f64 = 1.37484847732;
    const M: f64 = 0.004869491981;
    const GAMMA: f64 = -0.13161671474;
    const OMEGA1: f64 = -0.159596235474;
    const C00: f64 = 0.141189631152;
    const C10: f64 = 0.0809701286525;
    const C01: f64 = -0.281528535557;
    const C20: f64 = -0.178251207466;
    const C11: f64 = 0.15384112876;
    const C02: f64 = 0.106959469314;
    const D0: f64 = 0.0759196200467;
    const D1: f64 = -0.0217762490699;

    let chi2 = chi * chi;
    let psi2 = psi * psi;

    let c_poly =
        C00 + C10 * chi2 + C01 * psi2 + C20 * chi2 * chi2 + C11 * chi2 * psi2 + C02 * psi2 * psi2;

    let d_poly = D0 + D1 * chi2;

    let term1 = chi * GAMMA_STAR + chi * chi2 * (1.0 - GAMMA_STAR);
    let term2 = chi * psi2 * (1.0 - chi2) * (GAMMA + (M - GAMMA) * chi2 + (1.0 - psi2) * c_poly);
    let term3 = chi * chi2 * (1.0 - chi2) * (OMEGA1 - (1.0 - chi2) * d_poly);

    term1 + term2 + term3
}

pub(crate) fn deproject_csc(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (face, _phi_c, _theta_c, x_rel, y_rel) =
        quadcube_face_from_xy(inter.x_deg(), inter.y_deg());

    let x_norm = x_rel / 45.0;
    let y_norm = y_rel / 45.0;

    if x_norm.abs() > 1.0 || y_norm.abs() > 1.0 {
        return Err(WcsError::out_of_bounds(
            "CSC deprojection: point outside cube face",
        ));
    }

    let chi = csc_inverse_poly(x_norm, y_norm);
    let psi = csc_inverse_poly(y_norm, x_norm);

    let zeta = 1.0 / (1.0 + chi * chi + psi * psi).sqrt();
    let xi = chi * zeta;
    let eta = psi * zeta;

    let (l, m, n) = face_coords_to_direction_cosines(face, xi, eta, zeta);

    let theta = n.asin();
    let phi = m.atan2(l);

    Ok(native_coord_from_radians(phi, theta))
}

fn csc_inverse_poly(x: f64, y: f64) -> f64 {
    // CSC inverse polynomial coefficients from Paper II (Calabretta & Greisen 2002)
    // Section 5.6.2, page 24. P[i][j] corresponds to P_ij where the polynomial is:
    // f(X,Y) = X + X(1 - X²) Σᵢⱼ PᵢⱼX²ⁱY²ʲ
    // The matrix is organized as P[i][j] where i is X² power and j is Y² power.
    // Note: Paper II lists coefficients as P_ij where first subscript is X power.
    const P: [[f64; 7]; 7] = [
        // i=0: P_00, P_01, P_02, P_03, P_04, P_05, P_06
        [
            -0.27292696,
            -0.02819452,
            0.27058160,
            -0.60441560,
            0.93412077,
            -0.63915306,
            0.14381585,
        ],
        // i=1: P_10, P_11, P_12, P_13, P_14, P_15
        [
            -0.07629969,
            -0.01471565,
            -0.56800938,
            1.50880086,
            -1.41601920,
            0.52032238,
            0.0,
        ],
        // i=2: P_20, P_21, P_22, P_23, P_24
        [
            -0.22797056,
            0.48051509,
            0.30803317,
            -0.93678576,
            0.33887446,
            0.0,
            0.0,
        ],
        // i=3: P_30, P_31, P_32, P_33
        [
            0.54852384,
            -1.74114454,
            0.98938102,
            0.08693841,
            0.0,
            0.0,
            0.0,
        ],
        // i=4: P_40, P_41, P_42
        [-0.62930065, 1.71547508, -0.83180469, 0.0, 0.0, 0.0, 0.0],
        // i=5: P_50, P_51
        [0.25795794, -0.53022337, 0.0, 0.0, 0.0, 0.0, 0.0],
        // i=6: P_60
        [0.02584375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    let x2 = x * x;
    let y2 = y * y;

    let mut sum = 0.0;
    let mut x_pow = 1.0; // X^(2i) starts at X^0 = 1
    for (i, p_row) in P.iter().enumerate() {
        let mut y_pow = 1.0; // Y^(2j) starts at Y^0 = 1
        for p_val in p_row.iter().take(7 - i) {
            sum += p_val * x_pow * y_pow;
            y_pow *= y2;
        }
        x_pow *= x2;
    }

    x + x * (1.0 - x2) * sum
}

pub(crate) fn project_qsc(native: NativeCoord) -> WcsResult<IntermediateCoord> {
    let phi = native.phi().radians();
    let theta = native.theta().radians();

    let face = select_quadcube_face(phi, theta);

    if face.zeta <= 0.0 {
        return Err(WcsError::singularity(
            "QSC projection: point on back of cube face",
        ));
    }

    let xi_abs = face.xi.abs();
    let eta_abs = face.eta.abs();

    let (u, v) = if xi_abs >= eta_abs {
        let omega = if xi_abs > 1e-10 {
            face.eta / face.xi
        } else {
            0.0
        };
        let omega2 = omega * omega;

        let s = if face.xi >= 0.0 { 1.0 } else { -1.0 };
        let denom = 1.0 - 1.0 / 2.0_f64.sqrt() + omega2;
        let u_val = s * 45.0 * ((1.0 - face.zeta) / denom).sqrt();

        let v_val = if u_val.abs() > 1e-10 {
            let atan_omega = omega.atan() * RAD_TO_DEG;
            let asin_term = (omega / (2.0 * (1.0 + omega2)).sqrt()).asin() * RAD_TO_DEG;
            (u_val / 15.0) * (atan_omega - asin_term)
        } else {
            0.0
        };

        (u_val, v_val)
    } else {
        let omega = if eta_abs > 1e-10 {
            face.xi / face.eta
        } else {
            0.0
        };
        let omega2 = omega * omega;

        let s = if face.eta >= 0.0 { 1.0 } else { -1.0 };
        let denom = 1.0 - 1.0 / 2.0_f64.sqrt() + omega2;
        let u_val = s * 45.0 * ((1.0 - face.zeta) / denom).sqrt();

        let v_val = if u_val.abs() > 1e-10 {
            let atan_omega = omega.atan() * RAD_TO_DEG;
            let asin_term = (omega / (2.0 * (1.0 + omega2)).sqrt()).asin() * RAD_TO_DEG;
            (u_val / 15.0) * (atan_omega - asin_term)
        } else {
            0.0
        };

        (v_val, u_val)
    };

    let x = face.phi_c * RAD_TO_DEG + u;
    let y = face.theta_c * RAD_TO_DEG + v;

    Ok(IntermediateCoord::new(x, y))
}

pub(crate) fn deproject_qsc(inter: IntermediateCoord) -> WcsResult<NativeCoord> {
    let (face, _phi_c, _theta_c, x_rel, y_rel) =
        quadcube_face_from_xy(inter.x_deg(), inter.y_deg());

    let x_abs = x_rel.abs();
    let y_abs = y_rel.abs();

    if x_abs > 45.0 || y_abs > 45.0 {
        return Err(WcsError::out_of_bounds(
            "QSC deprojection: point outside cube face",
        ));
    }

    let (xi, eta, zeta) = if x_abs >= y_abs {
        let u = x_rel;
        let v = y_rel;

        let omega = qsc_inverse_omega(u, v);
        let omega2 = omega * omega;
        let denom = 1.0 - 1.0 / 2.0_f64.sqrt() + omega2;
        let zeta_val = 1.0 - (u / 45.0) * (u / 45.0) * denom;

        let factor = ((1.0 - zeta_val * zeta_val) / (1.0 + omega2)).sqrt();
        let xi_val = factor;
        let eta_val = omega * factor;

        let xi_signed = if x_rel >= 0.0 { xi_val } else { -xi_val };
        let eta_signed = if x_rel >= 0.0 { eta_val } else { -eta_val };

        (xi_signed, eta_signed, zeta_val)
    } else {
        let u = y_rel;
        let v = x_rel;

        let omega = qsc_inverse_omega(u, v);
        let omega2 = omega * omega;
        let denom = 1.0 - 1.0 / 2.0_f64.sqrt() + omega2;
        let zeta_val = 1.0 - (u / 45.0) * (u / 45.0) * denom;

        let factor = ((1.0 - zeta_val * zeta_val) / (1.0 + omega2)).sqrt();
        let eta_val = factor;
        let xi_val = omega * factor;

        let xi_signed = if y_rel >= 0.0 { xi_val } else { -xi_val };
        let eta_signed = if y_rel >= 0.0 { eta_val } else { -eta_val };

        (xi_signed, eta_signed, zeta_val)
    };

    let (l, m, n) = face_coords_to_direction_cosines(face, xi, eta, zeta);

    let theta = n.asin();
    let phi = m.atan2(l);

    Ok(native_coord_from_radians(phi, theta))
}

fn qsc_inverse_omega(u: f64, v: f64) -> f64 {
    if u.abs() < 1e-10 {
        return 0.0;
    }

    let ratio = 15.0 * v / u;
    let ratio_rad = ratio * DEG_TO_RAD;

    let cos_r = ratio_rad.cos();
    let denom = cos_r - 1.0 / 2.0_f64.sqrt();

    if denom.abs() < 1e-10 {
        return ratio_rad.sin() * 10.0;
    }

    ratio_rad.sin() / denom
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Projection;
    use celestial_core::assert_ulp_lt;
    use celestial_core::Angle;

    #[test]
    fn test_tsc_reference_point() {
        let proj = Projection::tsc();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_tsc_native_reference() {
        let proj = Projection::tsc();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_tsc_roundtrip() {
        let proj = Projection::tsc();
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_tsc_roundtrip_face_0() {
        let proj = Projection::tsc();
        let original = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_tsc_roundtrip_face_1() {
        let proj = Projection::tsc();
        let original = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_tsc_roundtrip_face_2() {
        let proj = Projection::tsc();
        let original = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_tsc_roundtrip_face_3() {
        let proj = Projection::tsc();
        let original = NativeCoord::new(Angle::from_degrees(180.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert!(
            (original.phi().degrees().abs() - recovered.phi().degrees().abs()).abs() < 1e-8,
            "phi mismatch: {} vs {}",
            original.phi().degrees(),
            recovered.phi().degrees()
        );
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_tsc_roundtrip_face_4() {
        let proj = Projection::tsc();
        let original = NativeCoord::new(Angle::from_degrees(-90.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_tsc_roundtrip_face_5() {
        let proj = Projection::tsc();
        let original = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 10);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 10);
    }

    #[test]
    fn test_tsc_roundtrip_various_angles() {
        let proj = Projection::tsc();
        for phi_deg in [-135.0, -90.0, -45.0, 0.0, 45.0, 90.0, 135.0] {
            for theta_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert!(
                    (original.phi().degrees() - recovered.phi().degrees()).abs() < 1e-8
                        || (original.phi().degrees().abs() - 180.0).abs() < 1e-8,
                    "phi mismatch at ({}, {}): {} vs {}",
                    phi_deg,
                    theta_deg,
                    original.phi().degrees(),
                    recovered.phi().degrees()
                );
                assert!(
                    (original.theta().degrees() - recovered.theta().degrees()).abs() < 1e-8,
                    "theta mismatch at ({}, {}): {} vs {}",
                    phi_deg,
                    theta_deg,
                    original.theta().degrees(),
                    recovered.theta().degrees()
                );
            }
        }
    }

    #[test]
    fn test_tsc_deproject_origin() {
        let proj = Projection::tsc();
        let inter = IntermediateCoord::new(0.0, 0.0);
        let result = proj.deproject(inter).unwrap();
        assert_eq!(result.phi().degrees(), 0.0);
        assert_eq!(result.theta().degrees(), 0.0);
    }

    #[test]
    fn test_tsc_pole() {
        let proj = Projection::tsc();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert_ulp_lt!(inter.y_deg(), 90.0, 2);
    }

    #[test]
    fn test_tsc_south_pole() {
        let proj = Projection::tsc();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-90.0));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert_ulp_lt!(inter.y_deg(), -90.0, 2);
    }

    #[test]
    fn test_csc_reference_point() {
        let proj = Projection::csc();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_csc_native_reference() {
        let proj = Projection::csc();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_csc_roundtrip() {
        let proj = Projection::csc();
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert!(
            (original.phi().degrees() - recovered.phi().degrees()).abs() < 0.01,
            "phi mismatch: {} vs {}",
            original.phi().degrees(),
            recovered.phi().degrees()
        );
        assert!(
            (original.theta().degrees() - recovered.theta().degrees()).abs() < 0.01,
            "theta mismatch: {} vs {}",
            original.theta().degrees(),
            recovered.theta().degrees()
        );
    }

    #[test]
    fn test_csc_roundtrip_face_0() {
        let proj = Projection::csc();
        let original = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert!(
            (original.phi().degrees() - recovered.phi().degrees()).abs() < 0.01,
            "phi mismatch: {} vs {}",
            original.phi().degrees(),
            recovered.phi().degrees()
        );
        assert!(
            (original.theta().degrees() - recovered.theta().degrees()).abs() < 0.01,
            "theta mismatch: {} vs {}",
            original.theta().degrees(),
            recovered.theta().degrees()
        );
    }

    #[test]
    fn test_csc_roundtrip_face_1() {
        let proj = Projection::csc();
        let original = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert!(
            (original.phi().degrees() - recovered.phi().degrees()).abs() < 0.01,
            "phi mismatch: {} vs {}",
            original.phi().degrees(),
            recovered.phi().degrees()
        );
        assert!(
            (original.theta().degrees() - recovered.theta().degrees()).abs() < 0.01,
            "theta mismatch: {} vs {}",
            original.theta().degrees(),
            recovered.theta().degrees()
        );
    }

    #[test]
    fn test_csc_roundtrip_face_2() {
        let proj = Projection::csc();
        let original = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert!(
            (original.phi().degrees() - recovered.phi().degrees()).abs() < 0.01,
            "phi mismatch: {} vs {}",
            original.phi().degrees(),
            recovered.phi().degrees()
        );
        assert!(
            (original.theta().degrees() - recovered.theta().degrees()).abs() < 0.01,
            "theta mismatch: {} vs {}",
            original.theta().degrees(),
            recovered.theta().degrees()
        );
    }

    #[test]
    fn test_csc_roundtrip_face_4() {
        let proj = Projection::csc();
        let original = NativeCoord::new(Angle::from_degrees(-90.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert!(
            (original.phi().degrees() - recovered.phi().degrees()).abs() < 0.01,
            "phi mismatch: {} vs {}",
            original.phi().degrees(),
            recovered.phi().degrees()
        );
        assert!(
            (original.theta().degrees() - recovered.theta().degrees()).abs() < 0.01,
            "theta mismatch: {} vs {}",
            original.theta().degrees(),
            recovered.theta().degrees()
        );
    }

    #[test]
    fn test_csc_roundtrip_face_5() {
        let proj = Projection::csc();
        let original = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert!(
            (original.phi().degrees() - recovered.phi().degrees()).abs() < 0.01,
            "phi mismatch: {} vs {}",
            original.phi().degrees(),
            recovered.phi().degrees()
        );
        assert!(
            (original.theta().degrees() - recovered.theta().degrees()).abs() < 0.01,
            "theta mismatch: {} vs {}",
            original.theta().degrees(),
            recovered.theta().degrees()
        );
    }

    #[test]
    fn test_csc_roundtrip_various_angles() {
        let proj = Projection::csc();
        for phi_deg in [-90.0, -45.0, 0.0, 45.0, 90.0] {
            for theta_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert!(
                    (original.phi().degrees() - recovered.phi().degrees()).abs() < 0.01,
                    "phi mismatch at ({}, {}): {} vs {}",
                    phi_deg,
                    theta_deg,
                    original.phi().degrees(),
                    recovered.phi().degrees()
                );
                assert!(
                    (original.theta().degrees() - recovered.theta().degrees()).abs() < 0.01,
                    "theta mismatch at ({}, {}): {} vs {}",
                    phi_deg,
                    theta_deg,
                    original.theta().degrees(),
                    recovered.theta().degrees()
                );
            }
        }
    }

    #[test]
    fn test_csc_deproject_origin() {
        let proj = Projection::csc();
        let inter = IntermediateCoord::new(0.0, 0.0);
        let result = proj.deproject(inter).unwrap();
        assert_eq!(result.phi().degrees(), 0.0);
        assert_eq!(result.theta().degrees(), 0.0);
    }

    #[test]
    fn test_csc_pole() {
        let proj = Projection::csc();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert_ulp_lt!(inter.y_deg(), 90.0, 2);
    }

    #[test]
    fn test_qsc_reference_point() {
        let proj = Projection::qsc();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert!(inter.y_deg().abs() < 1e-10);
    }

    #[test]
    fn test_qsc_native_reference() {
        let proj = Projection::qsc();
        let (phi0, theta0) = proj.native_reference();
        assert_eq!(phi0, 0.0);
        assert_eq!(theta0, 0.0);
    }

    #[test]
    fn test_qsc_roundtrip() {
        let proj = Projection::qsc();
        let original = NativeCoord::new(Angle::from_degrees(30.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 30);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 30);
    }

    #[test]
    fn test_qsc_roundtrip_face_0() {
        let proj = Projection::qsc();
        let original = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 30);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 30);
    }

    #[test]
    fn test_qsc_roundtrip_face_1() {
        let proj = Projection::qsc();
        let original = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 30);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 30);
    }

    #[test]
    fn test_qsc_roundtrip_face_2() {
        let proj = Projection::qsc();
        let original = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 30);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 30);
    }

    #[test]
    fn test_qsc_roundtrip_face_4() {
        let proj = Projection::qsc();
        let original = NativeCoord::new(Angle::from_degrees(-90.0), Angle::from_degrees(30.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 30);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 30);
    }

    #[test]
    fn test_qsc_roundtrip_face_5() {
        let proj = Projection::qsc();
        let original = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-60.0));
        let inter = proj.project(original).unwrap();
        let recovered = proj.deproject(inter).unwrap();
        assert_ulp_lt!(original.phi().degrees(), recovered.phi().degrees(), 30);
        assert_ulp_lt!(original.theta().degrees(), recovered.theta().degrees(), 30);
    }

    #[test]
    fn test_qsc_roundtrip_various_angles() {
        let proj = Projection::qsc();
        for phi_deg in [-90.0, -45.0, 0.0, 45.0, 90.0] {
            for theta_deg in [-60.0, -30.0, 0.0, 30.0, 60.0] {
                let original =
                    NativeCoord::new(Angle::from_degrees(phi_deg), Angle::from_degrees(theta_deg));
                let inter = proj.project(original).unwrap();
                let recovered = proj.deproject(inter).unwrap();
                assert!(
                    (original.phi().degrees() - recovered.phi().degrees()).abs() < 1e-6,
                    "phi mismatch at ({}, {}): {} vs {}",
                    phi_deg,
                    theta_deg,
                    original.phi().degrees(),
                    recovered.phi().degrees()
                );
                assert!(
                    (original.theta().degrees() - recovered.theta().degrees()).abs() < 1e-6,
                    "theta mismatch at ({}, {}): {} vs {}",
                    phi_deg,
                    theta_deg,
                    original.theta().degrees(),
                    recovered.theta().degrees()
                );
            }
        }
    }

    #[test]
    fn test_qsc_deproject_origin() {
        let proj = Projection::qsc();
        let inter = IntermediateCoord::new(0.0, 0.0);
        let result = proj.deproject(inter).unwrap();
        assert_eq!(result.phi().degrees(), 0.0);
        assert_eq!(result.theta().degrees(), 0.0);
    }

    #[test]
    fn test_qsc_pole() {
        let proj = Projection::qsc();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(90.0));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert_ulp_lt!(inter.y_deg(), 90.0, 2);
    }

    #[test]
    fn test_qsc_south_pole() {
        let proj = Projection::qsc();
        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-90.0));
        let inter = proj.project(native).unwrap();
        assert!(inter.x_deg().abs() < 1e-10);
        assert_ulp_lt!(inter.y_deg(), -90.0, 2);
    }

    #[test]
    fn test_quadcube_projections_native_reference() {
        let projections = [Projection::tsc(), Projection::csc(), Projection::qsc()];

        for proj in projections {
            let (phi0, theta0) = proj.native_reference();
            assert_eq!(phi0, 0.0);
            assert_eq!(theta0, 0.0);
        }
    }

    #[test]
    fn test_quadcube_projections_reference_maps_to_origin() {
        let projections = [Projection::tsc(), Projection::csc(), Projection::qsc()];

        for proj in &projections {
            let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(0.0));
            let inter = proj.project(native).unwrap();
            assert!(inter.x_deg().abs() < 1e-10, "x not zero for {:?}", proj);
            assert!(inter.y_deg().abs() < 1e-10, "y not zero for {:?}", proj);
        }
    }

    #[test]
    fn test_quadcube_face_selection() {
        let proj = Projection::tsc();

        let face_0 = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(60.0));
        let inter_0 = proj.project(face_0).unwrap();
        assert!(inter_0.y_deg() > 45.0);

        let face_1 = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(30.0));
        let inter_1 = proj.project(face_1).unwrap();
        assert!(inter_1.x_deg().abs() < 45.0 && inter_1.y_deg().abs() < 45.0);

        let face_2 = NativeCoord::new(Angle::from_degrees(90.0), Angle::from_degrees(30.0));
        let inter_2 = proj.project(face_2).unwrap();
        assert!(inter_2.x_deg() > 45.0 && inter_2.x_deg() < 135.0);

        let face_5 = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(-60.0));
        let inter_5 = proj.project(face_5).unwrap();
        assert!(inter_5.y_deg() < -45.0);
    }

    #[test]
    fn test_tsc_vs_tan_at_face_center() {
        let tsc = Projection::tsc();
        let tan = Projection::tan();

        let native = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(60.0));
        let tsc_inter = tsc.project(native).unwrap();

        let native_for_tan = NativeCoord::new(Angle::from_degrees(0.0), Angle::from_degrees(60.0));
        let tan_inter = tan.project(native_for_tan).unwrap();

        assert!(tsc_inter.y_deg() > 45.0);
        assert!(tan_inter.y_deg() < 0.0);
    }
}
