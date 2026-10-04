use super::*;
use crate::constants::HALF_PI;
use crate::errors::{AstroError, MathErrorKind};
use crate::matrix::Vector3;
use std::ops::Mul;

// f64 π/2 sits just below the true value, so its cosine is tiny but not zero.
const COS_HALF_PI: f64 = 6.123233995736766e-17;

fn unchecked(elements: [[f64; 3]; 3]) -> RotationMatrix3 {
    RotationMatrix3 { elements }
}

fn identity_with(row: usize, col: usize, value: f64) -> RotationMatrix3 {
    let mut elements = *RotationMatrix3::identity().elements();
    elements[row][col] = value;
    unchecked(elements)
}

fn from_array_error_kind(elements: [[f64; 3]; 3]) -> MathErrorKind {
    match RotationMatrix3::from_array(elements) {
        Err(AstroError::MathError { kind, .. }) => kind,
        other => panic!("expected a math error, got {other:?}"),
    }
}

// Expected values are an ERFA pnm06a output.
const ERFA_NPB: [[f64; 3]; 3] = [
    [
        0.9999989440476104,
        -0.0013328817612400115,
        -0.0005790767434730085,
    ],
    [
        0.0013328582543089545,
        0.9999991109044506,
        -4.097782710401556e-5,
    ],
    [
        0.0005791308472168153,
        4.0205956615939944e-5,
        0.9999998314954572,
    ],
];

#[test]
fn test_identity() {
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(RotationMatrix3::identity().elements(), &identity);
}

#[test]
fn test_from_array_accepts_erfa_matrix() {
    let m = RotationMatrix3::from_array(ERFA_NPB).unwrap();
    assert_eq!(m.elements(), &ERFA_NPB);
}

#[test]
fn test_from_array_rejects_non_rotations() {
    use MathErrorKind::{InvalidInput, NotFinite};
    let scaled = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]];
    let reflection = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]];
    let sheared = [[1.0, 1e-11, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(from_array_error_kind(scaled), InvalidInput);
    assert_eq!(from_array_error_kind(reflection), InvalidInput);
    assert_eq!(from_array_error_kind(sheared), InvalidInput);
    assert_eq!(from_array_error_kind([[0.0; 3]; 3]), InvalidInput);
    assert_eq!(from_array_error_kind([[f64::NAN; 3]; 3]), NotFinite);
    let mut infinite = *RotationMatrix3::identity().elements();
    infinite[1][2] = f64::INFINITY;
    assert_eq!(from_array_error_kind(infinite), NotFinite);
}

#[cfg(feature = "serde")]
#[test]
fn test_deserialize_validates() {
    let m = RotationMatrix3::from_array(ERFA_NPB).unwrap();
    let json = serde_json::to_string(&m).unwrap();
    assert_eq!(serde_json::from_str::<RotationMatrix3>(&json).unwrap(), m);
    let scaled = r#"{"elements": [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]}"#;
    assert!(serde_json::from_str::<RotationMatrix3>(scaled).is_err());
}

#[test]
fn test_rotate_z() {
    // ERFA convention: Rz(+psi) rotates anticlockwise looking from +z toward origin
    // This means [1,0,0] -> [cos(psi), -sin(psi), 0]
    // At 90°: [1,0,0] -> [0, -1, 0]
    let mut m = RotationMatrix3::identity();
    m.rotate_z(HALF_PI);
    let result = m.apply_to_vector([1.0, 0.0, 0.0]);
    assert_eq!(result, [COS_HALF_PI, -1.0, 0.0]);
}

#[test]
fn test_rotate_x() {
    // ERFA convention: Rx(+phi) rotates anticlockwise looking from +x toward origin
    // This means [0,1,0] -> [0, cos(phi), -sin(phi)]
    // At 90°: [0,1,0] -> [0, 0, -1]
    let mut m = RotationMatrix3::identity();
    m.rotate_x(HALF_PI);
    let result = m.apply_to_vector([0.0, 1.0, 0.0]);
    assert_eq!(result, [0.0, COS_HALF_PI, -1.0]);
}

#[test]
fn test_rotate_y() {
    // ERFA convention: Ry(+theta) rotates anticlockwise looking from +y toward origin
    // This means [0,0,1] -> [-sin(theta), 0, cos(theta)]
    // At 90°: [0,0,1] -> [-1, 0, 0]
    let mut m = RotationMatrix3::identity();
    m.rotate_y(HALF_PI);
    let result = m.apply_to_vector([0.0, 0.0, 1.0]);
    assert_eq!(result, [-1.0, 0.0, COS_HALF_PI]);
}

// ERFA rxp starts each row sum at 0.0, so negative-zero products sum to +0.0.
// The sign matters downstream: atan2(-0.0, -1.0) is -π, atan2(0.0, -1.0) is π.
#[test]
fn test_apply_to_vector_matches_erfa_rxp_signed_zero() {
    let out = RotationMatrix3::identity().apply_to_vector([-1.0, -0.0, -0.0]);
    assert_eq!(out.map(f64::to_bits), [-1.0, 0.0, 0.0].map(f64::to_bits));
}

#[test]
fn test_is_rotation_matrix_valid() {
    let mut m = RotationMatrix3::identity();
    m.rotate_z(0.5);
    assert!(m.is_rotation_matrix(1e-14));
}

#[test]
fn test_is_rotation_matrix_bad_determinant() {
    let m = identity_with(0, 0, 2.0);
    assert!(!m.is_rotation_matrix(1e-15));
}

#[test]
fn test_is_rotation_matrix_not_orthogonal() {
    let m = identity_with(0, 1, 0.1);
    assert!(!m.is_rotation_matrix(1e-15));
}

#[test]
fn test_is_rotation_matrix_rejects_nan() {
    let all_nan = unchecked([[f64::NAN; 3]; 3]);
    assert!(!all_nan.is_rotation_matrix(1e-14));
    let one_nan = identity_with(2, 2, f64::NAN);
    assert!(!one_nan.is_rotation_matrix(1e-14));
    assert!(!RotationMatrix3::identity().is_rotation_matrix(f64::NAN));
}

#[test]
fn test_max_difference_propagates_nan() {
    let identity = RotationMatrix3::identity();
    let first = identity_with(0, 0, f64::NAN);
    assert!(identity.max_difference(&first).is_nan());
    let last = identity_with(2, 2, f64::NAN);
    assert!(identity.max_difference(&last).is_nan());
}

// Expected values are eraS2c -> eraRxp -> eraC2s with ERFA's sin/cos/atan2 bound to
// Rust libm.
#[test]
fn test_transform_spherical_matches_erfa_near_pole() {
    let out = RotationMatrix3::identity().transform_spherical(0.3, HALF_PI - 1e-9);
    assert_eq!(out, (0.29999999999999993, 1.5707963257948965));
}

#[test]
fn test_transform_spherical_matches_erfa() {
    let mut m = RotationMatrix3::identity();
    m.rotate_z(0.4);
    m.rotate_x(0.3);
    assert_eq!(
        m.transform_spherical(2.0, -0.7),
        (1.6121309352387712, -0.9998216488211331)
    );
    assert_eq!(
        m.transform_spherical(0.25, 0.25),
        (-0.06796450854960329, 0.282901650627672)
    );
}

#[test]
fn test_transform_spherical_identity() {
    let m = RotationMatrix3::identity();
    assert_eq!(m.transform_spherical(1.0, 0.5), (1.0, 0.5));
}

#[test]
fn test_transform_spherical_rotation() {
    // ERFA Rz rotates in opposite direction to naive expectation
    // Rz(+90°) takes RA=0 to RA=-90° (or equivalently RA=270°=-HALF_PI)
    let mut m = RotationMatrix3::identity();
    m.rotate_z(HALF_PI);
    assert_eq!(m.transform_spherical(0.0, 0.0), (-HALF_PI, 0.0));
}

#[test]
fn test_transform_spherical_zero_norm() {
    let zero_matrix = unchecked([[0.0; 3]; 3]);
    let (_, dec) = zero_matrix.transform_spherical(0.0, 0.0);
    assert!(dec.is_finite());
}

// Each Mul overload is called by name so that every impl runs; written as operators,
// clippy::op_ref would flag the borrowed forms as redundant.
#[test]
fn test_mul_matrix_matrix() {
    let mut a = RotationMatrix3::identity();
    a.rotate_x(0.1);
    let mut b = RotationMatrix3::identity();
    b.rotate_y(0.2);

    let expected = a.multiply(&b);
    assert_eq!(Mul::mul(a, b), expected);
    assert_eq!(Mul::mul(a, &b), expected);
    assert_eq!(Mul::mul(&a, b), expected);
    assert_eq!(Mul::mul(&a, &b), expected);
}

#[test]
fn test_mul_matrix_vector() {
    let m = RotationMatrix3::identity();
    let v = Vector3::new(1.0, 2.0, 3.0);
    assert_eq!(Mul::mul(m, v), v);
    assert_eq!(Mul::mul(&m, v), v);
}

#[test]
fn test_display() {
    let mut m = RotationMatrix3::identity();
    m.rotate_z(0.1);
    let s = format!("{}", m);
    assert!(s.contains("RotationMatrix3:"));
    assert!(s.contains("["));
}

#[test]
fn test_max_difference() {
    let a = RotationMatrix3::identity();
    let b = identity_with(0, 1, 0.1);
    assert_eq!(a.max_difference(&b), 0.1);
}

#[test]
fn test_elements() {
    let m = RotationMatrix3::identity();
    let e = m.elements();
    assert_eq!(e[0][0], 1.0);
    assert_eq!(e[1][1], 1.0);
}
