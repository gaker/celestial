use crate::constants::{ARCSEC_TO_RAD, J2000_JD};
use crate::matrix::RotationMatrix3;
use crate::precession::iau2006::fw_matrix;
use crate::precession::{FukushimaWilliamsAngles, PrecessionIAU2006};

#[test]
fn test_compute_rejects_epochs_outside_the_model_range() {
    let p = PrecessionIAU2006::new();
    let beyond = 20.5 * crate::constants::DAYS_PER_JULIAN_CENTURY;
    assert!(p.compute(J2000_JD, beyond).is_err());
    assert!(p.compute(J2000_JD, -beyond).is_err());
    assert!(p.compute(f64::NAN, 0.0).is_err());
}

#[test]
fn test_compute_returns_rotation_matrices() {
    let p = PrecessionIAU2006::new();
    let result = p
        .compute(J2000_JD, 0.5 * crate::constants::DAYS_PER_JULIAN_CENTURY)
        .unwrap();
    assert!(result.bias_matrix.is_rotation_matrix(1e-14));
    assert!(result.precession_matrix.is_rotation_matrix(1e-14));
    assert!(result.bias_precession_matrix.is_rotation_matrix(1e-14));
}

#[test]
fn test_bias_matrix_is_constant() {
    let p = PrecessionIAU2006::new();
    let r1 = p.compute(J2000_JD, 0.0).unwrap();
    let r2 = p
        .compute(J2000_JD, crate::constants::DAYS_PER_JULIAN_CENTURY)
        .unwrap();
    assert_eq!(r1.bias_matrix, r2.bias_matrix);
}

// Expected values are ERFA bp06 outputs with Rust libm for sin/cos. At J2000.0 the
// precession matrix is B·Bᵀ, which rounding keeps from being exactly the identity.
#[test]
fn test_precession_at_j2000_matches_erfa_bp06() {
    let p = PrecessionIAU2006::new();
    let result = p.compute(J2000_JD, 0.0).unwrap();
    let expected = RotationMatrix3::from_array([
        [
            0.9999999999999997,
            2.917805702402595e-23,
            1.3234889800848443e-23,
        ],
        [
            2.917805702402595e-23,
            0.9999999999999999,
            -4.034804386554417e-17,
        ],
        [1.3234889800848443e-23, -4.034804386554417e-17, 1.0],
    ])
    .unwrap();
    assert_eq!(result.precession_matrix, expected);
}

#[test]
fn test_precession_changes_with_time() {
    let p = PrecessionIAU2006::new();
    let r0 = p.compute(J2000_JD, 0.0).unwrap();
    let r1 = p
        .compute(J2000_JD, crate::constants::DAYS_PER_JULIAN_CENTURY)
        .unwrap();
    let identity = RotationMatrix3::identity();
    assert!(
        r1.precession_matrix.max_difference(&identity) > 1e-4,
        "precession matrix should differ from identity after 1 century"
    );
    assert_ne!(r0.precession_matrix, r1.precession_matrix);
}

#[test]
fn test_combined_matrix_consistency() {
    let p = PrecessionIAU2006::new();
    let result = p
        .compute(J2000_JD, 0.5 * crate::constants::DAYS_PER_JULIAN_CENTURY)
        .unwrap();
    let bias_inverse = result.bias_matrix.transpose();
    let expected = result.bias_precession_matrix.multiply(&bias_inverse);
    assert_eq!(result.precession_matrix, expected);
}

// Expected values are ERFA pfw06 outputs at J2000.0.
#[test]
fn test_fukushima_williams_angles_match_erfa_pfw06_at_j2000() {
    let p = PrecessionIAU2006::new();
    assert_eq!(
        p.fukushima_williams_angles(0.0),
        FukushimaWilliamsAngles {
            gamma_bar: -2.5660218513765524e-7,
            phi_bar: 0.4090926336600278,
            psi_bar: -2.0253091528350866e-7,
            epsilon_a: 0.4090926006005829,
        }
    );
}

#[test]
fn test_fukushima_williams_angles_change_with_time() {
    let p = PrecessionIAU2006::new();
    let fw0 = p.fukushima_williams_angles(0.0);
    let fw1 = p.fukushima_williams_angles(1.0);
    assert_ne!(fw0.gamma_bar, fw1.gamma_bar);
    assert_ne!(fw0.phi_bar, fw1.phi_bar);
    assert_ne!(fw0.psi_bar, fw1.psi_bar);
    assert_ne!(fw0.epsilon_a, fw1.epsilon_a);
}

#[test]
fn test_fw_angles_to_matrix_returns_rotation() {
    let p = PrecessionIAU2006::new();
    let matrix = fw_matrix(&p.fukushima_williams_angles(0.5));
    assert!(matrix.is_rotation_matrix(1e-14));
}

#[test]
fn test_npb_matrix_returns_rotation() {
    let p = PrecessionIAU2006::new();
    let matrix = p.npb_matrix_iau2006a(0.5, 0.001 * ARCSEC_TO_RAD, 0.0005 * ARCSEC_TO_RAD);
    assert!(matrix.is_rotation_matrix(1e-14));
}

#[test]
fn test_npb_matrix_with_zero_nutation() {
    let p = PrecessionIAU2006::new();
    let fw = fw_matrix(&p.fukushima_williams_angles(0.5));
    let npb_matrix = p.npb_matrix_iau2006a(0.5, 0.0, 0.0);
    assert_eq!(npb_matrix, fw);
}

#[test]
fn test_two_part_date_equivalence() {
    let p = PrecessionIAU2006::new();
    let r1 = p.compute(J2000_JD, 1000.0).unwrap();
    let r2 = p.compute(J2000_JD + 500.0, 500.0).unwrap();
    assert_eq!(r1.bias_precession_matrix, r2.bias_precession_matrix);
}
