use crate::constants::{DAYS_PER_JULIAN_CENTURY, J2000_JD};
use crate::matrix::RotationMatrix3;
use crate::precession::{PrecessionIAU2000, PrecessionResult};

#[test]
fn test_compute_rejects_epochs_outside_the_model_range() {
    let p = PrecessionIAU2000::new();
    assert!(p.compute(J2000_JD, 20.0 * DAYS_PER_JULIAN_CENTURY).is_ok());
    assert!(p.compute(J2000_JD, 20.5 * DAYS_PER_JULIAN_CENTURY).is_err());
    assert!(p
        .compute(J2000_JD, -20.5 * DAYS_PER_JULIAN_CENTURY)
        .is_err());
    assert!(p.compute(f64::NAN, 0.0).is_err());
}

#[test]
fn test_compute_returns_rotation_matrices() {
    let p = PrecessionIAU2000::new();
    let result = p.compute(J2000_JD, 0.5 * DAYS_PER_JULIAN_CENTURY).unwrap();
    assert!(result.bias_matrix.is_rotation_matrix(1e-14));
    assert!(result.precession_matrix.is_rotation_matrix(1e-14));
    assert!(result.bias_precession_matrix.is_rotation_matrix(1e-14));
}

#[test]
fn test_bias_matrix_is_constant() {
    let p = PrecessionIAU2000::new();
    let r1 = p.compute(J2000_JD, 0.0).unwrap();
    let r2 = p.compute(J2000_JD, DAYS_PER_JULIAN_CENTURY).unwrap();
    assert_eq!(r1.bias_matrix, r2.bias_matrix);
}

#[test]
fn test_precession_at_j2000_is_identity() {
    let p = PrecessionIAU2000::new();
    let result = p.compute(J2000_JD, 0.0).unwrap();
    assert_eq!(result.precession_matrix, RotationMatrix3::identity());
}

#[test]
fn test_precession_changes_with_time() {
    let p = PrecessionIAU2000::new();
    let r0 = p.compute(J2000_JD, 0.0).unwrap();
    let r1 = p.compute(J2000_JD, DAYS_PER_JULIAN_CENTURY).unwrap();
    let identity = RotationMatrix3::identity();
    assert!(
        r1.precession_matrix.max_difference(&identity) > 1e-4,
        "precession matrix should differ from identity after 1 century"
    );
    assert_ne!(r0.precession_matrix, r1.precession_matrix);
}

// Expected values are ERFA bp00 outputs with Rust libm for sin/cos. The literals printed
// in t_erfa_c.c are not bit-reproducible by ERFA itself.
fn bp00_result(jd1: f64, jd2: f64) -> PrecessionResult {
    PrecessionIAU2000::new().compute(jd1, jd2).unwrap()
}

fn erfa_bp00_result() -> PrecessionResult {
    bp00_result(2400000.5, 50123.9999)
}

const ERFA_BP00_RBP: [(f64, f64, [[f64; 3]; 3]); 2] = [
    (
        2451545.0,
        36525.0,
        [
            [
                0.9997026830133443,
                -0.02236501297096165,
                -0.009713483964539135,
            ],
            [
                0.022365014311294465,
                0.9997498658978478,
                -0.00010849926694337283,
            ],
            [
                0.009713480878461423,
                -0.00010877519961119612,
                0.9999528171154775,
            ],
        ],
    ),
    (
        2451545.0,
        -73050.0,
        [
            [0.9988124884471957, 0.04467592415373793, 0.01943385507123954],
            [
                -0.04467592289443627,
                0.9990014380266931,
                -0.00043443541478397854,
            ],
            [
                -0.019433857966211224,
                -0.0004343058929944033,
                0.9998110504204984,
            ],
        ],
    ),
];

#[test]
fn test_bias_precession_matrix_matches_erfa_bp00_across_epochs() {
    for (jd1, jd2, expected) in ERFA_BP00_RBP {
        let got = bp00_result(jd1, jd2).bias_precession_matrix;
        assert_eq!(
            got,
            RotationMatrix3::from_array(expected).unwrap(),
            "{jd1} + {jd2}"
        );
    }
}

#[test]
fn test_bias_matrix_matches_erfa_bp00() {
    let expected = RotationMatrix3::from_array([
        [
            0.9999999999999942,
            -7.078279744199198e-8,
            8.056217146976134e-8,
        ],
        [
            7.078279477857338e-8,
            0.9999999999999969,
            3.3060414542221364e-8,
        ],
        [
            -8.056217380986972e-8,
            -3.306040883980552e-8,
            0.9999999999999962,
        ],
    ])
    .unwrap();
    assert_eq!(erfa_bp00_result().bias_matrix, expected);
}

#[test]
fn test_precession_matrix_matches_erfa_bp00() {
    let expected = RotationMatrix3::from_array([
        [
            0.9999995504864049,
            0.0008696113836207071,
            0.0003778928813389328,
        ],
        [
            -0.0008696113818227252,
            0.9999996218879366,
            -1.6906792630454113e-7,
        ],
        [
            -0.0003778928854764689,
            -1.595521004205125e-7,
            0.9999999285984683,
        ],
    ])
    .unwrap();
    assert_eq!(erfa_bp00_result().precession_matrix, expected);
}

#[test]
fn test_bias_precession_matrix_matches_erfa_bp00() {
    let expected = RotationMatrix3::from_array([
        [
            0.9999995505175088,
            0.000869540588361787,
            0.00037797347222390015,
        ],
        [
            -0.000869540599041085,
            0.9999996219494925,
            -1.3607758204411515e-7,
        ],
        [
            -0.0003779734476558179,
            -1.925857585841862e-7,
            0.9999999285680153,
        ],
    ])
    .unwrap();
    assert_eq!(erfa_bp00_result().bias_precession_matrix, expected);
}

#[test]
fn test_combined_matrix_consistency() {
    let p = PrecessionIAU2000::new();
    let result = p.compute(J2000_JD, 0.5 * DAYS_PER_JULIAN_CENTURY).unwrap();
    let expected = result.precession_matrix.multiply(&result.bias_matrix);
    assert_eq!(result.bias_precession_matrix, expected);
}
