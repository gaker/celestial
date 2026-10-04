use super::*;
use crate::constants::HALF_PI;

// f64 π/2 sits just below the true value, so its cosine is tiny but not zero.
const COS_HALF_PI: f64 = 6.123233995736766e-17;

#[test]
fn test_vector3_construction() {
    let v = Vector3::new(1.0, 2.0, 3.0);
    assert_eq!(v.x, 1.0);
    assert_eq!(v.y, 2.0);
    assert_eq!(v.z, 3.0);

    let zeros = Vector3::zeros();
    assert_eq!(zeros.x, 0.0);
    assert_eq!(zeros.y, 0.0);
    assert_eq!(zeros.z, 0.0);

    let x_axis = Vector3::x_axis();
    assert_eq!(x_axis, Vector3::new(1.0, 0.0, 0.0));

    let from_array = Vector3::from_array([4.0, 5.0, 6.0]);
    assert_eq!(from_array, Vector3::new(4.0, 5.0, 6.0));
}

#[test]
fn test_vector3_magnitude() {
    let v = Vector3::new(3.0, 4.0, 0.0);
    assert_eq!(v.magnitude(), 5.0);
    assert_eq!(v.magnitude_squared(), 25.0);

    let unit = v.normalize().unwrap();
    assert_eq!(unit.magnitude(), 1.0);
    assert_eq!(unit, Vector3::new(0.6000000000000001, 0.8, 0.0));
}

#[test]
fn test_vector3_arithmetic() {
    let a = Vector3::new(1.0, 2.0, 3.0);
    let b = Vector3::new(4.0, 5.0, 6.0);

    let sum = a + b;
    assert_eq!(sum, Vector3::new(5.0, 7.0, 9.0));

    let diff = b - a;
    assert_eq!(diff, Vector3::new(3.0, 3.0, 3.0));

    let scaled = a * 2.0;
    assert_eq!(scaled, Vector3::new(2.0, 4.0, 6.0));

    let scaled2 = 3.0 * a;
    assert_eq!(scaled2, Vector3::new(3.0, 6.0, 9.0));

    let divided = a / 2.0;
    assert_eq!(divided, Vector3::new(0.5, 1.0, 1.5));

    let negated = -a;
    assert_eq!(negated, Vector3::new(-1.0, -2.0, -3.0));
}

#[test]
fn test_vector3_dot_cross() {
    let a = Vector3::new(1.0, 0.0, 0.0);
    let b = Vector3::new(0.0, 1.0, 0.0);

    assert_eq!(a.dot(&b), 0.0);

    let c = a.cross(&b);
    assert_eq!(c, Vector3::new(0.0, 0.0, 1.0));

    let d = Vector3::new(1.0, 2.0, 3.0);
    let e = Vector3::new(4.0, 5.0, 6.0);
    assert_eq!(d.dot(&e), 32.0);
}

#[test]
fn test_vector3_spherical_conversion() {
    let v1 = Vector3::from_spherical(0.0, 0.0);
    assert_eq!(v1, Vector3::new(1.0, 0.0, 0.0));
    assert_eq!(v1.to_spherical(), (0.0, 0.0));

    let v2 = Vector3::from_spherical(HALF_PI, 0.0);
    assert_eq!(v2, Vector3::new(COS_HALF_PI, 1.0, 0.0));

    let v3 = Vector3::from_spherical(0.0, HALF_PI);
    assert_eq!(v3, Vector3::new(COS_HALF_PI, 0.0, 1.0));
}

#[test]
fn test_axis_constructors() {
    let y_axis = Vector3::y_axis();
    assert_eq!(y_axis, Vector3::new(0.0, 1.0, 0.0));

    let z_axis = Vector3::z_axis();
    assert_eq!(z_axis, Vector3::new(0.0, 0.0, 1.0));
}

#[test]
fn test_get_set_methods() {
    let mut v = Vector3::new(1.0, 2.0, 3.0);

    assert_eq!(v.get(0).unwrap(), 1.0);
    assert_eq!(v.get(1).unwrap(), 2.0);
    assert_eq!(v.get(2).unwrap(), 3.0);

    v.set(0, 10.0).unwrap();
    v.set(1, 20.0).unwrap();
    v.set(2, 30.0).unwrap();
    assert_eq!(v, Vector3::new(10.0, 20.0, 30.0));
}

#[test]
fn test_get_error() {
    let v = Vector3::new(1.0, 2.0, 3.0);
    let result = v.get(3);
    assert!(result.is_err());

    if let Err(err) = result {
        assert!(err.to_string().contains("index 3 out of bounds"));
    }
}

#[test]
fn test_set_error() {
    let mut v = Vector3::new(1.0, 2.0, 3.0);
    let result = v.set(5, 42.0);
    assert!(result.is_err());

    if let Err(err) = result {
        assert!(err.to_string().contains("index 5 out of bounds"));
    }
}

#[test]
fn test_normalize_rejects_zero_and_non_finite_length() {
    use MathErrorKind::{DivisionByZero, NotFinite};
    let kind = |v: Vector3| match v.normalize() {
        Err(AstroError::MathError { kind, .. }) => kind,
        other => panic!("expected a math error, got {other:?}"),
    };
    assert_eq!(kind(Vector3::zeros()), DivisionByZero);
    assert_eq!(kind(Vector3::new(1e-200, 0.0, 0.0)), DivisionByZero);
    assert_eq!(kind(Vector3::new(f64::NAN, 0.0, 0.0)), NotFinite);
    assert_eq!(kind(Vector3::new(0.0, f64::INFINITY, 0.0)), NotFinite);
    assert_eq!(kind(Vector3::new(0.0, 0.0, 1e200)), NotFinite);
}

// Expected values are ERFA pn outputs, which multiply by 1/|v| instead of dividing.
#[test]
fn test_normalize_matches_erfa_pn() {
    let v = Vector3::new(
        1.439176483123266,
        -0.42279464657465304,
        -0.07764850380203292,
    );
    let expected = Vector3::new(
        0.9581713955307932,
        -0.28148718470734346,
        -0.051696630761648946,
    );
    assert_eq!(v.normalize().unwrap(), expected);

    let v = Vector3::new(0.41912109110457063, 1.4010298681494708, 1.4635196977209102);
    let expected = Vector3::new(0.20257948231024275, 0.6771787710394523, 0.7073828280361426);
    assert_eq!(v.normalize().unwrap(), expected);
}

#[test]
fn test_to_array() {
    let v = Vector3::new(1.5, 2.5, 3.5);
    let arr = v.to_array();
    assert_eq!(arr, [1.5, 2.5, 3.5]);
}

#[test]
fn test_div_assign_operator() {
    let mut v = Vector3::new(10.0, 20.0, 30.0);
    v /= 2.0;
    assert_eq!(v, Vector3::new(5.0, 10.0, 15.0));
}

#[test]
fn test_display_formatting() {
    let v = Vector3::new(1.234567890, -2.345678901, 3.456789012);
    let display_output = format!("{}", v);

    assert!(display_output.contains("Vector3("));
    assert!(display_output.contains("1.234567890"));
    assert!(display_output.contains("-2.345678901"));
    assert!(display_output.contains("3.456789012"));
    assert!(display_output.ends_with(")"));
}

#[test]
fn test_to_spherical_north_pole() {
    let north_pole = Vector3::new(0.0, 0.0, 1.0);
    let (theta, phi) = north_pole.to_spherical();

    assert_eq!(theta, 0.0);
    assert_eq!(phi, HALF_PI);
}

#[test]
fn test_to_spherical_south_pole() {
    let south_pole = Vector3::new(0.0, 0.0, -1.0);
    let (theta, phi) = south_pole.to_spherical();

    assert_eq!(theta, 0.0);
    assert_eq!(phi, -HALF_PI);
}

#[test]
fn test_to_spherical_zero_z() {
    let on_equator = Vector3::new(1.0, 0.0, 0.0);
    let (theta, phi) = on_equator.to_spherical();

    assert_eq!(theta, 0.0);
    assert_eq!(phi, 0.0);
}

#[test]
fn test_to_spherical_zero_vector() {
    let zero = Vector3::zeros();
    let (theta, phi) = zero.to_spherical();

    assert_eq!(theta, 0.0);
    assert_eq!(phi, 0.0);
}

#[test]
fn test_spherical_roundtrip_at_poles() {
    let north_pole = Vector3::new(0.0, 0.0, 1.0);
    let (theta, phi) = north_pole.to_spherical();
    let roundtrip = Vector3::from_spherical(theta, phi);
    assert_eq!(roundtrip, Vector3::new(COS_HALF_PI, 0.0, 1.0));

    let south_pole = Vector3::new(0.0, 0.0, -1.0);
    let (theta, phi) = south_pole.to_spherical();
    let roundtrip = Vector3::from_spherical(theta, phi);
    assert_eq!(roundtrip, Vector3::new(COS_HALF_PI, 0.0, -1.0));
}
