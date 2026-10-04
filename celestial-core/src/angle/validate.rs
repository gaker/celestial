use super::Angle;
use crate::constants::{HALF_PI, PI};
use crate::errors::{AstroError, MathErrorKind};

pub fn validate_right_ascension(angle: Angle) -> Result<Angle, AstroError> {
    angle.normalized()
}

/// Validates declination angle.
///
/// - `beyond_pole = false`: standard range [-90°, +90°]
/// - `beyond_pole = true`: extended range [-180°, +180°] for GEM pier-flipped observations
///
/// The extended range supports the beyond-the-pole convention where German Equatorial
/// Mounts use Dec values from 90° to 180° for pier-flipped observations.
pub fn validate_declination(angle: Angle, beyond_pole: bool) -> Result<Angle, AstroError> {
    let limit = if beyond_pole {
        Limit::BEYOND_POLE
    } else {
        Limit::POLE
    };
    check_within(angle, limit, "validate_declination", "Declination")
}

pub fn validate_latitude(angle: Angle) -> Result<Angle, AstroError> {
    check_within(angle, Limit::POLE, "validate_latitude", "Latitude")
}

struct Limit {
    rad: f64,
    range: &'static str,
}

impl Limit {
    const POLE: Self = Self {
        rad: HALF_PI,
        range: "[-90°, +90°]",
    };
    const BEYOND_POLE: Self = Self {
        rad: PI,
        range: "[-180°, +180°]",
    };
}

fn check_within(angle: Angle, limit: Limit, op: &str, name: &str) -> Result<Angle, AstroError> {
    use MathErrorKind::{NotFinite, OutOfRange};
    let rad = angle.radians();
    if !rad.is_finite() {
        let reason = format!("{name} not finite");
        return Err(AstroError::math_error(op, NotFinite, &reason));
    }
    if (-limit.rad..=limit.rad).contains(&rad) {
        return Ok(angle);
    }
    let degrees = angle.degrees();
    let reason = format!("{name} {degrees:.2}° out of range {}", limit.range);
    Err(AstroError::math_error(op, OutOfRange, &reason))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::TWOPI;

    #[test]
    fn test_validate_right_ascension_valid() {
        let in_range = Angle::from_radians(0.75);
        assert_eq!(validate_right_ascension(in_range).unwrap().radians(), 0.75);

        let negative = Angle::from_radians(-HALF_PI);
        let wrapped = validate_right_ascension(negative).unwrap();
        assert_eq!(wrapped.radians(), TWOPI - HALF_PI);
    }

    #[test]
    fn test_validate_right_ascension_not_finite() {
        let angle = Angle::from_radians(f64::NAN);
        let result = validate_right_ascension(angle);
        assert!(result.is_err());
        if let Err(AstroError::MathError { kind, .. }) = result {
            assert_eq!(kind, MathErrorKind::NotFinite);
        } else {
            panic!("Expected MathError with NotFinite");
        }
    }

    #[test]
    fn test_validate_right_ascension_infinite() {
        let angle = Angle::from_radians(f64::INFINITY);
        let result = validate_right_ascension(angle);
        assert!(result.is_err());
        if let Err(AstroError::MathError { kind, .. }) = result {
            assert_eq!(kind, MathErrorKind::NotFinite);
        } else {
            panic!("Expected MathError with NotFinite");
        }
    }

    fn message(result: Result<Angle, AstroError>) -> String {
        result.unwrap_err().to_string()
    }

    fn accepted(result: Result<Angle, AstroError>) -> f64 {
        result.unwrap().radians()
    }

    #[test]
    fn test_validate_declination_accepts_closed_range() {
        for rad in [-HALF_PI, 0.75, HALF_PI] {
            let angle = Angle::from_radians(rad);
            assert_eq!(accepted(validate_declination(angle, false)), rad);
        }
        for rad in [-PI, 2.0, PI] {
            let angle = Angle::from_radians(rad);
            assert_eq!(accepted(validate_declination(angle, true)), rad);
        }
    }

    #[test]
    fn test_validate_declination_errors() {
        assert_eq!(
            message(validate_declination(Angle::from_degrees(95.0), false)),
            "Math error in validate_declination (out of range): Declination 95.00° out of range [-90°, +90°]"
        );
        assert_eq!(
            message(validate_declination(Angle::from_degrees(-185.0), true)),
            "Math error in validate_declination (out of range): Declination -185.00° out of range [-180°, +180°]"
        );
        assert_eq!(
            message(validate_declination(Angle::from_radians(f64::NAN), false)),
            "Math error in validate_declination (not finite): Declination not finite"
        );
    }

    #[test]
    fn test_validate_latitude_accepts_closed_range() {
        for rad in [-HALF_PI, 0.75, HALF_PI] {
            assert_eq!(accepted(validate_latitude(Angle::from_radians(rad))), rad);
        }
    }

    #[test]
    fn test_validate_latitude_errors_name_latitude() {
        assert_eq!(
            message(validate_latitude(Angle::from_degrees(95.0))),
            "Math error in validate_latitude (out of range): Latitude 95.00° out of range [-90°, +90°]"
        );
        assert_eq!(
            message(validate_latitude(Angle::from_radians(f64::INFINITY))),
            "Math error in validate_latitude (not finite): Latitude not finite"
        );
    }
}
