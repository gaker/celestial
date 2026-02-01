use super::core::Angle;
use crate::constants::{HALF_PI, PI};
use crate::{AstroError, MathErrorKind};

pub fn validate_right_ascension(angle: Angle) -> Result<Angle, AstroError> {
    let rad = angle.radians();
    if rad.is_finite() {
        let normalized = super::normalize::wrap_0_2pi(rad);
        return Ok(Angle::from_radians(normalized));
    }

    Err(AstroError::math_error(
        "validate_right_ascension",
        MathErrorKind::NotFinite,
        "RA Not Finite",
    ))
}

/// Validates declination angle.
///
/// - `beyond_pole = false`: standard range [-90°, +90°]
/// - `beyond_pole = true`: extended range [-180°, +180°] for GEM pier-flipped observations
///
/// The extended range supports the TPOINT beyond-the-pole convention where German Equatorial
/// Mounts use Dec values from 90° to 180° for pier-flipped observations.
pub fn validate_declination(angle: Angle, beyond_pole: bool) -> Result<Angle, AstroError> {
    let rad = angle.radians();
    if !rad.is_finite() {
        return Err(AstroError::math_error(
            "validate_declination",
            MathErrorKind::NotFinite,
            "Dec not Finite",
        ));
    }

    let (limit, range_desc) = if beyond_pole {
        (PI, "[-180°, +180°]")
    } else {
        (HALF_PI, "[-90°, +90°]")
    };

    if (-limit..=limit).contains(&rad) {
        return Ok(angle);
    }

    Err(AstroError::math_error(
        "validate_declination",
        MathErrorKind::OutOfRange,
        &format!("Dec {:.2}° out of range {}", angle.degrees(), range_desc),
    ))
}

pub fn validate_latitude(angle: Angle) -> Result<Angle, AstroError> {
    validate_declination(angle, false)
}

pub fn validate_longitude(angle: Angle, normalize: bool) -> Result<Angle, AstroError> {
    let rad = angle.radians();
    if !rad.is_finite() {
        return Err(AstroError::math_error(
            "validate_longitude",
            MathErrorKind::NotFinite,
            "Lon not finite",
        ));
    }

    if normalize {
        let normalized = super::normalize::wrap_0_2pi(rad);
        return Ok(Angle::from_radians(normalized));
    }

    if (-PI..=PI).contains(&rad) {
        return Ok(angle);
    }

    Err(AstroError::math_error(
        "validate_longitude",
        MathErrorKind::OutOfRange,
        "Lon out of range",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::TWOPI;

    #[test]
    fn test_validate_right_ascension_valid() {
        let angle = Angle::from_degrees(45.0);
        let result = validate_right_ascension(angle);
        assert!(result.is_ok());
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

    #[test]
    fn test_validate_declination() {
        // Valid standard range
        assert!(validate_declination(Angle::from_degrees(45.0), false).is_ok());
        assert!(validate_declination(Angle::from_degrees(-90.0), false).is_ok());

        // Out of standard range
        assert!(validate_declination(Angle::from_degrees(95.0), false).is_err());

        // Valid with beyond_pole
        assert!(validate_declination(Angle::from_degrees(120.0), true).is_ok());

        // Not finite
        assert!(validate_declination(Angle::from_radians(f64::NAN), false).is_err());
    }

    #[test]
    fn test_validate_latitude_delegates_to_declination() {
        let valid = Angle::from_degrees(45.0);
        assert!(validate_latitude(valid).is_ok());

        let invalid = Angle::from_degrees(95.0);
        assert!(validate_latitude(invalid).is_err());
    }

    #[test]
    fn test_validate_longitude_valid() {
        let angle = Angle::from_degrees(45.0);
        let result = validate_longitude(angle, false);
        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_longitude_not_finite() {
        let angle = Angle::from_radians(f64::NAN);
        let result = validate_longitude(angle, false);
        assert!(result.is_err());
        if let Err(AstroError::MathError { kind, .. }) = result {
            assert_eq!(kind, MathErrorKind::NotFinite);
        } else {
            panic!("Expected MathError with NotFinite");
        }
    }

    #[test]
    fn test_validate_longitude_normalized() {
        let angle = Angle::from_degrees(370.0);
        let result = validate_longitude(angle, true);
        assert!(result.is_ok());
        let normalized = result.unwrap();
        assert!(normalized.radians() >= 0.0 && normalized.radians() < TWOPI);
    }

    #[test]
    fn test_validate_longitude_out_of_range() {
        let angle = Angle::from_degrees(190.0);
        let result = validate_longitude(angle, false);
        assert!(result.is_err());
        if let Err(AstroError::MathError { kind, .. }) = result {
            assert_eq!(kind, MathErrorKind::OutOfRange);
        } else {
            panic!("Expected MathError with OutOfRange");
        }
    }
}
