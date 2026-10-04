use crate::constants::{HALF_PI, PI};
use crate::errors::{AstroError, AstroResult, MathErrorKind};

const LATITUDE_RANGE: &str = "Latitude outside valid range [-90°, +90°]";
const LONGITUDE_RANGE: &str = "Longitude outside valid range [-180°, +180°]";
const HEIGHT_RANGE: &str = "Height outside reasonable range [-12000, 100000] meters";
const HEIGHT_RANGE_M: std::ops::RangeInclusive<f64> = -12000.0..=100000.0;

pub(super) fn validate(latitude: f64, longitude: f64, height: f64) -> AstroResult<()> {
    finite(latitude, "Latitude must be finite")?;
    finite(longitude, "Longitude must be finite")?;
    finite(height, "Height must be finite")?;
    within(libm::fabs(latitude) <= HALF_PI, LATITUDE_RANGE)?;
    within(libm::fabs(longitude) <= PI, LONGITUDE_RANGE)?;
    within(HEIGHT_RANGE_M.contains(&height), HEIGHT_RANGE)
}

fn finite(value: f64, message: &str) -> AstroResult<()> {
    require(value.is_finite(), MathErrorKind::NotFinite, message)
}

fn within(condition: bool, message: &str) -> AstroResult<()> {
    require(condition, MathErrorKind::OutOfRange, message)
}

fn require(condition: bool, kind: MathErrorKind, message: &str) -> AstroResult<()> {
    if condition {
        return Ok(());
    }
    Err(AstroError::math_error("location_validation", kind, message))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::location::Location;

    fn error_kind(result: AstroResult<Location>) -> MathErrorKind {
        match result {
            Err(AstroError::MathError { kind, .. }) => kind,
            other => panic!("expected a math error, got {other:?}"),
        }
    }

    #[test]
    fn test_validation_errors_distinguish_non_finite_from_out_of_range() {
        use MathErrorKind::{NotFinite, OutOfRange};
        assert_eq!(error_kind(Location::new(f64::NAN, 0.0, 0.0)), NotFinite);
        assert_eq!(
            error_kind(Location::new(0.0, 0.0, f64::INFINITY)),
            NotFinite
        );
        assert_eq!(error_kind(Location::new(2.0, 0.0, 0.0)), OutOfRange);
        assert_eq!(error_kind(Location::new(0.0, 0.0, 1.0e6)), OutOfRange);
        assert_eq!(
            error_kind(Location::from_degrees(f64::NAN, 0.0, 0.0)),
            NotFinite
        );
        assert_eq!(
            error_kind(Location::from_degrees(0.0, 360.5, 0.0)),
            OutOfRange
        );
    }

    fn message(result: AstroResult<Location>) -> String {
        result.unwrap_err().to_string()
    }

    const NOT_FINITE: &str = "Math error in location_validation (not finite):";
    const OUT_OF_RANGE: &str = "Math error in location_validation (out of range):";
    const HEIGHT_REASON: &str = "Height outside reasonable range [-12000, 100000] meters";

    #[test]
    fn test_new_validation_messages() {
        use crate::constants::{PI, TWOPI};
        let cases = [
            ((f64::NAN, 0.0, 0.0), NOT_FINITE, "Latitude must be finite"),
            (
                (0.0, f64::INFINITY, 0.0),
                NOT_FINITE,
                "Longitude must be finite",
            ),
            ((0.0, 0.0, f64::NAN), NOT_FINITE, "Height must be finite"),
            (
                (-PI, 0.0, 0.0),
                OUT_OF_RANGE,
                "Latitude outside valid range [-90°, +90°]",
            ),
            (
                (0.0, TWOPI, 0.0),
                OUT_OF_RANGE,
                "Longitude outside valid range [-180°, +180°]",
            ),
            ((0.0, 0.0, -20000.0), OUT_OF_RANGE, HEIGHT_REASON),
            ((0.0, 0.0, 200000.0), OUT_OF_RANGE, HEIGHT_REASON),
        ];
        for ((lat, lon, height), kind, reason) in cases {
            let result = Location::new(lat, lon, height);
            assert_eq!(message(result), format!("{kind} {reason}"));
        }
    }

    #[test]
    fn test_from_degrees_validation_messages() {
        let cases = [
            ((f64::NAN, 0.0), NOT_FINITE, "Latitude must be finite"),
            ((0.0, f64::INFINITY), NOT_FINITE, "Longitude must be finite"),
            (
                (95.0, 0.0),
                OUT_OF_RANGE,
                "Latitude outside valid range [-90°, +90°]",
            ),
            (
                (-95.0, 0.0),
                OUT_OF_RANGE,
                "Latitude outside valid range [-90°, +90°]",
            ),
            (
                (0.0, -180.5),
                OUT_OF_RANGE,
                "Longitude outside valid range [-180°, +180°]",
            ),
            (
                (0.0, 360.5),
                OUT_OF_RANGE,
                "Longitude outside valid range [-180°, +180°]",
            ),
        ];
        for ((lat, lon), kind, reason) in cases {
            let result = Location::from_degrees(lat, lon, 0.0);
            assert_eq!(message(result), format!("{kind} {reason}"));
        }
    }
}
