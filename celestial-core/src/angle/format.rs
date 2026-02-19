//! Angle formatting and lightweight parsing for astronomical coordinates.
//!
//! This module provides formatters for displaying angles in astronomical notation
//! and a simple parser for common angle formats. For more flexible parsing with support
//! for verbose formats like "12 hours 30 minutes", see the [`parse`](super::parse) module.
//!
//! # Formatting Conventions
//!
//! Astronomy uses two primary sexagesimal (base-60) notations:
//!
//! ## Degrees-Minutes-Seconds (DMS)
//!
//! Used for declination, latitude, altitude, and general angular measurements.
//! - Format: `+DD° MM' SS.ss"`
//! - Sign is always shown (+ or -)
//! - 1 degree = 60 arcminutes = 3600 arcseconds
//!
//! ## Hours-Minutes-Seconds (HMS)
//!
//! Used for right ascension and hour angles.
//! - Format: `HHʰ MMᵐ SS.ssˢ`
//! - Always positive; negative angles wrap to [0, 24h)
//! - 24 hours = 360 degrees, so 1 hour = 15 degrees
//!
//! # Formatting Examples
//!
//! ```
//! use celestial_core::Angle;
//! use celestial_core::angle::{DmsFmt, HmsFmt};
//!
//! // Declination of Vega: +38° 47' 01"
//! let dec = Angle::from_degrees(38.783611);
//! let dms = DmsFmt { frac_digits: 0 };
//! assert_eq!(dms.fmt(dec), "+38° 47' 1\"");
//!
//! // Right ascension of Vega: 18h 36m 56s
//! let ra = Angle::from_hours(18.615556);
//! let hms = HmsFmt { frac_digits: 0 };
//! assert_eq!(hms.fmt(ra), "18ʰ 36ᵐ 56ˢ");
//!
//! // With fractional seconds
//! let hms_precise = HmsFmt { frac_digits: 2 };
//! assert_eq!(hms_precise.fmt(ra), "18ʰ 36ᵐ 56.00ˢ");
//! ```
//!
//! # Parsing Examples
//!
//! The [`parse_angle`] function handles common formats:
//!
//! ```
//! use celestial_core::angle::parse_angle;
//!
//! // HMS formats (tries HMS first, then DMS)
//! let ra = parse_angle("12h30m15s").unwrap();
//! assert!((ra.angle.hours() - 12.504166666666666).abs() < 1e-10);
//!
//! // Also accepts Unicode superscript notation
//! let ra2 = parse_angle("12ʰ30ᵐ15ˢ").unwrap();
//!
//! // Colon-separated (interpreted as HMS by default)
//! let ra3 = parse_angle("12:30:15").unwrap();
//!
//! // DMS formats
//! let dec = parse_angle("45°30'15\"").unwrap();
//! assert!((dec.angle.degrees() - 45.504166666666666).abs() < 1e-10);
//! ```
//!
//! # Default Display
//!
//! The `Display` trait formats angles as decimal degrees with 6 decimal places:
//!
//! ```
//! use celestial_core::Angle;
//!
//! let a = Angle::from_degrees(45.123456789);
//! assert_eq!(format!("{}", a), "45.123457°");
//! ```
use super::Angle;
use core::fmt;

/// Formatter for degrees-minutes-seconds (DMS) notation.
///
/// DMS is the standard format for declination, latitude, altitude, and general
/// angular measurements in astronomy. The sign is always explicit.
///
/// # Fields
///
/// * `frac_digits` - Number of decimal places for the arcseconds component.
///   Use 0 for whole arcseconds, 2-3 for sub-arcsecond precision.
///
/// # Output Format
///
/// `±DD° MM' SS.ss"` where:
/// - Sign is always shown (+ or -)
/// - Degrees, arcminutes are whole numbers
/// - Arcseconds include decimals per `frac_digits`
///
/// # Example
///
/// ```
/// use celestial_core::Angle;
/// use celestial_core::angle::DmsFmt;
///
/// let dec = Angle::from_degrees(-23.4392);
///
/// // Whole arcseconds
/// let fmt0 = DmsFmt { frac_digits: 0 };
/// assert_eq!(fmt0.fmt(dec), "-23° 26' 21\"");
///
/// // Sub-arcsecond precision (typical for catalogs)
/// let fmt2 = DmsFmt { frac_digits: 2 };
/// assert_eq!(fmt2.fmt(dec), "-23° 26' 21.12\"");
/// ```
pub struct DmsFmt {
    pub frac_digits: u8,
}

/// Formatter for hours-minutes-seconds (HMS) notation.
///
/// HMS is the standard format for right ascension and hour angles in astronomy.
/// The output is always positive; negative angles are wrapped to [0, 24h).
///
/// # Fields
///
/// * `frac_digits` - Number of decimal places for the seconds component.
///   Use 0 for whole seconds, 2-3 for sub-second precision.
///
/// # Output Format
///
/// `HHʰ MMᵐ SS.ssˢ` where:
/// - Uses Unicode superscript characters (ʰ, ᵐ, ˢ)
/// - Hours, minutes are whole numbers
/// - Seconds include decimals per `frac_digits`
/// - Negative angles wrap: -1.5h becomes 22h 30m
///
/// # Example
///
/// ```
/// use celestial_core::Angle;
/// use celestial_core::angle::HmsFmt;
///
/// let ra = Angle::from_hours(14.5);  // 14h 30m 00s
///
/// let fmt = HmsFmt { frac_digits: 1 };
/// assert_eq!(fmt.fmt(ra), "14ʰ 30ᵐ 0.0ˢ");
///
/// // Negative angles wrap to positive
/// let neg = Angle::from_hours(-1.5);
/// assert_eq!(fmt.fmt(neg), "22ʰ 30ᵐ 0.0ˢ");
/// ```
pub struct HmsFmt {
    pub frac_digits: u8,
}

impl DmsFmt {
    /// Formats an angle as degrees-minutes-seconds.
    ///
    /// Decomposes the angle into integer degrees and arcminutes, with arcseconds
    /// shown to the precision specified by `frac_digits`.
    ///
    /// # Arguments
    ///
    /// * `a` - The angle to format
    ///
    /// # Returns
    ///
    /// A string in the format `±DD° MM' SS.ss"`.
    #[inline]
    pub fn fmt(&self, a: Angle) -> String {
        let sign = if a.degrees() < 0.0 { '-' } else { '+' };
        let mut d = a.degrees().abs();
        let deg = libm::trunc(d);
        d = (d - deg) * 60.0;
        let min = libm::trunc(d);
        let sec = (d - min) * 60.0;
        format!(
            "{sign}{deg:.0}° {min:.0}' {sec:.*}\"",
            self.frac_digits as usize
        )
    }
}

impl HmsFmt {
    /// Formats an angle as hours-minutes-seconds.
    ///
    /// Decomposes the angle into integer hours and minutes, with seconds
    /// shown to the precision specified by `frac_digits`. Negative angles
    /// are wrapped to the range [0, 24h) using Euclidean remainder.
    ///
    /// # Arguments
    ///
    /// * `a` - The angle to format
    ///
    /// # Returns
    ///
    /// A string in the format `HHʰ MMᵐ SS.ssˢ` using Unicode superscript markers.
    #[inline]
    pub fn fmt(&self, a: Angle) -> String {
        let mut h = a.hours();
        h = h.rem_euclid(24.0);
        let hh = libm::trunc(h);
        h = (h - hh) * 60.0;
        let mm = libm::trunc(h);
        let ss = (h - mm) * 60.0;
        format!("{hh:.0}ʰ {mm:.0}ᵐ {ss:.*}ˢ", self.frac_digits as usize)
    }
}

impl fmt::Display for Angle {
    /// Formats the angle as decimal degrees with 6 decimal places.
    ///
    /// This provides a simple, unambiguous representation suitable for debugging
    /// and data export. For astronomical notation, use [`DmsFmt`] or [`HmsFmt`].
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:.6}°", self.degrees())
    }
}

/// Result of parsing an angle string.
///
/// This struct wraps the parsed [`Angle`] and may be extended in the future
/// to include metadata about the parse (detected unit, original sign, etc.).
///
/// # Example
///
/// ```
/// use celestial_core::angle::parse_angle;
///
/// let parsed = parse_angle("12h30m15s").unwrap();
/// let angle = parsed.angle;  // Extract the Angle
/// ```
pub struct ParsedAngle {
    /// The parsed angle value.
    pub angle: Angle,
}

/// Parses an angle string, trying HMS format first, then DMS.
///
/// This is a lightweight parser for angle formats.
/// For more flexible parsing including verbose formats ("12 hours 30 minutes"),
/// see [`super::parse::parse_hms`] and [`super::parse::parse_dms`].
///
/// # Supported Formats
///
/// **HMS (hours-minutes-seconds):**
/// - `12h30m15s` or `12h30m15.5s`
/// - `12ʰ30ᵐ15ˢ` (Unicode superscripts)
/// - `12:30:15` (colon-separated)
/// - `12h` (hours only)
/// - `-12h30m15s` (negative)
///
/// **DMS (degrees-minutes-seconds):**
/// - `45°30'15"` or `45°30'15.5"`
/// - `45d30m15s`
/// - `45:30:15` (colon-separated, tried if HMS fails)
/// - `45°` (degrees only)
/// - `-45°30'15"` (negative)
///
/// # Ambiguity
///
/// Colon-separated values like `12:30:15` are tried as HMS first. If you need
/// to parse this as DMS explicitly, use [`parse_dms`](super::parse::parse_dms).
///
/// # Errors
///
/// Returns [`AstroError`](crate::AstroError) if:
/// - The string is empty or contains no valid components
/// - Minutes or seconds are outside [0, 60)
/// - Fractional hours/degrees are mixed with minutes/seconds (e.g., "12.5h30m")
/// - The string cannot be parsed as either HMS or DMS
///
/// # Example
///
/// ```
/// use celestial_core::angle::parse_angle;
///
/// // Right ascension
/// let ra = parse_angle("05h14m32.27s").unwrap();
/// assert!((ra.angle.hours() - 5.242297).abs() < 1e-5);
///
/// // Declination
/// let dec = parse_angle("-08°12'05.9\"").unwrap();
/// assert!((dec.angle.degrees() - (-8.201639)).abs() < 1e-5);
/// ```
pub fn parse_angle(s: &str) -> Result<ParsedAngle, crate::AstroError> {
    parse_hms(s).or_else(|_| parse_dms(s))
}

/// Parses an HMS (hours-minutes-seconds) string into an angle.
///
/// Accepts formats like: `12h30m15s`, `12ʰ30ᵐ15ˢ`, `12:30:15`, `12h`, `-12h30m15s`
fn parse_hms(s: &str) -> Result<ParsedAngle, crate::AstroError> {
    let s = s.trim();
    let sign = if s.starts_with('-') { -1.0 } else { 1.0 };
    let s = s.trim_start_matches(['+', '-']);

    let parts: Vec<&str> = s
        .split(['h', 'ʰ', 'm', 'ᵐ', 's', 'ˢ', ':'])
        .map(|p| p.trim())
        .filter(|p| !p.is_empty())
        .collect();

    if parts.is_empty() {
        return Err(crate::AstroError::math_error(
            "parse_hms",
            crate::errors::MathErrorKind::InvalidInput,
            "Empty string",
        ));
    }

    if parts.len() > 3 {
        return Err(crate::AstroError::math_error(
            "parse_hms",
            crate::errors::MathErrorKind::InvalidInput,
            "Too many components (max 3: hours, minutes, seconds)",
        ));
    }

    let h = parts[0].parse::<f64>().map_err(|_| {
        crate::AstroError::math_error(
            "parse_hms",
            crate::errors::MathErrorKind::InvalidInput,
            "Invalid hours",
        )
    })?;

    let m = if parts.len() > 1 {
        parts[1].parse::<f64>().map_err(|_| {
            crate::AstroError::math_error(
                "parse_hms",
                crate::errors::MathErrorKind::InvalidInput,
                "Invalid minutes",
            )
        })?
    } else {
        0.0
    };

    let sec = if parts.len() > 2 {
        parts[2].parse::<f64>().map_err(|_| {
            crate::AstroError::math_error(
                "parse_hms",
                crate::errors::MathErrorKind::InvalidInput,
                "Invalid seconds",
            )
        })?
    } else {
        0.0
    };

    if parts.len() > 1 && h - libm::trunc(h) != 0.0 {
        return Err(crate::AstroError::math_error(
            "parse_hms",
            crate::errors::MathErrorKind::InvalidInput,
            "Cannot mix fractional hours with minutes/seconds",
        ));
    }

    if !(0.0..60.0).contains(&m) {
        return Err(crate::AstroError::math_error(
            "parse_hms",
            crate::errors::MathErrorKind::InvalidInput,
            "Minutes must be in range [0, 60)",
        ));
    }

    if !(0.0..60.0).contains(&sec) {
        return Err(crate::AstroError::math_error(
            "parse_hms",
            crate::errors::MathErrorKind::InvalidInput,
            "Seconds must be in range [0, 60)",
        ));
    }

    Ok(ParsedAngle {
        angle: Angle::from_hours(sign * (h.abs() + m / 60.0 + sec / 3600.0)),
    })
}

/// Parses a DMS (degrees-minutes-seconds) string into an angle.
///
/// Accepts formats like: `45°30'15"`, `45d30m15s`, `45:30:15`, `45°`, `-45°30'15"`
fn parse_dms(s: &str) -> Result<ParsedAngle, crate::AstroError> {
    let s = s.trim();
    let sign = if s.starts_with('-') { -1.0 } else { 1.0 };
    let s = s.trim_start_matches(['+', '-']);

    let parts: Vec<&str> = s
        .split(['°', '\'', '"', ':', 'd', 'm', 's'])
        .map(|p| p.trim())
        .filter(|p| !p.is_empty())
        .collect();

    if parts.is_empty() {
        return Err(crate::AstroError::math_error(
            "parse_dms",
            crate::errors::MathErrorKind::InvalidInput,
            "Empty string",
        ));
    }

    if parts.len() > 3 {
        return Err(crate::AstroError::math_error(
            "parse_dms",
            crate::errors::MathErrorKind::InvalidInput,
            "Too many components (max 3: degrees, arcminutes, arcseconds)",
        ));
    }

    let deg = parts[0].parse::<f64>().map_err(|_| {
        crate::AstroError::math_error(
            "parse_dms",
            crate::errors::MathErrorKind::InvalidInput,
            "Invalid degrees",
        )
    })?;

    let min = if parts.len() > 1 {
        parts[1].parse::<f64>().map_err(|_| {
            crate::AstroError::math_error(
                "parse_dms",
                crate::errors::MathErrorKind::InvalidInput,
                "Invalid arcminutes",
            )
        })?
    } else {
        0.0
    };

    let sec = if parts.len() > 2 {
        parts[2].parse::<f64>().map_err(|_| {
            crate::AstroError::math_error(
                "parse_dms",
                crate::errors::MathErrorKind::InvalidInput,
                "Invalid arcseconds",
            )
        })?
    } else {
        0.0
    };

    if parts.len() > 1 && deg - libm::trunc(deg) != 0.0 {
        return Err(crate::AstroError::math_error(
            "parse_dms",
            crate::errors::MathErrorKind::InvalidInput,
            "Cannot mix fractional degrees with arcminutes/arcseconds",
        ));
    }

    if !(0.0..60.0).contains(&min) {
        return Err(crate::AstroError::math_error(
            "parse_dms",
            crate::errors::MathErrorKind::InvalidInput,
            "Arcminutes must be in range [0, 60)",
        ));
    }

    if !(0.0..60.0).contains(&sec) {
        return Err(crate::AstroError::math_error(
            "parse_dms",
            crate::errors::MathErrorKind::InvalidInput,
            "Arcseconds must be in range [0, 60)",
        ));
    }

    Ok(ParsedAngle {
        angle: Angle::from_degrees(sign * (deg.abs() + min / 60.0 + sec / 3600.0)),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hms_format_normal() {
        let a = Angle::from_hours(12.5);
        let fmt = HmsFmt { frac_digits: 2 };
        let result = fmt.fmt(a);
        assert!(result.contains("12ʰ"));
        assert!(result.contains("30ᵐ"));
    }

    #[test]
    fn test_hms_format_extreme_positive() {
        let a = Angle::from_degrees(720.0);
        let fmt = HmsFmt { frac_digits: 0 };
        let result = fmt.fmt(a);
        assert!(result.contains("0ʰ"));
    }

    #[test]
    fn test_hms_format_extreme_negative() {
        let a = Angle::from_degrees(-750.0);
        let fmt = HmsFmt { frac_digits: 0 };
        let result = fmt.fmt(a);
        assert!(result.contains("22ʰ"));
    }

    #[test]
    fn test_dms_format_negative_with_precision() {
        let a = Angle::from_degrees(-12.345678);
        let fmt = DmsFmt { frac_digits: 2 };
        let result = fmt.fmt(a);
        assert_eq!(result, "-12° 20' 44.44\"");
    }

    #[test]
    fn test_hms_format_wraps_negative_angle() {
        let a = Angle::from_hours(-1.5);
        let fmt = HmsFmt { frac_digits: 1 };
        let result = fmt.fmt(a);
        assert_eq!(result, "22ʰ 30ᵐ 0.0ˢ");
    }

    #[test]
    fn test_angle_display_precision() {
        let a = Angle::from_degrees(1.23456789);
        assert_eq!(format!("{a}"), "1.234568°");
    }

    #[test]
    fn test_parse_hms() {
        let result = parse_hms("12h30m15s").unwrap();
        assert!((result.angle.hours() - 12.504166666666666).abs() < 1e-10);
    }

    #[test]
    fn test_parse_hms_unicode() {
        let result = parse_hms("12ʰ30ᵐ15ˢ").unwrap();
        assert!((result.angle.hours() - 12.504166666666666).abs() < 1e-10);
    }

    #[test]
    fn test_parse_hms_colon() {
        let result = parse_hms("12:30:15").unwrap();
        assert!((result.angle.hours() - 12.504166666666666).abs() < 1e-10);
    }

    #[test]
    fn test_parse_hms_partial() {
        let result = parse_hms("12h").unwrap();
        assert_eq!(result.angle.hours(), 12.0);
    }

    #[test]
    fn test_parse_dms_positive() {
        let result = parse_dms("45°30'15\"").unwrap();
        assert!((result.angle.degrees() - 45.50416666666667).abs() < 1e-10);
    }

    #[test]
    fn test_parse_dms_negative() {
        let result = parse_dms("-45°30'15\"").unwrap();
        assert!((result.angle.degrees() + 45.50416666666667).abs() < 1e-10);
    }

    #[test]
    fn test_parse_dms_colon() {
        let result = parse_dms("45:30:15").unwrap();
        assert!((result.angle.degrees() - 45.50416666666667).abs() < 1e-10);
    }

    #[test]
    fn test_parse_dms_partial() {
        let result = parse_dms("45°").unwrap();
        assert_eq!(result.angle.degrees(), 45.0);
    }

    #[test]
    fn test_parse_angle_dispatch_hms() {
        let result = parse_angle("12h30m15s").unwrap();
        assert!((result.angle.hours() - 12.504166666666666).abs() < 1e-10);
    }

    #[test]
    fn test_parse_angle_dispatch_dms() {
        let result = parse_angle("45°30'15\"").unwrap();
        assert!((result.angle.degrees() - 45.50416666666667).abs() < 1e-10);
    }

    #[test]
    fn test_parse_hms_negative() {
        let result = parse_hms("-01:30:00").unwrap();
        assert_eq!(result.angle.hours(), -1.5);
    }

    #[test]
    fn test_parse_hms_negative_with_seconds() {
        let result = parse_hms("-12h30m45s").unwrap();
        assert_eq!(result.angle.hours(), -12.5125);
    }

    #[test]
    fn test_parse_hms_invalid_minutes() {
        let result = parse_hms("12h99m00s");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_hms_invalid_seconds() {
        let result = parse_hms("12h30m80s");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_dms_invalid_arcminutes() {
        let result = parse_dms("45°80'00\"");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_dms_invalid_arcseconds() {
        let result = parse_dms("45°30'99\"");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_dms_negative_with_minutes_seconds() {
        let result = parse_dms("-45°30'15\"").unwrap();
        assert_eq!(result.angle.degrees(), -45.50416666666667);
    }

    #[test]
    fn test_parse_hms_rejects_fractional_hours_with_minutes() {
        let result = parse_hms("12.5h30m");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_hms_rejects_empty_string() {
        let result = parse_hms("");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_hms_rejects_too_many_components() {
        let result = parse_hms("12:30:15:99");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_hms_accepts_fractional_hours_alone() {
        let result = parse_hms("12.5h").unwrap();
        assert_eq!(result.angle.hours(), 12.5);
    }

    #[test]
    fn test_parse_dms_rejects_fractional_degrees_with_arcminutes() {
        let result = parse_dms("45.5°30'");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_dms_rejects_empty_string() {
        let result = parse_dms("   ");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_dms_rejects_too_many_components() {
        let result = parse_dms("45:30:15:99");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_dms_accepts_fractional_degrees_alone() {
        let result = parse_dms("45.5°").unwrap();
        assert_eq!(result.angle.degrees(), 45.5);
    }

    #[test]
    fn test_parse_dms_invalid_degrees() {
        let result = parse_dms("abc°");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_angle_fails_for_unknown_format() {
        let result = parse_angle("not an angle");
        assert!(result.is_err());
    }
}
