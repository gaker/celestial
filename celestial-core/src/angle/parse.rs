//! Angle parsing from string representations.
//!
//! This module provides flexible parsing for angles in formats commonly used in astronomy:
//!
//! - **HMS (Hours-Minutes-Seconds)**: Used for Right Ascension. 1 hour = 15 degrees.
//! - **DMS (Degrees-Minutes-Seconds)**: Used for Declination, altitude, and general angles.
//! - **Decimal**: Plain numeric values with explicit unit conversion.
//!
//! # Format Support
//!
//! Both HMS and DMS accept multiple notations:
//!
//! ```text
//! Colon-separated:  12:34:56.789
//! Letter markers:   12h34m56.789s  or  45d30m15s
//! Verbose:          12 hours 34 minutes 56 seconds
//! Symbol notation:  45d 30' 15"  or  45d 30' 15''
//! ```
//!
//! Signs are only valid at the beginning: `-12:34:56` works, `12:-34:56` does not.
//!
//! # Usage Patterns
//!
//! Two traits provide parsing:
//!
//! - [`AngleUnits`]: Explicit unit conversion via methods like `.deg()`, `.hms()`
//! - [`ParseAngle`]: Auto-detection via `.to_angle()` (tries HMS, then DMS, then decimal degrees)
//!
//! ```
//! use celestial_core::angle::{AngleUnits, ParseAngle};
//!
//! // Explicit unit - you know what format you have
//! let ra = "12:34:56".hms().unwrap();      // Right ascension
//! let dec = "-45:30:15".dms().unwrap();    // Declination
//! let alt = "30.5".deg().unwrap();         // Altitude in degrees
//!
//! // Auto-detection - useful for user input
//! let angle = "12:34:56".to_angle().unwrap();  // Parsed as HMS (tries first)
//! let angle = "45d30m15s".to_angle().unwrap(); // Parsed as DMS
//! let angle = "45.5".to_angle().unwrap();      // Parsed as decimal degrees
//! ```
//!
//! # HMS vs DMS Ambiguity
//!
//! The colon format `12:34:56` is ambiguous - it could be HMS or DMS. When using
//! auto-detection (`.to_angle()`), HMS is tried first. If you know the intended
//! interpretation, use `.hms()` or `.dms()` explicitly.
//!
//! For Right Ascension, always use `.hms()`. For Declination, always use `.dms()`.

use super::Angle;
use crate::AstroError;
use once_cell::sync::Lazy;
use regex::Regex;

/// Parse strings as angles with explicit unit specification.
///
/// Implemented for `str`. Each method interprets the string value in its respective unit.
///
/// # Decimal Methods
///
/// - `deg()` - Parse as decimal degrees
/// - `rad()` - Parse as radians
/// - `hours()` - Parse as decimal hours (1h = 15 deg)
/// - `arcmin()` - Parse as arcminutes (60' = 1 deg)
/// - `arcsec()` - Parse as arcseconds (3600" = 1 deg)
///
/// # Sexagesimal Methods
///
/// - `hms()` - Parse hours:minutes:seconds format
/// - `dms()` - Parse degrees:minutes:seconds format
pub trait AngleUnits {
    /// Parse as decimal degrees.
    fn deg(&self) -> Result<Angle, AstroError>;
    /// Parse as radians.
    fn rad(&self) -> Result<Angle, AstroError>;
    /// Parse as decimal hours (1 hour = 15 degrees).
    fn hours(&self) -> Result<Angle, AstroError>;
    /// Parse as arcminutes (60 arcmin = 1 degree).
    fn arcmin(&self) -> Result<Angle, AstroError>;
    /// Parse as arcseconds (3600 arcsec = 1 degree).
    fn arcsec(&self) -> Result<Angle, AstroError>;
    /// Parse degrees-minutes-seconds format. See module docs for accepted formats.
    fn dms(&self) -> Result<Angle, AstroError>;
    /// Parse hours-minutes-seconds format. See module docs for accepted formats.
    fn hms(&self) -> Result<Angle, AstroError>;
}

impl AngleUnits for str {
    #[inline]
    fn deg(&self) -> Result<Angle, AstroError> {
        parse_decimal(self).map(Angle::from_degrees)
    }

    #[inline]
    fn rad(&self) -> Result<Angle, AstroError> {
        parse_decimal(self).map(Angle::from_radians)
    }

    #[inline]
    fn hours(&self) -> Result<Angle, AstroError> {
        parse_decimal(self).map(Angle::from_hours)
    }

    #[inline]
    fn arcmin(&self) -> Result<Angle, AstroError> {
        parse_decimal(self).map(|v| Angle::from_degrees(v / 60.0))
    }

    #[inline]
    fn arcsec(&self) -> Result<Angle, AstroError> {
        parse_decimal(self).map(|v| Angle::from_degrees(v / 3600.0))
    }

    #[inline]
    fn dms(&self) -> Result<Angle, AstroError> {
        parse_dms(self)
    }

    #[inline]
    fn hms(&self) -> Result<Angle, AstroError> {
        parse_hms(self)
    }
}

/// Auto-detect and parse angle format.
///
/// Tries formats in order: HMS, then DMS, then decimal degrees.
/// Use this for user input where format is unknown.
///
/// For coordinates with known semantics (RA vs Dec), prefer explicit
/// `.hms()` or `.dms()` via [`AngleUnits`].
pub trait ParseAngle {
    /// Parse angle, auto-detecting format.
    ///
    /// Detection order: HMS -> DMS -> decimal degrees.
    fn to_angle(&self) -> Result<Angle, AstroError>;
}

impl ParseAngle for str {
    fn to_angle(&self) -> Result<Angle, AstroError> {
        parse_hms(self)
            .or_else(|_| parse_dms(self))
            .or_else(|_| parse_decimal(self).map(Angle::from_degrees))
    }
}

fn parse_decimal(s: &str) -> Result<f64, AstroError> {
    s.trim().parse::<f64>().map_err(|_| {
        AstroError::calculation_error("parse_decimal", &format!("Cannot parse '{}' as number", s))
    })
}

static HMS_REGEX: Lazy<Regex> = Lazy::new(|| {
    Regex::new(
        r#"(?xi)
        ^\s*
        ([+-])?                          # optional sign
        (\d{1,3})                        # hours (1-3 digits)
        (?:                              # separator group
            [:hH\s]+|                    # colons, h/H, spaces
            h(?:ou)?r?s?\s*              # hour/hours variants
        )
        (\d{1,2})                        # minutes (1-2 digits)
        (?:                              # separator group
            [:mM\s']+|                   # colons, m/M, spaces, apostrophes
            m(?:in(?:ute)?s?)?\s*        # min/minute variants
        )
        (\d{1,2}(?:\.\d+)?)              # seconds with optional decimal
        (?:                              # optional trailing markers
            [sS\s"']+|                   # s/S, spaces, quotes
            s(?:ec(?:ond)?s?)?           # sec/second variants
        )?
        \s*$
        "#,
    )
    .unwrap()
});

static DMS_REGEX: Lazy<Regex> = Lazy::new(|| {
    Regex::new(
        r#"(?xi)
        ^\s*
        ([+-])?                          # optional sign
        (\d{1,3})                        # degrees (1-3 digits)
        (?:                              # separator group
            [dD\s:*]+|                   # d/D, colon, asterisk, spaces
            d(?:eg(?:ree)?s?)?\s*        # deg/degree variants
        )
        (\d{1,2})                        # minutes (1-2 digits)
        (?:                              # separator group
            ['mM\s:]+|                   # apostrophes, m/M, spaces, colon
            m(?:in(?:ute)?s?)?\s*|       # min/minute variants
            arc\s?m(?:in(?:ute)?s?)?\s*  # arcmin/arcminute
        )
        (\d{1,2}(?:\.\d+)?)              # seconds with optional decimal
        (?:                              # optional trailing markers
            ["'sS\s]+|                   # quotes, s/S, spaces
            s(?:ec(?:ond)?s?)?|          # sec/second variants
            arc\s?s(?:ec(?:ond)?s?)?     # arcsec/arcsecond
        )?
        \s*$
        "#,
    )
    .unwrap()
});

static COLON_REGEX: Lazy<Regex> =
    Lazy::new(|| Regex::new(r#"^\s*([+-])?(\d{1,4}):(\d{1,3}):(\d{1,3}(?:\.\d+)?)\s*$"#).unwrap());

/// Parse a string as hours-minutes-seconds.
///
/// Accepts formats like `12:34:56`, `12h34m56s`, `12 hours 34 min 56 sec`.
/// Returns the angle with the value interpreted as hours (1h = 15 degrees).
///
/// Use this for Right Ascension values. The result can exceed 24h if the input does.
pub fn parse_hms(s: &str) -> Result<Angle, AstroError> {
    let s = normalize_input(s);

    if let Some(caps) = COLON_REGEX.captures(&s) {
        return parse_hms_captures(caps, &s);
    }

    if let Some(caps) = HMS_REGEX.captures(&s) {
        return parse_hms_captures(caps, &s);
    }

    Err(AstroError::calculation_error(
        "parse_hms",
        &format!("Cannot parse '{}' as HMS format", s),
    ))
}

/// Parse a string as degrees-minutes-seconds.
///
/// Accepts formats like `45:30:15`, `45d30m15s`, `45* 30' 15"`, `45 deg 30 arcmin 15 arcsec`.
/// Returns the angle with the value interpreted as degrees.
///
/// Use this for Declination, altitude, azimuth, or general angular measurements.
pub fn parse_dms(s: &str) -> Result<Angle, AstroError> {
    let s = normalize_input(s);

    if let Some(caps) = COLON_REGEX.captures(&s) {
        return parse_dms_captures(caps, &s);
    }

    if let Some(caps) = DMS_REGEX.captures(&s) {
        return parse_dms_captures(caps, &s);
    }

    Err(AstroError::calculation_error(
        "parse_dms",
        &format!("Cannot parse '{}' as DMS format", s),
    ))
}

fn parse_hms_captures(caps: regex::Captures, _original: &str) -> Result<Angle, AstroError> {
    let sign = caps
        .get(1)
        .map_or(1.0, |m| if m.as_str() == "-" { -1.0 } else { 1.0 });
    let hours: f64 = caps[2].parse().unwrap();
    let minutes: f64 = caps[3].parse().unwrap();
    let seconds: f64 = caps[4].parse().unwrap();

    let total_hours = sign * (hours + minutes / 60.0 + seconds / 3600.0);
    Ok(Angle::from_hours(total_hours))
}

fn parse_dms_captures(caps: regex::Captures, _original: &str) -> Result<Angle, AstroError> {
    let sign = caps
        .get(1)
        .map_or(1.0, |m| if m.as_str() == "-" { -1.0 } else { 1.0 });
    let degrees: f64 = caps[2].parse().unwrap();
    let minutes: f64 = caps[3].parse().unwrap();
    let seconds: f64 = caps[4].parse().unwrap();

    let total_degrees = sign * (degrees + minutes / 60.0 + seconds / 3600.0);
    Ok(Angle::from_degrees(total_degrees))
}

fn normalize_input(s: &str) -> String {
    let mut result = s.trim().to_string();

    result = result.replace("degrees", "d");
    result = result.replace("degree", "d");
    result = result.replace("deg", "d");
    result = result.replace('*', "d");

    result = result.replace("arcminutes", "m");
    result = result.replace("arcminute", "m");
    result = result.replace("arcmin", "m");
    result = result.replace("minutes", "m");
    result = result.replace("minute", "m");
    result = result.replace("min", "m");

    result = result.replace("arcseconds", "s");
    result = result.replace("arcsecond", "s");
    result = result.replace("arcsec", "s");
    result = result.replace("seconds", "s");
    result = result.replace("second", "s");
    result = result.replace("sec", "s");
    result = result.replace("''", "\"");

    result = result.replace("hours", "h");
    result = result.replace("hour", "h");
    result = result.replace("hrs", "h");
    result = result.replace("hr", "h");

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::PI;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_decimal_parsing() {
        assert_eq!("45.5".deg().unwrap().degrees(), 45.5);
        assert_eq!(format!("{:?}", PI).rad().unwrap().radians(), PI);
        assert_eq!("12.5".hours().unwrap().hours(), 12.5);

        assert!(("60.0".arcmin().unwrap().degrees() - 1.0).abs() < EPSILON);
        assert!(("3600.0".arcsec().unwrap().degrees() - 1.0).abs() < EPSILON);

        assert_eq!("-45.5".deg().unwrap().degrees(), -45.5);
        assert_eq!("-12.5".hours().unwrap().hours(), -12.5);

        assert_eq!("  45.5  ".deg().unwrap().degrees(), 45.5);
        assert_eq!(format!("\t{:?}\n", PI).rad().unwrap().radians(), PI);
    }

    #[test]
    fn test_hms_colon_format() {
        let angle = "12:34:56".hms().unwrap();
        let expected_hours = 12.0 + 34.0 / 60.0 + 56.0 / 3600.0;
        assert!((angle.hours() - expected_hours).abs() < EPSILON);

        let angle = "12:34:56.789".hms().unwrap();
        let expected = 12.0 + 34.0 / 60.0 + 56.789 / 3600.0;
        assert!((angle.hours() - expected).abs() < EPSILON);

        let angle = "-5:30:45".hms().unwrap();
        let expected = -(5.0 + 30.0 / 60.0 + 45.0 / 3600.0);
        assert!((angle.hours() - expected).abs() < EPSILON);

        assert!("0:0:0".hms().unwrap().hours() < EPSILON);
        assert!("23:59:59.999".hms().is_ok());
    }

    #[test]
    fn test_dms_colon_format() {
        let angle = "45:30:15".dms().unwrap();
        let expected_deg = 45.0 + 30.0 / 60.0 + 15.0 / 3600.0;
        assert!((angle.degrees() - expected_deg).abs() < EPSILON);

        let angle = "+45:30:15.5".dms().unwrap();
        let expected = 45.0 + 30.0 / 60.0 + 15.5 / 3600.0;
        assert!((angle.degrees() - expected).abs() < EPSILON);

        let angle = "-90:30:0".dms().unwrap();
        let expected = -(90.0 + 30.0 / 60.0);
        assert!((angle.degrees() - expected).abs() < EPSILON);
    }

    #[test]
    fn test_hms_verbose_formats() {
        let angle = "12h34m56s".hms().unwrap();
        let expected = 12.0 + 34.0 / 60.0 + 56.0 / 3600.0;
        assert!((angle.hours() - expected).abs() < EPSILON);

        let angle = "12H34M56S".hms().unwrap();
        assert!((angle.hours() - expected).abs() < EPSILON);

        let angle = "12h 34m 56s".hms().unwrap();
        assert!((angle.hours() - expected).abs() < EPSILON);

        let angle = "12 hours 34 minutes 56 seconds".hms().unwrap();
        assert!((angle.hours() - expected).abs() < EPSILON);

        let angle = "12hr 34min 56sec".hms().unwrap();
        assert!((angle.hours() - expected).abs() < EPSILON);
    }

    #[test]
    fn test_dms_verbose_formats() {
        let angle = "45d30m15s".dms().unwrap();
        let expected = 45.0 + 30.0 / 60.0 + 15.0 / 3600.0;
        assert!((angle.degrees() - expected).abs() < EPSILON);

        let angle = "45*30m15s".dms().unwrap();
        assert!((angle.degrees() - expected).abs() < EPSILON);

        let angle = "45 degrees 30 minutes 15 seconds".dms().unwrap();
        assert!((angle.degrees() - expected).abs() < EPSILON);

        let angle = "45deg 30min 15sec".dms().unwrap();
        assert!((angle.degrees() - expected).abs() < EPSILON);

        let angle = "45d 30 arcmin 15 arcsec".dms().unwrap();
        assert!((angle.degrees() - expected).abs() < EPSILON);
    }

    #[test]
    fn test_quote_formats() {
        let angle = "45d 30' 15\"".dms().unwrap();
        let expected = 45.0 + 30.0 / 60.0 + 15.0 / 3600.0;
        assert!((angle.degrees() - expected).abs() < EPSILON);

        let angle = "45d 30' 15''".dms().unwrap();
        assert!((angle.degrees() - expected).abs() < EPSILON);
    }

    #[test]
    fn test_auto_detection() {
        let angle = "12:34:56".to_angle().unwrap();
        let expected_hours = 12.0 + 34.0 / 60.0 + 56.0 / 3600.0;
        assert!((angle.hours() - expected_hours).abs() < EPSILON);

        let angle = "45d30m15s".to_angle().unwrap();
        let expected_deg = 45.0 + 30.0 / 60.0 + 15.0 / 3600.0;
        assert!((angle.degrees() - expected_deg).abs() < EPSILON);

        let angle = "45.5".to_angle().unwrap();
        assert_eq!(angle.degrees(), 45.5);
    }

    #[test]
    fn test_edge_cases() {
        assert!("0:0:0".hms().unwrap().radians().abs() < EPSILON);
        assert!("0:0:0".dms().unwrap().radians().abs() < EPSILON);
        assert!("0".deg().unwrap().radians().abs() < EPSILON);

        assert!("359:59:59".dms().is_ok());
        assert!("999:59:59".hms().is_ok());

        let angle = "12:34:56.123456789".hms().unwrap();
        assert!(angle.hours() > 12.0);

        assert!("01:02:03".hms().is_ok());
        assert!("001:02:03".dms().is_ok());
    }

    #[test]
    fn test_error_cases() {
        assert!("not_a_number".deg().is_err());
        assert!("12:34".hms().is_err());
        assert!("12:34:".hms().is_err());
        assert!(":12:34".hms().is_err());

        assert!("".deg().is_err());
        assert!("   ".deg().is_err());
    }

    #[test]
    fn test_sign_handling() {
        assert!("+45:30:15".dms().unwrap().degrees() > 0.0);
        assert!("+12:34:56".hms().unwrap().hours() > 0.0);

        assert!("-45:30:15".dms().unwrap().degrees() < 0.0);
        assert!("-12:34:56".hms().unwrap().hours() < 0.0);

        assert!("45:-30:15".dms().is_err());
        assert!("12:34:-56".hms().is_err());
    }

    #[test]
    fn test_whitespace_tolerance() {
        assert!("  45:30:15  ".dms().is_ok());
        assert!("\t12:34:56\n".hms().is_ok());

        assert!("45 : 30 : 15".dms().is_ok());
        assert!("12 h 34 m 56 s".hms().is_ok());
    }

    #[test]
    fn test_precision_preservation() {
        let input_deg = 123.456789012345;
        let angle = format!("{}", input_deg).deg().unwrap();
        assert!((angle.degrees() - input_deg).abs() < 1e-12);

        let angle = "12:34:56.123456".hms().unwrap();
        let back_to_hms = angle.hours();
        let expected = 12.0 + 34.0 / 60.0 + 56.123456 / 3600.0;
        assert!((back_to_hms - expected).abs() < 1e-9);
    }
}
