//! Core angle type for astronomical calculations.
//!
//! This module provides [`Angle`], the fundamental angular measurement type used throughout
//! the astronomy library. Angles are stored internally as radians (f64) but can be constructed
//! from and converted to degrees, hours, arcminutes, and arcseconds.
//!
//! # Design Rationale
//!
//! **Why radians internally?** All trigonometric functions in Rust (and most languages) operate
//! on radians. Storing radians avoids repeated conversions during calculations. The degree-based
//! constructors and accessors provide ergonomic APIs for human-readable values.
//!
//! **Why associated constants?** [`Angle::PI`], [`Angle::HALF_PI`], and [`Angle::ZERO`] exist
//! because angles are not just numbers. While `std::f64::consts::PI` gives you a raw float,
//! `Angle::PI` gives you a typed angle. This prevents accidentally mixing raw radians with
//! Angles and catches unit errors at compile time.
//!
//! # Quick Start
//!
//! ```
//! use astro_core::Angle;
//!
//! // Construction - pick the unit that matches your data
//! let from_deg = Angle::from_degrees(45.0);
//! let from_rad = Angle::from_radians(0.785398);
//! let from_hrs = Angle::from_hours(3.0);  // 3h = 45 degrees
//! let from_arcsec = Angle::from_arcseconds(162000.0);  // 45 degrees
//!
//! // Conversion - get any unit you need
//! assert!((from_deg.radians() - 0.785398).abs() < 1e-5);
//! assert!((from_deg.hours() - 3.0).abs() < 1e-10);
//!
//! // Trigonometry - no conversion needed
//! let (sin, cos) = from_deg.sin_cos();
//! ```
//!
//! # Hour Angles
//!
//! Astronomy uses hours (0-24h) for right ascension. One hour equals 15 degrees:
//!
//! ```
//! use astro_core::Angle;
//!
//! let ra = Angle::from_hours(6.0);  // 6h RA
//! assert!((ra.degrees() - 90.0).abs() < 1e-10);
//! ```
//!
//! # Validation
//!
//! Angles can be validated for specific astronomical contexts:
//!
//! ```
//! use astro_core::Angle;
//!
//! let dec = Angle::from_degrees(45.0);
//! assert!(dec.validate_declination(false).is_ok());  // -90 to +90
//!
//! let bad_dec = Angle::from_degrees(100.0);
//! assert!(bad_dec.validate_declination(false).is_err());  // Out of range
//!
//! // Right ascension auto-normalizes to [0, 360)
//! let ra = Angle::from_degrees(400.0);
//! let normalized = ra.validate_right_ascension().unwrap();
//! assert!((normalized.degrees() - 40.0).abs() < 1e-10);
//! ```
//!
//! # Convenience Functions
//!
//! For terser code, use the free functions [`deg`], [`rad`], [`hours`], [`arcsec`], [`arcmin`]:
//!
//! ```
//! use astro_core::angle::{deg, hours, arcsec};
//!
//! let a = deg(45.0);
//! let b = hours(3.0);
//! let c = arcsec(162000.0);
//!
//! assert!((a.degrees() - b.degrees()).abs() < 1e-10);
//! ```
//!
//! # Arithmetic
//!
//! Angles support addition, subtraction, negation, and scalar multiplication/division:
//!
//! ```
//! use astro_core::Angle;
//!
//! let a = Angle::from_degrees(30.0);
//! let b = Angle::from_degrees(15.0);
//!
//! let sum = a + b;  // 45 degrees
//! let diff = a - b;  // 15 degrees
//! let scaled = a * 2.0;  // 60 degrees
//! let neg = -a;  // -30 degrees
//! ```

use crate::constants::{HALF_PI, PI};

/// An angular measurement stored as radians.
///
/// `Angle` is the primary type for representing angles throughout this library.
/// It stores the angle as a 64-bit float in radians and provides conversions to/from
/// other angular units commonly used in astronomy.
///
/// # Internal Representation
///
/// Angles are stored as radians (`f64`). This choice optimizes for:
/// - Direct use with trigonometric functions
/// - Precision in intermediate calculations
/// - Consistency with mathematical conventions
///
/// # Derives
///
/// - `Copy`, `Clone`: Angles are small (8 bytes) and cheap to copy
/// - `Debug`: Shows internal radian value
/// - `PartialEq`, `PartialOrd`: Compare angles directly (compares radian values)
///
/// Note: `Eq` and `Ord` are not implemented because f64 can be NaN.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Angle {
    rad: f64,
}

impl Angle {
    /// Zero angle (0 radians).
    pub const ZERO: Self = Self { rad: 0.0 };

    /// Pi radians (180 degrees). Useful for half-circle operations.
    pub const PI: Self = Self { rad: PI };

    /// Pi/2 radians (90 degrees). Useful for right angles and pole declinations.
    pub const HALF_PI: Self = Self { rad: HALF_PI };

    /// Creates an angle from radians.
    ///
    /// This is the only `const` constructor because radians are the internal representation.
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Angle;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// let angle = Angle::from_radians(FRAC_PI_4);
    /// assert!((angle.degrees() - 45.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub const fn from_radians(rad: f64) -> Self {
        Self { rad }
    }

    /// Creates an angle from degrees.
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Angle;
    ///
    /// let angle = Angle::from_degrees(180.0);
    /// assert!((angle.radians() - std::f64::consts::PI).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_degrees(deg: f64) -> Self {
        Self {
            rad: deg.to_radians(),
        }
    }

    /// Creates an angle from hours.
    ///
    /// In astronomy, right ascension is measured in hours where 24h = 360 degrees.
    /// Each hour equals 15 degrees.
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Angle;
    ///
    /// let ra = Angle::from_hours(6.0);  // 6h = 90 degrees
    /// assert!((ra.degrees() - 90.0).abs() < 1e-10);
    ///
    /// let ra_24h = Angle::from_hours(24.0);  // Full circle
    /// assert!((ra_24h.degrees() - 360.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn from_hours(h: f64) -> Self {
        Self {
            rad: (h * 15.0).to_radians(),
        }
    }

    /// Creates an angle from arcseconds.
    ///
    /// One arcsecond = 1/3600 of a degree. Commonly used for:
    /// - Parallax measurements
    /// - Proper motion
    /// - Small angular separations
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Angle;
    ///
    /// let angle = Angle::from_arcseconds(3600.0);  // 1 degree
    /// assert!((angle.degrees() - 1.0).abs() < 1e-10);
    ///
    /// // Proxima Centauri's parallax is about 0.77 arcseconds
    /// let parallax = Angle::from_arcseconds(0.77);
    /// ```
    #[inline]
    pub fn from_arcseconds(arcsec: f64) -> Self {
        Self {
            rad: (arcsec / 3600.0).to_radians(),
        }
    }

    /// Creates an angle from arcminutes.
    ///
    /// One arcminute = 1/60 of a degree. Commonly used for:
    /// - Field of view specifications
    /// - Object sizes (e.g., the Moon is about 31 arcminutes)
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Angle;
    ///
    /// let angle = Angle::from_arcminutes(60.0);  // 1 degree
    /// assert!((angle.degrees() - 1.0).abs() < 1e-10);
    ///
    /// // Full Moon's apparent diameter
    /// let moon_diameter = Angle::from_arcminutes(31.0);
    /// ```
    #[inline]
    pub fn from_arcminutes(arcmin: f64) -> Self {
        Self {
            rad: (arcmin / 60.0).to_radians(),
        }
    }

    /// Returns the angle in radians.
    ///
    /// This is the internal representation, so no conversion occurs.
    #[inline]
    pub fn radians(self) -> f64 {
        self.rad
    }

    /// Returns the angle in degrees.
    #[inline]
    pub fn degrees(self) -> f64 {
        self.rad.to_degrees()
    }

    /// Returns the angle in hours.
    ///
    /// Useful for right ascension where 24h = 360 degrees.
    #[inline]
    pub fn hours(self) -> f64 {
        self.degrees() / 15.0
    }

    /// Returns the angle in arcseconds.
    #[inline]
    pub fn arcseconds(self) -> f64 {
        self.degrees() * 3600.0
    }

    /// Returns the angle in arcminutes.
    #[inline]
    pub fn arcminutes(self) -> f64 {
        self.degrees() * 60.0
    }

    /// Returns the sine of the angle.
    #[inline]
    pub fn sin(self) -> f64 {
        self.rad.sin()
    }

    /// Returns the cosine of the angle.
    #[inline]
    pub fn cos(self) -> f64 {
        self.rad.cos()
    }

    /// Returns both sine and cosine of the angle.
    ///
    /// Convenience method when you need both values.
    ///
    /// # Returns
    ///
    /// A tuple `(sin, cos)`.
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Angle;
    ///
    /// let angle = Angle::from_degrees(30.0);
    /// let (sin, cos) = angle.sin_cos();
    /// assert!((sin - 0.5).abs() < 1e-10);
    /// assert!((cos - 0.866025).abs() < 1e-5);
    /// ```
    #[inline]
    pub fn sin_cos(self) -> (f64, f64) {
        self.rad.sin_cos()
    }

    /// Returns the tangent of the angle.
    #[inline]
    pub fn tan(self) -> f64 {
        self.rad.tan()
    }

    /// Returns the absolute value of the angle.
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Angle;
    ///
    /// let negative = Angle::from_degrees(-45.0);
    /// let absolute = negative.abs();
    /// assert!((absolute.degrees() - 45.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn abs(self) -> Self {
        Self {
            rad: self.rad.abs(),
        }
    }

    /// Wraps the angle to the range [-pi, +pi) (i.e., [-180, +180) degrees).
    ///
    /// Use this for longitude-like quantities or angular differences where
    /// you want the shortest arc representation.
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Angle;
    ///
    /// let angle = Angle::from_degrees(270.0);
    /// let wrapped = angle.wrapped();
    /// assert!((wrapped.degrees() - (-90.0)).abs() < 1e-10);
    ///
    /// let angle2 = Angle::from_degrees(-270.0);
    /// let wrapped2 = angle2.wrapped();
    /// assert!((wrapped2.degrees() - 90.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn wrapped(self) -> Self {
        use super::normalize::wrap_pm_pi;
        Self {
            rad: wrap_pm_pi(self.rad),
        }
    }

    /// Normalizes the angle to the range [0, 2*pi) (i.e., [0, 360) degrees).
    ///
    /// Use this for right ascension or any angle that should be non-negative.
    ///
    /// # Example
    ///
    /// ```
    /// use astro_core::Angle;
    ///
    /// let angle = Angle::from_degrees(-90.0);
    /// let normalized = angle.normalized();
    /// assert!((normalized.degrees() - 270.0).abs() < 1e-10);
    ///
    /// let angle2 = Angle::from_degrees(450.0);
    /// let normalized2 = angle2.normalized();
    /// assert!((normalized2.degrees() - 90.0).abs() < 1e-10);
    /// ```
    #[inline]
    pub fn normalized(self) -> Self {
        Self {
            rad: super::normalize::wrap_0_2pi(self.rad),
        }
    }

    /// Validates the angle as a longitude.
    ///
    /// If `normalize` is true, wraps to [0, 2*pi) and returns Ok.
    /// If `normalize` is false, requires the angle to be in [-pi, +pi] or returns Err.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError`](crate::AstroError) if:
    /// - The angle is not finite (NaN or infinity)
    /// - `normalize` is false and the angle is outside [-180, +180] degrees
    #[inline]
    pub fn validate_longitude(self, normalize: bool) -> Result<Self, crate::AstroError> {
        super::validate::validate_longitude(self, normalize)
    }

    /// Validates the angle as a geographic latitude.
    ///
    /// Latitude must be in [-90, +90] degrees ([-pi/2, +pi/2] radians).
    ///
    /// # Errors
    ///
    /// Returns [`AstroError`](crate::AstroError) if:
    /// - The angle is not finite (NaN or infinity)
    /// - The angle is outside [-90, +90] degrees
    #[inline]
    pub fn validate_latitude(self) -> Result<Self, crate::AstroError> {
        super::validate::validate_latitude(self)
    }

    /// Validates the angle as a declination.
    ///
    /// - `beyond_pole = false`: standard range [-90°, +90°]
    /// - `beyond_pole = true`: extended range [-180°, +180°] for GEM pier-flipped observations
    ///
    /// # Errors
    ///
    /// Returns [`AstroError`](crate::AstroError) if:
    /// - The angle is not finite (NaN or infinity)
    /// - The angle is outside the valid range
    #[inline]
    pub fn validate_declination(self, beyond_pole: bool) -> Result<Self, crate::AstroError> {
        super::validate::validate_declination(self, beyond_pole)
    }

    /// Validates the angle as a right ascension, normalizing to [0, 360) degrees.
    ///
    /// Unlike declination, right ascension is cyclic. This method accepts any finite angle
    /// and normalizes it to [0, 2*pi).
    ///
    /// # Errors
    ///
    /// Returns [`AstroError`](crate::AstroError) if the angle is not finite (NaN or infinity).
    #[inline]
    pub fn validate_right_ascension(self) -> Result<Self, crate::AstroError> {
        super::validate::validate_right_ascension(self)
    }
}

/// Creates an angle from radians. Shorthand for [`Angle::from_radians`].
///
/// # Example
///
/// ```
/// use astro_core::angle::rad;
/// use std::f64::consts::PI;
///
/// let angle = rad(PI);
/// assert!((angle.degrees() - 180.0).abs() < 1e-10);
/// ```
#[inline]
pub fn rad(v: f64) -> Angle {
    Angle::from_radians(v)
}

/// Creates an angle from degrees. Shorthand for [`Angle::from_degrees`].
///
/// # Example
///
/// ```
/// use astro_core::angle::deg;
///
/// let angle = deg(45.0);
/// assert!((angle.radians() - std::f64::consts::FRAC_PI_4).abs() < 1e-10);
/// ```
#[inline]
pub fn deg(v: f64) -> Angle {
    Angle::from_degrees(v)
}

/// Creates an angle from hours. Shorthand for [`Angle::from_hours`].
///
/// # Example
///
/// ```
/// use astro_core::angle::hours;
///
/// let ra = hours(6.0);  // 6h = 90 degrees
/// assert!((ra.degrees() - 90.0).abs() < 1e-10);
/// ```
#[inline]
pub fn hours(v: f64) -> Angle {
    Angle::from_hours(v)
}

/// Creates an angle from arcseconds. Shorthand for [`Angle::from_arcseconds`].
///
/// # Example
///
/// ```
/// use astro_core::angle::arcsec;
///
/// let angle = arcsec(3600.0);  // 1 degree
/// assert!((angle.degrees() - 1.0).abs() < 1e-10);
/// ```
#[inline]
pub fn arcsec(v: f64) -> Angle {
    Angle::from_degrees(v / 3600.0)
}

/// Creates an angle from arcminutes. Shorthand for [`Angle::from_arcminutes`].
///
/// # Example
///
/// ```
/// use astro_core::angle::arcmin;
///
/// let angle = arcmin(60.0);  // 1 degree
/// assert!((angle.degrees() - 1.0).abs() < 1e-10);
/// ```
#[inline]
pub fn arcmin(v: f64) -> Angle {
    Angle::from_degrees(v / 60.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_arcseconds() {
        let angle = Angle::from_arcseconds(3600.0);
        assert!((angle.degrees() - 1.0).abs() < 1e-20);
    }

    #[test]
    fn test_from_arcminutes() {
        let angle = Angle::from_arcminutes(60.0);
        assert!((angle.degrees() - 1.0).abs() < 1e-20);
    }

    #[test]
    fn test_arcseconds_getter() {
        let angle = Angle::from_degrees(1.0);
        assert!((angle.arcseconds() - 3600.0).abs() < 1e-20);
    }

    #[test]
    fn test_arcminutes_getter() {
        let angle = Angle::from_degrees(1.0);
        assert!((angle.arcminutes() - 60.0).abs() < 1e-20);
    }

    #[test]
    fn test_sin() {
        let angle = Angle::from_degrees(30.0);
        assert!((angle.sin() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_tan() {
        let angle = Angle::from_degrees(45.0);
        assert!((angle.tan() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_helper_functions() {
        let a = rad(crate::constants::PI);
        assert!((a.degrees() - 180.0).abs() < 1e-20);

        let b = deg(90.0);
        assert!((b.radians() - crate::constants::HALF_PI).abs() < 1e-20);

        let c = hours(12.0);
        assert!((c.degrees() - 180.0).abs() < 1e-20);

        let d = arcsec(3600.0);
        assert!((d.degrees() - 1.0).abs() < 1e-20);

        let e = arcmin(60.0);
        assert!((e.degrees() - 1.0).abs() < 1e-20);
    }
}
