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
#[must_use]
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
    /// use celestial_core::angle::Angle;
    /// use celestial_core::constants::QUARTER_PI;
    ///
    /// let angle = Angle::from_radians(QUARTER_PI);
    /// assert_eq!(angle.degrees(), 45.0);
    /// ```
    #[inline]
    pub const fn from_radians(rad: f64) -> Self {
        Self { rad }
    }

    /// Returns the angle in radians.
    ///
    /// This is the internal representation, so no conversion occurs.
    #[inline]
    pub fn radians(self) -> f64 {
        self.rad
    }

    /// Returns the sine of the angle.
    #[inline]
    pub fn sin(self) -> f64 {
        libm::sin(self.rad)
    }

    /// Returns the cosine of the angle.
    #[inline]
    pub fn cos(self) -> f64 {
        libm::cos(self.rad)
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
    /// use celestial_core::angle::Angle;
    ///
    /// let angle = Angle::from_degrees(30.0);
    /// let (sin, cos) = angle.sin_cos();
    /// assert_eq!(sin, 0.5);
    /// assert_eq!(cos, 0.8660254037844386);
    /// ```
    #[inline]
    pub fn sin_cos(self) -> (f64, f64) {
        libm::sincos(self.rad)
    }

    /// Returns the tangent of the angle.
    #[inline]
    pub fn tan(self) -> f64 {
        libm::tan(self.rad)
    }

    /// Returns the absolute value of the angle.
    ///
    /// # Example
    ///
    /// ```
    /// use celestial_core::angle::Angle;
    ///
    /// let negative = Angle::from_degrees(-45.0);
    /// let absolute = negative.abs();
    /// assert_eq!(absolute.degrees(), 45.0);
    /// ```
    #[inline]
    pub fn abs(self) -> Self {
        Self {
            rad: libm::fabs(self.rad),
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
    /// use celestial_core::angle::Angle;
    ///
    /// let angle = Angle::from_degrees(270.0);
    /// let wrapped = angle.wrapped()?;
    /// assert_eq!(wrapped.degrees(), -90.0);
    ///
    /// let angle2 = Angle::from_degrees(-270.0);
    /// let wrapped2 = angle2.wrapped()?;
    /// assert_eq!(wrapped2.degrees(), 90.0);
    /// # Ok::<(), celestial_core::errors::AstroError>(())
    /// ```
    #[inline]
    pub fn wrapped(self) -> crate::errors::AstroResult<Self> {
        super::normalize::wrap_pm_pi(self.rad).map(Self::from_radians)
    }

    /// Normalizes the angle to the range [0, 2*pi) (i.e., [0, 360) degrees).
    ///
    /// Use this for right ascension or any angle that should be non-negative.
    ///
    /// # Example
    ///
    /// ```
    /// use celestial_core::angle::Angle;
    ///
    /// let angle = Angle::from_degrees(-90.0);
    /// let normalized = angle.normalized()?;
    /// assert_eq!(normalized.degrees(), 270.0);
    ///
    /// let angle2 = Angle::from_degrees(450.0);
    /// let normalized2 = angle2.normalized()?;
    /// assert_eq!(normalized2.degrees(), 90.0);
    /// # Ok::<(), celestial_core::errors::AstroError>(())
    /// ```
    #[inline]
    pub fn normalized(self) -> crate::errors::AstroResult<Self> {
        super::normalize::wrap_0_2pi(self.rad).map(Self::from_radians)
    }

    /// Validates the angle as a geographic latitude.
    ///
    /// Latitude must be in [-90, +90] degrees ([-pi/2, +pi/2] radians).
    ///
    /// # Errors
    ///
    /// Returns [`AstroError`](crate::errors::AstroError) if:
    /// - The angle is not finite (NaN or infinity)
    /// - The angle is outside [-90, +90] degrees
    #[inline]
    pub fn validate_latitude(self) -> Result<Self, crate::errors::AstroError> {
        super::validate::validate_latitude(self)
    }

    /// Validates the angle as a declination.
    ///
    /// - `beyond_pole = false`: standard range [-90°, +90°]
    /// - `beyond_pole = true`: extended range [-180°, +180°] for GEM pier-flipped observations
    ///
    /// # Errors
    ///
    /// Returns [`AstroError`](crate::errors::AstroError) if:
    /// - The angle is not finite (NaN or infinity)
    /// - The angle is outside the valid range
    #[inline]
    pub fn validate_declination(
        self,
        beyond_pole: bool,
    ) -> Result<Self, crate::errors::AstroError> {
        super::validate::validate_declination(self, beyond_pole)
    }

    /// Validates the angle as a right ascension, normalizing to [0, 360) degrees.
    ///
    /// Unlike declination, right ascension is cyclic. This method accepts any finite angle
    /// and normalizes it to [0, 2*pi).
    ///
    /// # Errors
    ///
    /// Returns [`AstroError`](crate::errors::AstroError) if the angle is not finite (NaN or infinity).
    #[inline]
    pub fn validate_right_ascension(self) -> Result<Self, crate::errors::AstroError> {
        super::validate::validate_right_ascension(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wrapping_rejects_non_finite_angle() {
        for rad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(Angle::from_radians(rad).wrapped().is_err());
            assert!(Angle::from_radians(rad).normalized().is_err());
        }
    }

    #[test]
    fn test_sin() {
        let angle = Angle::from_degrees(30.0);
        assert_eq!(angle.sin(), 0.5);
    }

    #[test]
    fn test_cos() {
        assert_eq!(Angle::from_degrees(0.0).cos(), 1.0);
        // The nearest f64 to π/3 sits just above the true value; its cosine is two ulps under 0.5.
        assert_eq!(Angle::from_degrees(60.0).cos(), 0.4999999999999999);
    }

    #[test]
    fn test_validate_latitude_delegates() {
        let angle = Angle::from_radians(0.75);
        assert_eq!(angle.validate_latitude().unwrap().radians(), 0.75);

        let err = Angle::from_degrees(95.0).validate_latitude().unwrap_err();
        assert_eq!(
            err.to_string(),
            "Math error in validate_latitude (out of range): Latitude 95.00° out of range [-90°, +90°]"
        );
    }

    #[test]
    fn test_tan() {
        let angle = Angle::from_degrees(45.0);
        // f64 π/4 sits just below the true value, so its tangent rounds to one ulp under 1.
        assert_eq!(angle.tan(), 0.9999999999999999);
    }
}
