use super::Angle;

/// Creates an angle from radians. Shorthand for [`Angle::from_radians`].
///
/// # Example
///
/// ```
/// use celestial_core::angle::rad;
/// use celestial_core::constants::PI;
///
/// let angle = rad(PI);
/// assert_eq!(angle.degrees(), 180.0);
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
/// use celestial_core::angle::deg;
/// use celestial_core::constants::QUARTER_PI;
///
/// let angle = deg(45.0);
/// assert_eq!(angle.radians(), QUARTER_PI);
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
/// use celestial_core::angle::hours;
///
/// let ra = hours(6.0);  // 6h = 90 degrees
/// assert_eq!(ra.degrees(), 90.0);
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
/// use celestial_core::angle::arcsec;
///
/// let angle = arcsec(3600.0);  // 1 degree
/// assert_eq!(angle.degrees(), 1.0);
/// ```
#[inline]
pub fn arcsec(v: f64) -> Angle {
    Angle::from_arcseconds(v)
}

/// Creates an angle from arcminutes. Shorthand for [`Angle::from_arcminutes`].
///
/// # Example
///
/// ```
/// use celestial_core::angle::arcmin;
///
/// let angle = arcmin(60.0);  // 1 degree
/// assert_eq!(angle.degrees(), 1.0);
/// ```
#[inline]
pub fn arcmin(v: f64) -> Angle {
    Angle::from_arcminutes(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_helper_functions() {
        let a = rad(crate::constants::PI);
        assert_eq!(a.degrees(), 180.0);

        let b = deg(90.0);
        assert_eq!(b.radians(), crate::constants::HALF_PI);

        let c = hours(12.0);
        assert_eq!(c.degrees(), 180.0);

        let d = arcsec(3600.0);
        assert_eq!(d.degrees(), 1.0);

        let e = arcmin(60.0);
        assert_eq!(e.degrees(), 1.0);
    }

    #[test]
    fn test_unit_shorthands_match_constructors() {
        for v in [0.371, 1.0, 59.9, 1234.5678] {
            assert_eq!(arcsec(v), Angle::from_arcseconds(v), "{v}");
            assert_eq!(arcmin(v), Angle::from_arcminutes(v), "{v}");
        }
    }
}
