/// Pixel types the detector can consume.
///
/// Implemented for `u8`, `u16`, `i16`, `i32`, `f32`, `f64`. Two methods:
/// [`Pixel::to_f64`] for arithmetic, [`Pixel::saturation_ceiling`] for the saturation
/// threshold used during centroiding. Integer types return their native max; float
/// types return infinity, which falls through to a tied-maximum heuristic in
/// [`crate::detect::find_bright_stars`].
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::Pixel;
///
/// assert_eq!(u16::saturation_ceiling(), 65535.0);
/// assert_eq!(42u16.to_f64(), 42.0);
/// ```
pub trait Pixel: Copy + PartialOrd + Send + Sync + 'static {
    /// Converts the pixel value to `f64` for arithmetic.
    fn to_f64(self) -> f64;

    /// Native saturation threshold for this pixel type.
    ///
    /// Integer types return their MAX as `f64` — `u16::MAX`, `i32::MAX`, etc. Float
    /// types return `f64::INFINITY`, which signals the detector to estimate
    /// saturation from the data (tied-maximum heuristic).
    fn saturation_ceiling() -> f64;
}

impl Pixel for u8 {
    fn to_f64(self) -> f64 {
        self as f64
    }
    fn saturation_ceiling() -> f64 {
        Self::MAX as f64
    }
}

impl Pixel for u16 {
    fn to_f64(self) -> f64 {
        self as f64
    }
    fn saturation_ceiling() -> f64 {
        Self::MAX as f64
    }
}

impl Pixel for i16 {
    fn to_f64(self) -> f64 {
        self as f64
    }
    fn saturation_ceiling() -> f64 {
        Self::MAX as f64
    }
}

impl Pixel for i32 {
    fn to_f64(self) -> f64 {
        self as f64
    }
    fn saturation_ceiling() -> f64 {
        Self::MAX as f64
    }
}

impl Pixel for f32 {
    fn to_f64(self) -> f64 {
        self as f64
    }
    fn saturation_ceiling() -> f64 {
        f64::INFINITY
    }
}

impl Pixel for f64 {
    fn to_f64(self) -> f64 {
        self
    }
    fn saturation_ceiling() -> f64 {
        Self::INFINITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_types_ceilings_match_native_max() {
        assert_eq!(<u8 as Pixel>::saturation_ceiling(), 255.0);
        assert_eq!(<u16 as Pixel>::saturation_ceiling(), 65535.0);
        assert_eq!(<i16 as Pixel>::saturation_ceiling(), 32767.0);
        assert_eq!(<i32 as Pixel>::saturation_ceiling(), i32::MAX as f64);
    }

    #[test]
    fn float_types_ceilings_are_infinity() {
        assert_eq!(<f32 as Pixel>::saturation_ceiling(), f64::INFINITY);
        assert_eq!(<f64 as Pixel>::saturation_ceiling(), f64::INFINITY);
    }

    #[test]
    fn to_f64_preserves_value_for_every_type() {
        assert_eq!(Pixel::to_f64(42u8), 42.0);
        assert_eq!(Pixel::to_f64(40000u16), 40000.0);
        assert_eq!(Pixel::to_f64(-100i16), -100.0);
        assert_eq!(Pixel::to_f64(1_000_000i32), 1_000_000.0);
        assert_eq!(Pixel::to_f64(3.5f32), 3.5);
        assert_eq!(Pixel::to_f64(1.5f64), 1.5);
    }
}
