/// Pixel types the detector can consume.
///
/// Implemented for `u8`, `i8`, `u16`, `i16`, `u32`, `i32`, `u64`, `i64`, `f32`, `f64`.
/// [`Pixel::to_f64`] converts a sample for arithmetic. The saturation threshold is
/// estimated from the data for every type; see
/// [`crate::detect::DetectionParams::saturation_adu`].
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::pixel::Pixel;
///
/// assert_eq!(42u16.to_f64(), 42.0);
/// ```
pub trait Pixel: Copy + PartialOrd + Send + Sync + 'static {
    /// Converts the pixel value to `f64` for arithmetic.
    fn to_f64(self) -> f64;
}

impl Pixel for u8 {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl Pixel for i8 {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl Pixel for u16 {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl Pixel for i16 {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl Pixel for u32 {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl Pixel for i32 {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl Pixel for u64 {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl Pixel for i64 {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl Pixel for f32 {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl Pixel for f64 {
    fn to_f64(self) -> f64 {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_f64_preserves_value_for_every_type() {
        assert_eq!(Pixel::to_f64(42u8), 42.0);
        assert_eq!(Pixel::to_f64(40000u16), 40000.0);
        assert_eq!(Pixel::to_f64(-100i16), -100.0);
        assert_eq!(Pixel::to_f64(1_000_000i32), 1_000_000.0);
        assert_eq!(Pixel::to_f64(-100i8), -100.0);
        assert_eq!(Pixel::to_f64(3_000_000_000u32), 3_000_000_000.0);
        assert_eq!(Pixel::to_f64(1u64 << 60), 1152921504606846976.0);
        assert_eq!(Pixel::to_f64(-(1i64 << 60)), -1152921504606846976.0);
        assert_eq!(Pixel::to_f64(3.5f32), 3.5);
        assert_eq!(Pixel::to_f64(1.5f64), 1.5);
    }
}
