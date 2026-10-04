//! Arithmetic operators for [`Angle`].
//!
//! Implements standard math ops: `+`, `-`, `*`, `/`, and unary `-`.

use super::Angle;
use core::ops::*;

impl Add for Angle {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::from_radians(self.radians() + rhs.radians())
    }
}

impl Sub for Angle {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::from_radians(self.radians() - rhs.radians())
    }
}

impl Mul<f64> for Angle {
    type Output = Self;

    #[inline]
    fn mul(self, k: f64) -> Self {
        Self::from_radians(self.radians() * k)
    }
}

impl Div<f64> for Angle {
    type Output = Self;

    #[inline]
    fn div(self, k: f64) -> Self {
        Self::from_radians(self.radians() / k)
    }
}

impl Neg for Angle {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::from_radians(-self.radians())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_sub() {
        let a = Angle::from_radians(1.0);
        let b = Angle::from_radians(0.5);
        assert_eq!((a + b).radians(), 1.5);
        assert_eq!((a - b).radians(), 0.5);
    }

    #[test]
    fn test_mul_div() {
        let a = Angle::from_radians(1.0);
        assert_eq!((a * 2.0).radians(), 2.0);
        assert_eq!((a / 2.0).radians(), 0.5);
    }

    #[test]
    fn test_neg() {
        let a = Angle::from_radians(1.0);
        assert_eq!((-a).radians(), -1.0);
        assert_eq!((-(-a)).radians(), 1.0);
    }
}
