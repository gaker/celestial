//! Arithmetic operators for [`Angle`].
//!
//! Implements standard math ops: `+`, `-`, `*`, `/`, and unary `-`.

use super::core::Angle;
use core::ops::*;

/// Angle + Angle → Angle
impl Add for Angle {
    type Output = Angle;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Angle::from_radians(self.radians() + rhs.radians())
    }
}

/// Angle - Angle → Angle
impl Sub for Angle {
    type Output = Angle;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Angle::from_radians(self.radians() - rhs.radians())
    }
}

/// Angle * scalar → Angle
impl Mul<f64> for Angle {
    type Output = Angle;
    #[inline]
    fn mul(self, k: f64) -> Self {
        Angle::from_radians(self.radians() * k)
    }
}

/// Angle / scalar → Angle
impl Div<f64> for Angle {
    type Output = Angle;
    #[inline]
    fn div(self, k: f64) -> Self {
        Angle::from_radians(self.radians() / k)
    }
}

/// -Angle → Angle
impl Neg for Angle {
    type Output = Angle;
    #[inline]
    fn neg(self) -> Self {
        Angle::from_radians(-self.radians())
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
