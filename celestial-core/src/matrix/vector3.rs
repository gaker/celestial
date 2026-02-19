//! 3D Cartesian vectors for astronomical coordinate calculations.
//!
//! Vectors are the workhorses of celestial coordinate math. When you transform a star's
//! position between reference frames, compute the angle between two objects, or calculate
//! parallax corrections, you're working with 3D vectors under the hood.
//!
//! # Cartesian vs Spherical
//!
//! Astronomical positions are usually given as spherical coordinates (RA/Dec, Az/Alt,
//! longitude/latitude), but transformations are cleanest in Cartesian form. The typical
//! workflow is:
//!
//! 1. Convert spherical → Cartesian with [`from_spherical`](Vector3::from_spherical)
//! 2. Apply rotation matrices for frame transformations
//! 3. Convert back with [`to_spherical`](Vector3::to_spherical)
//!
//! ```
//! use celestial_core::Vector3;
//! use std::f64::consts::FRAC_PI_4;
//!
//! // A star at RA=45°, Dec=30° (in radians)
//! let ra = FRAC_PI_4;
//! let dec = FRAC_PI_4 / 1.5;  // ~30°
//!
//! let cartesian = Vector3::from_spherical(ra, dec);
//! // Now apply rotations, then convert back:
//! let (new_ra, new_dec) = cartesian.to_spherical();
//! ```
//!
//! # Unit Vectors and Direction
//!
//! For celestial positions on the unit sphere (where distance doesn't matter), vectors
//! are normalized to unit length. The [`normalize`](Vector3::normalize) method returns
//! a unit vector pointing in the same direction:
//!
//! ```
//! use celestial_core::Vector3;
//!
//! let v = Vector3::new(3.0, 4.0, 0.0);
//! let unit = v.normalize();
//! assert!((unit.magnitude() - 1.0).abs() < 1e-15);
//! ```
//!
//! # Dot and Cross Products
//!
//! These operations have direct astronomical applications:
//!
//! - **Dot product**: Compute the cosine of the angle between two directions.
//!   For unit vectors, `a.dot(&b)` equals `cos(θ)` where θ is the separation angle.
//!
//! - **Cross product**: Find the axis perpendicular to two directions.
//!   Useful for computing rotation axes and angular momentum vectors.
//!
//! ```
//! use celestial_core::Vector3;
//!
//! let a = Vector3::x_axis();  // Points along +X
//! let b = Vector3::y_axis();  // Points along +Y
//!
//! // Perpendicular: dot product is zero
//! assert_eq!(a.dot(&b), 0.0);
//!
//! // Cross product gives +Z axis (right-hand rule)
//! let c = a.cross(&b);
//! assert_eq!(c, Vector3::z_axis());
//! ```
//!
//! # Coordinate Conventions
//!
//! The spherical coordinate convention used here matches standard astronomical practice:
//! - **θ (theta)**: Azimuthal angle from +X axis toward +Y axis (like right ascension)
//! - **φ (phi)**: Elevation angle from XY plane (like declination)
//!
//! This differs from the physics convention where φ is the azimuthal angle and θ is
//! the polar angle from +Z.
use crate::{AstroError, AstroResult, MathErrorKind};
use std::fmt;

/// A 3D Cartesian vector for coordinate calculations.
///
/// Used throughout the library for position vectors, direction vectors, and as
/// intermediate representations during coordinate transformations.
///
/// # Fields
///
/// Components are public for direct access when performance matters:
/// - `x`: First component (toward vernal equinox in equatorial coordinates)
/// - `y`: Second component (90° east in equatorial coordinates)
/// - `z`: Third component (toward celestial pole in equatorial coordinates)
///
/// # Construction
///
/// ```
/// use celestial_core::Vector3;
///
/// // Direct construction
/// let v = Vector3::new(1.0, 2.0, 3.0);
///
/// // Unit vectors along axes
/// let x = Vector3::x_axis();
/// let y = Vector3::y_axis();
/// let z = Vector3::z_axis();
///
/// // From spherical coordinates (RA, Dec in radians)
/// let star = Vector3::from_spherical(0.5, 0.3);
///
/// // From an array
/// let v = Vector3::from_array([1.0, 2.0, 3.0]);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3 {
    /// Creates a new vector from x, y, z components.
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Returns the zero vector `[0, 0, 0]`.
    #[inline]
    pub fn zeros() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Returns the unit vector along the X axis `[1, 0, 0]`.
    ///
    /// In equatorial coordinates, this points toward the vernal equinox.
    #[inline]
    pub fn x_axis() -> Self {
        Self::new(1.0, 0.0, 0.0)
    }

    /// Returns the unit vector along the Y axis `[0, 1, 0]`.
    ///
    /// In equatorial coordinates, this is 90° east of the vernal equinox on the equator.
    #[inline]
    pub fn y_axis() -> Self {
        Self::new(0.0, 1.0, 0.0)
    }

    /// Returns the unit vector along the Z axis `[0, 0, 1]`.
    ///
    /// In equatorial coordinates, this points toward the north celestial pole.
    #[inline]
    pub fn z_axis() -> Self {
        Self::new(0.0, 0.0, 1.0)
    }

    /// Returns the component at the given index (0=x, 1=y, 2=z).
    ///
    /// Returns an error for indices outside 0-2. For unchecked access, use
    /// indexing syntax `v[i]` or the public fields directly.
    pub fn get(&self, index: usize) -> AstroResult<f64> {
        match index {
            0 => Ok(self.x),
            1 => Ok(self.y),
            2 => Ok(self.z),
            _ => Err(AstroError::math_error(
                "Vector3::get",
                MathErrorKind::InvalidInput,
                &format!("index {} out of bounds (valid range: 0-2)", index),
            )),
        }
    }

    /// Sets the component at the given index (0=x, 1=y, 2=z).
    ///
    /// Returns an error for indices outside 0-2. For unchecked access, use
    /// indexing syntax `v[i] = value` or the public fields directly.
    pub fn set(&mut self, index: usize, value: f64) -> AstroResult<()> {
        match index {
            0 => {
                self.x = value;
                Ok(())
            }
            1 => {
                self.y = value;
                Ok(())
            }
            2 => {
                self.z = value;
                Ok(())
            }
            _ => Err(AstroError::math_error(
                "Vector3::set",
                MathErrorKind::InvalidInput,
                &format!("index {} out of bounds (valid range: 0-2)", index),
            )),
        }
    }

    /// Returns the Euclidean length (L2 norm) of the vector.
    ///
    /// For a unit vector, this returns 1.0. For the zero vector, returns 0.0.
    #[inline]
    pub fn magnitude(&self) -> f64 {
        libm::sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
    }

    /// Returns the squared magnitude.
    ///
    /// Faster than [`magnitude`](Self::magnitude) when you only need to compare
    /// lengths or don't need the actual distance.
    #[inline]
    pub fn magnitude_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Returns a unit vector pointing in the same direction.
    ///
    /// If the vector has zero length, returns the zero vector unchanged (avoids NaN).
    ///
    /// ```
    /// use celestial_core::Vector3;
    ///
    /// let v = Vector3::new(3.0, 4.0, 0.0);
    /// let unit = v.normalize();
    /// assert!((unit.magnitude() - 1.0).abs() < 1e-15);
    /// assert_eq!(unit, Vector3::new(0.6, 0.8, 0.0));
    /// ```
    pub fn normalize(&self) -> Self {
        let mag = self.magnitude();
        if mag == 0.0 {
            *self
        } else {
            Self::new(self.x / mag, self.y / mag, self.z / mag)
        }
    }

    /// Computes the dot product (inner product) with another vector.
    ///
    /// For unit vectors, this equals the cosine of the angle between them:
    /// `a.dot(&b) = cos(θ)`. Use this to compute angular separation between
    /// celestial positions.
    ///
    /// ```
    /// use celestial_core::Vector3;
    ///
    /// let a = Vector3::x_axis();
    /// let b = Vector3::y_axis();
    /// assert_eq!(a.dot(&b), 0.0);  // Perpendicular
    ///
    /// let c = Vector3::new(1.0, 2.0, 3.0);
    /// let d = Vector3::new(4.0, 5.0, 6.0);
    /// assert_eq!(c.dot(&d), 32.0);  // 1*4 + 2*5 + 3*6
    /// ```
    #[inline]
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Computes the cross product with another vector.
    ///
    /// The result is perpendicular to both input vectors, with direction given
    /// by the right-hand rule. The magnitude equals `|a||b|sin(θ)`.
    ///
    /// ```
    /// use celestial_core::Vector3;
    ///
    /// let x = Vector3::x_axis();
    /// let y = Vector3::y_axis();
    /// let z = x.cross(&y);
    /// assert_eq!(z, Vector3::z_axis());  // X × Y = Z
    /// ```
    pub fn cross(&self, other: &Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Returns the components as a `[f64; 3]` array.
    #[inline]
    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Creates a vector from a `[f64; 3]` array.
    #[inline]
    pub fn from_array(arr: [f64; 3]) -> Self {
        Self::new(arr[0], arr[1], arr[2])
    }

    /// Creates a unit vector from spherical coordinates.
    ///
    /// - `ra`: Azimuthal angle from +X toward +Y (right ascension), in radians
    /// - `dec`: Elevation from XY plane (declination), in radians
    ///
    /// The result is always a unit vector (magnitude = 1).
    ///
    /// ```
    /// use celestial_core::Vector3;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// // RA=0, Dec=0 → points along +X
    /// let v = Vector3::from_spherical(0.0, 0.0);
    /// assert!((v.x - 1.0).abs() < 1e-15);
    ///
    /// // RA=90°, Dec=0 → points along +Y
    /// let v = Vector3::from_spherical(FRAC_PI_2, 0.0);
    /// assert!((v.y - 1.0).abs() < 1e-15);
    ///
    /// // RA=0, Dec=90° → points along +Z (north pole)
    /// let v = Vector3::from_spherical(0.0, FRAC_PI_2);
    /// assert!((v.z - 1.0).abs() < 1e-15);
    /// ```
    pub fn from_spherical(ra: f64, dec: f64) -> Self {
        let (sin_ra, cos_ra) = libm::sincos(ra);
        let (sin_dec, cos_dec) = libm::sincos(dec);
        Self::new(cos_dec * cos_ra, cos_dec * sin_ra, sin_dec)
    }

    /// Converts the vector to spherical coordinates (θ, φ).
    ///
    /// Returns `(theta, phi)` where:
    /// - `theta`: Azimuthal angle from +X toward +Y (like RA), in radians `(-π, π]`
    /// - `phi`: Elevation from XY plane (like Dec), in radians `[-π/2, π/2]`
    ///
    /// The vector does not need to be normalized; direction is preserved regardless
    /// of magnitude. For the zero vector, returns `(0.0, 0.0)`.
    ///
    /// ```
    /// use celestial_core::Vector3;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let v = Vector3::new(0.0, 0.0, 1.0);  // North pole
    /// let (theta, phi) = v.to_spherical();
    /// assert_eq!(theta, 0.0);
    /// assert_eq!(phi, FRAC_PI_2);
    /// ```
    pub fn to_spherical(&self) -> (f64, f64) {
        let d2 = self.x * self.x + self.y * self.y;

        let theta = if d2 == 0.0 {
            0.0
        } else {
            libm::atan2(self.y, self.x)
        };
        let phi = if self.z == 0.0 {
            0.0
        } else {
            libm::atan2(self.z, libm::sqrt(d2))
        };

        (theta, phi)
    }
}

/// Vector + Vector
impl std::ops::Add for Vector3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

/// Vector - Vector
impl std::ops::Sub for Vector3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

/// Vector * scalar
impl std::ops::Mul<f64> for Vector3 {
    type Output = Self;

    fn mul(self, scalar: f64) -> Self {
        Self::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

/// scalar * Vector
impl std::ops::Mul<Vector3> for f64 {
    type Output = Vector3;

    fn mul(self, vec: Vector3) -> Vector3 {
        vec * self
    }
}

/// Vector / scalar
impl std::ops::Div<f64> for Vector3 {
    type Output = Self;

    fn div(self, scalar: f64) -> Self {
        Self::new(self.x / scalar, self.y / scalar, self.z / scalar)
    }
}

/// Vector /= scalar
impl std::ops::DivAssign<f64> for Vector3 {
    fn div_assign(&mut self, scalar: f64) {
        self.x /= scalar;
        self.y /= scalar;
        self.z /= scalar;
    }
}

/// -Vector
impl std::ops::Neg for Vector3 {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

/// v[i] indexing (panics if i > 2)
impl std::ops::Index<usize> for Vector3 {
    type Output = f64;

    fn index(&self, index: usize) -> &f64 {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Vector3 index out of bounds: {}", index),
        }
    }
}

/// v[i] = value mutable indexing (panics if i > 2)
impl std::ops::IndexMut<usize> for Vector3 {
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Vector3 index out of bounds: {}", index),
        }
    }
}

impl fmt::Display for Vector3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Vector3({:.9}, {:.9}, {:.9})", self.x, self.y, self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vector3_construction() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert_eq!(v.x, 1.0);
        assert_eq!(v.y, 2.0);
        assert_eq!(v.z, 3.0);

        let zeros = Vector3::zeros();
        assert_eq!(zeros.x, 0.0);
        assert_eq!(zeros.y, 0.0);
        assert_eq!(zeros.z, 0.0);

        let x_axis = Vector3::x_axis();
        assert_eq!(x_axis, Vector3::new(1.0, 0.0, 0.0));

        let from_array = Vector3::from_array([4.0, 5.0, 6.0]);
        assert_eq!(from_array, Vector3::new(4.0, 5.0, 6.0));
    }

    #[test]
    fn test_vector3_magnitude() {
        let v = Vector3::new(3.0, 4.0, 0.0);
        assert_eq!(v.magnitude(), 5.0);
        assert_eq!(v.magnitude_squared(), 25.0);

        let unit = v.normalize();
        assert!((unit.magnitude() - 1.0).abs() < 1e-15);
        assert_eq!(unit, Vector3::new(0.6, 0.8, 0.0));
    }

    #[test]
    fn test_vector3_arithmetic() {
        let a = Vector3::new(1.0, 2.0, 3.0);
        let b = Vector3::new(4.0, 5.0, 6.0);

        let sum = a + b;
        assert_eq!(sum, Vector3::new(5.0, 7.0, 9.0));

        let diff = b - a;
        assert_eq!(diff, Vector3::new(3.0, 3.0, 3.0));

        let scaled = a * 2.0;
        assert_eq!(scaled, Vector3::new(2.0, 4.0, 6.0));

        let scaled2 = 3.0 * a;
        assert_eq!(scaled2, Vector3::new(3.0, 6.0, 9.0));

        let divided = a / 2.0;
        assert_eq!(divided, Vector3::new(0.5, 1.0, 1.5));

        let negated = -a;
        assert_eq!(negated, Vector3::new(-1.0, -2.0, -3.0));
    }

    #[test]
    fn test_vector3_dot_cross() {
        let a = Vector3::new(1.0, 0.0, 0.0);
        let b = Vector3::new(0.0, 1.0, 0.0);

        assert_eq!(a.dot(&b), 0.0);

        let c = a.cross(&b);
        assert_eq!(c, Vector3::new(0.0, 0.0, 1.0));

        let d = Vector3::new(1.0, 2.0, 3.0);
        let e = Vector3::new(4.0, 5.0, 6.0);
        assert_eq!(d.dot(&e), 32.0);
    }

    #[test]
    fn test_vector3_spherical_conversion() {
        let v1 = Vector3::from_spherical(0.0, 0.0);
        assert!((v1.x - 1.0).abs() < 1e-15);
        assert!(v1.y.abs() < 1e-15);
        assert!(v1.z.abs() < 1e-15);

        let (ra, dec) = v1.to_spherical();
        assert!(ra.abs() < 1e-15);
        assert!(dec.abs() < 1e-15);

        let v2 = Vector3::from_spherical(crate::constants::HALF_PI, 0.0);
        assert!(v2.x.abs() < 1e-15);
        assert!((v2.y - 1.0).abs() < 1e-15);
        assert!(v2.z.abs() < 1e-15);

        let v3 = Vector3::from_spherical(0.0, crate::constants::HALF_PI);
        assert!(v3.x.abs() < 1e-15);
        assert!(v3.y.abs() < 1e-15);
        assert!((v3.z - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_axis_constructors() {
        // Test y_axis and z_axis constructors
        let y_axis = Vector3::y_axis();
        assert_eq!(y_axis, Vector3::new(0.0, 1.0, 0.0));

        let z_axis = Vector3::z_axis();
        assert_eq!(z_axis, Vector3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn test_get_set_methods() {
        let mut v = Vector3::new(1.0, 2.0, 3.0);

        // Test get method
        assert_eq!(v.get(0).unwrap(), 1.0);
        assert_eq!(v.get(1).unwrap(), 2.0);
        assert_eq!(v.get(2).unwrap(), 3.0);

        // Test set method
        v.set(0, 10.0).unwrap();
        v.set(1, 20.0).unwrap();
        v.set(2, 30.0).unwrap();
        assert_eq!(v, Vector3::new(10.0, 20.0, 30.0));
    }

    #[test]
    fn test_get_error() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        let result = v.get(3);
        assert!(result.is_err());

        if let Err(err) = result {
            assert!(err.to_string().contains("index 3 out of bounds"));
        }
    }

    #[test]
    fn test_set_error() {
        let mut v = Vector3::new(1.0, 2.0, 3.0);
        let result = v.set(5, 42.0);
        assert!(result.is_err());

        if let Err(err) = result {
            assert!(err.to_string().contains("index 5 out of bounds"));
        }
    }

    #[test]
    fn test_normalize_zero_vector() {
        let zero = Vector3::zeros();
        let normalized = zero.normalize();
        assert_eq!(normalized, zero); // Zero vector normalizes to itself
    }

    #[test]
    fn test_to_array() {
        let v = Vector3::new(1.5, 2.5, 3.5);
        let arr = v.to_array();
        assert_eq!(arr, [1.5, 2.5, 3.5]);
    }

    #[test]
    fn test_div_assign_operator() {
        let mut v = Vector3::new(10.0, 20.0, 30.0);
        v /= 2.0;
        assert_eq!(v, Vector3::new(5.0, 10.0, 15.0));
    }

    #[test]
    fn test_indexing_operators() {
        let mut v = Vector3::new(1.0, 2.0, 3.0);

        // Test read indexing
        assert_eq!(v[0], 1.0);
        assert_eq!(v[1], 2.0);
        assert_eq!(v[2], 3.0);

        // Test write indexing
        v[0] = 10.0;
        v[1] = 20.0;
        v[2] = 30.0;
        assert_eq!(v, Vector3::new(10.0, 20.0, 30.0));
    }

    #[test]
    #[should_panic(expected = "Vector3 index out of bounds: 4")]
    fn test_index_panic() {
        let v = Vector3::new(1.0, 2.0, 3.0);
        let _ = v[4];
    }

    #[test]
    #[should_panic(expected = "Vector3 index out of bounds: 7")]
    fn test_index_mut_panic() {
        let mut v = Vector3::new(1.0, 2.0, 3.0);
        v[7] = 42.0;
    }

    #[test]
    fn test_display_formatting() {
        let v = Vector3::new(1.234567890, -2.345678901, 3.456789012);
        let display_output = format!("{}", v);

        // Check that it contains expected parts
        assert!(display_output.contains("Vector3("));
        assert!(display_output.contains("1.234567890"));
        assert!(display_output.contains("-2.345678901"));
        assert!(display_output.contains("3.456789012"));
        assert!(display_output.ends_with(")"));
    }

    #[test]
    fn test_to_spherical_north_pole() {
        let north_pole = Vector3::new(0.0, 0.0, 1.0);
        let (theta, phi) = north_pole.to_spherical();

        assert_eq!(theta, 0.0);
        assert_eq!(phi, crate::constants::HALF_PI);
    }

    #[test]
    fn test_to_spherical_south_pole() {
        let south_pole = Vector3::new(0.0, 0.0, -1.0);
        let (theta, phi) = south_pole.to_spherical();

        assert_eq!(theta, 0.0);
        assert_eq!(phi, -crate::constants::HALF_PI);
    }

    #[test]
    fn test_to_spherical_zero_z() {
        let on_equator = Vector3::new(1.0, 0.0, 0.0);
        let (theta, phi) = on_equator.to_spherical();

        assert_eq!(theta, 0.0);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_to_spherical_zero_vector() {
        let zero = Vector3::zeros();
        let (theta, phi) = zero.to_spherical();

        assert_eq!(theta, 0.0);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_spherical_roundtrip_at_poles() {
        let north_pole = Vector3::new(0.0, 0.0, 1.0);
        let (theta, phi) = north_pole.to_spherical();
        let roundtrip = Vector3::from_spherical(theta, phi);

        assert_eq!(roundtrip.z, north_pole.z);
        assert!(roundtrip.x.abs() < 1e-15, "x component: {}", roundtrip.x);
        assert!(roundtrip.y.abs() < 1e-15, "y component: {}", roundtrip.y);

        let south_pole = Vector3::new(0.0, 0.0, -1.0);
        let (theta, phi) = south_pole.to_spherical();
        let roundtrip = Vector3::from_spherical(theta, phi);

        assert_eq!(roundtrip.z, south_pole.z);
        assert!(roundtrip.x.abs() < 1e-15, "x component: {}", roundtrip.x);
        assert!(roundtrip.y.abs() < 1e-15, "y component: {}", roundtrip.y);
    }
}
