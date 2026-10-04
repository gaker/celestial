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
//! use celestial_core::matrix::Vector3;
//! use celestial_core::constants::QUARTER_PI;
//!
//! // A star at RA=45°, Dec=30° (in radians)
//! let ra = QUARTER_PI;
//! let dec = QUARTER_PI / 1.5;  // ~30°
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
//! use celestial_core::matrix::Vector3;
//!
//! let v = Vector3::new(3.0, 4.0, 0.0);
//! let unit = v.normalize()?;
//! assert_eq!(unit.magnitude(), 1.0);
//! # Ok::<(), celestial_core::errors::AstroError>(())
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
//! use celestial_core::matrix::Vector3;
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
mod geometry;
mod ops;
#[cfg(test)]
mod tests;

use crate::errors::{AstroError, AstroResult, MathErrorKind};
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
/// use celestial_core::matrix::Vector3;
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
#[must_use]
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
    /// Returns an error for indices outside 0-2.
    pub fn get(&self, index: usize) -> AstroResult<f64> {
        match index {
            0 => Ok(self.x),
            1 => Ok(self.y),
            2 => Ok(self.z),
            _ => Err(index_error("Vector3::get", index)),
        }
    }

    /// Sets the component at the given index (0=x, 1=y, 2=z).
    ///
    /// Returns an error for indices outside 0-2.
    pub fn set(&mut self, index: usize, value: f64) -> AstroResult<()> {
        let component = match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => return Err(index_error("Vector3::set", index)),
        };
        *component = value;
        Ok(())
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
}

fn index_error(operation: &str, index: usize) -> AstroError {
    AstroError::math_error(
        operation,
        MathErrorKind::InvalidInput,
        &format!("index {} out of bounds (valid range: 0-2)", index),
    )
}

impl fmt::Display for Vector3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Vector3({:.9}, {:.9}, {:.9})", self.x, self.y, self.z)
    }
}
