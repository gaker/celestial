//! 3x3 rotation matrices for astronomical coordinate transformations.
//!
//! Rotation matrices are the fundamental tool for transforming coordinates between
//! reference frames in astronomy. When you convert a star's position from ICRS to
//! galactic coordinates, or account for Earth's precession over centuries, or rotate
//! from equatorial to horizon coordinates for telescope pointing -- you're applying
//! rotation matrices.
//!
//! # The Role of Rotation Matrices in Astronomy
//!
//! A rotation matrix is a 3x3 orthogonal matrix with determinant +1. When applied to
//! a position vector, it rotates that vector while preserving its length. In astronomy,
//! we use rotation matrices for:
//!
//! - **Frame bias**: The small rotation between the FK5 catalog frame and ICRS
//! - **Precession**: Earth's axis traces a cone over ~26,000 years, requiring frame updates
//! - **Nutation**: Short-period oscillations of Earth's axis (18.6-year cycle and harmonics)
//! - **Earth rotation**: Converting between celestial and terrestrial reference frames
//! - **Coordinate system changes**: ICRS to galactic, equatorial to ecliptic, etc.
//!
//! # Composing Transformations
//!
//! Rotation matrices compose by multiplication. To apply rotation A, then rotation B,
//! you compute `B * A` (note the order -- the rightmost matrix acts first on the vector).
//!
//! ```
//! use celestial_core::constants::ARCSEC_TO_RAD;
//! use celestial_core::matrix::RotationMatrix3;
//!
//! // Build up a combined precession-nutation-bias transformation
//! let mut bias = RotationMatrix3::identity();
//! bias.rotate_x(-0.0068192 * ARCSEC_TO_RAD);  // Frame bias around X (arcsec → rad)
//! bias.rotate_z(0.041775 * ARCSEC_TO_RAD);    // Frame bias around Z (arcsec → rad)
//!
//! let mut precession = RotationMatrix3::identity();
//! precession.rotate_z(0.00385);  // Example precession angles in radians (~794″, ~423″)
//! precession.rotate_y(-0.00205);
//!
//! // Combined transformation: precession * bias
//! let combined = precession * bias;
//! ```
//!
//! For the full celestial-to-terrestrial transformation, the IAU defines the
//! complete chain as: `W * R * NPB` where NPB is the frame bias-precession-nutation
//! matrix, R is Earth rotation, and W is polar motion.
//!
//! # Rotation Conventions (ERFA-Compatible)
//!
//! This implementation follows the ERFA (Essential Routines for Fundamental Astronomy)
//! conventions. Rotations are defined as positive when counterclockwise when looking
//! from the positive axis toward the origin:
//!
//! - `rotate_x(phi)`: Rotation about the X-axis by angle phi (radians)
//! - `rotate_y(theta)`: Rotation about the Y-axis by angle theta (radians)
//! - `rotate_z(psi)`: Rotation about the Z-axis by angle psi (radians)
//!
//! This is the "passive" or "alias" convention where we rotate the coordinate frame
//! rather than the vector. A positive rotation of 90 degrees about Z takes the
//! vector `[1, 0, 0]` to `[0, -1, 0]`.
//!
//! # Storage Layout
//!
//! Elements are stored in row-major order as `[[f64; 3]; 3]`. The element at row `i`,
//! column `j` is `matrix.elements()[i][j]`. When the matrix multiplies a column
//! vector, the result is the standard matrix-vector product:
//!
//! ```text
//! | r00 r01 r02 |   | x |   | r00*x + r01*y + r02*z |
//! | r10 r11 r12 | * | y | = | r10*x + r11*y + r12*z |
//! | r20 r21 r22 |   | z |   | r20*x + r21*y + r22*z |
//! ```
//!
//! # Inverting Rotations
//!
//! For a proper rotation matrix, the inverse equals the transpose. This is much cheaper
//! to compute than a general matrix inverse and is numerically stable:
//!
//! ```
//! use celestial_core::matrix::RotationMatrix3;
//!
//! let mut m = RotationMatrix3::identity();
//! m.rotate_z(0.5);
//!
//! let m_inverse = m.transpose();
//!
//! // Verify: m * m_inverse should be identity
//! let product = m * m_inverse;
//! assert_eq!(product.elements()[0][0], 1.0);
//! ```
//!
//! # Spherical Coordinate Transformations
//!
//! For the common case of transforming right ascension and declination (or longitude
//! and latitude), use [`transform_spherical`](RotationMatrix3::transform_spherical):
//!
//! ```
//! use celestial_core::matrix::RotationMatrix3;
//! use celestial_core::constants::PI;
//!
//! let mut frame_rotation = RotationMatrix3::identity();
//! frame_rotation.rotate_z(PI / 6.0);  // 30 degree rotation
//!
//! let (ra, dec) = (0.0, 0.0);  // On the celestial equator at RA=0
//! let (new_ra, new_dec) = frame_rotation.transform_spherical(ra, dec);
//! ```

mod ops;
mod rotations;
#[cfg(test)]
mod tests;
mod validation;

use std::fmt;

#[cfg(feature = "serde")]
use crate::errors::{AstroError, AstroResult};

/// A 3x3 rotation matrix for coordinate frame transformations.
///
/// This type represents proper rotation matrices (orthogonal with determinant +1).
/// All angles are in radians. The matrix uses row-major storage.
///
/// # Construction
///
/// ```
/// use celestial_core::matrix::RotationMatrix3;
///
/// // Start with identity and build up rotations
/// let mut m = RotationMatrix3::identity();
/// m.rotate_z(0.1);  // Rotate 0.1 radians about Z
/// m.rotate_x(0.05); // Then rotate 0.05 radians about X
///
/// // Or construct directly from elements
/// let m = RotationMatrix3::from_array([
///     [1.0, 0.0, 0.0],
///     [0.0, 1.0, 0.0],
///     [0.0, 0.0, 1.0],
/// ])?;
/// # Ok::<(), celestial_core::errors::AstroError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(try_from = "UncheckedRotationMatrix3"))]
#[must_use]
pub struct RotationMatrix3 {
    elements: [[f64; 3]; 3],
}

impl RotationMatrix3 {
    /// Creates the 3x3 identity matrix.
    ///
    /// The identity matrix leaves any vector unchanged when applied. It serves as
    /// the starting point for building up rotation sequences.
    ///
    /// ```
    /// use celestial_core::matrix::RotationMatrix3;
    ///
    /// let m = RotationMatrix3::identity();
    /// let v = [1.0, 2.0, 3.0];
    /// let result = m.apply_to_vector(v);
    /// assert_eq!(result, v);
    /// ```
    #[inline]
    pub fn identity() -> Self {
        Self {
            elements: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    /// Returns a reference to the underlying 3x3 array.
    ///
    /// Useful when you need direct access to all elements, for example when
    /// passing to external APIs or serialization.
    #[inline]
    pub fn elements(&self) -> &[[f64; 3]; 3] {
        &self.elements
    }
}

#[cfg(feature = "serde")]
#[derive(serde::Deserialize)]
struct UncheckedRotationMatrix3 {
    elements: [[f64; 3]; 3],
}

#[cfg(feature = "serde")]
impl TryFrom<UncheckedRotationMatrix3> for RotationMatrix3 {
    type Error = AstroError;

    fn try_from(raw: UncheckedRotationMatrix3) -> AstroResult<Self> {
        Self::from_array(raw.elements)
    }
}

impl fmt::Display for RotationMatrix3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "RotationMatrix3:")?;
        for row in &self.elements {
            writeln!(f, "  [{:12.9} {:12.9} {:12.9}]", row[0], row[1], row[2])?;
        }
        Ok(())
    }
}
