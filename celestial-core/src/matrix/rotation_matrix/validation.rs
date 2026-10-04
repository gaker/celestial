use super::RotationMatrix3;
use crate::errors::{AstroError, AstroResult, MathErrorKind};

// Rotations built in double precision are orthonormal to about 1e-16. 1e-12 still
// accepts matrices copied out to 13 significant digits, and rejects any scale or
// shear that would move a direction by more than about 0.2 microarcseconds.
const ROTATION_TOLERANCE: f64 = 1e-12;

impl RotationMatrix3 {
    /// Creates a rotation matrix from a 3x3 array of elements.
    ///
    /// The array is interpreted as row-major: `elements[i][j]` is row `i`, column `j`.
    ///
    /// Returns an error unless every element is finite and the matrix is a proper
    /// rotation: `M * M^T` within 1e-12 of identity and determinant within 1e-12 of +1.
    ///
    /// ```
    /// use celestial_core::matrix::RotationMatrix3;
    ///
    /// // A rotation by 90 degrees about Z
    /// let m = RotationMatrix3::from_array([
    ///     [0.0, 1.0, 0.0],
    ///     [-1.0, 0.0, 0.0],
    ///     [0.0, 0.0, 1.0],
    /// ])?;
    ///
    /// // A scaling matrix is not a rotation
    /// let scaled = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    /// assert!(RotationMatrix3::from_array(scaled).is_err());
    /// # Ok::<(), celestial_core::errors::AstroError>(())
    /// ```
    pub fn from_array(elements: [[f64; 3]; 3]) -> AstroResult<Self> {
        if !elements.iter().flatten().all(|e| e.is_finite()) {
            return Err(AstroError::math_error(
                "RotationMatrix3::from_array",
                MathErrorKind::NotFinite,
                "Matrix elements must be finite",
            ));
        }
        let matrix = Self { elements };
        if !matrix.is_rotation_matrix(ROTATION_TOLERANCE) {
            return Err(AstroError::math_error(
                "RotationMatrix3::from_array",
                MathErrorKind::InvalidInput,
                "Matrix is not a proper rotation",
            ));
        }
        Ok(matrix)
    }

    /// Computes the determinant of this matrix.
    ///
    /// For a proper rotation matrix, the determinant is always +1. A determinant
    /// of -1 indicates a reflection (improper rotation). Values far from +/-1
    /// indicate the matrix is not orthogonal.
    ///
    /// Used internally by [`is_rotation_matrix`](Self::is_rotation_matrix).
    pub fn determinant(&self) -> f64 {
        let m = &self.elements;

        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    /// Checks whether this matrix is a valid rotation matrix within a tolerance.
    ///
    /// A proper rotation matrix must satisfy two conditions:
    /// 1. Determinant equals +1 (not -1, which would be a reflection)
    /// 2. The matrix is orthogonal: `M * M^T = I`
    ///
    /// Due to floating-point arithmetic, these conditions are checked within
    /// the specified tolerance. [`from_array`](Self::from_array) runs this check
    /// at 1e-12; use this method to check a tighter bound.
    ///
    /// ```
    /// use celestial_core::matrix::RotationMatrix3;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(0.5);
    /// m.rotate_x(0.3);
    /// assert!(m.is_rotation_matrix(1e-14));
    /// ```
    pub fn is_rotation_matrix(&self, tolerance: f64) -> bool {
        let unit_determinant = libm::fabs(self.determinant() - 1.0) <= tolerance;
        let orthogonality_error = self
            .multiply(&self.transpose())
            .max_difference(&Self::identity());
        unit_determinant && orthogonality_error <= tolerance
    }

    /// Returns the maximum absolute difference between corresponding elements.
    ///
    /// Useful for comparing matrices, especially when testing against reference
    /// implementations like ERFA.
    ///
    /// ```
    /// use celestial_core::matrix::RotationMatrix3;
    ///
    /// let a = RotationMatrix3::identity();
    /// let mut b = RotationMatrix3::identity();
    /// b.rotate_z(0.001);
    ///
    /// // The largest change is sin(0.001) in the off-diagonal elements
    /// let diff = a.max_difference(&b);
    /// assert_eq!(diff, 0.0009999998333333417);
    /// ```
    pub fn max_difference(&self, other: &Self) -> f64 {
        let mut max_diff: f64 = 0.0;

        for i in 0..3 {
            for j in 0..3 {
                let diff = libm::fabs(self.elements[i][j] - other.elements[i][j]);
                if diff.is_nan() || diff > max_diff {
                    max_diff = diff;
                }
            }
        }

        max_diff
    }
}
