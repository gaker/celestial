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
//! use celestial_core::RotationMatrix3;
//!
//! // Build up a combined precession-nutation-bias transformation
//! let mut bias = RotationMatrix3::identity();
//! bias.rotate_x(-0.0068192);  // Frame bias around X
//! bias.rotate_z(0.041775);    // Frame bias around Z
//!
//! let mut precession = RotationMatrix3::identity();
//! precession.rotate_z(0.00385);  // Example precession angles
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
//! column `j` is accessed as `matrix[(i, j)]` or `matrix.get(i, j)`. When the matrix
//! multiplies a column vector, the result is the standard matrix-vector product:
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
//! use celestial_core::RotationMatrix3;
//!
//! let mut m = RotationMatrix3::identity();
//! m.rotate_z(0.5);
//!
//! let m_inverse = m.transpose();
//!
//! // Verify: m * m_inverse should be identity
//! let product = m * m_inverse;
//! assert!((product.get(0, 0) - 1.0).abs() < 1e-15);
//! ```
//!
//! # Spherical Coordinate Transformations
//!
//! For the common case of transforming right ascension and declination (or longitude
//! and latitude), use [`transform_spherical`](RotationMatrix3::transform_spherical):
//!
//! ```
//! use celestial_core::RotationMatrix3;
//! use std::f64::consts::PI;
//!
//! let mut frame_rotation = RotationMatrix3::identity();
//! frame_rotation.rotate_z(PI / 6.0);  // 30 degree rotation
//!
//! let (ra, dec) = (0.0, 0.0);  // On the celestial equator at RA=0
//! let (new_ra, new_dec) = frame_rotation.transform_spherical(ra, dec);
//! ```

use std::fmt;

/// A 3x3 rotation matrix for coordinate frame transformations.
///
/// This type represents proper rotation matrices (orthogonal with determinant +1).
/// All angles are in radians. The matrix uses row-major storage.
///
/// # Construction
///
/// ```
/// use celestial_core::RotationMatrix3;
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
/// ]);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
    /// use celestial_core::RotationMatrix3;
    ///
    /// let m = RotationMatrix3::identity();
    /// let v = [1.0, 2.0, 3.0];
    /// let result = m.apply_to_vector(v);
    /// assert_eq!(result, v);
    /// ```
    pub fn identity() -> Self {
        Self {
            elements: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    /// Creates a rotation matrix from a 3x3 array of elements.
    ///
    /// The array is interpreted as row-major: `elements[i][j]` is row `i`, column `j`.
    ///
    /// This does not validate that the matrix is a proper rotation. Use
    /// [`is_rotation_matrix`](Self::is_rotation_matrix) to check if needed.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    ///
    /// // A rotation by 90 degrees about Z
    /// let m = RotationMatrix3::from_array([
    ///     [0.0, 1.0, 0.0],
    ///     [-1.0, 0.0, 0.0],
    ///     [0.0, 0.0, 1.0],
    /// ]);
    /// ```
    pub fn from_array(elements: [[f64; 3]; 3]) -> Self {
        Self { elements }
    }

    /// Returns the element at the specified row and column.
    ///
    /// Indices are 0-based. Panics if `row >= 3` or `col >= 3`.
    ///
    /// You can also use indexing syntax: `matrix[(row, col)]`.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.elements[row][col]
    }

    /// Sets the element at the specified row and column.
    ///
    /// Indices are 0-based. Panics if `row >= 3` or `col >= 3`.
    ///
    /// You can also use indexing syntax: `matrix[(row, col)] = value`.
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        self.elements[row][col] = value;
    }

    /// Returns a reference to the underlying 3x3 array.
    ///
    /// Useful when you need direct access to all elements, for example when
    /// passing to external APIs or serialization.
    pub fn elements(&self) -> &[[f64; 3]; 3] {
        &self.elements
    }

    /// Applies a rotation about the X-axis to this matrix (in place).
    ///
    /// The rotation angle `phi` is in radians. Positive angles rotate counterclockwise
    /// when looking from the positive X-axis toward the origin (ERFA convention).
    ///
    /// This modifies `self` to become `Rx(phi) * self`, where `Rx` is the standard
    /// X-axis rotation matrix:
    ///
    /// ```text
    /// Rx(phi) = | 1    0         0      |
    ///           | 0    cos(phi)  sin(phi)|
    ///           | 0   -sin(phi)  cos(phi)|
    /// ```
    ///
    /// In astronomy, X-axis rotations appear in frame bias corrections and some
    /// nutation models.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_x(FRAC_PI_2);  // 90 degrees
    ///
    /// // [0, 1, 0] rotates to [0, 0, -1]
    /// let v = m.apply_to_vector([0.0, 1.0, 0.0]);
    /// assert!(v[0].abs() < 1e-15);
    /// assert!(v[1].abs() < 1e-15);
    /// assert!((v[2] + 1.0).abs() < 1e-15);
    /// ```
    pub fn rotate_x(&mut self, phi: f64) {
        let (s, c) = phi.sin_cos();

        let a10 = c * self.elements[1][0] + s * self.elements[2][0];
        let a11 = c * self.elements[1][1] + s * self.elements[2][1];
        let a12 = c * self.elements[1][2] + s * self.elements[2][2];
        let a20 = -s * self.elements[1][0] + c * self.elements[2][0];
        let a21 = -s * self.elements[1][1] + c * self.elements[2][1];
        let a22 = -s * self.elements[1][2] + c * self.elements[2][2];

        self.elements[1][0] = a10;
        self.elements[1][1] = a11;
        self.elements[1][2] = a12;
        self.elements[2][0] = a20;
        self.elements[2][1] = a21;
        self.elements[2][2] = a22;
    }

    /// Applies a rotation about the Z-axis to this matrix (in place).
    ///
    /// The rotation angle `psi` is in radians. Positive angles rotate counterclockwise
    /// when looking from the positive Z-axis toward the origin (ERFA convention).
    ///
    /// This modifies `self` to become `Rz(psi) * self`, where `Rz` is the standard
    /// Z-axis rotation matrix:
    ///
    /// ```text
    /// Rz(psi) = | cos(psi)  sin(psi)  0 |
    ///           |-sin(psi)  cos(psi)  0 |
    ///           |    0         0      1 |
    /// ```
    ///
    /// Z-axis rotations are ubiquitous in astronomy. Earth rotation about its axis,
    /// precession in right ascension, and rotations in longitude all use Rz.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(FRAC_PI_2);  // 90 degrees
    ///
    /// // [1, 0, 0] rotates to [0, -1, 0]
    /// let v = m.apply_to_vector([1.0, 0.0, 0.0]);
    /// assert!(v[0].abs() < 1e-15);
    /// assert!((v[1] + 1.0).abs() < 1e-15);
    /// assert!(v[2].abs() < 1e-15);
    /// ```
    pub fn rotate_z(&mut self, psi: f64) {
        let (s, c) = psi.sin_cos();

        let a00 = c * self.elements[0][0] + s * self.elements[1][0];
        let a01 = c * self.elements[0][1] + s * self.elements[1][1];
        let a02 = c * self.elements[0][2] + s * self.elements[1][2];
        let a10 = -s * self.elements[0][0] + c * self.elements[1][0];
        let a11 = -s * self.elements[0][1] + c * self.elements[1][1];
        let a12 = -s * self.elements[0][2] + c * self.elements[1][2];

        self.elements[0][0] = a00;
        self.elements[0][1] = a01;
        self.elements[0][2] = a02;
        self.elements[1][0] = a10;
        self.elements[1][1] = a11;
        self.elements[1][2] = a12;
    }

    /// Applies a rotation about the Y-axis to this matrix (in place).
    ///
    /// The rotation angle `theta` is in radians. Positive angles rotate counterclockwise
    /// when looking from the positive Y-axis toward the origin (ERFA convention).
    ///
    /// This modifies `self` to become `Ry(theta) * self`, where `Ry` is the standard
    /// Y-axis rotation matrix:
    ///
    /// ```text
    /// Ry(theta) = | cos(theta)  0  -sin(theta) |
    ///             |     0       1       0      |
    ///             | sin(theta)  0   cos(theta) |
    /// ```
    ///
    /// Y-axis rotations appear in obliquity of the ecliptic and some precession models.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_y(FRAC_PI_2);  // 90 degrees
    ///
    /// // [0, 0, 1] rotates to [-1, 0, 0]
    /// let v = m.apply_to_vector([0.0, 0.0, 1.0]);
    /// assert!((v[0] + 1.0).abs() < 1e-15);
    /// assert!(v[1].abs() < 1e-15);
    /// assert!(v[2].abs() < 1e-15);
    /// ```
    pub fn rotate_y(&mut self, theta: f64) {
        let (s, c) = theta.sin_cos();

        let a00 = c * self.elements[0][0] - s * self.elements[2][0];
        let a01 = c * self.elements[0][1] - s * self.elements[2][1];
        let a02 = c * self.elements[0][2] - s * self.elements[2][2];
        let a20 = s * self.elements[0][0] + c * self.elements[2][0];
        let a21 = s * self.elements[0][1] + c * self.elements[2][1];
        let a22 = s * self.elements[0][2] + c * self.elements[2][2];

        self.elements[0][0] = a00;
        self.elements[0][1] = a01;
        self.elements[0][2] = a02;
        self.elements[2][0] = a20;
        self.elements[2][1] = a21;
        self.elements[2][2] = a22;
    }

    /// Multiplies this matrix by another, returning the product.
    ///
    /// Matrix multiplication is not commutative: `A * B` is generally different
    /// from `B * A`. The result represents the composition of two rotations where
    /// `other` is applied first, then `self`.
    ///
    /// You can also use the `*` operator: `a * b` or `&a * &b`.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    ///
    /// let mut rx = RotationMatrix3::identity();
    /// rx.rotate_x(0.1);
    ///
    /// let mut rz = RotationMatrix3::identity();
    /// rz.rotate_z(0.2);
    ///
    /// // First rotate by X, then by Z
    /// let combined = rz.multiply(&rx);
    /// // Equivalent using operator:
    /// let combined_op = rz * rx;
    /// ```
    pub fn multiply(&self, other: &Self) -> Self {
        let mut result = [[0.0; 3]; 3];

        for (i, row) in result.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                for k in 0..3 {
                    *cell += self.elements[i][k] * other.elements[k][j];
                }
            }
        }

        Self::from_array(result)
    }

    /// Applies this rotation matrix to a 3D vector.
    ///
    /// Computes the standard matrix-vector product `M * v`. For coordinate
    /// transformations, this rotates the position vector from one frame to another.
    ///
    /// You can also use the `*` operator with [`Vector3`](super::Vector3):
    /// `matrix * vector`.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(std::f64::consts::FRAC_PI_2);  // 90 degrees
    ///
    /// let v = [1.0, 0.0, 0.0];
    /// let rotated = m.apply_to_vector(v);
    /// // Result is approximately [0, -1, 0]
    /// ```
    pub fn apply_to_vector(&self, vector: [f64; 3]) -> [f64; 3] {
        [
            self.elements[0][0] * vector[0]
                + self.elements[0][1] * vector[1]
                + self.elements[0][2] * vector[2],
            self.elements[1][0] * vector[0]
                + self.elements[1][1] * vector[1]
                + self.elements[1][2] * vector[2],
            self.elements[2][0] * vector[0]
                + self.elements[2][1] * vector[1]
                + self.elements[2][2] * vector[2],
        ]
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

    /// Returns the transpose of this matrix.
    ///
    /// For a rotation matrix, the transpose equals the inverse. This provides
    /// a numerically stable way to compute the reverse transformation without
    /// general matrix inversion.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(0.5);
    /// m.rotate_x(0.3);
    ///
    /// let m_inv = m.transpose();
    ///
    /// // Applying m then m_inv returns to the original
    /// let v = [1.0, 2.0, 3.0];
    /// let rotated = m.apply_to_vector(v);
    /// let restored = m_inv.apply_to_vector(rotated);
    ///
    /// assert!((restored[0] - v[0]).abs() < 1e-14);
    /// assert!((restored[1] - v[1]).abs() < 1e-14);
    /// assert!((restored[2] - v[2]).abs() < 1e-14);
    /// ```
    pub fn transpose(&self) -> Self {
        Self::from_array([
            [
                self.elements[0][0],
                self.elements[1][0],
                self.elements[2][0],
            ],
            [
                self.elements[0][1],
                self.elements[1][1],
                self.elements[2][1],
            ],
            [
                self.elements[0][2],
                self.elements[1][2],
                self.elements[2][2],
            ],
        ])
    }

    /// Checks whether this matrix is a valid rotation matrix within a tolerance.
    ///
    /// A proper rotation matrix must satisfy two conditions:
    /// 1. Determinant equals +1 (not -1, which would be a reflection)
    /// 2. The matrix is orthogonal: `M * M^T = I`
    ///
    /// Due to floating-point arithmetic, these conditions are checked within
    /// the specified tolerance.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(0.5);
    /// m.rotate_x(0.3);
    /// assert!(m.is_rotation_matrix(1e-14));
    ///
    /// // A scaling matrix is not a rotation
    /// let scaled = RotationMatrix3::from_array([
    ///     [2.0, 0.0, 0.0],
    ///     [0.0, 1.0, 0.0],
    ///     [0.0, 0.0, 1.0],
    /// ]);
    /// assert!(!scaled.is_rotation_matrix(1e-14));
    /// ```
    pub fn is_rotation_matrix(&self, tolerance: f64) -> bool {
        let det = self.determinant();
        if (det - 1.0).abs() > tolerance {
            return false;
        }

        let rt = self.transpose();
        let product = self.multiply(&rt);
        let identity = Self::identity();

        for i in 0..3 {
            for j in 0..3 {
                if (product.elements[i][j] - identity.elements[i][j]).abs() > tolerance {
                    return false;
                }
            }
        }

        true
    }

    /// Returns the maximum absolute difference between corresponding elements.
    ///
    /// Useful for comparing matrices, especially when testing against reference
    /// implementations like ERFA.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    ///
    /// let a = RotationMatrix3::identity();
    /// let b = RotationMatrix3::from_array([
    ///     [1.0, 0.001, 0.0],
    ///     [0.0, 1.0, 0.0],
    ///     [0.0, 0.0, 1.0],
    /// ]);
    ///
    /// let diff = a.max_difference(&b);
    /// assert!((diff - 0.001).abs() < 1e-15);
    /// ```
    pub fn max_difference(&self, other: &Self) -> f64 {
        let mut max_diff: f64 = 0.0;

        for i in 0..3 {
            for j in 0..3 {
                let diff = (self.elements[i][j] - other.elements[i][j]).abs();
                max_diff = max_diff.max(diff);
            }
        }

        max_diff
    }

    /// Transforms spherical coordinates (RA, Dec or longitude, latitude) through
    /// this rotation matrix.
    ///
    /// This is the common operation for coordinate frame transformations in astronomy.
    /// The input angles are in radians:
    /// - `ra`: Right ascension or longitude (azimuthal angle from X toward Y)
    /// - `dec`: Declination or latitude (elevation from the XY plane)
    ///
    /// Internally, this converts to a unit Cartesian vector, applies the rotation,
    /// and converts back to spherical coordinates.
    ///
    /// The output RA is in the range `(-pi, pi]`. The output Dec is in `[-pi/2, pi/2]`.
    ///
    /// ```
    /// use celestial_core::RotationMatrix3;
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// // Rotate the equatorial coordinate system by 45 degrees in RA
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(FRAC_PI_4);
    ///
    /// let (ra, dec) = (0.0, 0.0);  // Point on equator at RA=0
    /// let (new_ra, new_dec) = m.transform_spherical(ra, dec);
    ///
    /// // RA shifted by -45 degrees (rotation convention), Dec unchanged
    /// assert!((new_ra + FRAC_PI_4).abs() < 1e-14);
    /// assert!(new_dec.abs() < 1e-14);
    /// ```
    pub fn transform_spherical(&self, ra: f64, dec: f64) -> (f64, f64) {
        let (sin_ra, cos_ra) = ra.sin_cos();
        let (sin_dec, cos_dec) = dec.sin_cos();
        let vector = [cos_dec * cos_ra, cos_dec * sin_ra, sin_dec];

        let transformed = self.apply_to_vector(vector);

        let new_ra = transformed[1].atan2(transformed[0]);
        let norm = (transformed[0] * transformed[0]
            + transformed[1] * transformed[1]
            + transformed[2] * transformed[2])
            .sqrt();
        let z = if norm == 0.0 {
            0.0
        } else {
            (transformed[2] / norm).clamp(-1.0, 1.0)
        };
        let new_dec = z.asin();

        (new_ra, new_dec)
    }
}

impl std::ops::Mul for RotationMatrix3 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        self.multiply(&rhs)
    }
}

impl std::ops::Mul<&RotationMatrix3> for RotationMatrix3 {
    type Output = RotationMatrix3;

    fn mul(self, rhs: &RotationMatrix3) -> RotationMatrix3 {
        self.multiply(rhs)
    }
}

impl std::ops::Mul<RotationMatrix3> for &RotationMatrix3 {
    type Output = RotationMatrix3;

    fn mul(self, rhs: RotationMatrix3) -> RotationMatrix3 {
        self.multiply(&rhs)
    }
}

impl std::ops::Mul<&RotationMatrix3> for &RotationMatrix3 {
    type Output = RotationMatrix3;

    fn mul(self, rhs: &RotationMatrix3) -> RotationMatrix3 {
        self.multiply(rhs)
    }
}

impl std::ops::Index<(usize, usize)> for RotationMatrix3 {
    type Output = f64;

    fn index(&self, (row, col): (usize, usize)) -> &f64 {
        &self.elements[row][col]
    }
}

impl std::ops::IndexMut<(usize, usize)> for RotationMatrix3 {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut f64 {
        &mut self.elements[row][col]
    }
}

impl std::ops::Mul<super::Vector3> for RotationMatrix3 {
    type Output = super::Vector3;

    fn mul(self, vec: super::Vector3) -> super::Vector3 {
        let result = self.apply_to_vector([vec.x, vec.y, vec.z]);
        super::Vector3::from_array(result)
    }
}

impl std::ops::Mul<super::Vector3> for &RotationMatrix3 {
    type Output = super::Vector3;

    fn mul(self, vec: super::Vector3) -> super::Vector3 {
        let result = self.apply_to_vector([vec.x, vec.y, vec.z]);
        super::Vector3::from_array(result)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::HALF_PI;

    #[test]
    fn test_identity_and_get() {
        let m = RotationMatrix3::identity();
        assert_eq!(m.get(0, 0), 1.0);
        assert_eq!(m.get(1, 1), 1.0);
        assert_eq!(m.get(2, 2), 1.0);
        assert_eq!(m.get(0, 1), 0.0);
    }

    #[test]
    fn test_set() {
        let mut m = RotationMatrix3::identity();
        m.set(0, 1, 0.5);
        assert_eq!(m.get(0, 1), 0.5);
    }

    #[test]
    fn test_rotate_z() {
        // ERFA convention: Rz(+psi) rotates anticlockwise looking from +z toward origin
        // This means [1,0,0] -> [cos(psi), -sin(psi), 0]
        // At 90°: [1,0,0] -> [0, -1, 0]
        let mut m = RotationMatrix3::identity();
        m.rotate_z(HALF_PI);
        let result = m.apply_to_vector([1.0, 0.0, 0.0]);
        assert!(result[0].abs() < 1e-15);
        assert!((result[1] + 1.0).abs() < 1e-15);
        assert!(result[2].abs() < 1e-15);
    }

    #[test]
    fn test_rotate_x() {
        // ERFA convention: Rx(+phi) rotates anticlockwise looking from +x toward origin
        // This means [0,1,0] -> [0, cos(phi), -sin(phi)]
        // At 90°: [0,1,0] -> [0, 0, -1]
        let mut m = RotationMatrix3::identity();
        m.rotate_x(HALF_PI);
        let result = m.apply_to_vector([0.0, 1.0, 0.0]);
        assert!(result[0].abs() < 1e-15);
        assert!(result[1].abs() < 1e-15);
        assert!((result[2] + 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_rotate_y() {
        // ERFA convention: Ry(+theta) rotates anticlockwise looking from +y toward origin
        // This means [0,0,1] -> [-sin(theta), 0, cos(theta)]
        // At 90°: [0,0,1] -> [-1, 0, 0]
        let mut m = RotationMatrix3::identity();
        m.rotate_y(HALF_PI);
        let result = m.apply_to_vector([0.0, 0.0, 1.0]);
        assert!((result[0] + 1.0).abs() < 1e-15);
        assert!(result[1].abs() < 1e-15);
        assert!(result[2].abs() < 1e-15);
    }

    #[test]
    fn test_is_rotation_matrix_valid() {
        let mut m = RotationMatrix3::identity();
        m.rotate_z(0.5);
        assert!(m.is_rotation_matrix(1e-14));
    }

    #[test]
    fn test_is_rotation_matrix_bad_determinant() {
        let m = RotationMatrix3::from_array([[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert!(!m.is_rotation_matrix(1e-15));
    }

    #[test]
    fn test_is_rotation_matrix_not_orthogonal() {
        let m = RotationMatrix3::from_array([[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert!(!m.is_rotation_matrix(1e-15));
    }

    #[test]
    fn test_transform_spherical_identity() {
        let m = RotationMatrix3::identity();
        let (ra, dec) = (1.0, 0.5);
        let (new_ra, new_dec) = m.transform_spherical(ra, dec);
        assert!((new_ra - ra).abs() < 1e-14);
        assert!((new_dec - dec).abs() < 1e-14);
    }

    #[test]
    fn test_transform_spherical_rotation() {
        // ERFA Rz rotates in opposite direction to naive expectation
        // Rz(+90°) takes RA=0 to RA=-90° (or equivalently RA=270°=-HALF_PI)
        let mut m = RotationMatrix3::identity();
        m.rotate_z(HALF_PI);
        let (new_ra, new_dec) = m.transform_spherical(0.0, 0.0);
        assert!((new_ra + HALF_PI).abs() < 1e-14);
        assert!(new_dec.abs() < 1e-14);
    }

    #[test]
    fn test_transform_spherical_zero_norm() {
        let zero_matrix =
            RotationMatrix3::from_array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        let (_, dec) = zero_matrix.transform_spherical(0.0, 0.0);
        assert!(dec.is_finite());
    }

    #[test]
    fn test_mul_matrix_matrix() {
        let mut a = RotationMatrix3::identity();
        a.rotate_x(0.1);
        let mut b = RotationMatrix3::identity();
        b.rotate_y(0.2);

        let r1 = a * b;
        let r2 = a * &b;
        let r3 = &a * b;
        let r4 = &a * &b;

        assert_eq!(r1.get(0, 0), r2.get(0, 0));
        assert_eq!(r2.get(0, 0), r3.get(0, 0));
        assert_eq!(r3.get(0, 0), r4.get(0, 0));
    }

    #[test]
    fn test_index_operators() {
        let mut m = RotationMatrix3::identity();
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 0.0);
        m[(0, 1)] = 0.5;
        assert_eq!(m[(0, 1)], 0.5);
    }

    #[test]
    fn test_mul_matrix_vector() {
        use crate::Vector3;
        let m = RotationMatrix3::identity();
        let v = Vector3::new(1.0, 2.0, 3.0);
        let r1 = m * v;
        let r2 = &m * v;
        assert_eq!(r1, v);
        assert_eq!(r2, v);
    }

    #[test]
    fn test_display() {
        let mut m = RotationMatrix3::identity();
        m.rotate_z(0.1);
        let s = format!("{}", m);
        assert!(s.contains("RotationMatrix3:"));
        assert!(s.contains("["));
    }

    #[test]
    fn test_max_difference() {
        let a = RotationMatrix3::identity();
        let b = RotationMatrix3::from_array([[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert!((a.max_difference(&b) - 0.1).abs() < 1e-15);
    }

    #[test]
    fn test_elements() {
        let m = RotationMatrix3::identity();
        let e = m.elements();
        assert_eq!(e[0][0], 1.0);
        assert_eq!(e[1][1], 1.0);
    }
}
