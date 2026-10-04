use super::RotationMatrix3;
use crate::matrix::Vector3;

impl RotationMatrix3 {
    /// Multiplies this matrix by another, returning the product.
    ///
    /// Matrix multiplication is not commutative: `A * B` is generally different
    /// from `B * A`. The result represents the composition of two rotations where
    /// `other` is applied first, then `self`.
    ///
    /// You can also use the `*` operator: `a * b` or `&a * &b`.
    ///
    /// ```
    /// use celestial_core::matrix::RotationMatrix3;
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
    #[inline]
    pub fn multiply(&self, other: &Self) -> Self {
        let mut result = [[0.0; 3]; 3];

        for (i, row) in result.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                for k in 0..3 {
                    *cell += self.elements[i][k] * other.elements[k][j];
                }
            }
        }

        Self { elements: result }
    }

    /// Applies this rotation matrix to a 3D vector.
    ///
    /// Computes the standard matrix-vector product `M * v`. For coordinate
    /// transformations, this rotates the position vector from one frame to another.
    ///
    /// You can also use the `*` operator with [`Vector3`](crate::matrix::Vector3):
    /// `matrix * vector`.
    ///
    /// ```
    /// use celestial_core::matrix::RotationMatrix3;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(celestial_core::constants::HALF_PI);  // 90 degrees
    ///
    /// let v = [1.0, 0.0, 0.0];
    /// let rotated = m.apply_to_vector(v);
    /// // Result is approximately [0, -1, 0]
    /// ```
    #[inline]
    pub fn apply_to_vector(&self, vector: [f64; 3]) -> [f64; 3] {
        self.elements
            .map(|row| row.iter().zip(vector).fold(0.0, |sum, (m, v)| sum + m * v))
    }

    /// Returns the transpose of this matrix.
    ///
    /// For a rotation matrix, the transpose equals the inverse. This provides
    /// a numerically stable way to compute the reverse transformation without
    /// general matrix inversion.
    ///
    /// ```
    /// use celestial_core::matrix::RotationMatrix3;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(0.5);
    /// m.rotate_x(0.3);
    ///
    /// let m_inv = m.transpose();
    ///
    /// // Applying m then m_inv returns to the original, up to rounding in the last bit
    /// let v = [1.0, 2.0, 3.0];
    /// let rotated = m.apply_to_vector(v);
    /// let restored = m_inv.apply_to_vector(rotated);
    ///
    /// assert_eq!(restored, [1.0, 2.0, 2.9999999999999996]);
    /// ```
    #[inline]
    pub fn transpose(&self) -> Self {
        let elements = [
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
        ];
        Self { elements }
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
    /// use celestial_core::matrix::RotationMatrix3;
    /// use celestial_core::constants::QUARTER_PI;
    ///
    /// // Rotate the equatorial coordinate system by 45 degrees in RA
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(QUARTER_PI);
    ///
    /// let (ra, dec) = (0.0, 0.0);  // Point on equator at RA=0
    /// let (new_ra, new_dec) = m.transform_spherical(ra, dec);
    ///
    /// // RA shifted by -45 degrees (rotation convention), Dec unchanged
    /// assert_eq!(new_ra, -QUARTER_PI);
    /// assert_eq!(new_dec, 0.0);
    /// ```
    #[inline]
    pub fn transform_spherical(&self, ra: f64, dec: f64) -> (f64, f64) {
        (self * Vector3::from_spherical(ra, dec)).to_spherical()
    }
}

impl std::ops::Mul for RotationMatrix3 {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        self.multiply(&rhs)
    }
}

impl std::ops::Mul<&Self> for RotationMatrix3 {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: &Self) -> Self {
        self.multiply(rhs)
    }
}

impl std::ops::Mul<RotationMatrix3> for &RotationMatrix3 {
    type Output = RotationMatrix3;

    #[inline]
    fn mul(self, rhs: RotationMatrix3) -> RotationMatrix3 {
        self.multiply(&rhs)
    }
}

impl std::ops::Mul<&RotationMatrix3> for &RotationMatrix3 {
    type Output = RotationMatrix3;

    #[inline]
    fn mul(self, rhs: &RotationMatrix3) -> RotationMatrix3 {
        self.multiply(rhs)
    }
}

impl std::ops::Mul<Vector3> for RotationMatrix3 {
    type Output = Vector3;

    #[inline]
    fn mul(self, vec: Vector3) -> Vector3 {
        let result = self.apply_to_vector([vec.x, vec.y, vec.z]);
        Vector3::from_array(result)
    }
}

impl std::ops::Mul<Vector3> for &RotationMatrix3 {
    type Output = Vector3;

    #[inline]
    fn mul(self, vec: Vector3) -> Vector3 {
        let result = self.apply_to_vector([vec.x, vec.y, vec.z]);
        Vector3::from_array(result)
    }
}
