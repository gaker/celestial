use super::RotationMatrix3;

impl RotationMatrix3 {
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
    /// use celestial_core::matrix::RotationMatrix3;
    /// use celestial_core::constants::HALF_PI;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_x(HALF_PI);  // 90 degrees
    ///
    /// // [0, 1, 0] rotates to [0, 0, -1]. f64 π/2 sits just below the true value,
    /// // so its cosine leaves a tiny y component.
    /// let v = m.apply_to_vector([0.0, 1.0, 0.0]);
    /// assert_eq!(v, [0.0, 6.123233995736766e-17, -1.0]);
    /// ```
    #[inline]
    pub fn rotate_x(&mut self, phi: f64) {
        let (s, c) = libm::sincos(phi);

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
    /// use celestial_core::matrix::RotationMatrix3;
    /// use celestial_core::constants::HALF_PI;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_z(HALF_PI);  // 90 degrees
    ///
    /// // [1, 0, 0] rotates to [0, -1, 0]. f64 π/2 sits just below the true value,
    /// // so its cosine leaves a tiny x component.
    /// let v = m.apply_to_vector([1.0, 0.0, 0.0]);
    /// assert_eq!(v, [6.123233995736766e-17, -1.0, 0.0]);
    /// ```
    #[inline]
    pub fn rotate_z(&mut self, psi: f64) {
        let (s, c) = libm::sincos(psi);

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
    /// use celestial_core::matrix::RotationMatrix3;
    /// use celestial_core::constants::HALF_PI;
    ///
    /// let mut m = RotationMatrix3::identity();
    /// m.rotate_y(HALF_PI);  // 90 degrees
    ///
    /// // [0, 0, 1] rotates to [-1, 0, 0]. f64 π/2 sits just below the true value,
    /// // so its cosine leaves a tiny z component.
    /// let v = m.apply_to_vector([0.0, 0.0, 1.0]);
    /// assert_eq!(v, [-1.0, 0.0, 6.123233995736766e-17]);
    /// ```
    #[inline]
    pub fn rotate_y(&mut self, theta: f64) {
        let (s, c) = libm::sincos(theta);

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
}
