use super::Vector3;
use crate::errors::{AstroError, AstroResult, MathErrorKind};

impl Vector3 {
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
    /// Returns an error if the vector's length is zero or not finite.
    ///
    /// ```
    /// use celestial_core::matrix::Vector3;
    ///
    /// let v = Vector3::new(3.0, 4.0, 0.0);
    /// let unit = v.normalize()?;
    /// assert_eq!(unit.magnitude(), 1.0);
    /// assert_eq!(unit, Vector3::new(0.6000000000000001, 0.8, 0.0));
    /// # Ok::<(), celestial_core::errors::AstroError>(())
    /// ```
    #[inline]
    pub fn normalize(&self) -> AstroResult<Self> {
        let mag = self.magnitude();
        if mag == 0.0 {
            return Err(AstroError::math_error(
                "Vector3::normalize",
                MathErrorKind::DivisionByZero,
                "Cannot normalize a zero-length vector",
            ));
        }
        if !mag.is_finite() {
            return Err(AstroError::math_error(
                "Vector3::normalize",
                MathErrorKind::NotFinite,
                "Vector length is not finite",
            ));
        }
        Ok(*self * (1.0 / mag))
    }

    /// Computes the dot product (inner product) with another vector.
    ///
    /// For unit vectors, this equals the cosine of the angle between them:
    /// `a.dot(&b) = cos(θ)`. Use this to compute angular separation between
    /// celestial positions.
    ///
    /// ```
    /// use celestial_core::matrix::Vector3;
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
    /// use celestial_core::matrix::Vector3;
    ///
    /// let x = Vector3::x_axis();
    /// let y = Vector3::y_axis();
    /// let z = x.cross(&y);
    /// assert_eq!(z, Vector3::z_axis());  // X × Y = Z
    /// ```
    #[inline]
    pub fn cross(&self, other: &Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Creates a unit vector from spherical coordinates.
    ///
    /// - `ra`: Azimuthal angle from +X toward +Y (right ascension), in radians
    /// - `dec`: Elevation from XY plane (declination), in radians
    ///
    /// The result is always a unit vector (magnitude = 1).
    ///
    /// ```
    /// use celestial_core::matrix::Vector3;
    /// use celestial_core::constants::HALF_PI;
    ///
    /// // RA=0, Dec=0 → points along +X
    /// let v = Vector3::from_spherical(0.0, 0.0);
    /// assert_eq!(v, Vector3::new(1.0, 0.0, 0.0));
    ///
    /// // RA=90°, Dec=0 → points along +Y. f64 π/2 sits just below the true value,
    /// // so its cosine leaves a tiny x component.
    /// let v = Vector3::from_spherical(HALF_PI, 0.0);
    /// assert_eq!(v, Vector3::new(6.123233995736766e-17, 1.0, 0.0));
    ///
    /// // RA=0, Dec=90° → points along +Z (north pole)
    /// let v = Vector3::from_spherical(0.0, HALF_PI);
    /// assert_eq!(v, Vector3::new(6.123233995736766e-17, 0.0, 1.0));
    /// ```
    #[inline]
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
    /// use celestial_core::matrix::Vector3;
    /// use celestial_core::constants::HALF_PI;
    ///
    /// let v = Vector3::new(0.0, 0.0, 1.0);  // North pole
    /// let (theta, phi) = v.to_spherical();
    /// assert_eq!(theta, 0.0);
    /// assert_eq!(phi, HALF_PI);
    /// ```
    #[inline]
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
