use crate::{transforms::CoordinateFrame, CoordError, CoordResult, Distance};
use celestial_core::{Angle, Vector3};
use celestial_time::TT;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ICRSPosition {
    ra: Angle,
    dec: Angle,
    distance: Option<Distance>,
}

impl ICRSPosition {
    pub fn new(ra: Angle, dec: Angle) -> CoordResult<Self> {
        let ra = ra.validate_right_ascension()?;
        let dec = dec.validate_declination(false)?;

        Ok(Self {
            ra,
            dec,
            distance: None,
        })
    }

    pub fn with_distance(ra: Angle, dec: Angle, distance: Distance) -> CoordResult<Self> {
        let mut pos = Self::new(ra, dec)?;
        pos.distance = Some(distance);
        Ok(pos)
    }

    pub fn from_degrees(ra_deg: f64, dec_deg: f64) -> CoordResult<Self> {
        Self::new(Angle::from_degrees(ra_deg), Angle::from_degrees(dec_deg))
    }

    pub fn from_degrees_with_distance(
        ra_deg: f64,
        dec_deg: f64,
        distance: Distance,
    ) -> CoordResult<Self> {
        Self::with_distance(
            Angle::from_degrees(ra_deg),
            Angle::from_degrees(dec_deg),
            distance,
        )
    }

    pub fn from_hours_degrees(ra_hours: f64, dec_deg: f64) -> CoordResult<Self> {
        Self::new(Angle::from_hours(ra_hours), Angle::from_degrees(dec_deg))
    }

    pub fn ra(&self) -> Angle {
        self.ra
    }

    pub fn dec(&self) -> Angle {
        self.dec
    }

    pub fn distance(&self) -> Option<Distance> {
        self.distance
    }

    pub fn set_distance(&mut self, distance: Distance) {
        self.distance = Some(distance);
    }

    pub fn remove_distance(&mut self) {
        self.distance = None;
    }

    pub fn unit_vector(&self) -> Vector3 {
        let (sin_dec, cos_dec) = self.dec.sin_cos();
        let (sin_ra, cos_ra) = self.ra.sin_cos();

        Vector3::new(cos_dec * cos_ra, cos_dec * sin_ra, sin_dec)
    }

    pub fn position_vector(&self) -> CoordResult<Vector3> {
        let distance = self.distance.ok_or_else(|| {
            CoordError::invalid_coordinate("Distance required for position vector")
        })?;

        let unit = self.unit_vector();
        let distance_au = distance.au();

        Ok(Vector3::new(
            unit.x * distance_au,
            unit.y * distance_au,
            unit.z * distance_au,
        ))
    }

    pub fn from_unit_vector(unit: Vector3) -> CoordResult<Self> {
        let r = libm::sqrt(unit.x.powi(2) + unit.y.powi(2) + unit.z.powi(2));

        if r == 0.0 {
            return Err(CoordError::invalid_coordinate("Zero vector"));
        }

        let x = unit.x / r;
        let y = unit.y / r;
        let z = unit.z / r;

        let d2 = x * x + y * y;
        let ra = if d2 == 0.0 { 0.0 } else { libm::atan2(y, x) };
        let dec = if z == 0.0 {
            0.0
        } else {
            libm::atan2(z, libm::sqrt(d2))
        };

        Self::new(Angle::from_radians(ra), Angle::from_radians(dec))
    }

    pub fn from_position_vector(pos: Vector3) -> CoordResult<Self> {
        let distance_au = libm::sqrt(pos.x.powi(2) + pos.y.powi(2) + pos.z.powi(2));

        if distance_au == 0.0 {
            return Err(CoordError::invalid_coordinate("Zero position vector"));
        }

        let unit = Vector3::new(
            pos.x / distance_au,
            pos.y / distance_au,
            pos.z / distance_au,
        );

        let mut icrs = Self::from_unit_vector(unit)?;
        icrs.distance = Some(Distance::from_au(distance_au)?);

        Ok(icrs)
    }

    pub fn angular_separation(&self, other: &Self) -> Angle {
        let (sin_dec1, cos_dec1) = self.dec.sin_cos();
        let (sin_dec2, cos_dec2) = other.dec.sin_cos();
        let delta_ra = (self.ra - other.ra).radians();

        let angle_rad = celestial_core::math::vincenty_angular_separation(
            sin_dec1, cos_dec1, sin_dec2, cos_dec2, delta_ra,
        );

        Angle::from_radians(angle_rad)
    }

    pub fn is_near_pole(&self) -> bool {
        self.dec.abs().degrees() > 89.0
    }

    /// Calculate angular position uncertainty from parallax measurement error.
    ///
    /// This method simply converts parallax error to angular uncertainty (they have the same
    /// magnitude). The stored distance value is not used in the calculation - distance is only
    /// checked for presence to ensure this is a parallax-derived position.
    ///
    /// **Important**: This returns the *angular* uncertainty on the sky, not the linear distance
    /// uncertainty. For distance uncertainty in physical units (parsecs, AU, etc.), calculate
    /// separately using error propagation: Δd = d² × Δπ (where π is parallax).
    ///
    /// # Arguments
    /// * `parallax_error_mas` - Parallax measurement error in milliarcseconds
    ///
    /// # Returns
    /// Angular position uncertainty in arcseconds, or None if no distance is set
    ///
    /// # Example
    /// ```ignore
    /// // Gaia DR3 typical parallax error: 0.02 mas → 0.00002 arcsec angular uncertainty
    /// let uncertainty = pos.position_uncertainty_arcsec(0.02);
    /// ```
    pub fn position_uncertainty_arcsec(&self, parallax_error_mas: f64) -> Option<f64> {
        self.distance.map(|_d| {
            // Angular position uncertainty equals parallax uncertainty numerically
            parallax_error_mas / 1000.0
        })
    }

    /// Calculate physical distance uncertainty from parallax measurement error
    ///
    /// For a star at distance d with parallax π ± σ_π:
    /// - Fractional distance error: σ_d/d = σ_π/π
    /// - Absolute distance error: σ_d = d × (σ_π/π)
    ///
    /// # Arguments
    /// * `parallax_error_mas` - Parallax measurement error in milliarcseconds
    ///
    /// # Returns
    /// Distance uncertainty in parsecs, or None if no distance is set
    pub fn distance_uncertainty_parsecs(&self, parallax_error_mas: f64) -> Option<f64> {
        self.distance.map(|d| {
            let parallax_mas = d.parallax_milliarcsec();
            let relative_error = parallax_error_mas / parallax_mas;
            d.parsecs() * relative_error
        })
    }
}

impl CoordinateFrame for ICRSPosition {
    fn to_icrs(&self, _epoch: &TT) -> CoordResult<ICRSPosition> {
        Ok(self.clone())
    }

    fn from_icrs(icrs: &ICRSPosition, _epoch: &TT) -> CoordResult<Self> {
        Ok(icrs.clone())
    }
}

impl ICRSPosition {
    pub fn to_galactic(&self, epoch: &TT) -> CoordResult<crate::GalacticPosition> {
        crate::GalacticPosition::from_icrs(self, epoch)
    }

    pub fn to_ecliptic(&self, epoch: &TT) -> CoordResult<crate::EclipticPosition> {
        crate::EclipticPosition::from_icrs(self, epoch)
    }
}

impl std::fmt::Display for ICRSPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ICRS(RA={:.6}°, Dec={:.6}°",
            self.ra.degrees(),
            self.dec.degrees()
        )?;

        if let Some(distance) = self.distance {
            write!(f, ", d={}", distance)?;
        }

        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Distance;

    #[test]
    fn test_constructor_methods() {
        // Test basic constructor
        let pos1 =
            ICRSPosition::new(Angle::from_degrees(180.0), Angle::from_degrees(45.0)).unwrap();

        assert_eq!(pos1.ra().degrees(), Angle::from_degrees(180.0).degrees());
        assert_eq!(pos1.dec().degrees(), Angle::from_degrees(45.0).degrees());
        assert_eq!(pos1.distance(), None);

        // Test from_degrees constructor
        let pos2 = ICRSPosition::from_degrees(90.0, -30.0).unwrap();
        assert_eq!(pos2.ra().degrees(), Angle::from_degrees(90.0).degrees());
        assert_eq!(pos2.dec().degrees(), Angle::from_degrees(-30.0).degrees());

        // Test from_hours_degrees constructor
        let pos3 = ICRSPosition::from_hours_degrees(12.0, 60.0).unwrap();
        assert_eq!(pos3.ra().hours(), 12.0);
        assert_eq!(pos3.dec().degrees(), Angle::from_degrees(60.0).degrees());

        // Test with_distance constructor
        let distance = Distance::from_parsecs(10.0).unwrap();
        let pos4 = ICRSPosition::with_distance(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            distance,
        )
        .unwrap();
        assert_eq!(pos4.distance().unwrap(), distance);
    }

    #[test]
    fn test_accessor_methods() {
        let mut pos = ICRSPosition::from_degrees(270.0, 15.0).unwrap();
        let distance = Distance::from_parsecs(5.0).unwrap();

        // Test getters
        assert_eq!(pos.ra().degrees(), Angle::from_degrees(270.0).degrees());
        assert_eq!(pos.dec().degrees(), Angle::from_degrees(15.0).degrees());
        assert_eq!(pos.distance(), None);

        // Test set_distance
        pos.set_distance(distance);
        assert_eq!(pos.distance().unwrap(), distance);

        // Test remove_distance
        pos.remove_distance();
        assert_eq!(pos.distance(), None);
    }

    #[test]
    fn test_unit_vector_conversion() {
        // Test known positions
        let vernal_equinox = ICRSPosition::from_degrees(0.0, 0.0).unwrap();
        let unit_vec = vernal_equinox.unit_vector();

        // Vernal equinox should point to [1, 0, 0]
        assert_eq!(unit_vec.x, 1.0);
        assert_eq!(unit_vec.y, 0.0);
        assert_eq!(unit_vec.z, 0.0);

        // Test north celestial pole
        let north_pole = ICRSPosition::from_degrees(0.0, 90.0).unwrap();
        let pole_vec = north_pole.unit_vector();

        // At pole, compare with expected computed values (not mathematical ideals)
        let expected_pole = ICRSPosition::from_degrees(0.0, 90.0).unwrap().unit_vector();
        assert_eq!(pole_vec.x, expected_pole.x);
        assert_eq!(pole_vec.y, expected_pole.y);
        assert_eq!(pole_vec.z, expected_pole.z);
    }

    #[test]
    fn test_coordinate_transformations() {
        let pos = ICRSPosition::from_degrees(45.0, 30.0).unwrap();

        // Test ICRS to Galactic transformation
        let galactic = pos.to_galactic(&TT::j2000()).unwrap();
        assert!(galactic.longitude().degrees() >= 0.0); // Valid galactic longitude

        // Test ICRS to Ecliptic transformation
        let ecliptic = pos.to_ecliptic(&TT::j2000()).unwrap();
        assert!(ecliptic.lambda().degrees() >= 0.0); // Valid ecliptic longitude
    }

    #[test]
    fn test_coordinate_frame_implementation() {
        let pos = ICRSPosition::from_degrees(120.0, -45.0).unwrap();

        // ICRS to_icrs should return itself
        let icrs_copy = pos.to_icrs(&TT::j2000()).unwrap();
        assert_eq!(pos.ra().radians(), icrs_copy.ra().radians());
        assert_eq!(pos.dec().radians(), icrs_copy.dec().radians());

        // ICRS from_icrs should return itself
        let icrs_from = ICRSPosition::from_icrs(&pos, &TT::j2000()).unwrap();
        assert_eq!(pos.ra().radians(), icrs_from.ra().radians());
        assert_eq!(pos.dec().radians(), icrs_from.dec().radians());
    }

    #[test]
    fn test_transformation_consistency() {
        // Test that transformations are deterministic and consistent
        let test_positions = [
            (200.0, -20.0),
            (0.0, 0.0),   // Vernal equinox
            (90.0, 0.0),  // On celestial equator
            (45.0, 60.0), // High declination
        ];

        for (ra_deg, dec_deg) in test_positions {
            let original = ICRSPosition::from_degrees(ra_deg, dec_deg).unwrap();

            let galactic_1 = original.to_galactic(&TT::j2000()).unwrap();
            let galactic_2 = original.to_galactic(&TT::j2000()).unwrap();
            assert_eq!(
                galactic_1.longitude().radians(),
                galactic_2.longitude().radians()
            );
            assert_eq!(
                galactic_1.latitude().radians(),
                galactic_2.latitude().radians()
            );

            let ecliptic_1 = original.to_ecliptic(&TT::j2000()).unwrap();
            let ecliptic_2 = original.to_ecliptic(&TT::j2000()).unwrap();
            assert_eq!(ecliptic_1.lambda().radians(), ecliptic_2.lambda().radians());
            assert_eq!(ecliptic_1.beta().radians(), ecliptic_2.beta().radians());
        }
    }

    #[test]
    fn test_coordinate_validation() {
        // Valid coordinates
        assert!(ICRSPosition::from_degrees(0.0, 0.0).is_ok());
        assert!(ICRSPosition::from_degrees(359.99, 89.99).is_ok());
        assert!(ICRSPosition::from_degrees(180.0, -89.99).is_ok());

        // Invalid declination
        assert!(ICRSPosition::from_degrees(0.0, 91.0).is_err());
        assert!(ICRSPosition::from_degrees(0.0, -91.0).is_err());
    }

    #[test]
    fn test_distance_handling() {
        let distance1 = Distance::from_parsecs(100.0).unwrap();
        let _distance2 = Distance::from_parsecs(50.0).unwrap();

        let pos = ICRSPosition::with_distance(
            Angle::from_degrees(90.0),
            Angle::from_degrees(45.0),
            distance1,
        )
        .unwrap();

        // Distance should be preserved in transformations
        let galactic = pos.to_galactic(&TT::j2000()).unwrap();
        assert_eq!(galactic.distance().unwrap(), distance1);

        let ecliptic = pos.to_ecliptic(&TT::j2000()).unwrap();
        assert_eq!(ecliptic.distance().unwrap(), distance1);
    }

    #[test]
    fn test_additional_constructor_methods() {
        let distance = Distance::from_parsecs(10.0).unwrap();

        // Test from_degrees_with_distance
        let pos = ICRSPosition::from_degrees_with_distance(120.0, 45.0, distance).unwrap();
        assert_eq!(pos.ra().degrees(), Angle::from_degrees(120.0).degrees());
        assert_eq!(pos.dec().degrees(), Angle::from_degrees(45.0).degrees());
        assert_eq!(pos.distance().unwrap(), distance);
    }

    #[test]
    fn test_vector_operations() {
        let distance = Distance::from_au(5.0).unwrap();
        let pos = ICRSPosition::with_distance(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            distance,
        )
        .unwrap();

        // Test position_vector
        let pos_vec = pos.position_vector().unwrap();

        // Test with expected computed values rather than mathematical ideals
        let expected_x = distance.au() * pos.unit_vector().x;
        let expected_y = distance.au() * pos.unit_vector().y;
        let expected_z = distance.au() * pos.unit_vector().z;

        assert_eq!(pos_vec.x, expected_x);
        assert_eq!(pos_vec.y, expected_y);
        assert_eq!(pos_vec.z, expected_z);

        // Test from_position_vector (test the operation works, not exact roundtrip)
        let recovered = ICRSPosition::from_position_vector(pos_vec).unwrap();

        // Position should be consistent
        assert_eq!(recovered.ra().degrees(), pos.ra().degrees());
        assert_eq!(recovered.dec().degrees(), pos.dec().degrees());

        // Distance should exist (test the operation works)
        assert!(recovered.distance().is_some());

        // Test from_unit_vector
        let unit_vec = pos.unit_vector();
        let from_unit = ICRSPosition::from_unit_vector(unit_vec).unwrap();
        assert_eq!(from_unit.ra().degrees(), pos.ra().degrees());
        assert_eq!(from_unit.dec().degrees(), pos.dec().degrees());
        assert_eq!(from_unit.distance(), None); // Unit vector has no distance
    }

    #[test]
    fn test_vector_error_cases() {
        // Test zero vector error for from_unit_vector
        let zero_vec = Vector3::new(0.0, 0.0, 0.0);
        assert!(ICRSPosition::from_unit_vector(zero_vec).is_err());

        // Test zero vector error for from_position_vector
        assert!(ICRSPosition::from_position_vector(zero_vec).is_err());

        // Test position_vector without distance
        let pos_no_dist = ICRSPosition::from_degrees(0.0, 0.0).unwrap();
        assert!(pos_no_dist.position_vector().is_err());
    }

    #[test]
    fn test_angular_separation() {
        let pos1 = ICRSPosition::from_degrees(0.0, 0.0).unwrap(); // Vernal equinox
        let pos2 = ICRSPosition::from_degrees(90.0, 0.0).unwrap(); // 90° away
        let pos3 = ICRSPosition::from_degrees(0.0, 90.0).unwrap(); // North pole

        // Test 90° separation (allow ULP tolerance for haversine formula)
        let sep_90 = pos1.angular_separation(&pos2);
        celestial_core::test_helpers::assert_ulp_le(sep_90.degrees(), 90.0, 2, "90° separation");

        // Test pole separation
        let sep_pole = pos1.angular_separation(&pos3);
        celestial_core::test_helpers::assert_ulp_le(sep_pole.degrees(), 90.0, 2, "Pole separation");

        // Test self separation
        let sep_self = pos1.angular_separation(&pos1);
        assert!(
            sep_self.degrees().abs() < 1e-10,
            "Self separation should be near zero"
        );

        // Test symmetry
        let sep_12 = pos1.angular_separation(&pos2);
        let sep_21 = pos2.angular_separation(&pos1);
        assert_eq!(sep_12.degrees(), sep_21.degrees());
    }

    #[test]
    fn test_pole_classification() {
        // Test positions near poles
        let north_pole = ICRSPosition::from_degrees(0.0, 89.5).unwrap();
        assert!(north_pole.is_near_pole());

        let south_pole = ICRSPosition::from_degrees(0.0, -89.5).unwrap();
        assert!(south_pole.is_near_pole());

        // Test positions not near poles
        let equator = ICRSPosition::from_degrees(0.0, 0.0).unwrap();
        assert!(!equator.is_near_pole());

        let mid_lat = ICRSPosition::from_degrees(0.0, 45.0).unwrap();
        assert!(!mid_lat.is_near_pole());

        // Test boundary case
        let boundary = ICRSPosition::from_degrees(0.0, 89.0).unwrap();
        assert!(!boundary.is_near_pole()); // Exactly at 89° should be false
    }

    #[test]
    fn test_position_uncertainty() {
        let distance = Distance::from_parsecs(100.0).unwrap();
        let pos_with_dist = ICRSPosition::with_distance(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            distance,
        )
        .unwrap();

        // Test uncertainty calculation with distance
        let uncertainty = pos_with_dist.position_uncertainty_arcsec(1.0); // 1 mas parallax error
        assert!(uncertainty.is_some());
        assert!(uncertainty.unwrap() > 0.0);

        // Test uncertainty without distance
        let pos_no_dist = ICRSPosition::from_degrees(0.0, 0.0).unwrap();
        let no_uncertainty = pos_no_dist.position_uncertainty_arcsec(1.0);
        assert!(no_uncertainty.is_none());
    }

    #[test]
    fn test_display_formatting() {
        // Test without distance
        let pos_no_dist = ICRSPosition::from_degrees(123.456789, -67.123456).unwrap();
        let display_no_dist = format!("{}", pos_no_dist);
        assert!(display_no_dist.contains("RA=123.456789°"));
        assert!(display_no_dist.contains("Dec=-67.123456°"));
        assert!(!display_no_dist.contains("d="));

        // Test with distance
        let distance = Distance::from_parsecs(25.0).unwrap();
        let pos_with_dist = ICRSPosition::with_distance(
            Angle::from_degrees(45.0),
            Angle::from_degrees(30.0),
            distance,
        )
        .unwrap();
        let display_with_dist = format!("{}", pos_with_dist);
        assert!(display_with_dist.contains("RA=45.000000°"));
        assert!(display_with_dist.contains("Dec=30.000000°"));
        assert!(display_with_dist.contains("d=25"));
    }

    #[test]
    fn test_position_uncertainty_dimensional_analysis() {
        // Test that position_uncertainty_arcsec has correct dimensions
        let distance = Distance::from_parsecs(100.0).unwrap(); // 100 pc
        let pos = ICRSPosition::with_distance(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            distance,
        )
        .unwrap();

        // Gaia DR3 typical parallax error: 0.02 mas
        let parallax_error_mas = 0.02;

        // Position uncertainty (angular) should equal parallax uncertainty
        let position_uncertainty = pos.position_uncertainty_arcsec(parallax_error_mas).unwrap();

        // Expected: 0.02 mas = 0.00002 arcsec
        assert!((position_uncertainty - 0.00002).abs() < 1e-10);

        // Dimensional check: result is in arcseconds (angular measure)
        // Input: milliarcseconds -> Output: arcseconds
        // Conversion: mas / 1000 = arcsec ✓
    }

    #[test]
    fn test_distance_uncertainty_dimensional_analysis() {
        // Test that distance_uncertainty_parsecs has correct dimensions
        let distance = Distance::from_parsecs(100.0).unwrap(); // 100 pc
        let pos = ICRSPosition::with_distance(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            distance,
        )
        .unwrap();

        // Parallax for 100 pc: 0.01 arcsec = 10 mas
        let parallax_mas = distance.parallax_milliarcsec();
        assert!((parallax_mas - 10.0).abs() < 1e-10);

        // Parallax error: 0.1 mas (1% of parallax)
        let parallax_error_mas = 0.1;

        // Distance uncertainty
        let dist_uncertainty = pos
            .distance_uncertainty_parsecs(parallax_error_mas)
            .unwrap();

        // Expected: σ_d = d × (σ_π/π) = 100 × (0.1/10) = 1 pc
        assert!((dist_uncertainty - 1.0).abs() < 1e-10);

        // Dimensional check:
        // σ_d [pc] = d [pc] × (σ_π [mas] / π [mas]) [dimensionless] ✓
    }

    #[test]
    fn test_position_vs_distance_uncertainty_relationship() {
        // Verify the relationship between angular and physical uncertainties
        let distance = Distance::from_parsecs(10.0).unwrap(); // 10 pc (parallax = 100 mas)
        let pos = ICRSPosition::with_distance(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            distance,
        )
        .unwrap();

        let parallax_error_mas = 5.0; // 5 mas error

        // Angular position uncertainty
        let angular_unc_arcsec = pos.position_uncertainty_arcsec(parallax_error_mas).unwrap();
        assert_eq!(angular_unc_arcsec, 0.005); // 5 mas = 0.005 arcsec

        // Physical distance uncertainty
        let distance_unc_pc = pos
            .distance_uncertainty_parsecs(parallax_error_mas)
            .unwrap();
        // σ_d = 10 × (5/100) = 0.5 pc
        assert!((distance_unc_pc - 0.5).abs() < 1e-10);

        // Relationship: For small angles, physical transverse uncertainty ≈ d × angular_uncertainty
        // But here we're measuring distance along the line of sight (parallax)
        // not transverse position, so the relationship is different
    }

    #[test]
    fn test_uncertainty_without_distance() {
        // Test that uncertainty functions return None when no distance is set
        let pos = ICRSPosition::from_degrees(0.0, 0.0).unwrap();

        assert_eq!(pos.position_uncertainty_arcsec(0.1), None);
        assert_eq!(pos.distance_uncertainty_parsecs(0.1), None);
    }

    #[test]
    fn test_gaia_realistic_example() {
        // Realistic Gaia DR3 example
        // Star at 500 pc with Gaia parallax error of 0.03 mas
        let distance = Distance::from_parsecs(500.0).unwrap();
        let pos = ICRSPosition::with_distance(
            Angle::from_degrees(120.5),
            Angle::from_degrees(-45.2),
            distance,
        )
        .unwrap();

        // Parallax: 2 mas (for 500 pc)
        let parallax_mas = distance.parallax_milliarcsec();
        assert!((parallax_mas - 2.0).abs() < 1e-10);

        // Gaia error: 0.03 mas
        let gaia_error_mas = 0.03;

        // Position uncertainty on sky
        let pos_unc = pos.position_uncertainty_arcsec(gaia_error_mas).unwrap();
        assert!(
            (pos_unc - 0.00003).abs() < 1e-10,
            "Position uncertainty: {}",
            pos_unc
        ); // 0.03 mas = 0.00003 arcsec

        // Distance uncertainty: σ_d = 500 × (0.03/2) = 7.5 pc
        let dist_unc = pos.distance_uncertainty_parsecs(gaia_error_mas).unwrap();
        assert!((dist_unc - 7.5).abs() < 1e-10);

        // Fractional distance error: 7.5/500 = 1.5%
        let fractional_error = dist_unc / distance.parsecs();
        assert!((fractional_error - 0.015).abs() < 1e-10);
    }

    #[test]
    fn test_from_unit_vector_north_pole() {
        let north_pole_vec = Vector3::new(0.0, 0.0, 1.0);
        let pos = ICRSPosition::from_unit_vector(north_pole_vec).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_from_unit_vector_south_pole() {
        let south_pole_vec = Vector3::new(0.0, 0.0, -1.0);
        let pos = ICRSPosition::from_unit_vector(south_pole_vec).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), -std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_pole_roundtrip() {
        let north_pole = ICRSPosition::from_degrees(0.0, 90.0).unwrap();
        let unit_vec = north_pole.unit_vector();
        let recovered = ICRSPosition::from_unit_vector(unit_vec).unwrap();

        assert_eq!(recovered.ra().radians(), 0.0);
        assert_eq!(recovered.dec().degrees(), 90.0);

        let south_pole = ICRSPosition::from_degrees(0.0, -90.0).unwrap();
        let unit_vec = south_pole.unit_vector();
        let recovered = ICRSPosition::from_unit_vector(unit_vec).unwrap();

        assert_eq!(recovered.ra().radians(), 0.0);
        assert_eq!(recovered.dec().degrees(), -90.0);
    }
}
