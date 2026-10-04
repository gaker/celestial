use crate::distance::Distance;
use crate::errors::{CoordError, CoordResult};
use crate::frames::direction::spherical_angles;
use crate::frames::ecliptic::EclipticPosition;
use crate::frames::galactic::GalacticPosition;
use crate::transforms::CoordinateFrame;
use celestial_core::{angle::Angle, matrix::Vector3};
use celestial_time::scales::tt::TT;

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
        let (ra, dec) = spherical_angles(unit)?;
        Self::new(ra, dec)
    }

    pub fn from_position_vector(pos: Vector3) -> CoordResult<Self> {
        let mut icrs = Self::from_unit_vector(pos)?;
        icrs.distance = Some(Distance::from_au(pos.magnitude())?);

        Ok(icrs)
    }

    pub fn angular_separation(&self, other: &Self) -> Angle {
        Angle::from_radians(celestial_core::math::angular_separation(
            self.ra.radians(),
            self.dec.radians(),
            other.ra.radians(),
            other.dec.radians(),
        ))
    }

    pub fn is_near_pole(&self) -> bool {
        self.dec.abs().degrees() > 89.0
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
    pub fn to_galactic(&self, epoch: &TT) -> CoordResult<GalacticPosition> {
        GalacticPosition::from_icrs(self, epoch)
    }

    pub fn to_ecliptic(&self, epoch: &TT) -> CoordResult<EclipticPosition> {
        EclipticPosition::from_icrs(self, epoch)
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
    use crate::distance::Distance;

    #[test]
    fn test_constructor_methods() {
        let pos1 =
            ICRSPosition::new(Angle::from_degrees(180.0), Angle::from_degrees(45.0)).unwrap();

        assert_eq!(pos1.ra().degrees(), Angle::from_degrees(180.0).degrees());
        assert_eq!(pos1.dec().degrees(), Angle::from_degrees(45.0).degrees());
        assert_eq!(pos1.distance(), None);

        let pos2 = ICRSPosition::from_degrees(90.0, -30.0).unwrap();
        assert_eq!(pos2.ra().degrees(), Angle::from_degrees(90.0).degrees());
        assert_eq!(pos2.dec().degrees(), Angle::from_degrees(-30.0).degrees());

        let pos3 = ICRSPosition::from_hours_degrees(12.0, 60.0).unwrap();
        assert_eq!(pos3.ra().hours(), 12.0);
        assert_eq!(pos3.dec().degrees(), Angle::from_degrees(60.0).degrees());

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

        assert_eq!(pos.ra().degrees(), Angle::from_degrees(270.0).degrees());
        assert_eq!(pos.dec().degrees(), Angle::from_degrees(15.0).degrees());
        assert_eq!(pos.distance(), None);

        pos.set_distance(distance);
        assert_eq!(pos.distance().unwrap(), distance);

        pos.remove_distance();
        assert_eq!(pos.distance(), None);
    }

    #[test]
    fn test_unit_vector_conversion() {
        let vernal_equinox = ICRSPosition::from_degrees(0.0, 0.0).unwrap();
        let unit_vec = vernal_equinox.unit_vector();

        assert_eq!(unit_vec.x, 1.0);
        assert_eq!(unit_vec.y, 0.0);
        assert_eq!(unit_vec.z, 0.0);

        // x is the cosine of the double nearest pi/2, as eraS2c gives.
        let pole_vec = ICRSPosition::from_degrees(0.0, 90.0).unwrap().unit_vector();
        assert_eq!(
            [pole_vec.x, pole_vec.y, pole_vec.z],
            [6.123233995736766e-17, 0.0, 1.0]
        );
    }

    #[test]
    fn test_coordinate_transformations() {
        let pos = ICRSPosition::from_degrees(45.0, 30.0).unwrap();

        let galactic = pos.to_galactic(&TT::j2000()).unwrap();
        assert!(galactic.longitude().degrees() >= 0.0);

        let ecliptic = pos.to_ecliptic(&TT::j2000()).unwrap();
        assert!(ecliptic.lambda().degrees() >= 0.0);
    }

    #[test]
    fn test_coordinate_frame_implementation() {
        let pos = ICRSPosition::from_degrees(120.0, -45.0).unwrap();

        let icrs_copy = pos.to_icrs(&TT::j2000()).unwrap();
        assert_eq!(pos.ra().radians(), icrs_copy.ra().radians());
        assert_eq!(pos.dec().radians(), icrs_copy.dec().radians());

        let icrs_from = ICRSPosition::from_icrs(&pos, &TT::j2000()).unwrap();
        assert_eq!(pos.ra().radians(), icrs_from.ra().radians());
        assert_eq!(pos.dec().radians(), icrs_from.dec().radians());
    }

    #[test]
    fn test_transformation_consistency() {
        let test_positions = [(200.0, -20.0), (0.0, 0.0), (90.0, 0.0), (45.0, 60.0)];

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
        assert!(ICRSPosition::from_degrees(0.0, 0.0).is_ok());
        assert!(ICRSPosition::from_degrees(359.99, 89.99).is_ok());
        assert!(ICRSPosition::from_degrees(180.0, -89.99).is_ok());

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

        let galactic = pos.to_galactic(&TT::j2000()).unwrap();
        assert_eq!(galactic.distance().unwrap(), distance1);

        let ecliptic = pos.to_ecliptic(&TT::j2000()).unwrap();
        assert_eq!(ecliptic.distance().unwrap(), distance1);
    }

    #[test]
    fn test_additional_constructor_methods() {
        let distance = Distance::from_parsecs(10.0).unwrap();

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

        let pos_vec = pos.position_vector().unwrap();
        assert_eq!([pos_vec.x, pos_vec.y, pos_vec.z], [5.0, 0.0, 0.0]);

        let recovered = ICRSPosition::from_position_vector(pos_vec).unwrap();
        assert_eq!(recovered.ra().degrees(), pos.ra().degrees());
        assert_eq!(recovered.dec().degrees(), pos.dec().degrees());
        assert!(recovered.distance().is_some());

        let unit_vec = pos.unit_vector();
        let from_unit = ICRSPosition::from_unit_vector(unit_vec).unwrap();
        assert_eq!(from_unit.ra().degrees(), pos.ra().degrees());
        assert_eq!(from_unit.dec().degrees(), pos.dec().degrees());
        assert_eq!(from_unit.distance(), None);
    }

    #[test]
    fn test_vector_error_cases() {
        let zero_vec = Vector3::new(0.0, 0.0, 0.0);
        assert!(ICRSPosition::from_unit_vector(zero_vec).is_err());

        assert!(ICRSPosition::from_position_vector(zero_vec).is_err());

        let pos_no_dist = ICRSPosition::from_degrees(0.0, 0.0).unwrap();
        assert!(pos_no_dist.position_vector().is_err());
    }

    #[test]
    fn test_angular_separation() {
        let pos1 = ICRSPosition::from_degrees(0.0, 0.0).unwrap();
        let pos2 = ICRSPosition::from_degrees(90.0, 0.0).unwrap();
        let pos3 = ICRSPosition::from_degrees(0.0, 90.0).unwrap();

        assert_eq!(pos1.angular_separation(&pos2).degrees(), 90.0);
        assert_eq!(pos1.angular_separation(&pos3).degrees(), 90.0);
        assert_eq!(pos1.angular_separation(&pos1).degrees(), 0.0);

        let sep_12 = pos1.angular_separation(&pos2);
        let sep_21 = pos2.angular_separation(&pos1);
        assert_eq!(sep_12.degrees(), sep_21.degrees());
    }

    #[test]
    fn test_pole_classification() {
        let north_pole = ICRSPosition::from_degrees(0.0, 89.5).unwrap();
        assert!(north_pole.is_near_pole());

        let south_pole = ICRSPosition::from_degrees(0.0, -89.5).unwrap();
        assert!(south_pole.is_near_pole());

        let equator = ICRSPosition::from_degrees(0.0, 0.0).unwrap();
        assert!(!equator.is_near_pole());

        let mid_lat = ICRSPosition::from_degrees(0.0, 45.0).unwrap();
        assert!(!mid_lat.is_near_pole());

        let boundary = ICRSPosition::from_degrees(0.0, 89.0).unwrap();
        assert!(!boundary.is_near_pole());
    }

    #[test]
    fn test_display_formatting() {
        let pos_no_dist = ICRSPosition::from_degrees(123.456789, -67.123456).unwrap();
        let display_no_dist = format!("{}", pos_no_dist);
        assert!(display_no_dist.contains("RA=123.456789°"));
        assert!(display_no_dist.contains("Dec=-67.123456°"));
        assert!(!display_no_dist.contains("d="));

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
    fn test_distance_uncertainty_dimensional_analysis() {
        let distance = Distance::from_parsecs(100.0).unwrap();
        let pos = ICRSPosition::with_distance(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            distance,
        )
        .unwrap();

        assert_eq!(distance.parallax_milliarcsec(), 10.0);

        // σ_d = d × (σ_π/π) = 100 pc × (0.1 mas / 10 mas) = 1 pc
        let dist_uncertainty = pos.distance_uncertainty_parsecs(0.1).unwrap();
        assert_eq!(dist_uncertainty, 1.0);
    }

    #[test]
    fn test_uncertainty_without_distance() {
        let pos = ICRSPosition::from_degrees(0.0, 0.0).unwrap();

        assert_eq!(pos.distance_uncertainty_parsecs(0.1), None);
    }

    #[test]
    fn test_gaia_realistic_example() {
        // A star at 500 pc with a Gaia DR3-like parallax error of 0.03 mas.
        let distance = Distance::from_parsecs(500.0).unwrap();
        let pos = ICRSPosition::with_distance(
            Angle::from_degrees(120.5),
            Angle::from_degrees(-45.2),
            distance,
        )
        .unwrap();

        assert_eq!(distance.parallax_milliarcsec(), 2.0);

        // σ_d = 500 pc × (0.03 mas / 2 mas) = 7.5 pc, 1.5% of the distance
        let dist_unc = pos.distance_uncertainty_parsecs(0.03).unwrap();
        assert_eq!(dist_unc, 7.5);
        assert_eq!(dist_unc / distance.parsecs(), 0.015);
    }

    #[test]
    fn test_from_unit_vector_north_pole() {
        let north_pole_vec = Vector3::new(0.0, 0.0, 1.0);
        let pos = ICRSPosition::from_unit_vector(north_pole_vec).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), celestial_core::constants::HALF_PI);
    }

    #[test]
    fn test_from_unit_vector_south_pole() {
        let south_pole_vec = Vector3::new(0.0, 0.0, -1.0);
        let pos = ICRSPosition::from_unit_vector(south_pole_vec).unwrap();

        assert_eq!(pos.ra().radians(), 0.0);
        assert_eq!(pos.dec().radians(), -celestial_core::constants::HALF_PI);
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
