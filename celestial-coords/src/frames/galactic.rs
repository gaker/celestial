use super::direction::spherical_angles;
use crate::constants::ICRS_TO_GALACTIC;
use crate::distance::Distance;
use crate::errors::CoordResult;
use crate::frames::icrs::ICRSPosition;
use crate::transforms::CoordinateFrame;
use celestial_core::angle::Angle;
use celestial_core::matrix::{RotationMatrix3, Vector3};
use celestial_time::scales::tt::TT;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GalacticPosition {
    l: Angle,
    b: Angle,
    distance: Option<Distance>,
}

impl GalacticPosition {
    pub fn new(l: Angle, b: Angle) -> CoordResult<Self> {
        let l = l.normalized()?;
        let b = b.validate_latitude()?;

        Ok(Self {
            l,
            b,
            distance: None,
        })
    }

    pub fn with_distance(l: Angle, b: Angle, distance: Distance) -> CoordResult<Self> {
        let mut pos = Self::new(l, b)?;
        pos.distance = Some(distance);
        Ok(pos)
    }

    pub fn from_degrees(l_deg: f64, b_deg: f64) -> CoordResult<Self> {
        Self::new(Angle::from_degrees(l_deg), Angle::from_degrees(b_deg))
    }

    pub fn longitude(&self) -> Angle {
        self.l
    }

    pub fn latitude(&self) -> Angle {
        self.b
    }

    pub fn distance(&self) -> Option<Distance> {
        self.distance
    }

    pub fn set_distance(&mut self, distance: Distance) {
        self.distance = Some(distance);
    }

    pub fn galactic_center() -> Self {
        Self {
            l: Angle::ZERO,
            b: Angle::ZERO,
            distance: None,
        }
    }

    pub fn galactic_anticenter() -> Self {
        Self {
            l: Angle::PI,
            b: Angle::ZERO,
            distance: None,
        }
    }

    pub fn north_galactic_pole() -> Self {
        Self {
            l: Angle::ZERO,
            b: Angle::HALF_PI,
            distance: None,
        }
    }

    pub fn south_galactic_pole() -> Self {
        Self {
            l: Angle::ZERO,
            b: -Angle::HALF_PI,
            distance: None,
        }
    }

    pub fn is_near_galactic_plane(&self) -> bool {
        self.b.abs().degrees() < 10.0
    }

    pub fn is_in_galactic_bulge(&self) -> bool {
        self.b.abs().degrees() < 10.0 && (self.l.degrees() < 30.0 || self.l.degrees() > 330.0)
    }

    pub fn is_near_galactic_pole(&self) -> bool {
        self.b.abs().degrees() > 80.0
    }

    pub fn angular_distance_from_gc(&self) -> Angle {
        let gc = Self::galactic_center();
        self.angular_separation(&gc)
    }

    pub fn angular_separation(&self, other: &Self) -> Angle {
        Angle::from_radians(celestial_core::math::angular_separation(
            self.l.radians(),
            self.b.radians(),
            other.l.radians(),
            other.b.radians(),
        ))
    }
}

impl CoordinateFrame for GalacticPosition {
    fn to_icrs(&self, _epoch: &TT) -> CoordResult<ICRSPosition> {
        let galactic = Vector3::from_spherical(self.l.radians(), self.b.radians());
        let mut icrs = ICRSPosition::from_unit_vector(icrs_to_galactic()?.transpose() * galactic)?;
        if let Some(distance) = self.distance {
            icrs.set_distance(distance);
        }
        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, _epoch: &TT) -> CoordResult<Self> {
        let direction = Vector3::from_spherical(icrs.ra().radians(), icrs.dec().radians());
        let (l, b) = spherical_angles(icrs_to_galactic()? * direction)?;
        let mut galactic = Self::new(l, b)?;
        galactic.distance = icrs.distance();
        Ok(galactic)
    }
}

fn icrs_to_galactic() -> CoordResult<RotationMatrix3> {
    Ok(RotationMatrix3::from_array(ICRS_TO_GALACTIC)?)
}

impl std::fmt::Display for GalacticPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Galactic(l={:.6}°, b={:.6}°",
            self.l.degrees(),
            self.b.degrees()
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

    #[test]
    fn test_galactic_creation() {
        let pos = GalacticPosition::from_degrees(45.0, 30.0).unwrap();
        assert_eq!(pos.longitude(), Angle::from_degrees(45.0));
        assert_eq!(pos.latitude(), Angle::from_degrees(30.0));
        assert!(pos.distance().is_none());
    }

    #[test]
    fn test_galactic_validation() {
        assert!(GalacticPosition::from_degrees(0.0, 0.0).is_ok());
        assert!(GalacticPosition::from_degrees(359.999, 89.999).is_ok());

        // eraAnp of 380 degrees: 20 degrees, less the rounding of 380 degrees in radians.
        let pos = GalacticPosition::from_degrees(380.0, 45.0).unwrap();
        assert_eq!(pos.longitude().radians(), 0.3490658503988664);

        assert!(GalacticPosition::from_degrees(0.0, 95.0).is_err());
        assert!(GalacticPosition::from_degrees(0.0, -95.0).is_err());
    }

    #[test]
    fn test_special_positions() {
        let gc = GalacticPosition::galactic_center();
        assert_eq!(gc.longitude().degrees(), 0.0);
        assert_eq!(gc.latitude().degrees(), 0.0);

        let gac = GalacticPosition::galactic_anticenter();
        assert_eq!(gac.longitude().degrees(), 180.0);
        assert_eq!(gac.latitude().degrees(), 0.0);

        let ngp = GalacticPosition::north_galactic_pole();
        assert_eq!(ngp.latitude().degrees(), 90.0);

        let sgp = GalacticPosition::south_galactic_pole();
        assert_eq!(sgp.latitude().degrees(), -90.0);
    }

    #[test]
    fn test_galactic_regions() {
        let plane_pos = GalacticPosition::from_degrees(45.0, 5.0).unwrap();
        assert!(plane_pos.is_near_galactic_plane());
        assert!(!plane_pos.is_near_galactic_pole());

        let bulge_pos = GalacticPosition::from_degrees(5.0, 5.0).unwrap();
        assert!(bulge_pos.is_in_galactic_bulge());

        let pole_pos = GalacticPosition::from_degrees(0.0, 85.0).unwrap();
        assert!(pole_pos.is_near_galactic_pole());
        assert!(!pole_pos.is_near_galactic_plane());
    }

    #[test]
    fn test_angular_separation() {
        let pos1 = GalacticPosition::from_degrees(0.0, 0.0).unwrap();
        let pos2 = GalacticPosition::from_degrees(90.0, 0.0).unwrap();

        assert_eq!(pos1.angular_separation(&pos2).degrees(), 90.0);
        assert_eq!(pos2.angular_distance_from_gc().degrees(), 90.0);
    }

    #[test]
    fn test_coordinate_transformations() {
        let epoch = TT::j2000();
        let gal_pos = GalacticPosition::from_degrees(45.0, 30.0).unwrap();

        // eraG2icrs, then eraIcrs2g on its output.
        let icrs = gal_pos.to_icrs(&epoch).unwrap();
        assert_eq!(
            [icrs.ra().radians(), icrs.dec().radians()],
            [4.5324553575768585, 0.3996934769842426]
        );
        let back = GalacticPosition::from_icrs(&icrs, &epoch).unwrap();
        assert_eq!(
            [back.longitude().radians(), back.latitude().radians()],
            [0.7853981633974482, 0.5235987755982991]
        );
    }

    #[test]
    fn test_with_distance() {
        let distance = Distance::from_parsecs(100.0).unwrap();
        let pos = GalacticPosition::with_distance(
            Angle::from_degrees(45.0),
            Angle::from_degrees(30.0),
            distance,
        )
        .unwrap();

        assert_eq!(pos.distance().unwrap().parsecs(), 100.0);

        let epoch = TT::j2000();
        let icrs = pos.to_icrs(&epoch).unwrap();
        assert_eq!(icrs.distance().unwrap().parsecs(), 100.0);
    }
}
