use crate::{
    constants::GALACTIC_TO_ICRS, transforms::CoordinateFrame, CoordResult, Distance, ICRSPosition,
};
use astro_core::Angle;
use astro_time::TT;

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
        let l = l.validate_longitude(true)?;
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
        let (sin_b1, cos_b1) = self.b.sin_cos();
        let (sin_b2, cos_b2) = other.b.sin_cos();
        let delta_l = (self.l - other.l).radians();

        let angle_rad =
            astro_core::math::vincenty_angular_separation(sin_b1, cos_b1, sin_b2, cos_b2, delta_l);

        Angle::from_radians(angle_rad)
    }
}

impl CoordinateFrame for GalacticPosition {
    fn to_icrs(&self, _epoch: &TT) -> CoordResult<ICRSPosition> {
        let (sin_b, cos_b) = self.b.sin_cos();
        let (sin_l, cos_l) = self.l.sin_cos();
        let gal_cartesian = [cos_l * cos_b, sin_l * cos_b, sin_b];

        // Matrix multiplication: icrs = M^T * gal (transpose because matrix is stored as columns)
        // This works correctly because GALACTIC_TO_ICRS is orthonormal (rotation matrix).
        let icrs_cartesian = [
            GALACTIC_TO_ICRS[0][0] * gal_cartesian[0]
                + GALACTIC_TO_ICRS[1][0] * gal_cartesian[1]
                + GALACTIC_TO_ICRS[2][0] * gal_cartesian[2],
            GALACTIC_TO_ICRS[0][1] * gal_cartesian[0]
                + GALACTIC_TO_ICRS[1][1] * gal_cartesian[1]
                + GALACTIC_TO_ICRS[2][1] * gal_cartesian[2],
            GALACTIC_TO_ICRS[0][2] * gal_cartesian[0]
                + GALACTIC_TO_ICRS[1][2] * gal_cartesian[1]
                + GALACTIC_TO_ICRS[2][2] * gal_cartesian[2],
        ];

        let d2 = icrs_cartesian[0] * icrs_cartesian[0] + icrs_cartesian[1] * icrs_cartesian[1];
        let ra = if d2 != 0.0 {
            icrs_cartesian[1].atan2(icrs_cartesian[0])
        } else {
            0.0
        };
        let dec = if d2 != 0.0 || icrs_cartesian[2] != 0.0 {
            icrs_cartesian[2].atan2(d2.sqrt())
        } else {
            0.0
        };

        let mut icrs = ICRSPosition::new(Angle::from_radians(ra), Angle::from_radians(dec))?;

        if let Some(distance) = self.distance {
            icrs.set_distance(distance);
        }

        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, _epoch: &TT) -> CoordResult<Self> {
        let (sin_dec, cos_dec) = icrs.dec().sin_cos();
        let (sin_ra, cos_ra) = icrs.ra().sin_cos();
        let icrs_cartesian = [cos_ra * cos_dec, sin_ra * cos_dec, sin_dec];

        // Matrix multiplication: gal = M * icrs (standard row-major access)
        // For orthonormal matrices, M^T = M^(-1), so this is the inverse of to_icrs.
        let gal_cartesian = [
            GALACTIC_TO_ICRS[0][0] * icrs_cartesian[0]
                + GALACTIC_TO_ICRS[0][1] * icrs_cartesian[1]
                + GALACTIC_TO_ICRS[0][2] * icrs_cartesian[2],
            GALACTIC_TO_ICRS[1][0] * icrs_cartesian[0]
                + GALACTIC_TO_ICRS[1][1] * icrs_cartesian[1]
                + GALACTIC_TO_ICRS[1][2] * icrs_cartesian[2],
            GALACTIC_TO_ICRS[2][0] * icrs_cartesian[0]
                + GALACTIC_TO_ICRS[2][1] * icrs_cartesian[1]
                + GALACTIC_TO_ICRS[2][2] * icrs_cartesian[2],
        ];

        let d2 = gal_cartesian[0] * gal_cartesian[0] + gal_cartesian[1] * gal_cartesian[1];
        let l = if d2 != 0.0 {
            gal_cartesian[1].atan2(gal_cartesian[0])
        } else {
            0.0
        };
        let b = if d2 != 0.0 || gal_cartesian[2] != 0.0 {
            gal_cartesian[2].atan2(d2.sqrt())
        } else {
            0.0
        };

        let mut galactic = Self::new(Angle::from_radians(l), Angle::from_radians(b))?;

        if let Some(distance) = icrs.distance() {
            galactic.set_distance(distance);
        }

        Ok(galactic)
    }
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
        assert!((pos.longitude().degrees() - 45.0).abs() < 1e-12);
        assert!((pos.latitude().degrees() - 30.0).abs() < 1e-12);
        assert!(pos.distance().is_none());
    }

    #[test]
    fn test_galactic_validation() {
        // Valid coordinates
        assert!(GalacticPosition::from_degrees(0.0, 0.0).is_ok());
        assert!(GalacticPosition::from_degrees(359.999, 89.999).is_ok());

        // Longitude gets normalized
        let pos = GalacticPosition::from_degrees(380.0, 45.0).unwrap();
        assert!((pos.longitude().degrees() - 20.0).abs() < 1e-12);

        // Invalid latitude
        assert!(GalacticPosition::from_degrees(0.0, 95.0).is_err());
        assert!(GalacticPosition::from_degrees(0.0, -95.0).is_err());
    }

    #[test]
    fn test_special_positions() {
        let gc = GalacticPosition::galactic_center();
        assert_eq!(gc.longitude().degrees(), 0.0);
        assert_eq!(gc.latitude().degrees(), 0.0);

        let gac = GalacticPosition::galactic_anticenter();
        assert!((gac.longitude().degrees() - 180.0).abs() < 1e-12);
        assert_eq!(gac.latitude().degrees(), 0.0);

        let ngp = GalacticPosition::north_galactic_pole();
        assert!((ngp.latitude().degrees() - 90.0).abs() < 1e-12);

        let sgp = GalacticPosition::south_galactic_pole();
        assert!((sgp.latitude().degrees() - (-90.0)).abs() < 1e-12);
    }

    #[test]
    fn test_galactic_regions() {
        // Galactic plane
        let plane_pos = GalacticPosition::from_degrees(45.0, 5.0).unwrap();
        assert!(plane_pos.is_near_galactic_plane());
        assert!(!plane_pos.is_near_galactic_pole());

        // Galactic bulge
        let bulge_pos = GalacticPosition::from_degrees(5.0, 5.0).unwrap();
        assert!(bulge_pos.is_in_galactic_bulge());

        // Galactic pole
        let pole_pos = GalacticPosition::from_degrees(0.0, 85.0).unwrap();
        assert!(pole_pos.is_near_galactic_pole());
        assert!(!pole_pos.is_near_galactic_plane());
    }

    #[test]
    fn test_angular_separation() {
        let pos1 = GalacticPosition::from_degrees(0.0, 0.0).unwrap();
        let pos2 = GalacticPosition::from_degrees(90.0, 0.0).unwrap();

        let sep = pos1.angular_separation(&pos2);
        // Should be approximately 90 degrees
        assert!((sep.degrees() - 90.0).abs() < 1.0); // Allow 1° tolerance for approximation

        // Distance from galactic center
        let gc_dist = pos2.angular_distance_from_gc();
        assert!((gc_dist.degrees() - 90.0).abs() < 1.0);
    }

    #[test]
    fn test_coordinate_transformations() {
        let epoch = TT::j2000();
        let gal_pos = GalacticPosition::from_degrees(45.0, 30.0).unwrap();

        // Test Galactic -> ICRS -> Galactic round trip
        let icrs = gal_pos.to_icrs(&epoch).unwrap();
        let gal_recovered = GalacticPosition::from_icrs(&icrs, &epoch).unwrap();

        assert!(
            (gal_recovered.longitude().degrees() - gal_pos.longitude().degrees()).abs() < 1e-12
        );
        assert!((gal_recovered.latitude().degrees() - gal_pos.latitude().degrees()).abs() < 1e-12);
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

        // Test coordinate transformation preserves distance
        let epoch = TT::j2000();
        let icrs = pos.to_icrs(&epoch).unwrap();
        assert_eq!(icrs.distance().unwrap().parsecs(), 100.0);
    }
}
