use crate::distance::Distance;
use crate::errors::CoordResult;
use crate::frames::icrs::ICRSPosition;
use crate::lunar;
use crate::transforms::CoordinateFrame;
use celestial_core::angle::wrap_0_2pi;
use celestial_core::angle::Angle;
use celestial_time::scales::tt::TT;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SelenographicPosition {
    latitude: Angle,
    longitude: Angle,
    radius: Option<Distance>,
}

impl SelenographicPosition {
    pub fn new(latitude: Angle, longitude: Angle) -> CoordResult<Self> {
        let latitude = latitude.validate_latitude()?;
        let longitude = longitude.normalized()?;

        Ok(Self {
            latitude,
            longitude,
            radius: None,
        })
    }

    pub fn with_radius(latitude: Angle, longitude: Angle, radius: Distance) -> CoordResult<Self> {
        let mut pos = Self::new(latitude, longitude)?;
        pos.radius = Some(radius);
        Ok(pos)
    }

    pub fn from_degrees(lat_deg: f64, lon_deg: f64) -> CoordResult<Self> {
        Self::new(Angle::from_degrees(lat_deg), Angle::from_degrees(lon_deg))
    }

    pub fn latitude(&self) -> Angle {
        self.latitude
    }

    pub fn longitude(&self) -> Angle {
        self.longitude
    }

    pub fn radius(&self) -> Option<Distance> {
        self.radius
    }

    pub(crate) fn set_radius(&mut self, radius: Distance) {
        self.radius = Some(radius);
    }

    pub fn sub_earth_point(epoch: &TT) -> CoordResult<Self> {
        let (lon, lat) = lunar::compute_sub_earth_point(epoch)?;
        Self::new(lat, lon)
    }

    pub fn nearside_center() -> Self {
        Self {
            latitude: Angle::ZERO,
            longitude: Angle::ZERO,
            radius: None,
        }
    }

    pub fn farside_center() -> Self {
        Self {
            latitude: Angle::ZERO,
            longitude: Angle::PI,
            radius: None,
        }
    }

    pub fn north_pole() -> Self {
        Self {
            latitude: Angle::HALF_PI,
            longitude: Angle::ZERO,
            radius: None,
        }
    }

    pub fn south_pole() -> Self {
        Self {
            latitude: -Angle::HALF_PI,
            longitude: Angle::ZERO,
            radius: None,
        }
    }

    pub fn angular_separation(&self, other: &Self) -> Angle {
        Angle::from_radians(celestial_core::math::angular_separation(
            self.longitude.radians(),
            self.latitude.radians(),
            other.longitude.radians(),
            other.latitude.radians(),
        ))
    }

    pub fn is_visible_from_earth(&self, epoch: &TT) -> CoordResult<bool> {
        let sub_earth = Self::sub_earth_point(epoch)?;
        Ok(self.angular_separation(&sub_earth).degrees() < 90.0)
    }
}

impl CoordinateFrame for SelenographicPosition {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition> {
        let m = lunar::icrs_to_selenographic(epoch)?.transpose();
        let (ra, dec) = m.transform_spherical(self.longitude.radians(), self.latitude.radians());

        let mut icrs = ICRSPosition::new(
            Angle::from_radians(wrap_0_2pi(ra)?),
            Angle::from_radians(dec),
        )?;

        if let Some(radius) = self.radius {
            icrs.set_distance(radius);
        }
        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let m = lunar::icrs_to_selenographic(epoch)?;
        let (lon, lat) = m.transform_spherical(icrs.ra().radians(), icrs.dec().radians());

        let mut pos = Self::new(
            Angle::from_radians(lat),
            Angle::from_radians(wrap_0_2pi(lon)?),
        )?;

        if let Some(dist) = icrs.distance() {
            pos.set_radius(dist);
        }
        Ok(pos)
    }
}

impl std::fmt::Display for SelenographicPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Selenographic(lat={:.6}°, lon={:.6}°",
            self.latitude.degrees(),
            self.longitude.degrees()
        )?;

        if let Some(radius) = self.radius {
            write!(f, ", r={}", radius)?;
        }

        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::rounded;
    use celestial_core::matrix::Vector3;
    use celestial_time::julian::JulianDate;

    #[test]
    fn test_selenographic_creation() {
        let pos = SelenographicPosition::from_degrees(45.0, 30.0).unwrap();
        assert_eq!(pos.latitude(), Angle::from_degrees(45.0));
        assert_eq!(pos.longitude(), Angle::from_degrees(30.0));
        assert!(pos.radius().is_none());
    }

    #[test]
    fn test_selenographic_validation() {
        assert!(SelenographicPosition::from_degrees(0.0, 0.0).is_ok());
        assert!(SelenographicPosition::from_degrees(90.0, 180.0).is_ok());
        assert!(SelenographicPosition::from_degrees(-90.0, 359.0).is_ok());

        assert!(SelenographicPosition::from_degrees(95.0, 0.0).is_err());
        assert!(SelenographicPosition::from_degrees(-95.0, 0.0).is_err());
    }

    #[test]
    fn test_special_positions() {
        let nearside = SelenographicPosition::nearside_center();
        assert_eq!(nearside.latitude().degrees(), 0.0);
        assert_eq!(nearside.longitude().degrees(), 0.0);

        let farside = SelenographicPosition::farside_center();
        assert_eq!(farside.latitude().degrees(), 0.0);
        assert_eq!(farside.longitude().degrees(), 180.0);

        let north_pole = SelenographicPosition::north_pole();
        assert_eq!(north_pole.latitude().degrees(), 90.0);

        let south_pole = SelenographicPosition::south_pole();
        assert_eq!(south_pole.latitude().degrees(), -90.0);
    }

    #[test]
    fn test_angular_separation() {
        let nearside = SelenographicPosition::nearside_center();
        let farside = SelenographicPosition::farside_center();

        assert_eq!(nearside.angular_separation(&farside).degrees(), 180.0);

        let north = SelenographicPosition::north_pole();
        assert_eq!(nearside.angular_separation(&north).degrees(), 90.0);
    }

    #[test]
    fn test_visibility_from_earth_farside() {
        let epoch = TT::j2000();

        let farside = SelenographicPosition::farside_center();
        assert!(!farside.is_visible_from_earth(&epoch).unwrap());
    }

    #[test]
    fn test_sub_earth_point() {
        let epoch = TT::j2000();
        let sub_earth = SelenographicPosition::sub_earth_point(&epoch).unwrap();

        assert!(
            libm::fabs(sub_earth.latitude().degrees()) <= 7.5,
            "Sub-earth latitude = {}",
            sub_earth.latitude().degrees()
        );
        assert!(
            sub_earth.longitude().degrees() >= 0.0 && sub_earth.longitude().degrees() < 360.0,
            "Sub-earth longitude = {}",
            sub_earth.longitude().degrees()
        );
    }

    // eraMoon98 at Meeus' example 53.a; the Earth lies the opposite way.
    fn earth_from_moon() -> ICRSPosition {
        let moon = Vector3::new(
            -0.0016853456834658809,
            0.0016977178353056598,
            0.0005848854873104547,
        );
        let (ra, dec) = (-moon).to_spherical();
        ICRSPosition::new(
            Angle::from_radians(wrap_0_2pi(ra).unwrap()),
            Angle::from_radians(dec),
        )
        .unwrap()
    }

    fn meeus_53a() -> TT {
        TT::from_julian_date(JulianDate::new(2448724.5, 0.0))
    }

    #[test]
    fn test_earth_direction_maps_to_meeus_total_libration() {
        let p = SelenographicPosition::from_icrs(&earth_from_moon(), &meeus_53a()).unwrap();
        let printed = [p.latitude(), p.longitude()].map(|a| rounded(a.degrees(), 2));
        assert_eq!(printed, [4.20, 358.77]);
    }

    #[test]
    fn test_sub_earth_point_is_meeus_total_libration() {
        let p = SelenographicPosition::sub_earth_point(&meeus_53a()).unwrap();
        let printed = [p.latitude(), p.longitude()].map(|a| rounded(a.degrees(), 2));
        assert_eq!(printed, [4.20, 358.77]);
    }

    #[test]
    fn test_non_finite_epoch_is_err() {
        let epoch = TT::from_julian_date(JulianDate::new(f64::NAN, 0.0));
        let position = SelenographicPosition::from_degrees(10.0, 20.0).unwrap();
        let icrs = ICRSPosition::from_degrees(10.0, 20.0).unwrap();
        assert!(position.to_icrs(&epoch).is_err());
        assert!(SelenographicPosition::from_icrs(&icrs, &epoch).is_err());
        assert!(SelenographicPosition::sub_earth_point(&epoch).is_err());
        assert!(position.is_visible_from_earth(&epoch).is_err());
    }

    #[test]
    fn test_coordinate_frame_to_icrs() {
        let epoch = TT::j2000();
        let original = SelenographicPosition::from_degrees(0.0, 0.0).unwrap();

        let icrs = original.to_icrs(&epoch).unwrap();

        assert!(icrs.ra().degrees() >= 0.0 && icrs.ra().degrees() < 360.0);
        assert!(icrs.dec().degrees() >= -90.0 && icrs.dec().degrees() <= 90.0);
    }

    #[test]
    fn test_coordinate_frame_roundtrip() {
        let epoch = TT::j2000();
        let test_cases = [
            (0.0, 0.0),
            (5.0, 30.0),
            (-5.0, 90.0),
            (3.0, 180.0),
            (-3.0, 270.0),
        ];

        for (lat, lon) in test_cases {
            let original = SelenographicPosition::from_degrees(lat, lon).unwrap();
            let icrs = original.to_icrs(&epoch).unwrap();
            let recovered = SelenographicPosition::from_icrs(&icrs, &epoch).unwrap();

            // Two rotations and their inverses, so a few units in the last place of a radian.
            let error = original.angular_separation(&recovered).radians();
            assert!(error < 1e-15, "({lat}, {lon}): {error:e} rad");
        }
    }

    #[test]
    fn test_with_radius() {
        let radius = Distance::from_au(0.00257).unwrap();
        let pos = SelenographicPosition::with_radius(
            Angle::from_degrees(0.0),
            Angle::from_degrees(0.0),
            radius,
        )
        .unwrap();

        assert!(pos.radius().is_some());
        assert_eq!(pos.radius().unwrap(), radius);
    }

    #[test]
    fn test_display_formatting() {
        let pos = SelenographicPosition::from_degrees(45.123456, 30.654321).unwrap();
        let display = format!("{}", pos);
        assert!(display.contains("45.123456"));
        assert!(display.contains("30.654321"));
        assert!(display.contains("Selenographic"));
    }
}
