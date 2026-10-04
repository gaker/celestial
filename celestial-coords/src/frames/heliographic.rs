use crate::distance::Distance;
use crate::errors::CoordResult;
use crate::frames::icrs::ICRSPosition;
use crate::solar;
use crate::transforms::CoordinateFrame;
use celestial_core::angle::wrap_0_2pi;
use celestial_core::angle::Angle;
use celestial_time::scales::tt::TT;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HeliographicStonyhurst {
    latitude: Angle,
    longitude: Angle,
    radius: Option<Distance>,
}

impl HeliographicStonyhurst {
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

    pub fn to_carrington(&self, epoch: &TT) -> CoordResult<HeliographicCarrington> {
        let l0 = solar::compute_l0(epoch)?;
        let carrington_lon = self.longitude + l0;
        let normalized_lon = Angle::from_radians(wrap_0_2pi(carrington_lon.radians())?);

        let mut carr = HeliographicCarrington::new(self.latitude, normalized_lon)?;
        if let Some(r) = self.radius {
            carr.set_radius(r);
        }
        Ok(carr)
    }

    pub fn disk_center(epoch: &TT) -> CoordResult<Self> {
        Ok(Self {
            latitude: solar::compute_b0(epoch)?,
            longitude: Angle::ZERO,
            radius: None,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HeliographicCarrington {
    latitude: Angle,
    longitude: Angle,
    radius: Option<Distance>,
}

impl HeliographicCarrington {
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

    pub fn to_stonyhurst(&self, epoch: &TT) -> CoordResult<HeliographicStonyhurst> {
        let l0 = solar::compute_l0(epoch)?;
        let stonyhurst_lon = self.longitude - l0;
        let normalized_lon = Angle::from_radians(wrap_0_2pi(stonyhurst_lon.radians())?);

        let mut stony = HeliographicStonyhurst::new(self.latitude, normalized_lon)?;
        if let Some(r) = self.radius {
            stony.set_radius(r);
        }
        Ok(stony)
    }
}

impl CoordinateFrame for HeliographicStonyhurst {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition> {
        self.to_carrington(epoch)?.to_icrs(epoch)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        HeliographicCarrington::from_icrs(icrs, epoch)?.to_stonyhurst(epoch)
    }
}

impl CoordinateFrame for HeliographicCarrington {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition> {
        let m = solar::icrs_to_carrington(epoch)?;
        let (ra, dec) = m
            .transpose()
            .transform_spherical(self.longitude.radians(), self.latitude.radians());

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
        let m = solar::icrs_to_carrington(epoch)?;
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

impl std::fmt::Display for HeliographicStonyhurst {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "HeliographicStonyhurst(lat={:.6}°, lon={:.6}°",
            self.latitude.degrees(),
            self.longitude.degrees()
        )?;

        if let Some(radius) = self.radius {
            write!(f, ", r={}", radius)?;
        }

        write!(f, ")")
    }
}

impl std::fmt::Display for HeliographicCarrington {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "HeliographicCarrington(lat={:.6}°, lon={:.6}°",
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
    use crate::aberration::compute_earth_state;
    use crate::test_support::rounded;
    use celestial_core::test_helpers::assert_ulp_le;
    use celestial_time::julian::JulianDate;

    #[test]
    fn test_stonyhurst_creation() {
        let pos = HeliographicStonyhurst::from_degrees(45.0, 30.0).unwrap();
        assert_eq!(pos.latitude(), Angle::from_degrees(45.0));
        assert_eq!(pos.longitude(), Angle::from_degrees(30.0));
        assert!(pos.radius().is_none());
    }

    #[test]
    fn test_carrington_creation() {
        let pos = HeliographicCarrington::from_degrees(-30.0, 180.0).unwrap();
        assert_eq!(pos.latitude(), Angle::from_degrees(-30.0));
        assert_eq!(pos.longitude(), Angle::from_degrees(180.0));
        assert!(pos.radius().is_none());
    }

    #[test]
    fn test_stonyhurst_validation() {
        assert!(HeliographicStonyhurst::from_degrees(0.0, 0.0).is_ok());
        assert!(HeliographicStonyhurst::from_degrees(90.0, 180.0).is_ok());
        assert!(HeliographicStonyhurst::from_degrees(-90.0, 359.0).is_ok());

        assert!(HeliographicStonyhurst::from_degrees(95.0, 0.0).is_err());
        assert!(HeliographicStonyhurst::from_degrees(-95.0, 0.0).is_err());
    }

    #[test]
    fn test_stonyhurst_to_carrington_differs_by_l0() {
        let epoch = TT::j2000();
        let stonyhurst = HeliographicStonyhurst::from_degrees(15.0, 45.0).unwrap();
        let carrington = stonyhurst.to_carrington(&epoch).unwrap();

        assert_eq!(
            stonyhurst.latitude().degrees(),
            carrington.latitude().degrees()
        );

        let l0 = solar::compute_l0(&epoch).unwrap();
        assert_eq!(
            carrington.longitude().radians(),
            wrap_0_2pi((stonyhurst.longitude() + l0).radians()).unwrap()
        );
    }

    #[test]
    fn test_carrington_to_stonyhurst_roundtrip() {
        let epoch = TT::j2000();
        let original = HeliographicCarrington::from_degrees(30.0, 120.0).unwrap();
        let stonyhurst = original.to_stonyhurst(&epoch).unwrap();
        let roundtrip = stonyhurst.to_carrington(&epoch).unwrap();

        assert_eq!(roundtrip, original);
    }

    #[test]
    fn test_disk_center() {
        let epoch = TT::j2000();
        let center = HeliographicStonyhurst::disk_center(&epoch).unwrap();

        assert_eq!(center.latitude(), solar::compute_b0(&epoch).unwrap());
        assert_eq!(center.longitude().degrees(), 0.0);
    }

    #[test]
    fn test_solar_pole_is_iau_pole_at_j2000() {
        let epoch = TT::j2000();
        let poles = [
            HeliographicStonyhurst::from_degrees(90.0, 0.0)
                .unwrap()
                .to_icrs(&epoch)
                .unwrap(),
            HeliographicCarrington::from_degrees(90.0, 0.0)
                .unwrap()
                .to_icrs(&epoch)
                .unwrap(),
        ];
        for pole in poles {
            let printed = [pole.ra(), pole.dec()].map(|a| rounded(a.degrees(), 2));
            assert_eq!(printed, [286.13, 63.87]);
        }
    }

    #[test]
    fn test_disk_centre_is_earth_direction_displaced_by_aberration() {
        let epoch = TT::from_julian_date(JulianDate::new(2448908.50068, 0.0));
        let earth = compute_earth_state(&epoch).unwrap().heliocentric_position;
        let (ra, dec) = earth.to_spherical();
        let earth_direction = ICRSPosition::new(
            Angle::from_radians(wrap_0_2pi(ra).unwrap()),
            Angle::from_radians(dec),
        )
        .unwrap();

        let centre = HeliographicStonyhurst::disk_center(&epoch)
            .unwrap()
            .to_icrs(&epoch)
            .unwrap();
        let separation = centre.angular_separation(&earth_direction).arcseconds();
        assert_eq!(
            rounded(separation, 2),
            rounded(20.4898 / earth.magnitude(), 2)
        );
    }

    #[test]
    fn test_disk_centre_maps_back_to_b0_and_l0() {
        let epoch = TT::from_julian_date(JulianDate::new(2448908.50068, 0.0));
        let centre = HeliographicStonyhurst::disk_center(&epoch)
            .unwrap()
            .to_carrington(&epoch)
            .unwrap();
        let printed = [centre.latitude(), centre.longitude()].map(|a| rounded(a.degrees(), 2));
        assert_eq!(printed, [5.99, 238.63]);
    }

    #[test]
    fn test_non_finite_epoch_is_err() {
        let epoch = TT::from_julian_date(JulianDate::new(f64::NAN, 0.0));
        let icrs = ICRSPosition::from_degrees(10.0, 20.0).unwrap();
        let stonyhurst = HeliographicStonyhurst::from_degrees(10.0, 20.0).unwrap();
        let carrington = HeliographicCarrington::from_degrees(10.0, 20.0).unwrap();
        assert!(stonyhurst.to_icrs(&epoch).is_err());
        assert!(carrington.to_icrs(&epoch).is_err());
        assert!(HeliographicStonyhurst::from_icrs(&icrs, &epoch).is_err());
        assert!(HeliographicCarrington::from_icrs(&icrs, &epoch).is_err());
    }

    #[test]
    fn test_coordinate_frame_roundtrip() {
        let epoch = TT::j2000();
        let test_cases = [
            (20.0, 30.0),
            (0.0, 0.0),
            (45.0, 90.0),
            (-7.0, 180.0),
            (7.0, 270.0),
        ];

        for (lat, lon) in test_cases {
            let original = HeliographicStonyhurst::from_degrees(lat, lon).unwrap();
            let icrs = original.to_icrs(&epoch).unwrap();
            let recovered = HeliographicStonyhurst::from_icrs(&icrs, &epoch).unwrap();

            let ctx = format!("({lat}, {lon})");
            assert_ulp_le(
                recovered.latitude().radians(),
                original.latitude().radians(),
                2,
                &ctx,
            );
            assert_ulp_le(
                recovered.longitude().radians(),
                original.longitude().radians(),
                2,
                &ctx,
            );
        }
    }

    #[test]
    fn test_with_radius() {
        let radius = Distance::from_au(0.00465047).unwrap();
        let pos = HeliographicStonyhurst::with_radius(
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
        let pos = HeliographicStonyhurst::from_degrees(45.123456, 30.654321).unwrap();
        let display = format!("{}", pos);
        assert!(display.contains("45.123456"));
        assert!(display.contains("30.654321"));
        assert!(display.contains("HeliographicStonyhurst"));
    }
}
