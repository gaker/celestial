use crate::{solar, transforms::CoordinateFrame, CoordResult, Distance, ICRSPosition};
use celestial_core::constants::HALF_PI;
use celestial_core::matrix::RotationMatrix3;
use celestial_core::utils::normalize_angle_to_positive;
use celestial_core::Angle;
use celestial_time::TT;

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
        let longitude = longitude.validate_longitude(true)?;

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

    pub fn set_radius(&mut self, radius: Distance) {
        self.radius = Some(radius);
    }

    pub fn to_carrington(&self, epoch: &TT) -> CoordResult<HeliographicCarrington> {
        let l0 = solar::compute_l0(epoch);
        let carrington_lon = self.longitude + l0;
        let normalized_lon =
            Angle::from_radians(normalize_angle_to_positive(carrington_lon.radians()));

        let mut carr = HeliographicCarrington::new(self.latitude, normalized_lon)?;
        if let Some(r) = self.radius {
            carr.set_radius(r);
        }
        Ok(carr)
    }

    pub fn disk_center(epoch: &TT) -> Self {
        let orientation = solar::compute_solar_orientation(epoch);
        Self {
            latitude: orientation.b0,
            longitude: Angle::ZERO,
            radius: None,
        }
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
        let longitude = longitude.validate_longitude(true)?;

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

    pub fn set_radius(&mut self, radius: Distance) {
        self.radius = Some(radius);
    }

    pub fn to_stonyhurst(&self, epoch: &TT) -> CoordResult<HeliographicStonyhurst> {
        let l0 = solar::compute_l0(epoch);
        let stonyhurst_lon = self.longitude - l0;
        let normalized_lon =
            Angle::from_radians(normalize_angle_to_positive(stonyhurst_lon.radians()));

        let mut stony = HeliographicStonyhurst::new(self.latitude, normalized_lon)?;
        if let Some(r) = self.radius {
            stony.set_radius(r);
        }
        Ok(stony)
    }

    pub fn carrington_rotation_number(epoch: &TT) -> f64 {
        const CARRINGTON_EPOCH_JD: f64 = 2398220.0;
        const CARRINGTON_PERIOD_DAYS: f64 = 25.38;

        let jd = epoch.to_julian_date();
        let d = jd.jd1() + jd.jd2() - CARRINGTON_EPOCH_JD;
        d / CARRINGTON_PERIOD_DAYS
    }
}

fn heliographic_to_icrs_matrix(epoch: &TT) -> CoordResult<RotationMatrix3> {
    let orientation = solar::compute_solar_orientation(epoch);
    let b0 = orientation.b0.radians();
    let p = orientation.p.radians();

    let sun_icrs = solar::get_sun_icrs(epoch)?;
    let sun_ra = sun_icrs.ra().radians();
    let sun_dec = sun_icrs.dec().radians();

    let mut m = RotationMatrix3::identity();
    m.rotate_y(-b0);
    m.rotate_z(p);
    m.rotate_y(sun_dec - HALF_PI);
    m.rotate_z(-sun_ra);
    Ok(m)
}

impl CoordinateFrame for HeliographicStonyhurst {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition> {
        let m = heliographic_to_icrs_matrix(epoch)?;
        let (ra, dec) = m
            .transpose()
            .transform_spherical(self.longitude.radians(), self.latitude.radians());

        let mut icrs = ICRSPosition::new(
            Angle::from_radians(normalize_angle_to_positive(ra)),
            Angle::from_radians(dec),
        )?;

        if let Some(radius) = self.radius {
            icrs.set_distance(radius);
        }
        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let m = heliographic_to_icrs_matrix(epoch)?;
        let (lon, lat) = m.transform_spherical(icrs.ra().radians(), icrs.dec().radians());

        let mut pos = Self::new(
            Angle::from_radians(lat),
            Angle::from_radians(normalize_angle_to_positive(lon)),
        )?;

        if let Some(dist) = icrs.distance() {
            pos.set_radius(dist);
        }
        Ok(pos)
    }
}

impl CoordinateFrame for HeliographicCarrington {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition> {
        let stonyhurst = self.to_stonyhurst(epoch)?;
        stonyhurst.to_icrs(epoch)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let stonyhurst = HeliographicStonyhurst::from_icrs(icrs, epoch)?;
        stonyhurst.to_carrington(epoch)
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

    #[test]
    fn test_stonyhurst_creation() {
        let pos = HeliographicStonyhurst::from_degrees(45.0, 30.0).unwrap();
        assert!((pos.latitude().degrees() - 45.0).abs() < 1e-12);
        assert!((pos.longitude().degrees() - 30.0).abs() < 1e-12);
        assert!(pos.radius().is_none());
    }

    #[test]
    fn test_carrington_creation() {
        let pos = HeliographicCarrington::from_degrees(-30.0, 180.0).unwrap();
        assert!((pos.latitude().degrees() - (-30.0)).abs() < 1e-12);
        assert!((pos.longitude().degrees() - 180.0).abs() < 1e-12);
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

        let l0 = solar::compute_l0(&epoch);
        let expected_carr_lon =
            normalize_angle_to_positive((stonyhurst.longitude() + l0).radians())
                * celestial_core::constants::RAD_TO_DEG;

        assert!((carrington.longitude().degrees() - expected_carr_lon).abs() < 1e-10);
    }

    #[test]
    fn test_carrington_to_stonyhurst_roundtrip() {
        let epoch = TT::j2000();
        let original = HeliographicCarrington::from_degrees(30.0, 120.0).unwrap();
        let stonyhurst = original.to_stonyhurst(&epoch).unwrap();
        let roundtrip = stonyhurst.to_carrington(&epoch).unwrap();

        assert!((original.latitude().degrees() - roundtrip.latitude().degrees()).abs() < 1e-10);
        assert!((original.longitude().degrees() - roundtrip.longitude().degrees()).abs() < 1e-10);
    }

    #[test]
    fn test_disk_center() {
        let epoch = TT::j2000();
        let center = HeliographicStonyhurst::disk_center(&epoch);

        let b0 = solar::compute_b0(&epoch);
        assert!((center.latitude().degrees() - b0.degrees()).abs() < 1e-12);
        assert_eq!(center.longitude().degrees(), 0.0);
    }

    #[test]
    fn test_carrington_rotation_number() {
        let epoch = TT::j2000();
        let rotation = HeliographicCarrington::carrington_rotation_number(&epoch);

        assert!(
            rotation > 1900.0 && rotation < 2200.0,
            "Carrington rotation number at J2000 = {} should be reasonable",
            rotation
        );
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

            let lat_err = (original.latitude().degrees() - recovered.latitude().degrees()).abs();
            let lon_diff = (original.longitude().radians() - recovered.longitude().radians()).abs();
            let lon_err = if lon_diff > std::f64::consts::PI {
                std::f64::consts::TAU - lon_diff
            } else {
                lon_diff
            } * celestial_core::constants::RAD_TO_DEG;

            assert!(
                lat_err < 1.0 / 3600.0,
                "({}, {}): Latitude error {:.6} arcsec",
                lat,
                lon,
                lat_err * 3600.0,
            );
            assert!(
                lon_err < 1.0 / 3600.0,
                "({}, {}): Longitude error {:.6} arcsec",
                lat,
                lon,
                lon_err * 3600.0,
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
