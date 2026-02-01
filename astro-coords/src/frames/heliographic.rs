use crate::{solar, transforms::CoordinateFrame, CoordResult, Distance, ICRSPosition};
use astro_core::utils::normalize_angle_to_positive;
use astro_core::Angle;
use astro_time::TT;

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

impl CoordinateFrame for HeliographicStonyhurst {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition> {
        let orientation = solar::compute_solar_orientation(epoch);
        let b0 = orientation.b0.radians();
        let p = orientation.p.radians();

        let lat = self.latitude.radians();
        let lon = self.longitude.radians();
        let (sin_lat, cos_lat) = lat.sin_cos();
        let (sin_lon, cos_lon) = lon.sin_cos();

        let x_helio = cos_lat * cos_lon;
        let y_helio = cos_lat * sin_lon;
        let z_helio = sin_lat;

        let (sin_b0, cos_b0) = b0.sin_cos();
        let _x_proj = x_helio;
        let y_proj = y_helio * cos_b0 - z_helio * sin_b0;
        let z_proj = y_helio * sin_b0 + z_helio * cos_b0;

        let (sin_p, cos_p) = p.sin_cos();
        let y_disk = y_proj * cos_p + z_proj * sin_p;
        let z_disk = -y_proj * sin_p + z_proj * cos_p;

        let sun_icrs = solar::get_sun_icrs(epoch)?;
        let sun_ra = sun_icrs.ra().radians();
        let sun_dec = sun_icrs.dec().radians();

        let scale = 0.00465047; // ~1 solar radius in AU
        let offset_x = -y_disk * scale;
        let offset_y = z_disk * scale;

        let (_sin_dec, cos_dec) = sun_dec.sin_cos();
        let (_sin_ra, _cos_ra) = sun_ra.sin_cos();

        let ra_offset = offset_x / cos_dec;
        let dec_offset = offset_y;

        let new_ra = sun_ra + ra_offset;
        let new_dec = sun_dec + dec_offset;

        let mut icrs = ICRSPosition::new(
            Angle::from_radians(normalize_angle_to_positive(new_ra)),
            Angle::from_radians(new_dec),
        )?;

        if let Some(radius) = self.radius {
            icrs.set_distance(radius);
        }

        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let orientation = solar::compute_solar_orientation(epoch);
        let b0 = orientation.b0.radians();
        let p = orientation.p.radians();

        let sun_icrs = solar::get_sun_icrs(epoch)?;
        let sun_ra = sun_icrs.ra().radians();
        let sun_dec = sun_icrs.dec().radians();

        let ra_offset = icrs.ra().radians() - sun_ra;
        let dec_offset = icrs.dec().radians() - sun_dec;

        let scale = 0.00465047;
        let y_disk = -ra_offset * sun_dec.cos() / scale;
        let z_disk = dec_offset / scale;

        let (sin_p, cos_p) = p.sin_cos();
        let y_proj = y_disk * cos_p - z_disk * sin_p;
        let z_proj = y_disk * sin_p + z_disk * cos_p;

        let (sin_b0, cos_b0) = b0.sin_cos();
        let y_helio = y_proj * cos_b0 + z_proj * sin_b0;
        let z_helio = -y_proj * sin_b0 + z_proj * cos_b0;

        let x_helio = (1.0 - y_helio.powi(2) - z_helio.powi(2)).sqrt().max(0.0);

        let lat = z_helio.asin();
        let lon = y_helio.atan2(x_helio);

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
                * astro_core::constants::RAD_TO_DEG;

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
        let original = HeliographicStonyhurst::from_degrees(20.0, 30.0).unwrap();

        let icrs = original.to_icrs(&epoch).unwrap();
        let recovered = HeliographicStonyhurst::from_icrs(&icrs, &epoch).unwrap();

        assert!(
            (original.latitude().degrees() - recovered.latitude().degrees()).abs() < 5.0,
            "Latitude mismatch: {} vs {}",
            original.latitude().degrees(),
            recovered.latitude().degrees()
        );
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
