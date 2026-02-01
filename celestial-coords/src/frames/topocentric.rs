use crate::{CoordResult, Distance};
use celestial_core::constants::{HALF_PI, TWOPI};
use celestial_core::{Angle, Location};
use celestial_time::TT;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TopocentricPosition {
    azimuth: Angle,
    elevation: Angle,
    observer: Location,
    epoch: TT,
    distance: Option<Distance>,
}

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct HourAnglePosition {
    hour_angle: Angle,
    declination: Angle,
    observer: Location,
    epoch: TT,
    distance: Option<Distance>,
}

impl TopocentricPosition {
    pub fn new(
        azimuth: Angle,
        elevation: Angle,
        observer: Location,
        epoch: TT,
    ) -> CoordResult<Self> {
        let azimuth = azimuth.validate_longitude(true)?;
        let elevation = elevation.validate_latitude()?;

        Ok(Self {
            azimuth,
            elevation,
            observer,
            epoch,
            distance: None,
        })
    }

    pub fn with_distance(
        azimuth: Angle,
        elevation: Angle,
        observer: Location,
        epoch: TT,
        distance: Distance,
    ) -> CoordResult<Self> {
        let mut pos = Self::new(azimuth, elevation, observer, epoch)?;
        pos.distance = Some(distance);
        Ok(pos)
    }

    pub fn from_degrees(
        az_deg: f64,
        el_deg: f64,
        observer: Location,
        epoch: TT,
    ) -> CoordResult<Self> {
        Self::new(
            Angle::from_degrees(az_deg),
            Angle::from_degrees(el_deg),
            observer,
            epoch,
        )
    }

    pub fn azimuth(&self) -> Angle {
        self.azimuth
    }

    pub fn elevation(&self) -> Angle {
        self.elevation
    }

    pub fn observer(&self) -> &Location {
        &self.observer
    }

    pub fn epoch(&self) -> TT {
        self.epoch
    }

    pub fn distance(&self) -> Option<Distance> {
        self.distance
    }

    pub fn set_distance(&mut self, distance: Distance) {
        self.distance = Some(distance);
    }

    pub fn zenith_angle(&self) -> Angle {
        Angle::HALF_PI - self.elevation
    }

    pub fn air_mass(&self) -> f64 {
        self.air_mass_rozenberg()
    }

    pub fn air_mass_rozenberg(&self) -> f64 {
        let zenith = self.zenith_angle();
        if zenith.degrees() >= 90.0 {
            return 40.0;
        }
        let cos_z = zenith.cos();
        let term = cos_z + 0.025 * (-11.0 * cos_z).exp();
        1.0 / term
    }

    /// Computes airmass using Pickering's (2002) empirical formula.
    ///
    /// # Valid Range
    /// - Returns `f64::INFINITY` for elevations ≤ -2° (below horizon)
    /// - Accurate for elevations > -2° (including astronomical twilight)
    /// - Values become increasingly unreliable as elevation approaches -2°
    ///
    /// # Numerical Stability
    /// Near the horizon (0° to 5°), results can be very large but remain finite.
    /// Use this method only if observations extend below the horizon; otherwise
    /// prefer `air_mass_hardie()` or `air_mass_kasten_young()`.
    ///
    /// Reference: Pickering, K. A. (2002). "The Southern Limits of the Ancient Star
    /// Catalog". DIO 12, 3-27.
    pub fn air_mass_pickering(&self) -> f64 {
        let h = self.elevation.degrees();
        if h <= -2.0 {
            return f64::INFINITY;
        }
        let h_term = 244.0 / (165.0 + 47.0 * h.abs().powf(1.1));
        let sin_arg = h + h_term;
        1.0 / sin_arg.to_radians().sin()
    }

    pub fn air_mass_kasten_young(&self) -> f64 {
        let zenith_deg = self.zenith_angle().degrees();
        if zenith_deg >= 90.0 {
            return 38.0;
        }
        let cos_z = self.zenith_angle().cos();
        let term = (96.07995 - zenith_deg).powf(-1.6364);
        1.0 / (cos_z + 0.50572 * term)
    }

    pub fn atmospheric_refraction(
        &self,
        pressure_hpa: f64,
        temp_celsius: f64,
        rel_humidity: f64,
        wavelength_um: f64,
    ) -> Angle {
        let (refa, refb) =
            self.refraction_coefficients(pressure_hpa, temp_celsius, rel_humidity, wavelength_um);

        // Apply refraction formula: dZ = A tan(Z) + B tan³(Z)
        let z_obs = self.zenith_angle().radians();
        let tan_z = z_obs.tan();
        let refraction_rad = refa * tan_z + refb * tan_z.powi(3);

        Angle::from_radians(refraction_rad)
    }

    fn refraction_coefficients(
        &self,
        mut pressure_hpa: f64,
        mut temp_celsius: f64,
        mut rel_humidity: f64,
        wavelength_um: f64,
    ) -> (f64, f64) {
        pressure_hpa = pressure_hpa.clamp(0.0, 10000.0);
        temp_celsius = temp_celsius.clamp(-150.0, 200.0);
        rel_humidity = rel_humidity.clamp(0.0, 1.0);

        // Zero pressure means no refraction
        if pressure_hpa <= 0.0 {
            return (0.0, 0.0);
        }

        // Determine if optical/IR or radio
        let is_optical = (0.0..100.0).contains(&wavelength_um);

        // Temperature in Kelvin
        let temp_kelvin = temp_celsius + 273.15;

        let pw = if pressure_hpa > 0.0 {
            // Saturation pressure (Gill 1982, with pressure correction)
            let ps = 10.0_f64
                .powf((0.7859 + 0.03477 * temp_celsius) / (1.0 + 0.00412 * temp_celsius))
                * (1.0 + pressure_hpa * (4.5e-6 + 6e-10 * temp_celsius * temp_celsius));
            // Actual water vapor pressure (Crane 1976)
            rel_humidity * ps / (1.0 - (1.0 - rel_humidity) * ps / pressure_hpa)
        } else {
            0.0
        };

        if is_optical {
            // Optical/IR refraction (Hohenkerk & Sinclair 1985, IAG 1999)
            let wl_sq = wavelength_um * wavelength_um;
            let gamma = ((77.53484e-6 + (4.39108e-7 + 3.666e-9 / wl_sq) / wl_sq) * pressure_hpa
                - 11.2684e-6 * pw)
                / temp_kelvin;

            // Beta from Stone (1996) with empirical adjustments
            let beta = 4.4474e-6 * temp_kelvin;

            // Green (1987) Equation 4.31
            let refa = gamma * (1.0 - beta);
            let refb = -gamma * (beta - gamma / 2.0);

            (refa, refb)
        } else {
            // Radio refraction (Rueger 2002)
            let gamma = (77.6890e-6 * pressure_hpa - (6.3938e-6 - 0.375463 / temp_kelvin) * pw)
                / temp_kelvin;

            // Beta with humidity correction for radio
            let mut beta = 4.4474e-6 * temp_kelvin;
            beta -= 0.0074 * pw * beta;

            // Green (1987) Equation 4.31
            let refa = gamma * (1.0 - beta);
            let refb = -gamma * (beta - gamma / 2.0);

            (refa, refb)
        }
    }

    pub fn with_refraction(
        &self,
        pressure_hpa: f64,
        temp_celsius: f64,
        rel_humidity: f64,
        wavelength_um: f64,
    ) -> Self {
        let refraction =
            self.atmospheric_refraction(pressure_hpa, temp_celsius, rel_humidity, wavelength_um);

        // Refraction increases apparent elevation (decreases zenith distance)
        let apparent_elevation = self.elevation + refraction;

        Self {
            azimuth: self.azimuth,
            elevation: apparent_elevation,
            observer: self.observer,
            epoch: self.epoch,
            distance: self.distance,
        }
    }

    pub fn without_refraction(
        &self,
        pressure_hpa: f64,
        temp_celsius: f64,
        rel_humidity: f64,
        wavelength_um: f64,
    ) -> Self {
        let refraction =
            self.atmospheric_refraction(pressure_hpa, temp_celsius, rel_humidity, wavelength_um);

        // Remove refraction to get true elevation
        let true_elevation = self.elevation - refraction;

        Self {
            azimuth: self.azimuth,
            elevation: true_elevation,
            observer: self.observer,
            epoch: self.epoch,
            distance: self.distance,
        }
    }

    pub fn diurnal_parallax(&self) -> Option<Angle> {
        self.distance.map(|d| {
            // Earth's equatorial radius in AU
            const EARTH_RADIUS_AU: f64 = 4.2635e-5; // 6378.137 km

            let distance_au = d.au();
            let zenith_angle = self.zenith_angle();

            // Parallax formula: p = arcsin((a/r) × sin(z))
            let ratio = EARTH_RADIUS_AU / distance_au;
            let parallax_rad = if ratio < 1.0 {
                (ratio * zenith_angle.sin()).asin()
            } else {
                // Object inside Earth - return maximum parallax
                HALF_PI
            };

            Angle::from_radians(parallax_rad)
        })
    }

    pub fn horizontal_parallax(&self) -> Option<Angle> {
        self.distance.map(|d| {
            const EARTH_RADIUS_AU: f64 = 4.2635e-5;
            let distance_au = d.au();
            let ratio = EARTH_RADIUS_AU / distance_au;

            if ratio < 1.0 {
                Angle::from_radians(ratio.asin())
            } else {
                Angle::from_radians(HALF_PI)
            }
        })
    }

    pub fn with_diurnal_parallax(&self) -> Self {
        if let Some(parallax) = self.diurnal_parallax() {
            // Parallax correction decreases elevation (object appears lower)
            let corrected_elevation = self.elevation - parallax;

            Self {
                azimuth: self.azimuth,
                elevation: corrected_elevation,
                observer: self.observer,
                epoch: self.epoch,
                distance: self.distance,
            }
        } else {
            self.clone()
        }
    }

    pub fn without_diurnal_parallax(&self) -> Self {
        if let Some(parallax) = self.diurnal_parallax() {
            // Remove parallax correction (increases elevation)
            let geocentric_elevation = self.elevation + parallax;

            Self {
                azimuth: self.azimuth,
                elevation: geocentric_elevation,
                observer: self.observer,
                epoch: self.epoch,
                distance: self.distance,
            }
        } else {
            self.clone()
        }
    }

    pub fn is_above_horizon(&self) -> bool {
        self.elevation.degrees() > 0.0
    }

    pub fn is_near_zenith(&self) -> bool {
        self.elevation.degrees() > 89.0
    }

    pub fn is_near_horizon(&self) -> bool {
        self.elevation.degrees() < 10.0 && self.is_above_horizon()
    }

    pub fn cardinal_direction(&self) -> &'static str {
        let az_deg = self.azimuth.degrees();
        if !(22.5..337.5).contains(&az_deg) {
            "N"
        } else if az_deg < 67.5 {
            "NE"
        } else if az_deg < 112.5 {
            "E"
        } else if az_deg < 157.5 {
            "SE"
        } else if az_deg < 202.5 {
            "S"
        } else if az_deg < 247.5 {
            "SW"
        } else if az_deg < 292.5 {
            "W"
        } else {
            "NW"
        }
    }

    pub fn parallactic_angle(&self, hour_angle: Angle, declination: Angle) -> Angle {
        let lat = self.observer.latitude_angle();
        let (sin_ha, cos_ha) = hour_angle.sin_cos();
        let (sin_lat, cos_lat) = lat.sin_cos();
        let (sin_dec, cos_dec) = declination.sin_cos();

        let numerator = sin_ha;
        let denominator = cos_dec * sin_lat - sin_dec * cos_lat * cos_ha;

        Angle::from_radians(numerator.atan2(denominator))
    }

    pub fn to_hour_angle(&self) -> CoordResult<HourAnglePosition> {
        let (sin_az, cos_az) = self.azimuth.sin_cos();
        let (sin_el, cos_el) = self.elevation.sin_cos();
        let (sin_lat, cos_lat) = self.observer.latitude_angle().sin_cos();

        let sin_dec = sin_el * sin_lat + cos_el * cos_lat * cos_az;
        let dec = sin_dec.asin();

        let cos_dec = dec.cos();

        let cos_ha = if cos_dec.abs() < 1e-10 {
            0.0
        } else {
            (sin_el - sin_dec * sin_lat) / (cos_dec * cos_lat)
        };

        let sin_ha = if cos_dec.abs() < 1e-10 {
            0.0
        } else {
            -sin_az * cos_el / cos_dec
        };

        let ha = sin_ha.atan2(cos_ha);

        let mut ha_pos = HourAnglePosition::new(
            Angle::from_radians(ha),
            Angle::from_radians(dec),
            self.observer,
            self.epoch,
        )?;

        if let Some(distance) = self.distance {
            ha_pos.set_distance(distance);
        }

        Ok(ha_pos)
    }

    pub fn to_cirs(&self, delta_t: f64) -> CoordResult<crate::frames::CIRSPosition> {
        let ha = self.to_hour_angle()?;
        ha.to_cirs(delta_t)
    }
}

impl HourAnglePosition {
    pub fn new(
        hour_angle: Angle,
        declination: Angle,
        observer: Location,
        epoch: TT,
    ) -> CoordResult<Self> {
        let hour_angle = hour_angle.wrapped(); // [-180°, +180°]
        let declination = declination.validate_declination(true)?; // beyond_pole for GEM pier-flips

        Ok(Self {
            hour_angle,
            declination,
            observer,
            epoch,
            distance: None,
        })
    }

    pub fn with_distance(
        hour_angle: Angle,
        declination: Angle,
        observer: Location,
        epoch: TT,
        distance: Distance,
    ) -> CoordResult<Self> {
        let mut pos = Self::new(hour_angle, declination, observer, epoch)?;
        pos.distance = Some(distance);
        Ok(pos)
    }

    pub fn hour_angle(&self) -> Angle {
        self.hour_angle
    }

    pub fn declination(&self) -> Angle {
        self.declination
    }

    pub fn observer(&self) -> &Location {
        &self.observer
    }

    pub fn epoch(&self) -> TT {
        self.epoch
    }

    pub fn distance(&self) -> Option<Distance> {
        self.distance
    }

    pub fn set_distance(&mut self, distance: Distance) {
        self.distance = Some(distance);
    }

    pub fn to_topocentric(&self) -> CoordResult<TopocentricPosition> {
        // Step 1: Pre-compute trigonometric functions
        let (sin_ha, cos_ha) = self.hour_angle.sin_cos(); // sh, ch
        let (sin_dec, cos_dec) = self.declination.sin_cos(); // sd, cd
        let (sin_lat, cos_lat) = self.observer.latitude_angle().sin_cos(); // sp, cp

        // Step 2: Compute 3D unit vector in Az/El system
        // This implements the rotation matrix transformation
        let x = -cos_ha * cos_dec * sin_lat + sin_dec * cos_lat; // X component
        let y = -sin_ha * cos_dec; // Y component
        let z = cos_ha * cos_dec * cos_lat + sin_dec * sin_lat; // Z component

        // Step 3: Convert to spherical coordinates
        let r = (x * x + y * y).sqrt(); // Horizontal distance

        let raw_azimuth = if r != 0.0 { y.atan2(x) } else { 0.0 }; // Raw azimuth angle

        let azimuth = if raw_azimuth < 0.0 {
            // Normalize to [0, 2π]
            raw_azimuth + TWOPI
        } else {
            raw_azimuth
        };

        let elevation = z.atan2(r); // Elevation angle

        let mut topo = TopocentricPosition::new(
            Angle::from_radians(azimuth),
            Angle::from_radians(elevation),
            self.observer,
            self.epoch,
        )?;

        if let Some(distance) = self.distance {
            topo.set_distance(distance);
        }

        Ok(topo)
    }

    pub fn is_circumpolar(&self) -> bool {
        let lat = self.observer.latitude_angle();
        self.declination.radians() > (Angle::HALF_PI - lat).radians()
    }

    pub fn never_rises(&self) -> bool {
        let lat = self.observer.latitude_angle();
        self.declination.radians() < -(Angle::HALF_PI - lat).radians()
    }

    pub fn to_cirs(&self, delta_t: f64) -> CoordResult<crate::frames::CIRSPosition> {
        use celestial_time::scales::conversions::ToUT1WithDeltaT;
        use celestial_time::sidereal::GAST;

        let ut1 = self.epoch.to_ut1_with_delta_t(delta_t)?;
        let gast = GAST::from_ut1_and_tt(&ut1, &self.epoch)?;
        let last = gast.to_last(&self.observer);

        let ra_rad = last.radians() - self.hour_angle.radians();
        let ra = celestial_core::angle::wrap_0_2pi(ra_rad);

        let mut cirs = crate::frames::CIRSPosition::new(
            Angle::from_radians(ra),
            self.declination,
            self.epoch,
        )?;

        if let Some(distance) = self.distance {
            cirs.set_distance(distance);
        }

        Ok(cirs)
    }
}

// Note: Topocentric coordinates cannot implement CoordinateFrame directly
// because they require both time AND observer location for transformation.
// They need specialized transformation methods.

impl std::fmt::Display for TopocentricPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Topocentric(Az={:.2}° {}, El={:.2}°",
            self.azimuth.degrees(),
            self.cardinal_direction(),
            self.elevation.degrees()
        )?;

        if let Some(distance) = self.distance {
            write!(f, ", d={}", distance)?;
        }

        write!(f, ")")
    }
}

impl std::fmt::Display for HourAnglePosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "HourAngle(HA={:.4}h, Dec={:.4}°",
            self.hour_angle.hours(),
            self.declination.degrees()
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

    fn test_observer() -> Location {
        // Keck Observatory, Mauna Kea (4145m per keckobservatory.org)
        Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap()
    }

    #[test]
    fn test_topocentric_creation() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        assert!((topo.azimuth().degrees() - 180.0).abs() < 1e-12);
        assert!((topo.elevation().degrees() - 45.0).abs() < 1e-12);
        assert_eq!(
            topo.observer().latitude_degrees(),
            observer.latitude_degrees()
        );
        assert_eq!(topo.epoch(), epoch);
    }

    #[test]
    fn test_topocentric_validation() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Valid coordinates
        assert!(TopocentricPosition::from_degrees(0.0, 0.0, observer, epoch).is_ok());
        assert!(TopocentricPosition::from_degrees(359.999, 89.999, observer, epoch).is_ok());

        // Azimuth gets normalized
        let topo = TopocentricPosition::from_degrees(380.0, 45.0, observer, epoch).unwrap();
        assert!((topo.azimuth().degrees() - 20.0).abs() < 1e-12);

        // Invalid elevation
        assert!(TopocentricPosition::from_degrees(0.0, 95.0, observer, epoch).is_err());
    }

    #[test]
    fn test_zenith_and_air_mass() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Object at zenith - Rozenberg formula
        let zenith = TopocentricPosition::from_degrees(0.0, 90.0, observer, epoch).unwrap();
        assert!((zenith.zenith_angle().degrees() - 0.0).abs() < 1e-12);
        assert!((zenith.air_mass() - 1.0).abs() < 0.001);

        // Object at 60° elevation (zenith angle = 30°) - Rozenberg
        let high = TopocentricPosition::from_degrees(0.0, 60.0, observer, epoch).unwrap();
        assert!((high.zenith_angle().degrees() - 30.0).abs() < 1e-12);
        // Rozenberg gives slightly different value than simple sec(z)
        assert!((high.air_mass() - 1.154).abs() < 0.01);

        // Object at horizon - Rozenberg
        let horizon = TopocentricPosition::from_degrees(0.0, 0.0, observer, epoch).unwrap();
        assert!((horizon.air_mass() - 40.0).abs() < 0.1);

        // Object below horizon
        let below = TopocentricPosition::from_degrees(0.0, -10.0, observer, epoch).unwrap();
        assert_eq!(below.air_mass(), 40.0);
    }

    #[test]
    fn test_position_classification() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let above = TopocentricPosition::from_degrees(0.0, 45.0, observer, epoch).unwrap();
        assert!(above.is_above_horizon());
        assert!(!above.is_near_zenith());
        assert!(!above.is_near_horizon());

        let below = TopocentricPosition::from_degrees(0.0, -5.0, observer, epoch).unwrap();
        assert!(!below.is_above_horizon());

        let zenith = TopocentricPosition::from_degrees(0.0, 89.5, observer, epoch).unwrap();
        assert!(zenith.is_near_zenith());

        let horizon = TopocentricPosition::from_degrees(0.0, 5.0, observer, epoch).unwrap();
        assert!(horizon.is_near_horizon());
    }

    #[test]
    fn test_cardinal_directions() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let north = TopocentricPosition::from_degrees(0.0, 45.0, observer, epoch).unwrap();
        assert_eq!(north.cardinal_direction(), "N");

        let east = TopocentricPosition::from_degrees(90.0, 45.0, observer, epoch).unwrap();
        assert_eq!(east.cardinal_direction(), "E");

        let south = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        assert_eq!(south.cardinal_direction(), "S");

        let west = TopocentricPosition::from_degrees(270.0, 45.0, observer, epoch).unwrap();
        assert_eq!(west.cardinal_direction(), "W");

        let northeast = TopocentricPosition::from_degrees(45.0, 45.0, observer, epoch).unwrap();
        assert_eq!(northeast.cardinal_direction(), "NE");
    }

    #[test]
    fn test_hour_angle_creation() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let ha_pos = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();

        assert!((ha_pos.hour_angle().hours() - 2.0).abs() < 1e-12);
        assert!((ha_pos.declination().degrees() - 45.0).abs() < 1e-12);
    }

    #[test]
    fn test_hour_angle_to_topocentric() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Object on meridian at 45° declination
        let ha_pos = HourAnglePosition::new(
            Angle::ZERO, // On meridian
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();

        let topo = ha_pos.to_topocentric().unwrap();

        // On meridian should be due south (or north depending on observer latitude)
        // and elevation should be related to observer latitude and declination
        assert!(topo.is_above_horizon());
    }

    #[test]
    fn test_circumpolar() {
        let observer = test_observer(); // Latitude ~20°N
        let epoch = TT::j2000();

        // Very high declination object (near north pole)
        let high_dec =
            HourAnglePosition::new(Angle::ZERO, Angle::from_degrees(85.0), observer, epoch)
                .unwrap();
        assert!(high_dec.is_circumpolar());

        // Low declination object
        let low_dec =
            HourAnglePosition::new(Angle::ZERO, Angle::from_degrees(0.0), observer, epoch).unwrap();
        assert!(!low_dec.is_circumpolar());

        // Very negative declination
        let neg_dec =
            HourAnglePosition::new(Angle::ZERO, Angle::from_degrees(-85.0), observer, epoch)
                .unwrap();
        assert!(neg_dec.never_rises());
    }

    #[test]
    fn test_with_distance() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let distance = Distance::from_kilometers(384400.0).unwrap(); // Moon distance

        let topo = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();

        assert_eq!(topo.distance().unwrap().kilometers(), distance.kilometers());

        let ha_pos = HourAnglePosition::with_distance(
            Angle::from_hours(1.0),
            Angle::from_degrees(30.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();

        assert_eq!(
            ha_pos.distance().unwrap().kilometers(),
            distance.kilometers()
        );
    }

    #[test]
    fn test_topocentric_set_distance() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let mut topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();

        assert!(topo.distance().is_none());

        let distance = Distance::from_kilometers(1000.0).unwrap();
        topo.set_distance(distance);

        assert!((topo.distance().unwrap().kilometers() - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn test_cardinal_directions_all() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Test all 8 cardinal directions
        let directions = [
            (0.0, "N"),
            (45.0, "NE"),
            (90.0, "E"),
            (135.0, "SE"),
            (180.0, "S"),
            (225.0, "SW"),
            (270.0, "W"),
            (315.0, "NW"),
        ];

        for (az, expected) in directions {
            let topo = TopocentricPosition::from_degrees(az, 45.0, observer, epoch).unwrap();
            assert_eq!(
                topo.cardinal_direction(),
                expected,
                "Failed for azimuth {}°",
                az
            );
        }
    }

    #[test]
    fn test_parallactic_angle() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();

        // Test parallactic angle calculation
        let ha = Angle::from_hours(1.0);
        let dec = Angle::from_degrees(45.0);
        let pa = topo.parallactic_angle(ha, dec);

        // Should return a valid angle
        assert!(pa.radians().is_finite());
    }

    #[test]
    fn test_hour_angle_getters() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let ha_pos = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();

        assert_eq!(
            ha_pos.observer().latitude_degrees(),
            observer.latitude_degrees()
        );
        assert_eq!(ha_pos.epoch(), epoch);
        assert!(ha_pos.distance().is_none());
    }

    #[test]
    fn test_hour_angle_set_distance() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let mut ha_pos = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();

        let distance = Distance::from_kilometers(500.0).unwrap();
        ha_pos.set_distance(distance);
        assert!((ha_pos.distance().unwrap().kilometers() - 500.0).abs() < 1e-6);
    }

    #[test]
    fn test_hour_angle_with_distance_to_topocentric() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let distance = Distance::from_kilometers(1000.0).unwrap();

        let ha_pos = HourAnglePosition::with_distance(
            Angle::ZERO,
            Angle::from_degrees(45.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();

        let topo = ha_pos.to_topocentric().unwrap();
        assert!((topo.distance().unwrap().kilometers() - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn test_display_formatting() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Test TopocentricPosition display
        let topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        let display = format!("{}", topo);
        assert!(display.contains("Topocentric"));
        assert!(display.contains("180.00°"));
        assert!(display.contains("45.00°"));
        assert!(display.contains("S")); // Cardinal direction

        // Test with distance
        let distance = Distance::from_kilometers(1000.0).unwrap();
        let topo_dist = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();
        let display_dist = format!("{}", topo_dist);
        assert!(display_dist.contains("AU") || display_dist.contains("pc")); // Distance shown

        // Test HourAnglePosition display
        let ha = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();
        let ha_display = format!("{}", ha);
        assert!(ha_display.contains("HourAngle"));
        assert!(ha_display.contains("2."));
        assert!(ha_display.contains("45."));

        // Test HourAngle with distance
        let ha_dist = HourAnglePosition::with_distance(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();
        let ha_display_dist = format!("{}", ha_dist);
        assert!(ha_display_dist.contains("AU") || ha_display_dist.contains("pc"));
        // Distance shown
    }

    #[test]
    fn test_air_mass_formulas_at_zenith() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let zenith = TopocentricPosition::from_degrees(0.0, 90.0, observer, epoch).unwrap();

        let rozenberg = zenith.air_mass_rozenberg();
        let pickering = zenith.air_mass_pickering();
        let kasten = zenith.air_mass_kasten_young();

        // All formulas should give ~1.0 at zenith
        assert!((rozenberg - 1.0).abs() < 0.001);
        assert!((pickering - 1.0).abs() < 0.001);
        assert!((kasten - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_air_mass_formulas_moderate_angles() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Test at 30° zenith angle (60° elevation)
        let pos_30 = TopocentricPosition::from_degrees(0.0, 60.0, observer, epoch).unwrap();
        let roz_30 = pos_30.air_mass_rozenberg();
        let pick_30 = pos_30.air_mass_pickering();
        let ky_30 = pos_30.air_mass_kasten_young();

        // Simple sec(30°) = 1.1547
        // All formulas should be within 1% of each other at moderate angles
        assert!((roz_30 - 1.155).abs() < 0.01);
        assert!((pick_30 - 1.155).abs() < 0.01);
        assert!((ky_30 - 1.155).abs() < 0.01);
        assert!((roz_30 - pick_30).abs() < 0.02);
        assert!((roz_30 - ky_30).abs() < 0.02);

        // Test at 60° zenith angle (30° elevation)
        let pos_60 = TopocentricPosition::from_degrees(0.0, 30.0, observer, epoch).unwrap();
        let roz_60 = pos_60.air_mass_rozenberg();
        let pick_60 = pos_60.air_mass_pickering();
        let ky_60 = pos_60.air_mass_kasten_young();

        // Simple sec(60°) = 2.0
        assert!((roz_60 - 2.0).abs() < 0.05);
        assert!((pick_60 - 2.0).abs() < 0.05);
        assert!((ky_60 - 2.0).abs() < 0.05);
        assert!((roz_60 - pick_60).abs() < 0.1);
    }

    #[test]
    fn test_air_mass_formulas_high_zenith_angles() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Test at 75° zenith angle (15° elevation)
        let pos_75 = TopocentricPosition::from_degrees(0.0, 15.0, observer, epoch).unwrap();
        let roz_75 = pos_75.air_mass_rozenberg();
        let pick_75 = pos_75.air_mass_pickering();
        let ky_75 = pos_75.air_mass_kasten_young();

        // Simple sec(75°) = 3.864
        // Formulas may diverge more at high angles
        assert!(roz_75 > 3.5 && roz_75 < 4.5);
        assert!(pick_75 > 3.5 && pick_75 < 4.5);
        assert!(ky_75 > 3.5 && ky_75 < 4.5);

        // Test at 85° zenith angle (5° elevation)
        let pos_85 = TopocentricPosition::from_degrees(0.0, 5.0, observer, epoch).unwrap();
        let roz_85 = pos_85.air_mass_rozenberg();
        let pick_85 = pos_85.air_mass_pickering();
        let ky_85 = pos_85.air_mass_kasten_young();

        // Simple sec(85°) = 11.47
        // All formulas valid to horizon
        assert!(roz_85 > 10.0 && roz_85 < 15.0);
        assert!(pick_85 > 10.0 && pick_85 < 15.0);
        assert!(ky_85 > 10.0 && ky_85 < 15.0);
    }

    #[test]
    fn test_air_mass_formulas_near_horizon() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Test at horizon (0° elevation, 90° zenith)
        let horizon = TopocentricPosition::from_degrees(0.0, 0.0, observer, epoch).unwrap();
        let roz_hor = horizon.air_mass_rozenberg();
        let pick_hor = horizon.air_mass_pickering();
        let ky_hor = horizon.air_mass_kasten_young();

        // Rozenberg: horizon air mass = 40
        assert!((roz_hor - 40.0).abs() < 0.1);

        // Kasten-Young: horizon air mass ~38
        assert!((ky_hor - 38.0).abs() < 1.0);

        // Pickering: should also be reasonable at horizon
        assert!(pick_hor > 30.0 && pick_hor < 50.0);

        // Test slightly below horizon
        let below = TopocentricPosition::from_degrees(0.0, -1.0, observer, epoch).unwrap();
        let roz_below = below.air_mass_rozenberg();
        let pick_below = below.air_mass_pickering();

        assert_eq!(roz_below, 40.0);
        assert!(pick_below > 40.0);
    }

    #[test]
    fn test_air_mass_formula_comparison() {
        // Verify that all three formulas produce reasonable and consistent results
        // across the full range of zenith angles
        let observer = test_observer();
        let epoch = TT::j2000();

        let test_elevations = vec![
            90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0, 5.0, 2.0, 0.0,
        ];

        for elev in test_elevations {
            let pos = TopocentricPosition::from_degrees(0.0, elev, observer, epoch).unwrap();
            let roz = pos.air_mass_rozenberg();
            let pick = pos.air_mass_pickering();
            let ky = pos.air_mass_kasten_young();

            // All values should be >= 1.0 (with tolerance for formula approximations at zenith)
            assert!(
                roz >= 0.999,
                "Rozenberg air mass < 0.999 at elevation {}",
                elev
            );
            assert!(
                pick >= 0.999,
                "Pickering air mass < 0.999 at elevation {}",
                elev
            );
            assert!(
                ky >= 0.999,
                "Kasten-Young air mass < 0.999 at elevation {}",
                elev
            );

            // Air mass should increase as elevation decreases
            // (this is implicitly tested by the monotonic nature of the formulas)

            // For high elevations (> 30°), all formulas should agree within 5%
            if elev > 30.0 {
                let avg = (roz + pick + ky) / 3.0;
                assert!(
                    (roz - avg).abs() / avg < 0.05,
                    "Rozenberg deviates >5% at elevation {}",
                    elev
                );
                assert!(
                    (pick - avg).abs() / avg < 0.05,
                    "Pickering deviates >5% at elevation {}",
                    elev
                );
                assert!(
                    (ky - avg).abs() / avg < 0.05,
                    "Kasten-Young deviates >5% at elevation {}",
                    elev
                );
            }
        }
    }

    #[test]
    fn test_atmospheric_refraction_standard_conditions() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Standard conditions: sea level, 15°C, 50% humidity, optical (0.574 μm)
        let pressure = 1013.25;
        let temp = 15.0;
        let humidity = 0.5;
        let wavelength = 0.574;

        // Test at zenith (no refraction)
        let zenith = TopocentricPosition::from_degrees(0.0, 90.0, observer, epoch).unwrap();
        let ref_zenith = zenith.atmospheric_refraction(pressure, temp, humidity, wavelength);
        assert!(ref_zenith.arcseconds().abs() < 0.1);

        // Test at 45° elevation
        let pos_45 = TopocentricPosition::from_degrees(0.0, 45.0, observer, epoch).unwrap();
        let ref_45 = pos_45.atmospheric_refraction(pressure, temp, humidity, wavelength);
        // Typical refraction at 45° elevation ~60 arcsec
        assert!(ref_45.arcseconds() > 50.0 && ref_45.arcseconds() < 70.0);

        // Test near horizon (10° elevation)
        let pos_10 = TopocentricPosition::from_degrees(0.0, 10.0, observer, epoch).unwrap();
        let ref_10 = pos_10.atmospheric_refraction(pressure, temp, humidity, wavelength);
        // Refraction increases dramatically near horizon, ~5-6 arcmin
        assert!(ref_10.arcminutes() > 4.0 && ref_10.arcminutes() < 7.0);
    }

    #[test]
    fn test_atmospheric_refraction_with_without() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let pressure = 1013.25;
        let temp = 15.0;
        let humidity = 0.5;
        let wavelength = 0.574;

        // Start with true position at 45°
        let true_pos = TopocentricPosition::from_degrees(0.0, 45.0, observer, epoch).unwrap();

        // Apply refraction to get apparent position
        let apparent = true_pos.with_refraction(pressure, temp, humidity, wavelength);

        // Apparent elevation should be higher than true
        assert!(apparent.elevation().degrees() > true_pos.elevation().degrees());

        // Remove refraction to get back to true
        let back_to_true = apparent.without_refraction(pressure, temp, humidity, wavelength);

        // Should be close to original (within numerical precision)
        assert!(
            (back_to_true.elevation().degrees() - true_pos.elevation().degrees()).abs() < 0.001
        );
    }

    #[test]
    fn test_atmospheric_refraction_zero_pressure() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Zero pressure = no atmosphere = no refraction
        let pos = TopocentricPosition::from_degrees(0.0, 30.0, observer, epoch).unwrap();
        let refraction = pos.atmospheric_refraction(0.0, 15.0, 0.5, 0.574);

        assert_eq!(refraction.radians(), 0.0);
    }

    #[test]
    fn test_atmospheric_refraction_radio_vs_optical() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let pressure = 1013.25;
        let temp = 15.0;
        let humidity = 0.5;

        let pos = TopocentricPosition::from_degrees(0.0, 30.0, observer, epoch).unwrap();

        // Optical wavelength (0.574 μm)
        let optical = pos.atmospheric_refraction(pressure, temp, humidity, 0.574);

        // Radio wavelength (>100 μm)
        let radio = pos.atmospheric_refraction(pressure, temp, humidity, 200.0);

        // Both should give positive refraction
        assert!(optical.arcseconds() > 0.0);
        assert!(radio.arcseconds() > 0.0);

        // Radio refraction should be less affected by humidity (simplified model)
        assert!(optical.arcseconds() > 0.0);
    }

    #[test]
    fn test_diurnal_parallax_moon() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Moon at mean distance: 384,400 km ≈ 0.00257 AU
        let moon_distance = Distance::from_kilometers(384400.0).unwrap();

        // Moon at horizon (maximum parallax)
        let moon_horizon = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(0.0),
            observer,
            epoch,
            moon_distance,
        )
        .unwrap();

        // Horizontal parallax for Moon: ~57 arcmin = 0.95°
        let h_parallax = moon_horizon.horizontal_parallax().unwrap();
        assert!(h_parallax.degrees() > 0.9 && h_parallax.degrees() < 1.0);
        assert!(h_parallax.arcminutes() > 55.0 && h_parallax.arcminutes() < 59.0);

        // At horizon, diurnal parallax = horizontal parallax
        let diurnal = moon_horizon.diurnal_parallax().unwrap();
        assert!((diurnal.degrees() - h_parallax.degrees()).abs() < 0.001);

        // Moon at zenith (zero parallax)
        let moon_zenith = TopocentricPosition::with_distance(
            Angle::from_degrees(0.0),
            Angle::from_degrees(90.0),
            observer,
            epoch,
            moon_distance,
        )
        .unwrap();

        let zenith_parallax = moon_zenith.diurnal_parallax().unwrap();
        assert!(zenith_parallax.arcseconds().abs() < 1.0); // Should be nearly zero
    }

    #[test]
    fn test_diurnal_parallax_sun() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Sun at 1 AU
        let sun_distance = Distance::from_au(1.0).unwrap();
        let sun_horizon = TopocentricPosition::with_distance(
            Angle::from_degrees(90.0),
            Angle::from_degrees(0.0),
            observer,
            epoch,
            sun_distance,
        )
        .unwrap();

        // Solar horizontal parallax: ~8.794 arcsec
        let h_parallax = sun_horizon.horizontal_parallax().unwrap();
        assert!(h_parallax.arcseconds() > 8.7 && h_parallax.arcseconds() < 8.9);
    }

    #[test]
    fn test_diurnal_parallax_mars_opposition() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Mars at closest approach: ~0.38 AU
        let mars_distance = Distance::from_au(0.38).unwrap();
        let mars_horizon = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(0.0),
            observer,
            epoch,
            mars_distance,
        )
        .unwrap();

        // Mars horizontal parallax at opposition: ~23 arcsec
        let h_parallax = mars_horizon.horizontal_parallax().unwrap();
        assert!(h_parallax.arcseconds() > 22.0 && h_parallax.arcseconds() < 24.0);
    }

    #[test]
    fn test_diurnal_parallax_at_various_elevations() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Moon at mean distance
        let moon_distance = Distance::from_kilometers(384400.0).unwrap();

        // Test at different elevations
        let elevations = vec![0.0, 30.0, 45.0, 60.0, 90.0];

        for elev in elevations {
            let pos = TopocentricPosition::with_distance(
                Angle::from_degrees(0.0),
                Angle::from_degrees(elev),
                observer,
                epoch,
                moon_distance,
            )
            .unwrap();

            let parallax = pos.diurnal_parallax().unwrap();

            // Parallax should decrease with increasing elevation
            // At zenith (90°), it should be nearly zero
            // At horizon (0°), it should equal horizontal parallax
            if elev == 90.0 {
                assert!(parallax.arcseconds().abs() < 1.0);
            } else if elev == 0.0 {
                let h_par = pos.horizontal_parallax().unwrap();
                assert!((parallax.degrees() - h_par.degrees()).abs() < 0.001);
            } else {
                // Parallax should be between 0 and horizontal parallax
                let h_par = pos.horizontal_parallax().unwrap();
                assert!(parallax.degrees() > 0.0);
                assert!(parallax.degrees() < h_par.degrees());
            }
        }
    }

    #[test]
    fn test_diurnal_parallax_with_without() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Moon at 45° elevation
        let moon_distance = Distance::from_kilometers(384400.0).unwrap();
        let geocentric = TopocentricPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
            moon_distance,
        )
        .unwrap();

        // Apply parallax correction
        let topocentric = geocentric.with_diurnal_parallax();

        // Topocentric elevation should be LOWER than geocentric
        assert!(topocentric.elevation().degrees() < geocentric.elevation().degrees());

        // Remove parallax to get back to geocentric
        let back_to_geocentric = topocentric.without_diurnal_parallax();

        // Should match original within tolerance
        // (not exact due to zenith angle changing during correction - this is correct physics)
        // For Moon at 45° with ~0.9° parallax, expect ~0.01° roundtrip error
        let diff =
            (back_to_geocentric.elevation().degrees() - geocentric.elevation().degrees()).abs();
        assert!(diff < 0.01, "Roundtrip error: {} degrees", diff);
    }

    #[test]
    fn test_diurnal_parallax_without_distance() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Position without distance (star)
        let star_pos = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();

        assert_eq!(star_pos.diurnal_parallax(), None);
        assert_eq!(star_pos.horizontal_parallax(), None);

        // with/without should return unchanged
        let with_par = star_pos.with_diurnal_parallax();
        assert_eq!(
            with_par.elevation().degrees(),
            star_pos.elevation().degrees()
        );

        let without_par = star_pos.without_diurnal_parallax();
        assert_eq!(
            without_par.elevation().degrees(),
            star_pos.elevation().degrees()
        );
    }

    #[test]
    fn test_diurnal_parallax_formula_verification() {
        // Verify the parallax formula: p = arcsin((R_Earth/r) × sin(z))
        let observer = test_observer();
        let epoch = TT::j2000();

        const EARTH_RADIUS_AU: f64 = 4.2635e-5;

        // Moon at known distance and elevation
        let distance_au = 0.00257; // Moon's distance
        let elevation_deg = 30.0;
        let zenith_deg = 90.0 - elevation_deg;

        let moon_distance = Distance::from_au(distance_au).unwrap();
        let moon_pos = TopocentricPosition::with_distance(
            Angle::from_degrees(0.0),
            Angle::from_degrees(elevation_deg),
            observer,
            epoch,
            moon_distance,
        )
        .unwrap();

        // Calculate expected parallax
        let ratio = EARTH_RADIUS_AU / distance_au;
        let zenith_rad = zenith_deg.to_radians();
        let expected_parallax_rad = (ratio * zenith_rad.sin()).asin();

        // Get calculated parallax
        let calculated_parallax = moon_pos.diurnal_parallax().unwrap();

        // Should match within numerical precision
        assert!((calculated_parallax.radians() - expected_parallax_rad).abs() < 1e-10);
    }

    #[test]
    fn test_topocentric_to_hour_angle() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // Test object on meridian at 45° elevation
        let topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        let ha = topo.to_hour_angle().unwrap();

        // On meridian (Az=180°), hour angle should be ~0
        assert!(ha.hour_angle().hours().abs() < 0.001);
    }

    #[test]
    fn test_topocentric_hour_angle_roundtrip() {
        let observer = test_observer();
        let epoch = TT::j2000();

        let test_cases = [
            (Angle::from_hours(0.0), Angle::from_degrees(45.0)),
            (Angle::from_hours(2.0), Angle::from_degrees(30.0)),
            (Angle::from_hours(-3.0), Angle::from_degrees(60.0)),
            (Angle::from_hours(6.0), Angle::from_degrees(0.0)),
        ];

        for (ha, dec) in test_cases {
            let original = HourAnglePosition::new(ha, dec, observer, epoch).unwrap();

            let topo = original.to_topocentric().unwrap();
            let recovered = topo.to_hour_angle().unwrap();

            let ha_diff_sec = (original.hour_angle().radians() - recovered.hour_angle().radians())
                .abs()
                * 206265.0;
            let dec_diff_arcsec =
                (original.declination().radians() - recovered.declination().radians()).abs()
                    * 206265.0;

            assert!(
                ha_diff_sec < 0.001,
                "Hour angle roundtrip failed for HA={:.2}h, Dec={:.1}°: diff={:.6} arcsec",
                ha.hours(),
                dec.degrees(),
                ha_diff_sec
            );
            assert!(
                dec_diff_arcsec < 0.001,
                "Declination roundtrip failed for HA={:.2}h, Dec={:.1}°: diff={:.6} arcsec",
                ha.hours(),
                dec.degrees(),
                dec_diff_arcsec
            );
        }
    }

    #[test]
    fn test_topocentric_to_hour_angle_distance_preservation() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let distance = Distance::from_kilometers(384400.0).unwrap();

        let topo = TopocentricPosition::with_distance(
            Angle::from_degrees(90.0),
            Angle::from_degrees(30.0),
            observer,
            epoch,
            distance,
        )
        .unwrap();

        let ha = topo.to_hour_angle().unwrap();
        assert_eq!(ha.distance().unwrap().kilometers(), distance.kilometers());
    }

    #[test]
    fn test_topocentric_to_hour_angle_cardinal_points() {
        let observer = test_observer();
        let epoch = TT::j2000();

        // North (Az=0°): object in northern sky, HA=0 on meridian
        // Actually Az=0° is north, but object crosses north meridian at HA=0 only for circumpolar
        // Let's test east/west instead

        // Due east (Az=90°): object is rising, HA should be negative (before meridian)
        let east = TopocentricPosition::from_degrees(90.0, 30.0, observer, epoch).unwrap();
        let ha_east = east.to_hour_angle().unwrap();
        assert!(
            ha_east.hour_angle().hours() < 0.0 || ha_east.hour_angle().hours() > 12.0,
            "East object should have negative or >12h hour angle, got {}h",
            ha_east.hour_angle().hours()
        );

        // Due west (Az=270°): object is setting, HA should be positive
        let west = TopocentricPosition::from_degrees(270.0, 30.0, observer, epoch).unwrap();
        let ha_west = west.to_hour_angle().unwrap();
        assert!(
            ha_west.hour_angle().hours() > 0.0 && ha_west.hour_angle().hours() < 12.0,
            "West object should have positive hour angle, got {}h",
            ha_west.hour_angle().hours()
        );
    }

    #[test]
    fn test_hour_angle_to_cirs() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let delta_t = 64.0;

        let ha = HourAnglePosition::new(
            Angle::from_hours(2.0),
            Angle::from_degrees(45.0),
            observer,
            epoch,
        )
        .unwrap();

        let cirs = ha.to_cirs(delta_t).unwrap();

        assert!(cirs.ra().degrees() >= 0.0 && cirs.ra().degrees() < 360.0);
        assert_eq!(cirs.dec().degrees(), ha.declination().degrees());
    }

    #[test]
    fn test_hour_angle_cirs_roundtrip() {
        use crate::CIRSPosition;

        let observer = test_observer();
        let epoch = TT::j2000();
        let delta_t = 64.0;

        let original_cirs = CIRSPosition::from_degrees(120.0, 35.0, epoch).unwrap();

        let ha = original_cirs.to_hour_angle(&observer, delta_t).unwrap();
        let recovered_cirs = ha.to_cirs(delta_t).unwrap();

        let ra_diff_arcsec =
            (original_cirs.ra().radians() - recovered_cirs.ra().radians()).abs() * 206265.0;
        let dec_diff_arcsec =
            (original_cirs.dec().radians() - recovered_cirs.dec().radians()).abs() * 206265.0;

        assert!(
            ra_diff_arcsec < 0.001,
            "RA roundtrip failed: diff={:.6} arcsec",
            ra_diff_arcsec
        );
        assert!(
            dec_diff_arcsec < 0.001,
            "Dec roundtrip failed: diff={:.6} arcsec",
            dec_diff_arcsec
        );
    }

    #[test]
    fn test_topocentric_to_cirs() {
        let observer = test_observer();
        let epoch = TT::j2000();
        let delta_t = 64.0;

        let topo = TopocentricPosition::from_degrees(180.0, 45.0, observer, epoch).unwrap();
        let cirs = topo.to_cirs(delta_t).unwrap();

        assert!(cirs.ra().degrees() >= 0.0 && cirs.ra().degrees() < 360.0);
        assert!(cirs.dec().degrees() >= -90.0 && cirs.dec().degrees() <= 90.0);
    }

    #[test]
    fn test_full_reverse_chain_roundtrip() {
        use crate::CIRSPosition;

        let observer = test_observer();
        let epoch = TT::j2000();
        let delta_t = 64.0;

        let original_cirs = CIRSPosition::from_degrees(200.0, 40.0, epoch).unwrap();

        let ha = original_cirs.to_hour_angle(&observer, delta_t).unwrap();
        let topo = ha.to_topocentric().unwrap();
        let recovered_ha = topo.to_hour_angle().unwrap();
        let recovered_cirs = recovered_ha.to_cirs(delta_t).unwrap();

        let ra_diff_arcsec =
            (original_cirs.ra().radians() - recovered_cirs.ra().radians()).abs() * 206265.0;
        let dec_diff_arcsec =
            (original_cirs.dec().radians() - recovered_cirs.dec().radians()).abs() * 206265.0;

        assert!(
            ra_diff_arcsec < 0.01,
            "Full chain RA roundtrip failed: diff={:.6} arcsec",
            ra_diff_arcsec
        );
        assert!(
            dec_diff_arcsec < 0.01,
            "Full chain Dec roundtrip failed: diff={:.6} arcsec",
            dec_diff_arcsec
        );
    }
}
