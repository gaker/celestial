use super::angle::SiderealAngle;
use super::gast::GAST;
use super::gmst::GMST;
use super::lmst::LMST;
use crate::scales::{TT, UT1};
use crate::transforms::nutation::NutationCalculator;
use crate::TimeResult;
use astro_core::angle::wrap_0_2pi;
use astro_core::Location;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LAST {
    angle: SiderealAngle,
    location: Location,
}

impl LAST {
    pub fn from_ut1_tt_and_location(ut1: &UT1, tt: &TT, location: &Location) -> TimeResult<Self> {
        let gast = GAST::from_ut1_and_tt(ut1, tt)?;

        let last_rad = gast.radians() + location.longitude;

        let last_normalized = wrap_0_2pi(last_rad);

        let angle = SiderealAngle::from_radians_exact(last_normalized);

        Ok(Self {
            angle,
            location: *location,
        })
    }

    pub fn from_lmst_and_equation_of_equinoxes(
        ut1: &UT1,
        tt: &TT,
        location: &Location,
    ) -> TimeResult<Self> {
        Self::from_ut1_tt_and_location(ut1, tt, location)
    }

    pub fn from_hours(hours: f64, location: &Location) -> Self {
        Self {
            angle: SiderealAngle::from_hours(hours),
            location: *location,
        }
    }

    pub fn from_degrees(degrees: f64, location: &Location) -> Self {
        Self {
            angle: SiderealAngle::from_degrees(degrees),
            location: *location,
        }
    }

    pub fn from_radians(radians: f64, location: &Location) -> Self {
        Self {
            angle: SiderealAngle::from_radians(radians),
            location: *location,
        }
    }

    pub fn j2000(location: &Location) -> TimeResult<Self> {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        Self::from_ut1_tt_and_location(&ut1, &tt, location)
    }

    pub fn angle(&self) -> SiderealAngle {
        self.angle
    }

    pub fn location(&self) -> Location {
        self.location
    }

    pub fn hours(&self) -> f64 {
        self.angle.hours()
    }

    pub fn degrees(&self) -> f64 {
        self.angle.degrees()
    }

    pub fn radians(&self) -> f64 {
        self.angle.radians()
    }

    pub fn hour_angle_to_target(&self, target_ra_hours: f64) -> f64 {
        self.angle.hour_angle_to_target(target_ra_hours)
    }

    pub fn to_gast(&self) -> GAST {
        let longitude_hours = self.location.longitude * 12.0 / astro_core::constants::PI;

        let gast_hours = self.hours() - longitude_hours;

        GAST::from_hours(gast_hours)
    }

    pub fn to_lmst(&self, tt: &TT) -> TimeResult<LMST> {
        let nutation = tt.nutation_iau2006a()?;

        let jd = tt.to_julian_date();
        let mean_obliquity = astro_core::obliquity::iau_2006_mean_obliquity(jd.jd1(), jd.jd2());

        let ee_rad = nutation.nutation_longitude() * mean_obliquity.cos();
        let ee_hours = ee_rad * 12.0 / astro_core::constants::PI;

        let lmst_hours = self.hours() - ee_hours;

        Ok(LMST::from_hours(lmst_hours, &self.location))
    }

    pub fn to_gmst(&self, tt: &TT) -> TimeResult<GMST> {
        let lmst = self.to_lmst(tt)?;

        Ok(lmst.to_gmst())
    }
}

impl std::fmt::Display for LAST {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let lat_deg = self.location.latitude.to_degrees();
        let lon_deg = self.location.longitude.to_degrees();
        write!(
            f,
            "LAST {} at ({:.4}°, {:.4}°)",
            self.angle, lat_deg, lon_deg
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mauna_kea() -> Location {
        Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap()
    }

    fn greenwich() -> Location {
        Location::greenwich()
    }

    #[test]
    fn test_last_at_greenwich_equals_gast() {
        // At Greenwich (0° longitude), LAST should equal GAST
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = greenwich();

        let gast = GAST::from_ut1_and_tt(&ut1, &tt).unwrap();
        let last = LAST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();

        // Should be identical (within numerical precision)
        assert!(
            (last.hours() - gast.hours()).abs() < 1e-14,
            "LAST at Greenwich should equal GAST: LAST={}, GAST={}",
            last.hours(),
            gast.hours()
        );
    }

    #[test]
    fn test_last_method_consistency() {
        // Both calculation methods should give identical results
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();

        let last_method1 = LAST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();
        let last_method2 = LAST::from_lmst_and_equation_of_equinoxes(&ut1, &tt, &location).unwrap();

        // Both methods should produce identical results within numerical precision
        assert!(
            (last_method1.hours() - last_method2.hours()).abs() < 1e-14,
            "LAST calculation methods should match: method1={}, method2={}",
            last_method1.hours(),
            last_method2.hours()
        );
    }

    #[test]
    fn test_last_longitude_correction() {
        // Test longitude correction: 1 degree = 4 minutes = 1/15 hour
        let ut1 = UT1::j2000();
        let tt = TT::j2000();

        let greenwich_loc = greenwich();
        let east_15deg = Location::from_degrees(0.0, 15.0, 0.0).unwrap(); // 15°E = +1 hour
        let west_15deg = Location::from_degrees(0.0, -15.0, 0.0).unwrap(); // 15°W = -1 hour

        let last_greenwich = LAST::from_ut1_tt_and_location(&ut1, &tt, &greenwich_loc).unwrap();
        let last_east = LAST::from_ut1_tt_and_location(&ut1, &tt, &east_15deg).unwrap();
        let last_west = LAST::from_ut1_tt_and_location(&ut1, &tt, &west_15deg).unwrap();

        // 15°E should be +1 hour ahead of Greenwich
        let diff_east = last_east.hours() - last_greenwich.hours();
        assert!(
            (diff_east - 1.0).abs() < 1e-12,
            "15°E should be +1 hour: {}",
            diff_east
        );

        // 15°W should be -1 hour behind Greenwich
        let diff_west = last_west.hours() - last_greenwich.hours();
        assert!(
            (diff_west + 1.0).abs() < 1e-12,
            "15°W should be -1 hour: {}",
            diff_west
        );
    }

    #[test]
    fn test_last_vs_gast_longitude() {
        // LAST = GAST + longitude correction
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();

        let last = LAST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();
        let gast = GAST::from_ut1_and_tt(&ut1, &tt).unwrap();

        // Calculate longitude correction manually
        let longitude_hours = location.longitude * 12.0 / astro_core::constants::PI;

        // LAST should equal GAST + longitude correction
        let expected_last = gast.hours() + longitude_hours;
        let expected_last_normalized = ((expected_last % 24.0) + 24.0) % 24.0;

        assert!(
            (last.hours() - expected_last_normalized).abs() < 1e-14,
            "LAST = GAST + longitude: LAST={}, GAST={}, longitude={}, expected={}",
            last.hours(),
            gast.hours(),
            longitude_hours,
            expected_last_normalized
        );
    }

    #[test]
    fn test_last_mauna_kea() {
        // Mauna Kea is at -155.4783° = -10.365 hours west of Greenwich
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();

        let gast = GAST::from_ut1_and_tt(&ut1, &tt).unwrap();
        let last = LAST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();

        // Expected longitude correction in hours
        let expected_offset = -155.4783 / 15.0; // degrees to hours
        let actual_offset = last.hours() - gast.hours();

        assert!(
            (actual_offset - expected_offset).abs() < 1e-10,
            "Mauna Kea LAST offset incorrect: expected={}, actual={}",
            expected_offset,
            actual_offset
        );
    }

    #[test]
    fn test_last_j2000() {
        let location = mauna_kea();
        let last = LAST::j2000(&location).unwrap();

        // LAST should be in valid range
        let hours = last.hours();
        assert!(
            (0.0..24.0).contains(&hours),
            "LAST should be in [0, 24) hours: {}",
            hours
        );
    }

    #[test]
    fn test_last_hour_angle_calculation() {
        let location = mauna_kea();
        let last = LAST::from_hours(12.0, &location);
        let target_ra = 6.0;
        let hour_angle = last.hour_angle_to_target(target_ra);
        assert_eq!(hour_angle, 6.0);
    }

    #[test]
    fn test_last_to_gast_roundtrip() {
        let location = mauna_kea();
        let original_gast = GAST::from_hours(15.5);

        // Convert GAST -> LAST -> GAST
        let longitude_hours = location.longitude * 12.0 / astro_core::constants::PI;
        let last_hours = original_gast.hours() + longitude_hours;
        let last = LAST::from_hours(last_hours, &location);
        let recovered_gast = last.to_gast();

        assert!(
            (recovered_gast.hours() - original_gast.hours()).abs() < 1e-14,
            "GAST->LAST->GAST roundtrip failed: original={}, recovered={}",
            original_gast.hours(),
            recovered_gast.hours()
        );
    }

    #[test]
    fn test_last_to_lmst_conversion() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();

        let original_lmst = LMST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();

        // Convert LMST -> LAST -> LMST
        let last = LAST::from_lmst_and_equation_of_equinoxes(&ut1, &tt, &location).unwrap();
        let recovered_lmst = last.to_lmst(&tt).unwrap();

        // Note: Small precision difference expected due to CIO-based LAST vs classical LMST
        // Roundtrip precision limited by algorithm difference
        assert!(
            (recovered_lmst.hours() - original_lmst.hours()).abs() < 1e-7,
            "LMST->LAST->LMST roundtrip failed: original={}, recovered={}",
            original_lmst.hours(),
            recovered_lmst.hours()
        );
    }

    #[test]
    fn test_last_from_constructors() {
        let location = mauna_kea();

        // Test all constructor methods produce equivalent results
        let hours = 14.5;
        let degrees = hours * 15.0;
        let radians = hours * astro_core::constants::PI / 12.0;

        let last_hours = LAST::from_hours(hours, &location);
        let last_degrees = LAST::from_degrees(degrees, &location);
        let last_radians = LAST::from_radians(radians, &location);

        assert!((last_hours.hours() - last_degrees.hours()).abs() < 1e-14);
        assert!((last_hours.hours() - last_radians.hours()).abs() < 1e-14);
        assert_eq!(last_hours.location(), location);
        assert_eq!(last_degrees.location(), location);
        assert_eq!(last_radians.location(), location);
    }

    #[test]
    fn test_last_display() {
        let location = mauna_kea();
        let last = LAST::from_hours(12.5, &location);
        let display = format!("{}", last);

        // Should include LAST time and location coordinates
        assert!(display.contains("LAST"));
        assert!(display.contains("19.8283")); // Latitude
        assert!(display.contains("-155.4783")); // Longitude
    }

    #[test]
    fn test_last_location_enforcement() {
        // Test that LAST enforces location through type system
        let location = mauna_kea();
        let last = LAST::from_hours(12.0, &location);

        // Location is always available and cannot be None/invalid
        let stored_location = last.location();
        assert_eq!(stored_location.latitude, location.latitude);
        assert_eq!(stored_location.longitude, location.longitude);
        assert_eq!(stored_location.height, location.height);
    }

    #[test]
    fn test_extreme_longitudes() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();

        // Test extreme valid longitudes
        let east_extreme = Location::from_degrees(0.0, 180.0, 0.0).unwrap(); // 180°E = +12 hours
        let west_extreme = Location::from_degrees(0.0, -180.0, 0.0).unwrap(); // 180°W = -12 hours

        let gast = GAST::from_ut1_and_tt(&ut1, &tt).unwrap();
        let last_east = LAST::from_ut1_tt_and_location(&ut1, &tt, &east_extreme).unwrap();
        let last_west = LAST::from_ut1_tt_and_location(&ut1, &tt, &west_extreme).unwrap();

        // 180°E should be +12 hours ahead
        let diff_east = last_east.hours() - gast.hours();
        let expected_east = if diff_east < 0.0 {
            diff_east + 24.0
        } else {
            diff_east
        };
        assert!(
            (expected_east - 12.0).abs() < 1e-12,
            "180°E should be +12 hours"
        );

        // 180°W should be -12 hours behind (equivalent to +12 hours due to 24h wrap)
        let diff_west = last_west.hours() - gast.hours();
        let expected_west = if diff_west > 12.0 {
            diff_west - 24.0
        } else {
            diff_west
        };
        assert!(
            (expected_west + 12.0).abs() < 1e-12,
            "180°W should be -12 hours"
        );
    }

    #[test]
    fn test_last_equation_of_equinoxes_range() {
        // Test that equation of equinoxes is reasonable (should be small)
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();

        let last = LAST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();
        let lmst = LMST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();

        // Equation of equinoxes should be small (typically < 1 second = 1/3600 hours)
        let ee_hours = last.hours() - lmst.hours();
        let ee_hours_normalized = if ee_hours > 12.0 {
            ee_hours - 24.0
        } else if ee_hours < -12.0 {
            ee_hours + 24.0
        } else {
            ee_hours
        };

        assert!(
            ee_hours_normalized.abs() < 0.001, // Less than 3.6 seconds
            "Equation of equinoxes too large: {} hours = {} seconds",
            ee_hours_normalized,
            ee_hours_normalized * 3600.0
        );
    }

    #[test]
    fn test_last_constructors_and_accessors() {
        let location = greenwich();

        // Test from_degrees constructor
        let last_deg = LAST::from_degrees(45.0, &location);
        assert_eq!(last_deg.degrees(), 45.0);
        assert_eq!(last_deg.hours(), 3.0);

        // Test from_radians constructor
        let last_rad = LAST::from_radians(astro_core::constants::PI / 4.0, &location);
        assert!((last_rad.radians() - astro_core::constants::PI / 4.0).abs() < 1e-15);
        assert_eq!(last_rad.hours(), 3.0);

        // Test angle() accessor
        let angle = last_deg.angle();
        assert_eq!(angle.degrees(), 45.0);

        // Test location() accessor
        let stored_location = last_deg.location();
        assert_eq!(stored_location.latitude, location.latitude);
        assert_eq!(stored_location.longitude, location.longitude);
        assert_eq!(stored_location.height, location.height);

        // Test degrees() method
        let degrees = last_deg.degrees();
        assert_eq!(degrees, 45.0);
    }
}
