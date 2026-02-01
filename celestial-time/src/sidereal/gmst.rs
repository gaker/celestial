use super::angle::SiderealAngle;
use crate::scales::{TT, UT1};
use crate::JulianDate;
use crate::TimeResult;
use celestial_core::angle::wrap_0_2pi;
use celestial_core::constants::{J2000_JD, TWOPI};
use celestial_core::math::fmod;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GMST(SiderealAngle);

impl GMST {
    pub fn from_ut1_and_tt(ut1: &UT1, tt: &TT) -> TimeResult<Self> {
        let gmst_rad = calculate_gmst_iau2006(ut1, tt)?;
        let angle = SiderealAngle::from_radians_exact(gmst_rad);
        Ok(Self(angle))
    }

    pub fn from_hours(hours: f64) -> Self {
        Self(SiderealAngle::from_hours(hours))
    }

    pub fn from_degrees(degrees: f64) -> Self {
        Self(SiderealAngle::from_degrees(degrees))
    }

    pub fn from_radians(radians: f64) -> Self {
        Self(SiderealAngle::from_radians(radians))
    }

    pub fn j2000() -> TimeResult<Self> {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        Self::from_ut1_and_tt(&ut1, &tt)
    }

    pub fn angle(&self) -> SiderealAngle {
        self.0
    }

    pub fn hours(&self) -> f64 {
        self.0.hours()
    }

    pub fn degrees(&self) -> f64 {
        self.0.degrees()
    }

    pub fn radians(&self) -> f64 {
        self.0.radians()
    }

    pub fn hour_angle_to_target(&self, target_ra_hours: f64) -> f64 {
        self.0.hour_angle_to_target(target_ra_hours)
    }

    pub fn to_lmst(&self, location: &celestial_core::Location) -> crate::sidereal::LMST {
        use super::angle::SiderealAngle;
        use crate::sidereal::LMST;

        let gmst_rad = self.radians();
        let lmst_rad = gmst_rad + location.longitude;

        let lmst_normalized = wrap_0_2pi(lmst_rad);
        let angle = SiderealAngle::from_radians_exact(lmst_normalized);

        LMST::from_radians(angle.radians(), location)
    }
}

impl std::fmt::Display for GMST {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GMST {}", self.0)
    }
}

fn calculate_gmst_iau2006(ut1: &UT1, tt: &TT) -> TimeResult<f64> {
    let ut1_jd = ut1.to_julian_date();
    let tt_jd = tt.to_julian_date();

    let JulianDate {
        jd1: ut1_jd1,
        jd2: ut1_jd2,
    } = ut1_jd;
    let JulianDate {
        jd1: tt_jd1,
        jd2: tt_jd2,
    } = tt_jd;

    let t = ((tt_jd1 - J2000_JD) + tt_jd2) / celestial_core::constants::DAYS_PER_JULIAN_CENTURY;

    let era = calculate_era00(ut1_jd1, ut1_jd2)?;

    // IAU 2006 polynomial correction (Horner's method for precision)
    let polynomial_arcsec = 0.014506
        + t * (4612.156534
            + t * (1.3915817 + t * (-0.00000044 + t * (-0.000029956 + t * (-0.0000000368)))));

    let gmst = era + polynomial_arcsec * celestial_core::constants::ARCSEC_TO_RAD;

    Ok(wrap_0_2pi(gmst))
}

fn calculate_era00(ut1_jd1: f64, ut1_jd2: f64) -> TimeResult<f64> {
    let (d1, d2) = if ut1_jd1 < ut1_jd2 {
        (ut1_jd1, ut1_jd2)
    } else {
        (ut1_jd2, ut1_jd1)
    };

    let t = d1 + (d2 - J2000_JD);

    if t.is_infinite() || t.is_nan() || t.abs() > 1e12 {
        return Err(crate::TimeError::CalculationError(format!(
            "Time value out of valid range: {} days from J2000",
            t
        )));
    }

    let f = fmod(d1, 1.0) + fmod(d2, 1.0);

    let rotation_term = 0.00273781191135448 * t;

    if rotation_term.is_infinite() || rotation_term.is_nan() {
        return Err(crate::TimeError::CalculationError(
            "Earth rotation calculation overflow".to_string(),
        ));
    }

    let theta = TWOPI * (f + 0.7790572732640 + rotation_term);

    Ok(wrap_0_2pi(theta))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gmst_j2000() {
        let gmst = GMST::j2000().unwrap();

        let hours = gmst.hours();
        assert!(
            (0.0..24.0).contains(&hours),
            "GMST should be in [0, 24) hours: {}",
            hours
        );
        assert!(
            hours > 18.0 && hours < 19.0,
            "GMST at J2000.0 should be ~18.7 hours: {}",
            hours
        );
    }

    #[test]
    fn test_gmst_from_ut1_and_tt() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let gmst = GMST::from_ut1_and_tt(&ut1, &tt).unwrap();

        let gmst_j2000 = GMST::j2000().unwrap();
        assert!((gmst.hours() - gmst_j2000.hours()).abs() < 1e-10);
    }

    #[test]
    fn test_hour_angle_calculation() {
        let gmst = GMST::from_hours(12.0);
        let target_ra = 6.0;
        let hour_angle = gmst.hour_angle_to_target(target_ra);
        assert_eq!(hour_angle, 6.0);
    }

    #[test]
    fn test_gmst_constructors_and_accessors() {
        let gmst_deg = GMST::from_degrees(180.0);
        assert_eq!(gmst_deg.degrees(), 180.0);
        assert_eq!(gmst_deg.hours(), 12.0);

        let gmst_rad = GMST::from_radians(celestial_core::constants::PI);
        assert!((gmst_rad.radians() - celestial_core::constants::PI).abs() < 1e-15);
        assert_eq!(gmst_rad.hours(), 12.0);

        let angle = gmst_deg.angle();
        assert_eq!(angle.degrees(), 180.0);

        let degrees = gmst_deg.degrees();
        assert_eq!(degrees, 180.0);

        let display_str = format!("{}", gmst_deg);
        assert!(display_str.contains("GMST"));
        assert!(display_str.contains("12.000000h"));
    }

    #[test]
    fn test_overflow_protection() {
        use crate::scales::{TT, UT1};
        use crate::JulianDate;

        let extreme_jd = JulianDate::new(J2000_JD + 1e13, 0.0);
        let extreme_ut1 = UT1::from_julian_date(extreme_jd);
        let extreme_tt = TT::from_julian_date(extreme_jd);

        let result = GMST::from_ut1_and_tt(&extreme_ut1, &extreme_tt);
        assert!(result.is_err(), "Expected overflow protection to trigger");

        if let Err(err) = result {
            let error_message = format!("{}", err);
            assert!(
                error_message.contains("Time value out of valid range"),
                "Expected overflow error message, got: {}",
                error_message
            );
        }

        let infinite_jd = JulianDate::new(f64::INFINITY, 0.0);
        let infinite_ut1 = UT1::from_julian_date(infinite_jd);
        let infinite_tt = TT::from_julian_date(infinite_jd);

        let result = GMST::from_ut1_and_tt(&infinite_ut1, &infinite_tt);
        assert!(
            result.is_err(),
            "Expected infinite value protection to trigger"
        );

        let nan_jd = JulianDate::new(f64::NAN, 0.0);
        let nan_ut1 = UT1::from_julian_date(nan_jd);
        let nan_tt = TT::from_julian_date(nan_jd);

        let result = GMST::from_ut1_and_tt(&nan_ut1, &nan_tt);
        assert!(result.is_err(), "Expected NaN value protection to trigger");
    }
}
