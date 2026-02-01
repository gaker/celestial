use super::angle::SiderealAngle;
use crate::scales::{TT, UT1};
use crate::sidereal::LAST;
use crate::transforms::earth_rotation_angle;
use crate::transforms::nutation::NutationCalculator;
use crate::TimeResult;
use astro_core::angle::wrap_0_2pi;
use astro_core::cio::CioSolution;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GAST(SiderealAngle);

impl GAST {
    pub fn from_ut1_and_tt(ut1: &UT1, tt: &TT) -> TimeResult<Self> {
        let gast_rad = calculate_gast_iau2006a(ut1, tt)?;
        let angle = SiderealAngle::from_radians_exact(gast_rad);
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

    pub fn to_last(&self, location: &astro_core::Location) -> crate::sidereal::LAST {
        let gast_rad = self.radians();
        let last_rad = gast_rad + location.longitude;

        let last_normalized = wrap_0_2pi(last_rad);
        let angle = SiderealAngle::from_radians_exact(last_normalized);

        LAST::from_radians(angle.radians(), location)
    }
}

impl std::fmt::Display for GAST {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GAST {}", self.0)
    }
}

fn calculate_gast_iau2006a(ut1: &UT1, tt: &TT) -> TimeResult<f64> {
    let era = earth_rotation_angle(&ut1.to_julian_date())?;

    let tt_centuries = tt_to_centuries(tt)?;
    let npb_matrix = calculate_npb_matrix_iau2006a(tt)?;

    let cio_solution = CioSolution::calculate(&npb_matrix, tt_centuries).map_err(|e| {
        crate::TimeError::CalculationError(format!("CIO calculation failed: {}", e))
    })?;

    let gast = era - cio_solution.equation_of_origins;

    Ok(wrap_0_2pi(gast))
}

fn calculate_npb_matrix_iau2006a(tt: &TT) -> TimeResult<astro_core::RotationMatrix3> {
    let tt_jd = tt.to_julian_date();
    let t = astro_core::utils::jd_to_centuries(tt_jd.jd1(), tt_jd.jd2());

    let nutation_result = tt.nutation_iau2006a()?;
    let dpsi = nutation_result.nutation_longitude();
    let deps = nutation_result.nutation_obliquity();

    let precession_calc = astro_core::precession::PrecessionIAU2006::new();
    let npb_matrix = precession_calc.npb_matrix_iau2006a(t, dpsi, deps);

    Ok(npb_matrix)
}

fn tt_to_centuries(tt: &TT) -> TimeResult<f64> {
    crate::transforms::nutation::tt_to_centuries(tt)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gast_j2000() {
        let gast = GAST::j2000().unwrap();

        let hours = gast.hours();
        assert!(
            (0.0..24.0).contains(&hours),
            "GAST should be in [0, 24) hours: {}",
            hours
        );
        assert!(
            hours > 18.0 && hours < 19.0,
            "GAST at J2000.0 should be ~18.7 hours: {}",
            hours
        );
    }

    #[test]
    fn test_gast_from_ut1_and_tt() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let gast = GAST::from_ut1_and_tt(&ut1, &tt).unwrap();

        let gast_j2000 = GAST::j2000().unwrap();
        assert!((gast.hours() - gast_j2000.hours()).abs() < 1e-10);
    }

    #[test]
    fn test_hour_angle_calculation() {
        let gast = GAST::from_hours(12.0);
        let target_ra = 6.0;
        let hour_angle = gast.hour_angle_to_target(target_ra);
        assert_eq!(hour_angle, 6.0);
    }

    #[test]
    fn test_gast_constructors_and_accessors() {
        let gast_deg = GAST::from_degrees(270.0);
        assert_eq!(gast_deg.degrees(), 270.0);
        assert_eq!(gast_deg.hours(), 18.0);

        let gast_rad = GAST::from_radians(astro_core::constants::HALF_PI);
        assert!((gast_rad.radians() - astro_core::constants::HALF_PI).abs() < 1e-15);
        assert_eq!(gast_rad.hours(), 6.0);

        let angle = gast_deg.angle();
        assert_eq!(angle.degrees(), 270.0);

        let degrees = gast_deg.degrees();
        assert_eq!(degrees, 270.0);

        let display_str = format!("{}", gast_deg);
        assert!(display_str.contains("GAST"));
        assert!(display_str.contains("18.000000h"));
    }
}
