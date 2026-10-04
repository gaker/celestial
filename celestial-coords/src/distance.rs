use crate::errors::{CoordError, CoordResult};
use celestial_core::angle::Angle;
use celestial_core::constants::{
    ARCSEC_PER_RAD, AU_KM, AU_M, DAYS_PER_JULIAN_YEAR, MILLIARCSEC_TO_RAD, SECONDS_PER_DAY_F64,
    SPEED_OF_LIGHT_M_PER_S,
};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// IAU 2015 B2 defines the parsec as 648000/π au. The light year is c times a Julian year.
const AU_PER_PARSEC: f64 = ARCSEC_PER_RAD;
const KM_PER_PARSEC: f64 = ARCSEC_PER_RAD * AU_KM;
const LIGHT_YEAR_M: f64 = SPEED_OF_LIGHT_M_PER_S * SECONDS_PER_DAY_F64 * DAYS_PER_JULIAN_YEAR;
const LIGHT_YEARS_PER_PARSEC: f64 = ARCSEC_PER_RAD * (AU_M / LIGHT_YEAR_M);

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Distance {
    parsecs: f64,
}

impl Distance {
    /// Creates a Distance from parsecs.
    ///
    /// # Valid Range
    /// Must be positive and finite (0 < parsecs < ∞)
    ///
    /// # Errors
    /// Returns `CoordError::InvalidDistance` if value is ≤0, infinite, or NaN.
    pub fn from_parsecs(parsecs: f64) -> CoordResult<Self> {
        if !parsecs.is_finite() || parsecs <= 0.0 {
            return Err(CoordError::invalid_distance(format!(
                "Distance must be positive and finite, got {}",
                parsecs
            )));
        }
        Ok(Self { parsecs })
    }

    /// Creates a Distance from light-years.
    ///
    /// # Valid Range
    /// Must be positive and finite (0 < ly < ∞), and not so small that it underflows to 0 pc
    pub fn from_light_years(ly: f64) -> CoordResult<Self> {
        Self::from_parsecs(ly / LIGHT_YEARS_PER_PARSEC)
    }

    /// Creates a Distance from astronomical units.
    ///
    /// # Valid Range
    /// Must be positive and finite (0 < au < ∞), and not so small that it underflows to 0 pc
    pub fn from_au(au: f64) -> CoordResult<Self> {
        Self::from_parsecs(au / AU_PER_PARSEC)
    }

    /// Creates a Distance from kilometers.
    ///
    /// # Valid Range
    /// Must be positive and finite (0 < km < ∞), and not so small that it underflows to 0 pc
    pub fn from_kilometers(km: f64) -> CoordResult<Self> {
        Self::from_parsecs(km / KM_PER_PARSEC)
    }

    /// Creates a Distance from parallax in arcseconds.
    ///
    /// # Valid Range
    /// Must be positive and finite (0 < parallax_arcsec < ∞), and not so small that the distance
    /// overflows to ∞ pc
    ///
    /// # Note
    /// Distance (parsecs) = 1 / parallax (arcsec)
    pub fn from_parallax_arcsec(parallax_arcsec: f64) -> CoordResult<Self> {
        if !parallax_arcsec.is_finite() || parallax_arcsec <= 0.0 {
            return Err(CoordError::invalid_distance(format!(
                "Parallax must be positive and finite, got {} arcsec",
                parallax_arcsec
            )));
        }
        Self::from_parsecs(1.0 / parallax_arcsec)
    }

    pub fn from_parallax_milliarcsec(parallax_mas: f64) -> CoordResult<Self> {
        Self::from_parallax_arcsec(parallax_mas / 1000.0)
    }

    pub fn from_parallax_angle(parallax: Angle) -> CoordResult<Self> {
        Self::from_parallax_arcsec(parallax.arcseconds())
    }

    pub fn parsecs(self) -> f64 {
        self.parsecs
    }

    pub fn light_years(self) -> f64 {
        self.parsecs * LIGHT_YEARS_PER_PARSEC
    }

    pub fn au(self) -> f64 {
        self.parsecs * AU_PER_PARSEC
    }

    pub fn kilometers(self) -> f64 {
        self.parsecs * KM_PER_PARSEC
    }

    pub fn parallax_arcsec(self) -> f64 {
        1.0 / self.parsecs
    }

    pub fn parallax_milliarcsec(self) -> f64 {
        self.parallax_arcsec() * 1000.0
    }

    pub fn parallax_angle(self) -> Angle {
        Angle::from_arcseconds(self.parallax_arcsec())
    }

    pub fn distance_modulus(self) -> f64 {
        5.0 * libm::log10(self.parsecs) - 5.0
    }

    pub fn from_distance_modulus(dm: f64) -> CoordResult<Self> {
        let parsecs = libm::pow(10.0, (dm + 5.0) / 5.0);
        Self::from_parsecs(parsecs)
    }

    pub fn is_galactic(self) -> bool {
        self.parsecs < 100_000.0
    }

    pub fn is_local_group(self) -> bool {
        self.parsecs < 2_000_000.0
    }

    pub fn parallax_uncertainty_mas(self, relative_error: f64) -> f64 {
        let parallax_mas = self.parallax_milliarcsec();
        parallax_mas * relative_error
    }

    pub fn proper_motion_distance_au(self, pm_mas_per_year: f64, dt_years: f64) -> f64 {
        let angular_distance_rad = pm_mas_per_year * MILLIARCSEC_TO_RAD * dt_years;
        self.au() * angular_distance_rad
    }
}

impl std::ops::Add for Distance {
    type Output = CoordResult<Self>;

    fn add(self, other: Self) -> Self::Output {
        Self::from_parsecs(self.parsecs + other.parsecs)
    }
}

impl std::ops::Sub for Distance {
    type Output = CoordResult<Self>;

    fn sub(self, other: Self) -> Self::Output {
        Self::from_parsecs(self.parsecs - other.parsecs)
    }
}

impl std::ops::Mul<f64> for Distance {
    type Output = CoordResult<Self>;

    fn mul(self, factor: f64) -> Self::Output {
        Self::from_parsecs(self.parsecs * factor)
    }
}

impl std::ops::Div<f64> for Distance {
    type Output = CoordResult<Self>;

    fn div(self, divisor: f64) -> Self::Output {
        Self::from_parsecs(self.parsecs / divisor)
    }
}

impl PartialOrd for Distance {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.parsecs.partial_cmp(&other.parsecs)
    }
}

impl std::fmt::Display for Distance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.parsecs < 1e-3 {
            write!(f, "{:.3} AU", self.au())
        } else if self.parsecs < 1000.0 {
            write!(f, "{:.3} pc", self.parsecs)
        } else if self.parsecs < 1e6 {
            write!(f, "{:.3} kpc", self.parsecs / 1000.0)
        } else {
            write!(f, "{:.3} Mpc", self.parsecs / 1e6)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::rounded;

    #[test]
    fn test_distance_creation() {
        let d1 = Distance::from_parsecs(10.0).unwrap();
        assert_eq!(d1.parsecs(), 10.0);

        let d2 = Distance::from_parallax_arcsec(0.1).unwrap();
        assert_eq!(d2.parsecs(), 10.0);

        assert!(Distance::from_parsecs(-1.0).is_err());
        assert!(Distance::from_parsecs(0.0).is_err());
        assert!(Distance::from_parallax_arcsec(0.0).is_err());
    }

    // Expected values are the defining ratios rounded once to f64.
    #[test]
    fn test_parsec_is_648000_over_pi_au() {
        assert_eq!(
            Distance::from_parsecs(1.0).unwrap().au(),
            206264.80624709636
        );
        assert_eq!(
            Distance::from_au(206264.80624709636).unwrap().parsecs(),
            1.0
        );
    }

    #[test]
    fn test_light_year_is_c_times_julian_year() {
        assert_eq!(
            Distance::from_parsecs(1.0).unwrap().light_years(),
            3.2615637771674337
        );
        assert_eq!(
            Distance::from_light_years(3.2615637771674337)
                .unwrap()
                .parsecs(),
            1.0
        );
        assert_eq!(
            Distance::from_light_years(1.0).unwrap().parsecs(),
            0.30660139378555057
        );
    }

    #[test]
    fn test_kilometres_round_trip() {
        assert_eq!(
            Distance::from_parsecs(1.0).unwrap().kilometers(),
            30856775814913.67
        );
        assert_eq!(
            Distance::from_kilometers(30856775814913.67)
                .unwrap()
                .parsecs(),
            1.0
        );
    }

    #[test]
    fn test_parallax_angle() {
        let angle = Angle::from_arcseconds(0.1);
        let d = Distance::from_parallax_angle(angle).unwrap();
        assert_eq!(d.parsecs(), 10.0);
    }

    #[test]
    fn test_parallax_uncertainty_mas() {
        let d = Distance::from_parsecs(100.0).unwrap();
        let unc = d.parallax_uncertainty_mas(0.01);
        assert_eq!(unc, 0.1);
    }

    #[test]
    fn test_partial_ord() {
        let d1 = Distance::from_parsecs(10.0).unwrap();
        let d2 = Distance::from_parsecs(20.0).unwrap();
        assert!(d1 < d2);
    }

    #[test]
    fn test_parallax_calculations() {
        let proxima = Distance::from_parallax_arcsec(0.7687).unwrap();
        assert_eq!(rounded(proxima.parsecs(), 4), 1.3009);

        let distance = Distance::from_parallax_milliarcsec(768.7).unwrap();
        assert_eq!(distance, proxima);
    }

    #[test]
    fn test_distance_modulus() {
        let distance = Distance::from_parsecs(10.0).unwrap();
        let dm = distance.distance_modulus();
        assert_eq!(dm, 0.0);

        let recovered = Distance::from_distance_modulus(dm).unwrap();
        assert_eq!(recovered.parsecs(), 10.0);
    }

    #[test]
    fn test_distance_scales() {
        let galactic = Distance::from_parsecs(1000.0).unwrap();
        assert!(galactic.is_galactic());
        assert!(galactic.is_local_group());

        let extragalactic = Distance::from_parsecs(10_000_000.0).unwrap();
        assert!(!extragalactic.is_galactic());
        assert!(!extragalactic.is_local_group());
    }

    #[test]
    fn test_proper_motion_distance() {
        // At 1 pc an arcsecond spans 1 au, by the definition of the parsec.
        let distance = Distance::from_parsecs(1.0).unwrap();
        assert_eq!(distance.proper_motion_distance_au(1.0, 1.0), 0.001);
    }

    #[test]
    fn test_arithmetic_operations() {
        let d1 = Distance::from_parsecs(10.0).unwrap();
        let d2 = Distance::from_parsecs(5.0).unwrap();

        let sum = (d1 + d2).unwrap();
        assert_eq!(sum.parsecs(), 15.0);

        let diff = (d1 - d2).unwrap();
        assert_eq!(diff.parsecs(), 5.0);

        let doubled = (d1 * 2.0).unwrap();
        assert_eq!(doubled.parsecs(), 20.0);

        let halved = (d1 / 2.0).unwrap();
        assert_eq!(halved.parsecs(), 5.0);
    }

    #[test]
    fn test_display() {
        let close = Distance::from_au(1.0).unwrap();
        assert!(close.to_string().contains("AU"));

        let nearby = Distance::from_parsecs(10.0).unwrap();
        assert!(nearby.to_string().contains("pc"));

        let distant = Distance::from_parsecs(10000.0).unwrap();
        assert!(distant.to_string().contains("kpc"));

        let very_distant = Distance::from_parsecs(10_000_000.0).unwrap();
        assert!(very_distant.to_string().contains("Mpc"));
    }
}
