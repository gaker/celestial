use crate::{CoordError, CoordResult};
use celestial_core::Angle;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

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
    /// Must be positive and finite (0 < ly < ∞)
    pub fn from_light_years(ly: f64) -> CoordResult<Self> {
        const LY_TO_PC: f64 = 0.3066013937;
        Self::from_parsecs(ly * LY_TO_PC)
    }

    /// Creates a Distance from astronomical units.
    ///
    /// # Valid Range
    /// Must be positive and finite (0 < au < ∞)
    pub fn from_au(au: f64) -> CoordResult<Self> {
        const AU_TO_PC: f64 = 4.84813681109536e-6;
        Self::from_parsecs(au * AU_TO_PC)
    }

    /// Creates a Distance from kilometers.
    ///
    /// # Valid Range
    /// Must be positive and finite (0 < km < ∞)
    pub fn from_kilometers(km: f64) -> CoordResult<Self> {
        const KM_TO_PC: f64 = 3.24077929e-14;
        Self::from_parsecs(km * KM_TO_PC)
    }

    /// Creates a Distance from parallax in arcseconds.
    ///
    /// # Valid Range
    /// Must be positive and finite (0 < parallax_arcsec < ∞)
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
        const PC_TO_LY: f64 = 3.2615637769;
        self.parsecs * PC_TO_LY
    }

    pub fn au(self) -> f64 {
        const PC_TO_AU: f64 = 206264.806247096;
        self.parsecs * PC_TO_AU
    }

    pub fn kilometers(self) -> f64 {
        #[allow(clippy::excessive_precision)]
        const PC_TO_KM: f64 = 3.0856775814913673e13;
        self.parsecs * PC_TO_KM
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
        let parsecs = 10.0_f64.powf((dm + 5.0) / 5.0);
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
        let pm_rad_per_year =
            pm_mas_per_year * 1e-3 * (celestial_core::constants::PI / (180.0 * 3600.0));
        let angular_distance_rad = pm_rad_per_year * dt_years;
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

    #[test]
    fn test_from_light_years() {
        let d = Distance::from_light_years(1.0).unwrap();
        assert!((d.parsecs() - 0.3066013937).abs() < 1e-9);
    }

    #[test]
    fn test_parallax_angle() {
        let angle = Angle::from_arcseconds(0.1);
        let d = Distance::from_parallax_angle(angle).unwrap();
        assert!((d.parsecs() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_parallax_uncertainty_mas() {
        let d = Distance::from_parsecs(100.0).unwrap();
        let unc = d.parallax_uncertainty_mas(0.01);
        assert!((unc - 0.1).abs() < 1e-6);
    }

    #[test]
    fn test_partial_ord() {
        let d1 = Distance::from_parsecs(10.0).unwrap();
        let d2 = Distance::from_parsecs(20.0).unwrap();
        assert!(d1 < d2);
    }

    #[test]
    fn test_unit_conversions() {
        let distance = Distance::from_parsecs(1.0).unwrap();

        #[allow(clippy::excessive_precision)]
        {
            assert!((distance.light_years() - 3.261_563_776_9).abs() < 1e-9);
            assert!((distance.au() - 206264.806_247_096).abs() < 1e-6);
            assert!((distance.kilometers() - 3.085_677_581_491_367_3e13).abs() < 1e6);
        }
    }

    #[test]
    fn test_parallax_calculations() {
        let proxima = Distance::from_parallax_arcsec(0.7687).unwrap();
        assert!((proxima.parsecs() - 1.3009).abs() < 0.001);

        let distance = Distance::from_parallax_milliarcsec(768.7).unwrap();
        assert!((distance.parsecs() - 1.3009).abs() < 0.001);
    }

    #[test]
    fn test_distance_modulus() {
        let distance = Distance::from_parsecs(10.0).unwrap();
        let dm = distance.distance_modulus();
        assert!((dm - 0.0).abs() < 1e-12);

        let recovered = Distance::from_distance_modulus(dm).unwrap();
        assert!((recovered.parsecs() - 10.0).abs() < 1e-12);
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
        let distance = Distance::from_parsecs(1.0).unwrap();

        let linear_dist = distance.proper_motion_distance_au(1.0, 1.0);

        assert!(linear_dist > 0.0);
        assert!(linear_dist < 10.0);
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
