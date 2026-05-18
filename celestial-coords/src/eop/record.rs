use crate::{CoordError, CoordResult};
use celestial_core::constants::{
    ARCSEC_TO_RAD, DAYS_PER_JULIAN_CENTURY, J2000_JD, MILLIARCSEC_TO_RAD, MJD_ZERO_POINT,
};
use celestial_time::{transforms::earth_rotation_angle, JulianDate};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EopRecord {
    pub mjd: f64,

    pub x_p_encoded: i32,

    pub y_p_encoded: i32,

    pub ut1_utc_encoded: i32,

    pub lod_encoded: i32,

    pub dx_encoded: Option<i16>,

    pub dy_encoded: Option<i16>,

    pub xrt_encoded: Option<i32>,

    pub yrt_encoded: Option<i32>,

    pub flags: EopFlags,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EopFlags {
    pub source: EopSource,

    pub quality: EopQuality,

    pub has_polar_motion: bool,

    pub has_ut1_utc: bool,

    pub has_cip_offsets: bool,

    pub has_pole_rates: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EopSource {
    IersC04,

    IersFinals,

    IersPrediction,

    UserData,

    Interpolated,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EopQuality {
    HighPrecision,

    Standard,

    LowPrecision,

    Predicted,
}

fn validate_range(value: f64, limit: f64, label: &str, unit: &str) -> CoordResult<()> {
    if value.abs() > limit {
        return Err(CoordError::invalid_coordinate(format!(
            "{} out of range: {} {}", label, value, unit,
        )));
    }
    Ok(())
}

impl EopRecord {
    const ARCSEC_TO_UNITS: f64 = 10_000_000.0;
    const SEC_TO_UNITS: f64 = 10_000_000.0;
    const MILLIARCSEC_TO_UNITS: f64 = 10_000.0;

    pub fn new(
        mjd: f64,
        x_p_arcsec: f64,
        y_p_arcsec: f64,
        ut1_utc_sec: f64,
        lod_sec: f64,
    ) -> CoordResult<Self> {
        validate_range(x_p_arcsec, 6.0, "X polar motion", "arcsec")?;
        validate_range(y_p_arcsec, 6.0, "Y polar motion", "arcsec")?;
        validate_range(ut1_utc_sec, 1.0, "UT1-UTC", "sec")?;
        validate_range(lod_sec, 0.01, "LOD", "sec")?;
        Ok(Self {
            mjd,
            x_p_encoded: libm::round(x_p_arcsec * Self::ARCSEC_TO_UNITS) as i32,
            y_p_encoded: libm::round(y_p_arcsec * Self::ARCSEC_TO_UNITS) as i32,
            ut1_utc_encoded: libm::round(ut1_utc_sec * Self::SEC_TO_UNITS) as i32,
            lod_encoded: libm::round(lod_sec * Self::SEC_TO_UNITS) as i32,
            dx_encoded: None,
            dy_encoded: None,
            xrt_encoded: None,
            yrt_encoded: None,
            flags: EopFlags::default(),
        })
    }

    pub fn with_cip_offsets(
        mut self,
        dx_milliarcsec: f64,
        dy_milliarcsec: f64,
    ) -> CoordResult<Self> {
        validate_range(dx_milliarcsec, 1000.0, "CIP dX", "mas")?;
        validate_range(dy_milliarcsec, 1000.0, "CIP dY", "mas")?;
        self.dx_encoded = Some(libm::round(dx_milliarcsec * Self::MILLIARCSEC_TO_UNITS) as i16);
        self.dy_encoded = Some(libm::round(dy_milliarcsec * Self::MILLIARCSEC_TO_UNITS) as i16);
        self.flags.has_cip_offsets = true;
        Ok(self)
    }

    pub fn with_pole_rates(
        mut self,
        xrt_arcsec_per_day: f64,
        yrt_arcsec_per_day: f64,
    ) -> CoordResult<Self> {
        validate_range(xrt_arcsec_per_day, 1.0, "Pole rate xrt", "arcsec/day")?;
        validate_range(yrt_arcsec_per_day, 1.0, "Pole rate yrt", "arcsec/day")?;
        self.xrt_encoded = Some(libm::round(xrt_arcsec_per_day * Self::ARCSEC_TO_UNITS) as i32);
        self.yrt_encoded = Some(libm::round(yrt_arcsec_per_day * Self::ARCSEC_TO_UNITS) as i32);
        self.flags.has_pole_rates = true;
        Ok(self)
    }

    pub fn with_flags(mut self, flags: EopFlags) -> Self {
        self.flags = flags;
        self
    }

    pub fn to_parameters(&self) -> EopParameters {
        EopParameters {
            mjd: self.mjd,
            x_p: self.x_p_encoded as f64 / Self::ARCSEC_TO_UNITS,
            y_p: self.y_p_encoded as f64 / Self::ARCSEC_TO_UNITS,
            ut1_utc: self.ut1_utc_encoded as f64 / Self::SEC_TO_UNITS,
            lod: self.lod_encoded as f64 / Self::SEC_TO_UNITS,
            dx: self
                .dx_encoded
                .map(|v| v as f64 / Self::MILLIARCSEC_TO_UNITS),
            dy: self
                .dy_encoded
                .map(|v| v as f64 / Self::MILLIARCSEC_TO_UNITS),
            xrt: self.xrt_encoded.map(|v| v as f64 / Self::ARCSEC_TO_UNITS),
            yrt: self.yrt_encoded.map(|v| v as f64 / Self::ARCSEC_TO_UNITS),
            s_prime: 0.0,
            flags: self.flags,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct EopParameters {
    pub mjd: f64,

    pub x_p: f64,

    pub y_p: f64,

    pub ut1_utc: f64,

    pub lod: f64,

    pub dx: Option<f64>,

    pub dy: Option<f64>,

    pub xrt: Option<f64>,

    pub yrt: Option<f64>,

    pub s_prime: f64,

    pub flags: EopFlags,
}

impl EopParameters {
    /// Computes the TIO locator s' using the IAU 2000 linear approximation.
    ///
    /// s' ≈ -47 µas/century × t, where t is Julian centuries from J2000.0.
    /// Result is in radians (matching IAU SOFA convention).
    ///
    /// Note: For maximum precision, use `compute_s_prime_jd()` with 2-part JD.
    pub fn compute_s_prime(&mut self) {
        self.compute_s_prime_jd(self.mjd + MJD_ZERO_POINT, 0.0);
    }

    /// Computes s' from 2-part TT Julian Date for maximum precision.
    ///
    /// The 2-part JD allows for precision-preserving arithmetic when
    /// computing time intervals from J2000.0.
    pub fn compute_s_prime_jd(&mut self, tt1: f64, tt2: f64) {
        let t = ((tt1 - J2000_JD) + tt2) / DAYS_PER_JULIAN_CENTURY;
        self.s_prime = -47e-6 * t * ARCSEC_TO_RAD;
    }

    pub fn compute_era(&self) -> CoordResult<f64> {
        let ut1_jd = self.mjd
            + MJD_ZERO_POINT
            + self.ut1_utc / celestial_core::constants::SECONDS_PER_DAY_F64;

        let ut1_jd1 = libm::floor(ut1_jd);
        let ut1_jd2 = ut1_jd - ut1_jd1;

        let jd = JulianDate::new(ut1_jd1, ut1_jd2);
        earth_rotation_angle(&jd)
            .map_err(|e| CoordError::external_library("Earth Rotation Angle", &e.to_string()))
    }

    /// Returns CIP X coordinate corrected by dX offset (if available), in radians.
    pub fn corrected_cip_x(&self, x_iau: f64) -> f64 {
        x_iau + self.dx.unwrap_or(0.0) * MILLIARCSEC_TO_RAD
    }

    /// Returns CIP Y coordinate corrected by dY offset (if available), in radians.
    pub fn corrected_cip_y(&self, y_iau: f64) -> f64 {
        y_iau + self.dy.unwrap_or(0.0) * MILLIARCSEC_TO_RAD
    }
}

impl Default for EopFlags {
    fn default() -> Self {
        Self {
            source: EopSource::UserData,
            quality: EopQuality::Standard,
            has_polar_motion: true,
            has_ut1_utc: true,
            has_cip_offsets: false,
            has_pole_rates: false,
        }
    }
}

impl std::fmt::Display for EopParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "EOP(MJD={:.1}, xp={:.6}\", yp={:.6}\", UT1-UTC={:.7}s)",
            self.mjd, self.x_p, self.y_p, self.ut1_utc
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eop_record_encoding() {
        let record = EopRecord::new(
            59945.0,   // MJD 2023-01-01
            0.123456,  // x_p arcsec
            0.234567,  // y_p arcsec
            0.0123456, // UT1-UTC sec
            0.0012345, // LOD sec
        )
        .unwrap();

        let params = record.to_parameters();

        // Check precision preservation (should be within 0.1 µas for angles)
        assert!((params.x_p - 0.123456).abs() == 0.0);
        assert!((params.y_p - 0.234567).abs() == 0.0);
        assert!((params.ut1_utc - 0.0123456).abs() == 0.0);
        assert!((params.lod - 0.0012345).abs() == 0.0);
        assert_eq!(params.mjd, 59945.0);
    }

    #[test]
    fn test_cip_offsets() {
        let record = EopRecord::new(59945.0, 0.1, 0.2, 0.01, 0.001)
            .unwrap()
            .with_cip_offsets(0.5, -0.3) // mas
            .unwrap();

        let params = record.to_parameters();

        assert_eq!(params.dx, Some(0.5));
        assert_eq!(params.dy, Some(-0.3));
        assert!(params.flags.has_cip_offsets);
    }

    #[test]
    fn test_parameter_display() {
        let mut params = EopParameters {
            mjd: 59945.0,
            x_p: 0.123456,
            y_p: 0.234567,
            ut1_utc: 0.0123456,
            lod: 0.001,
            dx: None,
            dy: None,
            xrt: None,
            yrt: None,
            s_prime: 0.0,
            flags: EopFlags::default(),
        };

        params.compute_s_prime();

        let display = format!("{}", params);
        assert!(display.contains("MJD=59945.0"));
        assert!(display.contains("xp=0.123456"));
        assert!(display.contains("yp=0.234567"));
    }

    #[test]
    fn test_validation_x_polar_motion_out_of_range() {
        let result = EopRecord::new(59945.0, 6.1, 0.2, 0.01, 0.001);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("X polar motion out of range"));
    }

    #[test]
    fn test_validation_y_polar_motion_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, -6.1, 0.01, 0.001);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Y polar motion out of range"));
    }

    #[test]
    fn test_validation_ut1_utc_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 1.1, 0.001);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("UT1-UTC out of range"));
    }

    #[test]
    fn test_validation_lod_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01, 0.011);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("LOD out of range"));
    }

    #[test]
    fn test_validation_cip_dx_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01, 0.001)
            .unwrap()
            .with_cip_offsets(1001.0, 0.0);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("CIP dX out of range"));
    }

    #[test]
    fn test_validation_cip_dy_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01, 0.001)
            .unwrap()
            .with_cip_offsets(0.0, -1001.0);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("CIP dY out of range"));
    }

    #[test]
    fn with_pole_rates_round_trips_through_to_parameters() {
        let record = EopRecord::new(59945.0, 0.1, 0.2, 0.01, 0.001)
            .unwrap()
            .with_pole_rates(0.0001, -0.0002)
            .unwrap();
        let params = record.to_parameters();
        assert!((params.xrt.unwrap() - 0.0001).abs() < 1e-9);
        assert!((params.yrt.unwrap() - (-0.0002)).abs() < 1e-9);
        assert!(params.flags.has_pole_rates);
    }

    #[test]
    fn with_pole_rates_xrt_out_of_range_errors() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01, 0.001)
            .unwrap()
            .with_pole_rates(1.5, 0.0);
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Pole rate xrt out of range"));
    }

    #[test]
    fn with_pole_rates_yrt_out_of_range_errors() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01, 0.001)
            .unwrap()
            .with_pole_rates(0.0, -1.5);
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Pole rate yrt out of range"));
    }

    #[test]
    fn with_flags_overrides_default_flags() {
        let custom = EopFlags {
            source: EopSource::IersC04,
            quality: EopQuality::Predicted,
            has_polar_motion: false,
            has_ut1_utc: false,
            has_cip_offsets: true,
            has_pole_rates: true,
        };
        let record = EopRecord::new(59945.0, 0.1, 0.2, 0.01, 0.001)
            .unwrap()
            .with_flags(custom);
        assert_eq!(record.flags.source, EopSource::IersC04);
        assert_eq!(record.flags.quality, EopQuality::Predicted);
        assert!(!record.flags.has_polar_motion);
    }

    #[test]
    fn eop_flags_default_values() {
        let f = EopFlags::default();
        assert_eq!(f.source, EopSource::UserData);
        assert_eq!(f.quality, EopQuality::Standard);
        assert!(f.has_polar_motion);
        assert!(f.has_ut1_utc);
        assert!(!f.has_cip_offsets);
        assert!(!f.has_pole_rates);
    }

    fn sample_params(mjd: f64) -> EopParameters {
        EopParameters {
            mjd,
            x_p: 0.0, y_p: 0.0,
            ut1_utc: 0.0, lod: 0.0,
            dx: None, dy: None, xrt: None, yrt: None,
            s_prime: 0.0,
            flags: EopFlags::default(),
        }
    }

    #[test]
    fn compute_s_prime_jd_matches_compute_s_prime_at_j2000() {
        // At J2000 (MJD 51544.5), t=0 → s_prime should be 0.
        let mut params = sample_params(51544.5);
        params.compute_s_prime();
        assert!(params.s_prime.abs() < 1e-15);

        // compute_s_prime_jd with the same epoch via 2-part JD must agree.
        let mut params2 = sample_params(51544.5);
        params2.compute_s_prime_jd(2451545.0, 0.0);
        assert!((params.s_prime - params2.s_prime).abs() < 1e-15);
    }

    #[test]
    fn compute_s_prime_jd_is_linear_in_centuries() {
        // 1 century after J2000 → s' = -47e-6 * ARCSEC_TO_RAD radians.
        let mut params = sample_params(0.0);
        params.compute_s_prime_jd(2451545.0 + 36525.0, 0.0);
        let expected = -47e-6 * ARCSEC_TO_RAD;
        assert!((params.s_prime - expected).abs() < 1e-15);
    }

    #[test]
    fn compute_era_returns_radians_in_range() {
        // ERA wraps in [0, 2π); just verify it returns something in range
        // and is sensitive to mjd.
        let params_a = sample_params(60000.0);
        let params_b = sample_params(60000.5);
        let era_a = params_a.compute_era().unwrap();
        let era_b = params_b.compute_era().unwrap();
        assert!((0.0..celestial_core::constants::TWOPI).contains(&era_a));
        assert!((0.0..celestial_core::constants::TWOPI).contains(&era_b));
        assert!(era_a != era_b);
    }

    #[test]
    fn corrected_cip_x_adds_dx_offset_when_present() {
        let mut params = sample_params(60000.0);
        params.dx = Some(1.0); // 1 mas
        let corrected = params.corrected_cip_x(0.5);
        // Expected: 0.5 + 1.0 * MILLIARCSEC_TO_RAD
        let expected = 0.5 + MILLIARCSEC_TO_RAD;
        assert!((corrected - expected).abs() < 1e-15);
    }

    #[test]
    fn corrected_cip_y_falls_back_to_zero_when_dy_is_none() {
        let params = sample_params(60000.0);
        assert_eq!(params.dy, None);
        // No dy → corrected == input.
        assert_eq!(params.corrected_cip_y(0.7), 0.7);
    }

    #[test]
    fn corrected_cip_x_passthrough_when_dx_is_none() {
        let params = sample_params(60000.0);
        assert_eq!(params.dx, None);
        assert_eq!(params.corrected_cip_x(-0.3), -0.3);
    }
}
