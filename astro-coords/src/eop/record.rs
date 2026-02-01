use crate::{CoordError, CoordResult};
use astro_core::constants::{
    ARCSEC_TO_RAD, DAYS_PER_JULIAN_CENTURY, J2000_JD, MILLIARCSEC_TO_RAD, MJD_ZERO_POINT,
};
use astro_time::{transforms::earth_rotation_angle, JulianDate};

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
        if x_p_arcsec.abs() > 6.0 {
            return Err(CoordError::invalid_coordinate(format!(
                "X polar motion out of range: {} arcsec",
                x_p_arcsec
            )));
        }

        if y_p_arcsec.abs() > 6.0 {
            return Err(CoordError::invalid_coordinate(format!(
                "Y polar motion out of range: {} arcsec",
                y_p_arcsec
            )));
        }

        if ut1_utc_sec.abs() > 1.0 {
            return Err(CoordError::invalid_coordinate(format!(
                "UT1-UTC out of range: {} sec",
                ut1_utc_sec
            )));
        }

        if lod_sec.abs() > 0.01 {
            return Err(CoordError::invalid_coordinate(format!(
                "LOD out of range: {} sec",
                lod_sec
            )));
        }

        Ok(Self {
            mjd,
            x_p_encoded: (x_p_arcsec * Self::ARCSEC_TO_UNITS).round() as i32,
            y_p_encoded: (y_p_arcsec * Self::ARCSEC_TO_UNITS).round() as i32,
            ut1_utc_encoded: (ut1_utc_sec * Self::SEC_TO_UNITS).round() as i32,
            lod_encoded: (lod_sec * Self::SEC_TO_UNITS).round() as i32,
            dx_encoded: None,
            dy_encoded: None,
            flags: EopFlags::default(),
        })
    }

    pub fn with_cip_offsets(
        mut self,
        dx_milliarcsec: f64,
        dy_milliarcsec: f64,
    ) -> CoordResult<Self> {
        if dx_milliarcsec.abs() > 1000.0 {
            return Err(CoordError::invalid_coordinate(format!(
                "CIP dX out of range: {} mas",
                dx_milliarcsec
            )));
        }

        if dy_milliarcsec.abs() > 1000.0 {
            return Err(CoordError::invalid_coordinate(format!(
                "CIP dY out of range: {} mas",
                dy_milliarcsec
            )));
        }

        self.dx_encoded = Some((dx_milliarcsec * Self::MILLIARCSEC_TO_UNITS).round() as i16);
        self.dy_encoded = Some((dy_milliarcsec * Self::MILLIARCSEC_TO_UNITS).round() as i16);
        self.flags.has_cip_offsets = true;

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
        let ut1_jd =
            self.mjd + MJD_ZERO_POINT + self.ut1_utc / astro_core::constants::SECONDS_PER_DAY_F64;

        let ut1_jd1 = ut1_jd.floor();
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
}
