use crate::errors::{CoordError, CoordResult};
use celestial_core::constants::MILLIARCSEC_TO_RAD;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EopRecord {
    pub(crate) mjd: f64,

    pub(crate) x_p_encoded: i32,

    pub(crate) y_p_encoded: i32,

    pub(crate) ut1_utc_encoded: i32,

    pub(crate) lod_encoded: Option<i32>,

    pub(crate) dx_encoded: Option<i32>,

    pub(crate) dy_encoded: Option<i32>,

    pub(crate) xrt_encoded: Option<i32>,

    pub(crate) yrt_encoded: Option<i32>,

    pub(crate) flags: EopFlags,
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

    UserData,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EopQuality {
    HighPrecision,

    Standard,

    Predicted,
}

fn validate_range(value: f64, limit: f64, label: &str, unit: &str) -> CoordResult<()> {
    if !value.is_finite() || libm::fabs(value) > limit {
        return Err(CoordError::invalid_coordinate(format!(
            "{} out of range: {} {}",
            label, value, unit,
        )));
    }
    Ok(())
}

impl EopRecord {
    const ARCSEC_TO_UNITS: f64 = 10_000_000.0;
    const SEC_TO_UNITS: f64 = 10_000_000.0;
    const MILLIARCSEC_TO_UNITS: f64 = 10_000.0;

    pub fn new(mjd: f64, x_p_arcsec: f64, y_p_arcsec: f64, ut1_utc_sec: f64) -> CoordResult<Self> {
        if !mjd.is_finite() {
            return Err(CoordError::invalid_coordinate(format!(
                "EOP record MJD {} is not finite",
                mjd
            )));
        }
        validate_range(x_p_arcsec, 6.0, "X polar motion", "arcsec")?;
        validate_range(y_p_arcsec, 6.0, "Y polar motion", "arcsec")?;
        validate_range(ut1_utc_sec, 1.0, "UT1-UTC", "sec")?;
        Ok(Self {
            mjd,
            x_p_encoded: libm::round(x_p_arcsec * Self::ARCSEC_TO_UNITS) as i32,
            y_p_encoded: libm::round(y_p_arcsec * Self::ARCSEC_TO_UNITS) as i32,
            ut1_utc_encoded: libm::round(ut1_utc_sec * Self::SEC_TO_UNITS) as i32,
            lod_encoded: None,
            dx_encoded: None,
            dy_encoded: None,
            xrt_encoded: None,
            yrt_encoded: None,
            flags: EopFlags::default(),
        })
    }

    pub fn with_lod(mut self, lod_sec: f64) -> CoordResult<Self> {
        validate_range(lod_sec, 0.01, "LOD", "sec")?;
        self.lod_encoded = Some(libm::round(lod_sec * Self::SEC_TO_UNITS) as i32);
        Ok(self)
    }

    pub fn with_cip_offsets(
        mut self,
        dx_milliarcsec: f64,
        dy_milliarcsec: f64,
    ) -> CoordResult<Self> {
        validate_range(dx_milliarcsec, 1000.0, "CIP dX", "mas")?;
        validate_range(dy_milliarcsec, 1000.0, "CIP dY", "mas")?;
        self.dx_encoded = Some(libm::round(dx_milliarcsec * Self::MILLIARCSEC_TO_UNITS) as i32);
        self.dy_encoded = Some(libm::round(dy_milliarcsec * Self::MILLIARCSEC_TO_UNITS) as i32);
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

    pub(crate) fn with_flags(mut self, flags: EopFlags) -> Self {
        self.flags = flags;
        self
    }

    pub fn to_parameters(&self) -> EopParameters {
        EopParameters {
            mjd: self.mjd,
            x_p: self.x_p_encoded as f64 / Self::ARCSEC_TO_UNITS,
            y_p: self.y_p_encoded as f64 / Self::ARCSEC_TO_UNITS,
            ut1_utc: self.ut1_utc_encoded as f64 / Self::SEC_TO_UNITS,
            lod: self.lod_encoded.map(|v| v as f64 / Self::SEC_TO_UNITS),
            dx: self
                .dx_encoded
                .map(|v| v as f64 / Self::MILLIARCSEC_TO_UNITS),
            dy: self
                .dy_encoded
                .map(|v| v as f64 / Self::MILLIARCSEC_TO_UNITS),
            xrt: self.xrt_encoded.map(|v| v as f64 / Self::ARCSEC_TO_UNITS),
            yrt: self.yrt_encoded.map(|v| v as f64 / Self::ARCSEC_TO_UNITS),
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

    pub lod: Option<f64>,

    pub dx: Option<f64>,

    pub dy: Option<f64>,

    pub xrt: Option<f64>,

    pub yrt: Option<f64>,

    pub flags: EopFlags,
}

impl EopParameters {
    /// Returns CIP X coordinate corrected by dX offset (if available), in radians.
    pub(crate) fn corrected_cip_x(&self, x_iau: f64) -> f64 {
        x_iau + self.dx.unwrap_or(0.0) * MILLIARCSEC_TO_RAD
    }

    /// Returns CIP Y coordinate corrected by dY offset (if available), in radians.
    pub(crate) fn corrected_cip_y(&self, y_iau: f64) -> f64 {
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
        )
        .unwrap()
        .with_lod(0.0012345)
        .unwrap();

        let params = record.to_parameters();

        // Stored as whole multiples of 1e-7, so values given to seven decimals or fewer come
        // back exactly.
        assert_eq!(params.x_p, 0.123456);
        assert_eq!(params.y_p, 0.234567);
        assert_eq!(params.ut1_utc, 0.0123456);
        assert_eq!(params.lod, Some(0.0012345));
        assert_eq!(params.mjd, 59945.0);
    }

    #[test]
    fn test_cip_offsets() {
        let record = EopRecord::new(59945.0, 0.1, 0.2, 0.01)
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
        let params = EopParameters {
            mjd: 59945.0,
            x_p: 0.123456,
            y_p: 0.234567,
            ut1_utc: 0.0123456,
            lod: Some(0.001),
            dx: None,
            dy: None,
            xrt: None,
            yrt: None,
            flags: EopFlags::default(),
        };

        let display = format!("{}", params);
        assert!(display.contains("MJD=59945.0"));
        assert!(display.contains("xp=0.123456"));
        assert!(display.contains("yp=0.234567"));
    }

    #[test]
    fn test_validation_x_polar_motion_out_of_range() {
        let result = EopRecord::new(59945.0, 6.1, 0.2, 0.01);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("X polar motion out of range"));
    }

    #[test]
    fn test_validation_y_polar_motion_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, -6.1, 0.01);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Y polar motion out of range"));
    }

    #[test]
    fn test_validation_ut1_utc_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 1.1);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("UT1-UTC out of range"));
    }

    #[test]
    fn test_validation_lod_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01)
            .unwrap()
            .with_lod(0.011);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("LOD out of range"));
    }

    #[test]
    fn test_validation_cip_dx_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01)
            .unwrap()
            .with_cip_offsets(1001.0, 0.0);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("CIP dX out of range"));
    }

    #[test]
    fn test_validation_cip_dy_out_of_range() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01)
            .unwrap()
            .with_cip_offsets(0.0, -1001.0);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("CIP dY out of range"));
    }

    #[test]
    fn with_pole_rates_round_trips_through_to_parameters() {
        let record = EopRecord::new(59945.0, 0.1, 0.2, 0.01)
            .unwrap()
            .with_pole_rates(0.0001, -0.0002)
            .unwrap();
        let params = record.to_parameters();
        assert_eq!((params.xrt, params.yrt), (Some(0.0001), Some(-0.0002)));
        assert!(params.flags.has_pole_rates);
    }

    #[test]
    fn with_pole_rates_xrt_out_of_range_errors() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01)
            .unwrap()
            .with_pole_rates(1.5, 0.0);
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Pole rate xrt out of range"));
    }

    #[test]
    fn with_pole_rates_yrt_out_of_range_errors() {
        let result = EopRecord::new(59945.0, 0.1, 0.2, 0.01)
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
        let record = EopRecord::new(59945.0, 0.1, 0.2, 0.01)
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
            x_p: 0.0,
            y_p: 0.0,
            ut1_utc: 0.0,
            lod: None,
            dx: None,
            dy: None,
            xrt: None,
            yrt: None,
            flags: EopFlags::default(),
        }
    }

    #[test]
    fn corrected_cip_x_adds_dx_offset_when_present() {
        let mut params = sample_params(60000.0);
        params.dx = Some(1.0);
        assert_eq!(params.corrected_cip_x(0.5), 0.5 + MILLIARCSEC_TO_RAD);
    }

    #[test]
    fn corrected_cip_y_falls_back_to_zero_when_dy_is_none() {
        let params = sample_params(60000.0);
        assert_eq!(params.dy, None);
        assert_eq!(params.corrected_cip_y(0.7), 0.7);
    }

    #[test]
    fn corrected_cip_x_passthrough_when_dx_is_none() {
        let params = sample_params(60000.0);
        assert_eq!(params.dx, None);
        assert_eq!(params.corrected_cip_x(-0.3), -0.3);
    }

    #[test]
    fn non_finite_inputs_are_rejected() {
        assert!(EopRecord::new(f64::NAN, 0.1, 0.2, 0.01).is_err());
        assert!(EopRecord::new(59945.0, f64::NAN, 0.2, 0.01).is_err());
        assert!(EopRecord::new(59945.0, 0.1, 0.2, f64::INFINITY).is_err());
        let record = EopRecord::new(59945.0, 0.1, 0.2, 0.01).unwrap();
        assert!(record.clone().with_lod(f64::NAN).is_err());
        assert!(record.with_cip_offsets(f64::NAN, 0.0).is_err());
    }

    // MJD 45700 in C04 has dX = 3.287 mas, past what an i16 at 1e-4 mas can hold.
    #[test]
    fn large_cip_offsets_round_trip() {
        let params = EopRecord::new(45700.0, 0.1, 0.2, 0.01)
            .unwrap()
            .with_cip_offsets(8.661, -900.0)
            .unwrap()
            .to_parameters();
        assert_eq!((params.dx, params.dy), (Some(8.661), Some(-900.0)));
    }
}
