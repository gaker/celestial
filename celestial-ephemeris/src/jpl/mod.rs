mod chain;
mod chebyshev;
mod daf;
pub mod spk;
#[cfg(test)]
mod tests;
mod type2;

use celestial_core::errors::{AstroError, MathErrorKind};

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SpkError {
    #[error("SPK read failed: {0}")]
    Io(#[from] std::io::Error),
    #[error("invalid SPK file: {0}")]
    InvalidFormat(String),
    #[error("invalid SPK data: {0}")]
    InvalidData(String),
    #[error("TDB epoch JD {jd} is not a finite number of seconds from J2000")]
    InvalidEpoch { jd: f64 },
    #[error("no SPK segments connect body {body} to {center} at TDB JD {jd}")]
    SegmentNotFound { body: i32, center: i32, jd: f64 },
    #[error("SPK data type {0} is not supported; only type 2 is")]
    UnsupportedType(i32),
    #[error("SPK frame {0} is not supported; only frame 1 (J2000) is")]
    UnsupportedFrame(i32),
}

impl From<SpkError> for AstroError {
    fn from(error: SpkError) -> Self {
        let reason = error.to_string();
        match error {
            SpkError::Io(_) => AstroError::data_error("SPK", "read", &reason),
            SpkError::InvalidFormat(_) => AstroError::data_error("SPK", "open", &reason),
            SpkError::InvalidEpoch { .. } => {
                AstroError::math_error("SPK", MathErrorKind::NotFinite, &reason)
            }
            _ => AstroError::data_error("SPK", "evaluate", &reason),
        }
    }
}

pub mod bodies {
    pub const SOLAR_SYSTEM_BARYCENTER: i32 = 0;
    pub const MERCURY_BARYCENTER: i32 = 1;
    pub const VENUS_BARYCENTER: i32 = 2;
    pub const EARTH_MOON_BARYCENTER: i32 = 3;
    pub const MARS_BARYCENTER: i32 = 4;
    pub const JUPITER_BARYCENTER: i32 = 5;
    pub const SATURN_BARYCENTER: i32 = 6;
    pub const URANUS_BARYCENTER: i32 = 7;
    pub const NEPTUNE_BARYCENTER: i32 = 8;
    pub const PLUTO_BARYCENTER: i32 = 9;
    pub const SUN: i32 = 10;
    pub const MERCURY: i32 = 199;
    pub const VENUS: i32 = 299;
    pub const MOON: i32 = 301;
    pub const EARTH: i32 = 399;
}
