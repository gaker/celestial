pub mod constants;
pub mod julian;
pub mod parsing;
pub mod scales;
pub mod sidereal;
pub mod transforms;

pub use julian::JulianDate;
pub use scales::{GPS, TAI, TCB, TCG, TDB, TT, UT1, UTC};

pub use scales::{
    gps_from_calendar, tai_from_calendar, tcb_from_calendar, tcg_from_calendar, tdb_from_calendar,
    tt_from_calendar, ut1_from_calendar, utc_from_calendar,
};

pub use scales::conversions::{
    TcbToTdb, TdbToTcb, ToGPS, ToTAI, ToTAIWithOffset, ToTCB, ToTCG, ToTCGFromTCB, ToTDB, ToTT,
    ToTTFromTDB, ToTTWithDeltaT, ToUT1, ToUT1WithDUT1, ToUT1WithDeltaT, ToUT1WithOffset, ToUTC,
    ToUTCViaTAI, ToUTCWithDUT1,
};
pub use sidereal::{ObservatoryContext, SiderealAngle, GAST, GMST, LAST, LMST};
pub use transforms::{NutationCalculator, NutationModel, NutationResult};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub type TimeResult<T> = Result<T, TimeError>;

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum TimeError {
    InvalidDate,
    ConversionError(String),
    ParseError(String),
    CalculationError(String),
    InvalidEpoch(String),
}

impl std::fmt::Display for TimeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TimeError::InvalidDate => write!(f, "Invalid date"),
            TimeError::ConversionError(msg) => write!(f, "Conversion error: {}", msg),
            TimeError::ParseError(msg) => write!(f, "Parse error: {}", msg),
            TimeError::CalculationError(msg) => write!(f, "Calculation error: {}", msg),
            TimeError::InvalidEpoch(msg) => write!(f, "Invalid epoch: {}", msg),
        }
    }
}

impl std::error::Error for TimeError {}

impl From<astro_core::AstroError> for TimeError {
    fn from(err: astro_core::AstroError) -> Self {
        TimeError::CalculationError(err.to_string())
    }
}
