pub mod constants;
pub mod julian;
pub(crate) mod parsing;
pub mod scales;
pub mod sidereal;
pub mod transforms;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use thiserror::Error;

pub type TimeResult<T> = Result<T, TimeError>;

#[derive(Error, Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum TimeError {
    #[error("Invalid date: {0}")]
    InvalidDate(String),
    #[error("Conversion error: {0}")]
    ConversionError(String),
    #[error("Parse error: {0}")]
    ParseError(String),
    #[error("Calculation error: {0}")]
    CalculationError(String),
    #[error("Invalid epoch: {0}")]
    InvalidEpoch(String),
}

impl From<celestial_core::errors::AstroError> for TimeError {
    fn from(err: celestial_core::errors::AstroError) -> Self {
        Self::CalculationError(err.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        for (error, shown) in [
            (TimeError::InvalidDate("a".into()), "Invalid date: a"),
            (
                TimeError::ConversionError("b".into()),
                "Conversion error: b",
            ),
            (TimeError::ParseError("c".into()), "Parse error: c"),
            (
                TimeError::CalculationError("d".into()),
                "Calculation error: d",
            ),
            (TimeError::InvalidEpoch("e".into()), "Invalid epoch: e"),
        ] {
            assert_eq!(error.to_string(), shown);
        }
    }
}
