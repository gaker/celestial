use celestial_core::AstroError;
use thiserror::Error;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub type CoordResult<T> = Result<T, CoordError>;

#[derive(Debug, Error)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum CoordError {
    #[error("Invalid coordinate: {message}")]
    InvalidCoordinate { message: String },

    #[error("Epoch conversion failed: {source}")]
    EpochError {
        #[from]
        source: celestial_time::TimeError,
    },

    #[error("Core astronomical calculation failed: {message}")]
    CoreError { message: String },

    #[error("Invalid distance: {message}")]
    InvalidDistance { message: String },

    #[error("Observer location required for topocentric coordinates")]
    MissingObserver,

    #[error("Coordinate operation not supported: {message}")]
    UnsupportedOperation { message: String },

    #[error("Data parsing failed: {message}")]
    ParsingError { message: String },

    #[error("Data not available: {message}")]
    DataUnavailable { message: String },

    /// Errors from external libraries (filesystem, network, etc.)
    ///
    /// This is deliberately unstructured (just a string) since external error types vary widely.
    /// If richer context is needed for specific external errors, add dedicated variants.
    #[error("External error: {message}")]
    ExternalError { message: String },
}

impl CoordError {
    pub fn invalid_coordinate(message: impl Into<String>) -> Self {
        Self::InvalidCoordinate {
            message: message.into(),
        }
    }

    pub fn invalid_distance(message: impl Into<String>) -> Self {
        Self::InvalidDistance {
            message: message.into(),
        }
    }

    pub fn unsupported_operation(message: impl Into<String>) -> Self {
        Self::UnsupportedOperation {
            message: message.into(),
        }
    }

    pub fn parsing_error(message: impl Into<String>) -> Self {
        Self::ParsingError {
            message: message.into(),
        }
    }

    pub fn data_unavailable(message: impl Into<String>) -> Self {
        Self::DataUnavailable {
            message: message.into(),
        }
    }

    pub fn external_library(operation: &str, error: &str) -> Self {
        Self::ExternalError {
            message: format!("{}: {}", operation, error),
        }
    }

    pub fn from_core(error: AstroError) -> Self {
        Self::CoreError {
            message: error.to_string(),
        }
    }
}

impl From<AstroError> for CoordError {
    fn from(error: AstroError) -> Self {
        Self::from_core(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unsupported_operation() {
        let err = CoordError::unsupported_operation("test op");
        assert!(err.to_string().contains("test op"));
    }

    #[test]
    fn test_parsing_error() {
        let err = CoordError::parsing_error("parse fail");
        assert!(err.to_string().contains("parse fail"));
    }
}
