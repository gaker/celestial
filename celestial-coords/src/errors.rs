use celestial_core::errors::AstroError;
use thiserror::Error;

pub type CoordResult<T> = Result<T, CoordError>;

#[derive(Debug, Error)]
pub enum CoordError {
    #[error("Invalid coordinate: {message}")]
    InvalidCoordinate { message: String },

    #[error("Epoch conversion failed: {0}")]
    EpochError(#[from] celestial_time::TimeError),

    #[error("Core astronomical calculation failed: {0}")]
    CoreError(#[from] AstroError),

    #[error("Invalid distance: {message}")]
    InvalidDistance { message: String },

    #[error("Coordinate operation not supported: {message}")]
    UnsupportedOperation { message: String },

    #[error("Data parsing failed: {message}")]
    ParsingError { message: String },

    #[error("Data not available: {message}")]
    DataUnavailable { message: String },

    #[error("{context}: {source}")]
    Io {
        context: String,
        source: std::io::Error,
    },
}

impl CoordError {
    pub(crate) fn invalid_coordinate(message: impl Into<String>) -> Self {
        Self::InvalidCoordinate {
            message: message.into(),
        }
    }

    pub(crate) fn invalid_distance(message: impl Into<String>) -> Self {
        Self::InvalidDistance {
            message: message.into(),
        }
    }

    pub fn unsupported_operation(message: impl Into<String>) -> Self {
        Self::UnsupportedOperation {
            message: message.into(),
        }
    }

    pub(crate) fn parsing_error(message: impl Into<String>) -> Self {
        Self::ParsingError {
            message: message.into(),
        }
    }

    pub(crate) fn data_unavailable(message: impl Into<String>) -> Self {
        Self::DataUnavailable {
            message: message.into(),
        }
    }

    pub(crate) fn io(context: impl Into<String>, source: std::io::Error) -> Self {
        Self::Io {
            context: context.into(),
            source,
        }
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

    #[test]
    fn test_core_error_keeps_its_source() {
        let err = CoordError::from(AstroError::calculation_error("test", "failed"));
        let source = std::error::Error::source(&err);
        assert!(source.is_some_and(|s| s.is::<AstroError>()), "{err:?}");
    }
}
