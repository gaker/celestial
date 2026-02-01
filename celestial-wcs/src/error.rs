use thiserror::Error;

pub type WcsResult<T> = Result<T, WcsError>;

#[derive(Debug, Error)]
pub enum WcsError {
    #[error("Missing required WCS keyword: {keyword}")]
    MissingKeyword { keyword: String },

    #[error("Invalid WCS keyword '{keyword}': {message}")]
    InvalidKeyword { keyword: String, message: String },

    #[error("Unsupported projection: {code}")]
    UnsupportedProjection { code: String },

    #[error("Singularity in transformation: {message}")]
    Singularity { message: String },

    #[error("Coordinate out of bounds: {message}")]
    OutOfBounds { message: String },

    #[error("Invalid parameter: {message}")]
    InvalidParameter { message: String },

    #[error("Convergence failure: {message}")]
    ConvergenceFailure { message: String },

    #[error("Non-invertible matrix (determinant = {determinant})")]
    NonInvertibleMatrix { determinant: f64 },

    #[error("Coordinate error: {source}")]
    CoordinateError {
        #[from]
        source: celestial_coords::CoordError,
    },
}

impl WcsError {
    pub fn missing_keyword(keyword: impl Into<String>) -> Self {
        Self::MissingKeyword {
            keyword: keyword.into(),
        }
    }

    pub fn invalid_keyword(keyword: impl Into<String>, message: impl Into<String>) -> Self {
        Self::InvalidKeyword {
            keyword: keyword.into(),
            message: message.into(),
        }
    }

    pub fn unsupported_projection(code: impl Into<String>) -> Self {
        Self::UnsupportedProjection { code: code.into() }
    }

    pub fn singularity(message: impl Into<String>) -> Self {
        Self::Singularity {
            message: message.into(),
        }
    }

    pub fn out_of_bounds(message: impl Into<String>) -> Self {
        Self::OutOfBounds {
            message: message.into(),
        }
    }

    pub fn invalid_parameter(message: impl Into<String>) -> Self {
        Self::InvalidParameter {
            message: message.into(),
        }
    }

    pub fn convergence_failure(message: impl Into<String>) -> Self {
        Self::ConvergenceFailure {
            message: message.into(),
        }
    }

    pub fn non_invertible_matrix(determinant: f64) -> Self {
        Self::NonInvertibleMatrix { determinant }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_missing_keyword() {
        let err = WcsError::missing_keyword("CRPIX1");
        assert!(err.to_string().contains("CRPIX1"));
    }

    #[test]
    fn test_invalid_keyword() {
        let err = WcsError::invalid_keyword("CTYPE1", "unrecognized projection");
        assert!(err.to_string().contains("CTYPE1"));
        assert!(err.to_string().contains("unrecognized projection"));
    }

    #[test]
    fn test_unsupported_projection() {
        let err = WcsError::unsupported_projection("XYZ");
        assert!(err.to_string().contains("XYZ"));
    }

    #[test]
    fn test_singularity() {
        let err = WcsError::singularity("pole crossing");
        assert!(err.to_string().contains("pole crossing"));
    }

    #[test]
    fn test_out_of_bounds() {
        let err = WcsError::out_of_bounds("declination exceeds 90 degrees");
        assert!(err.to_string().contains("declination exceeds 90 degrees"));
    }

    #[test]
    fn test_non_invertible_matrix() {
        let err = WcsError::non_invertible_matrix(0.0);
        assert!(err.to_string().contains("0"));
    }
}
