//! Error types for astronomical calculations.
//!
//! This module provides a unified error type [`AstroError`] that covers the failure
//! modes encountered in astronomical computations: invalid dates, numerical issues,
//! data access problems, and calculation failures.
//!
//! # Error Categories
//!
//! | Variant | Use Case | Recoverable? |
//! |---------|----------|--------------|
//! | [`InvalidDate`](AstroError::InvalidDate) | Calendar validation failures | No |
//! | [`MathError`](AstroError::MathError) | Non-finite, out-of-range or invalid input, division by zero | No |
//! | [`DataError`](AstroError::DataError) | File I/O, network, parsing | Yes |
//! | [`CalculationError`](AstroError::CalculationError) | Algorithm failures | No |
//!
//! # Usage
//!
//! Most functions return [`AstroResult<T>`], which is `Result<T, AstroError>`.
//! Use the constructor methods for consistent error creation:
//!
//! ```
//! use celestial_core::errors::{AstroError, MathErrorKind};
//!
//! fn safe_divide(a: f64, b: f64) -> Result<f64, AstroError> {
//!     if b == 0.0 {
//!         return Err(AstroError::math_error(
//!             "safe_divide",
//!             MathErrorKind::DivisionByZero,
//!             "divisor is zero",
//!         ));
//!     }
//!     Ok(a / b)
//! }
//!
//! assert_eq!(safe_divide(6.0, 3.0)?, 2.0);
//! assert!(safe_divide(1.0, 0.0).is_err());
//! # Ok::<(), AstroError>(())
//! ```

use thiserror::Error;

/// Classification of mathematical errors.
///
/// Used with [`AstroError::MathError`] to distinguish between different
/// numerical failure modes.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
#[non_exhaustive]
pub enum MathErrorKind {
    /// Attempted division by zero or near-zero value.
    DivisionByZero,
    /// Input value is invalid for the operation.
    InvalidInput,
    /// Result is NaN or infinity.
    NotFinite,
    /// Value outside valid domain (e.g., latitude > 90°).
    OutOfRange,
}

/// Unified error type for astronomical calculations.
///
/// Covers calendar validation, numerical issues, data access, and algorithmic
/// failures. Use the constructor methods
/// ([`invalid_date`](Self::invalid_date), [`math_error`](Self::math_error), etc.)
/// for consistent error creation.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum AstroError {
    /// Invalid calendar date (e.g., February 30, month 13).
    #[error("Invalid date {year}-{month:02}-{day:02}: {message}")]
    InvalidDate {
        year: i32,
        month: i32,
        day: i32,
        message: String,
    },

    /// Numerical computation failure.
    #[error("Math error in {operation} ({kind}): {message}")]
    MathError {
        operation: String,
        kind: MathErrorKind,
        message: String,
    },

    /// Data access failure (file I/O, network, parsing).
    ///
    /// This is the only recoverable error variant — retry or fallback may succeed.
    #[error("Data error ({file_type} - {operation}): {message}")]
    DataError {
        file_type: String,
        operation: String,
        message: String,
    },

    /// Algorithm or calculation failure.
    #[error("Calculation error in {context}: {message}")]
    CalculationError { context: String, message: String },
}

impl core::fmt::Display for MathErrorKind {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::DivisionByZero => "division by zero",
            Self::InvalidInput => "invalid input",
            Self::NotFinite => "not finite",
            Self::OutOfRange => "out of range",
        })
    }
}

/// Convenience alias for `Result<T, AstroError>`.
pub type AstroResult<T> = Result<T, AstroError>;

impl AstroError {
    /// Creates an [`InvalidDate`](Self::InvalidDate) error.
    pub fn invalid_date(year: i32, month: i32, day: i32, reason: &str) -> Self {
        Self::InvalidDate {
            year,
            month,
            day,
            message: reason.to_string(),
        }
    }

    /// Creates a [`MathError`](Self::MathError) with the given kind.
    pub fn math_error(operation: &str, kind: MathErrorKind, reason: &str) -> Self {
        Self::MathError {
            operation: operation.to_string(),
            kind,
            message: reason.to_string(),
        }
    }

    /// Creates a [`DataError`](Self::DataError) (the only recoverable variant).
    pub fn data_error(file_type: &str, operation: &str, reason: &str) -> Self {
        Self::DataError {
            file_type: file_type.to_string(),
            operation: operation.to_string(),
            message: reason.to_string(),
        }
    }

    /// Creates a [`CalculationError`](Self::CalculationError).
    pub fn calculation_error(context: &str, reason: &str) -> Self {
        Self::CalculationError {
            context: context.to_string(),
            message: reason.to_string(),
        }
    }

    /// Returns `true` if retrying or using a fallback might succeed.
    ///
    /// Only [`DataError`](Self::DataError) is recoverable (network retry, alternate source).
    pub fn is_recoverable(&self) -> bool {
        matches!(self, Self::DataError { .. })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_invalid_date_error() {
        let err = AstroError::invalid_date(2000, 13, 1, "month out of range");
        assert_eq!(
            err.to_string(),
            "Invalid date 2000-13-01: month out of range"
        );
    }

    #[test]
    fn test_math_error_with_kind() {
        let err = AstroError::math_error(
            "vector normalization",
            MathErrorKind::DivisionByZero,
            "zero-length vector",
        );
        assert_eq!(
            err.to_string(),
            "Math error in vector normalization (division by zero): zero-length vector"
        );
    }

    #[test]
    fn test_data_error() {
        let err = AstroError::data_error("IERS Bulletin A", "download", "network timeout");
        assert!(err
            .to_string()
            .contains("Data error (IERS Bulletin A - download)"));
    }

    #[test]
    fn test_calculation_error() {
        let err = AstroError::calculation_error("orbit propagation", "insufficient data");
        assert!(err
            .to_string()
            .contains("Calculation error in orbit propagation"));
    }

    #[test]
    fn test_recoverable_errors() {
        assert!(AstroError::data_error("catalog", "download", "timeout").is_recoverable());
        assert!(!AstroError::invalid_date(2000, 13, 1, "bad month").is_recoverable());
    }

    #[test]
    fn test_send_sync() {
        fn _assert_send<T: Send>() {}
        fn _assert_sync<T: Sync>() {}
        _assert_send::<AstroError>();
        _assert_sync::<AstroError>();
    }

    #[test]
    fn test_non_recoverable_errors() {
        let math_err =
            AstroError::math_error("calculation", MathErrorKind::OutOfRange, "value too large");
        assert!(!math_err.is_recoverable());

        let calc_err = AstroError::calculation_error("orbit_propagation", "insufficient data");
        assert!(!calc_err.is_recoverable());
    }
}
