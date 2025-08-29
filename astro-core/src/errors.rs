//! Shared error types for astronomical calculation libraries.
//!
//! This module provides standardized error handling across the astro-core workspace.
//! Each error type captures specific context that helps with debugging and error
//! recovery in astronomical calculations.
//!
//! # Quick Start
//!
//! ```rust
//! use astro_core::{AstroError, AstroResult, MathErrorKind};
//!
//! fn calculate_julian_date(year: i32, month: i32, day: i32) -> AstroResult<f64> {
//!     if month < 1 || month > 12 {
//!         return Err(AstroError::invalid_date(year, month, day, "month out of range"));
//!     }
//!     
//!     // Perform calculation...
//!     Ok(2451545.0)
//! }
//! ```
//!
//! # Error Types
//!
//! - [`AstroError::InvalidDate`] - Calendar date validation failures
//! - [`AstroError::MathError`] - Numerical computation problems  
//! - [`AstroError::ExternalLibraryError`] - C library function failures
//! - [`AstroError::DataError`] - External data file issues
//! - [`AstroError::CalculationError`] - General calculation failures

use thiserror::Error;

/// Types of mathematical errors that can occur during numerical computations.
///
/// Used with [`AstroError::MathError`] to provide specific context about
/// what went wrong during a mathematical operation.
#[derive(Debug, Clone, PartialEq)]
pub enum MathErrorKind {
    /// Numeric overflow (result too large to represent)
    Overflow,
    /// Numeric underflow (result too small to represent)
    Underflow,
    /// Significant precision loss detected
    PrecisionLoss,
    /// Division by zero attempted
    DivisionByZero,
    /// Input value outside valid domain
    InvalidInput,
}

/// Standard error types for astronomical calculations.
///
/// Provides consistent error handling across all astro-core crates. Each variant
/// captures specific context to help with debugging and recovery.
#[derive(Error, Debug)]
pub enum AstroError {
    /// Calendar date validation failed.
    ///
    /// Use when date components are invalid (month 13, day 32, etc.).
    /// Common in time conversions, coordinate calculations, and ephemeris lookups.
    ///
    /// # Recovery
    ///
    /// Not recoverable - fix the input data.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::AstroError;
    /// let err = AstroError::invalid_date(2023, 13, 1, "month out of range");
    /// println!("{}", err); // "Invalid date 2023-13-01: month out of range"
    /// ```
    #[error("Invalid date {year}-{month:02}-{day:02}: {message}")]
    InvalidDate {
        /// Year component
        year: i32,
        /// Month component (1-12)
        month: i32,
        /// Day component (1-31)
        day: i32,
        /// Description of what's wrong
        message: String,
    },

    /// Numerical computation failed.
    ///
    /// Use for overflow, underflow, precision loss, division by zero, or invalid
    /// mathematical inputs. Includes the operation name and specific error kind.
    ///
    /// # Recovery
    ///
    /// Usually not recoverable - indicates fundamental input or algorithm issues.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::{AstroError, MathErrorKind};
    /// let err = AstroError::math_error(
    ///     "julian_date_conversion",
    ///     MathErrorKind::Overflow,
    ///     "year too large"
    /// );
    /// ```
    #[error("Math error in {operation} ({kind:?}): {message}")]
    MathError {
        /// Name of the operation that failed
        operation: String,
        /// Type of mathematical error
        kind: MathErrorKind,
        /// Detailed error description
        message: String,
    },

    /// External C library function returned an error.
    ///
    /// Use when calling ERFA, SOFA, or other C libraries. The status code comes
    /// from the library's return value.
    ///
    /// # Recovery
    ///
    /// Usually not recoverable - indicates bad inputs or library configuration.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::AstroError;
    /// let err = AstroError::external_library_error(
    ///     "eraCal2jd",
    ///     -1,
    ///     "bad calendar date"
    /// );
    /// ```
    #[error("External library error in {function}: status {status_code} - {message}")]
    ExternalLibraryError {
        /// Name of the C function that failed
        function: String,
        /// Library's error status code
        status_code: i32,
        /// Error message from library or wrapper
        message: String,
    },

    /// External data file operation failed.
    ///
    /// Use for downloads, parsing, or reading of astronomical data files like
    /// IERS bulletins, star catalogs, or ephemeris files.
    ///
    /// # Recovery
    ///
    /// Often recoverable - can retry download or fall back to cached data.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::AstroError;
    /// let err = AstroError::data_error(
    ///     "IERS Bulletin A",
    ///     "download",
    ///     "network timeout"
    /// );
    /// assert!(err.is_recoverable());
    /// ```
    #[error("Data error ({file_type} - {operation}): {message}")]
    DataError {
        /// Type of data file (e.g., "IERS Bulletin A", "Hipparcos catalog")
        file_type: String,
        /// Operation that failed (e.g., "download", "parse", "interpolate")
        operation: String,
        /// Specific error description
        message: String,
    },

    /// General calculation failed.
    ///
    /// Fallback for errors that don't fit other categories. Prefer specific
    /// error types when possible.
    ///
    /// # Recovery
    ///
    /// Usually not recoverable without additional context.
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::AstroError;
    /// let err = AstroError::calculation_error(
    ///     "orbit_propagation",
    ///     "insufficient observational data"
    /// );
    /// ```
    #[error("Calculation error in {context}: {message}")]
    CalculationError {
        /// Context where the error occurred
        context: String,
        /// Description of the failure
        message: String,
    },
}

// thiserror handles Display and Error trait implementations

/// Result type alias for astronomical operations.
///
/// Equivalent to `Result<T, AstroError>`. Use this as the return type for
/// functions that can fail with astronomical errors.
///
/// # Example
///
/// ```rust
/// # use astro_core::AstroResult;
/// fn parse_coordinates(ra_str: &str, dec_str: &str) -> AstroResult<(f64, f64)> {
///     // Parse and validate coordinates...
///     Ok((15.25, -22.5))
/// }
/// ```
pub type AstroResult<T> = Result<T, AstroError>;

/// Convenience constructors for common error patterns.
impl AstroError {
    /// Create an invalid date error.
    ///
    /// # Arguments
    ///
    /// * `year` - Year component
    /// * `month` - Month component (1-12)
    /// * `day` - Day component (1-31)
    /// * `reason` - Why the date is invalid
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::AstroError;
    /// let err = AstroError::invalid_date(2023, 13, 1, "month out of range");
    /// assert_eq!(err.to_string(), "Invalid date 2023-13-01: month out of range");
    /// ```
    pub fn invalid_date(year: i32, month: i32, day: i32, reason: &str) -> Self {
        Self::InvalidDate {
            year,
            month,
            day,
            message: reason.to_string(),
        }
    }

    /// Create a mathematical computation error.
    ///
    /// # Arguments
    ///
    /// * `operation` - Name of the operation that failed
    /// * `kind` - Type of mathematical error
    /// * `reason` - Specific description of the failure
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::{AstroError, MathErrorKind};
    /// let err = AstroError::math_error(
    ///     "sqrt_calculation",
    ///     MathErrorKind::InvalidInput,
    ///     "negative input"
    /// );
    /// ```
    pub fn math_error(operation: &str, kind: MathErrorKind, reason: &str) -> Self {
        Self::MathError {
            operation: operation.to_string(),
            kind,
            message: reason.to_string(),
        }
    }

    /// Create an external library function error.
    ///
    /// Use when C library functions return error status codes.
    ///
    /// # Arguments
    ///
    /// * `function` - Name of the C function that failed
    /// * `status_code` - Error code returned by the library
    /// * `message` - Error description (from library docs or wrapper)
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::AstroError;
    /// let err = AstroError::external_library_error(
    ///     "eraCal2jd",
    ///     -1,
    ///     "invalid calendar date"
    /// );
    /// ```
    pub fn external_library_error(function: &str, status_code: i32, message: &str) -> Self {
        Self::ExternalLibraryError {
            function: function.to_string(),
            status_code,
            message: message.to_string(),
        }
    }

    /// Create a data file operation error.
    ///
    /// Use for external data file failures (downloads, parsing, etc.).
    ///
    /// # Arguments
    ///
    /// * `file_type` - Type of data file (e.g., "IERS Bulletin A")
    /// * `operation` - Operation that failed (e.g., "download", "parse")
    /// * `reason` - Specific failure description
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::AstroError;
    /// let err = AstroError::data_error(
    ///     "star catalog",
    ///     "download",
    ///     "HTTP 404 - file not found"
    /// );
    /// assert!(err.is_recoverable());
    /// ```
    pub fn data_error(file_type: &str, operation: &str, reason: &str) -> Self {
        Self::DataError {
            file_type: file_type.to_string(),
            operation: operation.to_string(),
            message: reason.to_string(),
        }
    }

    /// Create a general calculation error.
    ///
    /// Fallback for errors that don't fit specific categories.
    /// Prefer more specific error types when possible.
    ///
    /// # Arguments
    ///
    /// * `context` - Where the error occurred
    /// * `reason` - What went wrong
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::AstroError;
    /// let err = AstroError::calculation_error(
    ///     "asteroid_orbit_fitting",
    ///     "insufficient observations"
    /// );
    /// ```
    pub fn calculation_error(context: &str, reason: &str) -> Self {
        Self::CalculationError {
            context: context.to_string(),
            message: reason.to_string(),
        }
    }

    /// Check if this error might be recoverable with retry logic.
    ///
    /// Returns `true` for errors that might succeed on retry (network timeouts,
    /// temporary file locks, etc.). Returns `false` for errors that indicate
    /// fundamental problems with input data or algorithm logic.
    ///
    /// # Recovery Strategies
    ///
    /// - `DataError`: Retry download, use cached data, or try alternative sources
    /// - All others: Fix input data or code logic - retrying won't help
    ///
    /// # Example
    ///
    /// ```rust
    /// # use astro_core::AstroError;
    /// let network_error = AstroError::data_error("IERS data", "download", "timeout");
    /// let date_error = AstroError::invalid_date(2023, 13, 1, "invalid month");
    ///
    /// assert!(network_error.is_recoverable());  // Can retry download
    /// assert!(!date_error.is_recoverable());    // Month 13 will never be valid
    /// ```
    pub fn is_recoverable(&self) -> bool {
        match self {
            Self::DataError { .. } => true,    // Could retry download/parsing
            Self::InvalidDate { .. } => false, // Can't fix bad input
            _ => false,
        }
    }
}

// Note: No automatic From implementations since errors should be explicit

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
            "nanosecond addition",
            MathErrorKind::Overflow,
            "value too large",
        );
        assert!(err.to_string().contains("Math error"));
        assert!(err.to_string().contains("Overflow"));
    }

    #[test]
    fn test_external_library_error() {
        let err = AstroError::external_library_error("libraryFunction", -2, "bad input");
        assert!(err.to_string().contains("bad input"));
        assert!(err.to_string().contains("libraryFunction"));
        assert!(err.to_string().contains("status -2"));
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
        // Compile-time check that AstroError implements Send + Sync
        fn _assert_send<T: Send>() {}
        fn _assert_sync<T: Sync>() {}
        _assert_send::<AstroError>();
        _assert_sync::<AstroError>();
    }

    #[test]
    fn test_non_recoverable_errors() {
        // Test that MathError is not recoverable
        let math_err = AstroError::math_error(
            "calculation",
            MathErrorKind::Overflow,
            "value too large"
        );
        assert!(!math_err.is_recoverable());

        // Test that ExternalLibraryError is not recoverable
        let lib_err = AstroError::external_library_error(
            "erfa_function",
            -1,
            "bad input"
        );
        assert!(!lib_err.is_recoverable());

        // Test that CalculationError is not recoverable
        let calc_err = AstroError::calculation_error(
            "orbit_propagation",
            "insufficient data"
        );
        assert!(!calc_err.is_recoverable());
    }
}