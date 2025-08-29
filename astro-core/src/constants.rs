
//! Physical and mathematical constants for astronomical calculations
//!
//! These constants are derived from internationally accepted standards and are
//! suitable for precision astronomical computations.

// WGS84 Earth ellipsoid parameters (EPSG:4326)
// Reference: NIMA TR8350.2, "Department of Defense World Geodetic System 1984"

/// Semi-major axis of WGS84 ellipsoid in kilometers
/// 
/// Value: 6378.137 km exactly (defining parameter)
pub const WGS84_A: f64 = 6378.137;

/// First eccentricity squared for WGS84 ellipsoid
/// 
/// Pre-computed from flattening f = 1/298.257223563:
/// e² = 2f - f² = 6.694379990141316×10⁻³
/// 
/// Used in geodetic coordinate transformations. Precision maintained
/// to 15 decimal places for nanosecond-accuracy time corrections.
pub const WGS84_E2: f64 = 6.694379990141316e-3;

// Time unit conversions

/// Nanoseconds in one second (10⁹)
/// 
/// Used for high-precision time arithmetic where nanosecond resolution
/// is required (e.g., pulsar timing, precision astrometry).
pub const NANOSECONDS_PER_SECOND: u32 = 1_000_000_000;

/// Seconds in one mean solar day
/// 
/// SI definition: exactly 86,400 seconds per day.
/// Does not account for leap seconds in UTC.
pub const SECONDS_PER_DAY: u32 = 86_400;

/// Seconds in one mean solar day (floating point)
/// 
/// Used in time scale conversions requiring fractional day arithmetic.
/// Maintains full precision for sub-second calculations.
pub const SECONDS_PER_DAY_F64: f64 = 86400.0;

/// Hours in one day
pub const HOURS_PER_DAY: f64 = 24.0;

/// Minutes in one day  
pub const MINUTES_PER_DAY: f64 = 1440.0;