//! Types for representing nutation computation results.
//!
//! This module provides the primary types for working with nutation:
//!
//! - [`NutationResult`]: The computed nutation angles (Δψ and Δε)
//! - [`NutationModel`]: A facade for selecting and using nutation models
//!
//! Nutation values are returned in radians and represent corrections to the
//! mean pole position due to the gravitational influence of the Moon and Sun
//! on Earth's equatorial bulge.

use super::iau2000a::NutationIAU2000A;
use crate::errors::AstroResult;

/// The result of a nutation computation, containing the nutation in longitude
/// and nutation in obliquity.
///
/// Both values are expressed in radians and represent small corrections
/// (typically on the order of arcseconds) to the mean celestial pole position.
///
/// # Coordinate System
///
/// - `delta_psi` (Δψ): Nutation in longitude, measured along the ecliptic
/// - `delta_eps` (Δε): Nutation in obliquity, measured perpendicular to the ecliptic
///
/// These corrections are applied to obtain the true (apparent) pole position
/// from the mean pole position at a given epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct NutationResult {
    /// Nutation in longitude (Δψ) in radians.
    ///
    /// This is the east-west oscillation of the celestial pole along the
    /// ecliptic. Positive values indicate eastward displacement.
    pub delta_psi: f64,

    /// Nutation in obliquity (Δε) in radians.
    ///
    /// This is the north-south oscillation of the celestial pole perpendicular
    /// to the ecliptic. Positive values indicate an increase in the obliquity
    /// of the ecliptic.
    pub delta_eps: f64,
}

/// A facade for computing nutation using a selected IAU model.
///
/// This type provides a unified interface for nutation computation,
/// abstracting over the underlying model implementation. Currently
/// supports IAU 2000A, with IAU 2000A as the default.
///
/// # Usage
///
/// ```ignore
/// let model = NutationModel::iau2000a();
/// let result = model.compute(2451545.0, 0.0)?;
/// // result.delta_psi and result.delta_eps are in radians
/// ```
///
/// # Two-Part Julian Date
///
/// The `compute` method accepts a two-part Julian Date for enhanced
/// precision. The date is interpreted as `jd1 + jd2`. A common convention
/// is `jd1 = 2451545.0` (J2000.0) and `jd2 = days since J2000.0`.
#[derive(Debug, Clone)]
pub struct NutationModel {
    /// The underlying nutation calculator implementation.
    calculator: NutationIAU2000A,
}

impl NutationModel {
    /// Creates a new `NutationModel` using the IAU 2000A nutation model.
    ///
    /// IAU 2000A is the full-precision model with 1365 terms, providing
    /// sub-milliarcsecond accuracy. For applications where lower precision
    /// is acceptable, consider using IAU 2000B (77 terms, ~1 mas accuracy).
    pub fn iau2000a() -> Self {
        Self {
            calculator: NutationIAU2000A::new(),
        }
    }

    /// Computes nutation angles for the given two-part Julian Date.
    ///
    /// # Arguments
    ///
    /// * `jd1` - First part of the Julian Date (typically the integer part or J2000.0)
    /// * `jd2` - Second part of the Julian Date (typically the fractional part)
    ///
    /// # Returns
    ///
    /// A [`NutationResult`] containing `delta_psi` and `delta_eps` in radians.
    ///
    /// # Errors
    ///
    /// Returns an error if the underlying computation encounters invalid inputs
    /// or numerical issues.
    pub fn compute(&self, jd1: f64, jd2: f64) -> AstroResult<NutationResult> {
        self.calculator.compute(jd1, jd2)
    }
}

impl Default for NutationModel {
    /// Returns the default nutation model (IAU 2000A).
    fn default() -> Self {
        Self::iau2000a()
    }
}
