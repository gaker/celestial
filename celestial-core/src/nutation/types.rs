//! Types for representing nutation computation results.
//!
//! [`NutationResult`] holds the computed nutation angles (Δψ and Δε).
//!
//! Nutation values are returned in radians and represent corrections to the
//! mean pole position due to the gravitational influence of the Moon and Sun
//! on Earth's equatorial bulge.

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
#[derive(Debug, Clone, Copy, PartialEq)]
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
