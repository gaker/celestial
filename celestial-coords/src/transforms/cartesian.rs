use celestial_core::Vector3;

/// Trait for Cartesian coordinate frame transformations.
/// Unlike `CoordinateFrame` which handles spherical sky positions,
/// this handles 3D Cartesian vectors (x, y, z).
pub trait CartesianFrame: Sized {
    /// Transform to ICRS Cartesian coordinates
    fn to_icrs(&self) -> Vector3;

    /// Create from ICRS Cartesian coordinates
    fn from_icrs(icrs: &Vector3) -> Self;
}
