use crate::errors::{CoordError, CoordResult};
use celestial_core::{angle::Angle, matrix::Vector3};

// Longitude and latitude of a direction of any non-zero length, without
// normalising first. The longitude is in (-π, π]; constructors wrap it.
pub(super) fn spherical_angles(v: Vector3) -> CoordResult<(Angle, Angle)> {
    if !(v.x.is_finite() && v.y.is_finite() && v.z.is_finite()) {
        return Err(CoordError::invalid_coordinate(
            "Direction vector is not finite",
        ));
    }
    if v == Vector3::zeros() {
        return Err(CoordError::invalid_coordinate("Zero vector"));
    }
    let (lon, lat) = v.to_spherical();
    Ok((Angle::from_radians(lon), Angle::from_radians(lat)))
}
