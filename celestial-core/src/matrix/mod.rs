//! 3D rotation matrices and vectors for coordinate transformations.
//!
//! - [`RotationMatrix3`]: 3×3 orthogonal matrix for frame rotations
//! - [`Vector3`]: 3D Cartesian vector

mod rotation_matrix;
mod vector3;

pub use rotation_matrix::RotationMatrix3;
pub use vector3::Vector3;
