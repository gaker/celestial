//! Optional observation-context headers for display and reporting.
//!
//! Unlike [`crate::metadata`], this module holds header accessors that aren't required
//! to solve an image. They feed [`crate::fit_wcs::SolveDisplay`] and any reporting the
//! caller wants to do alongside the solve result.

mod date_obs;
mod geodetic;

pub use date_obs::date_obs_from_image;
pub use geodetic::{geodetic_from_image, Geodetic};
