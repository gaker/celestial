//! Optional observation-context headers for display and reporting.
//!
//! Unlike [`crate::metadata`], this module holds header accessors that aren't required
//! to solve an image. They feed [`crate::fit_wcs::display::SolveDisplay`] and any reporting the
//! caller wants to do alongside the solve result.

pub mod geodetic;
