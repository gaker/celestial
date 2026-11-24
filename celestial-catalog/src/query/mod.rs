//! Query interface for HEALPix-indexed star catalogs.
//!
//! - [`catalog`]: open a catalog file, access the header, read stars by pixel
//! - [`cone`]: cone search with magnitude filtering and proper-motion propagation
//! - [`healpix`]: coordinate-to-pixel conversion, pixel counts, RA reduction
//! - [`quad`]: TAN projection and the quad hashes the plate solver matches
//! - [`validate`]: the epoch check a catalog writer runs before it builds a file

pub mod catalog;
pub mod cone;
mod disc;
mod format;
pub mod healpix;
pub mod quad;
#[cfg(any(test, feature = "test-support"))]
pub mod test_support;
pub mod validate;
