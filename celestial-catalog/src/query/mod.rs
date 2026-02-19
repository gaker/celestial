//! Query interface for HEALPix-indexed star catalogs.
//!
//! Three submodules cover the full query surface:
//!
//! - [`catalog`] — open a catalog file, access the header, read stars by pixel
//! - [`cone`] — cone search with magnitude filtering and proper-motion propagation
//! - [`healpix`] — coordinate-to-pixel conversion, disc queries, angular separation

pub mod catalog;
pub mod cone;
pub mod healpix;

pub use catalog::{Catalog, CatalogHeader, StarRecord};
pub use cone::{cone_search, cone_search_at_epoch, ConeSearchParams, ConeSearchResult};
