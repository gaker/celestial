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
pub mod quad;

pub use catalog::{Catalog, CatalogHeader, StarRecord};
pub use cone::{cone_search, cone_search_at_epoch, ConeSearchParams, ConeSearchResult};
pub use quad::{
    combinations_4, discrete_hash_key, discrete_hash_key_neighbors, find_spine, get_quad_hashes,
    neighbor_quads, normalize_quad, normalize_quad_ordered, orient_spine, tan_deproject_star,
    tan_project_star, Combinations4, Quad, QuadStar,
};
