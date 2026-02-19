//! HEALPix-indexed star catalog for fast positional queries.
//!
//! Provides memory-mapped access to a binary catalog of ~37 million stars
//! (Gaia DR3 + Hipparcos) organized as a HEALPix spatial index. The catalog
//! file is memory-mapped on open — no parsing, no loading into RAM. Star
//! records are read as zero-copy `repr(C)` slices directly from the map.
//!
//! # Modules
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`query::catalog`] | [`Catalog`](query::Catalog) reader, [`StarRecord`](query::StarRecord), [`CatalogHeader`](query::CatalogHeader), flag constants |
//! | [`query::cone`] | [`cone_search`](query::cone_search), [`ConeSearchParams`](query::ConeSearchParams), proper-motion propagation |
//! | [`query::healpix`] | Pixel indexing ([`ang2pix_nest`](query::ang2pix_nest)), disc queries, angular separation |
//!
//! # Quick Start
//!
//! ```ignore
//! use celestial_catalog::query::{Catalog, cone_search, ConeSearchParams};
//!
//! let catalog = Catalog::open("catalog.bin")?;
//!
//! let results = cone_search(&catalog, &ConeSearchParams {
//!     ra_deg: 83.633,
//!     dec_deg: -5.375,
//!     radius_deg: 0.5,
//!     max_mag: Some(14.0),
//!     max_results: Some(50),
//!     epoch: None,
//! });
//! ```
//!
//! # Binary Format
//!
//! The catalog file has three sections: a 64-byte header, a pixel offset table
//! (`npix × 16` bytes), and contiguous star data (`total_stars × 56` bytes).
//! Stars are grouped by HEALPix pixel, enabling spatial queries that touch
//! only the relevant pages.
//!
//! # Features
//!
//! - **`cli`** — Enables the `forge` and `query-catalog` binaries for building
//!   and querying catalogs from the command line.

pub mod query;
