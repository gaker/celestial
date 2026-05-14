//! Extract solver-ready metadata from FITS/XISF image headers.
//!
//! [`metadata_from_image`] reads position hint, epoch, and plate scale from an
//! [`Image`] and returns an [`ImageMetadata`] ready to feed the solver. Individual
//! pieces are also exposed: [`hint_from_image`], [`jd_from_image`],
//! [`plate_scale_from_image`].
//!
//! For images without FITS headers, construct [`ImageMetadata`] directly and pass it
//! to [`crate::Solver::with_metadata`].
//!
//! [`Image`]: celestial_images::formats::Image

mod error;
mod hint;
mod image_metadata;
mod plate_scale;
mod time;

pub use error::MetadataError;
pub use hint::hint_from_image;
pub use image_metadata::{metadata_from_image, ImageMetadata, MetadataOverrides};
pub use plate_scale::{plate_scale_from_image, PlateScale};
pub use time::jd_from_image;
