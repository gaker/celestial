//! Extract solver-ready metadata from FITS/XISF image headers.
//!
//! [`metadata_from_image`] reads position hint, epoch, and plate scale from an
//! [`Image`] and returns an [`ImageMetadata`] ready to feed the solver. Individual
//! pieces are also exposed: [`hint_from_image`], [`jd_from_image`],
//! [`plate_scale_from_image`].
//!
//! For images without FITS headers, construct [`ImageMetadata`] directly and pass it
//! to [`crate::solve::Solver::with_metadata`].
//!
//! [`Image`]: celestial_images::formats::image::Image
//! [`ImageMetadata`]: image_metadata::ImageMetadata
//! [`metadata_from_image`]: image_metadata::metadata_from_image
//! [`hint_from_image`]: hint::hint_from_image
//! [`jd_from_image`]: time::jd_from_image
//! [`plate_scale_from_image`]: plate_scale::plate_scale_from_image

pub mod error;
mod header_wcs;
pub mod hint;
pub mod image_metadata;
pub mod plate_scale;
pub mod time;
