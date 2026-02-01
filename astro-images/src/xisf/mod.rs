pub mod errors;
pub mod header;
pub mod reader;
pub mod writer;

pub use errors::{Result, XisfError};
pub use header::{
    ColorSpace, DataLocation, ImageInfo, PixelStorage, SampleFormat, XisfCompression, XisfHeader,
    XisfProperty, XisfPropertyValue,
};
pub use reader::XisfFile;
pub use writer::{wcs_to_xisf_properties, XisfDataType, XisfWriter};
