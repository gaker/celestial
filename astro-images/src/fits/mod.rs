pub mod compression;
pub mod data;
pub mod errors;
pub mod hdu;
pub mod header;
pub mod image;
pub mod io;
pub mod util;
pub mod wcs;

pub use compression::{CompressionAlgorithm, CompressionParams};
pub use data::array::TableValue;
pub use errors::{FitsError, Result};
pub use hdu::{
    AsciiTableHdu, AsciiTableRowIterator, BinaryTableHdu, BinaryTableRowIterator, Hdu, ImageHdu,
    PrimaryHdu,
};
pub use image::{FitsImage, ImageKind};
pub use io::{FitsFile, FitsReader, FitsWriter};
pub use wcs::WcsInfo;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fits_module_exports() {
        fn _type_checks() {
            fn _fits_error(_: FitsError) {}
            fn _result(_: Result<()>) {}
        }
    }
}
