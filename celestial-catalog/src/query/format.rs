//! On-disk layout of a catalog file: the magic and version, the section sizes, the pixel table
//! entry and the little-endian byte helpers that read and write header and record fields.

use super::catalog::StarRecord;

// Catalog::stars_in_pixel and Catalog::pixel_table reinterpret the little-endian file bytes as
// StarRecord and PixelEntry values in place.
#[cfg(target_endian = "big")]
compile_error!("celestial-catalog reads its little-endian catalog records in place");

pub(super) const CATALOG_MAGIC: &[u8; 4] = b"CCAT";
// Version 1 files were built with a HEALPix function that filed stars from RA 270 to 360 in
// the equatorial belt under the wrong base face. Version 2 uses the standard NEST numbering.
pub(super) const CATALOG_VERSION: u32 = 2;
pub(super) const HEADER_SIZE: usize = 64;
pub(super) const PIXEL_ENTRY_SIZE: usize = 16;
pub(super) const STAR_RECORD_SIZE: usize = 56;

// One entry of the pixel table: the byte offset of the pixel's first record from the start
// of the star data, and the number of records.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub(super) struct PixelEntry {
    pub(super) offset: u64,
    pub(super) count: u32,
    _reserved: u32,
}

const _: () = assert!(std::mem::size_of::<PixelEntry>() == PIXEL_ENTRY_SIZE);
const _: () = assert!(std::mem::size_of::<StarRecord>() == STAR_RECORD_SIZE);

// Every section starts on an 8-byte boundary of the file, and records follow each other
// without gaps, so a map that starts 8-byte aligned holds aligned records.
const _: () = assert!(
    std::mem::align_of::<StarRecord>() == 8
        && std::mem::align_of::<PixelEntry>() == 8
        && HEADER_SIZE.is_multiple_of(8)
        && PIXEL_ENTRY_SIZE.is_multiple_of(8)
        && STAR_RECORD_SIZE.is_multiple_of(8)
);

// Concatenates the fields into an N-byte array and zero-fills the rest. The callers pass
// fixed-size fields whose total length is at most N.
pub(super) fn concat_le<const N: usize>(fields: &[&[u8]]) -> [u8; N] {
    let mut out = [0u8; N];
    for (slot, byte) in out.iter_mut().zip(fields.iter().copied().flatten()) {
        *slot = *byte;
    }
    out
}

// Reads the N bytes at `at`. The callers read fixed offsets inside arrays of known length.
pub(super) fn le_field<const N: usize>(bytes: &[u8], at: usize) -> [u8; N] {
    std::array::from_fn(|i| bytes[at + i])
}
