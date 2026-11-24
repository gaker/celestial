//! Header, file length and pixel table checks run once by `Catalog::open`. `check_epoch` is
//! public so that a writer can reject an epoch before it builds a catalog file.

use celestial_core::constants::MAX_CENTURIES_FROM_J2000;

use super::catalog::CatalogHeader;
use super::format::{
    PixelEntry, CATALOG_MAGIC, CATALOG_VERSION, HEADER_SIZE, PIXEL_ENTRY_SIZE, STAR_RECORD_SIZE,
};
use super::healpix::nside_and_npix;
use crate::error::{CatalogError, CatalogResult};

pub(super) fn check_magic_and_version(magic: [u8; 4], version: u32) -> CatalogResult<()> {
    if &magic != CATALOG_MAGIC {
        return Err(CatalogError::invalid_header(format!(
            "Invalid catalog magic: expected {CATALOG_MAGIC:?}, got {magic:?}"
        )));
    }
    if version == 1 {
        return Err(CatalogError::invalid_header(
            "Catalog version 1 uses the old pixel numbering, which this library no longer reads; \
             rebuild the catalog with `forge build-index`",
        ));
    }
    if version != CATALOG_VERSION {
        return Err(CatalogError::invalid_header(format!(
            "Unsupported catalog version: expected {CATALOG_VERSION}, got {version}"
        )));
    }
    Ok(())
}

pub(super) fn check_header(header: &CatalogHeader) -> CatalogResult<()> {
    check_healpix(header)?;
    check_epoch(header.epoch)?;
    if !header.mag_limit.is_finite() {
        return Err(CatalogError::invalid_header(format!(
            "mag_limit {} is not finite",
            header.mag_limit
        )));
    }
    Ok(())
}

fn check_healpix(header: &CatalogHeader) -> CatalogResult<()> {
    let (nside, npix) = nside_and_npix(header.order).map_err(|e| match e {
        CatalogError::InvalidInput { message } => CatalogError::InvalidHeader { message },
        other => other,
    })?;
    if u64::from(header.nside) != nside {
        return Err(CatalogError::invalid_header(format!(
            "Inconsistent nside: order {} implies nside {nside}, got {}",
            header.order, header.nside
        )));
    }
    if header.npix != npix {
        return Err(CatalogError::invalid_header(format!(
            "Inconsistent npix: nside {} implies npix {npix}, got {}",
            header.nside, header.npix
        )));
    }
    Ok(())
}

// Treats an epoch more than 20 centuries from J2000, the span of the workspace's IAU models,
// as a corrupt header. Real catalog epochs are within a few decades of J2000.
pub fn check_epoch(epoch: f64) -> CatalogResult<()> {
    let span_years = MAX_CENTURIES_FROM_J2000 * 100.0;
    if !epoch.is_finite() || libm::fabs(epoch - 2000.0) > span_years {
        return Err(CatalogError::invalid_header(format!(
            "Catalog epoch J{} is not within {span_years} years of J2000",
            readable(epoch)
        )));
    }
    Ok(())
}

// Formats values outside [1e-3, 1e9) in exponent form, because the plain form of 1e300 has
// 301 digits.
fn readable(value: f64) -> String {
    let size = libm::fabs(value);
    if value == 0.0 || !value.is_finite() || (1e-3..1e9).contains(&size) {
        format!("{value}")
    } else {
        format!("{value:e}")
    }
}

// Where the sections of a checked catalog file start, and how many entries each holds.
#[derive(Debug, Clone, Copy)]
pub(super) struct Layout {
    pub(super) pixels: usize,
    pub(super) records_start: usize,
    pub(super) stars: usize,
}

impl Layout {
    // Accepts a file whose length is exactly the header, npix table entries and total_stars
    // records. The sum is formed in u128, where u64 counts times 56 cannot overflow.
    pub(super) fn new(header: &CatalogHeader, file_len: usize) -> CatalogResult<Self> {
        let table_len = u128::from(header.npix) * PIXEL_ENTRY_SIZE as u128;
        let records_len = u128::from(header.total_stars) * STAR_RECORD_SIZE as u128;
        let expected = HEADER_SIZE as u128 + table_len + records_len;
        if expected != file_len as u128 {
            return Err(CatalogError::invalid_header(format!(
                "file length {file_len} bytes does not match the {expected} bytes implied by \
                 npix {} and total_stars {}",
                header.npix, header.total_stars
            )));
        }
        // Each value is at most file_len, so it fits in usize.
        Ok(Self {
            pixels: header.npix as usize,
            records_start: HEADER_SIZE + table_len as usize,
            stars: header.total_stars as usize,
        })
    }
}

// Accepts a table whose entries cover the star records in pixel order: each entry starts
// where the previous one ended, and the last one ends at total_stars records.
pub(super) fn check_pixel_table(entries: &[PixelEntry], total_stars: u64) -> CatalogResult<()> {
    let records_len = total_stars
        .checked_mul(STAR_RECORD_SIZE as u64)
        .ok_or_else(|| {
            CatalogError::invalid_header(format!("total_stars {total_stars} records overflow u64"))
        })?;
    let mut next_offset = 0u64;
    for (pixel, entry) in entries.iter().enumerate() {
        next_offset = check_entry(pixel, entry, next_offset, records_len)?;
    }
    if next_offset != records_len {
        return Err(CatalogError::corrupt_table(format!(
            "pixel counts sum to {} records, but the header's total_stars is {total_stars}",
            next_offset / STAR_RECORD_SIZE as u64
        )));
    }
    Ok(())
}

// Returns the offset where the next pixel's records must start.
fn check_entry(
    pixel: usize,
    entry: &PixelEntry,
    offset: u64,
    records_len: u64,
) -> CatalogResult<u64> {
    let corrupt =
        |problem: String| CatalogError::corrupt_table(format!("pixel {pixel}: {problem}"));
    let (record_size, found) = (STAR_RECORD_SIZE as u64, entry.offset);
    if !found.is_multiple_of(record_size) {
        let problem = format!("offset {found} is not a multiple of {record_size}, the record size");
        return Err(corrupt(problem));
    }
    if found != offset {
        let problem = format!("offset {found} should be {offset}, where earlier pixels end");
        return Err(corrupt(problem));
    }
    let end = u64::from(entry.count)
        .checked_mul(record_size)
        .and_then(|len| len.checked_add(offset));
    end.filter(|end| *end <= records_len).ok_or_else(|| {
        let count = entry.count;
        corrupt(format!(
            "{count} records at offset {offset} run past the {records_len}-byte star data"
        ))
    })
}
