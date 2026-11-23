use super::FitsFile;
use crate::fits::block::align_to_block;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::parser::Header;
use std::io::{Read, Seek};

impl<R: Read + Seek> FitsFile<R> {
    // The declared size is checked against the file before anything is
    // allocated. Only the zero padding after the data may be missing,
    // which astropy also accepts.
    pub(super) fn data_size_within(
        &self,
        header: &Header,
        hdu_index: usize,
        available: u64,
    ) -> Result<usize> {
        let data_bytes = data_unit_bytes(header)?;
        if data_bytes > available {
            return Err(FitsError::InvalidFormat(format!(
                "HDU {} declares {} data bytes, but the file holds {} after its header",
                hdu_index, data_bytes, available
            )));
        }
        self.calculate_data_size(header)
    }

    pub(super) fn calculate_data_size(&self, header: &Header) -> Result<usize> {
        let size = align_to_block(data_unit_bytes(header)?);
        usize::try_from(size).map_err(|_| too_large())
    }
}

// The data unit before padding, as a length in memory.
pub(crate) fn data_unit_len(header: &Header) -> Result<usize> {
    usize::try_from(data_unit_bytes(header)?).map_err(|_| too_large())
}

// FITS 4.0 §4.4.1.1 and eq. 2: |BITPIX|/8 × GCOUNT × (PCOUNT + NAXIS1 ×
// ... × NAXISn) bytes, before padding. Random groups (§6) set NAXIS1 = 0
// and leave it out of the product. Missing PCOUNT and GCOUNT read as 0
// and 1, as astropy reads them.
fn data_unit_bytes(header: &Header) -> Result<u64> {
    let bitpix = required_integer(header, "BITPIX")?;
    let bitpix = i32::try_from(bitpix).map_err(|_| invalid_value("BITPIX", bitpix))?;
    let bytes_per_value = crate::core::types::BitPix::from_value(bitpix)
        .ok_or(FitsError::InvalidBitPix(bitpix))?
        .bytes_per_pixel() as u64;

    let axes = axis_lengths(header)?;
    let random_groups = header.get("GROUPS").as_bool() == Some(true) && axes.first() == Some(&0);
    let array_values = match axes.split_first() {
        None => Some(0),
        Some((_, rest)) if random_groups => rest.iter().try_fold(1u64, |n, &a| n.checked_mul(a)),
        Some(_) => axes.iter().try_fold(1u64, |n, &a| n.checked_mul(a)),
    };
    let pcount = optional_count(header, "PCOUNT", 0)?;
    let gcount = optional_count(header, "GCOUNT", 1)?;

    array_values
        .and_then(|n| n.checked_add(pcount))
        .and_then(|n| n.checked_mul(gcount))
        .and_then(|n| n.checked_mul(bytes_per_value))
        .ok_or_else(too_large)
}

fn too_large() -> FitsError {
    FitsError::InvalidFormat("Data dimensions too large".to_string())
}

fn invalid_value(keyword: &str, value: impl std::fmt::Display) -> FitsError {
    FitsError::InvalidKeywordValue {
        keyword: keyword.to_string(),
        value: value.to_string(),
    }
}

fn required_integer(header: &Header, name: &str) -> Result<i64> {
    let keyword = header
        .get_keyword(name)
        .ok_or_else(|| FitsError::KeywordNotFound {
            keyword: name.to_string(),
        })?;
    let value = keyword
        .value
        .as_ref()
        .ok_or_else(|| invalid_value(name, "no value"))?;
    value.as_integer().ok_or_else(|| invalid_value(name, value))
}

// FITS 4.0 §4.4.1.1: 0 <= NAXIS <= 999, and NAXISn is a non-negative integer
// for every n up to NAXIS.
fn axis_lengths(header: &Header) -> Result<Vec<u64>> {
    let naxis = required_integer(header, "NAXIS")?;
    if !(0..=999).contains(&naxis) {
        return Err(invalid_value("NAXIS", naxis));
    }
    (1..=naxis)
        .map(|n| {
            let name = format!("NAXIS{}", n);
            let length = required_integer(header, &name)?;
            u64::try_from(length).map_err(|_| invalid_value(&name, length))
        })
        .collect()
}

fn optional_count(header: &Header, name: &str, default: u64) -> Result<u64> {
    if header.get_keyword(name).is_none() {
        return Ok(default);
    }
    let count = required_integer(header, name)?;
    u64::try_from(count).map_err(|_| invalid_value(name, count))
}
