use super::{validate_fits_block_alignment, FitsFile, HduInfo};
use crate::fits::block::{align_to_block, CARD_SIZE, FITS_BLOCK_SIZE};
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::parser::Header;
use std::io::{Read, Seek, SeekFrom};

impl<R: Read + Seek> FitsFile<R> {
    // The scan stops at the first HDU that does not parse and keeps the HDUs
    // before it, as astropy does. `scan_error` returns the error. A
    // malformed primary HDU fails the open. EXTEND is only advisory (FITS 4.0
    // §4.4.2.1), so it does not end the scan.
    pub(super) fn scan_hdus(&mut self) -> Result<()> {
        let file_len = self.reader.seek(SeekFrom::End(0))?;
        let mut position = 0u64;

        loop {
            let hdu_info = match self.scan_single_hdu(position, self.hdus.len(), file_len) {
                Ok(hdu_info) => hdu_info,
                Err(e) if self.hdus.is_empty() => return Err(e),
                Err(e) => {
                    self.scan_error = Some(e);
                    return Ok(());
                }
            };
            position = Self::calculate_next_position(&hdu_info);
            self.hdus.push(hdu_info);

            if self.check_end_of_file(position)? {
                return Ok(());
            }
        }
    }

    pub(super) fn scan_single_hdu(
        &mut self,
        position: u64,
        hdu_index: usize,
        file_len: u64,
    ) -> Result<HduInfo> {
        let header_size = self.determine_header_size(position)?;
        validate_fits_block_alignment(header_size as u64, "header")?;
        let header = self.read_header_at(position, header_size)?;
        check_first_card(&header, hdu_index)?;
        check_structural_keywords_unique(&header, hdu_index)?;

        let data_start = position + header_size as u64;
        validate_fits_block_alignment(data_start, "data start position")?;
        let available = file_len.saturating_sub(data_start);

        Ok(HduInfo {
            index: hdu_index,
            header_start: position,
            header_size,
            data_start,
            data_size: self.data_size_within(&header, hdu_index, available)?,
        })
    }

    pub(super) fn calculate_next_position(hdu_info: &HduInfo) -> u64 {
        let position = hdu_info.data_start + hdu_info.data_size as u64;
        align_to_block(position)
    }

    // The scan ends when no bytes remain, or when every remaining byte is
    // NUL. astropy reads trailing NUL bytes as extra padding. Any other
    // trailing bytes must parse as an extension.
    pub(super) fn check_end_of_file(&mut self, position: u64) -> Result<bool> {
        let file_len = self.reader.seek(SeekFrom::End(0))?;
        if position >= file_len {
            return Ok(true);
        }
        self.reader.seek(SeekFrom::Start(position))?;

        let mut block = [0u8; FITS_BLOCK_SIZE];
        loop {
            let read = match self.reader.read(&mut block) {
                Ok(read) => read,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e.into()),
            };
            if read == 0 {
                return Ok(true);
            }
            if block[..read].iter().any(|&b| b != 0) {
                return Ok(false);
            }
        }
    }

    pub(super) fn determine_header_size(&mut self, start_position: u64) -> Result<usize> {
        let current_pos = self.reader.stream_position()?;
        self.reader.seek(SeekFrom::Start(start_position))?;
        let size = self.count_header_bytes(start_position);
        // The position is restored after a failed scan too, so the caller
        // reads from where it was. The scan error takes precedence.
        let restored = self.reader.seek(SeekFrom::Start(current_pos));
        let size = size?;
        restored?;
        Ok(size)
    }

    // Reads whole blocks from the current position up to the one holding END.
    fn count_header_bytes(&mut self, start_position: u64) -> Result<usize> {
        const MAX_HEADER_BLOCKS: usize = 1000;

        let mut block_buffer = vec![0u8; FITS_BLOCK_SIZE];
        for blocks_read in 0..MAX_HEADER_BLOCKS {
            self.reader.read_exact(&mut block_buffer).map_err(|e| {
                FitsError::InvalidFormat(format!(
                    "Header at byte {} has no END card in its {} complete blocks: {}",
                    start_position, blocks_read, e
                ))
            })?;
            let mut cards = block_buffer.as_chunks::<CARD_SIZE>().0.iter();
            if cards.any(|card| card[..8].trim_ascii() == b"END") {
                return Ok((blocks_read + 1) * FITS_BLOCK_SIZE);
            }
        }
        Err(FitsError::InvalidFormat(format!(
            "Header exceeds maximum size of {} blocks ({} bytes)",
            MAX_HEADER_BLOCKS,
            MAX_HEADER_BLOCKS * FITS_BLOCK_SIZE
        )))
    }
}

// FITS 4.0 §4.4.1: the primary header starts with SIMPLE and an extension
// header with XTENSION. astropy and fitsverify reject other files. The order
// of the other mandatory keywords is not checked, as astropy accepts any.
fn check_first_card(header: &Header, hdu_index: usize) -> Result<()> {
    let expected = if hdu_index == 0 { "SIMPLE" } else { "XTENSION" };
    match header.keywords().first() {
        Some(keyword) if keyword.name == expected => Ok(()),
        _ => Err(FitsError::InvalidFormat(format!(
            "HDU {} header does not start with {}",
            hdu_index, expected
        ))),
    }
}

// A mandatory keyword of FITS 4.0 §4.4.1 given twice leaves the HDU layout
// ambiguous, for example when a damaged END card merges two headers.
// fitsverify reports it as an error, so the header is malformed even when
// the values agree. Other keywords may repeat. Lookups return the first
// copy, and data and WCS readers apply the last, as astropy and cfitsio do.
fn check_structural_keywords_unique(header: &Header, hdu_index: usize) -> Result<()> {
    let mut seen = std::collections::HashSet::new();
    for keyword in header.keywords() {
        let name = keyword.name.as_str();
        if is_structural(name) && !seen.insert(name) {
            return Err(FitsError::InvalidFormat(format!(
                "HDU {} header repeats the keyword {}",
                hdu_index, name
            )));
        }
    }
    Ok(())
}

fn is_structural(name: &str) -> bool {
    const EXACT: [&str; 7] = [
        "SIMPLE", "XTENSION", "BITPIX", "NAXIS", "PCOUNT", "GCOUNT", "TFIELDS",
    ];
    const INDEXED: [&str; 2] = ["NAXIS", "TFORM"];
    let indexed = |prefix: &str| {
        name.strip_prefix(prefix)
            .is_some_and(|n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()))
    };
    EXACT.contains(&name) || INDEXED.iter().any(|p| indexed(p))
}
