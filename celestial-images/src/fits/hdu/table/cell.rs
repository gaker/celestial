use super::layout::Column;
use crate::fits::data::value::TableValue;
use crate::fits::errors::{FitsError, Result};
use crate::fits::hdu::HduType;
use crate::fits::header::parser::Header;
use std::io::{Read, Seek, SeekFrom};
use std::ops::Range;

// The cell format of one table type: text fields for TABLE, big-endian
// binary for BINTABLE. This module is private to the crate, so no
// other crate can implement it.
pub trait CellFormat: Sized + std::fmt::Debug {
    // FITS 4.0 §7.3.1 requires binary columns to fill NAXIS1 exactly, in
    // order. §7.2.4 lets ASCII fields sit anywhere in the row.
    const PACKED: bool;

    fn hdu_type() -> HduType;

    // Parses TFORMn of column `n` (1-based) and returns the column's byte
    // range in the row. `next` is the first byte after the previous column.
    fn place(tform: &str, header: &Header, n: usize, next: usize) -> Result<(Self, Range<usize>)>;

    fn decode<R: Read + Seek>(
        column: &Column<Self>,
        cell: &[u8],
        heap: &mut Heap<'_, R>,
    ) -> Result<TableValue>;
}

// The heap of a binary table (FITS 4.0 §7.3.5). ASCII tables have none.
pub struct Heap<'a, R> {
    reader: &'a mut R,
    start: u64,
    len: usize,
}

impl<'a, R> Heap<'a, R> {
    pub(super) fn new(reader: &'a mut R, start: u64, len: usize) -> Self {
        Self { reader, start, len }
    }
}

impl<R: Read + Seek> Heap<'_, R> {
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    // `range` must lie inside the heap, which the caller checks against
    // `len` before calling, so the allocation is bounded by PCOUNT.
    pub(crate) fn read(&mut self, range: Range<usize>) -> Result<Vec<u8>> {
        if range.start > range.end || range.end > self.len {
            return Err(FitsError::InvalidFormat(format!(
                "heap bytes {range:?} lie outside the {}-byte heap",
                self.len
            )));
        }
        self.reader
            .seek(SeekFrom::Start(self.start + range.start as u64))?;
        let mut bytes = vec![0u8; range.len()];
        self.reader.read_exact(&mut bytes)?;
        Ok(bytes)
    }
}
