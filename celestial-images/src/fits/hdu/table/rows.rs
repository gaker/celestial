use super::cell::CellFormat;
use super::TableHdu;
use crate::fits::data::value::TableValue;
use crate::fits::errors::Result;
use std::io::{Read, Seek};
use std::ops::Range;

// Reads one row per step, from blocks of rows read in order. A row that
// fails to decode yields its error, and the rows after it are still read.
pub struct TableRows<'a, R, F> {
    hdu: &'a TableHdu<F>,
    reader: &'a mut R,
    next: usize,
    block: Vec<u8>,
    loaded: Range<usize>,
}

impl<'a, R: Read + Seek, F: CellFormat> TableRows<'a, R, F> {
    pub(super) fn new(hdu: &'a TableHdu<F>, reader: &'a mut R) -> Self {
        Self {
            hdu,
            reader,
            next: 0,
            block: Vec::new(),
            loaded: 0..0,
        }
    }

    fn row(&mut self, row: usize) -> Result<Vec<TableValue>> {
        if !self.loaded.contains(&row) {
            self.loaded = 0..0;
            let rows = self.hdu.row_count();
            self.loaded = self
                .hdu
                .read_block(self.reader, row, rows, &mut self.block)?;
        }
        let bytes = self.hdu.row_in_block(&self.block, self.loaded.start, row)?;
        self.hdu.decode_row(row, bytes, self.reader)
    }
}

impl<R: Read + Seek, F: CellFormat> Iterator for TableRows<'_, R, F> {
    type Item = Result<Vec<TableValue>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.hdu.row_count() {
            return None;
        }
        let row = self.row(self.next);
        self.next += 1;
        Some(row)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.hdu.row_count().saturating_sub(self.next);
        (remaining, Some(remaining))
    }
}

impl<R: Read + Seek, F: CellFormat> ExactSizeIterator for TableRows<'_, R, F> {}
