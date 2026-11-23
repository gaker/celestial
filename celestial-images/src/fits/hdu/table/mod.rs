pub(crate) mod cell;
pub(crate) mod layout;
pub mod rows;
pub(crate) mod values;

use super::{ColumnInfo, HduTrait, HduType};
use crate::fits::data::value::TableValue;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::parser::Header;
use crate::fits::io::reader::HduInfo;
use cell::{CellFormat, Heap};
use layout::Layout;
use rows::TableRows;
use std::io::{Read, Seek, SeekFrom};
use std::ops::Range;

// Rows are read in blocks of about this many bytes.
const BLOCK_BYTES: usize = 1 << 16;

#[derive(Debug)]
pub struct TableHdu<F> {
    header: Header,
    info: HduInfo,
    layout: Layout<F>,
}

impl<F: CellFormat> TableHdu<F> {
    // The column layout is checked here, so a table whose columns do not
    // fit its rows is rejected before any cell is read.
    pub(crate) fn new(header: Header, info: HduInfo) -> Result<Self> {
        let layout = Layout::new(&header, info.data_size)?;
        Ok(Self {
            header,
            info,
            layout,
        })
    }

    pub fn number_of_fields(&self) -> Option<i64> {
        self.header.get("TFIELDS").as_i64()
    }

    pub fn number_of_rows(&self) -> Option<i64> {
        self.header.get("NAXIS2").as_i64()
    }

    pub fn extension_name(&self) -> Option<&str> {
        self.header.get("EXTNAME").as_str()
    }

    pub fn extension_version(&self) -> Option<i64> {
        self.header.get("EXTVER").as_i64()
    }

    pub fn column_count(&self) -> Result<usize> {
        Ok(self.layout.columns.len())
    }

    pub fn row_count(&self) -> usize {
        self.layout.rows
    }

    pub fn column_info(&self, column: usize) -> Result<ColumnInfo> {
        Ok(self.layout.column(column)?.info.clone())
    }

    pub fn all_column_info(&self) -> Result<Vec<ColumnInfo>> {
        Ok(self.layout.columns.iter().map(|c| c.info.clone()).collect())
    }

    // FITS 4.0 §7.2.2 and §7.3.2: names match without regard to case. An
    // exact match wins over a match that differs only in case.
    pub fn column_by_name(&self, name: &str) -> Result<usize> {
        self.layout
            .find(name)
            .ok_or_else(|| FitsError::InvalidFormat(format!("Column '{}' not found", name)))
    }

    // The stored bytes of one column, row after row.
    pub fn read_column_raw<R: Read + Seek>(
        &self,
        reader: &mut R,
        column: usize,
    ) -> Result<Vec<u8>> {
        let range = self.layout.column(column)?.range.clone();
        let mut out = Vec::with_capacity(range.len().saturating_mul(self.layout.rows));
        self.scan_rows(reader, |row, bytes, _| {
            out.extend_from_slice(cell(bytes, &range, row)?);
            Ok(())
        })?;
        Ok(out)
    }

    pub fn get_row<R: Read + Seek>(&self, reader: &mut R, row: usize) -> Result<Vec<TableValue>> {
        let mut block = Vec::new();
        let rows = self.read_block(reader, row, 1, &mut block)?;
        self.decode_row(rows.start, &block, reader)
    }

    pub fn get_column_values<R: Read + Seek>(
        &self,
        reader: &mut R,
        column: usize,
    ) -> Result<Vec<TableValue>> {
        let col = self.layout.column(column)?;
        let mut out = Vec::with_capacity(self.layout.rows);
        self.scan_rows(reader, |row, bytes, reader| {
            let cell = cell(bytes, &col.range, row)?;
            out.push(F::decode(col, cell, &mut self.heap(reader))?);
            Ok(())
        })?;
        Ok(out)
    }

    pub fn get_column_by_name<R: Read + Seek>(
        &self,
        reader: &mut R,
        name: &str,
    ) -> Result<Vec<TableValue>> {
        let column = self.column_by_name(name)?;
        self.get_column_values(reader, column)
    }

    pub fn iter_rows<'a, R: Read + Seek>(&'a self, reader: &'a mut R) -> TableRows<'a, R, F> {
        TableRows::new(self, reader)
    }

    // Calls `visit` with each row and its bytes, in order. Rows are read
    // a block at a time, so a column read seeks once per block. A heap
    // read in `visit` moves the reader, and the next block seeks back.
    fn scan_rows<R, V>(&self, reader: &mut R, mut visit: V) -> Result<()>
    where
        R: Read + Seek,
        V: FnMut(usize, &[u8], &mut R) -> Result<()>,
    {
        let mut block = Vec::new();
        let mut first = 0;
        while first < self.layout.rows {
            let rows = self.read_block(reader, first, self.layout.rows, &mut block)?;
            for row in rows.clone() {
                visit(row, self.row_in_block(&block, rows.start, row)?, reader)?;
            }
            first = rows.end;
        }
        Ok(())
    }

    // Reads up to `limit` rows from `first` into `block`, capped at
    // BLOCK_BYTES unless one row is longer, and returns the rows read.
    // The rows lie inside the data unit, which the file scan checked.
    pub(super) fn read_block<R: Read + Seek>(
        &self,
        reader: &mut R,
        first: usize,
        limit: usize,
        block: &mut Vec<u8>,
    ) -> Result<Range<usize>> {
        let (rows, row_bytes) = (self.layout.rows, self.layout.row_bytes);
        if first >= rows {
            return Err(FitsError::InvalidFormat(format!(
                "Row index {} out of range (0..{})",
                first, rows
            )));
        }
        let count = (BLOCK_BYTES / row_bytes.max(1))
            .clamp(1, limit)
            .min(rows - first);
        reader.seek(SeekFrom::Start(
            self.info.data_start + (first * row_bytes) as u64,
        ))?;
        block.resize(count * row_bytes, 0);
        reader.read_exact(block)?;
        Ok(first..first + count)
    }

    pub(super) fn row_in_block<'b>(
        &self,
        block: &'b [u8],
        first: usize,
        row: usize,
    ) -> Result<&'b [u8]> {
        let start = (row - first) * self.layout.row_bytes;
        block
            .get(start..start + self.layout.row_bytes)
            .ok_or_else(|| short_row(row))
    }

    // `bytes` holds one whole row.
    pub(super) fn decode_row<R: Read + Seek>(
        &self,
        row: usize,
        bytes: &[u8],
        reader: &mut R,
    ) -> Result<Vec<TableValue>> {
        let mut heap = self.heap(reader);
        self.layout
            .columns
            .iter()
            .map(|column| F::decode(column, cell(bytes, &column.range, row)?, &mut heap))
            .collect()
    }

    fn heap<'a, R>(&self, reader: &'a mut R) -> Heap<'a, R> {
        Heap::new(
            reader,
            self.info.data_start + self.layout.heap.start as u64,
            self.layout.heap.len(),
        )
    }
}

// `range` lies inside the row, which `Layout::new` checked.
fn cell<'b>(row_bytes: &'b [u8], range: &Range<usize>, row: usize) -> Result<&'b [u8]> {
    row_bytes.get(range.clone()).ok_or_else(|| short_row(row))
}

fn short_row(row: usize) -> FitsError {
    FitsError::InvalidFormat(format!("row {} is shorter than its columns", row + 1))
}

impl<F: CellFormat> HduTrait for TableHdu<F> {
    fn header(&self) -> &Header {
        &self.header
    }

    fn info(&self) -> &HduInfo {
        &self.info
    }

    fn hdu_type(&self) -> HduType {
        F::hdu_type()
    }
}
