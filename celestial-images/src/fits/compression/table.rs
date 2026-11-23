use crate::fits::errors::{FitsError, Result};
use crate::fits::hdu::binary_table::tform::{descriptor, heap_range, BinaryForm, Element, Number};
use crate::fits::hdu::table::layout::{Column, Layout};
use crate::fits::header::parser::Header;

// The rows and heap of a compressed-image table, laid out by the shared
// binary-table layout. Every read is checked against the data unit, so a
// bad descriptor returns an error.
pub(super) struct Table<'a> {
    data: &'a [u8],
    layout: Layout<BinaryForm>,
}

pub(super) struct HeapArray<'a> {
    pub(super) bytes: &'a [u8],
    pub(super) element: Element,
}

impl<'a> Table<'a> {
    pub(super) fn new(header: &Header, data: &'a [u8]) -> Result<Self> {
        Ok(Self {
            data,
            layout: Layout::new(header, data.len())?,
        })
    }

    pub(super) fn rows(&self) -> usize {
        self.layout.rows
    }

    // FITS 4.0 §10.1.3: column names are matched without regard to case.
    pub(super) fn column(&self, name: &str) -> Option<usize> {
        self.layout.find(name)
    }

    pub(super) fn is_array(&self, column: usize) -> bool {
        self.layout
            .columns
            .get(column)
            .is_some_and(|c| matches!(c.format, BinaryForm::Descriptor { .. }))
    }

    pub(super) fn array(&self, row: usize, column: usize) -> Result<HeapArray<'a>> {
        let col = self.layout.column(column)?;
        let BinaryForm::Descriptor { wide, element } = col.format else {
            return Err(invalid(format!(
                "column {} is not an array column",
                name(col)
            )));
        };
        let (count, offset) = descriptor(self.cell(row, col)?, wide).ok_or_else(|| {
            invalid(format!(
                "row {} of {} has a negative descriptor",
                row + 1,
                name(col)
            ))
        })?;
        let heap = self.data.get(self.layout.heap.clone()).unwrap_or_default();
        let bytes = heap_range(count, offset, element, heap.len())
            .and_then(|range| heap.get(range))
            .ok_or_else(|| {
                invalid(format!(
                    "row {} of {} points past the heap",
                    row + 1,
                    name(col)
                ))
            })?;
        Ok(HeapArray { bytes, element })
    }

    pub(super) fn real(&self, row: usize, column: usize) -> Result<f64> {
        Ok(match self.scalar(row, column)? {
            Number::U8(v) => f64::from(v),
            Number::I16(v) => f64::from(v),
            Number::I32(v) => f64::from(v),
            Number::I64(v) => v as f64,
            Number::F32(v) => f64::from(v),
            Number::F64(v) => v,
        })
    }

    pub(super) fn integer(&self, row: usize, column: usize) -> Result<i64> {
        match self.scalar(row, column)? {
            Number::U8(v) => Ok(i64::from(v)),
            Number::I16(v) => Ok(i64::from(v)),
            Number::I32(v) => Ok(i64::from(v)),
            Number::I64(v) => Ok(v),
            Number::F32(_) | Number::F64(_) => Err(invalid(format!(
                "column {} is not an integer column",
                name(self.layout.column(column)?)
            ))),
        }
    }

    fn scalar(&self, row: usize, column: usize) -> Result<Number> {
        let col = self.layout.column(column)?;
        match col.format {
            BinaryForm::Values { repeat: 1, element } if is_number(element) => {
                Number::read(element, self.cell(row, col)?)
                    .ok_or_else(|| invalid(format!("row {} is shorter than its columns", row + 1)))
            }
            _ => Err(invalid(format!(
                "column {} does not hold one number per row",
                name(col)
            ))),
        }
    }

    fn cell(&self, row: usize, col: &Column<BinaryForm>) -> Result<&'a [u8]> {
        if row >= self.layout.rows {
            return Err(invalid(format!(
                "row {} is past the {} table rows",
                row + 1,
                self.layout.rows
            )));
        }
        let start = row * self.layout.row_bytes;
        self.data
            .get(start + col.range.start..start + col.range.end)
            .ok_or_else(|| invalid(format!("row {} is shorter than its columns", row + 1)))
    }
}

fn is_number(element: Element) -> bool {
    matches!(
        element,
        Element::U8 | Element::I16 | Element::I32 | Element::I64 | Element::F32 | Element::F64
    )
}

fn name(column: &Column<BinaryForm>) -> &str {
    column.info.name.as_deref().unwrap_or("")
}

fn invalid(message: String) -> FitsError {
    FitsError::InvalidFormat(message)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::Keyword;

    fn header(columns: &[(&str, &str)], row_bytes: i64, rows: i64) -> Header {
        let mut h = Header::new();
        h.add_keyword(Keyword::integer("NAXIS1", row_bytes));
        h.add_keyword(Keyword::integer("NAXIS2", rows));
        h.add_keyword(Keyword::integer("TFIELDS", columns.len() as i64));
        for (i, (name, form)) in columns.iter().enumerate() {
            h.add_keyword(Keyword::string(format!("TTYPE{}", i + 1), name.to_string()));
            h.add_keyword(Keyword::string(format!("TFORM{}", i + 1), form.to_string()));
        }
        h
    }

    #[test]
    fn reads_descriptors_scalars_and_heap() {
        let mut h = header(
            &[
                ("COMPRESSED_DATA", "1PI"),
                ("ZSCALE", "1D"),
                ("ZBLANK", "1J"),
            ],
            20,
            1,
        );
        h.add_keyword(Keyword::integer("PCOUNT", 5));
        let mut data = Vec::new();
        data.extend_from_slice(&2i32.to_be_bytes());
        data.extend_from_slice(&1i32.to_be_bytes());
        data.extend_from_slice(&2.5f64.to_be_bytes());
        data.extend_from_slice(&(-7i32).to_be_bytes());
        data.extend_from_slice(&[9, 0, 1, 0, 2]);
        let t = Table::new(&h, &data).unwrap();
        let col = t.column("compressed_data").unwrap();
        let array = t.array(0, col).unwrap();
        assert_eq!(array.bytes, &[0, 1, 0, 2]);
        assert_eq!(array.element, Element::I16);
        assert_eq!(t.real(0, 1).unwrap(), 2.5);
        assert_eq!(t.integer(0, 2).unwrap(), -7);
        assert_eq!(t.real(0, 2).unwrap(), -7.0);
        assert!(t.integer(0, 1).is_err());
        assert!(t.array(1, col).is_err());
        assert!(t.real(0, col).is_err());
    }

    #[test]
    fn honours_theap() {
        let mut h = header(&[("COMPRESSED_DATA", "1PB")], 8, 1);
        h.add_keyword(Keyword::integer("THEAP", 10));
        h.add_keyword(Keyword::integer("PCOUNT", 3));
        let mut data = Vec::new();
        data.extend_from_slice(&1i32.to_be_bytes());
        data.extend_from_slice(&0i32.to_be_bytes());
        data.extend_from_slice(&[0, 0, 42]);
        let t = Table::new(&h, &data).unwrap();
        assert_eq!(t.array(0, 0).unwrap().bytes, &[42]);
    }

    #[test]
    fn rejects_bad_layouts() {
        let data = [0u8; 16];
        let wide = header(&[("COMPRESSED_DATA", "1PB"), ("ZSCALE", "1D")], 8, 1);
        assert!(Table::new(&wide, &data).is_err());
        let many_rows = header(&[("COMPRESSED_DATA", "1PB")], 8, 3);
        assert!(Table::new(&many_rows, &data).is_err());
        let mut theap = header(&[("COMPRESSED_DATA", "1PB")], 8, 1);
        theap.add_keyword(Keyword::integer("THEAP", 4));
        assert!(Table::new(&theap, &data).is_err());
        let negative = header(&[("COMPRESSED_DATA", "1PB")], -8, 1);
        assert!(Table::new(&negative, &data).is_err());
    }

    #[test]
    fn rejects_descriptors_outside_the_heap() {
        let mut h = header(&[("COMPRESSED_DATA", "1QB")], 16, 1);
        h.add_keyword(Keyword::integer("PCOUNT", 4));
        for (len, offset) in [(5i64, 0i64), (1, 4), (-1, 0), (0, -1), (i64::MAX, 0)] {
            let mut data = Vec::new();
            data.extend_from_slice(&len.to_be_bytes());
            data.extend_from_slice(&offset.to_be_bytes());
            data.extend_from_slice(&[1, 2, 3, 4]);
            let t = Table::new(&h, &data).unwrap();
            assert!(t.array(0, 0).is_err(), "{len} {offset}");
        }
    }
}
