use super::tform::{be, descriptor, heap_range, BinaryForm, Element, Number};
use crate::fits::data::scaling::{I8_OFFSET, U16_OFFSET, U32_OFFSET, U64_OFFSET};
use crate::fits::data::value::TableValue;
use crate::fits::errors::{FitsError, Result};
use crate::fits::hdu::table::cell::{CellFormat, Heap};
use crate::fits::hdu::table::layout::Column;
use crate::fits::hdu::table::values::{complex, integer, real, text};
use crate::fits::hdu::HduType;
use crate::fits::header::parser::Header;
use std::io::{Read, Seek};
use std::ops::Range;

impl CellFormat for BinaryForm {
    const PACKED: bool = true;

    fn hdu_type() -> HduType {
        HduType::BinaryTable
    }

    fn place(tform: &str, _: &Header, n: usize, next: usize) -> Result<(Self, Range<usize>)> {
        let (form, width) = BinaryForm::parse(tform)?;
        let end = next
            .checked_add(width)
            .ok_or_else(|| invalid(format!("column {n} ends past the largest row size")))?;
        Ok((form, next..end))
    }

    fn decode<R: Read + Seek>(
        column: &Column<Self>,
        cell: &[u8],
        heap: &mut Heap<'_, R>,
    ) -> Result<TableValue> {
        match column.format {
            BinaryForm::Values { repeat, element } => {
                cells(column, element, cell, repeat, Shape::Fixed)
            }
            BinaryForm::Descriptor { wide, element } => {
                let (count, bytes) = heap_array(column, cell, wide, element, heap)?;
                cells(column, element, &bytes, count, Shape::Array)
            }
        }
    }
}

// A fixed column with repeat 1 holds a scalar. Variable-length arrays are
// arrays even with one element.
#[derive(PartialEq)]
enum Shape {
    Fixed,
    Array,
}

// The element count and bytes of one variable-length array. Every
// descriptor is checked against the heap size before the bytes are read.
fn heap_array<R: Read + Seek>(
    column: &Column<BinaryForm>,
    cell: &[u8],
    wide: bool,
    element: Element,
    heap: &mut Heap<'_, R>,
) -> Result<(usize, Vec<u8>)> {
    let name = || column.info.name.clone().unwrap_or_default();
    let (count, offset) = descriptor(cell, wide)
        .ok_or_else(|| invalid(format!("column {} has a negative descriptor", name())))?;
    let range = heap_range(count, offset, element, heap.len()).ok_or_else(|| {
        invalid(format!(
            "column {} points to {count} elements at heap byte {offset}, past the {}-byte heap",
            name(),
            heap.len()
        ))
    })?;
    Ok((count, heap.read(range)?))
}

// `bytes` holds exactly `count` elements, which the column width or the
// heap range guarantees.
fn cells<F>(
    column: &Column<F>,
    element: Element,
    bytes: &[u8],
    count: usize,
    shape: Shape,
) -> Result<TableValue> {
    match element {
        Element::Char => Ok(text(bytes)),
        Element::Bit => Ok(TableValue::Bits(
            (0..count)
                .map(|i| bytes.get(i / 8).is_some_and(|b| b & (0x80 >> (i % 8)) != 0))
                .collect(),
        )),
        _ => {
            let mut values = bytes
                .chunks_exact(element.size())
                .map(|b| scalar(column, element, b));
            match (count, shape) {
                (1, Shape::Fixed) => values.next().unwrap_or(Ok(TableValue::Null)),
                _ => Ok(TableValue::Array(values.collect::<Result<_>>()?)),
            }
        }
    }
}

// One element from exactly `element.size()` bytes.
fn scalar<F>(column: &Column<F>, element: Element, b: &[u8]) -> Result<TableValue> {
    let short = || invalid(format!("{element:?} element of {} bytes", b.len()));
    Ok(match element {
        Element::Logical => logical(b.first().copied().ok_or_else(short)?),
        Element::C32 => {
            let re = f32::from_be_bytes(be(b).ok_or_else(short)?);
            let im = f32::from_be_bytes(be(b.get(4..).ok_or_else(short)?).ok_or_else(short)?);
            complex(
                column,
                f64::from(re),
                f64::from(im),
                TableValue::Complex32(re, im),
            )
        }
        Element::C64 => {
            let re = f64::from_be_bytes(be(b).ok_or_else(short)?);
            let im = f64::from_be_bytes(be(b.get(8..).ok_or_else(short)?).ok_or_else(short)?);
            complex(column, re, im, TableValue::Complex64(re, im))
        }
        Element::Bit | Element::Char => return Err(short()),
        _ => number(column, Number::read(element, b).ok_or_else(short)?),
    })
}

fn number<F>(column: &Column<F>, number: Number) -> TableValue {
    match number {
        Number::U8(v) => integer(
            column,
            v,
            TableValue::Byte,
            Some((I8_OFFSET, |v| TableValue::I8((v ^ 0x80) as i8))),
        ),
        Number::I16(v) => integer(
            column,
            v,
            TableValue::I16,
            Some((U16_OFFSET, |v| TableValue::U16((v as u16) ^ 0x8000))),
        ),
        Number::I32(v) => integer(
            column,
            v,
            TableValue::I32,
            Some((U32_OFFSET, |v| TableValue::U32((v as u32) ^ (1 << 31)))),
        ),
        Number::I64(v) => integer(
            column,
            v,
            TableValue::I64,
            Some((U64_OFFSET, |v| TableValue::U64((v as u64) ^ (1 << 63)))),
        ),
        Number::F32(v) => real(column, f64::from(v), TableValue::F32(v)),
        Number::F64(v) => real(column, v, TableValue::F64(v)),
    }
}

// FITS 4.0 §7.3.3.1: T is true, F is false and a zero byte is null.
fn logical(byte: u8) -> TableValue {
    match byte {
        b'T' => TableValue::Logical(true),
        b'F' => TableValue::Logical(false),
        // Other bytes, such as the 1 some writers store for true, have no
        // defined value, and an error would lose the rest of the row.
        _ => TableValue::Null,
    }
}

fn invalid(message: String) -> FitsError {
    FitsError::InvalidFormat(message)
}
