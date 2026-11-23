use super::cell::CellFormat;
use crate::fits::data::scaling::{integer_value, Linear};
use crate::fits::errors::{FitsError, Result};
use crate::fits::hdu::ColumnInfo;
use crate::fits::header::keywords::{Keyword, KeywordValue};
use crate::fits::header::parser::Header;
use std::ops::Range;

// One table column: its keywords, where its cells lie in a row, and how
// stored values become physical values.
#[derive(Debug)]
pub struct Column<F> {
    pub(crate) info: ColumnInfo,
    pub(crate) range: Range<usize>,
    pub(crate) format: F,
    pub(crate) scale: Linear,
    // TNULLn when it holds an integer. FITS 4.0 §7.3.2 compares it with
    // the stored value, before scaling.
    pub(crate) null: Option<i64>,
}

// The rows, columns and heap of a table data unit, checked against the
// header and the data size.
#[derive(Debug)]
pub(crate) struct Layout<F> {
    pub(crate) row_bytes: usize,
    pub(crate) rows: usize,
    // Byte offsets from the start of the data unit.
    pub(crate) heap: Range<usize>,
    pub(crate) columns: Vec<Column<F>>,
}

impl<F: CellFormat> Layout<F> {
    // `data_len` is the size of the data unit. The main table and PCOUNT
    // bytes after it must fit inside it.
    pub(crate) fn new(header: &Header, data_len: usize) -> Result<Self> {
        let row_bytes = count(header, "NAXIS1")?;
        let rows = count(header, "NAXIS2")?;
        let pcount = match header.get_keyword("PCOUNT") {
            None => 0,
            Some(_) => count(header, "PCOUNT")?,
        };
        let too_big = || {
            invalid(format!(
                "{rows} rows of {row_bytes} bytes and a {pcount}-byte heap exceed the {data_len}-byte data unit"
            ))
        };
        let main = row_bytes.checked_mul(rows).ok_or_else(too_big)?;
        let end = main
            .checked_add(pcount)
            .filter(|&end| end <= data_len)
            .ok_or_else(too_big)?;
        Ok(Self {
            row_bytes,
            rows,
            heap: heap(header, main, end)?,
            columns: columns(header, row_bytes)?,
        })
    }

    pub(crate) fn column(&self, index: usize) -> Result<&Column<F>> {
        self.columns.get(index).ok_or_else(|| {
            invalid(format!(
                "Column index {} out of range (0..{})",
                index,
                self.columns.len()
            ))
        })
    }

    pub(crate) fn find(&self, name: &str) -> Option<usize> {
        let named = |matches: &dyn Fn(&str) -> bool| {
            self.columns
                .iter()
                .position(|c| c.info.name.as_deref().is_some_and(matches))
        };
        named(&|n| n == name).or_else(|| named(&|n| n.eq_ignore_ascii_case(name)))
    }
}

// FITS 4.0 §7.3.5: the heap starts at THEAP, by default right after the
// main table, and ends PCOUNT bytes after the main table.
fn heap(header: &Header, main: usize, end: usize) -> Result<Range<usize>> {
    let start = match header.get_last_keyword("THEAP") {
        None => main,
        Some(_) => last_count(header, "THEAP")?,
    };
    if start < main || start > end {
        return Err(invalid(format!(
            "THEAP {start} lies outside the {main}..{end} bytes after the rows"
        )));
    }
    Ok(start..end)
}

fn columns<F: CellFormat>(header: &Header, row_bytes: usize) -> Result<Vec<Column<F>>> {
    let fields = count(header, "TFIELDS")?;
    if fields > 999 {
        return Err(invalid(format!("TFIELDS = {fields} is above 999")));
    }
    let mut out = Vec::with_capacity(fields);
    let mut next = 0;
    for n in 1..=fields {
        let column = column(header, n, next)?;
        if column.range.end > row_bytes {
            return Err(invalid(format!(
                "column {n} ends at byte {}, past NAXIS1 = {row_bytes}",
                column.range.end
            )));
        }
        next = column.range.end;
        out.push(column);
    }
    if F::PACKED && next != row_bytes {
        return Err(invalid(format!(
            "columns take {next} bytes, but NAXIS1 = {row_bytes}"
        )));
    }
    Ok(out)
}

fn column<F: CellFormat>(header: &Header, n: usize, next: usize) -> Result<Column<F>> {
    let key = format!("TFORM{n}");
    let tform = header
        .get(&key)
        .as_str()
        .ok_or(FitsError::KeywordNotFound { keyword: key })?;
    let (format, range) = F::place(tform, header, n, next)?;
    let null = header.get_last(&format!("TNULL{n}")).value();
    Ok(Column {
        info: info(header, n, tform, null),
        range,
        format,
        scale: Linear::from_header(header, &format!("TSCAL{n}"), &format!("TZERO{n}"))?,
        null: null.and_then(integer_value),
    })
}

fn info(header: &Header, n: usize, tform: &str, null: Option<&KeywordValue>) -> ColumnInfo {
    let text = |key: &str| {
        header
            .get(&format!("{key}{n}"))
            .as_str()
            .map(str::to_string)
    };
    let real = |key: &str| header.get_last(&format!("{key}{n}")).as_f64();
    let mut info = ColumnInfo::new(n - 1, tform.to_string());
    info.name = text("TTYPE");
    info.unit = text("TUNIT");
    info.null_value = null.map(|v| match v {
        KeywordValue::String(s) => s.clone(),
        other => other.to_string(),
    });
    info.scale = real("TSCAL");
    info.zero_offset = real("TZERO");
    info.display_format = text("TDISP");
    info
}

pub(crate) fn count(header: &Header, name: &str) -> Result<usize> {
    count_of(name, header.get_keyword(name))
}

// THEAP and TBCOLn take their last card, the one cfitsio reads the table
// with. The mandatory counts cannot repeat, since the scan rejects them.
pub(crate) fn last_count(header: &Header, name: &str) -> Result<usize> {
    count_of(name, header.get_last_keyword(name))
}

fn count_of(name: &str, keyword: Option<&Keyword>) -> Result<usize> {
    let keyword = keyword.ok_or_else(|| FitsError::KeywordNotFound {
        keyword: name.to_string(),
    })?;
    keyword
        .value
        .as_ref()
        .and_then(|v| v.as_integer())
        .and_then(|v| usize::try_from(v).ok())
        .ok_or_else(|| FitsError::InvalidKeywordValue {
            keyword: name.to_string(),
            value: format!("{:?}", keyword.value),
        })
}

fn invalid(message: String) -> FitsError {
    FitsError::InvalidFormat(message)
}
