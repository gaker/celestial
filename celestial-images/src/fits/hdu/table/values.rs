use super::layout::Column;
use crate::fits::data::value::TableValue;

// A Table 19 offset and the value of a stored integer shifted by it.
pub(crate) type Offset<T> = (i128, fn(T) -> TableValue);

// A stored integer and its physical value. FITS 4.0 §7.3.2: a value equal
// to TNULLn is undefined. With a scale of 1 and the Table 19 `offset`, the
// value keeps its width as `shifted`. Any other scaling gives
// TZEROn + TSCALn * stored as an f64.
pub(crate) fn integer<T: Copy + Into<i64>, F>(
    column: &Column<F>,
    stored: T,
    plain: fn(T) -> TableValue,
    offset: Option<Offset<T>>,
) -> TableValue {
    let value: i64 = stored.into();
    if column.null == Some(value) {
        return TableValue::Null;
    }
    if column.scale.is_identity() {
        return plain(stored);
    }
    match offset {
        Some((zero, shifted)) if column.scale.is_offset(zero) => shifted(stored),
        _ => TableValue::F64(column.scale.apply(value as f64)),
    }
}

// NaN stays a value: FITS 4.0 §7.3.3.1 makes it the null of float columns,
// and the caller can test for it.
pub(crate) fn real<F>(column: &Column<F>, stored: f64, plain: TableValue) -> TableValue {
    if column.scale.is_identity() {
        plain
    } else {
        TableValue::F64(column.scale.apply(stored))
    }
}

// FITS 4.0 §7.3.2: TZEROn offsets the real part only, and TSCALn scales
// both parts.
pub(crate) fn complex<F>(column: &Column<F>, re: f64, im: f64, plain: TableValue) -> TableValue {
    if column.scale.is_identity() {
        plain
    } else {
        TableValue::Complex64(column.scale.apply(re), column.scale.scale() * im)
    }
}

// Text ends at the first NUL (FITS 4.0 §7.3.3.1), and trailing spaces are
// padding.
pub(crate) fn text(bytes: &[u8]) -> TableValue {
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    let text = String::from_utf8_lossy(&bytes[..end]);
    TableValue::String(text.trim_end_matches(' ').to_string())
}
