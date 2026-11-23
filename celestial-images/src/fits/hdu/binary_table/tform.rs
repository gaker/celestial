use crate::fits::errors::{FitsError, Result};
use std::ops::Range;

// The element types of FITS 4.0 Table 18, other than the P and Q array
// descriptors.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Element {
    Logical,
    Bit,
    U8,
    I16,
    I32,
    I64,
    Char,
    F32,
    F64,
    C32,
    C64,
}

impl Element {
    // FITS 4.0 §7.3.1: type codes are upper case.
    fn from_code(code: u8) -> Result<Self> {
        Ok(match code {
            b'L' => Self::Logical,
            b'X' => Self::Bit,
            b'B' => Self::U8,
            b'I' => Self::I16,
            b'J' => Self::I32,
            b'K' => Self::I64,
            b'A' => Self::Char,
            b'E' => Self::F32,
            b'D' => Self::F64,
            b'C' => Self::C32,
            b'M' => Self::C64,
            _ => return Err(invalid(format!("unknown TFORM type '{}'", code as char))),
        })
    }

    // The bytes of one element. A bit takes one byte only as a whole
    // array, through `bytes`.
    pub(crate) fn size(self) -> usize {
        match self {
            Self::Logical | Self::Bit | Self::U8 | Self::Char => 1,
            Self::I16 => 2,
            Self::I32 | Self::F32 => 4,
            Self::I64 | Self::F64 | Self::C32 => 8,
            Self::C64 => 16,
        }
    }

    // The bytes of `count` elements. Bits are packed eight to a byte,
    // and the last byte is padded (Table 18).
    pub(crate) fn bytes(self, count: usize) -> Option<usize> {
        match self {
            Self::Bit => Some(count.div_ceil(8)),
            _ => count.checked_mul(self.size()),
        }
    }
}

// One element of a numeric column, read from its big-endian bytes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Number {
    U8(u8),
    I16(i16),
    I32(i32),
    I64(i64),
    F32(f32),
    F64(f64),
}

impl Number {
    // None when `element` is not a real or integer type, or when `bytes`
    // is shorter than one element.
    pub(crate) fn read(element: Element, bytes: &[u8]) -> Option<Self> {
        Some(match element {
            Element::U8 => Self::U8(*bytes.first()?),
            Element::I16 => Self::I16(i16::from_be_bytes(be(bytes)?)),
            Element::I32 => Self::I32(i32::from_be_bytes(be(bytes)?)),
            Element::I64 => Self::I64(i64::from_be_bytes(be(bytes)?)),
            Element::F32 => Self::F32(f32::from_be_bytes(be(bytes)?)),
            Element::F64 => Self::F64(f64::from_be_bytes(be(bytes)?)),
            _ => return None,
        })
    }
}

// The TFORMn of one BINTABLE column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinaryForm {
    Values { repeat: usize, element: Element },
    // A P (32-bit) or Q (64-bit, `wide`) array descriptor into the heap.
    Descriptor { wide: bool, element: Element },
}

impl BinaryForm {
    // FITS 4.0 §7.3.1 TFORMn = rTa: an optional repeat count, the type
    // code, and for P and Q the element type and an optional maximum
    // length. Returns the form and its width in the row.
    pub(crate) fn parse(tform: &str) -> Result<(Self, usize)> {
        let (repeat, code, array_type) = repeat_and_type(tform)?;
        if code == b'P' || code == b'Q' {
            return descriptor_form(tform, repeat, code == b'Q', array_type);
        }
        let element = Element::from_code(code)?;
        let width = element
            .bytes(repeat)
            .ok_or_else(|| invalid(format!("TFORM '{tform}' is too wide")))?;
        Ok((Self::Values { repeat, element }, width))
    }
}

// FITS 4.0 §7.3.5: the repeat count of a descriptor column is 0 or 1. A
// repeat of 0 is an empty field with no descriptor to read.
fn descriptor_form(
    tform: &str,
    repeat: usize,
    wide: bool,
    code: u8,
) -> Result<(BinaryForm, usize)> {
    let element = Element::from_code(code)?;
    match repeat {
        0 => Ok((BinaryForm::Values { repeat, element }, 0)),
        1 => Ok((
            BinaryForm::Descriptor { wide, element },
            if wide { 16 } else { 8 },
        )),
        _ => Err(invalid(format!(
            "TFORM '{tform}' must be 0 or 1 descriptors"
        ))),
    }
}

// The repeat count, 1 when absent, the type code, and the byte after it,
// which is the element type of a P or Q column.
fn repeat_and_type(tform: &str) -> Result<(usize, u8, u8)> {
    let text = tform.trim();
    let digits = text.bytes().take_while(u8::is_ascii_digit).count();
    let repeat = match digits {
        0 => 1,
        _ => text[..digits]
            .parse::<usize>()
            .map_err(|_| invalid(format!("TFORM '{tform}' has a bad repeat count")))?,
    };
    let mut rest = text.bytes().skip(digits);
    let code = rest
        .next()
        .ok_or_else(|| invalid(format!("TFORM '{tform}' has no type")))?;
    Ok((repeat, code, rest.next().unwrap_or(b' ')))
}

// The element count and heap offset of a P (8-byte) or Q (16-byte)
// descriptor. FITS 4.0 §7.3.5 leaves negative values undefined, so they
// give None.
pub(crate) fn descriptor(cell: &[u8], wide: bool) -> Option<(usize, usize)> {
    let (count, offset) = if wide {
        (
            i64::from_be_bytes(be(cell.get(..8)?)?),
            i64::from_be_bytes(be(cell.get(8..16)?)?),
        )
    } else {
        (
            i64::from(i32::from_be_bytes(be(cell.get(..4)?)?)),
            i64::from(i32::from_be_bytes(be(cell.get(4..8)?)?)),
        )
    };
    Some((usize::try_from(count).ok()?, usize::try_from(offset).ok()?))
}

// The heap bytes of `count` elements at `offset`, when they lie entirely
// inside a heap of `heap_len` bytes (FITS 4.0 §7.3.5).
pub(crate) fn heap_range(
    count: usize,
    offset: usize,
    element: Element,
    heap_len: usize,
) -> Option<Range<usize>> {
    let end = offset.checked_add(element.bytes(count)?)?;
    (end <= heap_len).then_some(offset..end)
}

// The first N bytes as an array, or None when there are fewer.
pub(crate) fn be<const N: usize>(bytes: &[u8]) -> Option<[u8; N]> {
    bytes.get(..N)?.try_into().ok()
}

fn invalid(message: String) -> FitsError {
    FitsError::InvalidFormat(message)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn width(tform: &str) -> usize {
        BinaryForm::parse(tform).unwrap().1
    }

    #[test]
    fn widths_follow_table_18() {
        assert_eq!(width("1PB(102)"), 8);
        assert_eq!(width("1QE(23)"), 16);
        assert_eq!(width("PJ"), 8);
        assert_eq!(width("0PB"), 0);
        assert_eq!(width("0QB(5)"), 0);
        assert_eq!(width("1D"), 8);
        assert_eq!(width("3E"), 12);
        assert_eq!(width("16X"), 2);
        assert_eq!(width("9X"), 2);
        assert_eq!(width("X"), 1);
        assert_eq!(width("0X"), 0);
        assert_eq!(width("20A"), 20);
        assert_eq!(width("2M"), 32);
        assert_eq!(width("3C"), 24);
        assert_eq!(width("L"), 1);
        assert_eq!(width("4K"), 32);
        assert_eq!(width("10Ab"), 10);
    }

    #[test]
    fn forms_name_their_element() {
        let form = |t: &str| BinaryForm::parse(t).unwrap().0;
        assert_eq!(
            form("3E"),
            BinaryForm::Values {
                repeat: 3,
                element: Element::F32
            }
        );
        assert_eq!(
            form("QJ(4)"),
            BinaryForm::Descriptor {
                wide: true,
                element: Element::I32
            }
        );
        assert_eq!(
            form("0PJ"),
            BinaryForm::Values {
                repeat: 0,
                element: Element::I32
            }
        );
    }

    #[test]
    fn bad_forms_are_errors() {
        for tform in [
            "2PB",
            "1Z",
            "",
            "42",
            "1PP",
            "1P",
            "j",
            "99999999999999999999999J",
        ] {
            assert!(BinaryForm::parse(tform).is_err(), "{tform}");
        }
        assert!(BinaryForm::parse(&format!("{}D", usize::MAX)).is_err());
    }

    #[test]
    fn descriptors_reject_negative_values_and_short_cells() {
        let mut p = 3i32.to_be_bytes().to_vec();
        p.extend_from_slice(&5i32.to_be_bytes());
        assert_eq!(descriptor(&p, false), Some((3, 5)));
        assert_eq!(descriptor(&p[..7], false), None);
        let mut q = (-1i64).to_be_bytes().to_vec();
        q.extend_from_slice(&0i64.to_be_bytes());
        assert_eq!(descriptor(&q, true), None);
    }

    #[test]
    fn heap_ranges_stay_inside_the_heap() {
        assert_eq!(heap_range(2, 4, Element::I32, 12), Some(4..12));
        assert_eq!(heap_range(2, 5, Element::I32, 12), None);
        assert_eq!(heap_range(17, 0, Element::Bit, 3), Some(0..3));
        assert_eq!(heap_range(usize::MAX, 0, Element::I64, 12), None);
        assert_eq!(heap_range(1, usize::MAX, Element::U8, 12), None);
        assert_eq!(heap_range(0, 12, Element::U8, 12), Some(12..12));
    }
}
