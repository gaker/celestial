use super::table::cell::{CellFormat, Heap};
use super::table::layout::{last_count, Column};
use super::table::rows::TableRows;
use super::table::values::{integer, real, text};
use super::table::TableHdu;
use super::HduType;
use crate::fits::data::value::TableValue;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::parser::Header;
use std::io::{Read, Seek};
use std::ops::Range;

pub type AsciiTableHdu = TableHdu<AsciiForm>;
pub type AsciiTableRowIterator<'a, R> = TableRows<'a, R, AsciiForm>;

// The field formats of FITS 4.0 Table 15. `decimals` is the d of Fw.d,
// Ew.d and Dw.d.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsciiForm {
    Char,
    Integer,
    Real { decimals: usize },
}

impl CellFormat for AsciiForm {
    const PACKED: bool = false;

    fn hdu_type() -> HduType {
        HduType::AsciiTable
    }

    // FITS 4.0 §7.2.1: field n starts at the 1-based column TBCOLn.
    fn place(tform: &str, header: &Header, n: usize, _: usize) -> Result<(Self, Range<usize>)> {
        let (form, width) = parse_tform(tform)?;
        let start = last_count(header, &format!("TBCOL{n}"))?
            .checked_sub(1)
            .ok_or_else(|| invalid(format!("TBCOL{n} = 0; fields start at column 1")))?;
        let end = start
            .checked_add(width)
            .ok_or_else(|| invalid(format!("field {n} ends past the largest row size")))?;
        Ok((form, start..end))
    }

    fn decode<R: Read + Seek>(
        column: &Column<Self>,
        cell: &[u8],
        _: &mut Heap<'_, R>,
    ) -> Result<TableValue> {
        if is_null(column, cell) {
            return Ok(TableValue::Null);
        }
        Ok(match column.format {
            AsciiForm::Char => text(cell),
            AsciiForm::Integer => integer(column, integer_field(cell)?, TableValue::I64, None),
            AsciiForm::Real { decimals } => {
                let v = real_field(cell, decimals)?;
                real(column, v, TableValue::F64(v))
            }
        })
    }
}

// FITS 4.0 §7.2.2: TNULLn is space filled to the field width.
fn is_null(column: &Column<AsciiForm>, cell: &[u8]) -> bool {
    column
        .info
        .null_value
        .as_ref()
        .is_some_and(|null| cell.trim_ascii_end() == null.as_bytes().trim_ascii_end())
}

// Aw and Iw, or Fw.d, Ew.d and Dw.d. A real format without .d reads as d = 0.
fn parse_tform(tform: &str) -> Result<(AsciiForm, usize)> {
    let bad = || invalid(format!("TFORM '{tform}' is not Aw, Iw, Fw.d, Ew.d or Dw.d"));
    let text = tform.trim();
    let (code, rest) = text.split_at_checked(1).ok_or_else(bad)?;
    let (width, decimals) = match rest.split_once('.') {
        Some((w, d)) => (w, Some(d)),
        None => (rest, None),
    };
    let width = digits(width).filter(|&w| w > 0).ok_or_else(bad)?;
    let form = match (code, decimals) {
        ("A", None) => AsciiForm::Char,
        ("I", None) => AsciiForm::Integer,
        ("F" | "E" | "D", None) => AsciiForm::Real { decimals: 0 },
        ("F" | "E" | "D", Some(d)) => AsciiForm::Real {
            decimals: digits(d).ok_or_else(bad)?,
        },
        _ => return Err(bad()),
    };
    Ok((form, width))
}

fn digits(text: &str) -> Option<usize> {
    if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    text.parse().ok()
}

// The field without the spaces around it.
fn field_text(cell: &[u8]) -> Result<&str> {
    let text = std::str::from_utf8(cell)
        .map_err(|_| invalid("numeric field holds bytes that are not ASCII".to_string()))?;
    Ok(text.trim_matches(' '))
}

// FITS 4.0 §7.2.5: an optional sign and decimal digits, with spaces
// before and after. A blank field is 0.
fn integer_field(cell: &[u8]) -> Result<i64> {
    let text = field_text(cell)?;
    if text.is_empty() {
        return Ok(0);
    }
    let unsigned = text.strip_prefix(['+', '-']).unwrap_or(text);
    if unsigned.is_empty() || !unsigned.bytes().all(|b| b.is_ascii_digit()) {
        return Err(invalid(format!("integer field holds '{text}'")));
    }
    text.parse()
        .map_err(|_| invalid(format!("integer field '{text}' does not fit 64 bits")))
}

// FITS 4.0 §7.2.5: [sign] digits [. digits] [exponent]. The exponent is E or
// D and a signed integer, or a sign and an integer. Without a decimal
// point, the last `decimals` digits are the fraction. A blank field is 0.
// The value is the decimal text rounded once to f64.
fn real_field(cell: &[u8], decimals: usize) -> Result<f64> {
    let text = field_text(cell)?;
    if text.is_empty() {
        return Ok(0.0);
    }
    let bad = || invalid(format!("real field holds '{text}'"));
    let (sign, unsigned) = match text.strip_prefix('-') {
        Some(rest) => ("-", rest),
        None => ("", text.strip_prefix('+').unwrap_or(text)),
    };
    let split = unsigned
        .find(|c: char| !c.is_ascii_digit() && c != '.')
        .unwrap_or(unsigned.len());
    let (mantissa, exponent) = unsigned.split_at(split);
    let (digits, shift) = mantissa_digits(mantissa, decimals).ok_or_else(bad)?;
    let exponent = exponent_value(exponent).ok_or_else(bad)?;
    let power = i64::try_from(shift)
        .ok()
        .and_then(|shift| exponent.checked_sub(shift))
        .ok_or_else(bad)?;
    format!("{sign}{digits}e{power}").parse().map_err(|_| bad())
}

// Returns the digits of a mantissa and how many of them follow the decimal
// point, or `None` when the mantissa has no digits or a second point.
fn mantissa_digits(mantissa: &str, decimals: usize) -> Option<(String, usize)> {
    let (whole, fraction, shift) = match mantissa.split_once('.') {
        Some((w, f)) => (w, f, f.len()),
        None => (mantissa, "", decimals),
    };
    let all_digits = |s: &str| s.bytes().all(|b| b.is_ascii_digit());
    if whole.len() + fraction.len() == 0 || !all_digits(fraction) {
        return None;
    }
    Some((format!("{whole}{fraction}"), shift))
}

fn exponent_value(text: &str) -> Option<i64> {
    let signed = match text.as_bytes().first() {
        None => return Some(0),
        Some(b'E' | b'D') => &text[1..],
        Some(b'+' | b'-') => text,
        Some(_) => return None,
    };
    let unsigned = signed.strip_prefix(['+', '-']).unwrap_or(signed);
    if unsigned.is_empty() || !unsigned.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    signed.parse().ok()
}

fn invalid(message: String) -> FitsError {
    FitsError::InvalidFormat(message)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::hdu::HduTrait;
    use crate::fits::header::keywords::{set_keyword, Keyword};
    use crate::fits::io::reader::HduInfo;
    use std::io::Cursor;

    // Fields given as (TFORM, TBCOL) with `extra` cards after them.
    fn table(
        fields: &[(&str, i64)],
        row_bytes: i64,
        rows: i64,
        extra: &[Keyword],
    ) -> AsciiTableHdu {
        let mut h = Header::new();
        h.add_keyword(Keyword::string("XTENSION", "TABLE"));
        h.add_keyword(Keyword::integer("BITPIX", 8));
        h.add_keyword(Keyword::integer("NAXIS", 2));
        h.add_keyword(Keyword::integer("NAXIS1", row_bytes));
        h.add_keyword(Keyword::integer("NAXIS2", rows));
        h.add_keyword(Keyword::integer("TFIELDS", fields.len() as i64));
        for (i, (form, col)) in fields.iter().enumerate() {
            h.add_keyword(Keyword::string(format!("TFORM{}", i + 1), form.to_string()));
            h.add_keyword(Keyword::integer(format!("TBCOL{}", i + 1), *col));
        }
        for card in extra {
            h.add_keyword(card.clone());
        }
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: (row_bytes * rows) as usize,
        };
        AsciiTableHdu::new(h, info).unwrap()
    }

    fn real(text: &str, decimals: usize) -> Result<f64> {
        real_field(text.as_bytes(), decimals)
    }

    #[test]
    fn tforms_follow_table_15() {
        assert_eq!(parse_tform("A12").unwrap(), (AsciiForm::Char, 12));
        assert_eq!(parse_tform("I6").unwrap(), (AsciiForm::Integer, 6));
        assert_eq!(
            parse_tform("F8.3").unwrap(),
            (AsciiForm::Real { decimals: 3 }, 8)
        );
        assert_eq!(
            parse_tform("E15.7").unwrap(),
            (AsciiForm::Real { decimals: 7 }, 15)
        );
        assert_eq!(
            parse_tform("D25.17").unwrap(),
            (AsciiForm::Real { decimals: 17 }, 25)
        );
        assert_eq!(
            parse_tform("E10").unwrap(),
            (AsciiForm::Real { decimals: 0 }, 10)
        );
        for bad in [
            "", "A", "A0", "I4.2", "X4", "a4", "F.3", "F8.x", "I-4", "B4",
        ] {
            assert!(parse_tform(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn integers_allow_spaces_and_a_sign() {
        assert_eq!(integer_field(b"  42 ").unwrap(), 42);
        assert_eq!(integer_field(b"+0003").unwrap(), 3);
        assert_eq!(integer_field(b"   -7").unwrap(), -7);
        assert_eq!(integer_field(b"     ").unwrap(), 0);
        for bad in [
            &b"4 2"[..],
            b"1.0",
            b"--1",
            b"+",
            b"x",
            b"99999999999999999999",
        ] {
            assert!(integer_field(bad).is_err());
        }
    }

    #[test]
    fn reals_follow_the_fortran_input_rules() {
        assert_eq!(real("1.5E3", 1).unwrap(), 1500.0);
        assert_eq!(real("2.5D2", 1).unwrap(), 250.0);
        assert_eq!(real("1.5-3", 1).unwrap(), 0.0015);
        assert_eq!(real("1.5+3", 1).unwrap(), 1500.0);
        assert_eq!(real("  .5 ", 2).unwrap(), 0.5);
        assert_eq!(real("-1.", 2).unwrap(), -1.0);
        assert_eq!(real("     ", 2).unwrap(), 0.0);
        assert_eq!(real("0.1", 0).unwrap(), 0.1);
        let long = "1.23456789012345678";
        assert_eq!(real(long, 0).unwrap(), long.parse::<f64>().unwrap());
    }

    #[test]
    fn reals_without_a_point_take_an_implicit_one() {
        assert_eq!(real("12345", 2).unwrap(), 123.45);
        assert_eq!(real("    1", 1).unwrap(), 0.1);
        assert_eq!(real("1E+02", 1).unwrap(), 10.0);
        assert_eq!(real("-5", 0).unwrap(), -5.0);
    }

    #[test]
    fn bad_reals_are_errors() {
        for bad in [
            "1 2", ".", "+", "1.2.3", "1E", "1E+", "1Q3", "abc", "1.5E3x", "NaN",
        ] {
            assert!(real(bad, 1).is_err(), "{bad}");
        }
    }

    // FITS 4.0 §7.2.5 names only the upper-case exponent letters.
    #[test]
    fn lower_case_exponents_are_errors() {
        for bad in ["1.5e3", "2.5d2", "1e-2"] {
            assert!(real(bad, 1).is_err(), "{bad}");
        }
    }

    #[test]
    fn fields_may_overlap_and_leave_gaps() {
        let t = table(&[("A3", 2), ("I2", 4), ("F3.1", 7)], 10, 1, &[]);
        let row = t
            .get_row(&mut Cursor::new(b"#ab12#123#".to_vec()), 0)
            .unwrap();
        assert_eq!(
            row,
            vec![
                TableValue::String("ab1".to_string()),
                TableValue::I64(12),
                TableValue::F64(12.3),
            ]
        );
        assert_eq!(t.hdu_type(), HduType::AsciiTable);
        assert_eq!(
            t.read_column_raw(&mut Cursor::new(b"#ab12#123#".to_vec()), 2)
                .unwrap(),
            b"123"
        );
    }

    #[test]
    fn scaled_fields_give_f64() {
        let extra = [Keyword::real("TSCAL1", 0.5), Keyword::real("TZERO2", 1.0)];
        let t = table(&[("I2", 1), ("E4.1", 3)], 6, 1, &extra);
        let row = t.get_row(&mut Cursor::new(b" 3 2.5".to_vec()), 0).unwrap();
        assert_eq!(row, vec![TableValue::F64(1.5), TableValue::F64(3.5)]);
    }

    #[test]
    fn bad_layouts_are_rejected() {
        // A card in `extra` replaces the card of the same name.
        let h = |fields: &[(&str, i64)], extra: &[Keyword]| {
            let mut cards = vec![
                Keyword::integer("NAXIS1", 8),
                Keyword::integer("NAXIS2", 1),
                Keyword::integer("TFIELDS", fields.len() as i64),
            ];
            for (i, (form, col)) in fields.iter().enumerate() {
                cards.push(Keyword::string(format!("TFORM{}", i + 1), form.to_string()));
                cards.push(Keyword::integer(format!("TBCOL{}", i + 1), *col));
            }
            for card in extra {
                set_keyword(&mut cards, card.clone());
            }
            let mut h = Header::new();
            for card in cards {
                h.add_keyword(card);
            }
            h
        };
        let info = HduInfo {
            index: 1,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 8,
        };
        for bad in [
            h(&[("I4", 6)], &[]),
            h(&[("I4", 0)], &[]),
            h(&[("I4", -1)], &[]),
            h(&[("Q4", 1)], &[]),
            h(&[("I4", 1)], &[Keyword::integer("TFIELDS", 2)]),
        ] {
            assert!(AsciiTableHdu::new(bad, info.clone()).is_err());
        }
        assert!(AsciiTableHdu::new(h(&[("I4", 5)], &[]), info).is_ok());
    }

    #[test]
    fn bad_field_text_is_an_error() {
        let t = table(&[("I3", 1)], 3, 1, &[]);
        assert!(t.get_row(&mut Cursor::new(b"1x2".to_vec()), 0).is_err());
        assert!(t
            .get_row(&mut Cursor::new(vec![0xff, b'1', b'2']), 0)
            .is_err());
    }
}
