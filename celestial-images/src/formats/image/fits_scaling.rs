use crate::core::errors::Result;
use crate::core::types::BitPix;
use crate::fits::data::scaling::{
    integer_value, Linear, I8_OFFSET, U16_OFFSET, U32_OFFSET, U64_OFFSET,
};
use crate::fits::header::keywords::{Keyword, KeywordValue};
use crate::fits::header::parser::Header;
use crate::fits::io::card::has_unwritable_name;
use crate::formats::pixel_data::PixelData;

// FITS 4.0 4.4.2.5: physical = BZERO + BSCALE * stored, and integer pixels
// equal to BLANK are undefined. `blank` is in stored units.
pub(super) struct Scaling {
    linear: Linear,
    pub(super) blank: Option<i64>,
}

impl Scaling {
    // A BLANK that is not an integer inside the stored type's range cannot
    // match any pixel, so it is ignored. BLANK never applies to float data.
    pub(super) fn from_header(header: &Header, bitpix: BitPix) -> Result<Self> {
        let blank = header
            .get_last("BLANK")
            .value()
            .and_then(integer_value)
            .filter(|&b| in_stored_range(b as i128, bitpix));
        Ok(Self {
            linear: Linear::from_header(header, "BSCALE", "BZERO")?,
            blank,
        })
    }

    pub(super) fn form(&self, bitpix: BitPix) -> Form {
        let offset = match bitpix {
            BitPix::U8 => Some(I8_OFFSET),
            BitPix::I16 => Some(U16_OFFSET),
            BitPix::I32 => Some(U32_OFFSET),
            BitPix::I64 => Some(U64_OFFSET),
            BitPix::F32 | BitPix::F64 => None,
        };
        if self.linear.is_identity() {
            Form::Plain
        } else if offset.is_some_and(|zero| self.linear.is_offset(zero)) {
            Form::Offset
        } else {
            Form::Physical
        }
    }

    pub(super) fn integer(&self, stored: i64) -> f64 {
        if Some(stored) == self.blank {
            f64::NAN
        } else {
            self.linear.apply(stored as f64)
        }
    }

    pub(super) fn real(&self, stored: f64) -> f64 {
        self.linear.apply(stored)
    }
}

// How stored integers map to in-memory pixels: unchanged, shifted by the
// Table 11 offset for their width, or scaled to physical f64 values.
pub(super) enum Form {
    Plain,
    Offset,
    Physical,
}

fn in_stored_range(stored: i128, bitpix: BitPix) -> bool {
    let (min, max) = match bitpix {
        BitPix::U8 => (0, u8::MAX as i128),
        BitPix::I16 => (i16::MIN as i128, i16::MAX as i128),
        BitPix::I32 => (i32::MIN as i128, i32::MAX as i128),
        BitPix::I64 => (i64::MIN as i128, i64::MAX as i128),
        BitPix::F32 | BitPix::F64 => return false,
    };
    (min..=max).contains(&stored)
}

pub(super) fn table11_bzero(pixels: &PixelData) -> Option<i128> {
    match pixels {
        PixelData::I8(_) => Some(I8_OFFSET),
        PixelData::U16(_) => Some(U16_OFFSET),
        PixelData::U32(_) => Some(U32_OFFSET),
        PixelData::U64(_) => Some(U64_OFFSET),
        _ => None,
    }
}

enum Card {
    Scaling,
    Blank,
    Other,
}

fn card(keyword: &Keyword) -> Card {
    let name = keyword.name.as_str();
    if name.eq_ignore_ascii_case("BSCALE") || name.eq_ignore_ascii_case("BZERO") {
        Card::Scaling
    } else if name.eq_ignore_ascii_case("BLANK") {
        Card::Blank
    } else {
        Card::Other
    }
}

fn with_integer(keyword: &Keyword, value: i64) -> Keyword {
    let mut out = keyword.clone();
    out.value = Some(KeywordValue::Integer(value));
    out
}

// In-memory keywords describe the in-memory pixels. BSCALE and BZERO are
// already applied, and BLANK holds the physical value of an undefined pixel.
// A u64 BLANK above i64::MAX has no integer card, so it is dropped.
pub(super) fn loaded_keywords(
    keywords: &[Keyword],
    pixels: &PixelData,
    stored_blank: Option<i64>,
) -> Vec<Keyword> {
    let physical = match pixels {
        PixelData::F32(_) | PixelData::F64(_) => None,
        _ => stored_blank
            .and_then(|b| i64::try_from(b as i128 + table11_bzero(pixels).unwrap_or(0)).ok()),
    };
    keywords
        .iter()
        .filter_map(|k| match card(k) {
            Card::Scaling => None,
            Card::Blank => physical.map(|v| with_integer(k, v)),
            Card::Other => Some(k.clone()),
        })
        .collect()
}

// The reverse of `loaded_keywords`: BLANK goes back to stored units, and
// I8, U16, U32 and U64 get their Table 11 BZERO. A BLANK that has no stored
// value for these pixels is dropped, the same as on load.
pub(super) fn keywords_for_save(keywords: &[Keyword], pixels: &PixelData) -> Vec<Keyword> {
    let bzero = table11_bzero(pixels);
    let stored_blank = |k: &Keyword| {
        let stored = k.value.as_ref().and_then(integer_value)? as i128 - bzero.unwrap_or(0);
        in_stored_range(stored, pixels.bitpix()).then_some(stored as i64)
    };
    // A keyword read from another format may have a name no FITS card holds.
    let mut out: Vec<Keyword> = keywords
        .iter()
        .filter(|k| !has_unwritable_name(k))
        .filter_map(|k| match card(k) {
            Card::Scaling => None,
            Card::Blank => stored_blank(k).map(|v| with_integer(k, v)),
            Card::Other => Some(k.clone()),
        })
        .collect();
    if let Some(bzero) = bzero {
        // Keyword integers are i64, so the u64 offset 2^63 is written as the
        // real 9223372036854776000.0, which parses back to exactly 2^63.
        let card = match i64::try_from(bzero) {
            Ok(v) => Keyword::integer("BZERO", v),
            Err(_) => Keyword::real("BZERO", bzero as f64),
        };
        out.push(card.with_comment("offset for the stored integer type"));
        out.push(Keyword::integer("BSCALE", 1).with_comment("default scaling factor"));
    }
    out
}
