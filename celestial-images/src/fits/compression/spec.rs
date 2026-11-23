use super::quantize::Method;
use super::rice::RiceParams;
use crate::core::types::BitPix;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::keywords::KeywordValue;
use crate::fits::header::parser::Header;

// FITS 4.0 §10.1.1 limits ZNAXIS to 999, the same as NAXIS.
const MAX_AXES: usize = 999;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Codec {
    Rice(RiceParams),
    Gzip1,
    Gzip2,
    NoCompress,
}

// The ZQUANTIZ keyword. `Lossless` is 'NONE', which fpack writes for floats
// compressed without quantization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Quantization {
    Lossless,
    Quantized(Method),
}

// The keywords of FITS 4.0 §10.1.1 and §10.2 that describe how the tiles of
// one compressed image are stored.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Spec {
    pub(super) codec: Codec,
    pub(super) bitpix: BitPix,
    pub(super) axes: Vec<usize>,
    pub(super) tile: Vec<usize>,
    pub(super) quantization: Quantization,
    pub(super) dither_seed: Option<i64>,
    pub(super) blank: Option<i64>,
}

impl Spec {
    // Unlike scaling and WCS cards, a repeated compression keyword takes its
    // first card, which is the one cfitsio and astropy decompress with.
    pub(super) fn from_header(header: &Header) -> Result<Self> {
        let axes = axes(header)?;
        Ok(Self {
            codec: codec(header)?,
            bitpix: bitpix(header)?,
            tile: tile(header, &axes)?,
            axes,
            quantization: quantization(header)?,
            dither_seed: optional_integer(header, "ZDITHER0")?,
            blank: optional_integer(header, "ZBLANK")?,
        })
    }
}

// True unless ZNAXIS or one of the ZNAXISn is 0, the compressed form of an
// HDU without pixels. Other bad values are left for `Spec` to report.
pub(crate) fn has_pixels(header: &Header) -> bool {
    match header.get("ZNAXIS").as_i64() {
        Some(0) => false,
        Some(n) => (1..=n.min(MAX_AXES as i64))
            .all(|i| header.get(&format!("ZNAXIS{i}")).as_i64() != Some(0)),
        None => true,
    }
}

// fpack and astropy accept these names, and both reject other spellings such
// as 'RICE' or 'rice_1'.
fn codec(header: &Header) -> Result<Codec> {
    let name = required(header, "ZCMPTYPE")?
        .as_string()
        .ok_or_else(|| bad_value(header, "ZCMPTYPE"))?
        .trim_end();
    match name {
        "RICE_1" | "RICE_ONE" => Ok(Codec::Rice(rice_params(header)?)),
        "GZIP_1" => Ok(Codec::Gzip1),
        "GZIP_2" => Ok(Codec::Gzip2),
        "NOCOMPRESS" => Ok(Codec::NoCompress),
        "HCOMPRESS_1" | "PLIO_1" => Err(FitsError::UnsupportedCompression(format!(
            "{name} tiles cannot be read"
        ))),
        _ => Err(FitsError::UnsupportedCompression(format!(
            "unknown ZCMPTYPE '{name}'"
        ))),
    }
}

// FITS 4.0 Table 37 names each parameter in ZNAMEi and gives its value in
// ZVALi, so the parameters are matched by name, not by position.
fn rice_params(header: &Header) -> Result<RiceParams> {
    let mut params = RiceParams::default();
    for i in 1..=MAX_AXES {
        let (name_key, value_key) = (format!("ZNAME{i}"), format!("ZVAL{i}"));
        let name = header.get_keyword(&name_key);
        let value = header.get_keyword(&value_key);
        let (name, value) = match (name, value) {
            (None, None) => break,
            (Some(n), Some(v)) => (n.value.as_ref(), v.value.as_ref()),
            _ => {
                return Err(invalid(format!(
                    "{name_key} and {value_key} must come as a pair"
                )))
            }
        };
        let name = name
            .and_then(|v| v.as_string())
            .ok_or_else(|| bad_value(header, &name_key))?;
        let value = value
            .and_then(|v| v.as_integer())
            .ok_or_else(|| bad_value(header, &value_key));
        match name.trim_end() {
            "BLOCKSIZE" => params.block_size = block_size(value?, &value_key)?,
            "BYTEPIX" => params.bytepix = bytepix(value?, &value_key)?,
            _ => {}
        }
    }
    Ok(params)
}

fn block_size(value: i64, key: &str) -> Result<usize> {
    match value {
        16 | 32 => Ok(value as usize),
        _ => Err(value_error(key, value)),
    }
}

// BYTEPIX 8 is a valid value that this reader does not decode, so it is
// reported as unsupported by the Rice decoder rather than as invalid here.
fn bytepix(value: i64, key: &str) -> Result<usize> {
    match value {
        1 | 2 | 4 | 8 => Ok(value as usize),
        _ => Err(value_error(key, value)),
    }
}

fn bitpix(header: &Header) -> Result<BitPix> {
    let value = required_integer(header, "ZBITPIX")?;
    i32::try_from(value)
        .ok()
        .and_then(BitPix::from_value)
        .ok_or(FitsError::InvalidBitPix(
            value.clamp(i32::MIN as i64, i32::MAX as i64) as i32,
        ))
}

fn axes(header: &Header) -> Result<Vec<usize>> {
    let count = required_integer(header, "ZNAXIS")?;
    if !(1..=MAX_AXES as i64).contains(&count) {
        return Err(value_error("ZNAXIS", count));
    }
    (1..=count)
        .map(|i| positive(header, &format!("ZNAXIS{i}"), None))
        .collect()
}

// FITS 4.0 §10.1.1: ZTILE1 defaults to ZNAXIS1 and the other ZTILEn to 1,
// so the default tiles are the image rows.
fn tile(header: &Header, axes: &[usize]) -> Result<Vec<usize>> {
    axes.iter()
        .enumerate()
        .map(|(i, &axis)| {
            let default = if i == 0 { axis } else { 1 };
            positive(header, &format!("ZTILE{}", i + 1), Some(default))
        })
        .collect()
}

fn quantization(header: &Header) -> Result<Quantization> {
    let Some(keyword) = header.get_keyword("ZQUANTIZ") else {
        return Ok(Quantization::Quantized(Method::NoDither));
    };
    let name = keyword
        .value
        .as_ref()
        .and_then(|v| v.as_string())
        .ok_or_else(|| bad_value(header, "ZQUANTIZ"))?;
    match name.trim_end() {
        "NO_DITHER" => Ok(Quantization::Quantized(Method::NoDither)),
        "SUBTRACTIVE_DITHER_1" => Ok(Quantization::Quantized(Method::SubtractiveDither1)),
        "SUBTRACTIVE_DITHER_2" => Ok(Quantization::Quantized(Method::SubtractiveDither2)),
        "NONE" => Ok(Quantization::Lossless),
        other => Err(FitsError::UnsupportedCompression(format!(
            "ZQUANTIZ = '{other}'"
        ))),
    }
}

fn positive(header: &Header, key: &str, default: Option<usize>) -> Result<usize> {
    let value = match (optional_integer(header, key)?, default) {
        (Some(v), _) => v,
        (None, Some(d)) => return Ok(d),
        (None, None) => required_integer(header, key)?,
    };
    usize::try_from(value)
        .ok()
        .filter(|&n| n > 0)
        .ok_or_else(|| value_error(key, value))
}

fn required<'a>(header: &'a Header, key: &str) -> Result<&'a KeywordValue> {
    header
        .get_keyword(key)
        .and_then(|k| k.value.as_ref())
        .ok_or_else(|| FitsError::KeywordNotFound {
            keyword: key.to_string(),
        })
}

fn required_integer(header: &Header, key: &str) -> Result<i64> {
    required(header, key)?
        .as_integer()
        .ok_or_else(|| bad_value(header, key))
}

fn optional_integer(header: &Header, key: &str) -> Result<Option<i64>> {
    match header.get_keyword(key) {
        None => Ok(None),
        Some(k) => k
            .value
            .as_ref()
            .and_then(|v| v.as_integer())
            .map(Some)
            .ok_or_else(|| bad_value(header, key)),
    }
}

fn bad_value(header: &Header, key: &str) -> FitsError {
    let value = header.get(key).value().map(|v| v.to_string());
    FitsError::InvalidKeywordValue {
        keyword: key.to_string(),
        value: value.unwrap_or_default(),
    }
}

fn value_error(key: &str, value: i64) -> FitsError {
    FitsError::InvalidKeywordValue {
        keyword: key.to_string(),
        value: value.to_string(),
    }
}

fn invalid(message: String) -> FitsError {
    FitsError::InvalidFormat(message)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::Keyword;

    fn header(cards: &[Keyword]) -> Header {
        let mut h = Header::new();
        for k in cards {
            h.add_keyword(k.clone());
        }
        h
    }

    fn base() -> Vec<Keyword> {
        vec![
            Keyword::string("ZCMPTYPE", "RICE_1  "),
            Keyword::integer("ZBITPIX", 16),
            Keyword::integer("ZNAXIS", 2),
            Keyword::integer("ZNAXIS1", 23),
            Keyword::integer("ZNAXIS2", 17),
        ]
    }

    fn with(extra: &[Keyword]) -> Header {
        let mut cards = base();
        cards.extend_from_slice(extra);
        header(&cards)
    }

    #[test]
    fn defaults_follow_the_standard() {
        let spec = Spec::from_header(&with(&[])).unwrap();
        assert_eq!(spec.codec, Codec::Rice(RiceParams::default()));
        assert_eq!(spec.bitpix, BitPix::I16);
        assert_eq!(spec.axes, vec![23, 17]);
        assert_eq!(spec.tile, vec![23, 1]);
        assert_eq!(spec.quantization, Quantization::Quantized(Method::NoDither));
        assert_eq!(spec.dither_seed, None);
        assert_eq!(spec.blank, None);
    }

    #[test]
    fn rice_parameters_are_matched_by_name() {
        let h = with(&[
            Keyword::string("ZNAME1", "BYTEPIX"),
            Keyword::integer("ZVAL1", 2),
            Keyword::string("ZNAME2", "OTHER"),
            Keyword::integer("ZVAL2", 99),
            Keyword::string("ZNAME3", "BLOCKSIZE"),
            Keyword::integer("ZVAL3", 16),
        ]);
        let params = RiceParams {
            block_size: 16,
            bytepix: 2,
        };
        assert_eq!(Spec::from_header(&h).unwrap().codec, Codec::Rice(params));
    }

    #[test]
    fn rice_parameters_are_checked() {
        for (name, value) in [("BLOCKSIZE", 64), ("BYTEPIX", 3), ("BLOCKSIZE", 0)] {
            let h = with(&[
                Keyword::string("ZNAME1", name),
                Keyword::integer("ZVAL1", value),
            ]);
            assert!(Spec::from_header(&h).is_err(), "{name} {value}");
        }
        let unpaired = with(&[Keyword::string("ZNAME1", "BYTEPIX")]);
        assert!(Spec::from_header(&unpaired).is_err());
    }

    #[test]
    fn codec_names() {
        for (name, codec) in [
            ("RICE_ONE", Codec::Rice(RiceParams::default())),
            ("GZIP_1", Codec::Gzip1),
            ("GZIP_2  ", Codec::Gzip2),
            ("NOCOMPRESS", Codec::NoCompress),
        ] {
            let mut cards = base();
            cards[0] = Keyword::string("ZCMPTYPE", name);
            assert_eq!(Spec::from_header(&header(&cards)).unwrap().codec, codec);
        }
        for name in ["HCOMPRESS_1", "PLIO_1", "RICE", "gzip_1"] {
            let mut cards = base();
            cards[0] = Keyword::string("ZCMPTYPE", name);
            assert!(matches!(
                Spec::from_header(&header(&cards)),
                Err(FitsError::UnsupportedCompression(_))
            ));
        }
    }

    #[test]
    fn quantization_names() {
        for (name, q) in [
            ("NONE", Quantization::Lossless),
            (
                "SUBTRACTIVE_DITHER_2",
                Quantization::Quantized(Method::SubtractiveDither2),
            ),
        ] {
            let h = with(&[Keyword::string("ZQUANTIZ", name)]);
            assert_eq!(Spec::from_header(&h).unwrap().quantization, q);
        }
        let h = with(&[Keyword::string("ZQUANTIZ", "DITHER")]);
        assert!(matches!(
            Spec::from_header(&h),
            Err(FitsError::UnsupportedCompression(_))
        ));
    }

    #[test]
    fn bad_structure_keywords_are_errors() {
        let replace = |i: usize, k: Keyword| {
            let mut cards = base();
            cards[i] = k;
            Spec::from_header(&header(&cards))
        };
        assert!(matches!(
            replace(1, Keyword::integer("ZBITPIX", 12)),
            Err(FitsError::InvalidBitPix(12))
        ));
        assert!(replace(2, Keyword::integer("ZNAXIS", 0)).is_err());
        assert!(replace(2, Keyword::integer("ZNAXIS", 1000)).is_err());
        assert!(replace(3, Keyword::integer("ZNAXIS1", 0)).is_err());
        assert!(replace(3, Keyword::real("ZNAXIS1", 2.5)).is_err());
        assert!(Spec::from_header(&with(&[Keyword::integer("ZTILE2", 0)])).is_err());
        assert!(Spec::from_header(&with(&[Keyword::string("ZDITHER0", "1")])).is_err());
    }

    #[test]
    fn repeated_compression_keywords_use_the_first_card() {
        let spec = Spec::from_header(&with(&[
            Keyword::integer("ZTILE1", 5),
            Keyword::integer("ZTILE1", 7),
            Keyword::integer("ZDITHER0", 1),
            Keyword::integer("ZDITHER0", 2),
        ]))
        .unwrap();
        assert_eq!(spec.tile, vec![5, 1]);
        assert_eq!(spec.dither_seed, Some(1));
    }

    #[test]
    fn has_pixels_detects_empty_images() {
        assert!(has_pixels(&with(&[])));
        let mut cards = base();
        cards[4] = Keyword::integer("ZNAXIS2", 0);
        assert!(!has_pixels(&header(&cards)));
        cards[2] = Keyword::integer("ZNAXIS", 0);
        assert!(!has_pixels(&header(&cards)));
    }
}
