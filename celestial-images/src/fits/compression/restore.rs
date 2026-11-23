use super::spec::Spec;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::keywords::{Keyword, KeywordValue};
use crate::fits::header::parser::Header;

// Keywords of FITS 4.0 §10.1.1 and §10.2 that describe the table or the
// compression, with no meaning for the restored image.
const COMPRESSION_KEYS: &[&str] = &[
    "ZIMAGE", "ZCMPTYPE", "ZBITPIX", "ZNAXIS", "ZMASKCMP", "ZQUANTIZ", "ZDITHER0", "ZSIMPLE",
    "ZTENSION", "ZEXTEND", "ZBLOCKED", "ZPCOUNT", "ZGCOUNT", "ZBLANK", "ZSCALE", "ZZERO",
];
const COMPRESSION_FAMILIES: &[&str] = &["ZNAXIS", "ZTILE", "ZNAME", "ZVAL"];
const TABLE_KEYS: &[&str] = &[
    "SIMPLE", "XTENSION", "BITPIX", "NAXIS", "PCOUNT", "GCOUNT", "TFIELDS", "THEAP", "EXTEND",
    "BLOCKED", "CHECKSUM", "DATASUM", "END",
];
const TABLE_FAMILIES: &[&str] = &[
    "NAXIS", "TTYPE", "TFORM", "TUNIT", "TNULL", "TSCAL", "TZERO", "TDISP", "TDIM",
];

// fpack and astropy name the table this way when the image had no EXTNAME.
const DEFAULT_EXTNAME: &str = "COMPRESSED_IMAGE";

// The header of the uncompressed image (FITS 4.0 §10.1.1): the mandatory
// cards rebuilt from the Z keywords, then every other card of the table
// except those that describe the table or the compression, then `blank` as
// a BLANK card.
pub(super) fn image_header(table: &Header, spec: &Spec, blank: Option<i64>) -> Result<Header> {
    let mut out = Header::new();
    let extension = is_extension(table)?;
    for keyword in mandatory_cards(table, spec, extension)? {
        out.add_keyword(keyword);
    }
    for keyword in table.iter().filter_map(other_card) {
        out.add_keyword(keyword);
    }
    if let Some(blank) = blank {
        out.add_keyword(Keyword::integer("BLANK", blank));
    }
    Ok(out)
}

// FITS 4.0 §4.4.1: SIMPLE or XTENSION, BITPIX, NAXIS and NAXISn, then the
// cards that close the mandatory set for a primary HDU or an extension.
fn mandatory_cards(table: &Header, spec: &Spec, extension: bool) -> Result<Vec<Keyword>> {
    let first = if extension {
        card(table, "ZTENSION", "XTENSION", "IMAGE".into())
    } else {
        card(table, "ZSIMPLE", "SIMPLE", true.into())
    };
    let bitpix = i64::from(spec.bitpix.value());
    let mut cards = vec![
        first,
        card(table, "ZBITPIX", "BITPIX", bitpix.into()),
        card(table, "ZNAXIS", "NAXIS", count(spec.axes.len())),
    ];
    cards.extend(spec.axes.iter().enumerate().map(|(i, &axis)| {
        let n = i + 1;
        card(
            table,
            &format!("ZNAXIS{n}"),
            &format!("NAXIS{n}"),
            count(axis),
        )
    }));
    if extension {
        cards.push(group_card(table, "ZPCOUNT", "PCOUNT", 0)?);
        cards.push(group_card(table, "ZGCOUNT", "GCOUNT", 1)?);
    } else {
        for (from, to) in [("ZEXTEND", "EXTEND"), ("ZBLOCKED", "BLOCKED")] {
            if let Some(value) = table.get(from).value() {
                cards.push(card(table, from, to, value.clone()));
            }
        }
    }
    Ok(cards)
}

// ZSIMPLE marks an image from a primary HDU and ZTENSION one from an
// extension, and §10.1.1 forbids both together. A table with neither is
// restored as an image extension.
fn is_extension(table: &Header) -> Result<bool> {
    let simple = table.get_keyword("ZSIMPLE");
    let tension = table.get_keyword("ZTENSION");
    match (simple, tension) {
        (Some(_), Some(_)) => Err(FitsError::InvalidFormat(
            "ZSIMPLE and ZTENSION are both present".to_string(),
        )),
        (Some(k), None) if k.value != Some(KeywordValue::Logical(true)) => {
            Err(bad_value("ZSIMPLE", &k.value))
        }
        (None, Some(k)) if !is_image(&k.value) => Err(bad_value("ZTENSION", &k.value)),
        (Some(_), None) => Ok(false),
        _ => Ok(true),
    }
}

fn is_image(value: &Option<KeywordValue>) -> bool {
    matches!(value, Some(KeywordValue::String(s)) if s.trim_end() == "IMAGE")
}

// A tile-compressed image has one group and no parameters, the only values
// an IMAGE extension allows (FITS 4.0 §7.1.1).
fn group_card(table: &Header, from: &str, to: &str, required: i64) -> Result<Keyword> {
    match table.get_keyword(from) {
        Some(k) if k.value != Some(KeywordValue::Integer(required)) => {
            Err(bad_value(from, &k.value))
        }
        _ => Ok(card(table, from, to, required.into())),
    }
}

// The value comes from `Spec`, which has checked it, and the comment from
// the Z keyword when there is one.
fn card(table: &Header, from: &str, to: &str, value: KeywordValue) -> Keyword {
    let keyword = Keyword::new(to.to_string()).with_value(value);
    match table.get_keyword(from).and_then(|k| k.comment.clone()) {
        Some(comment) => keyword.with_comment(comment),
        None => keyword,
    }
}

// Spec::from_header limits axis lengths to values read from i64 cards.
fn count(n: usize) -> KeywordValue {
    KeywordValue::Integer(n as i64)
}

// Names that the compressed-image table defines itself, so a copy from the
// image header would contradict the table.
pub(super) fn describes_table_or_compression(name: &str) -> bool {
    is_listed(name, COMPRESSION_KEYS, COMPRESSION_FAMILIES)
        || is_listed(name, TABLE_KEYS, TABLE_FAMILIES)
}

fn other_card(keyword: &Keyword) -> Option<Keyword> {
    let name = keyword.name.as_str();
    let renamed = match name {
        "ZHECKSUM" => "CHECKSUM",
        "ZDATASUM" => "DATASUM",
        _ if describes_table_or_compression(name) => return None,
        "EXTNAME" if is_default_extname(&keyword.value) => return None,
        _ => name,
    };
    let mut out = keyword.clone();
    out.name = renamed.to_string();
    Some(out)
}

fn is_listed(name: &str, keys: &[&str], families: &[&str]) -> bool {
    keys.contains(&name)
        || families.iter().any(|family| {
            name.strip_prefix(family)
                .is_some_and(|n| !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()))
        })
}

fn is_default_extname(value: &Option<KeywordValue>) -> bool {
    matches!(value, Some(KeywordValue::String(s)) if s.trim_end() == DEFAULT_EXTNAME)
}

fn bad_value(key: &str, value: &Option<KeywordValue>) -> FitsError {
    FitsError::InvalidKeywordValue {
        keyword: key.to_string(),
        value: value.as_ref().map(|v| v.to_string()).unwrap_or_default(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(extra: &[Keyword]) -> Header {
        let mut h = Header::new();
        for k in [
            Keyword::string("XTENSION", "BINTABLE"),
            Keyword::integer("BITPIX", 8),
            Keyword::integer("NAXIS", 2),
            Keyword::integer("NAXIS1", 8),
            Keyword::integer("NAXIS2", 3),
            Keyword::integer("PCOUNT", 99),
            Keyword::integer("GCOUNT", 1),
            Keyword::integer("TFIELDS", 1),
            Keyword::string("TTYPE1", "COMPRESSED_DATA"),
            Keyword::string("TFORM1", "1PB(9)"),
            Keyword::logical("ZIMAGE", true),
            Keyword::string("ZCMPTYPE", "RICE_1"),
            Keyword::string("ZNAME1", "BLOCKSIZE"),
            Keyword::integer("ZVAL1", 32),
            Keyword::integer("ZBITPIX", 16),
            Keyword::integer("ZNAXIS", 2),
            Keyword::integer("ZNAXIS1", 4).with_comment("width"),
            Keyword::integer("ZNAXIS2", 3),
            Keyword::integer("ZTILE1", 4),
            Keyword::string("CHECKSUM", "table"),
            Keyword::string("ZHECKSUM", "image"),
            Keyword::string("OBJECT", "M31"),
        ]
        .into_iter()
        .chain(extra.iter().cloned())
        {
            h.add_keyword(k);
        }
        h
    }

    fn restore(h: &Header) -> Result<Header> {
        image_header(h, &Spec::from_header(h)?, None)
    }

    fn names(h: &Header) -> Vec<&str> {
        h.iter().map(|k| k.name.as_str()).collect()
    }

    #[test]
    fn restores_a_primary_image() {
        let h = table(&[
            Keyword::logical("ZSIMPLE", true),
            Keyword::logical("ZEXTEND", true),
            Keyword::string("EXTNAME", "COMPRESSED_IMAGE"),
        ]);
        let out = restore(&h).unwrap();
        assert_eq!(
            names(&out),
            vec!["SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND", "CHECKSUM", "OBJECT"]
        );
        assert_eq!(out.get("BITPIX").as_i64(), Some(16));
        assert_eq!(out.get("NAXIS1").as_i64(), Some(4));
        assert_eq!(
            out.get_keyword("NAXIS1").unwrap().comment.as_deref(),
            Some("width")
        );
        assert_eq!(out.get("CHECKSUM").as_str(), Some("image"));
    }

    #[test]
    fn restores_an_image_extension() {
        let h = table(&[
            Keyword::string("ZTENSION", "IMAGE   "),
            Keyword::integer("ZPCOUNT", 0),
            Keyword::integer("ZGCOUNT", 1),
            Keyword::string("EXTNAME", "SCI"),
        ]);
        let out = restore(&h).unwrap();
        assert_eq!(
            names(&out),
            vec![
                "XTENSION", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "CHECKSUM",
                "OBJECT", "EXTNAME"
            ]
        );
        assert_eq!(out.get("XTENSION").as_str(), Some("IMAGE"));
        assert_eq!(out.get("PCOUNT").as_i64(), Some(0));
    }

    #[test]
    fn without_zsimple_or_ztension_restores_an_extension() {
        let out = restore(&table(&[])).unwrap();
        assert_eq!(out.get("XTENSION").as_str(), Some("IMAGE"));
        assert_eq!(out.get("GCOUNT").as_i64(), Some(1));
    }

    #[test]
    fn rejects_inconsistent_structure_keywords() {
        for extra in [
            vec![
                Keyword::logical("ZSIMPLE", true),
                Keyword::string("ZTENSION", "IMAGE"),
            ],
            vec![Keyword::logical("ZSIMPLE", false)],
            vec![Keyword::string("ZTENSION", "BINTABLE")],
            vec![
                Keyword::string("ZTENSION", "IMAGE"),
                Keyword::integer("ZPCOUNT", 4),
            ],
            vec![
                Keyword::string("ZTENSION", "IMAGE"),
                Keyword::integer("ZGCOUNT", 2),
            ],
        ] {
            assert!(restore(&table(&extra)).is_err(), "{extra:?}");
        }
    }

    #[test]
    fn keyword_families_need_digits() {
        assert!(is_listed("TTYPE12", &[], TABLE_FAMILIES));
        assert!(!is_listed("TTYPE", &[], TABLE_FAMILIES));
        assert!(!is_listed("TDIMX", &[], TABLE_FAMILIES));
        assert!(!is_listed("ZVALUE", &[], COMPRESSION_FAMILIES));
    }
}
