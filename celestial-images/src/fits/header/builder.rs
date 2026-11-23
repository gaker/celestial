use super::keywords::Keyword;
use crate::core::types::BitPix;

// The reference to the FITS paper that follows the mandatory cards.
const STANDARD_COMMENTS: [&str; 2] = [
    "FITS (Flexible Image Transport System) format is defined in 'Astronomy",
    "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H",
];

// FITS 4.0 §4.4.1.1: SIMPLE, BITPIX, NAXIS and NAXISn open a primary
// header. EXTEND = T (§4.4.2.1) allows extensions to follow. Callers check
// that every axis length fits an i64.
pub(crate) fn primary_cards(bitpix: BitPix, dimensions: &[usize]) -> Vec<Keyword> {
    let mut cards = vec![
        Keyword::logical("SIMPLE", true).with_comment("file does conform to FITS standard"),
        Keyword::integer("BITPIX", i64::from(bitpix.value()))
            .with_comment("number of bits per data pixel"),
        Keyword::integer("NAXIS", dimensions.len() as i64).with_comment("number of data axes"),
    ];
    cards.extend(dimensions.iter().enumerate().map(|(i, &dim)| {
        Keyword::integer(format!("NAXIS{}", i + 1), dim as i64)
            .with_comment(format!("length of data axis {}", i + 1))
    }));
    cards
        .push(Keyword::logical("EXTEND", true).with_comment("FITS dataset may contain extensions"));
    cards
}

// `primary_cards` and the standard comments. A header read from a file this
// crate wrote already holds the comments, so they are added only when
// `keywords` lacks them, and repeated saves do not repeat them.
pub(crate) fn primary_header_cards(
    bitpix: BitPix,
    dimensions: &[usize],
    keywords: &[Keyword],
) -> Vec<Keyword> {
    let mut cards = primary_cards(bitpix, dimensions);
    for text in STANDARD_COMMENTS {
        let present = keywords
            .iter()
            .any(|k| k.name == "COMMENT" && k.comment.as_deref() == Some(text));
        if !present {
            cards.push(Keyword::comment(text));
        }
    }
    cards
}
