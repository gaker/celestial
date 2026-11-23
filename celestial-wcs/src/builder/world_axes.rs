use crate::error::WcsResult;
use crate::header::KeywordProvider;

use super::parse::{is_celestial_pair, parse_ctype};

// The first index of these keywords numbers a world axis.
const WORLD_AXIS_KEYWORDS: [&str; 8] = [
    "CTYPE", "CRVAL", "CDELT", "CUNIT", "CROTA", "CD", "PC", "PV",
];

// Paper I numbers the world axes in header order, so a header may give
// latitude first. Reading through this view puts longitude first while
// CRPIXj and the second index of CDi_j and PCi_j keep numbering pixel axes.
pub(super) struct LongitudeFirst<'a, H> {
    header: &'a H,
    swapped: bool,
}

impl<'a, H: KeywordProvider> LongitudeFirst<'a, H> {
    pub(super) fn new(header: &'a H) -> WcsResult<Self> {
        let ctype1 = header.require_string("CTYPE1")?;
        let ctype2 = header.require_string("CTYPE2")?;
        let swapped = is_celestial_pair(parse_ctype(&ctype2)?.0, parse_ctype(&ctype1)?.0);
        Ok(Self { header, swapped })
    }

    // The pixel axis that shares its header number with each world axis.
    pub(super) fn paired_pixel_axes(&self) -> [usize; 2] {
        if self.swapped {
            [1, 0]
        } else {
            [0, 1]
        }
    }

    fn key(&self, key: &str) -> String {
        if self.swapped {
            swapped_key(key)
        } else {
            key.to_string()
        }
    }
}

impl<H: KeywordProvider> KeywordProvider for LongitudeFirst<'_, H> {
    fn get_string(&self, key: &str) -> Option<String> {
        self.header.get_string(&self.key(key))
    }

    fn get_float(&self, key: &str) -> Option<f64> {
        self.header.get_float(&self.key(key))
    }

    fn get_int(&self, key: &str) -> Option<i64> {
        self.header.get_int(&self.key(key))
    }
}

fn swapped_key(key: &str) -> String {
    let Some(at) = key.find(|c: char| c.is_ascii_digit()) else {
        return key.to_string();
    };
    let (name, rest) = key.split_at(at);
    let (axis, tail) = rest.split_at(1);
    let indexed = WORLD_AXIS_KEYWORDS.contains(&name) && (tail.is_empty() || tail.starts_with('_'));
    match axis {
        "1" if indexed => format!("{name}2{tail}"),
        "2" if indexed => format!("{name}1{tail}"),
        _ => key.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_swapped_key_renumbers_only_world_axes() {
        let cases = [
            ("CTYPE1", "CTYPE2"),
            ("CRVAL2", "CRVAL1"),
            ("CROTA1", "CROTA2"),
            ("CD1_2", "CD2_2"),
            ("PC2_1", "PC1_1"),
            ("PV1_29", "PV2_29"),
            ("CRPIX1", "CRPIX1"),
            ("A_2_0", "A_2_0"),
            ("LONPOLE", "LONPOLE"),
            ("CRVAL12", "CRVAL12"),
        ];
        for (key, expected) in cases {
            assert_eq!(swapped_key(key), expected, "{key}");
        }
    }
}
