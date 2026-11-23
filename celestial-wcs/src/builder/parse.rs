use std::collections::HashMap;

use crate::error::{WcsError, WcsResult};
use crate::header::KeywordProvider;
use crate::spherical::Projection;
use crate::trig::sincosd;

use super::build::MatrixSpec;
use super::units::{in_degrees, AngleUnit};
use super::world_axes::LongitudeFirst;

// Paper II §3: the coordinate type, padding dashes, then the projection code.
// The SIP convention appends a further "-SIP".
pub(super) fn parse_ctype(ctype: &str) -> WcsResult<(&str, &str, Option<&str>)> {
    let trimmed = ctype.trim();
    let invalid = |why: &str| WcsError::invalid_keyword("CTYPE", format!("{why}: '{ctype}'"));
    if trimmed.starts_with('-') {
        return Err(invalid("Invalid CTYPE format"));
    }
    let Some((prefix, rest)) = trimmed.split_once('-') else {
        return Err(invalid("Invalid CTYPE format (no dash separator)"));
    };
    let (code, suffix) = match rest.trim_start_matches('-').split_once('-') {
        Some((code, suffix)) => (code, Some(suffix)),
        None => (rest.trim_start_matches('-'), None),
    };
    if code.is_empty() {
        return Err(invalid("Missing projection code in CTYPE"));
    }
    if suffix.is_some_and(|s| s.is_empty() || s.contains('-')) {
        return Err(invalid("Malformed CTYPE suffix"));
    }
    Ok((prefix, code, suffix))
}

// Both axes must name the same projection and the same suffix.
pub(super) fn parse_ctype_pair<'a>(
    ctype1: &'a str,
    ctype2: &'a str,
) -> WcsResult<(&'a str, &'a str, Option<&'a str>)> {
    let (prefix1, code1, suffix1) = parse_ctype(ctype1)?;
    let (prefix2, code2, suffix2) = parse_ctype(ctype2)?;
    let mismatch = |message: String| Err(WcsError::invalid_keyword("CTYPE1/CTYPE2", message));
    if (code1, suffix1) != (code2, suffix2) {
        return mismatch(format!(
            "Mismatched projection codes: '{ctype1}' vs '{ctype2}'"
        ));
    }
    if !is_celestial_pair(prefix1, prefix2) {
        return mismatch(format!(
            "'{ctype1}' and '{ctype2}' are not a celestial pair"
        ));
    }
    Ok((prefix1, code1, suffix1))
}

// Paper II §3: RA/DEC, or xLON/xLAT and yzLN/yzLT with the same x or yz.
pub(super) fn is_celestial_pair(lon: &str, lat: &str) -> bool {
    let same = |lon_ending, lat_ending| {
        celestial_system(lon, lon_ending)
            .is_some_and(|system| celestial_system(lat, lat_ending) == Some(system))
    };
    (lon, lat) == ("RA", "DEC") || same("LON", "LAT") || same("LN", "LT")
}

// The left half of CTYPE is four characters.
fn celestial_system<'a>(prefix: &'a str, ending: &str) -> Option<&'a str> {
    prefix
        .strip_suffix(ending)
        .filter(|system| system.len() + ending.len() == 4)
}

// Paper I §2.4: an absent CRPIXj or CRVALi is 0, an absent CDELTi 1.
pub(super) fn axis_pair(header: &impl KeywordProvider, key: &str, default: f64) -> [f64; 2] {
    [1, 2].map(|i| header.get_float(&format!("{key}{i}")).unwrap_or(default))
}

// CDi_j and CDELTi are in the units of axis i; PC is a pure number.
pub(super) fn parse_matrix(
    header: &LongitudeFirst<'_, impl KeywordProvider>,
    units: [AngleUnit; 2],
) -> MatrixSpec {
    if let Some(cd) = parse_square(header, "CD", [[0.0; 2]; 2]) {
        return MatrixSpec::Cd([0, 1].map(|i| cd[i].map(|v| units[i].to_degrees(v))));
    }
    let cdelt = in_degrees(axis_pair(header, "CDELT", 1.0), units);
    let pixel = header.paired_pixel_axes();
    let unit = [0, 1].map(|i| [0, 1].map(|j| if pixel[i] == j { 1.0 } else { 0.0 }));

    // Paper II §6.2 bars CROTA2 beside PC, so any PC keyword means the header
    // uses the current formalism.
    if let Some(pc) = parse_square(header, "PC", unit) {
        return MatrixSpec::PcCdelt { pc, cdelt };
    }
    match header.get_float("CROTA2") {
        Some(rho) => MatrixSpec::Cd(crota_cd(cdelt, rho, pixel)),
        None => MatrixSpec::PcCdelt { pc: unit, cdelt },
    }
}

// None when the header has no element of the matrix; otherwise absent
// elements take their defaults.
fn parse_square(
    header: &impl KeywordProvider,
    prefix: &str,
    default: [[f64; 2]; 2],
) -> Option<[[f64; 2]; 2]> {
    let element = |i: usize, j: usize| header.get_float(&format!("{prefix}{}_{}", i + 1, j + 1));
    let found = [0, 1].map(|i| [0, 1].map(|j| element(i, j)));
    if found.iter().flatten().all(Option::is_none) {
        return None;
    }
    Some([0, 1].map(|i| [0, 1].map(|j| found[i][j].unwrap_or(default[i][j]))))
}

// Paper II eq. 189, the AIPS rotation of the CDELT-scaled axes. It equals
// CDELT times the PC of eq. 187 without rounding the ratio λ of eq. 188.
// Its columns are the pixel axes paired with longitude and latitude.
fn crota_cd(cdelt: [f64; 2], rho_deg: f64, pixel: [usize; 2]) -> [[f64; 2]; 2] {
    let (sin_rho, cos_rho) = sincosd(rho_deg);
    let rotated = [
        [cdelt[0] * cos_rho, -cdelt[1] * sin_rho],
        [cdelt[0] * sin_rho, cdelt[1] * cos_rho],
    ];
    [0, 1].map(|i| [0, 1].map(|j| rotated[i][pixel[j]]))
}

pub(super) fn parse_pv_params(header: &impl KeywordProvider) -> HashMap<(u8, u8), f64> {
    let mut pv_params = HashMap::new();

    for axis in 1..=2u8 {
        for index in 0..=Projection::PV_LAST_INDEX {
            let key = format!("PV{}_{}", axis, index);
            if let Some(value) = header.get_float(&key) {
                pv_params.insert((axis, index), value);
            }
        }
    }

    pv_params
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_ctype_valid_forms() {
        // Sweep covers: short and long axis names, single- vs triple-dash
        // separators, surrounding whitespace, and trailing dashes in the
        // prefix that must be trimmed.
        let cases: &[(&str, &str, &str)] = &[
            ("RA---TAN", "RA", "TAN"),
            ("DEC--TAN", "DEC", "TAN"),
            ("GLON-SIN", "GLON", "SIN"),
            ("GLAT-SIN", "GLAT", "SIN"),
            ("RA---ZEA", "RA", "ZEA"),
            ("GLON-TAN", "GLON", "TAN"),
            ("  RA---TAN  ", "RA", "TAN"),
        ];
        for &(input, expected_prefix, expected_proj) in cases {
            let (prefix, proj, suffix) = parse_ctype(input)
                .unwrap_or_else(|e| panic!("parse_ctype({:?}) failed: {}", input, e));
            assert_eq!(prefix, expected_prefix, "prefix for {:?}", input);
            assert_eq!(proj, expected_proj, "proj for {:?}", input);
            assert_eq!(suffix, None, "suffix for {:?}", input);
        }
    }

    #[test]
    fn test_parse_ctype_separates_the_sip_suffix() {
        assert_eq!(
            parse_ctype("RA---TAN-SIP").unwrap(),
            ("RA", "TAN", Some("SIP"))
        );
        assert_eq!(
            parse_ctype("DEC--TAN-SIP").unwrap(),
            ("DEC", "TAN", Some("SIP"))
        );
        for malformed in ["RA---TAN-", "RA---TAN-SIP-X"] {
            let err = parse_ctype(malformed).unwrap_err().to_string();
            assert!(err.contains("Malformed CTYPE suffix"), "{malformed}: {err}");
        }
    }

    #[test]
    fn test_parse_ctype_invalid_forms() {
        // No dash separator at all.
        let err = parse_ctype("RATAN").unwrap_err().to_string();
        assert!(err.contains("no dash separator"), "got: {}", err);

        // Trailing dash, no projection after it.
        let err = parse_ctype("RA---").unwrap_err().to_string();
        assert!(err.contains("Missing projection"), "got: {}", err);

        // Leading dash means an empty prefix.
        let err = parse_ctype("-TAN").unwrap_err().to_string();
        assert!(err.contains("Invalid CTYPE"), "got: {}", err);
    }

    #[test]
    fn test_parse_pv_params_reads_indices_up_to_29() {
        let mut header = crate::header::KeywordMap::new();
        header.set_float("PV2_29", 1e-9).set_float("PV2_30", 1.0);
        let pv = parse_pv_params(&header);
        assert_eq!(pv, HashMap::from([((2, 29), 1e-9)]));
    }
}
