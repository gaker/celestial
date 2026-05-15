use std::collections::HashMap;

use crate::error::{WcsError, WcsResult};
use crate::header::KeywordProvider;
use crate::spherical::Projection;

use super::build::MatrixSpec;


pub(super) fn create_projection_from_code(
    code: &str,
    pv_params: &HashMap<(u8, u8), f64>,
) -> WcsResult<Projection> {
    match code {
        "TAN" => Ok(Projection::tan()),
        "SIN" => {
            let xi = pv_params.get(&(2, 1)).copied().unwrap_or(0.0);
            let eta = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            if xi == 0.0 && eta == 0.0 {
                Ok(Projection::sin())
            } else {
                Ok(Projection::sin_with_params(xi, eta))
            }
        }
        "ARC" => Ok(Projection::arc()),
        "STG" => Ok(Projection::stg()),
        "ZEA" => Ok(Projection::zea()),
        "AZP" => {
            let mu = pv_params.get(&(2, 1)).copied().unwrap_or(0.0);
            let gamma = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            Ok(Projection::azp(mu, gamma))
        }
        "SZP" => {
            let mu = pv_params.get(&(2, 1)).copied().unwrap_or(0.0);
            let phi_c = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            let theta_c = pv_params.get(&(2, 3)).copied().unwrap_or(90.0);
            Ok(Projection::szp(mu, phi_c, theta_c))
        }
        "ZPN" => {
            let mut coeffs = Vec::new();
            for i in 0..=20 {
                if let Some(&val) = pv_params.get(&(2, i)) {
                    while coeffs.len() < i as usize {
                        coeffs.push(0.0);
                    }
                    coeffs.push(val);
                }
            }
            if coeffs.is_empty() {
                coeffs.push(0.0);
                coeffs.push(1.0);
            }
            Ok(Projection::zpn(coeffs))
        }
        "AIR" => {
            let theta_b = pv_params.get(&(2, 1)).copied().unwrap_or(90.0);
            Ok(Projection::air(theta_b))
        }
        "CAR" => Ok(Projection::car()),
        "MER" => Ok(Projection::mer()),
        "CEA" => {
            let lambda = pv_params.get(&(2, 1)).copied().unwrap_or(1.0);
            Ok(Projection::cea_with_lambda(lambda))
        }
        "CYP" => {
            let mu = pv_params.get(&(2, 1)).copied().unwrap_or(0.0);
            let lambda = pv_params.get(&(2, 2)).copied().unwrap_or(1.0);
            Ok(Projection::cyp(mu, lambda))
        }
        "SFL" => Ok(Projection::sfl()),
        "PAR" => Ok(Projection::par()),
        "MOL" => Ok(Projection::mol()),
        "AIT" => Ok(Projection::ait()),
        "COP" => {
            let theta_a = pv_params.get(&(2, 1)).copied().ok_or_else(|| {
                WcsError::missing_keyword("COP projection requires PV2_1 (theta_a)")
            })?;
            let eta = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            Ok(Projection::cop(theta_a, eta))
        }
        "COE" => {
            let theta_a = pv_params.get(&(2, 1)).copied().ok_or_else(|| {
                WcsError::missing_keyword("COE projection requires PV2_1 (theta_a)")
            })?;
            let eta = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            Ok(Projection::coe(theta_a, eta))
        }
        "COD" => {
            let theta_a = pv_params.get(&(2, 1)).copied().ok_or_else(|| {
                WcsError::missing_keyword("COD projection requires PV2_1 (theta_a)")
            })?;
            let eta = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            Ok(Projection::cod(theta_a, eta))
        }
        "COO" => {
            let theta_a = pv_params.get(&(2, 1)).copied().ok_or_else(|| {
                WcsError::missing_keyword("COO projection requires PV2_1 (theta_a)")
            })?;
            let eta = pv_params.get(&(2, 2)).copied().unwrap_or(0.0);
            Ok(Projection::coo(theta_a, eta))
        }
        "BON" => {
            let theta_1 = pv_params.get(&(2, 1)).ok_or_else(|| {
                WcsError::missing_keyword("BON projection requires PV2_1 (theta_1)")
            })?;
            Ok(Projection::bon(*theta_1))
        }
        "PCO" => Ok(Projection::pco()),
        "TSC" => Ok(Projection::tsc()),
        "CSC" => Ok(Projection::csc()),
        "QSC" => Ok(Projection::qsc()),
        _ => Err(WcsError::unsupported_projection(code)),
    }
}

pub(super) fn parse_ctype(ctype: &str) -> WcsResult<(&str, &str)> {
    let trimmed = ctype.trim();

    if let Some(dash_pos) = trimmed.rfind('-') {
        if dash_pos == 0 {
            return Err(WcsError::invalid_keyword(
                "CTYPE",
                format!("Invalid CTYPE format: '{}'", ctype),
            ));
        }

        let prefix_part = &trimmed[..dash_pos];
        let proj_part = &trimmed[dash_pos + 1..];

        let prefix = prefix_part.trim_end_matches('-');

        if proj_part.is_empty() {
            return Err(WcsError::invalid_keyword(
                "CTYPE",
                format!("Missing projection code in CTYPE: '{}'", ctype),
            ));
        }

        Ok((prefix, proj_part))
    } else {
        Err(WcsError::invalid_keyword(
            "CTYPE",
            format!("Invalid CTYPE format (no dash separator): '{}'", ctype),
        ))
    }
}

pub(super) fn parse_matrix(header: &impl KeywordProvider) -> WcsResult<MatrixSpec> {
    let cd11 = header.get_float("CD1_1");
    let cd12 = header.get_float("CD1_2");
    let cd21 = header.get_float("CD2_1");
    let cd22 = header.get_float("CD2_2");

    if cd11.is_some() || cd12.is_some() || cd21.is_some() || cd22.is_some() {
        let cd = [
            [cd11.unwrap_or(0.0), cd12.unwrap_or(0.0)],
            [cd21.unwrap_or(0.0), cd22.unwrap_or(0.0)],
        ];
        return Ok(MatrixSpec::Cd(cd));
    }

    let cdelt1 = header.get_float("CDELT1");
    let cdelt2 = header.get_float("CDELT2");

    if let (Some(c1), Some(c2)) = (cdelt1, cdelt2) {
        let pc11 = header.get_float("PC1_1").unwrap_or(1.0);
        let pc12 = header.get_float("PC1_2").unwrap_or(0.0);
        let pc21 = header.get_float("PC2_1").unwrap_or(0.0);
        let pc22 = header.get_float("PC2_2").unwrap_or(1.0);

        let pc = [[pc11, pc12], [pc21, pc22]];
        let cdelt = [c1, c2];

        return Ok(MatrixSpec::PcCdelt { pc, cdelt });
    }

    Err(WcsError::missing_keyword(
        "CD1_1 or CDELT1 (no transformation matrix found)",
    ))
}

pub(super) fn parse_pv_params(header: &impl KeywordProvider) -> HashMap<(u8, u8), f64> {
    let mut pv_params = HashMap::new();

    for axis in 1..=2u8 {
        for index in 0..=20u8 {
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
    use super::super::wcs::projection_to_code;

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
            let (prefix, proj) = parse_ctype(input).unwrap_or_else(|e| {
                panic!("parse_ctype({:?}) failed: {}", input, e)
            });
            assert_eq!(prefix, expected_prefix, "prefix for {:?}", input);
            assert_eq!(proj, expected_proj, "proj for {:?}", input);
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
    fn test_create_projection_from_code_round_trips_all_codes() {
        type PvEntry = ((u8, u8), f64);
        let with_params: &[(&str, &[PvEntry])] = &[
            ("SIN", &[]),
            ("ARC", &[]),
            ("STG", &[]),
            ("ZEA", &[]),
            ("AZP", &[((2, 1), 2.0), ((2, 2), 30.0)]),
            ("AZP", &[]),
            ("SZP", &[((2, 1), 2.0), ((2, 2), 45.0), ((2, 3), 60.0)]),
            ("SZP", &[]),
            ("ZPN", &[((2, 0), 0.0), ((2, 1), 1.0), ((2, 3), 0.1)]),
            ("ZPN", &[]),
            ("AIR", &[((2, 1), 45.0)]),
            ("AIR", &[]),
            ("CEA", &[((2, 1), 0.5)]),
            ("CEA", &[]),
            ("CYP", &[((2, 1), 1.0), ((2, 2), 2.0)]),
            ("COE", &[((2, 1), 45.0)]),
            ("COD", &[((2, 1), 30.0)]),
            ("COO", &[((2, 1), 60.0)]),
            ("BON", &[((2, 1), 45.0)]),
        ];
        for (code, params) in with_params {
            let pv: HashMap<(u8, u8), f64> = params.iter().copied().collect();
            let proj = create_projection_from_code(code, &pv)
                .unwrap_or_else(|e| panic!("{} failed to create: {}", code, e));
            assert_eq!(projection_to_code(&proj), *code, "code mismatch for {}", code);
        }
    }

    #[test]
    fn test_create_projection_missing_required_pv_param_errors() {
        for code in &["COE", "COD", "COO", "BON"] {
            let pv: HashMap<(u8, u8), f64> = HashMap::new();
            let err = create_projection_from_code(code, &pv).unwrap_err();
            assert!(
                err.to_string().contains("PV2_1"),
                "{} should report missing PV2_1, got: {}",
                code,
                err,
            );
        }
    }

}
