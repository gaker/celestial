use super::{MainTerm, ParseError, PertTerm};

const MAIN_WIDTH: usize = 99;
const PERT_WIDTH: usize = 93;
const PERT_HEADER_WORDS: [&str; 4] = ["PERTURBATIONS", "LONGITUDE", "LATITUDE", "DISTANCE"];

pub(super) fn parse_main_header(line: &str) -> Result<usize, ParseError> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 3 {
        return Err(ParseError::InvalidHeader(format!(
            "Main header too short: '{}'",
            line
        )));
    }
    parts[parts.len() - 1]
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid term count in: '{}'", line)))
}

pub(super) fn is_pert_header(line: &str) -> bool {
    PERT_HEADER_WORDS.iter().any(|word| line.contains(word))
}

pub(super) fn parse_pert_header(line: &str) -> Result<(usize, u8), ParseError> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 3 {
        return Err(ParseError::InvalidHeader(format!(
            "Pert header too short: '{}'",
            line
        )));
    }
    let term_count: usize = parts[parts.len() - 2]
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid term count in: '{}'", line)))?;
    let time_power: u8 = parts[parts.len() - 1]
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid time power in: '{}'", line)))?;
    Ok((term_count, time_power))
}

// The files are fixed-width, so a short line or a blank field is damage; reading
// it as zero would drop a term without a trace.
fn record(line: &str, width: usize) -> Result<&str, String> {
    if !line.is_ascii() || line.len() < width || !line[width..].trim().is_empty() {
        return Err(format!(
            "Line too short or too long ({} columns, expected {}): '{}'",
            line.len(),
            width,
            line
        ));
    }
    Ok(&line[..width])
}

fn int_fields<const N: usize>(line: &str, start: usize, name: &str) -> Result<[i32; N], String> {
    let mut out = [0; N];
    for (i, v) in out.iter_mut().enumerate() {
        let s = &line[start + i * 3..start + i * 3 + 3];
        *v = s
            .trim()
            .parse()
            .map_err(|_| format!("Invalid {}[{}]: '{}'", name, i, s))?;
    }
    Ok(out)
}

fn finite(s: &str) -> Option<f64> {
    s.trim().parse().ok().filter(|v: &f64| v.is_finite())
}

pub(super) fn parse_main_term(line: &str) -> Result<MainTerm, ParseError> {
    let line = record(line, MAIN_WIDTH).map_err(ParseError::InvalidMainTerm)?;
    let delaunay = int_fields(line, 0, "delaunay").map_err(ParseError::InvalidMainTerm)?;
    let mut coeffs = [0.0f64; 7];
    for (i, c) in coeffs.iter_mut().enumerate() {
        let start = if i == 0 { 14 } else { 15 + i * 12 };
        let end = if i == 0 { 27 } else { start + 12 };
        let s = &line[start..end];
        *c = finite(s)
            .ok_or_else(|| ParseError::InvalidMainTerm(format!("Invalid coeff[{}]: '{}'", i, s)))?;
    }
    Ok(MainTerm { delaunay, coeffs })
}

pub(super) fn parse_fortran_double(s: &str) -> Result<f64, ParseError> {
    let s = s.trim().replace('D', "E").replace('d', "e");
    finite(&s).ok_or_else(|| ParseError::InvalidPertTerm(format!("Invalid double: '{}'", s)))
}

pub(super) fn parse_pert_term(line: &str) -> Result<PertTerm, ParseError> {
    let line = record(line, PERT_WIDTH).map_err(ParseError::InvalidPertTerm)?;
    Ok(PertTerm {
        sin_coeff: parse_fortran_double(&line[5..25])?,
        cos_coeff: parse_fortran_double(&line[25..45])?,
        multipliers: int_fields(line, 45, "multiplier").map_err(ParseError::InvalidPertTerm)?,
    })
}
