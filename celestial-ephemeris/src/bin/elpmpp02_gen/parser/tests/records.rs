use super::*;

#[test]
fn headers() {
    assert_eq!(
        parse_main_header("MAIN PROBLEM. LONGITUDE.       1023").unwrap(),
        1023
    );
    let pert = "PERTURBATIONS. LONGITUDE.     11314         0";
    assert_eq!(parse_pert_header(pert).unwrap(), (11314, 0));
    assert_eq!(parse_pert_header("LATITUDE. 52 2").unwrap(), (52, 2));
    for line in [pert, "LONGITUDE 2 0", "LATITUDE 1 1", "DISTANCE 1 0"] {
        assert!(is_pert_header(line), "{line}");
    }
    assert!(!is_pert_header("MAIN PROBLEM. 1023"));
    assert!(!is_pert_header(PERT_LINE));
}

#[test]
fn header_errors() {
    let main = [
        ("AB", "Main header too short: 'AB'"),
        ("ELP MAIN abc", "Invalid term count in: 'ELP MAIN abc'"),
        ("ELP MAIN -1", "Invalid term count in: 'ELP MAIN -1'"),
    ];
    for (line, message) in main {
        let error = parse_main_header(line).unwrap_err();
        assert_eq!(error.to_string(), format!("Invalid header: {message}"));
    }
    let pert = [
        ("AB", "Pert header too short: 'AB'"),
        ("PERT LON abc 0", "Invalid term count in: 'PERT LON abc 0'"),
        (
            "PERT LON 100 xyz",
            "Invalid time power in: 'PERT LON 100 xyz'",
        ),
        (
            "PERT LON 100 256",
            "Invalid time power in: 'PERT LON 100 256'",
        ),
    ];
    for (line, message) in pert {
        let error = parse_pert_header(line).unwrap_err();
        assert_eq!(error.to_string(), format!("Invalid header: {message}"));
    }
}

#[test]
fn fortran_doubles() {
    let cases = [
        ("-0.1274921554086D+02", -0.1274921554086e2),
        (" 0.6368794709728D+01", 0.6368794709728e1),
        ("1.5d+00", 1.5),
        ("1.234e+05", 1.234e5),
    ];
    for (text, value) in cases {
        assert_eq!(parse_fortran_double(text).unwrap(), value, "{text}");
    }
    for (text, shown) in [
        ("not_a_number", "not_a_number"),
        ("NaN", "NaN"),
        ("1D+999", "1E+999"),
    ] {
        let error = parse_fortran_double(text).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!("Invalid pert term: Invalid double: '{shown}'")
        );
    }
}

#[test]
fn main_terms() {
    let term = parse_main_term(MAIN_LINE).unwrap();
    assert_eq!(term.delaunay, [0, 2, 0, 0]);
    assert_eq!(
        term.coeffs,
        [-411.60287, 168.48, -18433.81, -121.62, 0.40, -0.18, 0.00]
    );
    let trailing_blanks = format!("{MAIN_LINE}   ");
    assert_eq!(
        parse_main_term(&trailing_blanks).unwrap().coeffs,
        term.coeffs
    );
}

fn width_error(line: &str, expected: usize) -> String {
    format!(
        "Line too short or too long ({} columns, expected {expected}): '{line}'",
        line.len()
    )
}

#[test]
fn main_term_errors() {
    let short = &MAIN_LINE[..98];
    let long = format!("{MAIN_LINE} 1");
    let accented = format!("{}é", &MAIN_LINE[..98]);
    let cases = [
        ("short".to_string(), width_error("short", 99)),
        (short.to_string(), width_error(short, 99)),
        (long.clone(), width_error(&long, 99)),
        (accented.clone(), width_error(&accented, 99)),
        (
            with_field(MAIN_LINE, 0, "abc"),
            "Invalid delaunay[0]: 'abc'".into(),
        ),
        (
            with_field(MAIN_LINE, 9, " x "),
            "Invalid delaunay[3]: ' x '".into(),
        ),
        (
            with_field(MAIN_LINE, 14, "          inf"),
            "Invalid coeff[0]: '          inf'".into(),
        ),
        (
            with_field(MAIN_LINE, 51, "     abc.def"),
            "Invalid coeff[3]: '     abc.def'".into(),
        ),
        (
            format!("{:<99}", &MAIN_LINE[..27]),
            format!("Invalid coeff[1]: '{}'", " ".repeat(12)),
        ),
    ];
    for (line, message) in cases {
        let error = parse_main_term(&line).unwrap_err();
        assert_eq!(error.to_string(), format!("Invalid main term: {message}"));
    }
}

#[test]
fn pert_terms() {
    let term = parse_pert_term(PERT_LINE).unwrap();
    assert_eq!(
        (term.sin_coeff, term.cos_coeff),
        (-0.1274921554086e2, 0.6368794709728e1)
    );
    assert_eq!(
        term.multipliers,
        [0, 0, 1, 0, 0, -18, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    );
    let counting = "    1 0.1000000000000D+01 0.2000000000000D+01  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16";
    let expected: Vec<i32> = (1..=16).collect();
    assert_eq!(parse_pert_term(counting).unwrap().multipliers, *expected);
}

#[test]
fn pert_term_errors() {
    let short = &PERT_LINE[..92];
    let cases = [
        ("short".to_string(), width_error("short", 93)),
        (short.to_string(), width_error(short, 93)),
        (
            with_field(PERT_LINE, 60, " x "),
            "Invalid multiplier[5]: ' x '".into(),
        ),
        (
            format!("{:<93}", &PERT_LINE[..57]),
            "Invalid multiplier[4]: '   '".into(),
        ),
        (
            with_field(PERT_LINE, 5, "                 NaN"),
            "Invalid double: 'NaN'".into(),
        ),
        (
            with_field(PERT_LINE, 25, "                    "),
            "Invalid double: ''".into(),
        ),
    ];
    for (line, message) in cases {
        let error = parse_pert_term(&line).unwrap_err();
        assert_eq!(error.to_string(), format!("Invalid pert term: {message}"));
    }
}
