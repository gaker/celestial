use super::*;

#[test]
fn header_fields() {
    let line = " VSOP2013  3  1  0  32658    EARTH-MOON VARIABLE A   *T*00";
    let header = parse_header(line).unwrap();
    assert_eq!(header.planet, 3);
    assert_eq!(header.variable, Variable::A);
    assert_eq!(header.time_power, 0);
    assert_eq!(header.term_count, 32658);
    for (i, variable) in Variable::ALL.iter().enumerate() {
        let line = format!(" VSOP2013  5  {}  2  100    JUPITER VAR", i + 1);
        assert_eq!(parse_header(&line).unwrap().variable, *variable);
    }
}

#[test]
fn header_errors() {
    let cases = [
        (
            "SOME OTHER FORMAT",
            "Invalid header: Line doesn't start with ' VSOP2013': 'SOME OTHER FORMAT'",
        ),
        (
            " VSOP2013 1 2",
            "Invalid header: Not enough parts in header: ' VSOP2013 1 2'",
        ),
        (
            " VSOP2013  X  1  0  100    INVALID",
            "Invalid header: Invalid planet: 'X'",
        ),
        (
            " VSOP2013  3  Y  0  100    INVALID",
            "Invalid header: Invalid variable: 'Y'",
        ),
        (
            " VSOP2013  3  1  Z  100    INVALID",
            "Invalid header: Invalid time power: 'Z'",
        ),
        (
            " VSOP2013  3  1  0  XXX    INVALID",
            "Invalid header: Invalid term count: 'XXX'",
        ),
        (
            " VSOP2013  3  9  0  100    INVALID",
            "Invalid variable index: 9",
        ),
    ];
    for (line, message) in cases {
        assert_eq!(parse_header(line).unwrap_err().to_string(), message);
    }
}

#[test]
fn tokens() {
    let cases: [(&str, &[&str]); 7] = [
        ("1 2 3", &["1", "2", "3"]),
        ("-1 +2 -3", &["-1", "+2", "-3"]),
        ("1.5 -2.3 +0.7", &["1.5", "-2.3", "+0.7"]),
        ("-0.7736236063963646 -08", &["-0.7736236063963646", "-08"]),
        // VSOP files can put the exponent's sign right after the mantissa.
        ("1.234-05", &["1.234", "-05"]),
        ("9.87+03", &["9.87", "+03"]),
        ("1.0-01 2.0+02", &["1.0", "-01", "2.0", "+02"]),
    ];
    for (text, expected) in cases {
        assert_eq!(tokenize_numbers(text), expected, "{text}");
    }
    let tokens = tokenize_numbers(TERM_LINE);
    assert_eq!(tokens.len(), 22);
    assert_eq!((tokens[0], tokens[1], tokens[10]), ("2", "0", "-2"));
}

#[test]
fn term_fields() {
    let term = parse_term(TERM_LINE).unwrap();
    let mut multipliers = [0; 17];
    multipliers[2] = 2;
    multipliers[9] = -2;
    assert_eq!(term.multipliers, multipliers);
    assert_eq!(term.s_coeff, -7.736236063963646e-9);
    assert_eq!(term.c_coeff, 1.120495653357545e-5);
}

#[test]
fn term_errors() {
    let wide = "    2   9999999999999  0  2  0   0  0  0  0  0   -2   0   0   0      0   0  0  0 -0.77 -08  0.11 -04";
    let cases = [
        (
            "1 2 3 4 5",
            "Invalid term: Not enough tokens (5): '1 2 3 4 5'",
        ),
        (wide, "Invalid term: Invalid multiplier: '9999999999999'"),
    ];
    for (line, message) in cases {
        assert_eq!(parse_term(line).unwrap_err().to_string(), message);
    }
}

#[test]
fn fortran_floats() {
    let cases = [
        ("-0.7736236063963646", "-08", -0.7736236063963646e-8),
        ("0.1000001017641000", "+01", 0.1000001017641e1),
        ("0.1120495653357545", "-04", 0.1120495653357545e-4),
        ("-0.3107290842332163", "-04", -0.3107290842332163e-4),
        ("0.5970513434294535", "-05", 0.5970513434294535e-5),
    ];
    for (mantissa, exponent, expected) in cases {
        assert_eq!(parse_fortran_float(mantissa, exponent).unwrap(), expected);
    }
}

#[test]
fn fortran_float_errors() {
    let cases = [
        ("not_a_number", "-08", "Invalid mantissa: 'not_a_number'"),
        ("0.123", "abc", "Invalid exponent: 'abc'"),
        ("0.1", "+999", "Coefficient 0.1 +999 is not finite"),
        ("inf", "+00", "Invalid mantissa: 'inf'"),
        ("NaN", "+00", "Invalid mantissa: 'NaN'"),
    ];
    for (mantissa, exponent, message) in cases {
        let error = parse_fortran_float(mantissa, exponent).unwrap_err();
        assert_eq!(error.to_string(), format!("Invalid term: {message}"));
    }
}
