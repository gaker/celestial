use super::*;

mod files;
mod records;

const TERM_LINE: &str = "    2   0  0  2  0   0  0  0  0  0   -2   0   0   0      0   0  0  0 -0.7736236063963646 -08  0.1120495653357545 -04";

fn block(variable: Variable, time_power: u8, terms: usize) -> Vsop2013Block {
    let term = Vsop2013Term {
        multipliers: [0; 17],
        s_coeff: 1.0,
        c_coeff: 1.0,
    };
    Vsop2013Block {
        header: Vsop2013Header {
            planet: 3,
            variable,
            time_power,
            term_count: terms as u32,
        },
        terms: vec![term; terms],
    }
}

#[test]
fn variables_by_index_name_and_symbol() {
    let expected = [
        ("A", "A (semi-major axis)"),
        ("L", "Lambda (mean longitude)"),
        ("K", "K (e*cos(perihelion))"),
        ("H", "H (e*sin(perihelion))"),
        ("Q", "Q (sin(i/2)*cos(node))"),
        ("P", "P (sin(i/2)*sin(node))"),
    ];
    for (i, (variable, (symbol, name))) in Variable::ALL.iter().zip(expected).enumerate() {
        assert_eq!(Variable::from_index(i as u8 + 1), Some(*variable));
        assert_eq!(variable.to_string(), symbol);
        assert_eq!(variable.name(), name);
    }
    for index in [0, 7, 100] {
        assert_eq!(Variable::from_index(index), None);
    }
}

#[test]
fn parse_error_messages() {
    let io = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
    let cases = [
        (ParseError::IoError(io), "IO error: file not found"),
        (
            ParseError::InvalidHeader("bad header".to_string()),
            "Invalid header: bad header",
        ),
        (
            ParseError::InvalidTerm("bad term".to_string()),
            "Invalid term: bad term",
        ),
        (
            ParseError::InvalidVariable(99),
            "Invalid variable index: 99",
        ),
        (
            ParseError::MissingTerms {
                expected: 100,
                found: 50,
            },
            "Expected 100 terms, found 50",
        ),
    ];
    for (error, message) in cases {
        assert_eq!(error.to_string(), message);
    }
}

#[test]
fn io_errors_convert() {
    let io = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "access denied");
    let error: ParseError = io.into();
    assert_eq!(error.to_string(), "IO error: access denied");
}

#[test]
fn amplitude_is_the_hypotenuse() {
    let term = Vsop2013Term {
        multipliers: [0; 17],
        s_coeff: 3.0,
        c_coeff: 4.0,
    };
    assert_eq!(term.amplitude(), 5.0);
}

#[test]
fn totals_and_blocks_by_variable() {
    let vsop = Vsop2013File {
        planet: 3,
        blocks: vec![
            block(Variable::A, 0, 2),
            block(Variable::Lambda, 0, 1),
            block(Variable::A, 1, 1),
        ],
    };
    assert_eq!(vsop.total_terms(), 4);
    let powers: Vec<u8> = vsop
        .blocks_for_variable(Variable::A)
        .map(|b| b.header.time_power)
        .collect();
    assert_eq!(powers, [0, 1]);
    assert_eq!(vsop.blocks_for_variable(Variable::Lambda).count(), 1);
    assert_eq!(vsop.blocks_for_variable(Variable::K).count(), 0);
}

#[test]
fn planet_names() {
    let names = [
        "Mercury",
        "Venus",
        "Earth-Moon Barycenter",
        "Mars",
        "Jupiter",
        "Saturn",
        "Uranus",
        "Neptune",
        "Pluto",
    ];
    for (i, name) in names.iter().enumerate() {
        assert_eq!(planet_name(i as u8 + 1), *name);
    }
    for planet in [0, 10, 255] {
        assert_eq!(planet_name(planet), "Unknown");
    }
}
