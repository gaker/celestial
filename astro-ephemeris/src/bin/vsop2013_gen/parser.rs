use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Variable {
    A,
    Lambda,
    K,
    H,
    Q,
    P,
}

impl Variable {
    pub fn from_index(idx: u8) -> Option<Self> {
        match idx {
            1 => Some(Variable::A),
            2 => Some(Variable::Lambda),
            3 => Some(Variable::K),
            4 => Some(Variable::H),
            5 => Some(Variable::Q),
            6 => Some(Variable::P),
            _ => None,
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            Variable::A => "A (semi-major axis)",
            Variable::Lambda => "Lambda (mean longitude)",
            Variable::K => "K (e*cos(perihelion))",
            Variable::H => "H (e*sin(perihelion))",
            Variable::Q => "Q (sin(i/2)*cos(node))",
            Variable::P => "P (sin(i/2)*sin(node))",
        }
    }
}

impl fmt::Display for Variable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Variable::A => "A",
            Variable::Lambda => "L",
            Variable::K => "K",
            Variable::H => "H",
            Variable::Q => "Q",
            Variable::P => "P",
        };
        write!(f, "{}", name)
    }
}

#[derive(Debug, Clone)]
pub struct Vsop2013Header {
    pub planet: u8,
    pub variable: Variable,
    pub time_power: u8,
    pub term_count: u32,
}

#[derive(Debug, Clone)]
pub struct Vsop2013Term {
    pub multipliers: [i32; 17],
    pub s_coeff: f64,
    pub c_coeff: f64,
}

impl Vsop2013Term {
    pub fn amplitude(&self) -> f64 {
        (self.s_coeff.powi(2) + self.c_coeff.powi(2)).sqrt()
    }
}

#[derive(Debug, Clone)]
pub struct Vsop2013Block {
    pub header: Vsop2013Header,
    pub terms: Vec<Vsop2013Term>,
}

#[derive(Debug, Clone)]
pub struct Vsop2013File {
    pub planet: u8,
    pub blocks: Vec<Vsop2013Block>,
}

impl Vsop2013File {
    pub fn total_terms(&self) -> usize {
        self.blocks.iter().map(|b| b.terms.len()).sum()
    }

    pub fn blocks_for_variable(&self, var: Variable) -> impl Iterator<Item = &Vsop2013Block> {
        self.blocks.iter().filter(move |b| b.header.variable == var)
    }
}

#[derive(Debug)]
pub enum ParseError {
    IoError(std::io::Error),
    InvalidHeader(String),
    InvalidTerm(String),
    InvalidVariable(u8),
    MissingTerms { expected: u32, found: u32 },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::IoError(e) => write!(f, "IO error: {}", e),
            ParseError::InvalidHeader(s) => write!(f, "Invalid header: {}", s),
            ParseError::InvalidTerm(s) => write!(f, "Invalid term: {}", s),
            ParseError::InvalidVariable(v) => write!(f, "Invalid variable index: {}", v),
            ParseError::MissingTerms { expected, found } => {
                write!(f, "Expected {} terms, found {}", expected, found)
            }
        }
    }
}

impl std::error::Error for ParseError {}

impl From<std::io::Error> for ParseError {
    fn from(e: std::io::Error) -> Self {
        ParseError::IoError(e)
    }
}

fn parse_fortran_float(s: &str, exp: &str) -> Result<f64, ParseError> {
    let mantissa: f64 = s
        .trim()
        .parse()
        .map_err(|_| ParseError::InvalidTerm(format!("Invalid mantissa: '{}'", s)))?;
    let exponent: i32 = exp
        .trim()
        .parse()
        .map_err(|_| ParseError::InvalidTerm(format!("Invalid exponent: '{}'", exp)))?;
    Ok(mantissa * 10f64.powi(exponent))
}

fn parse_header(line: &str) -> Result<Vsop2013Header, ParseError> {
    if !line.starts_with(" VSOP2013") {
        return Err(ParseError::InvalidHeader(format!(
            "Line doesn't start with ' VSOP2013': '{}'",
            line
        )));
    }

    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 5 {
        return Err(ParseError::InvalidHeader(format!(
            "Not enough parts in header: '{}'",
            line
        )));
    }

    let planet: u8 = parts[1]
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid planet: '{}'", parts[1])))?;
    let var_idx: u8 = parts[2]
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid variable: '{}'", parts[2])))?;
    let time_power: u8 = parts[3]
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid time power: '{}'", parts[3])))?;
    let term_count: u32 = parts[4]
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid term count: '{}'", parts[4])))?;

    let variable = Variable::from_index(var_idx).ok_or(ParseError::InvalidVariable(var_idx))?;

    Ok(Vsop2013Header {
        planet,
        variable,
        time_power,
        term_count,
    })
}

fn tokenize_numbers(s: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let mut start = None;
    let bytes = s.as_bytes();

    for (i, &ch) in bytes.iter().enumerate() {
        let is_digit = ch.is_ascii_digit();
        let is_sign = ch == b'-' || ch == b'+';
        let is_dot = ch == b'.';

        match start {
            None => {
                if is_digit || is_sign || is_dot {
                    start = Some(i);
                }
            }
            Some(s_idx) => {
                if is_sign && i > 0 && bytes[i - 1].is_ascii_digit() {
                    tokens.push(&s[s_idx..i]);
                    start = Some(i);
                } else if !is_digit && !is_sign && !is_dot {
                    tokens.push(&s[s_idx..i]);
                    start = None;
                }
            }
        }
    }
    if let Some(s_idx) = start {
        tokens.push(&s[s_idx..]);
    }
    tokens
}

fn parse_term(line: &str) -> Result<Vsop2013Term, ParseError> {
    let tokens = tokenize_numbers(line);

    if tokens.len() < 22 {
        return Err(ParseError::InvalidTerm(format!(
            "Not enough tokens ({}): '{}'",
            tokens.len(),
            line
        )));
    }

    let mut multipliers = [0i32; 17];
    for i in 0..17 {
        multipliers[i] = tokens[i + 1].parse().map_err(|_| {
            ParseError::InvalidTerm(format!("Invalid multiplier: '{}'", tokens[i + 1]))
        })?;
    }

    let s_coeff = parse_fortran_float(tokens[18], tokens[19])?;
    let c_coeff = parse_fortran_float(tokens[20], tokens[21])?;

    Ok(Vsop2013Term {
        multipliers,
        s_coeff,
        c_coeff,
    })
}

pub fn parse_file(path: &Path) -> Result<Vsop2013File, ParseError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let lines = reader.lines();

    let mut blocks = Vec::new();
    let mut current_block: Option<Vsop2013Block> = None;
    let mut planet: Option<u8> = None;

    for line_result in lines {
        let line = line_result?;

        if line.starts_with(" VSOP2013") {
            if let Some(block) = current_block.take() {
                if block.terms.len() as u32 != block.header.term_count {
                    return Err(ParseError::MissingTerms {
                        expected: block.header.term_count,
                        found: block.terms.len() as u32,
                    });
                }
                blocks.push(block);
            }

            let header = parse_header(&line)?;
            if planet.is_none() {
                planet = Some(header.planet);
            }

            current_block = Some(Vsop2013Block {
                header,
                terms: Vec::new(),
            });
        } else if let Some(ref mut block) = current_block {
            let term = parse_term(&line)?;
            block.terms.push(term);
        }
    }

    if let Some(block) = current_block {
        if block.terms.len() as u32 != block.header.term_count {
            return Err(ParseError::MissingTerms {
                expected: block.header.term_count,
                found: block.terms.len() as u32,
            });
        }
        blocks.push(block);
    }

    Ok(Vsop2013File {
        planet: planet.unwrap_or(0),
        blocks,
    })
}

pub fn planet_name(planet: u8) -> &'static str {
    match planet {
        1 => "Mercury",
        2 => "Venus",
        3 => "Earth-Moon Barycenter",
        4 => "Mars",
        5 => "Jupiter",
        6 => "Saturn",
        7 => "Uranus",
        8 => "Neptune",
        9 => "Pluto",
        _ => "Unknown",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    const TEST_FILE: &str = "references/ephemeris/vsop2013/VSOP2013p3.dat";

    #[test]
    fn test_variable_from_index_all_valid() {
        assert_eq!(Variable::from_index(1), Some(Variable::A));
        assert_eq!(Variable::from_index(2), Some(Variable::Lambda));
        assert_eq!(Variable::from_index(3), Some(Variable::K));
        assert_eq!(Variable::from_index(4), Some(Variable::H));
        assert_eq!(Variable::from_index(5), Some(Variable::Q));
        assert_eq!(Variable::from_index(6), Some(Variable::P));
    }

    #[test]
    fn test_variable_from_index_invalid() {
        assert_eq!(Variable::from_index(0), None);
        assert_eq!(Variable::from_index(7), None);
        assert_eq!(Variable::from_index(100), None);
    }

    #[test]
    fn test_variable_name_all_variants() {
        assert_eq!(Variable::A.name(), "A (semi-major axis)");
        assert_eq!(Variable::Lambda.name(), "Lambda (mean longitude)");
        assert_eq!(Variable::K.name(), "K (e*cos(perihelion))");
        assert_eq!(Variable::H.name(), "H (e*sin(perihelion))");
        assert_eq!(Variable::Q.name(), "Q (sin(i/2)*cos(node))");
        assert_eq!(Variable::P.name(), "P (sin(i/2)*sin(node))");
    }

    #[test]
    fn test_variable_display_all() {
        assert_eq!(format!("{}", Variable::A), "A");
        assert_eq!(format!("{}", Variable::Lambda), "L");
        assert_eq!(format!("{}", Variable::K), "K");
        assert_eq!(format!("{}", Variable::H), "H");
        assert_eq!(format!("{}", Variable::Q), "Q");
        assert_eq!(format!("{}", Variable::P), "P");
    }

    #[test]
    fn test_parse_error_display_io_error() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
        let err = ParseError::IoError(io_err);
        let display = format!("{}", err);
        assert!(display.contains("IO error"));
        assert!(display.contains("file not found"));
    }

    #[test]
    fn test_parse_error_display_invalid_header() {
        let err = ParseError::InvalidHeader("bad header".to_string());
        let display = format!("{}", err);
        assert!(display.contains("Invalid header"));
        assert!(display.contains("bad header"));
    }

    #[test]
    fn test_parse_error_display_invalid_term() {
        let err = ParseError::InvalidTerm("bad term".to_string());
        let display = format!("{}", err);
        assert!(display.contains("Invalid term"));
        assert!(display.contains("bad term"));
    }

    #[test]
    fn test_parse_error_display_invalid_variable() {
        let err = ParseError::InvalidVariable(99);
        let display = format!("{}", err);
        assert!(display.contains("Invalid variable index"));
        assert!(display.contains("99"));
    }

    #[test]
    fn test_parse_error_display_missing_terms() {
        let err = ParseError::MissingTerms {
            expected: 100,
            found: 50,
        };
        let display = format!("{}", err);
        assert!(display.contains("Expected 100 terms"));
        assert!(display.contains("found 50"));
    }

    #[test]
    fn test_parse_error_from_io_error() {
        let io_err = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "access denied");
        let err: ParseError = io_err.into();
        match err {
            ParseError::IoError(e) => assert!(e.to_string().contains("access denied")),
            _ => panic!("Expected IoError variant"),
        }
    }

    #[test]
    fn test_parse_header() {
        let line = " VSOP2013  3  1  0  32658    EARTH-MOON VARIABLE A   *T*00";
        let header = parse_header(line).unwrap();
        assert_eq!(header.planet, 3);
        assert_eq!(header.variable, Variable::A);
        assert_eq!(header.time_power, 0);
        assert_eq!(header.term_count, 32658);
    }

    #[test]
    fn test_parse_header_all_variables() {
        for (idx, expected_var) in [
            (1, Variable::A),
            (2, Variable::Lambda),
            (3, Variable::K),
            (4, Variable::H),
            (5, Variable::Q),
            (6, Variable::P),
        ] {
            let line = format!(" VSOP2013  5  {}  2  100    JUPITER VAR", idx);
            let header = parse_header(&line).unwrap();
            assert_eq!(header.variable, expected_var);
        }
    }

    #[test]
    fn test_parse_header_error_not_vsop() {
        let line = "SOME OTHER FORMAT";
        let result = parse_header(line);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidHeader(s) => assert!(s.contains("doesn't start with")),
            _ => panic!("Expected InvalidHeader"),
        }
    }

    #[test]
    fn test_parse_header_error_not_enough_parts() {
        let line = " VSOP2013 1 2";
        let result = parse_header(line);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidHeader(s) => assert!(s.contains("Not enough parts")),
            _ => panic!("Expected InvalidHeader"),
        }
    }

    #[test]
    fn test_parse_header_error_invalid_planet() {
        let line = " VSOP2013  X  1  0  100    INVALID";
        let result = parse_header(line);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidHeader(s) => assert!(s.contains("Invalid planet")),
            _ => panic!("Expected InvalidHeader"),
        }
    }

    #[test]
    fn test_parse_header_error_invalid_variable_number() {
        let line = " VSOP2013  3  Y  0  100    INVALID";
        let result = parse_header(line);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidHeader(s) => assert!(s.contains("Invalid variable")),
            _ => panic!("Expected InvalidHeader"),
        }
    }

    #[test]
    fn test_parse_header_error_invalid_time_power() {
        let line = " VSOP2013  3  1  Z  100    INVALID";
        let result = parse_header(line);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidHeader(s) => assert!(s.contains("Invalid time power")),
            _ => panic!("Expected InvalidHeader"),
        }
    }

    #[test]
    fn test_parse_header_error_invalid_term_count() {
        let line = " VSOP2013  3  1  0  XXX    INVALID";
        let result = parse_header(line);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidHeader(s) => assert!(s.contains("Invalid term count")),
            _ => panic!("Expected InvalidHeader"),
        }
    }

    #[test]
    fn test_parse_header_error_invalid_variable_index() {
        let line = " VSOP2013  3  9  0  100    INVALID";
        let result = parse_header(line);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidVariable(idx) => assert_eq!(idx, 9),
            _ => panic!("Expected InvalidVariable"),
        }
    }

    #[test]
    fn test_tokenize_numbers_basic() {
        let tokens = tokenize_numbers("1 2 3");
        assert_eq!(tokens, vec!["1", "2", "3"]);
    }

    #[test]
    fn test_tokenize_numbers_with_signs() {
        let tokens = tokenize_numbers("-1 +2 -3");
        assert_eq!(tokens, vec!["-1", "+2", "-3"]);
    }

    #[test]
    fn test_tokenize_numbers_with_decimals() {
        let tokens = tokenize_numbers("1.5 -2.3 +0.7");
        assert_eq!(tokens, vec!["1.5", "-2.3", "+0.7"]);
    }

    #[test]
    fn test_tokenize_numbers_fortran_style() {
        let tokens = tokenize_numbers("-0.7736236063963646 -08");
        assert_eq!(tokens, vec!["-0.7736236063963646", "-08"]);
    }

    #[test]
    fn test_tokenize_numbers_digit_followed_by_sign() {
        // This tests lines 197-199: when a sign follows directly after a digit (no space)
        // This pattern appears in VSOP files where exponent immediately follows mantissa
        let tokens = tokenize_numbers("1.234-05");
        assert_eq!(tokens, vec!["1.234", "-05"]);

        let tokens = tokenize_numbers("9.87+03");
        assert_eq!(tokens, vec!["9.87", "+03"]);

        // Multiple occurrences
        let tokens = tokenize_numbers("1.0-01 2.0+02");
        assert_eq!(tokens, vec!["1.0", "-01", "2.0", "+02"]);
    }

    #[test]
    fn test_tokenize_numbers_vsop_term_line() {
        let line = "    2   0  0  2  0   0  0  0  0  0   -2   0   0   0      0   0  0  0 -0.7736236063963646 -08  0.1120495653357545 -04";
        let tokens = tokenize_numbers(line);
        assert!(tokens.len() >= 22);
        assert_eq!(tokens[0], "2");
        assert_eq!(tokens[1], "0");
        assert_eq!(tokens[10], "-2");
    }

    #[test]
    fn test_parse_term() {
        let line = "    2   0  0  2  0   0  0  0  0  0   -2   0   0   0      0   0  0  0 -0.7736236063963646 -08  0.1120495653357545 -04";
        let term = parse_term(line).unwrap();

        assert_eq!(term.multipliers[0], 0);
        assert_eq!(term.multipliers[1], 0);
        assert_eq!(term.multipliers[2], 2);
        assert_eq!(term.multipliers[9], -2);

        let expected_s = -7.736236063963646e-9;
        let expected_c = 1.120495653357545e-5;
        assert!((term.s_coeff - expected_s).abs() < 1e-20);
        assert!((term.c_coeff - expected_c).abs() < 1e-16);
    }

    #[test]
    fn test_parse_term_not_enough_tokens() {
        let line = "1 2 3 4 5";
        let result = parse_term(line);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidTerm(s) => assert!(s.contains("Not enough tokens")),
            _ => panic!("Expected InvalidTerm"),
        }
    }

    #[test]
    fn test_parse_term_invalid_multiplier() {
        // Use a value that's too large for i32 to cause a parse error
        let line = "    2   9999999999999  0  2  0   0  0  0  0  0   -2   0   0   0      0   0  0  0 -0.77 -08  0.11 -04";
        let result = parse_term(line);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidTerm(s) => assert!(s.contains("Invalid multiplier")),
            _ => panic!("Expected InvalidTerm"),
        }
    }

    #[test]
    fn test_parse_fortran_float() {
        assert!(
            (parse_fortran_float("-0.7736236063963646", "-08").unwrap() - (-7.736236063963646e-9))
                .abs()
                < 1e-20
        );
        assert!(
            (parse_fortran_float("0.1000001017641000", "+01").unwrap() - 1.000001017641).abs()
                < 1e-15
        );
        assert!(
            (parse_fortran_float("0.1120495653357545", "-04").unwrap() - 1.120495653357545e-5)
                .abs()
                < 1e-16
        );
    }

    #[test]
    fn test_parse_fortran_float_invalid_mantissa() {
        let result = parse_fortran_float("not_a_number", "-08");
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidTerm(s) => assert!(s.contains("Invalid mantissa")),
            _ => panic!("Expected InvalidTerm"),
        }
    }

    #[test]
    fn test_parse_fortran_float_invalid_exponent() {
        let result = parse_fortran_float("0.123", "abc");
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::InvalidTerm(s) => assert!(s.contains("Invalid exponent")),
            _ => panic!("Expected InvalidTerm"),
        }
    }

    #[test]
    fn test_amplitude() {
        let term = Vsop2013Term {
            multipliers: [0; 17],
            s_coeff: 3.0,
            c_coeff: 4.0,
        };
        assert!((term.amplitude() - 5.0).abs() < 1e-15);
    }

    #[test]
    fn test_vsop2013_file_total_terms() {
        let vsop = Vsop2013File {
            planet: 3,
            blocks: vec![
                Vsop2013Block {
                    header: Vsop2013Header {
                        planet: 3,
                        variable: Variable::A,
                        time_power: 0,
                        term_count: 2,
                    },
                    terms: vec![
                        Vsop2013Term {
                            multipliers: [0; 17],
                            s_coeff: 1.0,
                            c_coeff: 1.0,
                        },
                        Vsop2013Term {
                            multipliers: [0; 17],
                            s_coeff: 2.0,
                            c_coeff: 2.0,
                        },
                    ],
                },
                Vsop2013Block {
                    header: Vsop2013Header {
                        planet: 3,
                        variable: Variable::Lambda,
                        time_power: 0,
                        term_count: 1,
                    },
                    terms: vec![Vsop2013Term {
                        multipliers: [0; 17],
                        s_coeff: 3.0,
                        c_coeff: 3.0,
                    }],
                },
            ],
        };
        assert_eq!(vsop.total_terms(), 3);
    }

    #[test]
    fn test_vsop2013_file_blocks_for_variable() {
        let vsop = Vsop2013File {
            planet: 3,
            blocks: vec![
                Vsop2013Block {
                    header: Vsop2013Header {
                        planet: 3,
                        variable: Variable::A,
                        time_power: 0,
                        term_count: 1,
                    },
                    terms: vec![Vsop2013Term {
                        multipliers: [0; 17],
                        s_coeff: 1.0,
                        c_coeff: 1.0,
                    }],
                },
                Vsop2013Block {
                    header: Vsop2013Header {
                        planet: 3,
                        variable: Variable::Lambda,
                        time_power: 0,
                        term_count: 1,
                    },
                    terms: vec![Vsop2013Term {
                        multipliers: [0; 17],
                        s_coeff: 2.0,
                        c_coeff: 2.0,
                    }],
                },
                Vsop2013Block {
                    header: Vsop2013Header {
                        planet: 3,
                        variable: Variable::A,
                        time_power: 1,
                        term_count: 1,
                    },
                    terms: vec![Vsop2013Term {
                        multipliers: [0; 17],
                        s_coeff: 3.0,
                        c_coeff: 3.0,
                    }],
                },
            ],
        };

        let a_blocks: Vec<_> = vsop.blocks_for_variable(Variable::A).collect();
        assert_eq!(a_blocks.len(), 2);
        assert_eq!(a_blocks[0].header.time_power, 0);
        assert_eq!(a_blocks[1].header.time_power, 1);

        let lambda_blocks: Vec<_> = vsop.blocks_for_variable(Variable::Lambda).collect();
        assert_eq!(lambda_blocks.len(), 1);
    }

    #[test]
    fn test_planet_name_all_planets() {
        assert_eq!(planet_name(1), "Mercury");
        assert_eq!(planet_name(2), "Venus");
        assert_eq!(planet_name(3), "Earth-Moon Barycenter");
        assert_eq!(planet_name(4), "Mars");
        assert_eq!(planet_name(5), "Jupiter");
        assert_eq!(planet_name(6), "Saturn");
        assert_eq!(planet_name(7), "Uranus");
        assert_eq!(planet_name(8), "Neptune");
        assert_eq!(planet_name(9), "Pluto");
    }

    #[test]
    fn test_planet_name_unknown() {
        assert_eq!(planet_name(0), "Unknown");
        assert_eq!(planet_name(10), "Unknown");
        assert_eq!(planet_name(255), "Unknown");
    }

    #[test]
    fn test_parse_file_with_mock_data() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("test_vsop.dat");

        let content = " VSOP2013  3  1  0  2    EARTH VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
    2   1  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.2000000000000000 +00  0.0000000000000000 +00
 VSOP2013  3  2  0  1    EARTH VAR LAMBDA T^0
    3   0  1  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.3000000000000000 +00  0.0000000000000000 +00
";
        std::fs::write(&file_path, content).unwrap();

        let result = parse_file(&file_path);
        assert!(result.is_ok(), "Parse failed: {:?}", result.err());

        let vsop = result.unwrap();
        assert_eq!(vsop.planet, 3);
        assert_eq!(vsop.blocks.len(), 2);

        let block1 = &vsop.blocks[0];
        assert_eq!(block1.header.variable, Variable::A);
        assert_eq!(block1.header.time_power, 0);
        assert_eq!(block1.header.term_count, 2);
        assert_eq!(block1.terms.len(), 2);

        let block2 = &vsop.blocks[1];
        assert_eq!(block2.header.variable, Variable::Lambda);
        assert_eq!(block2.terms.len(), 1);
    }

    #[test]
    fn test_parse_file_missing_terms_error() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("bad_vsop.dat");

        let content = " VSOP2013  3  1  0  5    EARTH VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
";
        std::fs::write(&file_path, content).unwrap();

        let result = parse_file(&file_path);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::MissingTerms { expected, found } => {
                assert_eq!(expected, 5);
                assert_eq!(found, 1);
            }
            e => panic!("Expected MissingTerms, got {:?}", e),
        }
    }

    #[test]
    fn test_parse_file_missing_terms_mid_file() {
        // This test exercises lines 256-258: MissingTerms error when it occurs
        // before another header (mid-file), not at end of file
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("bad_vsop_mid.dat");

        // First block claims 5 terms but only has 1, followed by another header
        let content = " VSOP2013  3  1  0  5    EARTH VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
 VSOP2013  3  2  0  1    EARTH VAR LAMBDA T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.2000000000000000 +00  0.0000000000000000 +00
";
        std::fs::write(&file_path, content).unwrap();

        let result = parse_file(&file_path);
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::MissingTerms { expected, found } => {
                assert_eq!(expected, 5);
                assert_eq!(found, 1);
            }
            e => panic!("Expected MissingTerms, got {:?}", e),
        }
    }

    #[test]
    fn test_parse_file_empty() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("empty.dat");
        std::fs::write(&file_path, "").unwrap();

        let result = parse_file(&file_path);
        assert!(result.is_ok());
        let vsop = result.unwrap();
        assert_eq!(vsop.planet, 0);
        assert!(vsop.blocks.is_empty());
    }

    #[test]
    fn test_parse_file_io_error() {
        let result = parse_file(Path::new("/nonexistent/path/file.dat"));
        assert!(result.is_err());
        match result.unwrap_err() {
            ParseError::IoError(_) => {}
            e => panic!("Expected IoError, got {:?}", e),
        }
    }

    #[test]
    #[ignore = "requires local VSOP2013 data files"]
    fn test_parse_file() {
        let path = Path::new(TEST_FILE);
        let result = parse_file(path);
        assert!(result.is_ok(), "Parse failed: {:?}", result.err());

        let vsop = result.unwrap();
        assert_eq!(vsop.planet, 3);
        assert!(!vsop.blocks.is_empty());

        let first_block = &vsop.blocks[0];
        assert_eq!(first_block.header.variable, Variable::A);
        assert_eq!(first_block.header.time_power, 0);
        assert_eq!(first_block.header.term_count, 32658);
        assert_eq!(first_block.terms.len(), 32658);

        let total = vsop.total_terms();
        assert!(total > 100_000, "Expected > 100k terms, got {}", total);
    }

    #[test]
    #[ignore = "requires local VSOP2013 data files"]
    fn test_blocks_for_variable() {
        let path = Path::new(TEST_FILE);
        let vsop = parse_file(path).unwrap();
        let a_blocks: Vec<_> = vsop.blocks_for_variable(Variable::A).collect();
        assert!(!a_blocks.is_empty());
        for block in &a_blocks {
            assert_eq!(block.header.variable, Variable::A);
        }
    }
}
