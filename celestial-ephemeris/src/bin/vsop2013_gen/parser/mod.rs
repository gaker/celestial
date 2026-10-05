use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Variable {
    A,
    Lambda,
    K,
    H,
    Q,
    P,
}

impl Variable {
    pub(crate) const ALL: [Variable; 6] = [
        Variable::A,
        Variable::Lambda,
        Variable::K,
        Variable::H,
        Variable::Q,
        Variable::P,
    ];

    pub(crate) fn from_index(idx: u8) -> Option<Self> {
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

    pub(crate) fn name(&self) -> &'static str {
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
pub(crate) struct Vsop2013Header {
    pub(crate) planet: u8,
    pub(crate) variable: Variable,
    pub(crate) time_power: u8,
    pub(crate) term_count: u32,
}

#[derive(Debug, Clone)]
pub(crate) struct Vsop2013Term {
    pub(crate) multipliers: [i32; 17],
    pub(crate) s_coeff: f64,
    pub(crate) c_coeff: f64,
}

impl Vsop2013Term {
    pub(crate) fn amplitude(&self) -> f64 {
        libm::sqrt(self.s_coeff * self.s_coeff + self.c_coeff * self.c_coeff)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct Vsop2013Block {
    pub(crate) header: Vsop2013Header,
    pub(crate) terms: Vec<Vsop2013Term>,
}

#[derive(Debug, Clone)]
pub(crate) struct Vsop2013File {
    pub(crate) planet: u8,
    pub(crate) blocks: Vec<Vsop2013Block>,
}

impl Vsop2013File {
    pub(crate) fn total_terms(&self) -> usize {
        self.blocks.iter().map(|b| b.terms.len()).sum()
    }

    pub(crate) fn blocks_for_variable(
        &self,
        var: Variable,
    ) -> impl Iterator<Item = &Vsop2013Block> {
        self.blocks.iter().filter(move |b| b.header.variable == var)
    }
}

#[derive(Debug)]
pub(crate) enum ParseError {
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

// Parsing the decimal text in one step rounds once; scaling a parsed
// mantissa by a power of ten rounds twice and can miss by several ULP.
fn parse_fortran_float(s: &str, exp: &str) -> Result<f64, ParseError> {
    let exponent: i32 = exp
        .trim()
        .parse()
        .map_err(|_| ParseError::InvalidTerm(format!("Invalid exponent: '{}'", exp)))?;
    let value: f64 = format!("{}e{}", s.trim(), exponent)
        .parse()
        .map_err(|_| ParseError::InvalidTerm(format!("Invalid mantissa: '{}'", s)))?;
    if !value.is_finite() {
        return Err(ParseError::InvalidTerm(format!(
            "Coefficient {} {} is not finite",
            s, exp
        )));
    }
    Ok(value)
}

fn parse_header(line: &str) -> Result<Vsop2013Header, ParseError> {
    let parts = header_parts(line)?;
    let var_idx = header_field(parts[2], "variable")?;
    Ok(Vsop2013Header {
        planet: header_field(parts[1], "planet")?,
        variable: Variable::from_index(var_idx).ok_or(ParseError::InvalidVariable(var_idx))?,
        time_power: header_field(parts[3], "time power")?,
        term_count: header_field(parts[4], "term count")?,
    })
}

fn header_parts(line: &str) -> Result<Vec<&str>, ParseError> {
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
    Ok(parts)
}

fn header_field<T: FromStr>(field: &str, name: &str) -> Result<T, ParseError> {
    field
        .parse()
        .map_err(|_| ParseError::InvalidHeader(format!("Invalid {}: '{}'", name, field)))
}

// A sign straight after a digit starts a new token: the files can run an
// exponent into its mantissa.
fn tokenize_numbers(s: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let (mut start, mut after_digit) = (None, false);
    for (i, ch) in s.bytes().enumerate() {
        let is_sign = ch == b'-' || ch == b'+';
        let numeric = ch.is_ascii_digit() || is_sign || ch == b'.';
        if let Some(from) = start.filter(|_| !numeric || (is_sign && after_digit)) {
            tokens.push(&s[from..i]);
            start = None;
        }
        if numeric && start.is_none() {
            start = Some(i);
        }
        after_digit = ch.is_ascii_digit();
    }
    tokens.extend(start.map(|from| &s[from..]));
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
    Ok(Vsop2013Term {
        multipliers: parse_multipliers(&tokens[1..18])?,
        s_coeff: parse_fortran_float(tokens[18], tokens[19])?,
        c_coeff: parse_fortran_float(tokens[20], tokens[21])?,
    })
}

fn parse_multipliers(tokens: &[&str]) -> Result<[i32; 17], ParseError> {
    let mut multipliers = [0i32; 17];
    for (m, token) in multipliers.iter_mut().zip(tokens) {
        *m = token
            .parse()
            .map_err(|_| ParseError::InvalidTerm(format!("Invalid multiplier: '{}'", token)))?;
    }
    Ok(multipliers)
}

pub(crate) fn parse_file(path: &Path) -> Result<Vsop2013File, ParseError> {
    let blocks = read_blocks(BufReader::new(File::open(path)?))?;
    let first = blocks.first().ok_or_else(|| {
        ParseError::InvalidHeader(format!("No VSOP2013 header in {}", path.display()))
    })?;
    Ok(Vsop2013File {
        planet: first.header.planet,
        blocks,
    })
}

fn read_blocks(reader: impl BufRead) -> Result<Vec<Vsop2013Block>, ParseError> {
    let mut blocks: Vec<Vsop2013Block> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with(" VSOP2013") {
            check_complete(blocks.last())?;
            let header = parse_header(&line)?;
            blocks.push(Vsop2013Block {
                header,
                terms: Vec::new(),
            });
        } else if let Some(block) = blocks.last_mut() {
            block.terms.push(parse_term(&line)?);
        }
    }
    check_complete(blocks.last())?;
    Ok(blocks)
}

fn check_complete(block: Option<&Vsop2013Block>) -> Result<(), ParseError> {
    match block {
        Some(b) if b.terms.len() as u32 != b.header.term_count => Err(ParseError::MissingTerms {
            expected: b.header.term_count,
            found: b.terms.len() as u32,
        }),
        _ => Ok(()),
    }
}

pub(crate) fn planet_name(planet: u8) -> &'static str {
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
mod tests;
