use super::*;
use tempfile::TempDir;

const TEST_FILE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../vsop2013/VSOP2013p3.dat");

fn parse_text(content: &str) -> Result<Vsop2013File, ParseError> {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("test.dat");
    std::fs::write(&path, content).unwrap();
    parse_file(&path)
}

#[test]
fn file_blocks() {
    let content = " VSOP2013  3  1  0  2    EARTH VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
    2   1  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.2000000000000000 +00  0.0000000000000000 +00
 VSOP2013  3  2  0  1    EARTH VAR LAMBDA T^0
    3   0  1  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.3000000000000000 +00  0.0000000000000000 +00
";
    let vsop = parse_text(content).unwrap();
    assert_eq!(vsop.planet, 3);
    assert_eq!(vsop.blocks.len(), 2);
    let first = &vsop.blocks[0].header;
    assert_eq!(
        (first.variable, first.time_power, first.term_count),
        (Variable::A, 0, 2)
    );
    assert_eq!(vsop.blocks[0].terms.len(), 2);
    assert_eq!(vsop.blocks[1].header.variable, Variable::Lambda);
    assert_eq!(vsop.blocks[1].terms.len(), 1);
}

#[test]
fn short_blocks_are_rejected() {
    let at_end = " VSOP2013  3  1  0  5    EARTH VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
";
    // A short block followed by another header, rather than by end of file.
    let mid_file = " VSOP2013  3  1  0  5    EARTH VAR A T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.1000000000000000 +00  0.0000000000000000 +00
 VSOP2013  3  2  0  1    EARTH VAR LAMBDA T^0
    1   0  0  0  0   0  0  0  0  0    0   0   0   0      0   0  0  0  0.2000000000000000 +00  0.0000000000000000 +00
";
    for content in [at_end, mid_file] {
        let error = parse_text(content).unwrap_err();
        assert_eq!(error.to_string(), "Expected 5 terms, found 1");
    }
}

#[test]
fn files_without_a_header_are_rejected() {
    for content in ["", "\n"] {
        let error = parse_text(content).unwrap_err().to_string();
        assert!(
            error.starts_with("Invalid header: No VSOP2013 header in "),
            "{error}"
        );
    }
}

#[test]
fn missing_file_is_an_io_error() {
    let error = parse_file(Path::new("/nonexistent/path/file.dat")).unwrap_err();
    assert!(matches!(error, ParseError::IoError(_)), "{error}");
}

#[test]
#[ignore = "requires local VSOP2013 data files"]
fn shipped_emb_file() {
    let vsop = parse_file(Path::new(TEST_FILE)).unwrap();
    assert_eq!(vsop.planet, 3);
    let first = &vsop.blocks[0];
    assert_eq!(first.header.variable, Variable::A);
    assert_eq!(first.header.time_power, 0);
    assert_eq!(first.header.term_count, 32658);
    assert_eq!(first.terms.len(), 32658);
    assert_eq!(vsop.total_terms(), 294426);
    assert!(vsop
        .blocks_for_variable(Variable::A)
        .all(|b| b.header.variable == Variable::A));
}
