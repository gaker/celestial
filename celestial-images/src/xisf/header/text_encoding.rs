use crate::xisf::errors::{Result, XisfError};

// The text encodings of inline and embedded data blocks (XISF 1.0 section
// 10.3): Base64 and Base16 as RFC 4648 defines them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TextEncoding {
    Base64,
    Hex,
}

impl TextEncoding {
    pub(crate) fn parse(name: &str) -> Result<Self> {
        match name {
            "base64" => Ok(Self::Base64),
            "hex" => Ok(Self::Hex),
            _ => Err(invalid(format!("Unknown data block encoding \"{name}\""))),
        }
    }

    // The location attribute of an inline block in this encoding.
    pub(crate) fn inline_location(self) -> &'static str {
        match self {
            Self::Base64 => "inline:base64",
            Self::Hex => "inline:hex",
        }
    }

    // Section 10.3 makes white space irrelevant in both encodings. Encoders
    // write lowercase Base16 digits; uppercase digits are accepted too.
    pub(crate) fn decode(self, text: &str) -> Result<Vec<u8>> {
        let symbols = || text.chars().filter(|&c| !is_xml_space(c));
        match self {
            Self::Base64 => decode_base64(symbols, text.len()),
            Self::Hex => decode_hex(symbols(), text.len()),
        }
    }
}

// Digits without white space, as in a checksum digest (section 10.5).
pub(super) fn decode_hex_digits(text: &str) -> Result<Vec<u8>> {
    decode_hex(text.chars(), text.len())
}

fn is_xml_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\r' | '\n')
}

// `text_len` bounds the output: two digits make one byte.
fn decode_hex(symbols: impl Iterator<Item = char>, text_len: usize) -> Result<Vec<u8>> {
    let mut output = Vec::with_capacity(text_len / 2);
    let mut high = None;
    for c in symbols {
        let digit = c
            .to_digit(16)
            .ok_or_else(|| invalid_character("hexadecimal", c))?;
        match high.take() {
            // Both digits are below 16, so the byte cannot truncate.
            Some(high) => output.push((high * 16 + digit) as u8),
            None => high = Some(digit),
        }
    }
    if high.is_some() {
        let digits = output.len() * 2 + 1;
        let message = format!("Hexadecimal text has an odd number of digits ({digits})");
        return Err(invalid(message));
    }
    Ok(output)
}

// Bits left over after the last byte are ignored, as Python's
// base64.b64decode does. `text_len` bounds the output: four symbols make at
// most three bytes.
fn decode_base64<I: Iterator<Item = char>>(
    symbols: impl Fn() -> I,
    text_len: usize,
) -> Result<Vec<u8>> {
    let data_symbols = count_data_symbols(symbols())?;
    let mut output = Vec::with_capacity(text_len / 4 * 3);
    let (mut bits, mut pending) = (0u32, 0u32);
    for c in symbols().take(data_symbols) {
        bits = (bits << 6) | base64_value(c)?;
        pending += 6;
        if pending >= 8 {
            pending -= 8;
            // The shift leaves 8 bits, so the byte cannot truncate.
            output.push((bits >> pending) as u8);
            bits &= (1 << pending) - 1;
        }
    }
    Ok(output)
}

// Returns the number of symbols before the padding. Padding is required,
// and only at the end.
fn count_data_symbols(symbols: impl Iterator<Item = char>) -> Result<usize> {
    let (count, padding) = symbols.fold((0usize, 0usize), |(count, run), c| {
        (count + 1, if c == '=' { run + 1 } else { 0 })
    });
    if !count.is_multiple_of(4) {
        let message = format!("Base64 text has {count} characters, which is not a multiple of 4");
        return Err(invalid(message));
    }
    if padding > 2 {
        let message =
            format!("Base64 text ends with {padding} '=' characters; at most 2 are allowed");
        return Err(invalid(message));
    }
    Ok(count - padding)
}

fn base64_value(c: char) -> Result<u32> {
    let value = match c {
        'A'..='Z' => u32::from(c) - u32::from('A'),
        'a'..='z' => u32::from(c) - u32::from('a') + 26,
        '0'..='9' => u32::from(c) - u32::from('0') + 52,
        '+' => 62,
        '/' => 63,
        _ => return Err(invalid_character("Base64", c)),
    };
    Ok(value)
}

fn invalid_character(encoding: &str, c: char) -> XisfError {
    invalid(format!("Invalid {encoding} character {c:?}"))
}

fn invalid(message: String) -> XisfError {
    XisfError::InvalidFormat(message)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base64(text: &str) -> Result<Vec<u8>> {
        TextEncoding::Base64.decode(text)
    }

    fn hex(text: &str) -> Result<Vec<u8>> {
        TextEncoding::Hex.decode(text)
    }

    fn message(result: Result<Vec<u8>>) -> String {
        result.unwrap_err().to_string()
    }

    // Expected values from Python's base64.b64decode(text, validate=True).
    #[test]
    fn base64_matches_python() {
        assert_eq!(base64("TWFu").unwrap(), b"Man");
        assert_eq!(base64("TWE=").unwrap(), b"Ma");
        assert_eq!(base64("TQ==").unwrap(), b"M");
        assert_eq!(base64(" T W\nE=\t").unwrap(), b"Ma");
        assert_eq!(base64("").unwrap(), b"");
        assert_eq!(base64("+/+/").unwrap(), vec![251, 255, 191]);
        assert_eq!(base64("TR==").unwrap(), b"M");
        assert_eq!(
            base64("AAAAAAAA8D8=").unwrap(),
            vec![0, 0, 0, 0, 0, 0, 240, 63]
        );
        // The example of XISF 1.0 section 10.3.
        assert_eq!(
            base64("VGhpcyBpcyBhIHRlc3QgLSBURVNUIC0gMTIzNDU2Nzg5MA==").unwrap(),
            b"This is a test - TEST - 1234567890"
        );
    }

    // Python's base64.b64decode(text, validate=True) raises binascii.Error
    // for each of these.
    #[test]
    fn base64_rejects_what_python_rejects() {
        let cases = [
            (
                "TQ=",
                "Base64 text has 3 characters, which is not a multiple of 4",
            ),
            (
                "TQ",
                "Base64 text has 2 characters, which is not a multiple of 4",
            ),
            (
                "T===",
                "Base64 text ends with 3 '=' characters; at most 2 are allowed",
            ),
            ("TQ==TQ==", "Invalid Base64 character '='"),
            ("TW!u", "Invalid Base64 character '!'"),
            ("TQ=A", "Invalid Base64 character '='"),
            ("=AAA", "Invalid Base64 character '='"),
            ("TW-u", "Invalid Base64 character '-'"),
            ("TWFé", "Invalid Base64 character 'é'"),
        ];
        for (text, expected) in cases {
            assert_eq!(
                message(base64(text)),
                format!("Invalid XISF format: {expected}"),
                "{text}"
            );
        }
    }

    // Expected values from Python's binascii.unhexlify.
    #[test]
    fn hex_matches_python() {
        assert_eq!(hex("00ff10").unwrap(), vec![0, 255, 16]);
        assert_eq!(hex("00FF10").unwrap(), vec![0, 255, 16]);
        assert_eq!(hex(" 00 ff\r\n10 ").unwrap(), vec![0, 255, 16]);
        assert_eq!(hex("").unwrap(), b"");
    }

    #[test]
    fn hex_rejects_what_python_rejects() {
        let cases = [
            ("abc", "Hexadecimal text has an odd number of digits (3)"),
            ("zz", "Invalid hexadecimal character 'z'"),
            ("0x10", "Invalid hexadecimal character 'x'"),
            ("+1", "Invalid hexadecimal character '+'"),
        ];
        for (text, expected) in cases {
            assert_eq!(
                message(hex(text)),
                format!("Invalid XISF format: {expected}"),
                "{text}"
            );
        }
    }

    #[test]
    fn hex_digits_refuse_white_space() {
        assert_eq!(decode_hex_digits("0aff").unwrap(), vec![10, 255]);
        assert_eq!(
            message(decode_hex_digits("0a ff")),
            "Invalid XISF format: Invalid hexadecimal character ' '"
        );
    }

    #[test]
    fn encoding_names_are_parsed() {
        assert_eq!(TextEncoding::parse("base64").unwrap(), TextEncoding::Base64);
        assert_eq!(TextEncoding::parse("hex").unwrap(), TextEncoding::Hex);
        assert_eq!(
            TextEncoding::parse("Base64").unwrap_err().to_string(),
            "Invalid XISF format: Unknown data block encoding \"Base64\""
        );
    }
}
