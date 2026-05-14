use super::BinaryTableHdu;
use crate::fits::{FitsError, Result};

impl BinaryTableHdu {
    pub(super) fn parse_binary_format(&self, format: &str) -> Result<(String, usize)> {
        if format.is_empty() {
            return Err(FitsError::InvalidFormat("Empty column format".to_string()));
        }

        if format.contains('P') || format.contains('Q') {
            return self.parse_variable_length_format(format);
        }

        self.parse_standard_format(format)
    }

    fn parse_variable_length_format(&self, format: &str) -> Result<(String, usize)> {
        let chars: Vec<char> = format.chars().collect();
        let mut repeat_str = String::new();
        let mut i = 0;

        while i < chars.len() && chars[i].is_ascii_digit() {
            repeat_str.push(chars[i]);
            i += 1;
        }

        if i >= chars.len() || (chars[i] != 'P' && chars[i] != 'Q') {
            return Err(FitsError::InvalidFormat(format!(
                "Variable-length format '{}' must contain P or Q",
                format
            )));
        }

        let descriptor = chars[i];
        i += 1;

        if i >= chars.len() {
            return Err(FitsError::InvalidFormat(format!(
                "Missing data type in format '{}'",
                format
            )));
        }

        let data_type = chars[i];
        let repeat = if repeat_str.is_empty() {
            1
        } else {
            repeat_str.parse().unwrap_or(1)
        };

        Ok((format!("{}{}", descriptor, data_type), repeat))
    }

    fn parse_standard_format(&self, format: &str) -> Result<(String, usize)> {
        let mut repeat_str = String::new();
        let mut type_str = String::new();
        let mut parsing_repeat = true;

        for ch in format.chars() {
            if ch.is_ascii_digit() && parsing_repeat {
                repeat_str.push(ch);
            } else {
                parsing_repeat = false;
                type_str.push(ch);
            }
        }

        if type_str.is_empty() {
            return Err(FitsError::InvalidFormat(format!(
                "Invalid FITS format '{}' - missing data type",
                format
            )));
        }

        let repeat = if repeat_str.is_empty() {
            1
        } else {
            repeat_str.parse().map_err(|_| {
                FitsError::InvalidFormat(format!("Invalid repeat count in format '{}'", format))
            })?
        };

        Ok((type_str, repeat))
    }

    pub(super) fn get_element_size(&self, data_type: &str) -> Result<usize> {
        match data_type.chars().next().unwrap_or('X') {
            'L' => Ok(1),
            'X' => Ok(1),
            'B' => Ok(1),
            'I' => Ok(2),
            'J' => Ok(4),
            'K' => Ok(8),
            'A' => Ok(1),
            'E' => Ok(4),
            'D' => Ok(8),
            'C' => Ok(8),
            'M' => Ok(16),
            'P' => Ok(8),
            'Q' => Ok(16),
            _ => Err(FitsError::InvalidFormat(format!(
                "Unknown binary table format: {}",
                data_type
            ))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::Header;
    use crate::fits::io::reader::HduInfo;

    fn dummy_hdu() -> BinaryTableHdu {
        let info = HduInfo {
            index: 0,
            header_start: 0,
            header_size: 0,
            data_start: 0,
            data_size: 0,
        };
        BinaryTableHdu::new(Header::new(), info)
    }

    #[test]
    fn parse_binary_format_empty_is_rejected() {
        let hdu = dummy_hdu();
        let err = hdu.parse_binary_format("").unwrap_err();
        assert!(matches!(err, FitsError::InvalidFormat(_)));
        assert!(err.to_string().contains("Empty column format"));
    }

    #[test]
    fn parse_standard_format_single_char_types_have_repeat_one() {
        let hdu = dummy_hdu();
        for ch in ['L', 'X', 'B', 'I', 'J', 'K', 'A', 'E', 'D', 'C', 'M'] {
            let input = ch.to_string();
            let (ty, rep) = hdu.parse_binary_format(&input).unwrap();
            assert_eq!(ty, input);
            assert_eq!(rep, 1);
        }
    }

    #[test]
    fn parse_standard_format_with_repeat_prefix() {
        let hdu = dummy_hdu();
        let (ty, rep) = hdu.parse_binary_format("5J").unwrap();
        assert_eq!(ty, "J");
        assert_eq!(rep, 5);

        let (ty, rep) = hdu.parse_binary_format("128A").unwrap();
        assert_eq!(ty, "A");
        assert_eq!(rep, 128);
    }

    #[test]
    fn parse_standard_format_preserves_trailing_characters_in_type() {
        let hdu = dummy_hdu();
        let (ty, rep) = hdu.parse_binary_format("10Ab").unwrap();
        assert_eq!(ty, "Ab");
        assert_eq!(rep, 10);
    }

    #[test]
    fn parse_standard_format_zero_repeat_is_allowed() {
        let hdu = dummy_hdu();
        let (ty, rep) = hdu.parse_binary_format("0J").unwrap();
        assert_eq!(ty, "J");
        assert_eq!(rep, 0);
    }

    #[test]
    fn parse_standard_format_only_digits_is_rejected() {
        let hdu = dummy_hdu();
        let err = hdu.parse_binary_format("42").unwrap_err();
        assert!(matches!(err, FitsError::InvalidFormat(_)));
        assert!(err.to_string().contains("missing data type"));
    }

    #[test]
    fn parse_standard_format_overflowing_repeat_is_rejected() {
        let hdu = dummy_hdu();
        let huge = format!("{}J", "9".repeat(40));
        let err = hdu.parse_binary_format(&huge).unwrap_err();
        assert!(matches!(err, FitsError::InvalidFormat(_)));
        assert!(err.to_string().contains("Invalid repeat count"));
    }

    #[test]
    fn parse_variable_length_format_p_without_repeat() {
        let hdu = dummy_hdu();
        let (ty, rep) = hdu.parse_binary_format("PJ").unwrap();
        assert_eq!(ty, "PJ");
        assert_eq!(rep, 1);
    }

    #[test]
    fn parse_variable_length_format_q_without_repeat() {
        let hdu = dummy_hdu();
        let (ty, rep) = hdu.parse_binary_format("QE").unwrap();
        assert_eq!(ty, "QE");
        assert_eq!(rep, 1);
    }

    #[test]
    fn parse_variable_length_format_with_repeat_prefix() {
        let hdu = dummy_hdu();
        let (ty, rep) = hdu.parse_binary_format("7PJ").unwrap();
        assert_eq!(ty, "PJ");
        assert_eq!(rep, 7);

        let (ty, rep) = hdu.parse_binary_format("100QE").unwrap();
        assert_eq!(ty, "QE");
        assert_eq!(rep, 100);
    }

    #[test]
    fn parse_variable_length_format_all_element_types() {
        let hdu = dummy_hdu();
        for descriptor in ['P', 'Q'] {
            for data in ['L', 'X', 'B', 'I', 'J', 'K', 'A', 'E', 'D', 'C', 'M'] {
                let input = format!("{descriptor}{data}");
                let (ty, rep) = hdu.parse_binary_format(&input).unwrap();
                assert_eq!(ty, input);
                assert_eq!(rep, 1);
            }
        }
    }

    #[test]
    fn parse_variable_length_format_missing_descriptor_is_rejected() {
        // Input contains 'P' or 'Q' (triggering the variable-length path) but it's not
        // positioned after the optional digit prefix.
        let hdu = dummy_hdu();
        let err = hdu.parse_binary_format("AP").unwrap_err();
        assert!(matches!(err, FitsError::InvalidFormat(_)));
        assert!(err.to_string().contains("must contain P or Q"));
    }

    #[test]
    fn parse_variable_length_format_missing_data_type_is_rejected() {
        let hdu = dummy_hdu();
        let err = hdu.parse_binary_format("P").unwrap_err();
        assert!(matches!(err, FitsError::InvalidFormat(_)));
        assert!(err.to_string().contains("Missing data type"));

        let err = hdu.parse_binary_format("10Q").unwrap_err();
        assert!(matches!(err, FitsError::InvalidFormat(_)));
        assert!(err.to_string().contains("Missing data type"));
    }

    #[test]
    fn parse_variable_length_format_extra_chars_after_data_type_are_ignored() {
        // Only the character immediately after P/Q is used as the data type;
        // anything past that is discarded by the current implementation.
        let hdu = dummy_hdu();
        let (ty, rep) = hdu.parse_binary_format("PJextra").unwrap();
        assert_eq!(ty, "PJ");
        assert_eq!(rep, 1);
    }

    #[test]
    fn parse_binary_format_dispatches_on_p_or_q_anywhere_in_string() {
        // The dispatcher uses `contains('P') || contains('Q')` — so even with P in a
        // non-leading position the variable-length path is chosen.
        let hdu = dummy_hdu();
        let err = hdu.parse_binary_format("JP").unwrap_err();
        assert!(err.to_string().contains("must contain P or Q"));
    }

    #[test]
    fn get_element_size_covers_every_known_type() {
        let hdu = dummy_hdu();
        let cases = [
            ("L", 1usize),
            ("X", 1),
            ("B", 1),
            ("I", 2),
            ("J", 4),
            ("K", 8),
            ("A", 1),
            ("E", 4),
            ("D", 8),
            ("C", 8),
            ("M", 16),
            ("P", 8),
            ("Q", 16),
        ];
        for (ty, expected) in cases {
            assert_eq!(hdu.get_element_size(ty).unwrap(), expected, "type={ty}");
        }
    }

    #[test]
    fn get_element_size_uses_only_first_character() {
        let hdu = dummy_hdu();
        assert_eq!(hdu.get_element_size("Jfoo").unwrap(), 4);
        assert_eq!(hdu.get_element_size("PJ").unwrap(), 8);
        assert_eq!(hdu.get_element_size("QE").unwrap(), 16);
    }

    #[test]
    fn get_element_size_empty_string_defaults_to_x() {
        let hdu = dummy_hdu();
        assert_eq!(hdu.get_element_size("").unwrap(), 1);
    }

    #[test]
    fn get_element_size_unknown_type_is_rejected() {
        let hdu = dummy_hdu();
        let err = hdu.get_element_size("Z").unwrap_err();
        assert!(matches!(err, FitsError::InvalidFormat(_)));
        assert!(err.to_string().contains("Unknown binary table format"));
        assert!(err.to_string().contains('Z'));
    }
}
