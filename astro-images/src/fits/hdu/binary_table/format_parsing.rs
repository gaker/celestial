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
