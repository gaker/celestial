const FITS_BLOCK_SIZE: usize = 2880;

pub fn ones_complement_sum(data: &[u8]) -> u32 {
    let mut sum: u64 = 0;
    let mut i = 0;

    while i + 4 <= data.len() {
        let word = u32::from_be_bytes([data[i], data[i + 1], data[i + 2], data[i + 3]]);
        sum += word as u64;
        i += 4;
    }

    if i < data.len() {
        let mut last_word = [0u8; 4];
        last_word[..data.len() - i].copy_from_slice(&data[i..]);
        sum += u32::from_be_bytes(last_word) as u64;
    }

    fold_to_ones_complement(sum)
}

fn fold_to_ones_complement(mut sum: u64) -> u32 {
    while sum > 0xFFFF_FFFF {
        sum = (sum & 0xFFFF_FFFF) + (sum >> 32);
    }
    sum as u32
}

pub fn calculate_datasum(data: &[u8]) -> u32 {
    if data.is_empty() {
        return 0;
    }
    ones_complement_sum(data)
}

pub fn verify_datasum(data: &[u8], expected_str: &str) -> bool {
    let expected = match expected_str.trim().parse::<u32>() {
        Ok(v) => v,
        Err(_) => return false,
    };
    calculate_datasum(data) == expected
}

pub fn encode_checksum(sum: u32) -> String {
    let bytes = sum.to_be_bytes();
    let mut encoded = [b'0'; 16];

    for (i, &byte) in bytes.iter().enumerate() {
        let quotient = byte / 4;
        let remainder = byte % 4;
        let ch = [
            quotient + remainder,
            quotient + (4 - remainder),
            quotient,
            quotient + 4,
        ];

        for (j, &c) in ch.iter().enumerate() {
            let pos = (4 * j + i) % 16;
            encoded[pos] = ascii_encode(c + b'0');
        }
    }

    String::from_utf8(encoded.to_vec()).unwrap()
}

fn ascii_encode(val: u8) -> u8 {
    if val > b'Z' {
        val + 10
    } else if val > b'9' {
        val + 7
    } else {
        val
    }
}

pub fn decode_checksum(encoded: &str) -> Option<u32> {
    if encoded.len() != 16 {
        return None;
    }

    let chars: &[u8] = encoded.as_bytes();
    let mut bytes = [0u8; 4];

    for (i, byte) in bytes.iter_mut().enumerate() {
        let mut ch = [0u8; 4];
        for (j, c) in ch.iter_mut().enumerate() {
            let pos = (4 * j + i) % 16;
            *c = ascii_decode(chars[pos])? - b'0';
        }

        let quotient = ch[2];
        let remainder = ch[0].wrapping_sub(quotient);
        *byte = quotient * 4 + remainder;
    }

    Some(u32::from_be_bytes(bytes))
}

fn ascii_decode(c: u8) -> Option<u8> {
    if c >= b'a' {
        Some(c - 10)
    } else if c >= b'A' {
        Some(c - 7)
    } else if c.is_ascii_digit() {
        Some(c)
    } else {
        None
    }
}

pub fn calculate_hdu_checksum(header_bytes: &[u8], data_sum: u32) -> u32 {
    let header_sum = ones_complement_sum(header_bytes);
    add_ones_complement(header_sum, data_sum)
}

fn add_ones_complement(a: u32, b: u32) -> u32 {
    let sum = a as u64 + b as u64;
    fold_to_ones_complement(sum)
}

pub fn checksum_complement(sum: u32) -> u32 {
    !sum
}

pub fn verify_hdu_checksum(header_bytes: &[u8], data_bytes: &[u8]) -> bool {
    let total_sum = ones_complement_sum(header_bytes);
    let data_sum = ones_complement_sum(data_bytes);
    let combined = add_ones_complement(total_sum, data_sum);
    combined == 0xFFFF_FFFF || combined == 0
}

pub fn format_datasum(sum: u32) -> String {
    sum.to_string()
}

pub fn create_checksum_card_value(header_bytes: &[u8], data_sum: u32) -> String {
    let header_sum = ones_complement_sum(header_bytes);
    let combined = add_ones_complement(header_sum, data_sum);
    let complement = checksum_complement(combined);
    encode_checksum(complement)
}

pub fn pad_to_block(data: &[u8]) -> Vec<u8> {
    let remainder = data.len() % FITS_BLOCK_SIZE;
    if remainder == 0 {
        return data.to_vec();
    }

    let padding_needed = FITS_BLOCK_SIZE - remainder;
    let mut padded = Vec::with_capacity(data.len() + padding_needed);
    padded.extend_from_slice(data);
    padded.resize(padded.len() + padding_needed, 0);
    padded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ones_complement_sum_empty() {
        assert_eq!(ones_complement_sum(&[]), 0);
    }

    #[test]
    fn ones_complement_sum_single_word() {
        let data = [0x00, 0x00, 0x00, 0x01];
        assert_eq!(ones_complement_sum(&data), 1);
    }

    #[test]
    fn ones_complement_sum_multiple_words() {
        let data = [0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02];
        assert_eq!(ones_complement_sum(&data), 3);
    }

    #[test]
    fn ones_complement_sum_with_carry() {
        let data = [0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00, 0x02];
        let sum = ones_complement_sum(&data);
        assert_eq!(sum, 2);
    }

    #[test]
    fn ones_complement_sum_partial_word() {
        let data = [0x00, 0x00, 0x00, 0x01, 0x02];
        let sum = ones_complement_sum(&data);
        assert_eq!(sum, 1 + 0x02000000);
    }

    #[test]
    fn encode_decode_roundtrip() {
        let test_values = [0u32, 1, 0xFFFF_FFFF, 0x12345678, 0xDEADBEEF];

        for value in test_values {
            let encoded = encode_checksum(value);
            assert_eq!(encoded.len(), 16);

            let decoded = decode_checksum(&encoded);
            assert_eq!(decoded, Some(value));
        }
    }

    #[test]
    fn decode_invalid_length() {
        assert_eq!(decode_checksum("short"), None);
        assert_eq!(decode_checksum("toolongstringhere!"), None);
    }

    #[test]
    fn decode_invalid_chars() {
        assert_eq!(decode_checksum("!@#$%^&*(){}[]<>"), None);
    }

    #[test]
    fn calculate_datasum_empty() {
        assert_eq!(calculate_datasum(&[]), 0);
    }

    #[test]
    fn calculate_datasum_nonzero() {
        let data = vec![0x01; 100];
        assert_ne!(calculate_datasum(&data), 0);
    }

    #[test]
    fn verify_datasum_valid() {
        let data = [0x00, 0x00, 0x00, 0x05];
        let sum = calculate_datasum(&data);
        assert!(verify_datasum(&data, &sum.to_string()));
    }

    #[test]
    fn verify_datasum_invalid() {
        let data = [0x00, 0x00, 0x00, 0x05];
        assert!(!verify_datasum(&data, "12345"));
    }

    #[test]
    fn verify_datasum_bad_parse() {
        let data = [0x00, 0x00, 0x00, 0x05];
        assert!(!verify_datasum(&data, "not_a_number"));
    }

    #[test]
    fn format_datasum_values() {
        assert_eq!(format_datasum(0), "0");
        assert_eq!(format_datasum(12345), "12345");
        assert_eq!(format_datasum(u32::MAX), "4294967295");
    }

    #[test]
    fn checksum_complement_properties() {
        assert_eq!(checksum_complement(0), 0xFFFF_FFFF);
        assert_eq!(checksum_complement(0xFFFF_FFFF), 0);
    }

    #[test]
    fn add_ones_complement_no_carry() {
        assert_eq!(add_ones_complement(1, 2), 3);
    }

    #[test]
    fn add_ones_complement_with_carry() {
        assert_eq!(add_ones_complement(0xFFFF_FFFF, 2), 2);
    }

    #[test]
    fn pad_to_block_already_aligned() {
        let data = vec![0u8; FITS_BLOCK_SIZE];
        let padded = pad_to_block(&data);
        assert_eq!(padded.len(), FITS_BLOCK_SIZE);
    }

    #[test]
    fn pad_to_block_needs_padding() {
        let data = vec![0u8; 100];
        let padded = pad_to_block(&data);
        assert_eq!(padded.len(), FITS_BLOCK_SIZE);
        assert_eq!(padded.len() % FITS_BLOCK_SIZE, 0);
    }

    #[test]
    fn pad_to_block_empty() {
        let data: Vec<u8> = vec![];
        let padded = pad_to_block(&data);
        assert!(padded.is_empty());
    }

    #[test]
    fn verify_hdu_checksum_all_zeros() {
        let header = vec![0u8; FITS_BLOCK_SIZE];
        let data = vec![0u8; FITS_BLOCK_SIZE];
        assert!(verify_hdu_checksum(&header, &data));
    }

    #[test]
    fn calculate_hdu_checksum_combines_correctly() {
        let header = [0x00, 0x00, 0x00, 0x05];
        let data_sum = 10u32;
        let result = calculate_hdu_checksum(&header, data_sum);
        assert_eq!(result, 15);
    }

    #[test]
    fn encoding_chars_valid_ascii() {
        for i in 0..=63u8 {
            let sum = (i as u32) << 28;
            let encoded = encode_checksum(sum);
            assert!(encoded.is_ascii());
        }
    }
}
