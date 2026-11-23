pub(crate) fn ones_complement_sum(data: &[u8]) -> u32 {
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

pub(crate) fn calculate_datasum(data: &[u8]) -> u32 {
    if data.is_empty() {
        return 0;
    }
    ones_complement_sum(data)
}

pub(crate) fn verify_datasum(data: &[u8], expected_str: &str) -> bool {
    let expected = match expected_str.trim().parse::<u32>() {
        Ok(v) => v,
        Err(_) => return false,
    };
    calculate_datasum(data) == expected
}

// Encoding from the FITS checksum convention: each byte becomes four
// characters offset by '0' whose sum is the byte, then the 16 characters
// rotate right by one so they align with 32-bit words at card column 12.
const CHECKSUM_OFFSET: u8 = 0x30;

pub(crate) fn encode_checksum(sum: u32) -> String {
    let mut unrotated = [0u8; 16];
    for (i, byte) in sum.to_be_bytes().into_iter().enumerate() {
        let quotient = byte / 4 + CHECKSUM_OFFSET;
        let mut ch = [quotient; 4];
        ch[0] += byte % 4;
        for pair in ch.as_chunks_mut::<2>().0 {
            // Moving one unit between the pair keeps the sum and avoids the
            // punctuation between '9' and 'A' and between 'Z' and 'a'.
            while is_checksum_punctuation(pair[0]) || is_checksum_punctuation(pair[1]) {
                pair[0] += 1;
                pair[1] -= 1;
            }
        }
        for (j, &c) in ch.iter().enumerate() {
            unrotated[4 * j + i] = c;
        }
    }
    (0..16).map(|i| unrotated[(i + 15) % 16] as char).collect()
}

fn is_checksum_punctuation(c: u8) -> bool {
    (b':'..=b'@').contains(&c) || (b'['..=b'`').contains(&c)
}

fn add_ones_complement(a: u32, b: u32) -> u32 {
    let sum = a as u64 + b as u64;
    fold_to_ones_complement(sum)
}

pub(crate) fn checksum_complement(sum: u32) -> u32 {
    !sum
}

pub(crate) fn verify_hdu_checksum(header_bytes: &[u8], data_bytes: &[u8]) -> bool {
    let total_sum = ones_complement_sum(header_bytes);
    let data_sum = ones_complement_sum(data_bytes);
    let combined = add_ones_complement(total_sum, data_sum);
    combined == 0xFFFF_FFFF || combined == 0
}

pub(crate) fn format_datasum(sum: u32) -> String {
    sum.to_string()
}

pub(crate) fn create_checksum_card_value(header_bytes: &[u8], data_sum: u32) -> String {
    let header_sum = ones_complement_sum(header_bytes);
    let combined = add_ones_complement(header_sum, data_sum);
    let complement = checksum_complement(combined);
    encode_checksum(complement)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::block::FITS_BLOCK_SIZE;

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
    fn verify_hdu_checksum_all_zeros() {
        let header = vec![0u8; FITS_BLOCK_SIZE];
        let data = vec![0u8; FITS_BLOCK_SIZE];
        assert!(verify_hdu_checksum(&header, &data));
    }

    // Values and CHECKSUM strings from files written by astropy 6.0.1 with
    // checksum=True. Each value is the complement of the HDU sum taken with
    // the CHECKSUM string replaced by sixteen '0' characters.
    const ASTROPY_CHECKSUMS: [(u32, &str); 4] = [
        (0x062F_CDF7, "m3Adp18cm1Acm15c"),
        (0xA87E_9F12, "4ZQZ6ZOW4ZOW4ZOW"),
        (0xC4B1_053B, "9aa2FaX19aa1CaW1"),
        (0xCA29_EDF0, "ldAllb4klbAklb3k"),
    ];

    #[test]
    fn encode_checksum_matches_astropy() {
        for (value, expected) in ASTROPY_CHECKSUMS {
            assert_eq!(encode_checksum(value), expected, "0x{:08X}", value);
        }
    }

    #[test]
    fn encoded_checksum_has_no_punctuation() {
        for value in [
            0u32,
            0xFFFF_FFFF,
            0x0A0A_0A0A,
            0x2B2B_2B2B,
            0xEDED_EDED,
            0x1234_5678,
        ] {
            let encoded = encode_checksum(value);
            assert!(
                encoded.bytes().all(|b| b.is_ascii_alphanumeric()),
                "0x{:08X} encoded as {}",
                value,
                encoded
            );
        }
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
