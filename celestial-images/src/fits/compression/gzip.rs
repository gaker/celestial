use crate::fits::errors::{FitsError, Result};
use flate2::read::GzDecoder;
use std::io::Read;

// DEFLATE expands at most 1032 to 1, so a tile larger than this cannot come
// from `len` compressed bytes.
pub(super) fn max_inflated(len: usize) -> usize {
    len.saturating_mul(1032)
}

// Inflates a GZIP_1 stream that must hold exactly `expected` bytes. Reading
// stops one byte past `expected`, so a corrupt stream cannot grow the buffer.
pub(super) fn inflate(bytes: &[u8], expected: usize) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    out.try_reserve_exact(expected)
        .map_err(|_| FitsError::InvalidFormat(format!("cannot allocate a {expected}-byte tile")))?;
    let limit = u64::try_from(expected)
        .ok()
        .and_then(|n| n.checked_add(1))
        .unwrap_or(u64::MAX);
    GzDecoder::new(bytes)
        .take(limit)
        .read_to_end(&mut out)
        .map_err(|e| FitsError::InvalidFormat(format!("GZIP tile does not inflate: {e}")))?;
    if out.len() != expected {
        return Err(FitsError::InvalidFormat(format!(
            "GZIP tile inflates to {} or more bytes where {} are expected",
            out.len(),
            expected
        )));
    }
    Ok(out)
}

// GZIP_2 (FITS 4.0 §10.4.2) stores the most significant byte of every value
// first, then the next byte of every value, and so on.
pub(super) fn unshuffle(bytes: &[u8], width: usize) -> Vec<u8> {
    let count = bytes.len() / width.max(1);
    if width <= 1 || count == 0 {
        return bytes.to_vec();
    }
    let mut out = vec![0u8; bytes.len()];
    for (byte, plane) in bytes.chunks_exact(count).take(width).enumerate() {
        for (value, &b) in plane.iter().enumerate() {
            out[value * width + byte] = b;
        }
    }
    out
}

// The inverse of `unshuffle`: byte `byte` of every value goes to plane `byte`.
pub(super) fn shuffle(bytes: &[u8], width: usize) -> Vec<u8> {
    let count = bytes.len() / width.max(1);
    if width <= 1 || count == 0 {
        return bytes.to_vec();
    }
    let mut out = vec![0u8; bytes.len()];
    for (value, chunk) in bytes.chunks_exact(width).enumerate() {
        for (byte, &b) in chunk.iter().enumerate() {
            out[byte * count + value] = b;
        }
    }
    out
}

pub(super) fn compress(data: &[u8]) -> Result<Vec<u8>> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder
        .write_all(data)
        .map_err(|e| FitsError::InvalidFormat(format!("GZIP compression failed: {}", e)))?;
    encoder
        .finish()
        .map_err(|e| FitsError::InvalidFormat(format!("GZIP compression failed: {}", e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inflate_returns_the_tile() {
        let data: Vec<u8> = (0..100).collect();
        let stream = compress(&data).unwrap();
        assert_eq!(inflate(&stream, 100).unwrap(), data);
    }

    #[test]
    fn inflate_rejects_short_and_long_tiles() {
        let stream = compress(&[7u8; 50]).unwrap();
        assert!(matches!(
            inflate(&stream, 100),
            Err(FitsError::InvalidFormat(_))
        ));
        assert!(matches!(
            inflate(&stream, 49),
            Err(FitsError::InvalidFormat(_))
        ));
    }

    #[test]
    fn inflate_rejects_bytes_that_are_not_gzip() {
        assert!(inflate(&[1, 2, 3, 4], 4).is_err());
    }

    // The example of FITS 4.0 §10.4.2: A1 A2 B1 B2 ... is stored as
    // A1 B1 C1 D1 E1 A2 B2 C2 D2 E2.
    #[test]
    fn unshuffle_reverses_the_standard_example() {
        let stored = [0xA1, 0xB1, 0xC1, 0xD1, 0xE1, 0xA2, 0xB2, 0xC2, 0xD2, 0xE2];
        assert_eq!(
            unshuffle(&stored, 2),
            vec![0xA1, 0xA2, 0xB1, 0xB2, 0xC1, 0xC2, 0xD1, 0xD2, 0xE1, 0xE2]
        );
        assert_eq!(unshuffle(&stored, 1), stored.to_vec());
        assert!(unshuffle(&[], 4).is_empty());
    }

    #[test]
    fn shuffle_writes_the_standard_example() {
        let values = [0xA1, 0xA2, 0xB1, 0xB2, 0xC1, 0xC2, 0xD1, 0xD2, 0xE1, 0xE2];
        assert_eq!(
            shuffle(&values, 2),
            vec![0xA1, 0xB1, 0xC1, 0xD1, 0xE1, 0xA2, 0xB2, 0xC2, 0xD2, 0xE2]
        );
        let wide: Vec<u8> = (0..24).collect();
        for width in [1, 4, 8] {
            assert_eq!(unshuffle(&shuffle(&wide, width), width), wide);
        }
    }

    #[test]
    fn max_inflated_bounds_a_constant_stream() {
        let data = vec![0u8; 1 << 20];
        let stream = compress(&data).unwrap();
        assert!(max_inflated(stream.len()) >= data.len());
    }
}
