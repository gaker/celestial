// FITS 4.0 §3.1: a file is a sequence of 2880-byte blocks, and a header is
// a sequence of 80-byte cards.
pub(crate) const FITS_BLOCK_SIZE: usize = 2880;
pub(crate) const CARD_SIZE: usize = 80;

// Pads `bytes` to a whole number of blocks. FITS 4.0 §3.3.2 and §3.4 fill
// a header with spaces and a data unit with zeros.
pub(crate) fn pad_to_block(bytes: &mut Vec<u8>, fill: u8) {
    bytes.resize(
        bytes.len().div_ceil(FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE,
        fill,
    );
}

// The start of the block after the one holding byte `position - 1`.
pub(crate) fn align_to_block(position: u64) -> u64 {
    let block = FITS_BLOCK_SIZE as u64;
    position.div_ceil(block) * block
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pad_to_block_fills_to_the_next_block() {
        let mut bytes = vec![1u8; 3];
        pad_to_block(&mut bytes, b' ');
        assert_eq!(bytes.len(), FITS_BLOCK_SIZE);
        assert_eq!(&bytes[..4], b"\x01\x01\x01 ");
        assert_eq!(bytes[FITS_BLOCK_SIZE - 1], b' ');
    }

    #[test]
    fn pad_to_block_keeps_whole_blocks_and_empty_input() {
        let mut bytes = vec![7u8; FITS_BLOCK_SIZE];
        pad_to_block(&mut bytes, 0);
        assert_eq!(bytes, vec![7u8; FITS_BLOCK_SIZE]);

        let mut empty = Vec::new();
        pad_to_block(&mut empty, 0);
        assert!(empty.is_empty());
    }

    #[test]
    fn align_to_block_rounds_up() {
        let cases = [(0, 0), (1, 2880), (2880, 2880), (2881, 5760), (5761, 8640)];
        for (position, aligned) in cases {
            assert_eq!(align_to_block(position), aligned);
        }
    }
}
