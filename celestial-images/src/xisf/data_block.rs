use crate::core::sample_order::to_planar;
use crate::core::types::ByteOrder;
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::compression::XisfCompression;
use crate::xisf::header::image_info::ImageInfo;
use crate::xisf::header::location::DataLocation;
use crate::xisf::header::sample_format::SampleFormat;
use flate2::read::ZlibDecoder;
use std::io::{ErrorKind, Read, Seek, SeekFrom};

// An LZ4 sequence writes at most 255 output bytes per input byte, through a
// match-length extension byte. A block that claims more cannot be valid.
const LZ4_MAX_EXPANSION: u64 = 255;

// Reads an image data block and returns its samples as little-endian bytes in
// the stored pixel order. The result length always equals the geometry size.
pub(super) fn read_image_block<R: Read + Seek>(
    reader: &mut R,
    stream_len: u64,
    info: &ImageInfo,
) -> Result<Vec<u8>> {
    if let Some(defect) = &info.block_defect {
        return Err(defect.to_error());
    }
    let expected = info.checked_data_size()?;
    let stored = read_stored(reader, stream_len, &info.location)?;
    // XISF 1.0 section 10.6.2: a block that fails its checksum is not
    // decompressed.
    if let Some(checksum) = &info.checksum {
        checksum.verify(&stored)?;
    }
    let mut data = decode(info, stored, expected)?;
    if let Some(item_size) = info.shuffle_item_size {
        data = unshuffle(&data, item_size);
    }
    if info.byte_order == ByteOrder::BigEndian {
        swap_to_little_endian(&mut data, info.sample_format.bytes_per_sample());
    }
    Ok(data)
}

pub(super) fn past_end_of_file(what: &str, stream_len: u64) -> XisfError {
    XisfError::Io(std::io::Error::new(
        ErrorKind::UnexpectedEof,
        format!("{what} extends past the end of the {stream_len}-byte file"),
    ))
}

// The bounds check runs before the allocation, so a forged size cannot
// request more memory than the file holds.
pub(super) fn read_stored<R: Read + Seek>(
    reader: &mut R,
    stream_len: u64,
    location: &DataLocation,
) -> Result<Vec<u8>> {
    if let Some(bytes) = location.header_bytes() {
        return Ok(bytes.to_vec());
    }
    if !location.ends_within(stream_len) {
        return Err(past_end_of_file("The data block", stream_len));
    }
    let len = usize::try_from(location.size)
        .map_err(|_| past_end_of_file("The data block", stream_len))?;
    reader.seek(SeekFrom::Start(location.offset))?;
    let mut buffer = vec![0u8; len];
    reader.read_exact(&mut buffer)?;
    Ok(buffer)
}

fn decode(info: &ImageInfo, stored: Vec<u8>, expected: usize) -> Result<Vec<u8>> {
    if !info.compression.is_compressed() {
        check_length(stored.len() as u64, expected, "The stored data block")?;
        return Ok(stored);
    }
    let declared = info.uncompressed_size.ok_or_else(|| {
        XisfError::InvalidFormat("The compression attribute has no uncompressed size".into())
    })?;
    check_length(declared, expected, "The declared uncompressed size")?;
    // Without a subblocks attribute the block is a single subblock (section 10.6.1).
    let whole_block = [(stored.len() as u64, declared)];
    let subblocks = match info.subblocks.as_slice() {
        [] => &whole_block[..],
        listed => listed,
    };
    decompress_subblocks(info.compression, &stored, subblocks, expected)
}

fn check_length(actual: u64, expected: usize, what: &str) -> Result<()> {
    if actual == expected as u64 {
        return Ok(());
    }
    Err(XisfError::InvalidFormat(format!(
        "{what} is {actual} bytes; the image geometry needs {expected}"
    )))
}

fn decompress_subblocks(
    codec: XisfCompression,
    stored: &[u8],
    subblocks: &[(u64, u64)],
    expected: usize,
) -> Result<Vec<u8>> {
    check_subblock_totals(stored.len(), subblocks, expected)?;
    let mut output = Vec::new();
    let mut rest = stored;
    for &(compressed, uncompressed) in subblocks {
        let (chunk, tail) = rest
            .split_at_checked(compressed as usize)
            .ok_or_else(|| subblock_error("lengths exceed the data block"))?;
        rest = tail;
        let decoded = decompress_subblock(codec, chunk, uncompressed as usize)?;
        if output.is_empty() {
            output = decoded;
        } else {
            output.extend_from_slice(&decoded);
        }
    }
    Ok(output)
}

// Both sums are bounded by in-memory lengths once they match, so the casts
// to usize in `decompress_subblocks` cannot truncate.
fn check_subblock_totals(stored: usize, subblocks: &[(u64, u64)], expected: usize) -> Result<()> {
    let total = |pick: fn(&(u64, u64)) -> u64| {
        subblocks
            .iter()
            .try_fold(0u64, |acc, pair| acc.checked_add(pick(pair)))
    };
    if total(|pair| pair.0) != Some(stored as u64) {
        return Err(subblock_error(
            "compressed lengths do not sum to the block size",
        ));
    }
    if total(|pair| pair.1) != Some(expected as u64) {
        return Err(subblock_error(
            "uncompressed lengths do not sum to the image size",
        ));
    }
    Ok(())
}

fn subblock_error(problem: &str) -> XisfError {
    XisfError::InvalidFormat(format!("Compression subblock {problem}"))
}

// Equal lengths mark an uncompressed subblock, which is copied unmodified
// (XISF 1.0 section 10.6.1).
fn decompress_subblock(codec: XisfCompression, chunk: &[u8], size: usize) -> Result<Vec<u8>> {
    if chunk.len() == size {
        return Ok(chunk.to_vec());
    }
    let decoded = match codec {
        XisfCompression::None => chunk.to_vec(),
        XisfCompression::Zlib => inflate(chunk, size)?,
        XisfCompression::Lz4 | XisfCompression::Lz4Hc => lz4_block(chunk, size)?,
        XisfCompression::Zstd => return Err(XisfError::UnsupportedCodec("zstd".into())),
    };
    if decoded.len() != size {
        return Err(subblock_error(&format!(
            "decompresses to {} bytes instead of {size}",
            decoded.len()
        )));
    }
    Ok(decoded)
}

// Reading one byte past `size` detects a stream that is too long without
// inflating all of it.
fn inflate(chunk: &[u8], size: usize) -> Result<Vec<u8>> {
    let mut output = Vec::new();
    ZlibDecoder::new(chunk)
        .take((size as u64).saturating_add(1))
        .read_to_end(&mut output)
        .map_err(|e| XisfError::InvalidFormat(format!("zlib decompression failed: {e}")))?;
    Ok(output)
}

// XISF stores raw LZ4 blocks without frame headers (section 10.6.6).
fn lz4_block(chunk: &[u8], size: usize) -> Result<Vec<u8>> {
    if size as u64 > (chunk.len() as u64).saturating_mul(LZ4_MAX_EXPANSION) {
        return Err(XisfError::InvalidFormat(format!(
            "A {}-byte LZ4 block cannot hold {size} bytes",
            chunk.len()
        )));
    }
    lz4_flex::block::decompress(chunk, size)
        .map_err(|e| XisfError::InvalidFormat(format!("LZ4 decompression failed: {e}")))
}

// The reverse transform of XISF 1.0 section 10.6.3. Byte j of item i is
// stored at j * items + i. Trailing bytes that do not fill an item are
// stored unshuffled after the shuffled bytes.
fn unshuffle(shuffled: &[u8], item_size: usize) -> Vec<u8> {
    let items = shuffled.len() / item_size;
    let mut output = shuffled.to_vec();
    if items == 0 {
        return output;
    }
    let planes = shuffled.chunks_exact(items).take(item_size);
    for (j, plane) in planes.enumerate() {
        for (item, &byte) in output.chunks_exact_mut(item_size).zip(plane) {
            item[j] = byte;
        }
    }
    output
}

fn swap_to_little_endian(data: &mut [u8], bytes_per_sample: usize) {
    for sample in data.chunks_exact_mut(bytes_per_sample) {
        sample.reverse();
    }
}

// Reorders the samples of an image with Normal storage to planar order.
pub(super) fn planar_samples(data: &[u8], info: &ImageInfo) -> Vec<u8> {
    let channels = info.num_channels();
    match info.sample_format {
        SampleFormat::UInt8 => to_planar(data, channels),
        SampleFormat::UInt16 => planar_items::<2>(data, channels),
        SampleFormat::UInt32 | SampleFormat::Float32 => planar_items::<4>(data, channels),
        SampleFormat::UInt64 | SampleFormat::Float64 => planar_items::<8>(data, channels),
    }
}

// Moves each N-byte sample as one item.
fn planar_items<const N: usize>(data: &[u8], channels: usize) -> Vec<u8> {
    to_planar(data.as_chunks::<N>().0, channels).into_flattened()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unshuffle_reverses_the_spec_example_layout() {
        // Items (a0 a1 a2)(b0 b1 b2) are stored as a0 b0 a1 b1 a2 b2.
        let shuffled = [0xA0, 0xB0, 0xA1, 0xB1, 0xA2, 0xB2];
        assert_eq!(
            unshuffle(&shuffled, 3),
            vec![0xA0, 0xA1, 0xA2, 0xB0, 0xB1, 0xB2]
        );
    }

    #[test]
    fn unshuffle_leaves_the_partial_trailing_item() {
        // Three 2-byte items and one trailing byte D0.
        let shuffled = [0xA0, 0xB0, 0xC0, 0xA1, 0xB1, 0xC1, 0xD0];
        assert_eq!(
            unshuffle(&shuffled, 2),
            vec![0xA0, 0xA1, 0xB0, 0xB1, 0xC0, 0xC1, 0xD0]
        );
    }

    #[test]
    fn unshuffle_is_identity_for_one_byte_items_and_short_blocks() {
        let data = [5u8, 6, 7];
        assert_eq!(unshuffle(&data, 1), data);
        assert_eq!(unshuffle(&data, 8), data);
        assert_eq!(unshuffle(&[], 4), Vec::<u8>::new());
    }

    #[test]
    fn swap_to_little_endian_reverses_each_sample() {
        let mut data = [1u8, 2, 3, 4, 5, 6, 7, 8];
        swap_to_little_endian(&mut data, 4);
        assert_eq!(data, [4, 3, 2, 1, 8, 7, 6, 5]);
        swap_to_little_endian(&mut data, 1);
        assert_eq!(data, [4, 3, 2, 1, 8, 7, 6, 5]);
    }

    #[test]
    fn lz4_block_rejects_impossible_expansion() {
        assert!(lz4_block(&[0x10, 7], 511).is_err());
        assert_eq!(lz4_block(&[0x10, 7], 1).unwrap(), vec![7]);
    }

    #[test]
    fn inflate_stops_one_byte_past_the_declared_size() {
        use flate2::write::ZlibEncoder;
        use std::io::Write;
        let mut encoder = ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(&[0u8; 1000]).unwrap();
        let block = encoder.finish().unwrap();
        assert_eq!(inflate(&block, 10).unwrap().len(), 11);
        assert_eq!(inflate(&block, 1000).unwrap().len(), 1000);
    }

    #[test]
    fn zero_length_subblocks_are_accepted() {
        let decoded =
            decompress_subblocks(XisfCompression::Zlib, &[1, 2], &[(0, 0), (2, 2)], 2).unwrap();
        assert_eq!(decoded, vec![1, 2]);
    }
}
