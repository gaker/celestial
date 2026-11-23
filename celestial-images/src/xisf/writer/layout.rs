use super::{ImageData, PropertyDataBlock, XisfWriter, HEADER_ALIGNMENT, XISF_SIGNATURE};
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::checksum::Checksum;
use crate::xisf::header::compression::XisfCompression;
use crate::xisf::header::location::DataLocation;
use byteorder::{LittleEndian, WriteBytesExt};
use flate2::write::ZlibEncoder;
use flate2::Compression;
use std::io::{Seek, Write};

pub(super) fn pad_to_alignment(data: &[u8], alignment: usize) -> Vec<u8> {
    let mut padded = data.to_vec();
    let remainder = padded.len() % alignment;
    if remainder != 0 {
        padded.resize(padded.len() + (alignment - remainder), 0);
    }
    padded
}

pub(super) fn align_to(size: usize, alignment: usize) -> usize {
    let remainder = size % alignment;
    if remainder == 0 {
        size
    } else {
        size + (alignment - remainder)
    }
}

// Largest input of the reference codecs: LZ4_MAX_INPUT_SIZE, and the 32-bit
// length of zlib. A larger block would need subblocks (XISF 1.0 section
// 10.6), so the writer stores it uncompressed instead.
const LZ4_MAX_INPUT: usize = 0x7E00_0000;
const ZLIB_MAX_INPUT: u64 = u32::MAX as u64;

// Returns the bytes to store and the codec they are stored with. A block
// that the codec does not make smaller is stored uncompressed. The writer
// has no zstd encoder, so zstd is an error.
pub(super) fn compress_data(
    compression: XisfCompression,
    data: &[u8],
) -> Result<(Vec<u8>, XisfCompression)> {
    if compression == XisfCompression::Zstd {
        return Err(XisfError::UnsupportedCodec(
            compression.as_str().to_string(),
        ));
    }
    Ok(match encode(compression, data) {
        Some(compressed) if compressed.len() < data.len() => (compressed, compression),
        _ => (data.to_vec(), XisfCompression::None),
    })
}

// None for no compression, or when the block is too large for the codec.
fn encode(compression: XisfCompression, data: &[u8]) -> Option<Vec<u8>> {
    match compression {
        XisfCompression::None | XisfCompression::Zstd => None,
        XisfCompression::Lz4 | XisfCompression::Lz4Hc if data.len() > LZ4_MAX_INPUT => None,
        XisfCompression::Zlib if data.len() as u64 > ZLIB_MAX_INPUT => None,
        // XISF stores raw LZ4 blocks; the size is in the compression attribute.
        XisfCompression::Lz4 | XisfCompression::Lz4Hc => Some(lz4_flex::block::compress(data)),
        XisfCompression::Zlib => zlib(data),
    }
}

// The encoder writes to a Vec, which does not fail, so an error here stores
// the block uncompressed.
fn zlib(data: &[u8]) -> Option<Vec<u8>> {
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).ok()?;
    encoder.finish().ok()
}

pub(super) fn build_property_blocks(
    properties: &[crate::xisf::header::property::XisfProperty],
) -> Vec<PropertyDataBlock> {
    let mut blocks = Vec::new();
    for (index, property) in properties.iter().enumerate() {
        if property.value.needs_data_block() {
            if let Some(data) = property.value.to_le_bytes() {
                blocks.push(PropertyDataBlock {
                    property_index: index,
                    data,
                    location: DataLocation::new(0, 0),
                    checksum: None,
                });
            }
        }
    }
    blocks
}

// Places the property blocks, then the image blocks, one after another
// from `first_offset`.
pub(super) fn assign_offsets(
    property_blocks: &mut [PropertyDataBlock],
    images: &mut [ImageData],
    first_offset: u64,
) {
    let mut offset = first_offset;
    for block in property_blocks.iter_mut() {
        let size = block.data.len() as u64;
        block.location = DataLocation::new(offset, size);
        offset += size;
    }
    for image in images.iter_mut() {
        image.info.location = DataLocation::new(offset, image.data.len() as u64);
        offset += image.data.len() as u64;
    }
}

// Attached blocks follow the signature, the 8-byte header length field and
// the padded XML header.
fn first_block_offset(padded_xml_size: usize) -> u64 {
    (XISF_SIGNATURE.len() + 8 + padded_xml_size) as u64
}

// XISF 1.0 section 10.5: the digest covers the block as stored, so a
// compressed block is hashed compressed (section 10.6.2).
pub(super) fn add_sha1_checksums(
    property_blocks: &mut [PropertyDataBlock],
    images: &mut [ImageData],
) {
    for block in property_blocks.iter_mut() {
        block.checksum = Some(Checksum::sha1(&block.data));
    }
    for image in images.iter_mut() {
        image.info.checksum = Some(Checksum::sha1(&image.data));
    }
}

pub(super) fn write_header_length<W: Write>(writer: &mut W, length: usize) -> Result<()> {
    let length = u32::try_from(length).map_err(|_| {
        XisfError::InvalidFormat(format!(
            "The XML header is {length} bytes; its length field holds at most {}",
            u32::MAX
        ))
    })?;
    writer.write_u32::<LittleEndian>(length)?;
    writer.write_u32::<LittleEndian>(0)?;
    Ok(())
}

pub(super) fn write_property_data<W: Write>(
    writer: &mut W,
    property_blocks: &[PropertyDataBlock],
) -> Result<()> {
    for block in property_blocks {
        writer.write_all(&block.data)?;
    }
    Ok(())
}

pub(super) fn write_image_data<W: Write>(writer: &mut W, images: &[ImageData]) -> Result<()> {
    for image in images {
        writer.write_all(&image.data)?;
    }
    Ok(())
}

// `created` is the XISF:CreationTime value.
pub(super) fn generate_final_xml<W: Write + Seek>(
    writer: &mut XisfWriter<W>,
    created: &str,
) -> Result<Vec<u8>> {
    assign_offsets(&mut writer.property_blocks, &mut writer.images, 0);
    let first_pass = super::xml::generate_xml_content(writer, created)?;
    let mut padded_size = align_to(first_pass.len(), HEADER_ALIGNMENT);

    loop {
        let first_offset = first_block_offset(padded_size);
        assign_offsets(
            &mut writer.property_blocks,
            &mut writer.images,
            first_offset,
        );
        let final_xml = super::xml::generate_xml_content(writer, created)?;
        let actual_padded = align_to(final_xml.len(), HEADER_ALIGNMENT);

        if actual_padded <= padded_size {
            return Ok(pad_to_alignment(&final_xml, HEADER_ALIGNMENT));
        }
        padded_size = actual_padded;
    }
}

pub(super) fn validate_data_length(data_len: usize, geometry: &[usize]) -> Result<()> {
    let expected_pixels = geometry
        .iter()
        .try_fold(1usize, |acc, &d| acc.checked_mul(d));
    if expected_pixels != Some(data_len) {
        return Err(XisfError::InvalidFormat(format!(
            "Data length {} does not match geometry {:?}",
            data_len, geometry
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xisf::header::image_info::ImageInfo;
    use crate::xisf::header::property::{XisfProperty, XisfPropertyValue};

    fn make_image(size: u64) -> ImageData {
        ImageData {
            info: ImageInfo {
                geometry: vec![10, 10],
                sample_format: crate::xisf::header::sample_format::SampleFormat::UInt8,
                bounds: (0.0, 255.0),
                color_space: crate::xisf::header::color_space::ColorSpace::Gray,
                pixel_storage: crate::xisf::header::pixel_storage::PixelStorage::Planar,
                location: DataLocation::new(0, 0),
                compression: XisfCompression::None,
                uncompressed_size: None,
                shuffle_item_size: None,
                subblocks: Vec::new(),
                byte_order: crate::core::types::ByteOrder::LittleEndian,
                checksum: None,
                block_defect: None,
            },
            data: vec![0u8; size as usize],
        }
    }

    fn make_block(property_index: usize, size: usize) -> PropertyDataBlock {
        PropertyDataBlock {
            property_index,
            data: vec![0u8; size],
            location: DataLocation::new(0, 0),
            checksum: None,
        }
    }

    #[test]
    fn pad_to_alignment_already_aligned_is_untouched() {
        let data = vec![1u8, 2, 3, 4, 5, 6, 7, 8];
        let padded = pad_to_alignment(&data, 4);
        assert_eq!(padded, data);
    }

    #[test]
    fn pad_to_alignment_needs_padding_appends_zeros() {
        let data = vec![1u8, 2, 3];
        let padded = pad_to_alignment(&data, 8);
        assert_eq!(padded.len(), 8);
        assert_eq!(&padded[..3], &data[..]);
        assert_eq!(&padded[3..], &[0, 0, 0, 0, 0]);
    }

    #[test]
    fn pad_to_alignment_empty_stays_empty() {
        let padded = pad_to_alignment(&[], 16);
        assert!(padded.is_empty());
    }

    #[test]
    fn align_to_zero_stays_zero() {
        assert_eq!(align_to(0, 16), 0);
    }

    #[test]
    fn align_to_already_aligned_is_unchanged() {
        assert_eq!(align_to(16, 16), 16);
        assert_eq!(align_to(32, 16), 32);
        assert_eq!(align_to(48, 16), 48);
    }

    #[test]
    fn align_to_rounds_up_to_next_multiple() {
        assert_eq!(align_to(1, 16), 16);
        assert_eq!(align_to(17, 16), 32);
        assert_eq!(align_to(31, 16), 32);
        assert_eq!(align_to(33, 16), 48);
    }

    #[test]
    fn compress_data_none_returns_passthrough() {
        let data = vec![1u8, 2, 3, 4];
        let (out, used) = compress_data(XisfCompression::None, &data).unwrap();
        assert_eq!(out, data);
        assert_eq!(used, XisfCompression::None);
    }

    #[test]
    fn compress_data_lz4_shrinks_compressible_data() {
        let data = vec![0u8; 4096];
        let (out, used) = compress_data(XisfCompression::Lz4, &data).unwrap();
        assert!(out.len() < data.len());
        assert_eq!(used, XisfCompression::Lz4);
    }

    #[test]
    fn compress_data_lz4_hc_preserves_variant_on_success() {
        let data = vec![42u8; 2048];
        let (out, used) = compress_data(XisfCompression::Lz4Hc, &data).unwrap();
        assert!(out.len() < data.len());
        assert_eq!(used, XisfCompression::Lz4Hc);
    }

    #[test]
    fn compress_data_lz4_falls_back_when_output_not_smaller() {
        let data = vec![7u8, 42, 99];
        let (out, used) = compress_data(XisfCompression::Lz4, &data).unwrap();
        assert_eq!(out, data);
        assert_eq!(used, XisfCompression::None);
    }

    #[test]
    fn compress_data_lz4_writes_a_raw_lz4_block() {
        let data: Vec<u8> = (0..4096u32).map(|i| (i % 7) as u8).collect();
        for codec in [XisfCompression::Lz4, XisfCompression::Lz4Hc] {
            let (out, used) = compress_data(codec, &data).unwrap();
            assert_eq!(used, codec);
            assert_eq!(lz4_flex::block::decompress(&out, data.len()).unwrap(), data);
        }
    }

    #[test]
    fn compress_data_zlib_shrinks_compressible_data() {
        let data = vec![5u8; 4096];
        let (out, used) = compress_data(XisfCompression::Zlib, &data).unwrap();
        assert!(out.len() < data.len());
        assert_eq!(used, XisfCompression::Zlib);
    }

    #[test]
    fn compress_data_zlib_falls_back_when_not_smaller() {
        let data = vec![1u8, 2];
        let (out, used) = compress_data(XisfCompression::Zlib, &data).unwrap();
        assert_eq!(out, data);
        assert_eq!(used, XisfCompression::None);
    }

    #[test]
    fn compress_data_rejects_zstd() {
        let data = vec![0u8; 1024];
        let err = compress_data(XisfCompression::Zstd, &data).unwrap_err();
        assert!(
            matches!(&err, XisfError::UnsupportedCodec(codec) if codec == "zstd"),
            "{err}"
        );
    }

    #[test]
    fn build_property_blocks_skips_scalars() {
        let props = vec![
            XisfProperty::string("A", "x"),
            XisfProperty::float64("B", 1.5),
            XisfProperty::int32("C", 3),
            XisfProperty::boolean("D", true),
        ];
        assert!(build_property_blocks(&props).is_empty());
    }

    #[test]
    fn build_property_blocks_captures_vectors_and_matrices() {
        let props = vec![
            XisfProperty::string("skip", "x"),
            XisfProperty::f64_vector("vec", vec![1.0, 2.0, 3.0]),
            XisfProperty::float64("scalar", 1.0),
            XisfProperty::f64_matrix("mat", 2, 2, vec![1.0, 0.0, 0.0, 1.0]),
        ];
        let blocks = build_property_blocks(&props);
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].property_index, 1);
        assert_eq!(blocks[0].data.len(), 24);
        assert_eq!(blocks[1].property_index, 3);
        assert_eq!(blocks[1].data.len(), 32);
    }

    #[test]
    fn build_property_blocks_skips_value_with_no_bytes() {
        let props: Vec<XisfProperty> = vec![XisfProperty::new(
            "id",
            XisfPropertyValue::F64Vector(vec![]),
        )];
        let blocks = build_property_blocks(&props);
        assert_eq!(blocks.len(), 1);
        assert!(blocks[0].data.is_empty());
    }

    #[test]
    fn assign_offsets_chains_across_blocks_and_images() {
        let mut blocks = vec![make_block(0, 16), make_block(1, 32)];
        let mut images = vec![make_image(100), make_image(200)];
        assign_offsets(&mut blocks, &mut images, 0);
        assert_eq!(blocks[0].location.offset, 0);
        assert_eq!(blocks[0].location.size, 16);
        assert_eq!(blocks[1].location.offset, 16);
        assert_eq!(blocks[1].location.size, 32);
        assert_eq!(images[0].info.location.offset, 48);
        assert_eq!(images[0].info.location.size, 100);
        assert_eq!(images[1].info.location.offset, 148);
        assert_eq!(images[1].info.location.size, 200);
    }

    #[test]
    fn first_block_offset_starts_after_signature_and_xml() {
        let mut blocks = vec![make_block(0, 10)];
        let mut images = vec![make_image(20)];
        let padded_xml_size = 64;
        let first_offset = first_block_offset(padded_xml_size);
        assign_offsets(&mut blocks, &mut images, first_offset);
        let header_size = XISF_SIGNATURE.len() + 8 + padded_xml_size;
        assert_eq!(blocks[0].location.offset, header_size as u64);
        assert_eq!(images[0].info.location.offset, (header_size + 10) as u64);
    }

    #[test]
    fn assign_offsets_with_no_blocks_places_image_at_header_end() {
        let mut blocks: Vec<PropertyDataBlock> = Vec::new();
        let mut images = vec![make_image(50)];
        assign_offsets(&mut blocks, &mut images, first_block_offset(32));
        let header_size = XISF_SIGNATURE.len() + 8 + 32;
        assert_eq!(images[0].info.location.offset, header_size as u64);
    }

    #[test]
    fn write_header_length_emits_eight_bytes() {
        let mut buf = Vec::new();
        write_header_length(&mut buf, 0x1234ABCD).unwrap();
        assert_eq!(buf.len(), 8);
        assert_eq!(&buf[..4], &0x1234ABCDu32.to_le_bytes());
        assert_eq!(&buf[4..], &[0u8; 4]);
    }

    #[test]
    fn write_header_length_rejects_a_length_past_u32() {
        let mut buf = Vec::new();
        write_header_length(&mut buf, u32::MAX as usize).unwrap();
        assert_eq!(&buf[..4], &[0xFF; 4]);

        let error = write_header_length(&mut buf, 1 << 32)
            .err()
            .map(|e| e.to_string());
        assert_eq!(
            error.as_deref(),
            Some(
                "Invalid XISF format: The XML header is 4294967296 bytes; \
                 its length field holds at most 4294967295"
            )
        );
    }

    #[test]
    fn write_property_data_concatenates_block_bytes() {
        let blocks = vec![
            PropertyDataBlock {
                property_index: 0,
                data: vec![1, 2, 3],
                location: DataLocation::new(0, 3),
                checksum: None,
            },
            PropertyDataBlock {
                property_index: 1,
                data: vec![4, 5],
                location: DataLocation::new(3, 2),
                checksum: None,
            },
        ];
        let mut buf = Vec::new();
        write_property_data(&mut buf, &blocks).unwrap();
        assert_eq!(buf, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn write_image_data_concatenates_image_bytes() {
        let images = vec![
            ImageData {
                info: make_image(0).info,
                data: vec![10, 20, 30],
            },
            ImageData {
                info: make_image(0).info,
                data: vec![40, 50],
            },
        ];
        let mut buf = Vec::new();
        write_image_data(&mut buf, &images).unwrap();
        assert_eq!(buf, vec![10, 20, 30, 40, 50]);
    }

    #[test]
    fn validate_data_length_accepts_matching_geometry() {
        assert!(validate_data_length(100, &[10, 10]).is_ok());
        assert!(validate_data_length(300, &[10, 10, 3]).is_ok());
        assert!(validate_data_length(0, &[0, 5]).is_ok());
    }

    #[test]
    fn validate_data_length_rejects_overflowing_geometry() {
        assert!(validate_data_length(0, &[usize::MAX, 2]).is_err());
        assert!(validate_data_length(0, &[usize::MAX / 2, 2, 3]).is_err());
    }

    #[test]
    fn validate_data_length_rejects_mismatch_with_details() {
        let err = validate_data_length(99, &[10, 10]).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("99"));
        assert!(msg.contains("[10, 10]"));

        let err = validate_data_length(50, &[10, 10, 3]).unwrap_err();
        assert!(matches!(err, XisfError::InvalidFormat(_)));
    }
}
