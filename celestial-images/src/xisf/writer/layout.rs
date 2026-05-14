use super::{ImageData, PropertyDataBlock, XisfWriter, HEADER_ALIGNMENT, XISF_SIGNATURE};
use crate::xisf::header::{DataLocation, XisfCompression};
use crate::xisf::{Result, XisfError};
use byteorder::{LittleEndian, WriteBytesExt};
use std::io::{Seek, Write};

pub(super) fn build_geometry(width: usize, height: usize, channels: usize) -> Vec<usize> {
    if channels > 1 {
        vec![width, height, channels]
    } else {
        vec![width, height]
    }
}

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

pub(super) fn compress_data(
    compression: XisfCompression,
    data: &[u8],
) -> (Vec<u8>, XisfCompression) {
    match compression {
        XisfCompression::None => (data.to_vec(), XisfCompression::None),
        XisfCompression::Lz4 | XisfCompression::Lz4Hc => {
            let compressed = lz4_flex::compress_prepend_size(data);
            if compressed.len() < data.len() {
                (compressed, compression)
            } else {
                (data.to_vec(), XisfCompression::None)
            }
        }
        XisfCompression::Zlib => {
            use flate2::write::ZlibEncoder;
            use flate2::Compression;
            let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
            if encoder.write_all(data).is_ok() {
                if let Ok(compressed) = encoder.finish() {
                    if compressed.len() < data.len() {
                        return (compressed, XisfCompression::Zlib);
                    }
                }
            }
            (data.to_vec(), XisfCompression::None)
        }
        XisfCompression::Zstd => (data.to_vec(), XisfCompression::None),
    }
}

pub(super) fn build_property_blocks(
    properties: &[crate::xisf::header::XisfProperty],
) -> Vec<PropertyDataBlock> {
    let mut blocks = Vec::new();
    for (index, property) in properties.iter().enumerate() {
        if property.value.needs_data_block() {
            if let Some(data) = property.value.to_le_bytes() {
                blocks.push(PropertyDataBlock {
                    property_index: index,
                    data,
                    location: DataLocation::new(0, 0),
                });
            }
        }
    }
    blocks
}

pub(super) fn set_placeholder_offsets(
    property_blocks: &mut [PropertyDataBlock],
    images: &mut [ImageData],
) {
    let mut offset = 0u64;
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

pub(super) fn calculate_final_offsets(
    property_blocks: &mut [PropertyDataBlock],
    images: &mut [ImageData],
    padded_xml_size: usize,
) {
    let header_size = XISF_SIGNATURE.len() + 8 + padded_xml_size;
    let mut offset = header_size as u64;
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

pub(super) fn write_header_length<W: Write>(writer: &mut W, length: u32) -> Result<()> {
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

pub(super) fn generate_final_xml<W: Write + Seek>(writer: &mut XisfWriter<W>) -> Result<Vec<u8>> {
    set_placeholder_offsets(&mut writer.property_blocks, &mut writer.images);
    let first_pass = super::xml::generate_xml_content(writer)?;
    let mut padded_size = align_to(first_pass.len(), HEADER_ALIGNMENT);

    loop {
        calculate_final_offsets(
            &mut writer.property_blocks,
            &mut writer.images,
            padded_size,
        );
        let final_xml = super::xml::generate_xml_content(writer)?;
        let actual_padded = align_to(final_xml.len(), HEADER_ALIGNMENT);

        if actual_padded <= padded_size {
            return Ok(pad_to_alignment(&final_xml, HEADER_ALIGNMENT));
        }
        padded_size = actual_padded;
    }
}

pub(super) fn validate_data_length(data_len: usize, geometry: &[usize]) -> Result<()> {
    let expected_pixels: usize = geometry.iter().product();
    if data_len != expected_pixels {
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
    use crate::xisf::header::{ImageInfo, XisfProperty, XisfPropertyValue};

    fn make_image(size: u64) -> ImageData {
        ImageData {
            info: ImageInfo {
                geometry: vec![10, 10],
                sample_format: crate::xisf::header::SampleFormat::UInt8,
                bounds: (0.0, 255.0),
                color_space: crate::xisf::header::ColorSpace::Gray,
                pixel_storage: crate::xisf::header::PixelStorage::Planar,
                location: DataLocation::new(0, 0),
                compression: XisfCompression::None,
                uncompressed_size: None,
            },
            data: vec![0u8; size as usize],
        }
    }

    fn make_block(property_index: usize, size: usize) -> PropertyDataBlock {
        PropertyDataBlock {
            property_index,
            data: vec![0u8; size],
            location: DataLocation::new(0, 0),
        }
    }

    #[test]
    fn build_geometry_single_channel_is_two_dim() {
        assert_eq!(build_geometry(1920, 1080, 1), vec![1920, 1080]);
        assert_eq!(build_geometry(100, 50, 0), vec![100, 50]);
    }

    #[test]
    fn build_geometry_multi_channel_is_three_dim() {
        assert_eq!(build_geometry(1920, 1080, 3), vec![1920, 1080, 3]);
        assert_eq!(build_geometry(64, 64, 4), vec![64, 64, 4]);
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
        let (out, used) = compress_data(XisfCompression::None, &data);
        assert_eq!(out, data);
        assert_eq!(used, XisfCompression::None);
    }

    #[test]
    fn compress_data_lz4_shrinks_compressible_data() {
        let data = vec![0u8; 4096];
        let (out, used) = compress_data(XisfCompression::Lz4, &data);
        assert!(out.len() < data.len());
        assert_eq!(used, XisfCompression::Lz4);
    }

    #[test]
    fn compress_data_lz4_hc_preserves_variant_on_success() {
        let data = vec![42u8; 2048];
        let (out, used) = compress_data(XisfCompression::Lz4Hc, &data);
        assert!(out.len() < data.len());
        assert_eq!(used, XisfCompression::Lz4Hc);
    }

    #[test]
    fn compress_data_lz4_falls_back_when_output_not_smaller() {
        let data = vec![7u8, 42, 99];
        let (out, used) = compress_data(XisfCompression::Lz4, &data);
        assert_eq!(out, data);
        assert_eq!(used, XisfCompression::None);
    }

    #[test]
    fn compress_data_zlib_shrinks_compressible_data() {
        let data = vec![5u8; 4096];
        let (out, used) = compress_data(XisfCompression::Zlib, &data);
        assert!(out.len() < data.len());
        assert_eq!(used, XisfCompression::Zlib);
    }

    #[test]
    fn compress_data_zlib_falls_back_when_not_smaller() {
        let data = vec![1u8, 2];
        let (out, used) = compress_data(XisfCompression::Zlib, &data);
        assert_eq!(out, data);
        assert_eq!(used, XisfCompression::None);
    }

    #[test]
    fn compress_data_zstd_always_falls_back() {
        let data = vec![0u8; 1024];
        let (out, used) = compress_data(XisfCompression::Zstd, &data);
        assert_eq!(out, data);
        assert_eq!(used, XisfCompression::None);
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
    fn set_placeholder_offsets_chains_across_blocks_and_images() {
        let mut blocks = vec![make_block(0, 16), make_block(1, 32)];
        let mut images = vec![make_image(100), make_image(200)];
        set_placeholder_offsets(&mut blocks, &mut images);
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
    fn calculate_final_offsets_starts_after_signature_and_xml() {
        let mut blocks = vec![make_block(0, 10)];
        let mut images = vec![make_image(20)];
        let padded_xml_size = 64;
        calculate_final_offsets(&mut blocks, &mut images, padded_xml_size);
        let header_size = XISF_SIGNATURE.len() + 8 + padded_xml_size;
        assert_eq!(blocks[0].location.offset, header_size as u64);
        assert_eq!(images[0].info.location.offset, (header_size + 10) as u64);
    }

    #[test]
    fn calculate_final_offsets_with_no_blocks_places_image_at_header_end() {
        let mut blocks: Vec<PropertyDataBlock> = Vec::new();
        let mut images = vec![make_image(50)];
        calculate_final_offsets(&mut blocks, &mut images, 32);
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
    fn write_property_data_concatenates_block_bytes() {
        let blocks = vec![
            PropertyDataBlock {
                property_index: 0,
                data: vec![1, 2, 3],
                location: DataLocation::new(0, 3),
            },
            PropertyDataBlock {
                property_index: 1,
                data: vec![4, 5],
                location: DataLocation::new(3, 2),
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
    fn validate_data_length_rejects_mismatch_with_details() {
        let err = validate_data_length(99, &[10, 10]).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("99"));
        assert!(msg.contains("[10, 10]"));

        let err = validate_data_length(50, &[10, 10, 3]).unwrap_err();
        assert!(matches!(err, XisfError::InvalidFormat(_)));
    }
}
