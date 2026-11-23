pub mod data_type;
mod layout;
mod metadata;
pub(crate) mod wcs;
mod xml;

use data_type::XisfDataType;

use crate::core::shape::plane_dimensions;
use crate::fits::header::keywords::{set_keyword, Keyword};
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::checksum::Checksum;
use crate::xisf::header::color_space::ColorSpace;
use crate::xisf::header::compression::XisfCompression;
use crate::xisf::header::image_info::ImageInfo;
use crate::xisf::header::location::DataLocation;
use crate::xisf::header::pixel_storage::PixelStorage;
use crate::xisf::header::property::XisfProperty;
use crate::xisf::header::sample_format::SampleFormat;
use std::fs::File;
use std::io::{BufWriter, Cursor, Seek, Write};
use std::path::Path;
use std::time::SystemTime;

const XISF_SIGNATURE: &[u8] = b"XISF0100";
const HEADER_ALIGNMENT: usize = 16;

pub struct XisfWriter<W> {
    pub(super) writer: W,
    pub(super) images: Vec<ImageData>,
    // The FITSKeyword and Property children of image 0.
    pub(super) keywords: Vec<Keyword>,
    pub(super) properties: Vec<XisfProperty>,
    pub(super) compression: XisfCompression,
    pub(super) property_blocks: Vec<PropertyDataBlock>,
    sha1_checksums: bool,
    // The XISF:CreationTime value. None writes the system time at `write`.
    // Tests fix it so the output is the same on every run.
    pub(crate) creation_time: Option<SystemTime>,
}

pub(super) struct ImageData {
    pub(super) info: ImageInfo,
    pub(super) data: Vec<u8>,
}

pub(super) struct PropertyDataBlock {
    pub(super) property_index: usize,
    pub(super) data: Vec<u8>,
    pub(super) location: DataLocation,
    pub(super) checksum: Option<Checksum>,
}

impl XisfWriter<BufWriter<File>> {
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path)?;
        Ok(Self::new(BufWriter::new(file)))
    }
}

impl XisfWriter<Cursor<Vec<u8>>> {
    pub fn write_to_vec(mut self) -> Result<Vec<u8>> {
        self.write_internal()?;
        Ok(self.writer.into_inner())
    }
}

impl<W: Write + Seek> XisfWriter<W> {
    pub fn new(writer: W) -> Self {
        Self {
            writer,
            images: Vec::new(),
            keywords: Vec::new(),
            properties: Vec::new(),
            compression: XisfCompression::None,
            property_blocks: Vec::new(),
            sha1_checksums: false,
            creation_time: None,
        }
    }

    pub fn compression(mut self, compression: XisfCompression) -> Self {
        self.compression = compression;
        self
    }

    // Gives every attached data block a SHA-1 checksum attribute (XISF 1.0
    // section 10.5), computed over the stored bytes.
    pub fn sha1_checksums(mut self, enabled: bool) -> Self {
        self.sha1_checksums = enabled;
        self
    }

    pub fn add_image<T: XisfDataType>(
        &mut self,
        data: &[T],
        width: usize,
        height: usize,
        channels: usize,
    ) -> Result<()> {
        let bounds = T::calculate_bounds(data);
        self.add_image_with_bounds(data, width, height, channels, bounds)
    }

    pub fn add_image_with_bounds<T: XisfDataType>(
        &mut self,
        data: &[T],
        width: usize,
        height: usize,
        channels: usize,
        bounds: (f64, f64),
    ) -> Result<()> {
        let geometry = plane_dimensions(width, height, channels);
        check_image(&geometry, bounds, data.len())?;
        let mut info = image_info(geometry, T::SAMPLE_FORMAT, bounds);
        let data = stored_block(self.compression, &T::to_le_bytes(data), &mut info)?;
        self.images.push(ImageData { info, data });
        Ok(())
    }

    // Keywords and properties added here belong to image 0. Replaces every
    // keyword of the same name with one; `add_keyword` appends.
    pub fn set_keyword(&mut self, name: &str, value: &str) {
        set_keyword(&mut self.keywords, Keyword::string(name, value));
    }

    pub fn add_keyword(&mut self, keyword: Keyword) {
        self.keywords.push(keyword);
    }

    pub fn add_property(&mut self, property: XisfProperty) {
        self.properties.push(property);
    }

    pub fn add_properties(&mut self, properties: impl IntoIterator<Item = XisfProperty>) {
        self.properties.extend(properties);
    }

    pub fn write(mut self) -> Result<()> {
        self.write_internal()?;
        self.writer.flush()?;
        Ok(())
    }

    fn write_internal(&mut self) -> Result<()> {
        // The keywords and properties are children of image 0, and their
        // data blocks would be written with no element to reference them.
        let has_metadata = !self.keywords.is_empty() || !self.properties.is_empty();
        if has_metadata && self.images.is_empty() {
            return Err(XisfError::InvalidFormat(
                "Keywords and properties are written with image 0, and the writer has no image"
                    .to_string(),
            ));
        }
        let created = self.creation_time.unwrap_or_else(SystemTime::now);
        let created = metadata::time_point(created)?;
        self.property_blocks = layout::build_property_blocks(&self.properties);
        if self.sha1_checksums {
            layout::add_sha1_checksums(&mut self.property_blocks, &mut self.images);
        }
        self.writer.write_all(XISF_SIGNATURE)?;

        let padded_xml = layout::generate_final_xml(self, &created)?;
        layout::write_header_length(&mut self.writer, padded_xml.len())?;
        self.writer.write_all(&padded_xml)?;
        layout::write_property_data(&mut self.writer, &self.property_blocks)?;
        layout::write_image_data(&mut self.writer, &self.images)?;
        Ok(())
    }
}

impl<W> XisfWriter<W> {
    // The keywords and properties written as children of image `index`.
    pub(super) fn image_metadata(&self, index: usize) -> (&[Keyword], &[XisfProperty]) {
        if index == 0 {
            (&self.keywords, &self.properties)
        } else {
            (&[], &[])
        }
    }
}

// XISF 1.0 geometry dimensions are at least 1; readers refuse 0.
fn check_image(geometry: &[usize], bounds: (f64, f64), data_len: usize) -> Result<()> {
    if geometry.contains(&0) {
        return Err(XisfError::InvalidFormat(format!(
            "Image geometry {geometry:?} has a zero dimension"
        )));
    }
    check_bounds(bounds)?;
    layout::validate_data_length(data_len, geometry)
}

// XISF 1.0 section 11.5.1 gives the representable range as lower:upper; a
// range that is not two finite numbers in order cannot describe the samples.
fn check_bounds((lower, upper): (f64, f64)) -> Result<()> {
    if lower.is_finite() && upper.is_finite() && lower <= upper {
        return Ok(());
    }
    Err(XisfError::InvalidFormat(format!(
        "Image bounds ({lower}, {upper}) are not two finite numbers with lower <= upper"
    )))
}

// Returns the bytes to store for `raw` and sets the compression and size
// attributes of `info` to match them.
fn stored_block(compression: XisfCompression, raw: &[u8], info: &mut ImageInfo) -> Result<Vec<u8>> {
    let (stored, used) = layout::compress_data(compression, raw)?;
    info.location = DataLocation::new(0, stored.len() as u64);
    info.compression = used;
    info.uncompressed_size = used.is_compressed().then_some(raw.len() as u64);
    Ok(stored)
}

fn image_info(geometry: Vec<usize>, sample_format: SampleFormat, bounds: (f64, f64)) -> ImageInfo {
    let color_space = if geometry.get(2) == Some(&3) {
        ColorSpace::Rgb
    } else {
        ColorSpace::Gray
    };
    ImageInfo {
        geometry,
        sample_format,
        bounds,
        color_space,
        pixel_storage: PixelStorage::Planar,
        location: DataLocation::new(0, 0),
        compression: XisfCompression::None,
        uncompressed_size: None,
        shuffle_item_size: None,
        subblocks: Vec::new(),
        byte_order: crate::core::types::ByteOrder::LittleEndian,
        checksum: None,
        block_defect: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::{Keyword, KeywordValue};
    use crate::xisf::errors::XisfError;
    use crate::xisf::header::property::XisfProperty;
    use crate::xisf::header::sample_format::SampleFormat;
    use crate::xisf::reader::XisfFile;
    use celestial_time::scales::utc::UTC;
    use std::io::Cursor;
    use std::time::{Duration, UNIX_EPOCH};

    fn roundtrip<T: XisfDataType>(
        data: &[T],
        width: usize,
        height: usize,
        channels: usize,
        compression: XisfCompression,
    ) -> XisfFile<Cursor<Vec<u8>>> {
        let mut writer = XisfWriter::new(Cursor::new(Vec::new())).compression(compression);
        writer.add_image(data, width, height, channels).unwrap();
        let bytes = writer.write_to_vec().unwrap();
        XisfFile::new(Cursor::new(bytes)).unwrap()
    }

    #[test]
    fn new_writer_starts_empty() {
        let w = XisfWriter::new(Cursor::new(Vec::new()));
        assert!(w.images.is_empty());
        assert!(w.keywords.is_empty());
        assert!(w.properties.is_empty());
        assert_eq!(w.compression, XisfCompression::None);
        assert!(w.property_blocks.is_empty());
    }

    #[test]
    fn compression_builder_sets_field() {
        let w = XisfWriter::new(Cursor::new(Vec::new())).compression(XisfCompression::Lz4);
        assert_eq!(w.compression, XisfCompression::Lz4);
    }

    #[test]
    fn add_image_length_mismatch_is_rejected() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        let err = w.add_image::<u8>(&[0, 1, 2], 2, 2, 1).unwrap_err();
        assert!(matches!(err, XisfError::InvalidFormat(_)));
        assert!(err.to_string().contains("Data length"));
    }

    #[test]
    fn roundtrip_u8_gray() {
        let data: Vec<u8> = (0..100).collect();
        let mut reader = roundtrip(&data, 10, 10, 1, XisfCompression::None);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.geometry, vec![10, 10]);
        assert!(matches!(info.sample_format, SampleFormat::UInt8));

        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<u8 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn roundtrip_u16_gray() {
        let data: Vec<u16> = (0u16..64).collect();
        let mut reader = roundtrip(&data, 8, 8, 1, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().sample_format,
            SampleFormat::UInt16
        ));
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<u16 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn roundtrip_u32_gray() {
        let data: Vec<u32> = (0u32..36).collect();
        let mut reader = roundtrip(&data, 6, 6, 1, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().sample_format,
            SampleFormat::UInt32
        ));
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<u32 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn roundtrip_f32_rgb_preserves_bounds() {
        let data: Vec<f32> = (0..300).map(|i| i as f32 / 300.0).collect();
        let mut reader = roundtrip(&data, 10, 10, 3, XisfCompression::None);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.geometry, vec![10, 10, 3]);
        assert!(matches!(info.sample_format, SampleFormat::Float32));
        assert!(matches!(info.color_space, ColorSpace::Rgb));
        let (min, max) = info.bounds;
        assert!(min >= 0.0 && max <= 1.0);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<f32 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn roundtrip_f64_gray() {
        let data: Vec<f64> = (0..16).map(|i| i as f64 * 0.25).collect();
        let mut reader = roundtrip(&data, 4, 4, 1, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().sample_format,
            SampleFormat::Float64
        ));
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(<f64 as XisfDataType>::from_le_bytes(&bytes), data);
    }

    #[test]
    fn single_channel_gets_gray_colorspace() {
        let data = vec![128u8; 25];
        let reader = roundtrip(&data, 5, 5, 1, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().color_space,
            ColorSpace::Gray
        ));
    }

    #[test]
    fn three_channel_gets_rgb_colorspace() {
        let data = vec![128u8; 75];
        let reader = roundtrip(&data, 5, 5, 3, XisfCompression::None);
        assert!(matches!(
            reader.image_info(0).unwrap().color_space,
            ColorSpace::Rgb
        ));
    }

    #[test]
    fn lz4_compressed_roundtrip_preserves_data() {
        let data = vec![42u8; 4096];
        let mut reader = roundtrip(&data, 64, 64, 1, XisfCompression::Lz4);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::Lz4);
        assert_eq!(info.uncompressed_size, Some(4096));
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, data);
    }

    #[test]
    fn lz4_hc_compressed_roundtrip_preserves_data() {
        let data = vec![7u8; 2048];
        let mut reader = roundtrip(&data, 32, 64, 1, XisfCompression::Lz4Hc);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::Lz4Hc);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, data);
    }

    #[test]
    fn zlib_compressed_roundtrip_preserves_data() {
        let data = vec![3u8; 4096];
        let mut reader = roundtrip(&data, 64, 64, 1, XisfCompression::Zlib);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::Zlib);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, data);
    }

    #[test]
    fn keywords_or_properties_without_an_image_are_rejected() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_property(XisfProperty::f64_vector("V", vec![1.0]));
        let err = w.write_to_vec().unwrap_err();
        assert!(matches!(err, XisfError::InvalidFormat(_)), "{err}");

        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_keyword(Keyword::string("OBJECT", "M31"));
        let err = w.write_to_vec().unwrap_err();
        assert!(matches!(err, XisfError::InvalidFormat(_)), "{err}");
    }

    #[test]
    fn a_writer_without_images_or_metadata_writes_an_empty_unit() {
        let bytes = XisfWriter::new(Cursor::new(Vec::new()))
            .write_to_vec()
            .unwrap();
        assert_eq!(XisfFile::new(Cursor::new(bytes)).unwrap().num_images(), 0);
    }

    #[test]
    fn zstd_is_rejected() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new())).compression(XisfCompression::Zstd);
        let err = w.add_image(&[5u8; 512], 16, 32, 1).unwrap_err();
        assert!(
            matches!(&err, XisfError::UnsupportedCodec(codec) if codec == "zstd"),
            "{err}"
        );
    }

    #[test]
    fn incompressible_lz4_falls_back_to_uncompressed() {
        let data: Vec<u8> = (0..255).collect();
        let mut reader = roundtrip(&data, 255, 1, 1, XisfCompression::Lz4);
        let info = reader.image_info(0).unwrap();
        assert_eq!(info.compression, XisfCompression::None);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, data);
    }

    #[test]
    fn add_image_with_bounds_uses_explicit_bounds() {
        let data: Vec<f32> = vec![0.5; 4];
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image_with_bounds(&data, 2, 2, 1, (-10.0, 10.0))
            .unwrap();
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.image_info(0).unwrap().bounds, (-10.0, 10.0));
    }

    #[test]
    fn set_keyword_appears_in_roundtrip() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.set_keyword("TELESCOP", "Hubble");
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        let kw = reader.get_keyword("TELESCOP").unwrap();
        assert_eq!(kw.value, Some(KeywordValue::String("Hubble".to_string())));
    }

    #[test]
    fn set_keyword_replaces_previous_value() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.set_keyword("OBJECT", "A");
        w.set_keyword("OBJECT", "B");
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        let objects: Vec<_> = reader
            .keywords()
            .iter()
            .filter(|k| k.name == "OBJECT")
            .map(|k| k.value.clone())
            .collect();
        assert_eq!(objects, [Some(KeywordValue::String("B".into()))]);
    }

    #[test]
    fn add_keyword_appears_in_roundtrip() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.add_keyword(Keyword::integer("EXPTIME", 30).with_comment("seconds"));
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        let kw = reader.get_keyword("EXPTIME").unwrap();
        assert_eq!(kw.value, Some(KeywordValue::Integer(30)));
    }

    #[test]
    fn keywords_and_properties_read_back_on_image_0_only() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.add_image::<u8>(&[1u8; 4], 2, 2, 1).unwrap();
        w.add_keyword(Keyword::string("OBJECT", "M31"));
        w.add_property(XisfProperty::f64_vector("V", vec![1.0, 2.0]));
        let bytes = w.write_to_vec().unwrap();
        let mut reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        let object = [Keyword::string("OBJECT", "M31")];
        assert_eq!(reader.keywords(), object);
        assert_eq!(reader.image_keywords(0), Some(&object[..]));
        assert_eq!(reader.image_keywords(1), Some(&[][..]));
        let vector = XisfProperty::f64_vector("V", vec![1.0, 2.0]);
        assert_eq!(reader.properties(0).unwrap(), [vector]);
        assert_eq!(reader.properties(1).unwrap(), []);
    }

    // XISF 1.0 section 10.5: the digest covers the stored bytes, which are
    // the compressed bytes of a compressed block (section 10.6.2).
    #[test]
    fn sha1_checksums_are_written_and_verified_by_the_reader() {
        for compression in [XisfCompression::None, XisfCompression::Zlib] {
            let data = vec![9u8; 64];
            let mut w = XisfWriter::new(Cursor::new(Vec::new()))
                .compression(compression)
                .sha1_checksums(true);
            w.add_image(&data, 8, 8, 1).unwrap();
            w.add_property(XisfProperty::f64_vector("V", vec![1.0, 2.0]));
            let bytes = w.write_to_vec().unwrap();

            let mut reader = XisfFile::new(Cursor::new(bytes.clone())).unwrap();
            let info = reader.image_info(0).unwrap().clone();
            assert_eq!(info.compression, compression);
            let start = info.location.offset as usize;
            let stored = &bytes[start..start + info.location.size as usize];
            assert_eq!(info.checksum, Some(Checksum::sha1(stored)));
            assert_eq!(reader.read_image_data_raw(0).unwrap(), data);
            let vector = XisfProperty::f64_vector("V", vec![1.0, 2.0]);
            assert_eq!(reader.properties(0).unwrap(), [vector]);

            let mut corrupted = bytes;
            corrupted[start] ^= 1;
            let mut reader = XisfFile::new(Cursor::new(corrupted)).unwrap();
            assert!(matches!(
                reader.read_image_data_raw(0),
                Err(XisfError::ChecksumMismatch { .. })
            ));
        }
    }

    #[test]
    fn sha1_checksums_cover_property_blocks() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new())).sha1_checksums(true);
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.add_property(XisfProperty::f64_vector("V", vec![1.0, 2.0]));
        let mut bytes = w.write_to_vec().unwrap();
        let header_end = 16 + u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]);
        // The property block is written first, right after the header.
        bytes[header_end as usize] ^= 1;
        let mut reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.properties(0).unwrap(), []);
    }

    #[test]
    fn checksums_are_off_by_default() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new())).compression(XisfCompression::Zlib);
        w.add_image(&[9u8; 64], 8, 8, 1).unwrap();
        w.add_property(XisfProperty::f64_vector("V", vec![1.0, 2.0]));
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes.clone())).unwrap();
        assert_eq!(reader.image_info(0).unwrap().checksum, None);
        let text = String::from_utf8_lossy(&bytes);
        assert!(!text.contains("checksum="));
        assert!(!text.contains("ChecksumAlgorithms"));
    }

    fn unix_now() -> i64 {
        let since = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap();
        i64::try_from(since.as_secs()).unwrap()
    }

    // Writes a file with one image, property A and keyword OBJECT, checks
    // that the Metadata element is the last child of the root element and
    // holds `extra` after the two properties XISF 1.0 section 11.4.1 makes
    // mandatory, and that the metadata does not change what the reader gives
    // for the image.
    fn check_metadata(w: XisfWriter<Cursor<Vec<u8>>>, extra: &str) {
        let mut w = w;
        // 64 equal bytes, which every codec makes smaller.
        w.add_image(&[9u8; 64], 8, 8, 1).unwrap();
        w.add_property(XisfProperty::int32("A", 1));
        w.add_keyword(Keyword::string("OBJECT", "M31"));
        let before = unix_now();
        let bytes = w.write_to_vec().unwrap();
        let after = unix_now();
        let text = String::from_utf8_lossy(&bytes);
        let time_attribute = "type=\"TimePoint\" value=\"";
        let time_start = text.find(time_attribute).unwrap() + time_attribute.len();
        let time = &text[time_start..time_start + 20];
        let metadata = format!(
            "  <Metadata>\n    \
             <Property id=\"XISF:CreationTime\" type=\"TimePoint\" value=\"{time}\"/>\n    \
             <Property id=\"XISF:CreatorApplication\" type=\"String\">celestial-images {}</Property>\n\
             {extra}  </Metadata>\n</xisf>",
            env!("CARGO_PKG_VERSION")
        );
        assert!(text.contains(&metadata), "{text}");
        assert_eq!(text.matches("<Metadata>").count(), 1);
        assert!(time.ends_with('Z'), "{time}");
        let (seconds, nanos) = time.parse::<UTC>().unwrap().to_unix().unwrap();
        assert_eq!(nanos, 0);
        assert!((before..=after).contains(&seconds), "{time}");

        let mut reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.properties(0).unwrap(), [XisfProperty::int32("A", 1)]);
        assert_eq!(reader.keywords(), [Keyword::string("OBJECT", "M31")]);
    }

    // XISF 1.0 section 11.4.2: XISF:CompressionCodecs lists the codecs of
    // the compressed data blocks.
    #[test]
    fn compression_codecs_are_listed_in_the_metadata() {
        check_metadata(
            XisfWriter::new(Cursor::new(Vec::new())).compression(XisfCompression::Zlib),
            "    <Property id=\"XISF:CompressionCodecs\" type=\"String\">zlib</Property>\n",
        );
        check_metadata(
            XisfWriter::new(Cursor::new(Vec::new()))
                .compression(XisfCompression::Lz4Hc)
                .sha1_checksums(true),
            concat!(
                "    <Property id=\"XISF:ChecksumAlgorithms\" type=\"String\">sha-1</Property>\n",
                "    <Property id=\"XISF:CompressionCodecs\" type=\"String\">lz4hc</Property>\n",
            ),
        );
    }

    // A block the codec does not make smaller is stored uncompressed, and
    // then no codec is listed.
    #[test]
    fn compression_codecs_are_left_out_without_compressed_blocks() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new())).compression(XisfCompression::Zlib);
        w.add_image(&[1u8, 2, 3, 4], 2, 2, 1).unwrap();
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes.clone())).unwrap();
        assert_eq!(
            reader.image_info(0).unwrap().compression,
            XisfCompression::None
        );
        assert!(!String::from_utf8_lossy(&bytes).contains("CompressionCodecs"));
    }

    #[test]
    fn a_fixed_creation_time_is_written() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.creation_time = Some(UNIX_EPOCH + Duration::from_secs(1_418_128_695));
        w.add_image(&[9u8; 4], 2, 2, 1).unwrap();
        let text = String::from_utf8_lossy(&w.write_to_vec().unwrap()).into_owned();
        let property =
            r#"<Property id="XISF:CreationTime" type="TimePoint" value="2014-12-09T12:38:15Z"/>"#;
        assert!(text.contains(property), "{text}");
    }

    #[test]
    fn a_creation_time_before_1970_is_an_error() {
        let mut out = Vec::new();
        let mut w = XisfWriter::new(Cursor::new(&mut out));
        w.creation_time = Some(UNIX_EPOCH - Duration::from_secs(1));
        w.add_image(&[9u8; 4], 2, 2, 1).unwrap();
        let err = w.write().unwrap_err();
        assert!(matches!(err, XisfError::CreationTime(_)), "{err}");
        assert_eq!(out, []);
    }

    // XISF 1.0 section 11.4: every header has one Metadata element.
    #[test]
    fn every_file_has_the_mandatory_metadata() {
        check_metadata(XisfWriter::new(Cursor::new(Vec::new())), "");
    }

    // XISF 1.0 section 11.4.2: the XISF:ChecksumAlgorithms property lists
    // the algorithms of the data block checksums.
    #[test]
    fn sha1_checksums_are_listed_in_the_metadata() {
        check_metadata(
            XisfWriter::new(Cursor::new(Vec::new())).sha1_checksums(true),
            "    <Property id=\"XISF:ChecksumAlgorithms\" type=\"String\">sha-1</Property>\n",
        );
    }

    #[test]
    fn add_property_and_add_properties_accumulate() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_property(XisfProperty::int32("A", 1));
        w.add_properties(vec![
            XisfProperty::int32("B", 2),
            XisfProperty::int32("C", 3),
        ]);
        assert_eq!(w.properties.len(), 3);
        assert_eq!(w.properties[0].id, "A");
        assert_eq!(w.properties[2].id, "C");
    }

    #[test]
    fn roundtrip_property_data_blocks_for_vector_and_matrix() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.add_property(XisfProperty::f64_vector("Vec", vec![1.0, 2.0, 3.0]));
        w.add_property(XisfProperty::f64_matrix(
            "Mat",
            2,
            2,
            vec![1.0, 0.0, 0.0, 1.0],
        ));
        w.add_property(XisfProperty::float64("Scalar", 0.5));
        w.add_property(XisfProperty::boolean("Flag", true));
        w.add_property(XisfProperty::int32("Count", 42));
        w.add_property(XisfProperty::string("Name", "value"));

        // Just verify write succeeds and produces a readable file.
        let bytes = w.write_to_vec().unwrap();
        assert!(bytes.starts_with(b"XISF0100"));
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.num_images(), 1);
    }

    #[test]
    fn write_to_vec_starts_with_xisf_signature() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        let bytes = w.write_to_vec().unwrap();
        assert!(bytes.starts_with(b"XISF0100"));
    }

    #[test]
    fn create_writes_to_filesystem() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.xisf");
        let mut w = XisfWriter::create(&path).unwrap();
        w.add_image::<u8>(&[0u8; 16], 4, 4, 1).unwrap();
        w.write().unwrap();

        let mut reader = XisfFile::new(std::fs::File::open(&path).unwrap()).unwrap();
        assert_eq!(reader.num_images(), 1);
        let bytes = reader.read_image_data_raw(0).unwrap();
        assert_eq!(bytes, vec![0u8; 16]);
    }

    #[test]
    fn write_flushes_underlying_writer() {
        // write() takes self by value and calls flush. Using Cursor keeps bytes in memory —
        // we just confirm the call path compiles and succeeds.
        let w = XisfWriter::new(Cursor::new(Vec::new()));
        assert!(w.write().is_ok());
    }

    #[test]
    fn multiple_images_get_consecutive_locations() {
        let mut w = XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[1u8; 16], 4, 4, 1).unwrap();
        w.add_image::<u8>(&[2u8; 16], 4, 4, 1).unwrap();
        let bytes = w.write_to_vec().unwrap();
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.num_images(), 2);
        let first = reader.image_info(0).unwrap().location.offset;
        let second = reader.image_info(1).unwrap().location.offset;
        assert!(second > first);
    }

    // XISF 1.0 geometry dimensions are at least 1, and a reader refuses 0.
    #[test]
    fn add_image_refuses_a_zero_dimension() {
        let mut writer = XisfWriter::new(Cursor::new(Vec::new()));
        for (width, height) in [(0, 0), (0, 4), (4, 0)] {
            assert!(matches!(
                writer.add_image::<u8>(&[], width, height, 1),
                Err(crate::xisf::errors::XisfError::InvalidFormat(_))
            ));
        }
    }

    #[test]
    fn writer_without_images_still_produces_valid_file() {
        let w = XisfWriter::new(Cursor::new(Vec::new()));
        let bytes = w.write_to_vec().unwrap();
        assert!(bytes.starts_with(b"XISF0100"));
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.num_images(), 0);
    }
}
