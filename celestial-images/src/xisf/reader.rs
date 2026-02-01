use crate::fits::header::{Keyword, KeywordValue};
use crate::xisf::header::{
    parse_geometry, ColorSpace, DataLocation, ImageInfo, PixelStorage, SampleFormat,
    XisfCompression, XisfHeader,
};
use crate::xisf::{Result, XisfError};
use byteorder::{LittleEndian, ReadBytesExt};
use quick_xml::events::{BytesEnd, BytesStart, Event};
use quick_xml::Reader;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

const XISF_SIGNATURE: &[u8] = b"XISF0100";

fn strip_fits_quotes(s: &str) -> String {
    let trimmed = s.trim();
    if trimmed.starts_with('\'') && trimmed.ends_with('\'') && trimmed.len() >= 2 {
        trimmed[1..trimmed.len() - 1].trim_end().to_string()
    } else {
        trimmed.to_string()
    }
}

fn parse_keyword_value(s: &str) -> Option<KeywordValue> {
    if s.is_empty() {
        return None;
    }

    // Boolean
    if s == "T" {
        return Some(KeywordValue::Logical(true));
    }
    if s == "F" {
        return Some(KeywordValue::Logical(false));
    }

    // Integer (try first since it's more specific)
    if let Ok(i) = s.parse::<i64>() {
        return Some(KeywordValue::Integer(i));
    }

    // Real (float)
    if let Ok(f) = s.parse::<f64>() {
        return Some(KeywordValue::Real(f));
    }

    // String (everything else)
    Some(KeywordValue::String(s.to_string()))
}

#[derive(Debug)]
pub struct XisfFile<R> {
    reader: R,
    header: XisfHeader,
}

impl XisfFile<File> {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Self::new(file)
    }
}

impl<R: Read + Seek> XisfFile<R> {
    pub fn new(mut reader: R) -> Result<Self> {
        let mut signature = [0u8; 8];
        reader.read_exact(&mut signature)?;

        if signature != XISF_SIGNATURE {
            return Err(XisfError::InvalidFormat(
                "Missing XISF signature".to_string(),
            ));
        }

        let header = Self::parse_header(&mut reader)?;

        Ok(Self { reader, header })
    }

    pub fn num_images(&self) -> usize {
        self.header.images.len()
    }

    pub fn image_info(&self, index: usize) -> Option<&ImageInfo> {
        self.header.images.get(index)
    }

    pub fn keywords(&self) -> &[Keyword] {
        &self.header.keywords
    }

    pub fn get_keyword(&self, key: &str) -> Option<&Keyword> {
        self.header.keywords.iter().find(|k| k.name == key)
    }

    #[deprecated(note = "Use keywords() instead")]
    pub fn fits_keywords(&self) -> &[Keyword] {
        &self.header.keywords
    }

    #[deprecated(note = "Use get_keyword() instead")]
    pub fn get_fits_keyword(&self, key: &str) -> Option<&Keyword> {
        self.header.keywords.iter().find(|k| k.name == key)
    }

    fn parse_header(reader: &mut R) -> Result<XisfHeader> {
        let header_length = reader.read_u32::<LittleEndian>()? as usize;
        let _reserved = reader.read_u32::<LittleEndian>()?;

        let mut xml_content = vec![0u8; header_length];
        reader.read_exact(&mut xml_content)?;

        let xml_end = xml_content
            .iter()
            .position(|&b| b == 0)
            .unwrap_or(xml_content.len());

        let xml_str = std::str::from_utf8(&xml_content[..xml_end])
            .map_err(|e| XisfError::XmlParse(format!("Invalid UTF-8: {}", e)))?;

        Self::parse_xml(xml_str)
    }

    fn parse_xml(xml: &str) -> Result<XisfHeader> {
        let mut reader = Reader::from_str(xml);
        reader.trim_text(true);

        let mut version = String::new();
        let mut images = Vec::new();
        let mut keywords = Vec::new();
        let mut current_image: Option<ImageInfo> = None;
        let mut buf = Vec::new();

        loop {
            match reader.read_event_into(&mut buf) {
                Ok(Event::Start(ref e)) => {
                    Self::handle_xml_element(e, &mut version, &mut current_image, &mut keywords)?
                }
                Ok(Event::Empty(ref e)) => {
                    Self::handle_xml_element(e, &mut version, &mut current_image, &mut keywords)?;
                    if e.name().as_ref() == b"Image" {
                        if let Some(img) = current_image.take() {
                            images.push(img);
                        }
                    }
                }
                Ok(Event::End(ref e)) => {
                    Self::handle_xml_end_element(e, &mut current_image, &mut images)
                }
                Ok(Event::Eof) => break,
                Err(e) => return Err(XisfError::XmlParse(e.to_string())),
                _ => {}
            }
            buf.clear();
        }

        Ok(XisfHeader {
            version,
            images,
            keywords,
        })
    }

    fn handle_xml_element(
        element: &BytesStart,
        version: &mut String,
        current_image: &mut Option<ImageInfo>,
        keywords: &mut Vec<Keyword>,
    ) -> Result<()> {
        match element.name().as_ref() {
            b"xisf" => Self::parse_xisf_attributes(element, version),
            b"Image" => {
                *current_image = Self::parse_image_element(element)?;
                Ok(())
            }
            b"FITSKeyword" => Self::parse_fits_keyword_element(element, keywords),
            _ => Ok(()),
        }
    }

    fn handle_xml_end_element(
        element: &BytesEnd,
        current_image: &mut Option<ImageInfo>,
        images: &mut Vec<ImageInfo>,
    ) {
        if element.name().as_ref() == b"Image" {
            if let Some(img) = current_image.take() {
                images.push(img);
            }
        }
    }

    fn parse_xisf_attributes(element: &BytesStart, version: &mut String) -> Result<()> {
        for attr in element.attributes() {
            let attr = attr.map_err(|e| XisfError::XmlParse(e.to_string()))?;
            if attr.key.as_ref() == b"version" {
                *version = String::from_utf8_lossy(&attr.value).to_string();
            }
        }
        Ok(())
    }

    fn parse_image_element(element: &BytesStart) -> Result<Option<ImageInfo>> {
        let mut geometry = Vec::new();
        let mut sample_format = SampleFormat::Float32;
        let mut bounds = (0.0, 1.0);
        let mut color_space = ColorSpace::Gray;
        let mut pixel_storage = PixelStorage::Planar;
        let mut location = None;
        let mut compression = XisfCompression::None;
        let mut uncompressed_size = None;

        for attr in element.attributes() {
            let attr = attr.map_err(|e| XisfError::XmlParse(e.to_string()))?;
            let key = attr.key.as_ref();
            let value = String::from_utf8_lossy(&attr.value);

            match key {
                b"geometry" => geometry = parse_geometry(&value)?,
                b"sampleFormat" => sample_format = SampleFormat::parse(&value)?,
                b"bounds" => bounds = Self::parse_bounds(&value)?,
                b"colorSpace" => color_space = ColorSpace::parse(&value),
                b"pixelStorage" => pixel_storage = PixelStorage::parse(&value),
                b"location" => location = Some(DataLocation::parse(&value)?),
                b"compression" => {
                    let (comp, size) = Self::parse_compression(&value);
                    compression = comp;
                    uncompressed_size = size;
                }
                _ => {}
            }
        }

        Ok(location.map(|loc| ImageInfo {
            geometry,
            sample_format,
            bounds,
            color_space,
            pixel_storage,
            location: loc,
            compression,
            uncompressed_size,
        }))
    }

    fn parse_bounds(value: &str) -> Result<(f64, f64)> {
        let parts: Vec<&str> = value.split(':').collect();
        if parts.len() == 2 {
            let lower = parts[0].parse().unwrap_or(0.0);
            let upper = parts[1].parse().unwrap_or(1.0);
            Ok((lower, upper))
        } else {
            Ok((0.0, 1.0))
        }
    }

    fn parse_compression(value: &str) -> (XisfCompression, Option<u64>) {
        // Format: "algorithm:uncompressed_size" e.g., "lz4:1048576"
        let parts: Vec<&str> = value.split(':').collect();
        let compression = XisfCompression::parse(parts.first().unwrap_or(&""));
        let uncompressed_size = parts.get(1).and_then(|s| s.parse().ok());
        (compression, uncompressed_size)
    }

    fn parse_fits_keyword_element(element: &BytesStart, keywords: &mut Vec<Keyword>) -> Result<()> {
        let mut name = String::new();
        let mut value_str = String::new();
        let mut comment = String::new();

        for attr in element.attributes() {
            let attr = attr.map_err(|e| XisfError::XmlParse(e.to_string()))?;
            let key = attr.key.as_ref();
            let val = String::from_utf8_lossy(&attr.value);

            match key {
                b"name" => name = val.to_string(),
                b"value" => value_str = strip_fits_quotes(&val),
                b"comment" => comment = val.to_string(),
                _ => {}
            }
        }

        if !name.is_empty() {
            let keyword_value = parse_keyword_value(&value_str);
            let mut kw = Keyword::new(name);
            if let Some(v) = keyword_value {
                kw = kw.with_value(v);
            }
            if !comment.is_empty() {
                kw = kw.with_comment(comment);
            }
            keywords.push(kw);
        }
        Ok(())
    }

    pub fn read_image_data_raw(&mut self, index: usize) -> Result<Vec<u8>> {
        let image_info = self.get_image_info_clone(index)?;
        let buffer = self.read_raw_bytes(&image_info)?;

        if image_info.pixel_storage == PixelStorage::Normal && image_info.num_channels() > 1 {
            let bytes_per_sample = image_info.sample_format.bytes_per_sample();
            Ok(deinterleave_normal_storage(
                &buffer,
                &image_info,
                bytes_per_sample,
            ))
        } else {
            Ok(buffer)
        }
    }

    fn get_image_info_clone(&self, index: usize) -> Result<ImageInfo> {
        self.header
            .images
            .get(index)
            .cloned()
            .ok_or_else(|| XisfError::InvalidFormat(format!("Image {} not found", index)))
    }

    fn read_raw_bytes(&mut self, info: &ImageInfo) -> Result<Vec<u8>> {
        self.reader.seek(SeekFrom::Start(info.location.offset))?;
        let mut buffer = vec![0u8; info.location.size as usize];
        self.reader.read_exact(&mut buffer)?;

        // Decompress if needed
        match info.compression {
            XisfCompression::None => Ok(buffer),
            XisfCompression::Lz4 | XisfCompression::Lz4Hc => {
                lz4_flex::decompress_size_prepended(&buffer).map_err(|e| {
                    XisfError::InvalidFormat(format!("LZ4 decompression failed: {}", e))
                })
            }
            XisfCompression::Zlib => {
                use flate2::read::ZlibDecoder;
                let mut decoder = ZlibDecoder::new(&buffer[..]);
                let mut decompressed = Vec::new();
                decoder.read_to_end(&mut decompressed)?;
                Ok(decompressed)
            }
            XisfCompression::Zstd => Err(XisfError::InvalidFormat(
                "Zstd decompression not yet implemented".to_string(),
            )),
        }
    }

    pub fn read_image_data_typed<T>(&mut self, index: usize) -> Result<Vec<T>>
    where
        T: crate::fits::data::array::DataArray,
    {
        let raw_bytes = self.read_image_data_raw(index)?;
        let byte_order = crate::core::ByteOrder::LittleEndian;
        T::from_bytes(&raw_bytes, byte_order).map_err(|e| XisfError::InvalidFormat(e.to_string()))
    }

    pub fn read_image_data<T>(&mut self, index: usize) -> Result<Vec<T>>
    where
        T: Clone + Default,
    {
        let image_info = self
            .header
            .images
            .get(index)
            .ok_or_else(|| XisfError::InvalidFormat(format!("Image {} not found", index)))?;

        let total_pixels: usize = image_info.geometry.iter().product();
        Ok(vec![T::default(); total_pixels])
    }

    pub fn header(&self) -> &XisfHeader {
        &self.header
    }
}

fn deinterleave_normal_storage(data: &[u8], info: &ImageInfo, bytes_per_sample: usize) -> Vec<u8> {
    let num_channels = info.num_channels();
    let pixels_per_channel = info.pixels_per_channel();
    let channel_size = pixels_per_channel * bytes_per_sample;
    let mut output = vec![0u8; data.len()];

    for pixel_idx in 0..pixels_per_channel {
        for channel in 0..num_channels {
            let src_offset = (pixel_idx * num_channels + channel) * bytes_per_sample;
            let dst_offset = channel * channel_size + pixel_idx * bytes_per_sample;
            copy_sample(&mut output, data, dst_offset, src_offset, bytes_per_sample);
        }
    }
    output
}

fn copy_sample(dst: &mut [u8], src: &[u8], dst_offset: usize, src_offset: usize, len: usize) {
    dst[dst_offset..dst_offset + len].copy_from_slice(&src[src_offset..src_offset + len]);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn create_valid_xisf_xml() -> String {
        r#"<?xml version="1.0" encoding="UTF-8"?>
<xisf version="1.0" xmlns="http://www.pixinsight.com/xisf" 
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" 
      xsi:schemaLocation="http://www.pixinsight.com/xisf http://pixinsight.com/xisf/xisf-1.0.xsd">
<Image geometry="1920:1080:3" sampleFormat="UInt16" bounds="0:65535" 
       colorSpace="RGB" location="attachment:1024:12441600" />
<FITSKeyword name="TELESCOP" value="Hubble" />
<FITSKeyword name="FILTER" value="F606W" />
</xisf>"#
            .to_string()
    }

    fn create_xisf_data_with_xml(xml: &str) -> Vec<u8> {
        use byteorder::{LittleEndian, WriteBytesExt};

        let xml_bytes = xml.as_bytes();
        let header_len = xml_bytes.len() as u32;

        let mut data = Vec::new();
        data.extend_from_slice(XISF_SIGNATURE);
        data.write_u32::<LittleEndian>(header_len).unwrap();
        data.write_u32::<LittleEndian>(0).unwrap();
        data.extend_from_slice(xml_bytes);
        data
    }

    fn create_minimal_xisf() -> Vec<u8> {
        let xml = r#"<?xml version="1.0"?>
<xisf version="1.0">
</xisf>"#;
        create_xisf_data_with_xml(xml)
    }

    #[test]
    fn xisf_file_new_valid_signature() {
        let data = create_minimal_xisf();
        let cursor = Cursor::new(data);
        let result = XisfFile::new(cursor);
        assert!(result.is_ok());
    }

    #[test]
    fn xisf_file_new_invalid_signature() {
        let mut data = Vec::new();
        data.extend_from_slice(b"INVALID0");
        data.extend_from_slice(b"<?xml version=\"1.0\"?><xisf version=\"1.0\"></xisf>");

        let cursor = Cursor::new(data);
        let result = XisfFile::new(cursor);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), XisfError::InvalidFormat(_)));
    }

    #[test]
    fn xisf_file_new_truncated_signature() {
        let data = vec![b'X', b'I', b'S', b'F'];
        let cursor = Cursor::new(data);
        let result = XisfFile::new(cursor);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), XisfError::Io(_)));
    }

    #[test]
    fn parse_xml_valid_complete() {
        let xml = create_valid_xisf_xml();
        let result = XisfFile::<std::io::Cursor<Vec<u8>>>::parse_xml(&xml);
        assert!(result.is_ok());

        let header = result.unwrap();
        assert_eq!(header.version, "1.0");
        assert_eq!(header.images.len(), 1);
        assert_eq!(header.keywords.len(), 2);
        let telescop = header.keywords.iter().find(|k| k.name == "TELESCOP");
        assert!(telescop.is_some());
        assert_eq!(
            telescop.unwrap().value,
            Some(KeywordValue::String("Hubble".to_string()))
        );
    }

    #[test]
    fn parse_xml_missing_version() {
        let xml = r#"<?xml version="1.0"?>
<xisf>
</xisf>"#;
        let result = XisfFile::<std::io::Cursor<Vec<u8>>>::parse_xml(xml);
        assert!(result.is_ok());

        let header = result.unwrap();
        assert_eq!(header.version, "");
    }

    #[test]
    fn parse_xml_multiple_images() {
        let xml = r#"<?xml version="1.0"?>
<xisf version="1.0">
<Image geometry="1920:1080" sampleFormat="UInt16" location="attachment:1024:4147200" />
<Image geometry="960:540:3" sampleFormat="Float32" location="attachment:4148224:6220800" />
</xisf>"#;

        let result = XisfFile::<std::io::Cursor<Vec<u8>>>::parse_xml(xml);
        assert!(result.is_ok());

        let header = result.unwrap();
        assert_eq!(header.images.len(), 2);

        let first_image = &header.images[0];
        assert_eq!(first_image.geometry, vec![1920, 1080]);
        assert!(matches!(first_image.sample_format, SampleFormat::UInt16));

        let second_image = &header.images[1];
        assert_eq!(second_image.geometry, vec![960, 540, 3]);
        assert!(matches!(second_image.sample_format, SampleFormat::Float32));
    }

    #[test]
    fn parse_xml_malformed() {
        let invalid_xmls = [
            "<?xml version=\"1.0\"?><xisf><unclosed>",
            "<?xml version=\"1.0\"?><xisf><Image></Image>",
            "not xml at all",
            "",
            "<?xml version=\"1.0\"?><different_root></different_root>",
        ];

        for invalid_xml in invalid_xmls {
            let result = XisfFile::<std::io::Cursor<Vec<u8>>>::parse_xml(invalid_xml);
            if let Err(e) = result {
                assert!(matches!(e, XisfError::XmlParse(_)));
            }
        }
    }

    #[test]
    fn parse_image_element_complete() {
        let xml = r#"<Image geometry="1920:1080:3" sampleFormat="UInt16" bounds="0:65535" colorSpace="RGB" location="attachment:1024:12441600" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result = XisfFile::<Cursor<Vec<u8>>>::parse_image_element(&element);
            assert!(result.is_ok());

            let image_info = result.unwrap().unwrap();
            assert_eq!(image_info.geometry, vec![1920, 1080, 3]);
            assert!(matches!(image_info.sample_format, SampleFormat::UInt16));
            assert_eq!(image_info.bounds, (0.0, 65535.0));
            assert!(matches!(image_info.color_space, ColorSpace::RGB));
            assert_eq!(image_info.location.offset, 1024);
            assert_eq!(image_info.location.size, 12441600);
        }
    }

    #[test]
    fn parse_image_element_minimal() {
        let xml = r#"<Image location="attachment:0:1024" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result = XisfFile::<Cursor<Vec<u8>>>::parse_image_element(&element);
            assert!(result.is_ok());

            let image_info = result.unwrap().unwrap();
            assert_eq!(image_info.geometry, Vec::<usize>::new());
            assert!(matches!(image_info.sample_format, SampleFormat::Float32));
            assert_eq!(image_info.bounds, (0.0, 1.0));
            assert!(matches!(image_info.color_space, ColorSpace::Gray));
        }
    }

    #[test]
    fn parse_image_element_missing_location() {
        let xml = r#"<Image geometry="1920:1080" sampleFormat="UInt16" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result = XisfFile::<Cursor<Vec<u8>>>::parse_image_element(&element);
            assert!(result.is_ok());
            assert!(result.unwrap().is_none());
        }
    }

    #[test]
    fn parse_bounds_valid() {
        assert_eq!(
            XisfFile::<Cursor<Vec<u8>>>::parse_bounds("0:65535").unwrap(),
            (0.0, 65535.0)
        );
        assert_eq!(
            XisfFile::<Cursor<Vec<u8>>>::parse_bounds("-1.5:1.5").unwrap(),
            (-1.5, 1.5)
        );
        assert_eq!(
            XisfFile::<Cursor<Vec<u8>>>::parse_bounds("0.0:1.0").unwrap(),
            (0.0, 1.0)
        );
    }

    #[test]
    fn parse_bounds_invalid() {
        assert_eq!(
            XisfFile::<Cursor<Vec<u8>>>::parse_bounds("not_a_number").unwrap(),
            (0.0, 1.0)
        );
        assert_eq!(
            XisfFile::<Cursor<Vec<u8>>>::parse_bounds("0:1:2").unwrap(),
            (0.0, 1.0)
        );
        assert_eq!(
            XisfFile::<Cursor<Vec<u8>>>::parse_bounds("").unwrap(),
            (0.0, 1.0)
        );
    }

    #[test]
    fn parse_bounds_partial_invalid() {
        assert_eq!(
            XisfFile::<Cursor<Vec<u8>>>::parse_bounds("5:invalid").unwrap(),
            (5.0, 1.0)
        );
        assert_eq!(
            XisfFile::<Cursor<Vec<u8>>>::parse_bounds("invalid:10").unwrap(),
            (0.0, 10.0)
        );
    }

    #[test]
    fn parse_fits_keyword_element() {
        let xml = r#"<FITSKeyword name="TELESCOP" value="Hubble Space Telescope" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result =
                XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
            assert!(result.is_ok());
            assert_eq!(keywords.len(), 1);
            assert_eq!(keywords[0].name, "TELESCOP");
            assert_eq!(
                keywords[0].value,
                Some(KeywordValue::String("Hubble Space Telescope".to_string()))
            );
        }
    }

    #[test]
    fn parse_fits_keyword_element_with_comment() {
        let xml =
            r#"<FITSKeyword name="EXPTIME" value="300.0" comment="Exposure time in seconds" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result =
                XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
            assert!(result.is_ok());
            assert_eq!(keywords.len(), 1);
            assert_eq!(keywords[0].name, "EXPTIME");
            assert_eq!(keywords[0].value, Some(KeywordValue::Real(300.0)));
            assert_eq!(
                keywords[0].comment,
                Some("Exposure time in seconds".to_string())
            );
        }
    }

    #[test]
    fn parse_fits_keyword_element_history() {
        let xml = r#"<FITSKeyword name="HISTORY" value="" comment="Dark frame subtracted" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result =
                XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
            assert!(result.is_ok());
            assert_eq!(keywords.len(), 1);
            assert_eq!(keywords[0].name, "HISTORY");
            assert_eq!(keywords[0].value, None);
            assert_eq!(
                keywords[0].comment,
                Some("Dark frame subtracted".to_string())
            );
        }
    }

    #[test]
    fn parse_fits_keyword_element_missing_attributes() {
        let xml = r#"<FITSKeyword name="FILTER" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result =
                XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
            assert!(result.is_ok());
            assert_eq!(keywords.len(), 1);
            assert_eq!(keywords[0].name, "FILTER");
            assert_eq!(keywords[0].value, None);
        }
    }

    #[test]
    fn parse_fits_keyword_element_empty_name() {
        let xml = r#"<FITSKeyword value="some_value" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result =
                XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
            assert!(result.is_ok());
            assert!(keywords.is_empty());
        }
    }

    #[test]
    fn parse_fits_keyword_strips_quotes() {
        // PixInsight writes FITS-style quoted strings in XISF
        let xml = r#"<FITSKeyword name="TELESCOP" value="'C14 / Hyperstar V3'" comment="The telescope" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result =
                XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
            assert!(result.is_ok());
            assert_eq!(keywords.len(), 1);
            assert_eq!(keywords[0].name, "TELESCOP");
            assert_eq!(
                keywords[0].value,
                Some(KeywordValue::String("C14 / Hyperstar V3".to_string()))
            );
            assert_eq!(keywords[0].comment, Some("The telescope".to_string()));
        }
    }

    #[test]
    fn parse_fits_keyword_strips_quotes_with_braces() {
        // UUID-style values from PixInsight
        let xml = r#"<FITSKeyword name="SBUUID" value="'{c1810877-52f1-4d56-8065-d1bb6c9c8a32}'" comment="UUID" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result =
                XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
            assert!(result.is_ok());
            assert_eq!(
                keywords[0].value,
                Some(KeywordValue::String(
                    "{c1810877-52f1-4d56-8065-d1bb6c9c8a32}".to_string()
                ))
            );
        }
    }

    #[test]
    fn parse_fits_keyword_numeric_not_stripped() {
        // Numeric values should not be affected
        let xml = r#"<FITSKeyword name="FOCALLEN" value="675." comment="Focal length" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result =
                XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
            assert!(result.is_ok());
            assert_eq!(keywords[0].value, Some(KeywordValue::Real(675.0)));
        }
    }

    #[test]
    fn xisf_file_methods() {
        let data = create_xisf_data_with_xml(&create_valid_xisf_xml());
        let cursor = Cursor::new(data);
        let xisf_file = XisfFile::new(cursor).unwrap();

        assert_eq!(xisf_file.num_images(), 1);
        assert!(xisf_file.image_info(0).is_some());
        assert!(xisf_file.image_info(999).is_none());

        let keywords = xisf_file.keywords();
        assert_eq!(keywords.len(), 2);
        let telescop = xisf_file.get_keyword("TELESCOP");
        assert!(telescop.is_some());
        assert_eq!(
            telescop.unwrap().value,
            Some(KeywordValue::String("Hubble".to_string()))
        );
        assert!(xisf_file.get_keyword("NONEXISTENT").is_none());
    }

    #[test]
    fn read_image_data_valid_index() {
        let data = create_xisf_data_with_xml(&create_valid_xisf_xml());
        let cursor = Cursor::new(data);
        let mut xisf_file = XisfFile::new(cursor).unwrap();

        let result: Result<Vec<f32>> = xisf_file.read_image_data(0);
        assert!(result.is_ok());

        let data = result.unwrap();
        assert_eq!(data.len(), 1920 * 1080 * 3);
        assert!(data.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn read_image_data_invalid_index() {
        let data = create_xisf_data_with_xml(&create_valid_xisf_xml());
        let cursor = Cursor::new(data);
        let mut xisf_file = XisfFile::new(cursor).unwrap();

        let result: Result<Vec<f32>> = xisf_file.read_image_data(999);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), XisfError::InvalidFormat(_)));
    }

    #[test]
    fn handle_xml_element_coverage() {
        let mut version = String::new();
        let mut current_image = None;
        let mut keywords = Vec::new();

        let xml = r#"<UnknownElement attr="value" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result = XisfFile::<Cursor<Vec<u8>>>::handle_xml_element(
                &element,
                &mut version,
                &mut current_image,
                &mut keywords,
            );
            assert!(result.is_ok());
            assert!(version.is_empty());
            assert!(current_image.is_none());
            assert!(keywords.is_empty());
        }
    }

    #[test]
    fn handle_xml_end_element() {
        let mut current_image = Some(ImageInfo {
            geometry: vec![100, 100],
            sample_format: SampleFormat::UInt8,
            bounds: (0.0, 255.0),
            color_space: ColorSpace::Gray,
            pixel_storage: PixelStorage::Planar,
            location: DataLocation {
                offset: 0,
                size: 10000,
            },
            compression: XisfCompression::None,
            uncompressed_size: None,
        });
        let mut images = Vec::new();

        let xml = r#"</Image>"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();

        if let Ok(Event::End(element)) = reader.read_event_into(&mut buf) {
            XisfFile::<Cursor<Vec<u8>>>::handle_xml_end_element(
                &element,
                &mut current_image,
                &mut images,
            );
            assert!(current_image.is_none());
            assert_eq!(images.len(), 1);
        }
    }

    #[test]
    fn parse_xisf_attributes() {
        let xml = r#"<xisf version="1.0" other_attr="ignored" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut version = String::new();

        if let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) {
            let result = XisfFile::<Cursor<Vec<u8>>>::parse_xisf_attributes(&element, &mut version);
            assert!(result.is_ok());
            assert_eq!(version, "1.0");
        }
    }

    #[test]
    fn parse_header_truncated_length() {
        use byteorder::{LittleEndian, WriteBytesExt};

        let mut data = Vec::new();
        data.extend_from_slice(XISF_SIGNATURE);
        data.write_u32::<LittleEndian>(1000).unwrap();
        data.write_u32::<LittleEndian>(0).unwrap();
        data.extend_from_slice(&[b'x'; 100]);

        let cursor = Cursor::new(data);
        let result = XisfFile::new(cursor);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), XisfError::Io(_)));
    }

    #[test]
    fn parse_header_xml_too_large() {
        let data = create_minimal_xisf();
        let cursor = Cursor::new(data);
        let result = XisfFile::new(cursor);
        assert!(result.is_ok());
    }

    #[test]
    fn parse_header_invalid_utf8_in_xml() {
        use byteorder::{LittleEndian, WriteBytesExt};

        let invalid_xml = b"<?xml version=\"1.0\"?><xisf\xff></xisf>";
        let mut data = Vec::new();
        data.extend_from_slice(XISF_SIGNATURE);
        data.write_u32::<LittleEndian>(invalid_xml.len() as u32)
            .unwrap();
        data.write_u32::<LittleEndian>(0).unwrap();
        data.extend_from_slice(invalid_xml);

        let cursor = Cursor::new(data);
        let result = XisfFile::new(cursor);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), XisfError::XmlParse(_)));
    }

    #[test]
    fn extreme_geometry_values() {
        let xml = r#"<?xml version="1.0"?>
<xisf version="1.0">
<Image geometry="65535:65535" sampleFormat="UInt8" location="attachment:0:4294836225" />
</xisf>"#;

        let result = XisfFile::<std::io::Cursor<Vec<u8>>>::parse_xml(xml);
        assert!(result.is_ok());

        let header = result.unwrap();
        assert_eq!(header.images.len(), 1);
        assert_eq!(header.images[0].geometry, vec![65535, 65535]);

        let total_pixels: u64 = header.images[0]
            .geometry
            .iter()
            .map(|&x| x as u64)
            .product();
        assert_eq!(total_pixels, 65535u64 * 65535u64);
    }

    fn create_xisf_with_u8_data(width: usize, height: usize, pixel_data: &[u8]) -> Vec<u8> {
        use byteorder::{LittleEndian, WriteBytesExt};

        let data_size = pixel_data.len();
        let header_block_size = 1024usize;
        let data_offset = 16 + header_block_size;

        let xml = format!(
            r#"<?xml version="1.0"?>
<xisf version="1.0">
<Image geometry="{}:{}" sampleFormat="UInt8" bounds="0:255"
       colorSpace="Gray" location="attachment:{}:{}" />
</xisf>"#,
            width, height, data_offset, data_size
        );

        let xml_bytes = xml.as_bytes();
        let mut data = Vec::new();
        data.extend_from_slice(XISF_SIGNATURE);
        data.write_u32::<LittleEndian>(header_block_size as u32)
            .unwrap();
        data.write_u32::<LittleEndian>(0).unwrap();
        data.extend_from_slice(xml_bytes);
        data.resize(16 + header_block_size, 0);
        data.extend_from_slice(pixel_data);
        data
    }

    fn create_xisf_with_u16_data(width: usize, height: usize, pixel_data: &[u16]) -> Vec<u8> {
        use byteorder::{LittleEndian, WriteBytesExt};

        let data_size = pixel_data.len() * 2;
        let header_block_size = 1024usize;
        let data_offset = 16 + header_block_size;

        let xml = format!(
            r#"<?xml version="1.0"?>
<xisf version="1.0">
<Image geometry="{}:{}" sampleFormat="UInt16" bounds="0:65535"
       colorSpace="Gray" location="attachment:{}:{}" />
</xisf>"#,
            width, height, data_offset, data_size
        );

        let xml_bytes = xml.as_bytes();
        let mut data = Vec::new();
        data.extend_from_slice(XISF_SIGNATURE);
        data.write_u32::<LittleEndian>(header_block_size as u32)
            .unwrap();
        data.write_u32::<LittleEndian>(0).unwrap();
        data.extend_from_slice(xml_bytes);
        data.resize(16 + header_block_size, 0);

        for &val in pixel_data {
            data.write_u16::<LittleEndian>(val).unwrap();
        }
        data
    }

    #[test]
    fn read_image_data_typed_u8() {
        let pixels = vec![10u8, 20, 30, 40, 50, 60];
        let xisf_data = create_xisf_with_u8_data(3, 2, &pixels);
        let cursor = Cursor::new(xisf_data);
        let mut xisf_file = XisfFile::new(cursor).unwrap();

        let result: Result<Vec<u8>> = xisf_file.read_image_data_typed(0);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), pixels);
    }

    #[test]
    fn read_image_data_typed_u16_as_i16() {
        let pixels = vec![1000u16, 2000, 3000, 4000];
        let xisf_data = create_xisf_with_u16_data(2, 2, &pixels);
        let cursor = Cursor::new(xisf_data);
        let mut xisf_file = XisfFile::new(cursor).unwrap();

        let result: Result<Vec<i16>> = xisf_file.read_image_data_typed(0);
        assert!(result.is_ok());
        let data = result.unwrap();
        assert_eq!(data.len(), 4);
        assert_eq!(data[0], 1000);
        assert_eq!(data[1], 2000);
        assert_eq!(data[2], 3000);
        assert_eq!(data[3], 4000);
    }

    #[test]
    fn read_image_data_typed_invalid_index() {
        let pixels = vec![1u8, 2, 3, 4];
        let xisf_data = create_xisf_with_u8_data(2, 2, &pixels);
        let cursor = Cursor::new(xisf_data);
        let mut xisf_file = XisfFile::new(cursor).unwrap();

        let result: Result<Vec<u8>> = xisf_file.read_image_data_typed(999);
        assert!(result.is_err());
    }

    #[test]
    fn deinterleave_normal_storage_rgb_u8() {
        let info = ImageInfo {
            geometry: vec![2, 2, 3],
            sample_format: SampleFormat::UInt8,
            bounds: (0.0, 255.0),
            color_space: ColorSpace::RGB,
            pixel_storage: PixelStorage::Normal,
            location: DataLocation {
                offset: 0,
                size: 12,
            },
            compression: XisfCompression::None,
            uncompressed_size: None,
        };

        let normal_data = vec![
            1, 2, 3, // pixel 0: R0, G0, B0
            4, 5, 6, // pixel 1: R1, G1, B1
            7, 8, 9, // pixel 2: R2, G2, B2
            10, 11, 12, // pixel 3: R3, G3, B3
        ];

        let planar = deinterleave_normal_storage(&normal_data, &info, 1);

        let expected_planar = vec![
            1, 4, 7, 10, // R channel
            2, 5, 8, 11, // G channel
            3, 6, 9, 12, // B channel
        ];
        assert_eq!(planar, expected_planar);
    }

    #[test]
    fn deinterleave_normal_storage_rgb_u16() {
        let info = ImageInfo {
            geometry: vec![2, 1, 3],
            sample_format: SampleFormat::UInt16,
            bounds: (0.0, 65535.0),
            color_space: ColorSpace::RGB,
            pixel_storage: PixelStorage::Normal,
            location: DataLocation {
                offset: 0,
                size: 12,
            },
            compression: XisfCompression::None,
            uncompressed_size: None,
        };

        let normal_data: Vec<u8> = vec![
            0x01, 0x00, // R0 = 1
            0x02, 0x00, // G0 = 2
            0x03, 0x00, // B0 = 3
            0x04, 0x00, // R1 = 4
            0x05, 0x00, // G1 = 5
            0x06, 0x00, // B1 = 6
        ];

        let planar = deinterleave_normal_storage(&normal_data, &info, 2);

        let expected_planar: Vec<u8> = vec![
            0x01, 0x00, 0x04, 0x00, // R channel
            0x02, 0x00, 0x05, 0x00, // G channel
            0x03, 0x00, 0x06, 0x00, // B channel
        ];
        assert_eq!(planar, expected_planar);
    }

    fn create_xisf_with_normal_storage_rgb(width: usize, height: usize, data: &[u8]) -> Vec<u8> {
        use byteorder::{LittleEndian, WriteBytesExt};

        let data_size = data.len();
        let header_block_size = 1024usize;
        let data_offset = 16 + header_block_size;

        let xml = format!(
            r#"<?xml version="1.0"?>
<xisf version="1.0">
<Image geometry="{}:{}:3" sampleFormat="UInt8" bounds="0:255"
       colorSpace="RGB" pixelStorage="Normal" location="attachment:{}:{}" />
</xisf>"#,
            width, height, data_offset, data_size
        );

        let xml_bytes = xml.as_bytes();
        let mut result = Vec::new();
        result.extend_from_slice(XISF_SIGNATURE);
        result
            .write_u32::<LittleEndian>(header_block_size as u32)
            .unwrap();
        result.write_u32::<LittleEndian>(0).unwrap();
        result.extend_from_slice(xml_bytes);
        result.resize(16 + header_block_size, 0);
        result.extend_from_slice(data);
        result
    }

    #[test]
    fn read_normal_storage_converts_to_planar() {
        let normal_data = vec![
            1, 2, 3, // pixel 0: R0, G0, B0
            4, 5, 6, // pixel 1: R1, G1, B1
            7, 8, 9, // pixel 2: R2, G2, B2
            10, 11, 12, // pixel 3: R3, G3, B3
        ];

        let xisf_data = create_xisf_with_normal_storage_rgb(2, 2, &normal_data);
        let cursor = Cursor::new(xisf_data);
        let mut xisf_file = XisfFile::new(cursor).unwrap();

        let info = xisf_file.image_info(0).unwrap();
        assert_eq!(info.pixel_storage, PixelStorage::Normal);

        let result = xisf_file.read_image_data_raw(0).unwrap();

        let expected_planar = vec![
            1, 4, 7, 10, // R channel
            2, 5, 8, 11, // G channel
            3, 6, 9, 12, // B channel
        ];
        assert_eq!(result, expected_planar);
    }
}
