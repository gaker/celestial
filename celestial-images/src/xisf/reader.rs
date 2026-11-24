use crate::fits::header::keywords::{find_keyword, Keyword};
use crate::fits::header::parser::Field;
use crate::xisf::data_block::{past_end_of_file, planar_samples, read_image_block};
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::image_info::{ImageInfo, XisfHeader};
use crate::xisf::header::property::XisfProperty;
use crate::xisf::header_text::HeaderScan;
use crate::xisf::image_element::parse_image_element;
use crate::xisf::keyword_reader::fits_keyword_element;
use crate::xisf::property_reader::image_property_specs;
use crate::xisf::property_reader::spec::{read_properties, PropertySpec};
use crate::xisf::writer::data_type::XisfDataType;
use byteorder::{LittleEndian, ReadBytesExt};
use quick_xml::events::attributes::Attribute;
use quick_xml::events::{BytesEnd, BytesStart, Event};
use quick_xml::{Reader, XmlVersion};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

const XISF_SIGNATURE: &[u8] = b"XISF0100";

#[derive(Debug)]
pub struct XisfFile<R> {
    reader: R,
    header: XisfHeader,
    stream_len: u64,
    properties: Vec<Vec<PropertySpec>>,
}

impl XisfFile<File> {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Self::new(file)
    }
}

impl<R: Read + Seek> XisfFile<R> {
    pub fn new(mut reader: R) -> Result<Self> {
        let stream_len = reader.seek(SeekFrom::End(0))?;
        reader.seek(SeekFrom::Start(0))?;
        let mut signature = [0u8; 8];
        reader.read_exact(&mut signature)?;

        if signature != XISF_SIGNATURE {
            return Err(XisfError::InvalidFormat(
                "Missing XISF signature".to_string(),
            ));
        }

        let xml = Self::read_xml(&mut reader, stream_len)?;
        let header = Self::parse_xml(&xml)?;
        let properties = image_property_specs(&xml)?;

        Ok(Self {
            reader,
            header,
            stream_len,
            properties,
        })
    }

    pub fn num_images(&self) -> usize {
        self.header.images.len()
    }

    pub fn image_info(&self, index: usize) -> Option<&ImageInfo> {
        self.header.images.get(index)
    }

    // The keywords of image 0, as `image_keywords` gives them, or the
    // FITSKeyword children of the root element when the file has no image.
    pub fn keywords(&self) -> &[Keyword] {
        match self.image_keywords(0) {
            Some(keywords) => keywords,
            None => &self.header.root_keywords,
        }
    }

    // The keywords XISF 1.0 section 11.6 associates with image `index`, as
    // `KeywordScan::finish` resolves them. None when the file has no image
    // `index`.
    pub fn image_keywords(&self, index: usize) -> Option<&[Keyword]> {
        self.header.keywords.get(index).map(Vec::as_slice)
    }

    pub fn get_keyword(&self, key: &str) -> Option<&Keyword> {
        let keywords = self.keywords();
        keywords.get(find_keyword(keywords, key)?)
    }

    pub fn get(&self, key: &str) -> Field<'_> {
        Field::from_value(self.get_keyword(key).and_then(|k| k.value.as_ref()))
    }

    // Reads the properties of image `index`, with vector and matrix values
    // from their data blocks.
    pub fn properties(&mut self, index: usize) -> Result<Vec<XisfProperty>> {
        let specs = self
            .properties
            .get(index)
            .ok_or_else(|| XisfError::InvalidFormat(format!("Image {} not found", index)))?;
        Ok(read_properties(&mut self.reader, self.stream_len, specs))
    }

    fn read_xml(reader: &mut R, stream_len: u64) -> Result<String> {
        let header_length = reader.read_u32::<LittleEndian>()?;
        let _reserved = reader.read_u32::<LittleEndian>()?;
        // Checked before allocating, so a forged length cannot request 4 GiB.
        if 16 + u64::from(header_length) > stream_len {
            return Err(past_end_of_file("The XML header", stream_len));
        }

        let mut xml_content = vec![0u8; header_length as usize];
        reader.read_exact(&mut xml_content)?;

        let xml_end = xml_content
            .iter()
            .position(|&b| b == 0)
            .unwrap_or(xml_content.len());

        let xml_str = std::str::from_utf8(&xml_content[..xml_end])
            .map_err(|e| XisfError::XmlParse(format!("Invalid UTF-8: {}", e)))?;

        Ok(xml_str.to_string())
    }

    fn parse_xml(xml: &str) -> Result<XisfHeader> {
        let mut reader = Reader::from_str(xml);
        reader.config_mut().trim_text(true);
        let mut scan = HeaderScan::new();
        loop {
            match reader.read_event() {
                Ok(Event::Start(e)) => Self::start_element(&e, &mut reader, &mut scan)?,
                Ok(Event::Empty(e)) => Self::empty_element(&e, &mut scan)?,
                Ok(Event::End(e)) => Self::end_element(&e, &mut scan),
                Ok(Event::Eof) => return scan.finish(xml.len()),
                Err(e) => return Err(XisfError::XmlParse(e.to_string())),
                _ => {}
            }
        }
    }

    fn start_element(
        element: &BytesStart,
        reader: &mut Reader<&[u8]>,
        scan: &mut HeaderScan,
    ) -> Result<()> {
        if scan.is_image_data(element) {
            return scan.read_data(element, reader);
        }
        Self::keyword_element(element, scan)?;
        scan.enter();
        let (version, image) = scan.targets();
        Self::handle_xml_element(element, version, image)?;
        if element.name().as_ref() == "Image" {
            scan.open_image(element, reader)?;
        }
        Ok(())
    }

    fn empty_element(element: &BytesStart, scan: &mut HeaderScan) -> Result<()> {
        if scan.is_image_data(element) {
            return scan.set_data(element, Some(String::new()));
        }
        Self::keyword_element(element, scan)?;
        let (version, image) = scan.targets();
        Self::handle_xml_element(element, version, image)?;
        if element.name().as_ref() == "Image" {
            scan.close_empty_image(element)?;
            scan.finish_image();
        }
        Ok(())
    }

    fn end_element(element: &BytesEnd, scan: &mut HeaderScan) {
        scan.leave();
        if element.name().as_ref() == "Image" {
            scan.close_image();
            scan.finish_image();
        }
    }

    fn handle_xml_element(
        element: &BytesStart,
        version: &mut String,
        current_image: &mut Option<ImageInfo>,
    ) -> Result<()> {
        match element.name().as_ref() {
            "xisf" => Self::parse_xisf_attributes(element, version),
            "Image" => {
                *current_image = Some(parse_image_element(element)?);
                Ok(())
            }
            _ => Ok(()),
        }
    }

    // XISF 1.0 section 11.6 allows a FITSKeyword element only as a child of
    // an Image element or of the root element. Elsewhere, as in a Metadata
    // or Thumbnail element, it is skipped.
    fn keyword_element(element: &BytesStart, scan: &mut HeaderScan) -> Result<()> {
        let (in_image, at_root) = (scan.in_image(), scan.at_root());
        let keywords = &mut scan.keywords;
        match element.name().as_ref() {
            "FITSKeyword" if in_image => {
                Self::parse_fits_keyword_element(element, &mut keywords.image)
            }
            "FITSKeyword" if at_root => keywords.add_root_keyword(element),
            "Reference" if in_image => keywords.add_reference(element),
            "Reference" => keywords.mark_referenced(element),
            _ => Ok(()),
        }
    }

    fn parse_xisf_attributes(element: &BytesStart, version: &mut String) -> Result<()> {
        for attr in element.attributes() {
            let attr = attr.map_err(|e| XisfError::XmlParse(e.to_string()))?;
            if attr.key.as_ref() == "version" {
                *version = attr.value.to_string();
            }
        }
        Ok(())
    }

    fn parse_fits_keyword_element(element: &BytesStart, keywords: &mut Vec<Keyword>) -> Result<()> {
        keywords.extend(fits_keyword_element(element)?);
        Ok(())
    }

    // Returns the samples as little-endian bytes in planar order.
    pub fn read_image_data_raw(&mut self, index: usize) -> Result<Vec<u8>> {
        let stored = self.read_stored_samples(index)?;
        let info = image_at(&self.header.images, index)?;
        if info.stores_interleaved() {
            return Ok(planar_samples(&stored, info));
        }
        Ok(stored)
    }

    // Returns the samples as little-endian bytes in stored order, which
    // `ImageInfo::stores_interleaved` gives.
    pub(crate) fn read_stored_samples(&mut self, index: usize) -> Result<Vec<u8>> {
        let info = image_at(&self.header.images, index)?;
        read_image_block(&mut self.reader, self.stream_len, info)
    }

    pub fn read_image_data_typed<T: XisfDataType>(&mut self, index: usize) -> Result<Vec<T>> {
        let stored = &image_at(&self.header.images, index)?.sample_format;
        if *stored != T::SAMPLE_FORMAT {
            return Err(XisfError::SampleFormatMismatch {
                stored: stored.as_str(),
                requested: T::SAMPLE_FORMAT.as_str(),
            });
        }
        let raw_bytes = self.read_image_data_raw(index)?;
        Ok(T::from_le_bytes(&raw_bytes))
    }
}

// Borrowed rather than cloned, since an inline or embedded image holds its
// whole decoded block.
fn image_at(images: &[ImageInfo], index: usize) -> Result<&ImageInfo> {
    images
        .get(index)
        .ok_or_else(|| XisfError::InvalidFormat(format!("Image {} not found", index)))
}

// The attribute value with the five predefined XML entities and character
// references decoded (XML 1.0 section 3.3.3).
pub(crate) fn attribute_text(attr: &Attribute) -> Result<String> {
    attr.normalized_value(XmlVersion::Implicit1_0)
        .map(|text| text.into_owned())
        .map_err(|e| XisfError::XmlParse(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::types::ByteOrder;
    use crate::fits::header::keywords::KeywordValue;
    use crate::xisf::header::block_attrs::BlockDefect;
    use crate::xisf::header::color_space::ColorSpace;
    use crate::xisf::header::compression::XisfCompression;
    use crate::xisf::header::location::DataLocation;
    use crate::xisf::header::pixel_storage::PixelStorage;
    use crate::xisf::header::sample_format::SampleFormat;
    use crate::xisf::image_element::parse_bounds;
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
        assert_eq!(header.keywords[0].len(), 2);
        let telescop = header.keywords[0].iter().find(|k| k.name == "TELESCOP");
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
<Image geometry="1920:1080:1" sampleFormat="UInt16" location="attachment:1024:4147200" />
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

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
        let result = parse_image_element(&element);
        assert!(result.is_ok());

        let image_info = result.unwrap();
        assert_eq!(image_info.geometry, vec![1920, 1080, 3]);
        assert!(matches!(image_info.sample_format, SampleFormat::UInt16));
        assert_eq!(image_info.bounds, (0.0, 65535.0));
        assert!(matches!(image_info.color_space, ColorSpace::Rgb));
        assert_eq!(image_info.location.offset, 1024);
        assert_eq!(image_info.location.size, 12441600);
    }

    #[test]
    fn parse_image_element_minimal() {
        let xml = r#"<Image location="attachment:0:1024" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
        let result = parse_image_element(&element);
        assert!(result.is_ok());

        let image_info = result.unwrap();
        assert_eq!(image_info.geometry, Vec::<usize>::new());
        assert!(matches!(image_info.sample_format, SampleFormat::Float32));
        assert_eq!(image_info.bounds, (0.0, 1.0));
        assert!(matches!(image_info.color_space, ColorSpace::Gray));
    }

    #[test]
    fn parse_image_element_missing_location() {
        let xml = r#"<Image geometry="1920:1080" sampleFormat="UInt16" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
        let result = parse_image_element(&element);
        assert_eq!(
            result.unwrap().block_defect,
            Some(BlockDefect::Invalid(
                "The Image element has no location attribute".to_string()
            ))
        );
    }

    #[test]
    fn parse_bounds_valid() {
        assert_eq!(parse_bounds("0:65535").unwrap(), (0.0, 65535.0));
        assert_eq!(parse_bounds("-1.5:1.5").unwrap(), (-1.5, 1.5));
        assert_eq!(parse_bounds("0.0:1.0").unwrap(), (0.0, 1.0));
    }

    #[test]
    fn parse_bounds_invalid() {
        for value in ["not_a_number", "0:1:2", ""] {
            assert!(matches!(
                parse_bounds(value),
                Err(XisfError::InvalidFormat(_))
            ));
        }
    }

    #[test]
    fn parse_bounds_partial_invalid() {
        for value in ["5:invalid", "invalid:10"] {
            assert!(matches!(
                parse_bounds(value),
                Err(XisfError::InvalidFormat(_))
            ));
        }
    }

    #[test]
    fn parse_fits_keyword_element() {
        let xml = r#"<FITSKeyword name="TELESCOP" value="Hubble Space Telescope" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
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

    #[test]
    fn parse_fits_keyword_element_with_comment() {
        let xml =
            r#"<FITSKeyword name="EXPTIME" value="300.0" comment="Exposure time in seconds" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
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

    #[test]
    fn parse_fits_keyword_element_history() {
        let xml = r#"<FITSKeyword name="HISTORY" value="" comment="Dark frame subtracted" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
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

    #[test]
    fn parse_fits_keyword_element_missing_attributes() {
        let xml = r#"<FITSKeyword name="FILTER" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
        let result =
            XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
        assert!(result.is_ok());
        assert_eq!(keywords.len(), 1);
        assert_eq!(keywords[0].name, "FILTER");
        assert_eq!(keywords[0].value, None);
    }

    #[test]
    fn parse_fits_keyword_element_empty_name() {
        let xml = r#"<FITSKeyword value="some_value" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
        let result =
            XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
        assert!(result.is_ok());
        assert!(keywords.is_empty());
    }

    #[test]
    fn parse_fits_keyword_strips_quotes() {
        // PixInsight writes FITS-style quoted strings in XISF
        let xml = r#"<FITSKeyword name="TELESCOP" value="'C14 / Hyperstar V3'" comment="The telescope" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
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

    #[test]
    fn parse_fits_keyword_strips_quotes_with_braces() {
        // UUID-style values from PixInsight
        let xml = r#"<FITSKeyword name="SBUUID" value="'{c1810877-52f1-4d56-8065-d1bb6c9c8a32}'" comment="UUID" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
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

    #[test]
    fn parse_fits_keyword_numeric_not_stripped() {
        // Numeric values should not be affected
        let xml = r#"<FITSKeyword name="FOCALLEN" value="675." comment="Focal length" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut keywords = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
        let result =
            XisfFile::<Cursor<Vec<u8>>>::parse_fits_keyword_element(&element, &mut keywords);
        assert!(result.is_ok());
        assert_eq!(keywords[0].value, Some(KeywordValue::Real(675.0)));
    }

    #[test]
    fn fits_keyword_attributes_decode_xml_entities() {
        let xml = r#"<xisf version="1.0"><Image geometry="1:1:1" sampleFormat="UInt8" location="attachment:0:1"><FITSKeyword name="OBSERVER" value="&apos;O&apos;&apos;Brien &amp; co&apos;" comment="&lt;a&gt; &quot;b&quot; &#65;&#x42;" /></Image></xisf>"#;
        let header = XisfFile::<Cursor<Vec<u8>>>::parse_xml(xml).unwrap();
        assert_eq!(
            header.keywords[0],
            vec![Keyword::string("OBSERVER", "O'Brien & co").with_comment("<a> \"b\" AB")]
        );
    }

    // The FITS reader maps bytes outside 0x20-0x7E to '?', and so does this.
    #[test]
    fn fits_keyword_text_outside_printable_ascii_becomes_question_marks() {
        let xml = "<xisf version=\"1.0\"><Image geometry=\"1:1:1\" sampleFormat=\"UInt8\" location=\"attachment:0:1\"><FITSKeyword name=\"T\u{c9}ST\" value=\"'a\u{e9}&#9;b'\" comment=\"c&#2;d\" /></Image></xisf>";
        let header = XisfFile::<Cursor<Vec<u8>>>::parse_xml(xml).unwrap();
        assert_eq!(
            header.keywords[0],
            vec![Keyword::string("T??ST", "a???b").with_comment("c?d")]
        );
    }

    #[test]
    fn fits_keyword_complex_value_is_read_as_complex() {
        let xml = r#"<xisf version="1.0"><Image geometry="1:1:1" sampleFormat="UInt8" location="attachment:0:1"><FITSKeyword name="CPLX" value="(1.5, -2.0E+01)" comment="" /></Image></xisf>"#;
        let header = XisfFile::<Cursor<Vec<u8>>>::parse_xml(xml).unwrap();
        assert_eq!(
            header.keywords[0][0].value,
            Some(KeywordValue::Complex(1.5, -20.0))
        );
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
    fn get_returns_typed_field() {
        let data = create_xisf_data_with_xml(&create_valid_xisf_xml());
        let cursor = Cursor::new(data);
        let xisf_file = XisfFile::new(cursor).unwrap();

        let telescop = xisf_file.get("TELESCOP");
        assert!(telescop.exists());
        assert_eq!(telescop.as_str(), Some("Hubble"));
        assert_eq!(telescop.as_i64(), None);
        assert_eq!(
            telescop.value(),
            Some(&KeywordValue::String("Hubble".to_string()))
        );

        let missing = xisf_file.get("NONEXISTENT");
        assert!(!missing.exists());
        assert_eq!(missing.as_str(), None);
        assert_eq!(missing.value(), None);
    }

    #[test]
    fn handle_xml_element_coverage() {
        let mut version = String::new();
        let mut current_image = None;

        let xml = r#"<UnknownElement attr="value" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
        let result = XisfFile::<Cursor<Vec<u8>>>::handle_xml_element(
            &element,
            &mut version,
            &mut current_image,
        );
        assert!(result.is_ok());
        assert!(version.is_empty());
        assert!(current_image.is_none());
    }

    #[test]
    fn image_end_tag_lists_the_open_image() {
        let mut scan = HeaderScan::new();
        scan.image = Some(ImageInfo {
            geometry: vec![100, 100],
            sample_format: SampleFormat::UInt8,
            bounds: (0.0, 255.0),
            color_space: ColorSpace::Gray,
            pixel_storage: PixelStorage::Planar,
            location: DataLocation::new(0, 10000),
            compression: XisfCompression::None,
            uncompressed_size: None,
            shuffle_item_size: None,
            subblocks: Vec::new(),
            byte_order: ByteOrder::LittleEndian,
            checksum: None,
            block_defect: None,
        });

        // The reader rejects an end tag with no open element, so the fixture
        // opens the Image element first.
        let xml = r#"<Image></Image>"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();

        let Ok(Event::Start(_)) = reader.read_event_into(&mut buf) else {
            panic!("expected a start tag");
        };
        let Ok(Event::End(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an end tag");
        };
        XisfFile::<Cursor<Vec<u8>>>::end_element(&element, &mut scan);
        assert!(scan.image.is_none());
        assert_eq!(scan.header.images.len(), 1);
        assert_eq!(scan.finish(0).unwrap().keywords.len(), 1);
    }

    #[test]
    fn parse_xisf_attributes() {
        let xml = r#"<xisf version="1.0" other_attr="ignored" />"#;
        let mut reader = Reader::from_str(xml);
        let mut buf = Vec::new();
        let mut version = String::new();

        let Ok(Event::Empty(element)) = reader.read_event_into(&mut buf) else {
            panic!("expected an empty element");
        };
        let result = XisfFile::<Cursor<Vec<u8>>>::parse_xisf_attributes(&element, &mut version);
        assert!(result.is_ok());
        assert_eq!(version, "1.0");
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
<Image geometry="65535:65535:1" sampleFormat="UInt8" location="attachment:0:4294836225" />
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
    fn read_image_data_typed_u16() {
        let pixels = vec![1000u16, 2000, 3000, 65535];
        let xisf_data = create_xisf_with_u16_data(2, 2, &pixels);
        let cursor = Cursor::new(xisf_data);
        let mut xisf_file = XisfFile::new(cursor).unwrap();

        let result: Result<Vec<u16>> = xisf_file.read_image_data_typed(0);
        assert_eq!(result.unwrap(), pixels);
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
    fn planar_samples_of_normal_storage_rgb_u8() {
        let info = ImageInfo {
            geometry: vec![2, 2, 3],
            sample_format: SampleFormat::UInt8,
            bounds: (0.0, 255.0),
            color_space: ColorSpace::Rgb,
            pixel_storage: PixelStorage::Normal,
            location: DataLocation::new(0, 12),
            compression: XisfCompression::None,
            uncompressed_size: None,
            shuffle_item_size: None,
            subblocks: Vec::new(),
            byte_order: ByteOrder::LittleEndian,
            checksum: None,
            block_defect: None,
        };

        let normal_data = vec![
            1, 2, 3, // pixel 0: R0, G0, B0
            4, 5, 6, // pixel 1: R1, G1, B1
            7, 8, 9, // pixel 2: R2, G2, B2
            10, 11, 12, // pixel 3: R3, G3, B3
        ];

        let planar = planar_samples(&normal_data, &info);

        let expected_planar = vec![
            1, 4, 7, 10, // R channel
            2, 5, 8, 11, // G channel
            3, 6, 9, 12, // B channel
        ];
        assert_eq!(planar, expected_planar);
    }

    #[test]
    fn planar_samples_of_normal_storage_rgb_u16() {
        let info = ImageInfo {
            geometry: vec![2, 1, 3],
            sample_format: SampleFormat::UInt16,
            bounds: (0.0, 65535.0),
            color_space: ColorSpace::Rgb,
            pixel_storage: PixelStorage::Normal,
            location: DataLocation::new(0, 12),
            compression: XisfCompression::None,
            uncompressed_size: None,
            shuffle_item_size: None,
            subblocks: Vec::new(),
            byte_order: ByteOrder::LittleEndian,
            checksum: None,
            block_defect: None,
        };

        let normal_data: Vec<u8> = vec![
            0x01, 0x00, // R0 = 1
            0x02, 0x00, // G0 = 2
            0x03, 0x00, // B0 = 3
            0x04, 0x00, // R1 = 4
            0x05, 0x00, // G1 = 5
            0x06, 0x00, // B1 = 6
        ];

        let planar = planar_samples(&normal_data, &info);

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

    fn xisf_with_block(attrs: &str, block: &[u8]) -> Vec<u8> {
        xisf_with_block_size(attrs, block, &block.len().to_string())
    }

    fn xisf_with_block_size(attrs: &str, block: &[u8], size: &str) -> Vec<u8> {
        let header_block_size = 1024usize;
        let xml = format!(
            r#"<?xml version="1.0"?><xisf version="1.0"><Image {attrs} location="attachment:{}:{size}" /></xisf>"#,
            16 + header_block_size
        );
        let mut data = create_xisf_data_with_xml(&xml);
        data[8..12].copy_from_slice(&(header_block_size as u32).to_le_bytes());
        data.resize(16 + header_block_size, 0);
        data.extend_from_slice(block);
        data
    }

    fn read_block(attrs: &str, block: &[u8]) -> Result<Vec<u8>> {
        XisfFile::new(Cursor::new(xisf_with_block(attrs, block)))?.read_image_data_raw(0)
    }

    fn zlib(data: &[u8]) -> Vec<u8> {
        use flate2::write::ZlibEncoder;
        use std::io::Write;
        let mut encoder = ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    }

    // The forward transform from the Shuffle routine in XISF 1.0 section 10.6.3.
    fn spec_shuffle(input: &[u8], item_size: usize) -> Vec<u8> {
        let items = input.len() / item_size;
        let mut output = Vec::with_capacity(input.len());
        for j in 0..item_size {
            for i in 0..items {
                output.push(input[i * item_size + j]);
            }
        }
        output.extend_from_slice(&input[items * item_size..]);
        output
    }

    // An LZ4 block holding one literal-only sequence, the shortest valid
    // block for data that has no matches.
    fn lz4_literal_block(data: &[u8]) -> Vec<u8> {
        let mut block = match data.len().checked_sub(15) {
            None => vec![(data.len() as u8) << 4],
            Some(extra) => vec![0xF0, u8::try_from(extra).unwrap()],
        };
        block.extend_from_slice(data);
        block
    }

    // Hand-encoded LZ4 block for [1, 2, 3, 4] repeated 8 times: 4 literals,
    // a match at offset 4 of length 4 + 15 + 4, then 5 closing literals.
    const LZ4_REPEAT_BLOCK: [u8; 14] = [0x4F, 1, 2, 3, 4, 0x04, 0x00, 0x04, 0x50, 4, 1, 2, 3, 4];

    fn repeat_1234() -> Vec<u8> {
        [1u8, 2, 3, 4].repeat(8)
    }

    const SHUFFLE_CASES: [(&str, usize); 5] = [
        ("UInt16", 2),
        ("UInt32", 4),
        ("Float32", 4),
        ("UInt64", 8),
        ("Float64", 8),
    ];

    fn distinct_bytes(len: usize) -> Vec<u8> {
        (0..len).map(|i| (i * 7 + 3) as u8).collect()
    }

    #[test]
    fn zlib_shuffled_u16_matches_review_example() {
        let raw: Vec<u8> = [0x0102u16, 0x0304, 0x0506, 0x0708]
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        let block = zlib(&[0x02, 0x04, 0x06, 0x08, 0x01, 0x03, 0x05, 0x07]);
        let attrs = r#"geometry="2:2:1" sampleFormat="UInt16" compression="zlib+sh:8:2""#;
        assert_eq!(read_block(attrs, &block).unwrap(), raw);
    }

    #[test]
    fn zlib_shuffled_unshuffles_each_item_size() {
        for (format, item) in SHUFFLE_CASES {
            let raw = distinct_bytes(4 * item);
            let block = zlib(&spec_shuffle(&raw, item));
            let attrs = format!(
                r#"geometry="2:2" sampleFormat="{format}" compression="zlib+sh:{}:{item}""#,
                raw.len()
            );
            assert_eq!(read_block(&attrs, &block).unwrap(), raw, "{format}");
        }
    }

    #[test]
    fn lz4_shuffled_unshuffles_each_item_size() {
        for codec in ["lz4+sh", "lz4hc+sh"] {
            for (format, item) in SHUFFLE_CASES {
                let raw = distinct_bytes(4 * item);
                let block = lz4_literal_block(&spec_shuffle(&raw, item));
                let attrs = format!(
                    r#"geometry="2:2" sampleFormat="{format}" compression="{codec}:{}:{item}""#,
                    raw.len()
                );
                assert_eq!(read_block(&attrs, &block).unwrap(), raw, "{codec} {format}");
            }
        }
    }

    #[test]
    fn lz4_shuffled_block_with_match_unshuffles() {
        let raw: Vec<u8> = [1u8, 1, 2, 2, 3, 3, 4, 4].repeat(4);
        let attrs = r#"geometry="4:4" sampleFormat="UInt16" compression="lz4+sh:32:2""#;
        assert_eq!(read_block(attrs, &LZ4_REPEAT_BLOCK).unwrap(), raw);
    }

    #[test]
    fn shuffle_keeps_trailing_partial_item_in_place() {
        let raw = distinct_bytes(8);
        let block = zlib(&spec_shuffle(&raw, 3));
        let attrs = r#"geometry="2:4" sampleFormat="UInt8" compression="zlib+sh:8:3""#;
        assert_eq!(read_block(attrs, &block).unwrap(), raw);
    }

    #[test]
    fn shuffled_codec_requires_positive_item_size() {
        let block = zlib(&[1, 0, 2, 0]);
        for compression in ["zlib+sh:4", "zlib+sh:4:0", "zlib+sh:4:x", "zlib:4:2"] {
            let attrs =
                format!(r#"geometry="2:1" sampleFormat="UInt16" compression="{compression}""#);
            assert!(read_block(&attrs, &block).is_err(), "{compression}");
        }
    }

    #[test]
    fn lz4_spec_block_decodes() {
        for codec in ["lz4", "lz4hc"] {
            let attrs = format!(r#"geometry="8:4" sampleFormat="UInt8" compression="{codec}:32""#);
            assert_eq!(
                read_block(&attrs, &LZ4_REPEAT_BLOCK).unwrap(),
                repeat_1234(),
                "{codec}"
            );
        }
    }

    #[test]
    fn lz4_block_must_decode_to_declared_size() {
        // A complete, valid block that decodes to 31 bytes, one short of the
        // declared 32.
        let short_block = [0x4F, 1, 2, 3, 4, 0x04, 0x00, 0x04, 0x40, 0x04, 1, 2, 3];
        let attrs = r#"geometry="8:4" sampleFormat="UInt8" compression="lz4:32""#;
        let err = read_block(attrs, &short_block).unwrap_err().to_string();
        assert!(
            err.contains("decompresses to 31 bytes instead of 32"),
            "{err}"
        );
    }

    #[test]
    fn byte_order_big_is_swapped_to_little() {
        let be: Vec<u8> = [1u16, 2, 3, 4]
            .iter()
            .flat_map(|v| v.to_be_bytes())
            .collect();
        let le: Vec<u8> = [1u16, 2, 3, 4]
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        let attrs = r#"geometry="2:2:1" sampleFormat="UInt16" byteOrder="big""#;
        assert_eq!(read_block(attrs, &be).unwrap(), le);
    }

    #[test]
    fn byte_order_big_applies_after_unshuffle() {
        let values = [1.5f64, -2.25, 1e300, f64::MIN_POSITIVE];
        let be: Vec<u8> = values.iter().flat_map(|v| v.to_be_bytes()).collect();
        let le: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        let block = zlib(&spec_shuffle(&be, 8));
        let attrs =
            r#"geometry="2:2" sampleFormat="Float64" byteOrder="big" compression="zlib+sh:32:8""#;
        assert_eq!(read_block(attrs, &block).unwrap(), le);
    }

    #[test]
    fn byte_order_little_and_unknown_values() {
        let le = vec![1u8, 0, 2, 0];
        let little = r#"geometry="2:1" sampleFormat="UInt16" byteOrder="little""#;
        assert_eq!(read_block(little, &le).unwrap(), le);
        let other = r#"geometry="2:1" sampleFormat="UInt16" byteOrder="middle""#;
        assert!(read_block(other, &le).is_err());
    }

    #[test]
    fn unknown_codec_is_an_error() {
        let data = vec![9u8, 9, 9, 9];
        for codec in ["brotli:4", "zstd2:4", "lz4-hc:4", "ZLIB:4", "none:4", ":4"] {
            let attrs = format!(r#"geometry="2:2:1" sampleFormat="UInt8" compression="{codec}""#);
            let err = read_block(&attrs, &data).unwrap_err();
            assert!(
                matches!(err, XisfError::UnsupportedCodec(_)),
                "{codec}: {err}"
            );
        }
    }

    #[test]
    fn zstd_is_reported_as_unsupported() {
        let attrs = r#"geometry="2:2:1" sampleFormat="UInt8" compression="zstd:4""#;
        let err = read_block(attrs, &[1, 2, 3]).unwrap_err();
        assert!(matches!(err, XisfError::UnsupportedCodec(_)), "{err}");
    }

    #[test]
    fn inflated_length_must_match_declared_size() {
        let attrs = r#"geometry="2:2:1" sampleFormat="UInt8" compression="zlib:4""#;
        assert!(read_block(attrs, &zlib(&[7, 7])).is_err());
    }

    #[test]
    fn declared_uncompressed_size_must_match_geometry() {
        let attrs = r#"geometry="2:2:1" sampleFormat="UInt8" compression="zlib:2""#;
        assert!(read_block(attrs, &zlib(&[7, 7])).is_err());
        let missing = r#"geometry="2:2:1" sampleFormat="UInt8" compression="zlib""#;
        assert!(read_block(missing, &zlib(&[7, 7, 7, 7])).is_err());
    }

    #[test]
    fn stored_length_must_match_geometry() {
        let mono = r#"geometry="2:2:1" sampleFormat="UInt8""#;
        assert!(read_block(mono, &[1, 2, 3]).is_err());
        assert!(read_block(mono, &[1, 2, 3, 4, 5]).is_err());
        let normal = r#"geometry="2:2:3" sampleFormat="UInt8" pixelStorage="Normal""#;
        assert!(read_block(normal, &[1, 2, 3]).is_err());
        let planar = r#"geometry="2:2:3" sampleFormat="UInt16""#;
        assert!(read_block(planar, &[0; 23]).is_err());
    }

    #[test]
    fn subblocks_are_decoded_and_concatenated() {
        let raw = distinct_bytes(16);
        let (first, second) = raw.split_at(10);
        let mut block = zlib(first);
        let first_len = block.len();
        block.extend_from_slice(&zlib(second));
        let attrs = format!(
            r#"geometry="4:4" sampleFormat="UInt8" compression="zlib:16" subblocks="{first_len},10:{},6""#,
            block.len() - first_len
        );
        assert_eq!(read_block(&attrs, &block).unwrap(), raw);
    }

    #[test]
    fn uncompressed_subblock_is_copied_and_unshuffled() {
        let raw = distinct_bytes(16);
        let shuffled = spec_shuffle(&raw, 2);
        let (first, second) = shuffled.split_at(10);
        let mut block = zlib(first);
        let first_len = block.len();
        block.extend_from_slice(second);
        let attrs = format!(
            r#"geometry="4:2" sampleFormat="UInt16" compression="zlib+sh:16:2" subblocks="{first_len},10:6,6""#
        );
        assert_eq!(read_block(&attrs, &block).unwrap(), raw);
    }

    #[test]
    fn whole_block_with_equal_sizes_is_stored_uncompressed() {
        let raw = distinct_bytes(8);
        let attrs = r#"geometry="2:2" sampleFormat="UInt16" compression="lz4+sh:8:2""#;
        assert_eq!(read_block(attrs, &spec_shuffle(&raw, 2)).unwrap(), raw);
    }

    #[test]
    fn subblock_sizes_are_validated() {
        let raw = distinct_bytes(16);
        let block = zlib(&raw);
        let c = block.len();
        for subblocks in [
            format!("{},16", c + 1),
            format!("{},15", c),
            format!("{},8:0,8", c),
            format!("{},17", c),
            format!("{c}"),
            format!("{c},16:"),
            "18446744073709551615,16:2,0".to_string(),
            String::new(),
        ] {
            let attrs = format!(
                r#"geometry="4:4" sampleFormat="UInt8" compression="zlib:16" subblocks="{subblocks}""#
            );
            assert!(read_block(&attrs, &block).is_err(), "{subblocks}");
        }
    }

    #[test]
    fn subblock_must_decompress_to_its_declared_size() {
        let raw = distinct_bytes(16);
        let mut block = zlib(&raw[..8]);
        let first_len = block.len();
        block.extend_from_slice(&zlib(&raw[8..]));
        let attrs = format!(
            r#"geometry="4:4" sampleFormat="UInt8" compression="zlib:16" subblocks="{first_len},9:{},7""#,
            block.len() - first_len
        );
        assert!(read_block(&attrs, &block).is_err());
    }

    #[test]
    fn huge_geometry_is_rejected_without_overflow() {
        let attrs =
            r#"geometry="4294967296:4294967296:3" sampleFormat="UInt16" pixelStorage="Normal""#;
        assert!(read_block(attrs, &[1; 6]).is_err());
    }

    #[test]
    fn location_past_end_of_file_is_rejected_before_allocating() {
        let attrs = r#"geometry="2:2:1" sampleFormat="UInt8""#;
        for size in ["18446744073709551615", "1099511627776", "5"] {
            let bytes = xisf_with_block_size(attrs, &[1; 4], size);
            let err = XisfFile::new(Cursor::new(bytes))
                .unwrap()
                .read_image_data_raw(0)
                .unwrap_err();
            assert!(err.to_string().contains("past the end"), "{size}: {err}");
        }
    }

    #[test]
    fn header_length_past_end_of_file_is_rejected_before_allocating() {
        let mut data = XISF_SIGNATURE.to_vec();
        data.extend_from_slice(&u32::MAX.to_le_bytes());
        data.extend_from_slice(&0u32.to_le_bytes());
        data.extend_from_slice(b"<xisf/>");
        let err = XisfFile::new(Cursor::new(data)).unwrap_err();
        assert!(err.to_string().contains("past the end"), "{err}");
    }

    #[test]
    fn lz4_declared_size_beyond_block_capacity_is_rejected() {
        let attrs =
            r#"geometry="65536:65536" sampleFormat="Float64" compression="lz4:34359738368""#;
        let err = read_block(attrs, &LZ4_REPEAT_BLOCK).unwrap_err();
        assert!(err.to_string().contains("cannot hold"), "{err}");
    }

    #[test]
    fn read_image_data_typed_returns_file_samples() {
        let samples = [1.5f32, 2.5, 3.5, 4.5];
        let bytes: Vec<u8> = samples.iter().flat_map(|v| v.to_le_bytes()).collect();
        let xisf = xisf_with_block(r#"geometry="2:2" sampleFormat="Float32""#, &bytes);
        let mut file = XisfFile::new(Cursor::new(xisf)).unwrap();
        let read: Vec<f32> = file.read_image_data_typed(0).unwrap();
        assert_eq!(read, samples);
    }

    #[test]
    fn unsupported_codec_leaves_the_rest_of_the_file_readable() {
        let xml = r#"<?xml version="1.0"?><xisf version="1.0">
<Image geometry="2:2" sampleFormat="UInt8" compression="brotli:4" location="attachment:1040:4" />
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:1044:4" />
<FITSKeyword name="OBJECT" value="'M31'" />
</xisf>"#;
        let mut data = create_xisf_data_with_xml(xml);
        data[8..12].copy_from_slice(&1024u32.to_le_bytes());
        data.resize(1040, 0);
        data.extend_from_slice(&[9, 9, 9, 9, 1, 2, 3, 4]);

        let mut file = XisfFile::new(Cursor::new(data)).unwrap();
        assert!(matches!(
            file.read_image_data_raw(0),
            Err(XisfError::UnsupportedCodec(_))
        ));
        assert_eq!(file.read_image_data_raw(1).unwrap(), vec![1, 2, 3, 4]);
        assert_eq!(file.get("OBJECT").as_str(), Some("M31"));
    }

    #[test]
    fn get_keyword_returns_first_of_repeated_name_and_accepts_hierarch_alias() {
        let xml = r#"<?xml version="1.0"?><xisf version="1.0">
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:1024:4">
<FITSKeyword name="OBJECT" value="'FIRST'" />
<FITSKeyword name="OBJECT" value="'SECOND'" />
<FITSKeyword name="HIERARCH ESO DET DIT" value="1" />
<FITSKeyword name="HIERARCH ESO DET DIT" value="2" />
</Image>
</xisf>"#;
        let mut data = create_xisf_data_with_xml(xml);
        data[8..12].copy_from_slice(&1008u32.to_le_bytes());
        data.resize(1024, 0);
        data.extend_from_slice(&[1, 2, 3, 4]);

        let file = XisfFile::new(Cursor::new(data)).unwrap();
        assert_eq!(file.get("OBJECT").as_str(), Some("FIRST"));
        assert_eq!(file.get("HIERARCH ESO DET DIT").as_i64(), Some(1));
        assert_eq!(file.get("ESO DET DIT").as_i64(), Some(1));
    }

    // XISF 1.0 section 11.6: a FITSKeyword child of an Image element belongs
    // to that image. One outside every Image element applies to every image.
    #[test]
    fn each_image_has_its_own_keywords() {
        let xml = r#"<?xml version="1.0"?><xisf version="1.0">
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:1024:4">
<FITSKeyword name="OBJECT" value="'A'" />
</Image>
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:1024:4">
<FITSKeyword name="OBJECT" value="'B'" />
</Image>
<FITSKeyword name="TELESCOP" value="'T'" />
</xisf>"#;
        let mut data = create_xisf_data_with_xml(xml);
        data.extend_from_slice(&[1, 2, 3, 4]);
        let file = XisfFile::new(Cursor::new(data)).unwrap();
        let telescop = Keyword::string("TELESCOP", "T");
        assert_eq!(
            file.keywords(),
            [Keyword::string("OBJECT", "A"), telescop.clone()]
        );
        assert_eq!(
            file.image_keywords(1),
            Some(&[Keyword::string("OBJECT", "B"), telescop][..])
        );
        assert_eq!(file.image_keywords(2), None);
    }

    fn keywords_of(xml: &str) -> XisfFile<Cursor<Vec<u8>>> {
        XisfFile::new(Cursor::new(create_xisf_data_with_xml(xml))).unwrap()
    }

    // A root keyword without a uid, or with a uid no image references,
    // applies to an image only when the image has no keyword of its name,
    // so the first and the last card of that name are the image's own.
    #[test]
    fn root_keywords_do_not_repeat_a_name_the_image_has() {
        let file = keywords_of(
            r#"<xisf version="1.0">
<FITSKeyword name="CRVAL1" value="20.0" />
<FITSKeyword name="EXPTIME" value="30" uid="UNUSED" />
<FITSKeyword name="HISTORY" value="" comment="root" />
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:0:4">
<FITSKeyword name="CRVAL1" value="10.0" />
<FITSKeyword name="HISTORY" value="" comment="image" />
</Image>
</xisf>"#,
        );
        assert_eq!(
            file.keywords(),
            [
                Keyword::real("CRVAL1", 10.0),
                Keyword::history("image"),
                Keyword::integer("EXPTIME", 30),
                Keyword::history("root"),
            ]
        );
    }

    // XISF 1.0 section 11.6: a Reference child of an Image element gives the
    // image the root FITSKeyword whose uid it names, at its position.
    #[test]
    fn a_referenced_root_keyword_belongs_only_to_the_referencing_images() {
        let file = keywords_of(
            r#"<xisf version="1.0">
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:0:4">
<FITSKeyword name="OBJECT" value="'A'" />
</Image>
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:0:4">
<Reference ref="KW1" />
<FITSKeyword name="OBJECT" value="'B'" />
</Image>
<FITSKeyword name="FILTER" value="'Ha'" uid="KW1" />
</xisf>"#,
        );
        assert_eq!(file.keywords(), [Keyword::string("OBJECT", "A")]);
        assert_eq!(
            file.image_keywords(1),
            Some(
                &[
                    Keyword::string("FILTER", "Ha"),
                    Keyword::string("OBJECT", "B")
                ][..]
            )
        );
    }

    // XISF 1.0 section 11.6 requires access to every FITSKeyword element.
    #[test]
    fn a_file_without_images_gives_the_root_keywords() {
        let file = keywords_of(
            r#"<xisf version="1.0">
<FITSKeyword name="OBJECT" value="'R'" />
<FITSKeyword name="TELESCOP" value="'T'" uid="K" />
</xisf>"#,
        );
        let root = [
            Keyword::string("OBJECT", "R"),
            Keyword::string("TELESCOP", "T"),
        ];
        assert_eq!(file.keywords(), root);
        assert_eq!(file.get("OBJECT").as_str(), Some("R"));
        assert_eq!(file.image_keywords(0), None);
    }

    // A FITSKeyword is a child of an Image element or of the root element
    // only (XISF 1.0 section 11.6).
    #[test]
    fn keywords_inside_metadata_or_thumbnail_are_skipped() {
        let file = keywords_of(
            r#"<xisf version="1.0">
<Metadata><FITSKeyword name="META" value="'m'" /></Metadata>
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:0:4">
<Thumbnail geometry="1:1:1" sampleFormat="UInt8" location="attachment:0:1">
<FITSKeyword name="THUMBKW" value="'t'" />
</Thumbnail>
<FITSKeyword name="OBJECT" value="'A'" />
</Image>
</xisf>"#,
        );
        assert_eq!(file.keywords(), [Keyword::string("OBJECT", "A")]);
    }

    // XISF 1.0 section 11.13: a Reference element names its uid wherever it
    // is, so FILTER is not a keyword for every image. Only a Reference child
    // of an Image element gives that image the keyword.
    #[test]
    fn a_reference_outside_an_image_marks_the_keyword_as_referenced() {
        let file = keywords_of(
            r#"<xisf version="1.0">
<FITSKeyword name="FILTER" value="'Ha'" uid="U" />
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:0:4">
<Thumbnail geometry="1:1:1" sampleFormat="UInt8" location="attachment:0:1">
<Reference ref="U" />
</Thumbnail>
<FITSKeyword name="OBJECT" value="'A'" />
</Image>
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:0:4" />
<Reference ref="U" />
</xisf>"#,
        );
        assert_eq!(file.keywords(), [Keyword::string("OBJECT", "A")]);
        assert_eq!(file.image_keywords(1), Some(&[][..]));
    }

    #[test]
    fn a_repeated_reference_attaches_the_keyword_once() {
        let file = keywords_of(
            r#"<xisf version="1.0">
<FITSKeyword name="FILTER" value="'Ha'" uid="U" />
<Image geometry="2:2" sampleFormat="UInt8" location="attachment:0:4">
<Reference ref="U" />
<Reference ref="U"></Reference>
</Image>
</xisf>"#,
        );
        assert_eq!(file.keywords(), [Keyword::string("FILTER", "Ha")]);
    }

    // A header of `padding` spaces plus `images` Image elements and `root`
    // unreferenced root keywords, so each image gets `root` copies.
    fn copies_header(images: usize, root: usize, length: usize) -> String {
        let mut xml = String::from(r#"<xisf version="1.0">"#);
        for i in 0..root {
            xml += &format!(r#"<FITSKeyword name="K{i}" value="1" />"#);
        }
        for _ in 0..images {
            xml += r#"<Image geometry="2:2" sampleFormat="UInt8" location="attachment:0:4" />"#;
        }
        let end = "</xisf>";
        assert!(xml.len() + end.len() <= length);
        xml += &" ".repeat(length - xml.len() - end.len());
        xml + end
    }

    // Each image holds a copy of each root keyword without a uid, so a
    // small header could make the reader copy keywords without bound. The
    // reader refuses a header that asks for more copies than it has bytes.
    #[test]
    fn root_keyword_copies_are_bounded_by_the_header_length() {
        let file = keywords_of(&copies_header(200, 200, 40_000));
        assert_eq!(file.image_keywords(199).map(<[Keyword]>::len), Some(200));

        let data = create_xisf_data_with_xml(&copies_header(200, 200, 39_999));
        let err = XisfFile::new(Cursor::new(data)).unwrap_err();
        assert!(matches!(err, XisfError::InvalidFormat(_)), "{err}");
        assert!(err.to_string().contains("200 images"), "{err}");
    }

    #[test]
    fn read_image_data_typed_rejects_other_sample_formats() {
        let xisf = create_xisf_with_u8_data(2, 2, &[1, 2, 3, 4]);
        let mut file = XisfFile::new(Cursor::new(xisf)).unwrap();
        let err = file.read_image_data_typed::<f32>(0).unwrap_err();
        assert!(
            matches!(err, XisfError::SampleFormatMismatch { .. }),
            "{err}"
        );
    }
}
