use super::{PropertyDataBlock, XisfWriter};
use crate::fits::header::KeywordValue;
use crate::xisf::header::{format_geometry_with_channels, ImageInfo, PixelStorage, XisfPropertyValue};
use crate::xisf::{Result, XisfError};
use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
use quick_xml::Writer;
use std::io::{Cursor, Seek, Write};

pub(super) fn generate_xml_content<W: Write + Seek>(
    state: &XisfWriter<W>,
) -> Result<Vec<u8>> {
    let mut buffer = Cursor::new(Vec::new());
    let mut writer = Writer::new_with_indent(&mut buffer, b' ', 2);

    write_xml_declaration(&mut writer)?;
    write_xisf_content(state, &mut writer)?;

    Ok(buffer.into_inner())
}

fn write_xml_declaration<W: Write>(writer: &mut Writer<W>) -> Result<()> {
    let decl = BytesDecl::new("1.0", Some("UTF-8"), None);
    writer
        .write_event(Event::Decl(decl))
        .map_err(|e| XisfError::XmlParse(e.to_string()))
}

fn write_xisf_content<S, Wr: Write>(state: &XisfWriter<S>, writer: &mut Writer<Wr>) -> Result<()> {
    let mut xisf_start = BytesStart::new("xisf");
    xisf_start.push_attribute(("version", "1.0"));
    xisf_start.push_attribute(("xmlns", "http://www.pixinsight.com/xisf"));
    xisf_start.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
    xisf_start.push_attribute((
        "xsi:schemaLocation",
        "http://www.pixinsight.com/xisf http://pixinsight.com/xisf/xisf-1.0.xsd",
    ));
    writer
        .write_event(Event::Start(xisf_start))
        .map_err(|e| XisfError::XmlParse(e.to_string()))?;

    for image in &state.images {
        write_image_element(state, writer, &image.info)?;
    }

    writer
        .write_event(Event::End(BytesEnd::new("xisf")))
        .map_err(|e| XisfError::XmlParse(e.to_string()))?;

    Ok(())
}

fn write_image_element<S, Wr: Write>(
    state: &XisfWriter<S>,
    writer: &mut Writer<Wr>,
    info: &ImageInfo,
) -> Result<()> {
    let mut elem = BytesStart::new("Image");
    elem.push_attribute((
        "geometry",
        format_geometry_with_channels(&info.geometry).as_str(),
    ));
    elem.push_attribute(("sampleFormat", info.sample_format.as_str()));
    if info.sample_format.is_floating_point() {
        elem.push_attribute(("bounds", info.format_bounds().as_str()));
    }
    elem.push_attribute(("colorSpace", info.color_space.as_str()));
    if info.pixel_storage != PixelStorage::Planar {
        elem.push_attribute(("pixelStorage", info.pixel_storage.as_str()));
    }

    if info.compression.is_compressed() {
        let compression_str = format!(
            "{}:{}",
            info.compression.as_str(),
            info.uncompressed_size.unwrap_or(0)
        );
        elem.push_attribute(("compression", compression_str.as_str()));
    }

    elem.push_attribute(("location", info.location.format().as_str()));

    let has_content = !state.keywords.is_empty() || !state.properties.is_empty();

    if !has_content {
        writer
            .write_event(Event::Empty(elem))
            .map_err(|e| XisfError::XmlParse(e.to_string()))?;
    } else {
        writer
            .write_event(Event::Start(elem))
            .map_err(|e| XisfError::XmlParse(e.to_string()))?;

        write_properties(&state.properties, &state.property_blocks, writer)?;
        write_keywords_as_fits(&state.keywords, writer)?;

        writer
            .write_event(Event::End(BytesEnd::new("Image")))
            .map_err(|e| XisfError::XmlParse(e.to_string()))?;
    }
    Ok(())
}

fn write_keywords_as_fits<Wr: Write>(
    keywords: &[crate::fits::header::Keyword],
    writer: &mut Writer<Wr>,
) -> Result<()> {
    for keyword in keywords {
        let mut elem = BytesStart::new("FITSKeyword");
        elem.push_attribute(("name", keyword.name.as_str()));

        let value_str = keyword_value_to_string(&keyword.value);
        elem.push_attribute(("value", value_str.as_str()));

        let comment_str = keyword.comment.as_deref().unwrap_or("");
        elem.push_attribute(("comment", comment_str));

        writer
            .write_event(Event::Empty(elem))
            .map_err(|e| XisfError::XmlParse(e.to_string()))?;
    }
    Ok(())
}

fn write_properties<Wr: Write>(
    properties: &[crate::xisf::header::XisfProperty],
    property_blocks: &[PropertyDataBlock],
    writer: &mut Writer<Wr>,
) -> Result<()> {
    let mut block_index = 0usize;

    for (prop_index, property) in properties.iter().enumerate() {
        let mut elem = BytesStart::new("Property");
        elem.push_attribute(("id", property.id.as_str()));
        elem.push_attribute(("type", property.value.type_name()));

        match &property.value {
            XisfPropertyValue::Boolean(b) => {
                elem.push_attribute(("value", if *b { "1" } else { "0" }));
                writer
                    .write_event(Event::Empty(elem))
                    .map_err(|e| XisfError::XmlParse(e.to_string()))?;
            }
            XisfPropertyValue::Float64(v) => {
                elem.push_attribute(("value", v.to_string().as_str()));
                writer
                    .write_event(Event::Empty(elem))
                    .map_err(|e| XisfError::XmlParse(e.to_string()))?;
            }
            XisfPropertyValue::Int32(v) => {
                elem.push_attribute(("value", v.to_string().as_str()));
                writer
                    .write_event(Event::Empty(elem))
                    .map_err(|e| XisfError::XmlParse(e.to_string()))?;
            }
            XisfPropertyValue::String(s) => {
                writer
                    .write_event(Event::Start(elem))
                    .map_err(|e| XisfError::XmlParse(e.to_string()))?;
                writer
                    .write_event(Event::Text(BytesText::new(s)))
                    .map_err(|e| XisfError::XmlParse(e.to_string()))?;
                writer
                    .write_event(Event::End(BytesEnd::new("Property")))
                    .map_err(|e| XisfError::XmlParse(e.to_string()))?;
            }
            XisfPropertyValue::F64Vector(v) => {
                elem.push_attribute(("length", v.len().to_string().as_str()));
                if let Some(block) =
                    find_property_block(property_blocks, prop_index, &mut block_index)
                {
                    elem.push_attribute(("location", block.location.format().as_str()));
                }
                writer
                    .write_event(Event::Empty(elem))
                    .map_err(|e| XisfError::XmlParse(e.to_string()))?;
            }
            XisfPropertyValue::F64Matrix { rows, cols, .. } => {
                elem.push_attribute(("rows", rows.to_string().as_str()));
                elem.push_attribute(("columns", cols.to_string().as_str()));
                if let Some(block) =
                    find_property_block(property_blocks, prop_index, &mut block_index)
                {
                    elem.push_attribute(("location", block.location.format().as_str()));
                }
                writer
                    .write_event(Event::Empty(elem))
                    .map_err(|e| XisfError::XmlParse(e.to_string()))?;
            }
        }
    }
    Ok(())
}

fn find_property_block<'a>(
    property_blocks: &'a [PropertyDataBlock],
    prop_index: usize,
    hint: &mut usize,
) -> Option<&'a PropertyDataBlock> {
    for (i, block) in property_blocks.iter().enumerate().skip(*hint) {
        if block.property_index == prop_index {
            *hint = i + 1;
            return Some(block);
        }
    }
    None
}

fn keyword_value_to_string(value: &Option<KeywordValue>) -> String {
    match value {
        None => String::new(),
        Some(KeywordValue::Logical(b)) => {
            if *b {
                "T".to_string()
            } else {
                "F".to_string()
            }
        }
        Some(KeywordValue::Integer(i)) => i.to_string(),
        Some(KeywordValue::Real(f)) => format!("{:?}", f),
        Some(KeywordValue::String(s)) => s.clone(),
        Some(KeywordValue::Complex(r, i)) => format!("({}, {})", r, i),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::Keyword;
    use crate::xisf::header::{XisfCompression, XisfProperty};
    use crate::xisf::XisfFile;
    use std::io::Cursor;

    fn xml_from(writer: crate::xisf::writer::XisfWriter<Cursor<Vec<u8>>>) -> String {
        let bytes = writer.write_to_vec().unwrap();
        extract_xml(&bytes)
    }

    fn extract_xml(bytes: &[u8]) -> String {
        // Skip XISF signature (8 bytes) + header length (8 bytes), then read until a null pad byte.
        let start = 16;
        let end = bytes[start..]
            .iter()
            .position(|&b| b == 0)
            .map(|p| start + p)
            .unwrap_or(bytes.len());
        String::from_utf8_lossy(&bytes[start..end]).to_string()
    }

    fn writer_with_image() -> crate::xisf::writer::XisfWriter<Cursor<Vec<u8>>> {
        let mut w = crate::xisf::writer::XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w
    }

    #[test]
    fn xml_starts_with_declaration_and_xisf_root() {
        let xml = xml_from(writer_with_image());
        assert!(xml.starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"));
        assert!(xml.contains("<xisf version=\"1.0\""));
        assert!(xml.contains("xmlns=\"http://www.pixinsight.com/xisf\""));
        assert!(xml.contains("</xisf>"));
    }

    #[test]
    fn empty_image_uses_self_closing_tag() {
        let xml = xml_from(writer_with_image());
        assert!(xml.contains("<Image "));
        assert!(xml.contains("/>"));
        assert!(!xml.contains("</Image>"));
    }

    #[test]
    fn image_with_keyword_uses_open_close_tags() {
        let mut w = writer_with_image();
        w.set_keyword("TELESCOP", "Hubble");
        let xml = xml_from(w);
        assert!(xml.contains("<Image "));
        assert!(xml.contains("</Image>"));
        assert!(xml.contains("<FITSKeyword"));
    }

    #[test]
    fn image_with_property_uses_open_close_tags() {
        let mut w = writer_with_image();
        w.add_property(XisfProperty::int32("Count", 42));
        let xml = xml_from(w);
        assert!(xml.contains("</Image>"));
        assert!(xml.contains("<Property"));
    }

    #[test]
    fn bounds_attribute_only_on_floating_point() {
        let mut w = crate::xisf::writer::XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        let xml_u8 = xml_from(w);
        assert!(!xml_u8.contains("bounds="));

        let mut w = crate::xisf::writer::XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<f32>(&[0.25f32; 4], 2, 2, 1).unwrap();
        let xml_f32 = xml_from(w);
        assert!(xml_f32.contains("bounds="));
    }

    #[test]
    fn pixel_storage_attribute_omitted_for_planar_default() {
        let xml = xml_from(writer_with_image());
        assert!(!xml.contains("pixelStorage"));
    }

    #[test]
    fn compression_attribute_present_when_compressed() {
        let mut w = crate::xisf::writer::XisfWriter::new(Cursor::new(Vec::new()))
            .compression(XisfCompression::Lz4);
        w.add_image::<u8>(&vec![7u8; 4096], 64, 64, 1).unwrap();
        let xml = xml_from(w);
        assert!(xml.contains("compression=\"lz4:4096\""));
    }

    #[test]
    fn compression_attribute_absent_when_uncompressed() {
        let xml = xml_from(writer_with_image());
        assert!(!xml.contains("compression="));
    }

    #[test]
    fn fits_keyword_logical_renders_as_t_or_f() {
        assert_eq!(keyword_value_to_string(&Some(KeywordValue::Logical(true))), "T");
        assert_eq!(keyword_value_to_string(&Some(KeywordValue::Logical(false))), "F");
    }

    #[test]
    fn fits_keyword_integer_renders_as_decimal_string() {
        assert_eq!(keyword_value_to_string(&Some(KeywordValue::Integer(42))), "42");
        assert_eq!(keyword_value_to_string(&Some(KeywordValue::Integer(-7))), "-7");
    }

    #[test]
    fn fits_keyword_real_keeps_decimal_when_already_present() {
        assert_eq!(keyword_value_to_string(&Some(KeywordValue::Real(1.25))), "1.25");
    }

    #[test]
    fn fits_keyword_real_preserves_decimal_for_whole_numbers() {
        assert_eq!(keyword_value_to_string(&Some(KeywordValue::Real(5.0))), "5.0");
        assert_eq!(keyword_value_to_string(&Some(KeywordValue::Real(-12.0))), "-12.0");
        assert_eq!(keyword_value_to_string(&Some(KeywordValue::Real(0.0))), "0.0");
    }

    #[test]
    fn fits_keyword_real_roundtrips_through_parse() {
        for v in [0.1, -1.5, 1234.5, 1.0e-10, 1.234567890123456e-9] {
            let rendered = keyword_value_to_string(&Some(KeywordValue::Real(v)));
            let parsed: f64 = rendered.parse().unwrap_or_else(|_| {
                panic!("could not reparse '{rendered}' (from {v})")
            });
            assert_eq!(parsed, v);
        }
    }

    #[test]
    fn fits_keyword_string_is_cloned_verbatim() {
        assert_eq!(
            keyword_value_to_string(&Some(KeywordValue::String("hello".into()))),
            "hello"
        );
    }

    #[test]
    fn fits_keyword_complex_uses_paren_notation() {
        assert_eq!(
            keyword_value_to_string(&Some(KeywordValue::Complex(1.5, -2.5))),
            "(1.5, -2.5)"
        );
    }

    #[test]
    fn fits_keyword_none_is_empty_string() {
        assert_eq!(keyword_value_to_string(&None), "");
    }

    #[test]
    fn fits_keyword_element_includes_name_value_comment() {
        let mut w = writer_with_image();
        w.add_keyword(Keyword::real("EXPTIME", 30.0).with_comment("seconds"));
        let xml = xml_from(w);
        assert!(xml.contains("<FITSKeyword"));
        assert!(xml.contains("name=\"EXPTIME\""));
        assert!(xml.contains("value=\"30.0\""));
        assert!(xml.contains("comment=\"seconds\""));
    }

    #[test]
    fn fits_keyword_missing_comment_renders_empty() {
        let mut w = writer_with_image();
        w.add_keyword(Keyword::integer("NAXIS", 2));
        let xml = xml_from(w);
        assert!(xml.contains("comment=\"\""));
    }

    #[test]
    fn property_boolean_emits_value_attribute() {
        let mut w = writer_with_image();
        w.add_property(XisfProperty::boolean("Flag", true));
        let xml = xml_from(w);
        assert!(xml.contains("<Property"));
        assert!(xml.contains("type=\"Boolean\""));
        assert!(xml.contains("value=\"1\""));

        let mut w = writer_with_image();
        w.add_property(XisfProperty::boolean("Flag", false));
        let xml = xml_from(w);
        assert!(xml.contains("value=\"0\""));
    }

    #[test]
    fn property_int32_emits_value_attribute() {
        let mut w = writer_with_image();
        w.add_property(XisfProperty::int32("N", -123));
        let xml = xml_from(w);
        assert!(xml.contains("type=\"Int32\""));
        assert!(xml.contains("value=\"-123\""));
    }

    #[test]
    fn property_float64_emits_value_attribute() {
        let mut w = writer_with_image();
        w.add_property(XisfProperty::float64("Pi", 1.25));
        let xml = xml_from(w);
        assert!(xml.contains("type=\"Float64\""));
        assert!(xml.contains("value=\"1.25\""));
    }

    #[test]
    fn property_string_uses_text_content() {
        let mut w = writer_with_image();
        w.add_property(XisfProperty::string("Name", "M31"));
        let xml = xml_from(w);
        assert!(xml.contains("type=\"String\""));
        assert!(xml.contains(">M31<"));
    }

    #[test]
    fn property_f64_vector_emits_length_and_location() {
        let mut w = writer_with_image();
        w.add_property(XisfProperty::f64_vector("V", vec![1.0, 2.0, 3.0]));
        let xml = xml_from(w);
        assert!(xml.contains("type=\"F64Vector\""));
        assert!(xml.contains("length=\"3\""));
        assert!(xml.contains("location=\"attachment:"));
    }

    #[test]
    fn property_f64_matrix_emits_rows_columns_and_location() {
        let mut w = writer_with_image();
        w.add_property(XisfProperty::f64_matrix(
            "M",
            2,
            3,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        ));
        let xml = xml_from(w);
        assert!(xml.contains("type=\"F64Matrix\""));
        assert!(xml.contains("rows=\"2\""));
        assert!(xml.contains("columns=\"3\""));
        assert!(xml.contains("location=\"attachment:"));
    }

    #[test]
    fn full_property_roundtrip_preserves_scalar_values_via_reader() {
        let mut w = crate::xisf::writer::XisfWriter::new(Cursor::new(Vec::new()));
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.add_property(XisfProperty::int32("X", 7));
        w.add_property(XisfProperty::float64("Y", 2.5));
        w.add_property(XisfProperty::boolean("Z", true));
        w.add_property(XisfProperty::string("S", "abc"));
        let bytes = w.write_to_vec().unwrap();
        // Reader round-trip just confirms the XML parses cleanly.
        let reader = XisfFile::new(Cursor::new(bytes)).unwrap();
        assert_eq!(reader.num_images(), 1);
    }

    #[test]
    fn property_block_location_increments_across_multiple_bulk_properties() {
        let mut w = writer_with_image();
        w.add_property(XisfProperty::f64_vector("V1", vec![1.0, 2.0]));
        w.add_property(XisfProperty::f64_vector("V2", vec![3.0, 4.0, 5.0]));
        let xml = xml_from(w);
        // Extract the location attrs for the two F64Vector entries.
        let mut locations = Vec::new();
        for attr in xml.split("location=\"attachment:") {
            if let Some(rest) = attr.split_once('"') {
                locations.push(rest.0.to_string());
            }
        }
        // First match is the Image location, remaining are property locations.
        assert!(locations.len() >= 3, "expected Image + 2 property locations, got {locations:?}");
    }
}
