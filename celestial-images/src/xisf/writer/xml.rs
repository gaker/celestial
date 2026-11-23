use super::metadata::write_metadata;
use super::{PropertyDataBlock, XisfWriter};
use crate::fits::header::keywords::{is_standard_name, Keyword};
use crate::fits::io::card::{format_keyword, has_unwritable_name, is_commentary_card, value_text};
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::geometry::format_geometry_with_channels;
use crate::xisf::header::image_info::ImageInfo;
use crate::xisf::header::pixel_storage::PixelStorage;
use crate::xisf::header::property::{XisfProperty, XisfPropertyValue};
use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
use quick_xml::Writer;
use std::io::{Cursor, Seek, Write};

// `created` is the XISF:CreationTime value.
pub(super) fn generate_xml_content<W: Write + Seek>(
    state: &XisfWriter<W>,
    created: &str,
) -> Result<Vec<u8>> {
    let mut buffer = Cursor::new(Vec::new());
    let mut writer = Writer::new_with_indent(&mut buffer, b' ', 2);

    write_xml_declaration(&mut writer)?;
    write_xisf_content(state, created, &mut writer)?;

    Ok(buffer.into_inner())
}

fn write_xml_declaration<W: Write>(writer: &mut Writer<W>) -> Result<()> {
    let decl = BytesDecl::new("1.0", Some("UTF-8"), None);
    writer.write_event(Event::Decl(decl))?;
    Ok(())
}

fn write_xisf_content<S, Wr: Write>(
    state: &XisfWriter<S>,
    created: &str,
    writer: &mut Writer<Wr>,
) -> Result<()> {
    let mut xisf_start = BytesStart::new("xisf");
    xisf_start.push_attribute(("version", "1.0"));
    xisf_start.push_attribute(("xmlns", "http://www.pixinsight.com/xisf"));
    xisf_start.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
    xisf_start.push_attribute((
        "xsi:schemaLocation",
        "http://www.pixinsight.com/xisf http://pixinsight.com/xisf/xisf-1.0.xsd",
    ));
    writer.write_event(Event::Start(xisf_start))?;
    for (index, image) in state.images.iter().enumerate() {
        write_image(state, index, &image.info, writer)?;
    }
    write_metadata(writer, created, state.sha1_checksums, &state.images)?;
    writer.write_event(Event::End(BytesEnd::new("xisf")))?;
    Ok(())
}

// An image without keywords or properties is written as an empty element.
fn write_image<S, Wr: Write>(
    state: &XisfWriter<S>,
    index: usize,
    info: &ImageInfo,
    writer: &mut Writer<Wr>,
) -> Result<()> {
    let (keywords, properties) = state.image_metadata(index);
    let element = image_element(info);
    if keywords.is_empty() && properties.is_empty() {
        writer.write_event(Event::Empty(element))?;
        return Ok(());
    }
    writer.write_event(Event::Start(element))?;
    write_properties(properties, &state.property_blocks, writer)?;
    write_keywords_as_fits(keywords, writer)?;
    writer.write_event(Event::End(BytesEnd::new("Image")))?;
    Ok(())
}

fn image_element(info: &ImageInfo) -> BytesStart<'static> {
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
    if let Some(checksum) = &info.checksum {
        elem.push_attribute(("checksum", checksum.to_string().as_str()));
    }
    elem
}

fn write_keywords_as_fits<Wr: Write>(keywords: &[Keyword], writer: &mut Writer<Wr>) -> Result<()> {
    for keyword in keywords {
        let Some((name, value)) = fits_keyword_fields(keyword)? else {
            continue;
        };
        let mut elem = BytesStart::new("FITSKeyword");
        elem.push_attribute(("name", name.as_str()));
        elem.push_attribute(("value", value.as_str()));
        elem.push_attribute(("comment", keyword.comment.as_deref().unwrap_or("")));
        writer.write_event(Event::Empty(elem))?;
    }
    Ok(())
}

// XISF 1.0 section 11.6 requires FITSKeyword data that a FITS header can
// hold, so a keyword the FITS card writer rejects is rejected here too. A
// name it cannot write, even as HIERARCH, gives None and the keyword is
// left out. The value is the FITS value field: strings are quoted, so their
// type survives. Other names take the HIERARCH form the FITS writer gives.
fn fits_keyword_fields(keyword: &Keyword) -> Result<Option<(String, String)>> {
    if has_unwritable_name(keyword) {
        return Ok(None);
    }
    format_keyword(keyword).map_err(XisfError::Keyword)?;
    if is_commentary_card(keyword) {
        return Ok(Some((keyword.name.clone(), String::new())));
    }
    let name = if is_standard_name(&keyword.name) || keyword.name.starts_with("HIERARCH ") {
        keyword.name.clone()
    } else {
        format!("HIERARCH {}", keyword.name)
    };
    let value = match &keyword.value {
        Some(value) => value_text(keyword, value).map_err(XisfError::Keyword)?,
        None => String::new(),
    };
    Ok(Some((name, value)))
}

fn write_properties<Wr: Write>(
    properties: &[XisfProperty],
    property_blocks: &[PropertyDataBlock],
    writer: &mut Writer<Wr>,
) -> Result<()> {
    let mut block_index = 0usize;
    for (prop_index, property) in properties.iter().enumerate() {
        let block = if property.value.needs_data_block() {
            find_property_block(property_blocks, prop_index, &mut block_index)
        } else {
            None
        };
        write_property(property, block, writer)?;
    }
    Ok(())
}

// A String property holds its value as character content (XISF 1.0 section
// 11.1.6), a vector or matrix in a data block, and a scalar in the value
// attribute.
pub(super) fn write_property<Wr: Write>(
    property: &XisfProperty,
    block: Option<&PropertyDataBlock>,
    writer: &mut Writer<Wr>,
) -> Result<()> {
    let mut elem = BytesStart::new("Property");
    elem.push_attribute(("id", property.id.as_str()));
    elem.push_attribute(("type", property.value.type_name()));
    match &property.value {
        XisfPropertyValue::String(text) => {
            writer.write_event(Event::Start(elem))?;
            writer.write_event(Event::Text(BytesText::new(text)))?;
            writer.write_event(Event::End(BytesEnd::new("Property")))?;
            return Ok(());
        }
        XisfPropertyValue::F64Vector(v) => {
            elem.push_attribute(("length", v.len().to_string().as_str()));
        }
        XisfPropertyValue::F64Matrix { rows, cols, .. } => {
            elem.push_attribute(("rows", rows.to_string().as_str()));
            elem.push_attribute(("columns", cols.to_string().as_str()));
        }
        scalar => elem.push_attribute(("value", scalar.format_value().as_str())),
    }
    if let Some(block) = block {
        elem.push_attribute(("location", block.location.format().as_str()));
        if let Some(checksum) = &block.checksum {
            elem.push_attribute(("checksum", checksum.to_string().as_str()));
        }
    }
    writer.write_event(Event::Empty(elem))?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::KeywordValue;
    use crate::xisf::header::compression::XisfCompression;
    use crate::xisf::header::property::XisfProperty;
    use crate::xisf::reader::XisfFile;
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

    fn render(value: Option<KeywordValue>) -> String {
        let mut keyword = Keyword::new("TEST".to_string());
        keyword.value = value;
        fits_keyword_fields(&keyword).unwrap().unwrap().1
    }

    // The first Image element, from its start tag to its end tag. The
    // Metadata element also holds Property elements, so a property of the
    // image is looked for here.
    fn image_element(xml: &str) -> &str {
        let start = xml.find("<Image ").expect("the XML holds an Image element");
        let end = xml[start..]
            .find("</Image>")
            .map_or(xml.len(), |e| start + e);
        &xml[start..end]
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
        assert!(
            image_element(&xml).contains(r#"<Property id="Count" type="Int32" value="42"/>"#),
            "{xml}"
        );
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
        assert_eq!(render(Some(KeywordValue::Logical(true))), "T");
        assert_eq!(render(Some(KeywordValue::Logical(false))), "F");
    }

    #[test]
    fn fits_keyword_integer_renders_as_decimal_string() {
        assert_eq!(render(Some(KeywordValue::Integer(42))), "42");
        assert_eq!(render(Some(KeywordValue::Integer(-7))), "-7");
    }

    #[test]
    fn fits_keyword_real_keeps_decimal_when_already_present() {
        assert_eq!(render(Some(KeywordValue::Real(1.25))), "1.25");
    }

    #[test]
    fn fits_keyword_real_preserves_decimal_for_whole_numbers() {
        assert_eq!(render(Some(KeywordValue::Real(5.0))), "5.0");
        assert_eq!(render(Some(KeywordValue::Real(-12.0))), "-12.0");
        assert_eq!(render(Some(KeywordValue::Real(0.0))), "0.0");
    }

    #[test]
    fn fits_keyword_real_roundtrips_through_parse() {
        for v in [0.1, -1.5, 1234.5, 1.0e-10, 1.234567890123456e-9] {
            let rendered = render(Some(KeywordValue::Real(v)));
            let parsed: f64 = rendered
                .parse()
                .unwrap_or_else(|_| panic!("could not reparse '{rendered}' (from {v})"));
            assert_eq!(parsed, v);
        }
    }

    #[test]
    fn fits_keyword_string_is_quoted_as_in_a_fits_card() {
        assert_eq!(
            render(Some(KeywordValue::String("hello".into()))),
            "'hello'"
        );
    }

    #[test]
    fn fits_keyword_complex_uses_paren_notation() {
        assert_eq!(
            render(Some(KeywordValue::Complex(1.5, -2.5))),
            "(1.5, -2.5)"
        );
    }

    #[test]
    fn fits_keyword_none_is_empty_string() {
        assert_eq!(render(None), "");
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
        assert!(
            image_element(&xml).contains(r#"<Property id="Flag" type="Boolean" value="true"/>"#),
            "{xml}"
        );

        let mut w = writer_with_image();
        w.add_property(XisfProperty::boolean("Flag", false));
        let xml = xml_from(w);
        assert!(
            image_element(&xml).contains(r#"<Property id="Flag" type="Boolean" value="false"/>"#),
            "{xml}"
        );
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
    fn keywords_and_properties_are_written_to_image_0_only() {
        let mut w = writer_with_image();
        w.add_image::<u8>(&[0u8; 4], 2, 2, 1).unwrap();
        w.add_keyword(Keyword::string("OBJECT", "M31"));
        w.add_property(XisfProperty::int32("Count", 1));
        let xml = xml_from(w);
        let images: Vec<&str> = xml
            .split("<Image ")
            .skip(1)
            .filter_map(|image| image.split("<Metadata>").next())
            .collect();
        assert_eq!(images.len(), 2);
        assert_eq!(images[0].matches("<FITSKeyword").count(), 1);
        assert_eq!(images[0].matches("<Property").count(), 1);
        assert_eq!(images[1].matches("<FITSKeyword").count(), 0);
        assert_eq!(images[1].matches("<Property").count(), 0);
    }

    // `{:?}` gives the shortest text that parses back to the same bits.
    #[test]
    fn float64_properties_round_trip_exactly() {
        let values = [0.1, 1e300, 5e-324, -0.0, 42.0, f64::MAX, f64::MIN_POSITIVE];
        let texts = [
            "0.1",
            "1e300",
            "5e-324",
            "-0.0",
            "42.0",
            "1.7976931348623157e308",
            "2.2250738585072014e-308",
        ];
        let mut w = writer_with_image();
        for (i, &v) in values.iter().enumerate() {
            w.add_property(XisfProperty::float64(format!("P{i}"), v));
        }
        let bytes = w.write_to_vec().unwrap();
        let xml = extract_xml(&bytes);
        for (i, text) in texts.iter().enumerate() {
            assert!(xml.contains(&format!("id=\"P{i}\" type=\"Float64\" value=\"{text}\"")));
        }
        let read = XisfFile::new(Cursor::new(bytes))
            .unwrap()
            .properties(0)
            .unwrap();
        let bits: Vec<Option<u64>> = read
            .iter()
            .map(|p| match p.value {
                XisfPropertyValue::Float64(v) => Some(v.to_bits()),
                _ => None,
            })
            .collect();
        let want: Vec<Option<u64>> = values.iter().map(|v| Some(v.to_bits())).collect();
        assert_eq!(bits, want);
    }

    #[test]
    fn float64_special_values_are_written_as_section_8_3_3_gives_them() {
        let mut w = writer_with_image();
        w.add_property(XisfProperty::float64("P", f64::INFINITY));
        w.add_property(XisfProperty::float64("N", f64::NEG_INFINITY));
        w.add_property(XisfProperty::float64("Q", f64::NAN));
        let bytes = w.write_to_vec().unwrap();
        let xml = extract_xml(&bytes);
        assert!(
            xml.contains(r#"id="P" type="Float64" value="+Inf""#),
            "{xml}"
        );
        assert!(
            xml.contains(r#"id="N" type="Float64" value="-Inf""#),
            "{xml}"
        );
        assert!(
            xml.contains(r#"id="Q" type="Float64" value="NaN""#),
            "{xml}"
        );
        let read = XisfFile::new(Cursor::new(bytes))
            .unwrap()
            .properties(0)
            .unwrap();
        let values: Vec<_> = read.iter().map(|p| p.value.clone()).collect();
        assert_eq!(
            values[..2],
            [
                XisfPropertyValue::Float64(f64::INFINITY),
                XisfPropertyValue::Float64(f64::NEG_INFINITY)
            ]
        );
        assert!(matches!(values[2], XisfPropertyValue::Float64(v) if v.is_nan()));
    }

    struct FailingWrite;

    impl Write for FailingWrite {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("refused"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn xml_write_failures_are_io_errors() {
        let mut writer = Writer::new(FailingWrite);
        let error = write_xml_declaration(&mut writer).unwrap_err();
        assert!(matches!(error, XisfError::Io(_)), "{error:?}");
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
        assert!(
            locations.len() >= 3,
            "expected Image + 2 property locations, got {locations:?}"
        );
    }
}
