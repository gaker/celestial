use super::Image;
use crate::core::errors::{ImageError, Result};
use crate::core::shape::plane_shape_with_axes;
use crate::fits::header::keywords::Keyword;
use crate::fits::io::writer::is_checksum_keyword;
use crate::formats::atomic_write::write_atomically;
use crate::formats::pixel_data::PixelData;
use crate::xisf::errors::XisfError;
use crate::xisf::header::property::XisfProperty;
use crate::xisf::header::sample_format::SampleFormat;
use crate::xisf::reader::XisfFile;
use crate::xisf::solution::{celestial_wcs, solution_keywords, SOLUTION_IDS};
use crate::xisf::writer::data_type::XisfDataType;
use crate::xisf::writer::wcs::solution_properties;
use crate::xisf::writer::XisfWriter;
use std::io::Cursor;
use std::path::Path;

impl Image {
    pub(super) fn open_xisf(path: &Path) -> Result<Self> {
        let mut xisf = XisfFile::open(path).map_err(ImageError::Xisf)?;
        let info = xisf.image_info(0).ok_or_else(|| {
            ImageError::Xisf(XisfError::InvalidFormat("No images in file".to_string()))
        })?;
        let dimensions = info.geometry.clone();
        let sample_format = info.sample_format.clone();
        let stored_interleaved = info.stores_interleaved();

        let stored = xisf.read_stored_samples(0).map_err(ImageError::Xisf)?;
        let keywords = describing_keywords(xisf.keywords());
        let properties = xisf.properties(0).map_err(ImageError::Xisf)?;
        let (keywords, xisf_properties) = with_solution_keywords(keywords, properties);

        let pixels = pixel_data(&sample_format, stored);
        Ok(Self::from_parts(
            pixels,
            dimensions,
            keywords,
            xisf_properties,
            stored_interleaved,
        ))
    }

    pub(super) fn save_xisf(&self, path: &Path) -> Result<()> {
        let (width, height, channels) =
            plane_shape_with_axes(&self.dimensions, ImageError::InvalidDimensions)?;
        let mut writer = XisfWriter::new(Cursor::new(Vec::new()));

        if channels == 3 {
            let mut planar_image = self.clone();
            planar_image.interleaved_to_planar();
            write_xisf_pixels(&mut writer, &planar_image.pixels, width, height, channels)?;
        } else {
            write_xisf_pixels(&mut writer, &self.pixels, width, height, channels)?;
        }

        write_xisf_file(writer, &self.keywords, &self.xisf_properties, path)
            .map_err(ImageError::Xisf)
    }
}

// Adds image 0's metadata to `writer`, which must already hold the image:
// `keywords` without the FITS structural cards, then the astrometric
// solution properties computed from `keywords` merged over `properties`.
// Writes the XISF unit to memory, then replaces the file at `path` with it,
// so an error leaves an existing file unchanged.
pub(crate) fn write_xisf_file(
    mut writer: XisfWriter<Cursor<Vec<u8>>>,
    keywords: &[Keyword],
    properties: &[XisfProperty],
    path: &Path,
) -> crate::xisf::errors::Result<()> {
    for keyword in describing_keywords(keywords) {
        writer.add_keyword(keyword);
    }
    let generated = solution_properties(keywords);
    writer.add_properties(merged_properties(generated, properties));
    let bytes = writer.write_to_vec()?;
    Ok(write_atomically(path, &bytes)?)
}

fn pixel_data(format: &SampleFormat, bytes: Vec<u8>) -> PixelData {
    match format {
        SampleFormat::UInt8 => PixelData::U8(bytes),
        SampleFormat::UInt16 => PixelData::U16(XisfDataType::from_le_bytes(&bytes)),
        SampleFormat::UInt32 => PixelData::U32(XisfDataType::from_le_bytes(&bytes)),
        SampleFormat::UInt64 => PixelData::U64(XisfDataType::from_le_bytes(&bytes)),
        SampleFormat::Float32 => PixelData::F32(XisfDataType::from_le_bytes(&bytes)),
        SampleFormat::Float64 => PixelData::F64(XisfDataType::from_le_bytes(&bytes)),
    }
}

// The structural cards describe a FITS data unit, not this image: XISF gives
// the geometry and sample format in the Image element. Checksums cover FITS
// bytes that an XISF file does not hold.
fn describing_keywords(keywords: &[Keyword]) -> Vec<Keyword> {
    let structural = |k: &Keyword| {
        k.is_mandatory()
            || matches!(k.name.as_str(), "BZERO" | "BSCALE")
            || is_checksum_keyword(&k.name)
    };
    keywords
        .iter()
        .filter(|k| !structural(k))
        .cloned()
        .collect()
}

// Properties derived from the keywords replace stored ones with the same id,
// so the file holds one solution, and it agrees with the keywords.
fn merged_properties(generated: Vec<XisfProperty>, stored: &[XisfProperty]) -> Vec<XisfProperty> {
    let replaced = |p: &XisfProperty| generated.iter().any(|g| g.id == p.id);
    let kept: Vec<XisfProperty> = stored.iter().filter(|p| !replaced(p)).cloned().collect();
    generated.into_iter().chain(kept).collect()
}

// Keywords that build a celestial WCS are the solution, and the solution
// properties are dropped so the image does not hold two. Without such
// keywords a usable solution becomes keywords. Keywords whose WCS does not
// build leave the properties as they are.
fn with_solution_keywords(
    mut keywords: Vec<Keyword>,
    mut properties: Vec<XisfProperty>,
) -> (Vec<Keyword>, Vec<XisfProperty>) {
    let derived = match celestial_wcs(&keywords) {
        Ok(Some(_)) => Some(Vec::new()),
        Ok(None) => solution_keywords(&properties),
        Err(_) => None,
    };
    if let Some(derived) = derived {
        keywords.extend(derived);
        properties.retain(|p| !SOLUTION_IDS.contains(&p.id.as_str()));
    }
    (keywords, properties)
}

fn write_xisf_pixels<W: std::io::Write + std::io::Seek>(
    writer: &mut XisfWriter<W>,
    pixels: &PixelData,
    width: usize,
    height: usize,
    channels: usize,
) -> Result<()> {
    let written = match pixels {
        PixelData::U8(data) => writer.add_image(data, width, height, channels),
        PixelData::U16(data) => writer.add_image(data, width, height, channels),
        PixelData::U32(data) => writer.add_image(data, width, height, channels),
        PixelData::U64(data) => writer.add_image(data, width, height, channels),
        PixelData::F32(data) => writer.add_image(data, width, height, channels),
        PixelData::F64(data) => writer.add_image(data, width, height, channels),
        PixelData::I8(_) => Err(signed_samples("Int8", true)),
        PixelData::I16(_) => Err(signed_samples("Int16", true)),
        PixelData::I32(_) => Err(signed_samples("Int32", true)),
        PixelData::I64(_) => Err(signed_samples("Int64", false)),
    };
    written.map_err(ImageError::Xisf)
}

// XISF 1.0 section 11.5.1, Table 11, has no signed integer sample format.
fn signed_samples(format: &str, has_float_form: bool) -> XisfError {
    let float = if has_float_form {
        " or PixelData::signed_to_float"
    } else {
        ""
    };
    XisfError::UnsupportedFormat(format!(
        "{format} pixels: XISF 1.0 has no signed integer sample format; \
         convert them with PixelData::signed_to_offset_unsigned{float}"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::{Keyword, KeywordValue};
    use crate::xisf::header::property::{XisfProperty, XisfPropertyValue};

    fn tmp_xisf() -> tempfile::NamedTempFile {
        tempfile::Builder::new().suffix(".xisf").tempfile().unwrap()
    }

    // An image with no axes holds no samples, so its width is 0 and the
    // writer refuses it, the way FITS refuses data under NAXIS = 0.
    #[test]
    fn roundtrip_u8_mono() {
        let tmp = tmp_xisf();
        let original: Vec<u8> = (0u8..64).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![8usize, 8]);
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    // The writer gives every image at least two dimensions, so a FITS-like
    // NAXIS = 1 image [n] is saved and read back as [n, 1].
    #[test]
    fn one_dimensional_image_is_saved_as_two_dimensional() {
        let tmp = tmp_xisf();
        let original = vec![0.25f32, 0.5, 0.75, 1.0];
        let img = Image::new(PixelData::F32(original.clone()), vec![4usize]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![4usize, 1]);
        assert_eq!(restored.pixels.as_f32().unwrap(), &original);
        let file = crate::xisf::reader::XisfFile::open(tmp.path()).unwrap();
        assert_eq!(file.image_info(0).unwrap().geometry, vec![4usize, 1]);
    }

    #[test]
    fn roundtrip_u16_mono() {
        let tmp = tmp_xisf();
        let original: Vec<u16> = (0u16..64).collect();
        let img = Image::new(PixelData::U16(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u16().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f32_mono() {
        let tmp = tmp_xisf();
        let original: Vec<f32> = (0..64).map(|i| i as f32 * 0.1).collect();
        let img = Image::new(PixelData::F32(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f32().unwrap(), &original);
    }

    #[test]
    fn roundtrip_f64_mono() {
        let tmp = tmp_xisf();
        let original: Vec<f64> = (0..64).map(|i| i as f64 * 0.01).collect();
        let img = Image::new(PixelData::F64(original.clone()), vec![8usize, 8]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_f64().unwrap(), &original);
    }

    #[test]
    fn rgb_roundtrip_preserves_channel_count() {
        let tmp = tmp_xisf();
        let original: Vec<u8> = (0u8..48).collect();
        let img = Image::new(PixelData::U8(original.clone()), vec![4usize, 4, 3]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.dimensions, vec![4usize, 4, 3]);
        assert_eq!(restored.channels(), 3);
        // Pixel layout in memory must be interleaved (RGB RGB ...) so
        // display code can pass the buffer straight to a texture upload.
        assert_eq!(restored.pixels.as_u8().unwrap(), &original);
    }

    #[test]
    fn open_reads_uint32_samples_as_u32() {
        let tmp = tmp_xisf();
        let samples = vec![0u32, 1, 2147483648, 4294967295];
        let mut writer = XisfWriter::create(tmp.path()).unwrap();
        writer.add_image(&samples, 4, 1, 1).unwrap();
        writer.write().unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u32(), Some(&samples));
    }

    #[test]
    fn roundtrip_u32_mono() {
        let tmp = tmp_xisf();
        let original = vec![0u32, 7, 2147483648, 4294967295];
        let img = Image::new(PixelData::U32(original.clone()), vec![2usize, 2]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u32(), Some(&original));
    }

    #[test]
    fn roundtrip_u64_mono() {
        let tmp = tmp_xisf();
        let original = vec![0u64, (1 << 53) + 1, 1 << 63, u64::MAX];
        let img = Image::new(PixelData::U64(original.clone()), vec![2usize, 2]);
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(restored.pixels.as_u64(), Some(&original));
    }

    #[test]
    fn save_rejects_i8() {
        let tmp = tmp_xisf();
        let img = Image::new(PixelData::I8(vec![-1, 0, 1, 2]), vec![2usize, 2]);
        assert!(matches!(
            img.save(tmp.path()),
            Err(ImageError::Xisf(
                crate::xisf::errors::XisfError::UnsupportedFormat(_)
            ))
        ));
    }

    #[test]
    fn save_reports_write_failures_as_xisf_io_errors() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("missing").join("out.xisf");
        let img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        assert!(matches!(
            img.save(&path),
            Err(ImageError::Xisf(crate::xisf::errors::XisfError::Io(_)))
        ));
    }

    #[test]
    fn failed_save_leaves_existing_file_untouched() {
        let tmp = tmp_xisf();
        std::fs::write(tmp.path(), b"keep").unwrap();
        for pixels in [
            PixelData::I8(vec![-1, 0, 1, 2]),
            PixelData::I64(vec![-1, 0, 1, 2]),
        ] {
            let img = Image::new(pixels, vec![2usize, 2]);
            assert!(matches!(
                img.save(tmp.path()),
                Err(ImageError::Xisf(
                    crate::xisf::errors::XisfError::UnsupportedFormat(_)
                ))
            ));
            assert_eq!(std::fs::read(tmp.path()).unwrap(), b"keep");
        }
    }

    #[test]
    fn keywords_roundtrip_through_xisf() {
        let tmp = tmp_xisf();
        let mut img = Image::new(PixelData::U8(vec![0; 16]), vec![4usize, 4]);
        img.set_keyword(Keyword::string("OBJECT", "NGC1234"));
        img.set_keyword(Keyword::real("EXPTIME", 60.0));
        img.save(tmp.path()).unwrap();

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(
            restored.get_keyword("OBJECT").unwrap().value,
            Some(crate::fits::header::keywords::KeywordValue::String(
                "NGC1234".to_string()
            ))
        );
    }

    #[test]
    fn open_rejects_planar_rgb_shorter_than_geometry() {
        let tmp = tmp_xisf();
        std::fs::write(tmp.path(), short_planar_rgb()).unwrap();
        assert!(matches!(
            Image::open(tmp.path()),
            Err(ImageError::Xisf(XisfError::InvalidFormat(_)))
        ));
    }

    // A 2x2 RGB UInt8 image that stores 3 of its 12 bytes.
    fn short_planar_rgb() -> Vec<u8> {
        let xml = r#"<xisf version="1.0"><Image geometry="2:2:3" sampleFormat="UInt8" location="attachment:128:3" /></xisf>"#;
        let mut bytes = b"XISF0100".to_vec();
        bytes.extend_from_slice(&112u32.to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(xml.as_bytes());
        bytes.resize(128, 0);
        bytes.extend_from_slice(&[1, 2, 3]);
        bytes
    }

    // XISF 1.0 Table 11 has no signed sample format, so D11 refuses them.
    #[test]
    fn i16_is_refused_and_converts_to_unsigned_on_disk() {
        let tmp = tmp_xisf();
        let img = Image::new(PixelData::I16(vec![0, 100, -50, -1]), vec![2usize, 2]);
        assert!(matches!(
            img.save(tmp.path()),
            Err(ImageError::Xisf(XisfError::UnsupportedFormat(_)))
        ));

        let converted = img.pixels.signed_to_offset_unsigned().unwrap();
        Image::new(converted, vec![2usize, 2])
            .save(tmp.path())
            .unwrap();
        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(
            restored.pixels.as_u16(),
            Some(&vec![32768, 32868, 32718, 32767])
        );
    }

    #[test]
    fn i32_is_refused_and_converts_to_unsigned_on_disk() {
        let tmp = tmp_xisf();
        let img = Image::new(PixelData::I32(vec![0, 100, -50, -1]), vec![2usize, 2]);
        assert!(matches!(
            img.save(tmp.path()),
            Err(ImageError::Xisf(XisfError::UnsupportedFormat(_)))
        ));

        let converted = img.pixels.signed_to_offset_unsigned().unwrap();
        Image::new(converted, vec![2usize, 2])
            .save(tmp.path())
            .unwrap();
        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(
            restored.pixels.as_u32(),
            Some(&vec![2147483648, 2147483748, 2147483598, 2147483647])
        );
    }

    // The XML header of an XISF file, as text.
    fn header_xml(path: &Path) -> String {
        let bytes = std::fs::read(path).unwrap();
        let len = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize;
        String::from_utf8_lossy(&bytes[16..16 + len]).to_string()
    }

    const SOLUTION_IDS: [&str; 6] = [
        "PCL:AstrometricSolution:ProjectionSystem",
        "PCL:AstrometricSolution:ReferenceCelestialCoordinates",
        "PCL:AstrometricSolution:ReferenceImageCoordinates",
        "PCL:AstrometricSolution:LinearTransformationMatrix",
        "PCL:AstrometricSolution:ReferenceNativeCoordinates",
        "PCL:AstrometricSolution:CelestialPoleNativeCoordinates",
    ];

    // CRPIX (512, 512), CRVAL (37.5, 30), CD diag(-1e-3, 1e-3), with the
    // integer-valued cards a FITS reader types as Integer.
    fn integer_valued_tan_image() -> Image {
        let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        img.keywords = vec![
            Keyword::string("CTYPE1", "RA---TAN"),
            Keyword::string("CTYPE2", "DEC--TAN"),
            Keyword::integer("CRPIX1", 512),
            Keyword::integer("CRPIX2", 512),
            Keyword::real("CRVAL1", 37.5),
            Keyword::integer("CRVAL2", 30),
            Keyword::real("CD1_1", -1e-3),
            Keyword::integer("CD1_2", 0),
            Keyword::integer("CD2_1", 0),
            Keyword::real("CD2_2", 1e-3),
        ];
        img
    }

    #[test]
    fn integer_valued_wcs_keywords_keep_the_astrometric_solution() {
        let tmp = tmp_xisf();
        integer_valued_tan_image().save(tmp.path()).unwrap();
        let xml = header_xml(tmp.path());
        for id in SOLUTION_IDS {
            assert!(xml.contains(&format!("id=\"{id}\"")), "{id} missing");
        }
    }

    #[test]
    fn save_without_wcs_writes_no_astrometric_properties() {
        let tmp = tmp_xisf();
        Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2])
            .save(tmp.path())
            .unwrap();
        assert!(!header_xml(tmp.path()).contains("AstrometricSolution"));
    }

    #[test]
    fn save_rejects_i16_and_i32_before_creating_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("signed.xisf");
        for pixels in [
            PixelData::I16(vec![-1, 0, 1, 2]),
            PixelData::I32(vec![-1, 0, 1, 2]),
        ] {
            let err = Image::new(pixels, vec![2usize, 2]).save(&path).unwrap_err();
            let message = err.to_string();
            assert!(
                matches!(err, ImageError::Xisf(XisfError::UnsupportedFormat(_))),
                "{message}"
            );
            assert!(message.contains("PixelData::signed_to_offset_unsigned"));
            assert!(message.contains("PixelData::signed_to_float"));
            assert!(!path.exists());
        }
    }

    #[test]
    fn save_rejects_i8_and_i64_naming_the_conversion_helpers() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("signed.xisf");
        let i8_err = Image::new(PixelData::I8(vec![-1, 0, 1, 2]), vec![2usize, 2])
            .save(&path)
            .unwrap_err()
            .to_string();
        assert!(i8_err.contains("PixelData::signed_to_offset_unsigned"));
        assert!(i8_err.contains("PixelData::signed_to_float"));
        let i64_err = Image::new(PixelData::I64(vec![-1, 0, 1, 2]), vec![2usize, 2])
            .save(&path)
            .unwrap_err()
            .to_string();
        assert!(i64_err.contains("PixelData::signed_to_offset_unsigned"));
        assert!(!path.exists());
    }

    #[test]
    fn keyword_values_keep_their_fits_type_through_xisf() {
        let tmp = tmp_xisf();
        let keywords = vec![
            Keyword::string("FILTER", "T"),
            Keyword::string("OBJECT", "1234"),
            Keyword::string("OBSERVER", "O'Brien"),
            Keyword::string("NOTE", "a < b & c > \"d\" 'e'").with_comment("<&\"'>"),
            Keyword::string("LEADING", "  padded"),
            Keyword::logical("FLAG", true),
            Keyword::integer("GAIN", -120),
            Keyword::real("EXPTIME", 300.0),
            Keyword::real("XPIXSZ", 3.76),
            Keyword::real("TINY", 1.0e-300),
            Keyword::new("CPLX".to_string()).with_value(KeywordValue::Complex(1.5, -2.0)),
        ];
        let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        img.keywords = keywords.clone();
        img.save(tmp.path()).unwrap();

        assert_eq!(Image::open(tmp.path()).unwrap().keywords, keywords);
    }

    // The FITS reader keeps HIERARCH names under "HIERARCH <name>", commentary
    // cards and undefined values. XISF must hand back the same keywords.
    #[test]
    fn fits_keywords_keep_hierarch_commentary_and_undefined_values_through_xisf() {
        let dir = tempfile::tempdir().unwrap();
        let fits_path = dir.path().join("in.fits");
        let xisf_path = dir.path().join("out.xisf");
        let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        img.keywords = vec![
            Keyword::real("ESO DET DIT", 1.5).with_comment("detector"),
            Keyword::string("HIERARCH ESO OBS NAME", "field 1"),
            Keyword::comment("  indented comment"),
            Keyword::history("dark subtracted"),
            Keyword::commentary("", "blank keyword text"),
            Keyword::new("UNDEF".to_string()).with_comment("no value"),
        ];
        img.save(&fits_path).unwrap();
        let from_fits = Image::open(&fits_path).unwrap();
        from_fits.save(&xisf_path).unwrap();

        let expected: Vec<Keyword> = from_fits
            .keywords
            .iter()
            .filter(|k| !k.is_mandatory())
            .cloned()
            .collect();
        assert!(expected.iter().any(|k| k.name == "HIERARCH ESO DET DIT"));
        assert_eq!(Image::open(&xisf_path).unwrap().keywords, expected);
    }

    #[test]
    fn structural_fits_keywords_are_not_written_to_xisf() {
        let tmp = tmp_xisf();
        let mut img = Image::new(PixelData::U16(vec![0; 4]), vec![2usize, 2]);
        img.keywords = vec![
            Keyword::logical("SIMPLE", true),
            Keyword::integer("BITPIX", 16),
            Keyword::integer("NAXIS", 2),
            Keyword::integer("NAXIS1", 2),
            Keyword::integer("NAXIS2", 2),
            Keyword::logical("EXTEND", true),
            Keyword::integer("BZERO", 32768),
            Keyword::integer("BSCALE", 1),
            Keyword::string("CHECKSUM", "hcHjjc9ghcEghc9g"),
            Keyword::string("DATASUM", "0"),
            Keyword::string("OBJECT", "M31"),
        ];
        img.save(tmp.path()).unwrap();
        let xml = header_xml(tmp.path());
        assert_eq!(xml.matches("<FITSKeyword").count(), 1, "{xml}");
        assert!(xml.contains("name=\"OBJECT\""));
    }

    #[test]
    fn structural_fits_keywords_in_an_xisf_file_are_not_loaded() {
        let tmp = tmp_xisf();
        let mut writer = XisfWriter::create(tmp.path()).unwrap();
        writer.add_image(&[0u16; 4], 2, 2, 1).unwrap();
        writer.add_keyword(Keyword::integer("BITPIX", 16));
        writer.add_keyword(Keyword::integer("BZERO", 32768));
        writer.add_keyword(Keyword::string("OBJECT", "M31"));
        writer.write().unwrap();

        let names: Vec<String> = Image::open(tmp.path())
            .unwrap()
            .keywords
            .into_iter()
            .map(|k| k.name)
            .collect();
        assert_eq!(names, vec!["OBJECT".to_string()]);
    }

    #[test]
    fn image_properties_are_read_back() {
        let tmp = tmp_xisf();
        let properties = vec![
            XisfProperty::float64("Test:Float", -1.25e-7),
            XisfProperty::int32("Test:Int", -42),
            XisfProperty::boolean("Test:Flag", true),
            XisfProperty::string("Test:Text", "  a < b & 'c' \"d\"  "),
            XisfProperty::f64_vector("Test:Vector", vec![1.5, -2.5, 3.25]),
            XisfProperty::f64_matrix("Test:Matrix", 2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        ];
        let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        img.xisf_properties = properties.clone();
        img.save(tmp.path()).unwrap();

        assert_eq!(Image::open(tmp.path()).unwrap().xisf_properties, properties);
    }

    // D7: ReferenceImageCoordinates is CRPIX - 0.5, so the reader adds 0.5.
    #[test]
    fn solution_properties_alone_give_wcs_keywords_on_open() {
        let tmp = tmp_xisf();
        let mut writer = XisfWriter::create(tmp.path()).unwrap();
        writer.add_image(&[0u8; 4], 2, 2, 1).unwrap();
        writer.add_properties([
            XisfProperty::string(SOLUTION_IDS[0], "Gnomonic"),
            XisfProperty::f64_vector(SOLUTION_IDS[1], vec![37.5, 30.0]),
            XisfProperty::f64_vector(SOLUTION_IDS[2], vec![511.5, 255.25]),
            XisfProperty::f64_matrix(SOLUTION_IDS[3], 2, 2, vec![-1e-3, 2e-6, 3e-6, 1e-3]),
        ]);
        writer.write().unwrap();

        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.get("CTYPE1").as_str(), Some("RA---TAN"));
        assert_eq!(img.get("CTYPE2").as_str(), Some("DEC--TAN"));
        assert_eq!(img.get("CRPIX1").as_f64(), Some(512.0));
        assert_eq!(img.get("CRPIX2").as_f64(), Some(255.75));
        assert_eq!(img.get("CRVAL1").as_f64(), Some(37.5));
        assert_eq!(img.get("CRVAL2").as_f64(), Some(30.0));
        assert_eq!(img.get("CD1_1").as_f64(), Some(-1e-3));
        assert_eq!(img.get("CD1_2").as_f64(), Some(2e-6));
        assert_eq!(img.get("CD2_1").as_f64(), Some(3e-6));
        assert_eq!(img.get("CD2_2").as_f64(), Some(1e-3));
    }

    #[test]
    fn keyword_wcs_replaces_stored_solution_properties_on_save_and_open() {
        let tmp = tmp_xisf();
        let mut img = integer_valued_tan_image();
        img.xisf_properties = vec![
            XisfProperty::string(SOLUTION_IDS[0], "Mercator"),
            XisfProperty::int32("Test:Kept", 3),
        ];
        img.save(tmp.path()).unwrap();
        let xml = header_xml(tmp.path());
        assert_eq!(xml.matches("ProjectionSystem").count(), 1, "{xml}");
        assert!(xml.contains(">Gnomonic<"));

        let restored = Image::open(tmp.path()).unwrap();
        assert_eq!(
            restored.xisf_properties,
            vec![XisfProperty::int32("Test:Kept", 3)]
        );
        assert_eq!(restored.get("CRPIX1").as_i64(), Some(512));
    }

    // XISF 1.0 section 11.5.3.7.2.2: an unknown name makes the solution
    // unavailable, so no keywords are made and the properties stay as stored.
    #[test]
    fn unusable_solution_properties_are_kept_verbatim_on_open() {
        let tmp = tmp_xisf();
        let properties = vec![
            XisfProperty::string(SOLUTION_IDS[0], "ZenithalEquidistant"),
            XisfProperty::f64_vector(SOLUTION_IDS[1], vec![37.5, 30.0]),
            XisfProperty::f64_vector(SOLUTION_IDS[2], vec![511.5, 255.25]),
            XisfProperty::f64_matrix(SOLUTION_IDS[3], 2, 2, vec![-1e-3, 0.0, 0.0, 1e-3]),
        ];
        let mut writer = XisfWriter::create(tmp.path()).unwrap();
        writer.add_image(&[0u8; 4], 2, 2, 1).unwrap();
        writer.add_properties(properties.clone());
        writer.write().unwrap();

        let img = Image::open(tmp.path()).unwrap();
        assert!(img.keywords.is_empty());
        assert_eq!(img.xisf_properties, properties);
    }

    // The FITSKeyword elements hold a WCS that Table 16 cannot name.
    #[test]
    fn a_wcs_without_an_xisf_projection_system_is_saved_as_keywords_only() {
        let tmp = tmp_xisf();
        let mut img = integer_valued_tan_image();
        img.set_keyword(Keyword::string("CTYPE1", "RA---ARC"));
        img.set_keyword(Keyword::string("CTYPE2", "DEC--ARC"));
        img.save(tmp.path()).unwrap();
        assert!(!header_xml(tmp.path()).contains("AstrometricSolution"));
        assert_eq!(Image::open(tmp.path()).unwrap().keywords, img.keywords);
    }

    // Saves an image with these keywords and reports whether the file holds
    // solution properties. The keywords must read back unchanged.
    fn saves_a_solution(keywords: Vec<Keyword>) -> bool {
        let tmp = tmp_xisf();
        let mut img = Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2]);
        img.keywords = keywords.clone();
        img.save(tmp.path()).unwrap();
        assert_eq!(Image::open(tmp.path()).unwrap().keywords, keywords);
        header_xml(tmp.path()).contains("AstrometricSolution")
    }

    fn tan_with(changes: Vec<Keyword>, removed: &[&str]) -> Vec<Keyword> {
        let mut keywords = integer_valued_tan_image().keywords;
        keywords.retain(|k| !removed.contains(&k.name.as_str()));
        for keyword in changes {
            keywords.retain(|k| k.name != keyword.name);
            keywords.push(keyword);
        }
        keywords
    }

    // IRAF and SCAMP write TNX and ZPX, which celestial-wcs does not build.
    #[test]
    fn wcs_keywords_that_do_not_build_are_saved_without_solution_properties() {
        let ctypes = |code: &str| {
            vec![
                Keyword::string("CTYPE1", format!("RA---{code}").as_str()),
                Keyword::string("CTYPE2", format!("DEC--{code}").as_str()),
            ]
        };
        let cd = ["CD1_1", "CD1_2", "CD2_1", "CD2_2"];
        let zero_cdelt = vec![Keyword::real("CDELT1", 0.0), Keyword::real("CDELT2", 0.0)];
        let mut impossible_pole = ctypes("CAR");
        impossible_pole.push(Keyword::real("LONPOLE", 90.0));
        let cases = [
            ("TNX", tan_with(ctypes("TNX"), &[])),
            ("ZPX", tan_with(ctypes("ZPX"), &[])),
            ("HPX", tan_with(ctypes("HPX"), &[])),
            ("SIP without A_ORDER", tan_with(ctypes("TAN-SIP"), &[])),
            ("CDELT 0", tan_with(zero_cdelt, &cd)),
            ("lone CTYPE1", tan_with(vec![], &["CTYPE2"])),
            (
                "LINEAR second axis",
                tan_with(vec![Keyword::string("CTYPE2", "LINEAR")], &[]),
            ),
            ("impossible LATPOLE", tan_with(impossible_pole, &[])),
        ];
        for (case, keywords) in cases {
            assert!(!saves_a_solution(keywords), "{case}");
        }
    }

    // Paper II: without RADESYS, an equinox before 1984 means FK4 and a later
    // one FK5. EPOCH is the older name for EQUINOX.
    #[test]
    fn equinox_without_radesys_selects_fk4_or_fk5() {
        let equinox = |name: &str, year: f64| tan_with(vec![Keyword::real(name, year)], &[]);
        assert!(!saves_a_solution(equinox("EQUINOX", 1950.0)));
        assert!(!saves_a_solution(equinox("EPOCH", 1950.0)));
        assert!(!saves_a_solution(equinox("EQUINOX", 1990.0)));
        assert!(saves_a_solution(equinox("EQUINOX", 2000.0)));
        assert!(saves_a_solution(equinox("EPOCH", 2000.0)));
    }

    // A WCS value of the wrong type, or a missing CRVAL, would be replaced by
    // a default, so no solution is written. The keywords stay.
    #[test]
    fn mistyped_or_missing_wcs_values_give_no_solution() {
        let cases = [
            tan_with(vec![Keyword::string("CRVAL1", "10.0")], &[]),
            tan_with(vec![Keyword::string("CRVAL2", "+20:00:00")], &[]),
            tan_with(vec![Keyword::string("CD1_1", "-0.001")], &[]),
            tan_with(vec![Keyword::real("CTYPE2", 1.0)], &[]),
            tan_with(vec![Keyword::string("EQUINOX", "J2000")], &[]),
            tan_with(vec![], &["CRVAL1"]),
            tan_with(vec![], &["CRVAL2"]),
        ];
        for keywords in cases {
            assert!(!saves_a_solution(keywords.clone()), "{keywords:?}");
        }
    }

    // XISF 1.0 section 7: one unreadable property leaves the rest available.
    #[test]
    fn a_malformed_property_leaves_the_image_readable() {
        let xml = concat!(
            r#"<xisf version="1.0"><Image geometry="2:2:1" sampleFormat="UInt8" location="attachment:4096:4">"#,
            r#"<Property id="bad" type="Boolean" value="yes"/>"#,
            r#"<Property id="good" type="Int32" value="7"/>"#,
            r#"<FITSKeyword name="OBJECT" value="'M31'" comment=""/></Image></xisf>"#
        );
        let mut bytes = b"XISF0100".to_vec();
        bytes.extend_from_slice(&(xml.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&[0; 4]);
        bytes.extend_from_slice(xml.as_bytes());
        bytes.resize(4096, 0);
        bytes.extend_from_slice(&[1, 2, 3, 4]);
        let tmp = tmp_xisf();
        std::fs::write(tmp.path(), bytes).unwrap();

        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_u8(), Some(&vec![1, 2, 3, 4]));
        assert_eq!(img.keywords, vec![Keyword::string("OBJECT", "M31")]);
        assert_eq!(img.xisf_properties, vec![XisfProperty::int32("good", 7)]);
    }

    // An XISF file whose FITSKeyword names include ones no FITS card can hold.
    fn open_xisf_with_unwritable_names(dir: &Path) -> Image {
        let long = "X".repeat(80);
        let keywords = format!(
            concat!(
                "<FITSKeyword name=\"A=B\" value=\"1\" comment=\"\"/>",
                "<FITSKeyword name=\"HIERARCH  PAD\" value=\"1\" comment=\"\"/>",
                "<FITSKeyword name=\"{long}\" value=\"1\" comment=\"\"/>",
                "<FITSKeyword name=\" \" value=\"3.0\" comment=\"c\"/>",
                "<FITSKeyword name=\"TAB&#9;\" value=\"'x&#9;y'\" comment=\"\u{e9}&#2;\"/>",
                "<FITSKeyword name=\"OBJECT\" value=\"'M31'\" comment=\"\"/>"
            ),
            long = long
        );
        let xml = format!(
            r#"<xisf version="1.0"><Image geometry="2:2:1" sampleFormat="UInt8" location="attachment:4096:4">{keywords}</Image></xisf>"#
        );
        let mut bytes = b"XISF0100".to_vec();
        bytes.extend_from_slice(&(xml.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&[0; 4]);
        bytes.extend_from_slice(xml.as_bytes());
        bytes.resize(4096, 0);
        bytes.extend_from_slice(&[1, 2, 3, 4]);
        std::fs::write(dir.join("in.xisf"), bytes).unwrap();
        Image::open(dir.join("in.xisf")).unwrap()
    }

    // The names a save writes for the image's keywords, without the cards
    // the FITS writer adds for the data unit and the FITS standard reference.
    fn saved_names(img: &Image, path: &Path) -> Vec<String> {
        let added = [
            "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND", "COMMENT",
        ];
        img.save(path).unwrap();
        Image::open(path)
            .unwrap()
            .keywords
            .into_iter()
            .map(|k| k.name)
            .filter(|name| !added.contains(&name.as_str()))
            .collect()
    }

    // Every keyword the reader accepts must save. A name that cannot be
    // written even as HIERARCH is left out.
    #[test]
    fn keywords_read_from_xisf_can_be_saved_again() {
        let dir = tempfile::tempdir().unwrap();
        let img = open_xisf_with_unwritable_names(dir.path());
        assert_eq!(
            saved_names(&img, &dir.path().join("out.xisf")),
            vec!["HIERARCH TAB?", "OBJECT"]
        );
    }

    #[test]
    fn keywords_read_from_xisf_can_be_saved_as_fits() {
        let dir = tempfile::tempdir().unwrap();
        let img = open_xisf_with_unwritable_names(dir.path());
        assert_eq!(
            saved_names(&img, &dir.path().join("out.fits")),
            vec!["HIERARCH TAB?", "OBJECT"]
        );
    }

    // A 1x1 UInt8 image at offset 4096, followed by 128 f64 values of 1.5
    // that `count` vector properties share.
    fn xisf_with_shared_block(count: usize) -> Vec<u8> {
        let property =
            r#"<Property id="v" type="F64Vector" length="128" location="attachment:4097:1024"/>"#;
        let xml = format!(
            r#"<xisf version="1.0"><Image geometry="1:1:1" sampleFormat="UInt8" location="attachment:4096:1">{}</Image></xisf>"#,
            property.repeat(count)
        );
        let mut file = b"XISF0100".to_vec();
        file.extend_from_slice(&(xml.len() as u32).to_le_bytes());
        file.extend_from_slice(&[0; 4]);
        file.extend_from_slice(xml.as_bytes());
        file.resize(4096, 0);
        file.push(7);
        for _ in 0..128 {
            file.extend_from_slice(&1.5f64.to_le_bytes());
        }
        file
    }

    // The file is 5121 bytes, so 5 copies of the 1024-byte block fit the
    // property budget. The image opens without the other 3.
    #[test]
    fn open_skips_properties_past_the_file_length_budget() {
        let tmp = tmp_xisf();
        std::fs::write(tmp.path(), xisf_with_shared_block(8)).unwrap();
        let img = Image::open(tmp.path()).unwrap();
        assert_eq!(img.pixels.as_u8(), Some(&vec![7]));
        let expected = XisfPropertyValue::F64Vector(vec![1.5; 128]);
        let values: Vec<&XisfPropertyValue> =
            img.xisf_properties.iter().map(|p| &p.value).collect();
        assert_eq!(values, vec![&expected; 5]);
    }

    #[test]
    fn save_refuses_an_empty_image_before_creating_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.xisf");
        let img = Image::new(PixelData::U8(vec![]), vec![0usize, 0]);
        assert!(matches!(
            img.save(&path),
            Err(ImageError::Xisf(XisfError::InvalidFormat(_)))
        ));
        assert!(!path.exists());
    }
}
