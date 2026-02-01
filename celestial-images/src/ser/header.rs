use crate::ser::{ColorId, Result, SerError};
use celestial_core::Location;
use byteorder::{ByteOrder, LittleEndian};
use std::ffi::CStr;

#[derive(Debug, Clone)]
pub struct SerHeader {
    pub file_id: String,
    pub lu_id: u32,
    pub color_id: ColorId,
    pub little_endian: bool,
    pub image_width: u32,
    pub image_height: u32,
    pub pixel_depth_per_plane: u32,
    pub frame_count: u32,
    pub observer: String,
    pub instrument: String,
    pub telescope: String,
    pub date_time: i64,
    pub date_time_utc: i64,
    pub location: Option<Location>,
}

impl SerHeader {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        if bytes.len() < 178 {
            return Err(SerError::InvalidHeader("Header too short".to_string()));
        }

        let file_id = Self::read_string(&bytes[0..14])?;
        if file_id != "LUCAM-RECORDER" {
            return Err(SerError::InvalidFileId { actual: file_id });
        }

        let lu_id = LittleEndian::read_u32(&bytes[14..18]);
        let color_id = ColorId::from_u32(LittleEndian::read_u32(&bytes[18..22]))?;
        let little_endian = LittleEndian::read_u32(&bytes[22..26]) != 0;
        let image_width = LittleEndian::read_u32(&bytes[26..30]);
        let image_height = LittleEndian::read_u32(&bytes[30..34]);
        let pixel_depth_per_plane = LittleEndian::read_u32(&bytes[34..38]);
        let frame_count = LittleEndian::read_u32(&bytes[38..42]);
        let observer = Self::read_string(&bytes[42..82])?;
        let instrument = Self::read_string(&bytes[82..122])?;
        let telescope = Self::read_string(&bytes[122..162])?;
        let date_time = LittleEndian::read_i64(&bytes[162..170]);
        let date_time_utc = LittleEndian::read_i64(&bytes[170..178]);

        crate::ser::types::SerFile::validate_dimensions(image_width, image_height)?;
        crate::ser::types::SerFile::validate_pixel_depth(pixel_depth_per_plane)?;

        Ok(Self {
            file_id,
            lu_id,
            color_id,
            little_endian,
            image_width,
            image_height,
            pixel_depth_per_plane,
            frame_count,
            observer,
            instrument,
            telescope,
            date_time,
            date_time_utc,
            location: None,
        })
    }

    pub fn to_bytes(&self) -> [u8; 178] {
        let mut bytes = [0u8; 178];

        Self::write_string(&mut bytes[0..14], &self.file_id);
        LittleEndian::write_u32(&mut bytes[14..18], self.lu_id);
        LittleEndian::write_u32(&mut bytes[18..22], self.color_id as u32);
        LittleEndian::write_u32(&mut bytes[22..26], if self.little_endian { 1 } else { 0 });
        LittleEndian::write_u32(&mut bytes[26..30], self.image_width);
        LittleEndian::write_u32(&mut bytes[30..34], self.image_height);
        LittleEndian::write_u32(&mut bytes[34..38], self.pixel_depth_per_plane);
        LittleEndian::write_u32(&mut bytes[38..42], self.frame_count);
        Self::write_string(&mut bytes[42..82], &self.observer);
        Self::write_string(&mut bytes[82..122], &self.instrument);
        Self::write_string(&mut bytes[122..162], &self.telescope);
        LittleEndian::write_i64(&mut bytes[162..170], self.date_time);
        LittleEndian::write_i64(&mut bytes[170..178], self.date_time_utc);

        bytes
    }

    pub fn bytes_per_pixel(&self) -> u32 {
        crate::ser::types::SerFile::calculate_bytes_per_pixel(
            self.pixel_depth_per_plane,
            self.color_id.planes(),
        )
    }

    pub fn frame_size(&self) -> u64 {
        crate::ser::types::SerFile::calculate_frame_size(
            self.image_width,
            self.image_height,
            self.bytes_per_pixel(),
        )
    }

    pub fn has_trailer(&self) -> bool {
        self.date_time > 0
    }

    fn read_string(bytes: &[u8]) -> Result<String> {
        let null_pos = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
        let valid_bytes = &bytes[..null_pos];

        match std::str::from_utf8(valid_bytes) {
            Ok(s) => Ok(s.to_string()),
            Err(_) => match CStr::from_bytes_with_nul(bytes) {
                Ok(c_str) => Ok(c_str.to_string_lossy().to_string()),
                Err(_) => {
                    let mut padded = bytes.to_vec();
                    padded.push(0);
                    match CStr::from_bytes_with_nul(&padded) {
                        Ok(c_str) => Ok(c_str.to_string_lossy().to_string()),
                        Err(_) => Err(SerError::InvalidHeader("Invalid string data".to_string())),
                    }
                }
            },
        }
    }

    fn write_string(buf: &mut [u8], s: &str) {
        let bytes = s.as_bytes();
        let copy_len = bytes.len().min(buf.len());
        buf[..copy_len].copy_from_slice(&bytes[..copy_len]);

        for byte in buf.iter_mut().skip(copy_len) {
            *byte = 0;
        }
    }
}

impl Default for SerHeader {
    fn default() -> Self {
        Self {
            file_id: "LUCAM-RECORDER".to_string(),
            lu_id: 0,
            color_id: ColorId::Mono,
            little_endian: true,
            image_width: 1,
            image_height: 1,
            pixel_depth_per_plane: 16,
            frame_count: 0,
            observer: String::new(),
            instrument: String::new(),
            telescope: String::new(),
            date_time: 0,
            date_time_utc: 0,
            location: None,
        }
    }
}

impl SerHeader {
    pub fn with_location(mut self, location: Location) -> Self {
        self.location = Some(location);
        self
    }

    pub fn set_location(&mut self, location: Location) {
        self.location = Some(location);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_valid_header_bytes() -> Vec<u8> {
        let mut bytes = vec![0u8; 178];

        bytes[0..14].copy_from_slice(b"LUCAM-RECORDER");
        bytes[14..18].copy_from_slice(&123u32.to_le_bytes());
        bytes[18..22].copy_from_slice(&8u32.to_le_bytes());
        bytes[22..26].copy_from_slice(&1u32.to_le_bytes());
        bytes[26..30].copy_from_slice(&640u32.to_le_bytes());
        bytes[30..34].copy_from_slice(&480u32.to_le_bytes());
        bytes[34..38].copy_from_slice(&16u32.to_le_bytes());
        bytes[38..42].copy_from_slice(&100u32.to_le_bytes());

        bytes[42..50].copy_from_slice(b"Observer");
        bytes[82..92].copy_from_slice(b"Instrument");
        bytes[122..131].copy_from_slice(b"Telescope");
        bytes[162..170].copy_from_slice(&1234567890i64.to_le_bytes());
        bytes[170..178].copy_from_slice(&1234567891i64.to_le_bytes());

        bytes
    }

    #[test]
    fn from_bytes_header_too_short() {
        let short_bytes = vec![0u8; 100];
        let result = SerHeader::from_bytes(&short_bytes);
        assert!(matches!(result, Err(SerError::InvalidHeader(_))));
    }

    #[test]
    fn from_bytes_invalid_file_id() {
        let mut bytes = vec![0u8; 178];
        bytes[0..14].copy_from_slice(b"INVALID-HEADER");
        bytes[26..30].copy_from_slice(&1u32.to_le_bytes());
        bytes[30..34].copy_from_slice(&1u32.to_le_bytes());
        bytes[34..38].copy_from_slice(&16u32.to_le_bytes());

        let result = SerHeader::from_bytes(&bytes);
        assert!(matches!(result, Err(SerError::InvalidFileId { .. })));
    }

    #[test]
    fn from_bytes_successful_parse() {
        let bytes = create_valid_header_bytes();
        let header = SerHeader::from_bytes(&bytes).unwrap();

        assert_eq!(header.file_id, "LUCAM-RECORDER");
        assert_eq!(header.lu_id, 123);
        assert_eq!(header.color_id, ColorId::BayerRggb);
        assert!(header.little_endian);
        assert_eq!(header.image_width, 640);
        assert_eq!(header.image_height, 480);
        assert_eq!(header.pixel_depth_per_plane, 16);
        assert_eq!(header.frame_count, 100);
        assert_eq!(header.observer, "Observer");
        assert_eq!(header.instrument, "Instrument");
        assert_eq!(header.telescope, "Telescope");
        assert_eq!(header.date_time, 1234567890);
        assert_eq!(header.date_time_utc, 1234567891);
        assert!(header.location.is_none());
    }

    #[test]
    fn to_bytes_serialization() {
        let header = SerHeader {
            lu_id: 456,
            color_id: ColorId::Rgb,
            little_endian: false,
            image_width: 1920,
            image_height: 1080,
            pixel_depth_per_plane: 8,
            frame_count: 200,
            observer: "Test Observer".to_string(),
            instrument: "Test Camera".to_string(),
            telescope: "Test Scope".to_string(),
            date_time: 9876543210,
            date_time_utc: 9876543211,
            ..Default::default()
        };

        let bytes = header.to_bytes();
        assert_eq!(bytes.len(), 178);

        let deserialized = SerHeader::from_bytes(&bytes).unwrap();
        assert_eq!(deserialized.lu_id, 456);
        assert_eq!(deserialized.color_id, ColorId::Rgb);
        assert!(!deserialized.little_endian);
        assert_eq!(deserialized.image_width, 1920);
        assert_eq!(deserialized.image_height, 1080);
        assert_eq!(deserialized.pixel_depth_per_plane, 8);
        assert_eq!(deserialized.frame_count, 200);
        assert_eq!(deserialized.observer, "Test Observer");
        assert_eq!(deserialized.instrument, "Test Camera");
        assert_eq!(deserialized.telescope, "Test Scope");
        assert_eq!(deserialized.date_time, 9876543210);
        assert_eq!(deserialized.date_time_utc, 9876543211);
    }

    #[test]
    fn has_trailer_logic() {
        let header_zero = SerHeader {
            date_time: 0,
            ..Default::default()
        };
        assert!(!header_zero.has_trailer());

        let header_positive = SerHeader {
            date_time: 1,
            ..Default::default()
        };
        assert!(header_positive.has_trailer());

        let header_negative = SerHeader {
            date_time: -1,
            ..Default::default()
        };
        assert!(!header_negative.has_trailer());
    }

    #[test]
    fn read_string_valid_utf8() {
        let test_string = b"Hello\0World";
        let result = SerHeader::read_string(test_string).unwrap();
        assert_eq!(result, "Hello");
    }

    #[test]
    fn read_string_invalid_utf8_with_null() {
        let invalid_utf8 = b"\xFF\xFE\xFD\0";
        let result = SerHeader::read_string(invalid_utf8).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn read_string_invalid_utf8_no_null() {
        let invalid_utf8 = b"\xFF\xFE\xFD";
        let result = SerHeader::read_string(invalid_utf8).unwrap();
        assert!(!result.is_empty());
    }

    #[test]
    fn read_string_completely_invalid() {
        let completely_invalid = vec![255; 100];
        let result = SerHeader::read_string(&completely_invalid);
        assert!(result.is_ok());
    }

    #[test]
    fn write_string_normal_case() {
        let mut buf = [0u8; 10];
        SerHeader::write_string(&mut buf, "Hello");
        assert_eq!(&buf[0..5], b"Hello");
        assert_eq!(&buf[5..], &[0, 0, 0, 0, 0]);
    }

    #[test]
    fn write_string_truncation() {
        let mut buf = [0u8; 5];
        SerHeader::write_string(&mut buf, "Hello World");
        assert_eq!(&buf, b"Hello");
    }

    #[test]
    fn write_string_zero_padding() {
        let mut buf = [255u8; 10];
        SerHeader::write_string(&mut buf, "Hi");
        assert_eq!(&buf[0..2], b"Hi");
        assert_eq!(&buf[2..], &[0, 0, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn with_location() {
        let location = Location::from_degrees(40.0, -74.0, 100.0).unwrap();
        let header = SerHeader::default().with_location(location);
        assert!(header.location.is_some());
    }

    #[test]
    fn set_location() {
        let location = Location::from_degrees(51.5, -0.1, 25.0).unwrap();
        let mut header = SerHeader::default();
        header.set_location(location);
        assert!(header.location.is_some());
    }

    #[test]
    fn bytes_per_pixel_calculation() {
        let header_mono = SerHeader {
            pixel_depth_per_plane: 8,
            color_id: ColorId::Mono,
            ..Default::default()
        };
        assert_eq!(header_mono.bytes_per_pixel(), 1);

        let header_rgb = SerHeader {
            pixel_depth_per_plane: 16,
            color_id: ColorId::Rgb,
            ..Default::default()
        };
        assert_eq!(header_rgb.bytes_per_pixel(), 6);
    }

    #[test]
    fn frame_size_calculation() {
        let header = SerHeader {
            image_width: 100,
            image_height: 100,
            pixel_depth_per_plane: 16,
            color_id: ColorId::Mono,
            ..Default::default()
        };
        assert_eq!(header.frame_size(), 20000);
    }

    #[test]
    fn debug_and_clone() {
        let header = SerHeader::default();
        let cloned = header.clone();
        assert_eq!(header.file_id, cloned.file_id);

        let debug_str = format!("{:?}", header);
        assert!(debug_str.contains("SerHeader"));
    }

    #[test]
    fn round_trip_serialization() {
        let original = create_valid_header_bytes();
        let parsed = SerHeader::from_bytes(&original).unwrap();
        let serialized = parsed.to_bytes();
        let reparsed = SerHeader::from_bytes(&serialized).unwrap();

        assert_eq!(parsed.file_id, reparsed.file_id);
        assert_eq!(parsed.image_width, reparsed.image_width);
        assert_eq!(parsed.image_height, reparsed.image_height);
        assert_eq!(parsed.observer, reparsed.observer);
    }

    #[test]
    fn from_bytes_invalid_dimensions() {
        let mut bytes = vec![0u8; 178];
        bytes[0..14].copy_from_slice(b"LUCAM-RECORDER");
        bytes[26..30].copy_from_slice(&0u32.to_le_bytes());
        bytes[30..34].copy_from_slice(&480u32.to_le_bytes());
        bytes[34..38].copy_from_slice(&16u32.to_le_bytes());

        let result = SerHeader::from_bytes(&bytes);
        assert!(result.is_err());
    }

    #[test]
    fn read_string_completely_invalid_cstr() {
        let invalid_bytes = [255u8; 4];
        let result = SerHeader::read_string(&invalid_bytes);
        assert!(result.is_ok());
    }
}
