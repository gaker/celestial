use crate::core::types::ByteOrder;
use crate::ser::error::{Result, SerError};
use crate::ser::timestamp::TICKS_MASK;
use crate::ser::types::ColorId;
use celestial_core::location::Location;

pub(crate) const HEADER_SIZE: usize = 178;
pub(crate) const TIMESTAMP_SIZE: u64 = 8;
const FILE_ID: &[u8; 14] = b"LUCAM-RECORDER";
// The spec stores width and height as Int32. File sizes are checked
// separately, with checked arithmetic, in `frame_size` and `trailer_end`.
const MAX_DIMENSION: u32 = i32::MAX as u32;
const MAX_PIXEL_DEPTH: u32 = 16;

#[derive(Debug, Clone)]
pub struct SerHeader {
    pub lu_id: u32,
    pub color_id: ColorId,
    pub byte_order: ByteOrder,
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
    pub(crate) fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let bytes: &[u8; HEADER_SIZE] = bytes
            .get(..HEADER_SIZE)
            .and_then(|b| b.try_into().ok())
            .ok_or_else(|| SerError::InvalidHeader("Header too short".to_string()))?;
        if &bytes[0..14] != FILE_ID {
            return Err(SerError::InvalidFileId {
                actual: read_string(&bytes[0..14]),
            });
        }

        let header = Self {
            lu_id: le_u32(bytes, 14),
            color_id: ColorId::from_u32(le_u32(bytes, 18))?,
            // Capture software writes 1 for big-endian and 0 for little-endian,
            // the reverse of the spec text. Other values read as little-endian.
            byte_order: if le_u32(bytes, 22) == 1 {
                ByteOrder::BigEndian
            } else {
                ByteOrder::LittleEndian
            },
            image_width: le_u32(bytes, 26),
            image_height: le_u32(bytes, 30),
            pixel_depth_per_plane: le_u32(bytes, 34),
            frame_count: le_u32(bytes, 38),
            observer: read_string(&bytes[42..82]),
            instrument: read_string(&bytes[82..122]),
            telescope: read_string(&bytes[122..162]),
            date_time: date_ticks(bytes, 162),
            date_time_utc: date_ticks(bytes, 170),
            location: None,
        };
        header.validate()?;
        Ok(header)
    }

    pub(crate) fn to_bytes(&self) -> [u8; HEADER_SIZE] {
        let mut bytes = [0u8; HEADER_SIZE];
        let byte_order_flag: u32 = match self.byte_order {
            ByteOrder::LittleEndian => 0,
            ByteOrder::BigEndian => 1,
        };

        bytes[0..14].copy_from_slice(FILE_ID);
        bytes[14..18].copy_from_slice(&self.lu_id.to_le_bytes());
        bytes[18..22].copy_from_slice(&(self.color_id as u32).to_le_bytes());
        bytes[22..26].copy_from_slice(&byte_order_flag.to_le_bytes());
        bytes[26..30].copy_from_slice(&self.image_width.to_le_bytes());
        bytes[30..34].copy_from_slice(&self.image_height.to_le_bytes());
        bytes[34..38].copy_from_slice(&self.pixel_depth_per_plane.to_le_bytes());
        bytes[38..42].copy_from_slice(&self.frame_count.to_le_bytes());
        write_string(&mut bytes[42..82], &self.observer);
        write_string(&mut bytes[82..122], &self.instrument);
        write_string(&mut bytes[122..162], &self.telescope);
        bytes[162..170].copy_from_slice(&self.date_time.to_le_bytes());
        bytes[170..178].copy_from_slice(&self.date_time_utc.to_le_bytes());

        bytes
    }

    pub fn bytes_per_pixel(&self) -> u32 {
        let bytes_per_plane = if self.pixel_depth_per_plane <= 8 {
            1
        } else {
            2
        };
        bytes_per_plane * self.color_id.planes()
    }

    pub fn frame_size(&self) -> Result<u64> {
        u64::from(self.image_width)
            .checked_mul(u64::from(self.image_height))
            .and_then(|pixels| pixels.checked_mul(u64::from(self.bytes_per_pixel())))
            .ok_or_else(|| self.size_overflow())
    }

    // The trailer starts here, so this is also the length of a file without one.
    pub(crate) fn frames_end(&self) -> Result<u64> {
        self.frame_offset(self.frame_count)
    }

    pub(crate) fn trailer_end(&self) -> Result<u64> {
        self.timestamp_offset(self.frame_count)
    }

    pub(crate) fn frame_offset(&self, frame_id: u32) -> Result<u64> {
        self.frame_size()?
            .checked_mul(u64::from(frame_id))
            .and_then(|frames| frames.checked_add(HEADER_SIZE as u64))
            .ok_or_else(|| self.size_overflow())
    }

    pub(crate) fn timestamp_offset(&self, frame_id: u32) -> Result<u64> {
        self.frames_end()?
            .checked_add(TIMESTAMP_SIZE * u64::from(frame_id))
            .ok_or_else(|| self.size_overflow())
    }

    pub(crate) fn validate(&self) -> Result<()> {
        let (width, height) = (self.image_width, self.image_height);
        if width == 0 || height == 0 || width > MAX_DIMENSION || height > MAX_DIMENSION {
            return Err(SerError::InvalidDimensions { width, height });
        }
        let depth = self.pixel_depth_per_plane;
        if depth == 0 || depth > MAX_PIXEL_DEPTH {
            return Err(SerError::InvalidPixelDepth(depth));
        }
        self.trailer_end().map(|_| ())
    }

    fn size_overflow(&self) -> SerError {
        SerError::InvalidHeader(format!(
            "{} frames of {}x{} {} at depth {} overflow a 64-bit file size",
            self.frame_count,
            self.image_width,
            self.image_height,
            self.color_id,
            self.pixel_depth_per_plane
        ))
    }

    pub fn with_location(mut self, location: Location) -> Self {
        self.location = Some(location);
        self
    }

    pub fn set_location(&mut self, location: Location) {
        self.location = Some(location);
    }
}

impl Default for SerHeader {
    fn default() -> Self {
        Self {
            lu_id: 0,
            color_id: ColorId::Mono,
            byte_order: ByteOrder::LittleEndian,
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

fn le_u32(bytes: &[u8; HEADER_SIZE], offset: usize) -> u32 {
    let mut field = [0u8; 4];
    field.copy_from_slice(&bytes[offset..offset + 4]);
    u32::from_le_bytes(field)
}

// Masking the Kind bits leaves 62 bits, so the result is never negative.
fn date_ticks(bytes: &[u8; HEADER_SIZE], offset: usize) -> i64 {
    let mut field = [0u8; 8];
    field.copy_from_slice(&bytes[offset..offset + 8]);
    (u64::from_le_bytes(field) & TICKS_MASK) as i64
}

// Header strings are read as Latin-1, where each byte is the code point of one
// char, so no byte is lost. The text ends at the first NUL, and trailing
// spaces are trimmed.
fn read_string(bytes: &[u8]) -> String {
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    let text: String = bytes[..end].iter().map(|&b| char::from(b)).collect();
    text.trim_end_matches(' ').to_string()
}

// Readers decode these fields in different 8-bit encodings, and a NUL ends
// the field, so only printable ASCII is written and any other char becomes
// '?'. A string longer than the field is cut, and the rest of the field is NUL.
fn write_string(field: &mut [u8], text: &str) {
    field.fill(0);
    for (byte, c) in field.iter_mut().zip(text.chars()) {
        *byte = u8::try_from(c)
            .ok()
            .filter(|b| (b' '..=b'~').contains(b))
            .unwrap_or(b'?');
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shown<T>(result: Result<T>) -> std::result::Result<T, String> {
        result.map_err(|e| e.to_string())
    }

    fn create_valid_header_bytes() -> Vec<u8> {
        let mut bytes = vec![0u8; HEADER_SIZE];

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
        let result = SerHeader::from_bytes(&[0u8; HEADER_SIZE - 1]);
        assert!(matches!(result, Err(SerError::InvalidHeader(m)) if m == "Header too short"));
    }

    #[test]
    fn from_bytes_invalid_file_id() {
        let mut bytes = create_valid_header_bytes();
        bytes[0..14].copy_from_slice(b"INVALID-HEADER");
        let result = SerHeader::from_bytes(&bytes);
        assert!(matches!(
            result,
            Err(SerError::InvalidFileId { actual }) if actual == "INVALID-HEADER"
        ));
    }

    #[test]
    fn from_bytes_successful_parse() {
        let bytes = create_valid_header_bytes();
        let header = SerHeader::from_bytes(&bytes).unwrap();

        assert_eq!(header.lu_id, 123);
        assert_eq!(header.color_id, ColorId::BayerRggb);
        assert_eq!(header.byte_order, ByteOrder::BigEndian);
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
    fn byte_order_flag_one_is_big_endian() {
        let mut bytes = create_valid_header_bytes();
        for (flag, order) in [
            (0u32, ByteOrder::LittleEndian),
            (1, ByteOrder::BigEndian),
            (2, ByteOrder::LittleEndian),
            (7, ByteOrder::LittleEndian),
            (u32::MAX, ByteOrder::LittleEndian),
        ] {
            bytes[22..26].copy_from_slice(&flag.to_le_bytes());
            let header = SerHeader::from_bytes(&bytes).unwrap();
            assert_eq!(header.byte_order, order, "flag {flag}");
        }
        let little = SerHeader::default();
        assert_eq!(little.to_bytes()[22..26], [0, 0, 0, 0]);
        let big = SerHeader {
            byte_order: ByteOrder::BigEndian,
            ..SerHeader::default()
        };
        assert_eq!(big.to_bytes()[22..26], [1, 0, 0, 0]);
    }

    #[test]
    fn dates_have_the_kind_bits_masked() {
        // 2024-06-01T12:34:56.1234567 in .NET ticks.
        let ticks: u64 = 638_528_420_961_234_567;
        let mut bytes = create_valid_header_bytes();
        bytes[162..170].copy_from_slice(&(ticks | 1 << 63).to_le_bytes());
        bytes[170..178].copy_from_slice(&(ticks | 1 << 62).to_le_bytes());
        let header = SerHeader::from_bytes(&bytes).unwrap();
        assert_eq!(header.date_time, ticks as i64);
        assert_eq!(header.date_time_utc, ticks as i64);
    }

    #[test]
    fn to_bytes_serialization() {
        let header = SerHeader {
            lu_id: 456,
            color_id: ColorId::Rgb,
            byte_order: ByteOrder::BigEndian,
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
        assert_eq!(&bytes[0..14], b"LUCAM-RECORDER");

        let deserialized = SerHeader::from_bytes(&bytes).unwrap();
        assert_eq!(deserialized.lu_id, 456);
        assert_eq!(deserialized.color_id, ColorId::Rgb);
        assert_eq!(deserialized.byte_order, ByteOrder::BigEndian);
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
    fn read_string_stops_at_the_first_nul() {
        assert_eq!(read_string(b"Hello\0World"), "Hello");
        assert_eq!(read_string(b"Hello"), "Hello");
        assert_eq!(read_string(b"\0Hello"), "");
    }

    #[test]
    fn read_string_trims_trailing_spaces() {
        assert_eq!(read_string(b"  Scope   \0\0"), "  Scope");
        assert_eq!(read_string(b"    "), "");
    }

    #[test]
    fn read_string_is_latin1() {
        assert_eq!(read_string(b"\xC5k\0\0"), "\u{C5}k");
        assert_eq!(read_string(b"\xFF\xFE\xFD"), "\u{FF}\u{FE}\u{FD}");
        let every_byte: Vec<u8> = (1..=255).collect();
        let text = read_string(&every_byte);
        let codes: Vec<u32> = text.chars().map(u32::from).collect();
        assert_eq!(codes, (1..=255).collect::<Vec<u32>>());
    }

    #[test]
    fn latin1_observer_in_a_header() {
        let mut bytes = create_valid_header_bytes();
        bytes[42..82].fill(0);
        bytes[42] = 0xC5;
        bytes[43] = b'k';
        let header = SerHeader::from_bytes(&bytes).unwrap();
        assert_eq!(header.observer, "\u{C5}k");
        assert_eq!(header.to_bytes()[42..44], [b'?', b'k']);
    }

    #[test]
    fn write_string_normal_case() {
        let mut buf = [0u8; 10];
        write_string(&mut buf, "Hello");
        assert_eq!(&buf, b"Hello\0\0\0\0\0");
    }

    #[test]
    fn write_string_truncation() {
        let mut buf = [0u8; 5];
        write_string(&mut buf, "Hello World");
        assert_eq!(&buf, b"Hello");
    }

    #[test]
    fn write_string_zero_padding() {
        let mut buf = [255u8; 10];
        write_string(&mut buf, "Hi");
        assert_eq!(&buf, b"Hi\0\0\0\0\0\0\0\0");
    }

    #[test]
    fn write_string_writes_printable_ascii() {
        let mut buf = [0u8; 8];
        write_string(&mut buf, "\u{C5}ngstr\u{F6}m");
        assert_eq!(&buf, b"?ngstr?m");
        write_string(&mut buf, "a\0b\x7F\tc~ ");
        assert_eq!(&buf, b"a?b??c~ ");
        write_string(&mut buf, "M\u{2083}\u{1F52D}x");
        assert_eq!(&buf, b"M??x\0\0\0\0");
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
        let cases = [
            (ColorId::Mono, 1, 1),
            (ColorId::Mono, 8, 1),
            (ColorId::BayerGrbg, 9, 2),
            (ColorId::Mono, 16, 2),
            (ColorId::Rgb, 8, 3),
            (ColorId::Bgr, 12, 6),
            (ColorId::Rgb, 16, 6),
        ];
        for (color_id, depth, bytes) in cases {
            let header = SerHeader {
                color_id,
                pixel_depth_per_plane: depth,
                ..Default::default()
            };
            assert_eq!(header.bytes_per_pixel(), bytes, "{color_id} depth {depth}");
        }
    }

    #[test]
    fn size_calculation() {
        let header = SerHeader {
            image_width: 100,
            image_height: 100,
            pixel_depth_per_plane: 16,
            color_id: ColorId::Mono,
            frame_count: 10,
            ..Default::default()
        };
        assert_eq!(shown(header.frame_size()), Ok(20000));
        assert_eq!(shown(header.frames_end()), Ok(178 + 10 * 20000));
        assert_eq!(shown(header.trailer_end()), Ok(178 + 10 * 20000 + 80));
    }

    #[test]
    fn overflowing_sizes_are_errors() {
        let huge = SerHeader {
            image_width: u32::MAX,
            image_height: u32::MAX,
            color_id: ColorId::Rgb,
            frame_count: 1,
            ..Default::default()
        };
        let message = "Invalid SER header: 1 frames of 4294967295x4294967295 RGB at depth 16 \
                       overflow a 64-bit file size";
        assert_eq!(shown(huge.frame_size()), Err(message.to_string()));
        assert_eq!(shown(huge.trailer_end()), Err(message.to_string()));

        let many = SerHeader {
            image_width: 65535,
            image_height: 65535,
            color_id: ColorId::Rgb,
            frame_count: u32::MAX,
            ..Default::default()
        };
        assert_eq!(shown(many.frame_size()), Ok(25_769_017_350));
        let message = "Invalid SER header: 4294967295 frames of 65535x65535 RGB at depth 16 \
                       overflow a 64-bit file size";
        assert_eq!(shown(many.frames_end()), Err(message.to_string()));
        assert_eq!(
            shown(SerHeader::from_bytes(&many.to_bytes()).map(|_| ())),
            Err(message.to_string())
        );
    }

    #[test]
    fn trailer_end_overflow_alone_is_an_error() {
        // frames_end fits in u64, and adding 8 bytes per frame does not.
        let header = SerHeader {
            image_width: 65535,
            image_height: 65535,
            frame_count: 2_147_549_185,
            ..Default::default()
        };
        assert_eq!(shown(header.frames_end()), Ok(18_446_744_069_414_453_428));
        assert_eq!(
            shown(header.trailer_end()),
            Err(
                "Invalid SER header: 2147549185 frames of 65535x65535 MONO at depth 16 \
                 overflow a 64-bit file size"
                    .to_string()
            )
        );
    }

    #[test]
    fn debug_and_clone() {
        let header = SerHeader::default();
        let cloned = header.clone();
        assert_eq!(header.to_bytes(), cloned.to_bytes());

        let debug_str = format!("{:?}", header);
        assert!(debug_str.contains("SerHeader"));
    }

    #[test]
    fn round_trip_serialization() {
        let original = create_valid_header_bytes();
        let parsed = SerHeader::from_bytes(&original).unwrap();
        assert_eq!(parsed.to_bytes().as_slice(), original.as_slice());
    }

    #[test]
    fn from_bytes_invalid_dimensions() {
        let too_big = 1u32 << 31;
        for (width, height) in [(0u32, 480u32), (640, 0), (too_big, 480), (640, too_big)] {
            let mut bytes = create_valid_header_bytes();
            bytes[26..30].copy_from_slice(&width.to_le_bytes());
            bytes[30..34].copy_from_slice(&height.to_le_bytes());
            let result = SerHeader::from_bytes(&bytes);
            assert!(
                matches!(result, Err(SerError::InvalidDimensions { width: w, height: h }) if (w, h) == (width, height)),
                "{width}x{height}"
            );
        }
    }

    // The spec stores width and height as Int32.
    #[test]
    fn from_bytes_accepts_dimensions_up_to_i32_max() {
        let max = i32::MAX as u32;
        for (width, height) in [(70_000u32, 480u32), (640, 70_000), (max, 1), (1, max)] {
            let mut bytes = create_valid_header_bytes();
            bytes[26..30].copy_from_slice(&width.to_le_bytes());
            bytes[30..34].copy_from_slice(&height.to_le_bytes());
            let header = SerHeader::from_bytes(&bytes).unwrap();
            assert_eq!((header.image_width, header.image_height), (width, height));
        }
    }

    #[test]
    fn from_bytes_invalid_pixel_depth() {
        for depth in [0u32, 17, 32] {
            let mut bytes = create_valid_header_bytes();
            bytes[34..38].copy_from_slice(&depth.to_le_bytes());
            let result = SerHeader::from_bytes(&bytes);
            assert!(
                matches!(result, Err(SerError::InvalidPixelDepth(d)) if d == depth),
                "{depth}"
            );
        }
    }
}
