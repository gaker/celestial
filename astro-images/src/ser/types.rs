use crate::ser::SerError;
use astro_core::constants::SECONDS_PER_DAY_F64;
use astro_time::{constants::UNIX_EPOCH_JD, TimeError, TimeResult, UTC};
use byteorder::{ByteOrder, LittleEndian};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum ColorId {
    Mono = 0,
    BayerRggb = 8,
    BayerGrbg = 9,
    BayerGbrg = 10,
    BayerBggr = 11,
    BayerCyym = 16,
    BayerYcmy = 17,
    BayerYmcy = 18,
    BayerMyyc = 19,
    Rgb = 100,
    Bgr = 101,
}

impl ColorId {
    pub fn from_u32(value: u32) -> Result<Self, SerError> {
        match value {
            0 => Ok(Self::Mono),
            8 => Ok(Self::BayerRggb),
            9 => Ok(Self::BayerGrbg),
            10 => Ok(Self::BayerGbrg),
            11 => Ok(Self::BayerBggr),
            16 => Ok(Self::BayerCyym),
            17 => Ok(Self::BayerYcmy),
            18 => Ok(Self::BayerYmcy),
            19 => Ok(Self::BayerMyyc),
            100 => Ok(Self::Rgb),
            101 => Ok(Self::Bgr),
            _ => Err(SerError::UnsupportedColorFormat(value)),
        }
    }

    pub fn planes(self) -> u32 {
        match self {
            Self::Mono
            | Self::BayerRggb
            | Self::BayerGrbg
            | Self::BayerGbrg
            | Self::BayerBggr
            | Self::BayerCyym
            | Self::BayerYcmy
            | Self::BayerYmcy
            | Self::BayerMyyc => 1,
            Self::Rgb | Self::Bgr => 3,
        }
    }

    pub fn is_bayer(self) -> bool {
        matches!(
            self,
            Self::BayerRggb
                | Self::BayerGrbg
                | Self::BayerGbrg
                | Self::BayerBggr
                | Self::BayerCyym
                | Self::BayerYcmy
                | Self::BayerYmcy
                | Self::BayerMyyc
        )
    }
}

impl fmt::Display for ColorId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::Mono => "MONO",
            Self::BayerRggb => "BAYER_RGGB",
            Self::BayerGrbg => "BAYER_GRBG",
            Self::BayerGbrg => "BAYER_GBRG",
            Self::BayerBggr => "BAYER_BGGR",
            Self::BayerCyym => "BAYER_CYYM",
            Self::BayerYcmy => "BAYER_YCMY",
            Self::BayerYmcy => "BAYER_YMCY",
            Self::BayerMyyc => "BAYER_MYYC",
            Self::Rgb => "RGB",
            Self::Bgr => "BGR",
        };
        write!(f, "{}", name)
    }
}

pub type SerFrameId = u32;

#[derive(Debug, Clone)]
pub struct SerTimestamp {
    pub ticks: u64,
    pub utc: Option<UTC>,
}

impl SerTimestamp {
    pub const TICKS_PER_SECOND: u64 = 10_000_000;
    pub const UNIX_EPOCH_TICKS: u64 = 621_355_968_000_000_000;
    pub const DOTNET_EPOCH_TICKS: u64 = 0;

    pub fn new(ticks: u64) -> Self {
        let utc = Self::ticks_to_utc(ticks).ok();
        Self { ticks, utc }
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        Self::new(LittleEndian::read_u64(bytes))
    }

    pub fn to_bytes(&self, buf: &mut [u8]) {
        LittleEndian::write_u64(buf, self.ticks);
    }

    pub fn from_astro_time(utc: UTC) -> TimeResult<Self> {
        let (unix_seconds, unix_nanos) = utc_to_unix_timestamp(&utc)?;
        let unix_total = unix_seconds as f64 + unix_nanos as f64 / 1e9;
        let ticks = Self::unix_seconds_to_ticks(unix_total)?;
        Ok(Self {
            ticks,
            utc: Some(utc),
        })
    }

    pub fn to_utc(&self) -> TimeResult<UTC> {
        if let Some(utc) = &self.utc {
            Ok(*utc)
        } else {
            Self::ticks_to_utc(self.ticks)
        }
    }

    pub fn to_unix_seconds(&self) -> f64 {
        if self.ticks >= Self::UNIX_EPOCH_TICKS {
            (self.ticks - Self::UNIX_EPOCH_TICKS) as f64 / Self::TICKS_PER_SECOND as f64
        } else {
            0.0
        }
    }

    pub fn from_unix_seconds(seconds: f64) -> TimeResult<Self> {
        let ticks = Self::unix_seconds_to_ticks(seconds)?;
        let unix_secs = seconds.trunc() as i64;
        let unix_nanos = ((seconds.fract() * 1e9) as u32).min(999_999_999);
        let utc = UTC::new(unix_secs, unix_nanos);
        Ok(Self {
            ticks,
            utc: Some(utc),
        })
    }

    fn unix_seconds_to_ticks(seconds: f64) -> TimeResult<u64> {
        if !seconds.is_finite() || seconds < 0.0 {
            return Err(TimeError::CalculationError(
                "Invalid unix timestamp".to_string(),
            ));
        }

        let ticks = (seconds * Self::TICKS_PER_SECOND as f64) as u64 + Self::UNIX_EPOCH_TICKS;
        Ok(ticks)
    }

    fn ticks_to_utc(ticks: u64) -> TimeResult<UTC> {
        if ticks < Self::UNIX_EPOCH_TICKS {
            return Err(TimeError::CalculationError(
                "Timestamp before Unix epoch".to_string(),
            ));
        }

        let delta_ticks = ticks - Self::UNIX_EPOCH_TICKS;
        let seconds_since_epoch = delta_ticks / Self::TICKS_PER_SECOND;
        if seconds_since_epoch > i64::MAX as u64 {
            return Err(TimeError::CalculationError(
                "Timestamp exceeds supported range".to_string(),
            ));
        }

        let fractional_ticks = delta_ticks % Self::TICKS_PER_SECOND;
        let nanos = (fractional_ticks * 100) as u32;

        Ok(UTC::new(seconds_since_epoch as i64, nanos))
    }

    pub fn precision_100ns() -> u64 {
        1
    }
}

fn utc_to_unix_timestamp(utc: &UTC) -> TimeResult<(i64, u32)> {
    let jd = utc.to_julian_date().to_f64();
    if jd < UNIX_EPOCH_JD {
        return Err(TimeError::CalculationError(
            "Timestamp before Unix epoch".to_string(),
        ));
    }

    let total_seconds = (jd - UNIX_EPOCH_JD) * SECONDS_PER_DAY_F64;
    let seconds = total_seconds.trunc() as i64;
    let nanos = ((total_seconds - seconds as f64) * 1e9) as u32;

    Ok((seconds, nanos.min(999_999_999)))
}

#[derive(Debug)]
pub struct SerFrame<'a> {
    pub id: SerFrameId,
    pub data: &'a [u8],
    pub timestamp: Option<SerTimestamp>,
}

impl<'a> SerFrame<'a> {
    pub fn new(id: SerFrameId, data: &'a [u8], timestamp: Option<SerTimestamp>) -> Self {
        Self {
            id,
            data,
            timestamp,
        }
    }
}

pub struct SerFile;

impl SerFile {
    pub const HEADER_SIZE: usize = 178;
    pub const TRAILER_ENTRY_SIZE: usize = 8;
    pub const FILE_ID: &'static str = "LUCAM-RECORDER";

    pub fn calculate_bytes_per_pixel(pixel_depth: u32, planes: u32) -> u32 {
        let bytes_per_plane = if pixel_depth <= 8 { 1 } else { 2 };
        bytes_per_plane * planes
    }

    pub fn calculate_frame_size(width: u32, height: u32, bytes_per_pixel: u32) -> u64 {
        width as u64 * height as u64 * bytes_per_pixel as u64
    }

    pub fn calculate_trailer_offset(
        frame_count: u32,
        width: u32,
        height: u32,
        bytes_per_pixel: u32,
    ) -> u64 {
        let header_size = Self::HEADER_SIZE as u64;
        let data_size =
            frame_count as u64 * Self::calculate_frame_size(width, height, bytes_per_pixel);
        header_size + data_size
    }

    pub fn validate_dimensions(width: u32, height: u32) -> Result<(), SerError> {
        if width == 0 || height == 0 || width > 65535 || height > 65535 {
            Err(SerError::InvalidDimensions { width, height })
        } else {
            Ok(())
        }
    }

    pub fn validate_pixel_depth(depth: u32) -> Result<(), SerError> {
        if depth == 0 || depth > 16 {
            Err(SerError::InvalidPixelDepth(depth))
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use astro_time::UTC;

    #[test]
    fn color_id_from_u32_all_variants() {
        assert_eq!(ColorId::from_u32(0).unwrap(), ColorId::Mono);
        assert_eq!(ColorId::from_u32(8).unwrap(), ColorId::BayerRggb);
        assert_eq!(ColorId::from_u32(9).unwrap(), ColorId::BayerGrbg);
        assert_eq!(ColorId::from_u32(10).unwrap(), ColorId::BayerGbrg);
        assert_eq!(ColorId::from_u32(11).unwrap(), ColorId::BayerBggr);
        assert_eq!(ColorId::from_u32(16).unwrap(), ColorId::BayerCyym);
        assert_eq!(ColorId::from_u32(17).unwrap(), ColorId::BayerYcmy);
        assert_eq!(ColorId::from_u32(18).unwrap(), ColorId::BayerYmcy);
        assert_eq!(ColorId::from_u32(19).unwrap(), ColorId::BayerMyyc);
        assert_eq!(ColorId::from_u32(100).unwrap(), ColorId::Rgb);
        assert_eq!(ColorId::from_u32(101).unwrap(), ColorId::Bgr);

        assert!(matches!(
            ColorId::from_u32(999),
            Err(SerError::UnsupportedColorFormat(999))
        ));
    }

    #[test]
    fn color_id_planes() {
        assert_eq!(ColorId::Rgb.planes(), 3);
        assert_eq!(ColorId::Bgr.planes(), 3);
        assert_eq!(ColorId::Mono.planes(), 1);
        assert_eq!(ColorId::BayerRggb.planes(), 1);
        assert_eq!(ColorId::BayerGrbg.planes(), 1);
        assert_eq!(ColorId::BayerGbrg.planes(), 1);
        assert_eq!(ColorId::BayerBggr.planes(), 1);
        assert_eq!(ColorId::BayerCyym.planes(), 1);
        assert_eq!(ColorId::BayerYcmy.planes(), 1);
        assert_eq!(ColorId::BayerYmcy.planes(), 1);
        assert_eq!(ColorId::BayerMyyc.planes(), 1);
    }

    #[test]
    fn color_id_is_bayer() {
        assert!(!ColorId::Mono.is_bayer());
        assert!(!ColorId::Rgb.is_bayer());
        assert!(!ColorId::Bgr.is_bayer());

        assert!(ColorId::BayerRggb.is_bayer());
        assert!(ColorId::BayerGrbg.is_bayer());
        assert!(ColorId::BayerGbrg.is_bayer());
        assert!(ColorId::BayerBggr.is_bayer());
        assert!(ColorId::BayerCyym.is_bayer());
        assert!(ColorId::BayerYcmy.is_bayer());
        assert!(ColorId::BayerYmcy.is_bayer());
        assert!(ColorId::BayerMyyc.is_bayer());
    }

    #[test]
    fn color_id_display_all_variants() {
        assert_eq!(ColorId::Mono.to_string(), "MONO");
        assert_eq!(ColorId::BayerRggb.to_string(), "BAYER_RGGB");
        assert_eq!(ColorId::BayerGrbg.to_string(), "BAYER_GRBG");
        assert_eq!(ColorId::BayerGbrg.to_string(), "BAYER_GBRG");
        assert_eq!(ColorId::BayerBggr.to_string(), "BAYER_BGGR");
        assert_eq!(ColorId::BayerCyym.to_string(), "BAYER_CYYM");
        assert_eq!(ColorId::BayerYcmy.to_string(), "BAYER_YCMY");
        assert_eq!(ColorId::BayerYmcy.to_string(), "BAYER_YMCY");
        assert_eq!(ColorId::BayerMyyc.to_string(), "BAYER_MYYC");
        assert_eq!(ColorId::Rgb.to_string(), "RGB");
        assert_eq!(ColorId::Bgr.to_string(), "BGR");
    }

    #[test]
    fn ser_timestamp_new() {
        // Lines 102-104: Test timestamp creation with ticks conversion
        let ts = SerTimestamp::new(1000000);
        assert_eq!(ts.ticks, 1000000);
        assert!(ts.utc.is_none());
    }

    #[test]
    fn ser_timestamp_from_bytes() {
        let bytes = [0xEF, 0xCD, 0xAB, 0x90, 0x78, 0x56, 0x34, 0x12]; // little-endian
        let ts = SerTimestamp::from_bytes(&bytes);
        assert_eq!(ts.ticks, 0x1234567890ABCDEF);
    }

    #[test]
    fn ser_timestamp_to_bytes() {
        let ts = SerTimestamp::new(0x1234567890ABCDEF);
        let mut buf = [0u8; 8];
        ts.to_bytes(&mut buf);
        assert_eq!(buf, [0xEF, 0xCD, 0xAB, 0x90, 0x78, 0x56, 0x34, 0x12]);
    }

    #[test]
    fn ser_timestamp_from_astro_time() {
        let utc = UTC::new(1672531200, 500_000_000);
        let ts = SerTimestamp::from_astro_time(utc).unwrap();
        assert!(ts.ticks > 0);
        assert!(ts.utc.is_some());
    }

    #[test]
    fn ser_timestamp_to_astro_time_cached() {
        let utc = UTC::new(1672531200, 0);
        let ts = SerTimestamp::from_astro_time(utc).unwrap();
        let recovered = ts.to_utc().unwrap();
        assert_eq!(recovered.to_julian_date(), utc.to_julian_date());
    }

    #[test]
    fn ser_timestamp_to_astro_time_compute() {
        let unix_1980_ticks = SerTimestamp::UNIX_EPOCH_TICKS
            + (8 * 365 * 24 * 60 * 60) as u64 * SerTimestamp::TICKS_PER_SECOND; // ~1980
        let mut ts = SerTimestamp::new(unix_1980_ticks);
        ts.utc = None; // Force computation
        let result = ts.to_utc();
        assert!(result.is_ok());
    }

    #[test]
    fn ser_timestamp_to_utc() {
        let utc = UTC::new(1672531200, 0);
        let ts = SerTimestamp::from_astro_time(utc).unwrap();
        let utc_result = ts.to_utc();
        assert!(utc_result.is_ok());
    }

    #[test]
    fn ser_timestamp_to_unix_seconds_valid() {
        let unix_time = 1672531200.5; // 2023-01-01 00:00:00.5 UTC
        let ts = SerTimestamp::from_unix_seconds(unix_time).unwrap();
        let recovered = ts.to_unix_seconds();
        assert!((recovered - unix_time).abs() < 1e-6);
    }

    #[test]
    fn ser_timestamp_to_unix_seconds_before_epoch() {
        let ts = SerTimestamp::new(100); // Way before unix epoch
        let unix_seconds = ts.to_unix_seconds();
        assert_eq!(unix_seconds, 0.0);
    }

    #[test]
    fn ser_timestamp_from_unix_seconds() {
        let unix_time = 1_672_531_200.123_456_7;
        let ts = SerTimestamp::from_unix_seconds(unix_time).unwrap();
        assert!(ts.ticks > SerTimestamp::UNIX_EPOCH_TICKS);
        assert!(ts.utc.is_some());
        let recovered = ts.to_unix_seconds();
        assert!((recovered - unix_time).abs() < 1e-6);
    }

    #[test]
    fn ser_timestamp_invalid_unix_seconds() {
        assert!(SerTimestamp::from_unix_seconds(-1.0).is_err());
        assert!(SerTimestamp::from_unix_seconds(f64::INFINITY).is_err());
        assert!(SerTimestamp::from_unix_seconds(f64::NAN).is_err());

        let valid = SerTimestamp::from_unix_seconds(1672531200.0).unwrap();
        assert!(valid.ticks > SerTimestamp::UNIX_EPOCH_TICKS);
    }

    #[test]
    fn ser_timestamp_ticks_to_astro_time_valid() {
        let unix_1980_ticks = SerTimestamp::UNIX_EPOCH_TICKS
            + (8 * 365 * 24 * 60 * 60) as u64 * SerTimestamp::TICKS_PER_SECOND;
        let utc = SerTimestamp::ticks_to_utc(unix_1980_ticks).unwrap();
        let (unix_secs, _) = utc_to_unix_timestamp(&utc).unwrap();
        assert!(unix_secs > 0);
    }

    #[test]
    fn ser_timestamp_ticks_to_astro_time_before_epoch() {
        let ticks = 1000; // Before unix epoch
        let result = SerTimestamp::ticks_to_utc(ticks);
        assert!(result.is_err());
    }

    #[test]
    fn ser_timestamp_ticks_conversion_edge_cases() {
        let base_1980 = SerTimestamp::UNIX_EPOCH_TICKS
            + (8 * 365 * 24 * 60 * 60) as u64 * SerTimestamp::TICKS_PER_SECOND;
        let ticks = base_1980 + 15_000_000; // Add 1.5 seconds
        let utc = SerTimestamp::ticks_to_utc(ticks).unwrap();
        let (unix_secs, unix_nanos) = utc_to_unix_timestamp(&utc).unwrap();
        assert!(unix_secs > 0);
        // Julian Date conversion has limited precision (~milliseconds)
        // so we just check the nanos are in the right ballpark
        assert!(unix_nanos > 400_000_000 && unix_nanos < 600_000_000);
    }

    #[test]
    fn ser_timestamp_precision_constant() {
        assert_eq!(SerTimestamp::precision_100ns(), 1);
    }

    #[test]
    fn ser_frame_new() {
        let data = [1, 2, 3, 4, 5];
        let ts = SerTimestamp::new(1000);
        let frame = SerFrame::new(42, &data, Some(ts));
        assert_eq!(frame.id, 42);
        assert_eq!(frame.data, &data);
        assert!(frame.timestamp.is_some());
    }

    #[test]
    fn ser_file_constants() {
        // Test all constants
        assert_eq!(SerFile::HEADER_SIZE, 178);
        assert_eq!(SerFile::TRAILER_ENTRY_SIZE, 8);
        assert_eq!(SerFile::FILE_ID, "LUCAM-RECORDER");
    }

    #[test]
    fn ser_file_calculate_bytes_per_pixel() {
        assert_eq!(SerFile::calculate_bytes_per_pixel(8, 1), 1);
        assert_eq!(SerFile::calculate_bytes_per_pixel(16, 1), 2);
        assert_eq!(SerFile::calculate_bytes_per_pixel(8, 3), 3);
        assert_eq!(SerFile::calculate_bytes_per_pixel(16, 3), 6);
    }

    #[test]
    fn ser_file_calculate_frame_size() {
        assert_eq!(SerFile::calculate_frame_size(100, 100, 2), 20000);
        assert_eq!(SerFile::calculate_frame_size(640, 480, 3), 921600);
    }

    #[test]
    fn ser_file_calculate_trailer_offset() {
        let offset = SerFile::calculate_trailer_offset(10, 100, 100, 2);
        let expected = 178 + (10 * 100 * 100 * 2); // header + data
        assert_eq!(offset, expected as u64);
    }

    #[test]
    fn ser_file_validate_dimensions_valid() {
        assert!(SerFile::validate_dimensions(640, 480).is_ok());
        assert!(SerFile::validate_dimensions(1, 1).is_ok());
        assert!(SerFile::validate_dimensions(65535, 65535).is_ok());
    }

    #[test]
    fn ser_file_validate_dimensions_invalid() {
        assert!(SerFile::validate_dimensions(0, 100).is_err());
        assert!(SerFile::validate_dimensions(100, 0).is_err());
        assert!(SerFile::validate_dimensions(65536, 100).is_err());
        assert!(SerFile::validate_dimensions(100, 65536).is_err());
    }

    #[test]
    fn ser_file_validate_pixel_depth_valid() {
        assert!(SerFile::validate_pixel_depth(1).is_ok());
        assert!(SerFile::validate_pixel_depth(8).is_ok());
        assert!(SerFile::validate_pixel_depth(16).is_ok());
    }

    #[test]
    fn ser_file_validate_pixel_depth_invalid() {
        assert!(SerFile::validate_pixel_depth(0).is_err());
        assert!(SerFile::validate_pixel_depth(17).is_err());
        assert!(SerFile::validate_pixel_depth(32).is_err());
    }

    #[test]
    fn ser_timestamp_constants() {
        assert_eq!(SerTimestamp::TICKS_PER_SECOND, 10_000_000);
        assert_eq!(SerTimestamp::UNIX_EPOCH_TICKS, 621_355_968_000_000_000);
        assert_eq!(SerTimestamp::DOTNET_EPOCH_TICKS, 0);
    }

    #[test]
    fn color_id_and_clone() {
        let color = ColorId::BayerRggb;
        let cloned = color;
        assert_eq!(color, cloned);

        let debug_str = format!("{:?}", color);
        assert!(debug_str.contains("BayerRggb"));
    }

    #[test]
    fn ser_timestamp_and_clone() {
        let ts = SerTimestamp::new(1000);
        let cloned = ts.clone();
        assert_eq!(ts.ticks, cloned.ticks);

        let debug_str = format!("{:?}", ts);
        assert!(debug_str.contains("SerTimestamp"));
    }
}
