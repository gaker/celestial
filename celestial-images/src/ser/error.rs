use celestial_time::TimeError;
use std::io;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SerError {
    #[error("Invalid SER header: {0}")]
    InvalidHeader(String),

    #[error("Invalid file ID, expected 'LUCAM-RECORDER', got {actual:?}")]
    InvalidFileId { actual: String },

    #[error("Unsupported color format: {0}")]
    UnsupportedColorFormat(u32),

    #[error("Invalid pixel depth: {0}")]
    InvalidPixelDepth(u32),

    #[error("Frame index {frame} out of bounds (total frames: {total})")]
    FrameOutOfBounds { frame: u32, total: u32 },

    #[error("Invalid frame dimensions: {width}x{height}")]
    InvalidDimensions { width: u32, height: u32 },

    #[error("Timestamp format error: {0}")]
    TimestampError(String),

    #[error("Time conversion failed: {0}")]
    Time(#[from] TimeError),

    #[error("Buffer size mismatch: expected {expected}, got {actual}")]
    BufferSizeMismatch { expected: usize, actual: usize },

    #[error("Write buffer of {megabytes} MiB is larger than the 1024 MiB limit")]
    BufferTooLarge { megabytes: usize },

    #[error("File truncated: expected {expected} bytes, got {actual}")]
    FileTruncated { expected: u64, actual: u64 },

    #[error(
        "Frame {frame} {}, but frame 0 {}",
        if *frame0_has_timestamp { "has no timestamp" } else { "has a timestamp" },
        if *frame0_has_timestamp { "has one" } else { "does not" }
    )]
    TimestampPresence {
        frame: u32,
        frame0_has_timestamp: bool,
    },

    #[error("A SER file holds at most {} frames", u32::MAX)]
    TooManyFrames,

    #[error("Pixel depth {depth} is stored in 8-bit samples, not 16-bit")]
    NotSixteenBit { depth: u32 },

    #[error("SER writer stopped after an earlier I/O error: {0}")]
    WriterFailed(String),

    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
}

pub type Result<T> = std::result::Result<T, SerError>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Error as IoError, ErrorKind};

    #[test]
    fn invalid_header_display() {
        let err = SerError::InvalidHeader("test message".to_string());
        let display = format!("{}", err);
        assert_eq!(display, "Invalid SER header: test message");
    }

    #[test]
    fn invalid_file_id_display() {
        let err = SerError::InvalidFileId {
            actual: "BADHEADER".to_string(),
        };
        let display = format!("{}", err);
        assert_eq!(
            display,
            "Invalid file ID, expected 'LUCAM-RECORDER', got \"BADHEADER\""
        );
    }

    #[test]
    fn unsupported_color_format_display() {
        let err = SerError::UnsupportedColorFormat(999);
        let display = format!("{}", err);
        assert_eq!(display, "Unsupported color format: 999");
    }

    #[test]
    fn invalid_pixel_depth_display() {
        let err = SerError::InvalidPixelDepth(32);
        let display = format!("{}", err);
        assert_eq!(display, "Invalid pixel depth: 32");
    }

    #[test]
    fn frame_out_of_bounds_display() {
        let err = SerError::FrameOutOfBounds {
            frame: 10,
            total: 5,
        };
        let display = format!("{}", err);
        assert_eq!(display, "Frame index 10 out of bounds (total frames: 5)");
    }

    #[test]
    fn invalid_dimensions_display() {
        let err = SerError::InvalidDimensions {
            width: 0,
            height: 100,
        };
        let display = format!("{}", err);
        assert_eq!(display, "Invalid frame dimensions: 0x100");
    }

    #[test]
    fn timestamp_error_display() {
        let err = SerError::TimestampError("invalid format".to_string());
        let display = format!("{}", err);
        assert_eq!(display, "Timestamp format error: invalid format");
    }

    #[test]
    fn new_variant_displays() {
        let cases = [
            (
                SerError::BufferTooLarge { megabytes: 7 },
                "Write buffer of 7 MiB is larger than the 1024 MiB limit",
            ),
            (
                SerError::WriterFailed("disk full".to_string()),
                "SER writer stopped after an earlier I/O error: disk full",
            ),
            (
                SerError::TimestampPresence {
                    frame: 2,
                    frame0_has_timestamp: true,
                },
                "Frame 2 has no timestamp, but frame 0 has one",
            ),
            (
                SerError::TimestampPresence {
                    frame: 1,
                    frame0_has_timestamp: false,
                },
                "Frame 1 has a timestamp, but frame 0 does not",
            ),
            (
                SerError::TooManyFrames,
                "A SER file holds at most 4294967295 frames",
            ),
            (
                SerError::NotSixteenBit { depth: 8 },
                "Pixel depth 8 is stored in 8-bit samples, not 16-bit",
            ),
            (
                SerError::Time(TimeError::ConversionError("x".to_string())),
                "Time conversion failed: Conversion error: x",
            ),
        ];
        for (err, shown) in cases {
            assert_eq!(err.to_string(), shown);
        }
    }

    #[test]
    fn buffer_size_mismatch_display() {
        let err = SerError::BufferSizeMismatch {
            expected: 1024,
            actual: 512,
        };
        let display = format!("{}", err);
        assert_eq!(display, "Buffer size mismatch: expected 1024, got 512");
    }

    #[test]
    fn file_truncated_display() {
        let err = SerError::FileTruncated {
            expected: 10000,
            actual: 5000,
        };
        let display = format!("{}", err);
        assert_eq!(display, "File truncated: expected 10000 bytes, got 5000");
    }

    #[test]
    fn io_error_conversion() {
        let io_err = IoError::new(ErrorKind::NotFound, "file not found");
        let ser_err = SerError::from(io_err);

        match ser_err {
            SerError::Io(ref inner) => {
                assert_eq!(inner.kind(), ErrorKind::NotFound);
                assert_eq!(inner.to_string(), "file not found");
            }
            _ => panic!("Expected IoError variant"),
        }

        let display = format!("{}", ser_err);
        assert_eq!(display, "I/O error: file not found");
    }

    #[test]
    fn io_error_automatic_conversion() {
        fn test_conversion() -> Result<()> {
            Err(IoError::new(ErrorKind::PermissionDenied, "access denied"))?
        }

        let result = test_conversion();
        match result.unwrap_err() {
            SerError::Io(ref inner) => {
                assert_eq!(inner.kind(), ErrorKind::PermissionDenied);
            }
            _ => panic!("Expected automatic IoError conversion"),
        }
    }

    #[test]
    fn error_format() {
        let err = SerError::InvalidHeader("debug test".to_string());
        let debug = format!("{:?}", err);
        assert!(debug.contains("InvalidHeader"));
        assert!(debug.contains("debug test"));
    }

    #[test]
    fn result_type_alias() {
        let success: Result<i32> = Ok(42);
        assert_eq!(success.ok(), Some(42));

        let failure: Result<i32> = Err(SerError::InvalidPixelDepth(0));
        assert!(failure.is_err());
    }

    #[test]
    fn error_is_send_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<SerError>();
    }

    #[test]
    fn error_source_chain() {
        let io_err = IoError::new(ErrorKind::BrokenPipe, "connection lost");
        let ser_err = SerError::from(io_err);

        let source = std::error::Error::source(&ser_err);
        assert!(source.is_some());
        assert_eq!(source.unwrap().to_string(), "connection lost");
    }
}
