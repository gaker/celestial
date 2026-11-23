#[derive(Debug, thiserror::Error)]
pub enum ImageError {
    #[error("Unsupported format")]
    UnsupportedFormat,

    #[error("Format detection failed: {0}")]
    FormatDetectionFailed(String),

    #[error("FITS error: {0}")]
    Fits(#[from] crate::fits::errors::FitsError),

    #[error("XISF error: {0}")]
    Xisf(#[from] crate::xisf::errors::XisfError),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Pixel conversion: {0}")]
    PixelConversion(String),

    #[error("Invalid dimensions: {0}")]
    InvalidDimensions(String),

    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("PNG error: {0}")]
    Png(String),

    #[error("TIFF error: {0}")]
    Tiff(String),
}

pub type Result<T> = std::result::Result<T, ImageError>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Error, ErrorKind};

    #[test]
    fn unsupported_format_error_display() {
        let error = ImageError::UnsupportedFormat;
        assert_eq!(error.to_string(), "Unsupported format");
    }

    #[test]
    fn format_detection_failed_error_display() {
        let error = ImageError::FormatDetectionFailed("No magic bytes found".to_string());
        assert_eq!(
            error.to_string(),
            "Format detection failed: No magic bytes found"
        );
    }

    #[test]
    fn invalid_dimensions_and_input_display() {
        let error = ImageError::InvalidDimensions("width must be nonzero".to_string());
        assert_eq!(
            error.to_string(),
            "Invalid dimensions: width must be nonzero"
        );
        let error = ImageError::InvalidInput("no finite samples".to_string());
        assert_eq!(error.to_string(), "Invalid input: no finite samples");
    }

    #[test]
    fn io_error_conversion() {
        let io_error = Error::new(ErrorKind::NotFound, "File not found");
        let image_error: ImageError = io_error.into();

        assert!(matches!(image_error, ImageError::Io(_)));
        assert!(image_error.to_string().contains("File not found"));
    }

    #[test]
    fn result_type_alias() {
        let success_result: Result<i32> = Ok(42);
        let error_result: Result<i32> = Err(ImageError::UnsupportedFormat);

        assert!(success_result.is_ok());
        if let Ok(value) = success_result {
            assert_eq!(value, 42);
        }
        assert!(error_result.is_err());
    }

    #[test]
    fn error_is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<ImageError>();
    }

    #[test]
    fn error_chain_with_multiple_levels() {
        let io_error = Error::new(ErrorKind::PermissionDenied, "Access denied");
        let image_error = ImageError::Io(io_error);

        let error_chain = image_error.to_string();
        assert!(error_chain.contains("I/O error"));
        assert!(error_chain.contains("Access denied"));
    }

    #[test]
    fn long_format_detection_message() {
        let long_message = "A".repeat(10000);
        let error = ImageError::FormatDetectionFailed(long_message.clone());
        assert!(error.to_string().contains(&long_message));
    }
}
