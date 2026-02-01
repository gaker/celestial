#[derive(Debug, thiserror::Error)]
pub enum XisfError {
    #[error("Invalid XISF format: {0}")]
    InvalidFormat(String),

    #[error("XML parsing error: {0}")]
    XmlParse(String),

    #[error("Missing required element: {0}")]
    MissingElement(String),

    #[error("Invalid geometry: {0}")]
    InvalidGeometry(String),

    #[error("Unsupported sample format: {0}")]
    UnsupportedFormat(String),

    #[error("Data not found at specified location")]
    DataNotFound,

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, XisfError>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Error, ErrorKind};

    #[test]
    fn invalid_format_error_display() {
        let error = XisfError::InvalidFormat("Bad signature".to_string());
        assert_eq!(error.to_string(), "Invalid XISF format: Bad signature");
    }

    #[test]
    fn xml_parse_error_display() {
        let error = XisfError::XmlParse("Unclosed tag".to_string());
        assert_eq!(error.to_string(), "XML parsing error: Unclosed tag");
    }

    #[test]
    fn missing_element_error_display() {
        let error = XisfError::MissingElement("Image".to_string());
        assert_eq!(error.to_string(), "Missing required element: Image");
    }

    #[test]
    fn invalid_geometry_error_display() {
        let error = XisfError::InvalidGeometry("Not a number".to_string());
        assert_eq!(error.to_string(), "Invalid geometry: Not a number");
    }

    #[test]
    fn unsupported_format_error_display() {
        let error = XisfError::UnsupportedFormat("Complex64".to_string());
        assert_eq!(error.to_string(), "Unsupported sample format: Complex64");
    }

    #[test]
    fn data_not_found_error_display() {
        let error = XisfError::DataNotFound;
        assert_eq!(error.to_string(), "Data not found at specified location");
    }

    #[test]
    fn io_error_conversion() {
        let io_error = Error::new(ErrorKind::PermissionDenied, "Access denied");
        let xisf_error: XisfError = io_error.into();

        assert!(matches!(xisf_error, XisfError::Io(_)));
        assert!(xisf_error.to_string().contains("Access denied"));
    }

    #[test]
    fn result_type_alias() {
        let success_result: Result<i32> = Ok(42);
        let error_result: Result<i32> = Err(XisfError::DataNotFound);

        assert!(success_result.is_ok());
        if let Ok(value) = success_result {
            assert_eq!(value, 42);
        }
        assert!(error_result.is_err());
    }

    #[test]
    fn error_is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<XisfError>();
    }

    #[test]
    fn error_chain_with_io_error() {
        let io_error = Error::new(ErrorKind::UnexpectedEof, "Truncated file");
        let xisf_error = XisfError::Io(io_error);

        let error_chain = xisf_error.to_string();
        assert!(error_chain.contains("I/O error"));
        assert!(error_chain.contains("Truncated file"));
    }

    #[test]
    fn long_error_messages() {
        let long_message = "A".repeat(10000);
        let error = XisfError::InvalidFormat(long_message.clone());
        assert!(error.to_string().contains(&long_message));
    }

    #[test]
    fn empty_error_messages() {
        let error = XisfError::InvalidFormat(String::new());
        assert_eq!(error.to_string(), "Invalid XISF format: ");
    }

    #[test]
    fn unicode_in_error_messages() {
        let unicode_message = "файл не найден 🚫";
        let error = XisfError::MissingElement(unicode_message.to_string());
        assert!(error.to_string().contains(unicode_message));
    }
}
