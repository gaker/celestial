mod columns;
pub(crate) mod decode;
pub(crate) mod encode;
mod gzip;
mod quantize;
mod restore;
mod rice;
pub(crate) mod spec;
mod table;
mod tiles;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressionAlgorithm {
    Gzip,
    Gzip2,
    Rice,
    HCompress,
    Plio,
}

impl CompressionAlgorithm {
    pub fn fits_name(&self) -> &'static str {
        match self {
            Self::Gzip => "GZIP_1",
            Self::Gzip2 => "GZIP_2",
            Self::Rice => "RICE_1",
            Self::HCompress => "HCOMPRESS_1",
            Self::Plio => "PLIO_1",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compression_algorithm_fits_name() {
        assert_eq!(CompressionAlgorithm::Gzip.fits_name(), "GZIP_1");
        assert_eq!(CompressionAlgorithm::Gzip2.fits_name(), "GZIP_2");
        assert_eq!(CompressionAlgorithm::Rice.fits_name(), "RICE_1");
        assert_eq!(CompressionAlgorithm::HCompress.fits_name(), "HCOMPRESS_1");
        assert_eq!(CompressionAlgorithm::Plio.fits_name(), "PLIO_1");
    }
}
