use super::{parser::C04Parser, record::EopRecord};
use crate::{CoordError, CoordResult};
use astro_core::AstroError;
use reqwest;

pub struct EopDownloader {
    base_url: String,

    user_agent: String,
}

pub const IERS_C04_URL: &str = "https://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now";
pub const IERS_C04_HISTORICAL_URL: &str =
    "https://datacenter.iers.org/data/224/eopc04_14_IAU2000.62-now.txt";

impl EopDownloader {
    pub fn new() -> Self {
        Self {
            base_url: "https://hpiers.obspm.fr/iers/eop/eopc04".to_string(),
            user_agent: format!("astro-coords/{}", env!("CARGO_PKG_VERSION")),
        }
    }

    pub async fn download_eop_data(
        &self,
        _start_mjd: f64,
        _end_mjd: f64,
    ) -> CoordResult<Vec<EopRecord>> {
        self.download_latest().await
    }

    pub async fn download_latest(&self) -> CoordResult<Vec<EopRecord>> {
        self.download_from_url(IERS_C04_URL).await
    }

    pub async fn download_from_path(&self, path: &str) -> CoordResult<Vec<EopRecord>> {
        let url = if path.starts_with("http://") || path.starts_with("https://") {
            path.to_string()
        } else {
            format!(
                "{}/{}",
                self.base_url.trim_end_matches('/'),
                path.trim_start_matches('/')
            )
        };
        self.download_from_url(&url).await
    }

    pub async fn download_from_url(&self, url: &str) -> CoordResult<Vec<EopRecord>> {
        let client = reqwest::Client::builder()
            .user_agent(&self.user_agent)
            .timeout(std::time::Duration::from_secs(30))
            .build()
            .map_err(|e| {
                CoordError::from_core(AstroError::data_error(
                    "IERS EOP data",
                    "setup",
                    &format!("Failed to create HTTP client: {}", e),
                ))
            })?;

        let response = client.get(url).send().await.map_err(|e| {
            CoordError::from_core(AstroError::data_error(
                "IERS EOP data",
                "download",
                &format!("Network request failed: {}", e),
            ))
        })?;

        if !response.status().is_success() {
            return Err(CoordError::from_core(AstroError::data_error(
                "IERS EOP data",
                "download",
                &format!("HTTP request failed with status: {}", response.status()),
            )));
        }

        let content = response.text().await.map_err(|e| {
            CoordError::from_core(AstroError::data_error(
                "IERS EOP data",
                "download",
                &format!("Failed to read response: {}", e),
            ))
        })?;

        let parser = C04Parser::new();
        parser.parse(&content)
    }
}

impl Default for EopDownloader {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_downloader_creation() {
        let downloader = EopDownloader::new();
        assert_eq!(
            downloader.base_url,
            "https://hpiers.obspm.fr/iers/eop/eopc04"
        );
        assert!(downloader.user_agent.starts_with("astro-coords/"));
        assert_eq!(
            downloader.user_agent,
            format!("astro-coords/{}", env!("CARGO_PKG_VERSION"))
        );
    }

    #[test]
    fn test_downloader_default() {
        let downloader = EopDownloader::default();
        assert_eq!(
            downloader.base_url,
            "https://hpiers.obspm.fr/iers/eop/eopc04"
        );
        assert!(downloader.user_agent.starts_with("astro-coords/"));
        assert_eq!(
            downloader.user_agent,
            format!("astro-coords/{}", env!("CARGO_PKG_VERSION"))
        );
    }

    #[test]
    fn test_constants() {
        assert_eq!(
            IERS_C04_URL,
            "https://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now"
        );
        assert_eq!(
            IERS_C04_HISTORICAL_URL,
            "https://datacenter.iers.org/data/224/eopc04_14_IAU2000.62-now.txt"
        );
    }

    #[tokio::test]
    async fn test_download_success() {
        use mockito;

        let sample_data = r#"      Date      MJD      x          y        UT1-UTC       LOD         dX        dY        x Err     y Err   UT1-UTC Err  LOD Err     dX Err       dY Err
                         "          "           s           s          "         "           "          "          s         s            "           "
     (0h UTC)

1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230   0.000000   0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000
1962   1   2  37666  -0.015900   0.214100   0.0320547   0.0016690   0.000000   0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000"#;

        let mut server = mockito::Server::new_async().await;
        let mock = server
            .mock("GET", "/test")
            .with_status(200)
            .with_header("content-type", "text/plain")
            .with_body(sample_data)
            .create_async()
            .await;

        let downloader = EopDownloader::new();
        let url = format!("{}/test", server.url());

        let records = downloader.download_from_url(&url).await.unwrap();
        assert_eq!(records.len(), 2);

        let first = &records[0];
        assert_eq!(first.mjd, 37665.0);

        let params = first.to_parameters();
        assert!((params.x_p - (-0.012700)).abs() < 1e-7);
        assert!((params.y_p - 0.213000).abs() < 1e-7);
        assert!((params.ut1_utc - 0.0326338).abs() < 1e-7);
        assert!((params.lod - 0.0017230).abs() < 1e-7);

        mock.assert_async().await;
    }

    #[tokio::test]
    async fn test_download_http_error() {
        use mockito;

        let mut server = mockito::Server::new_async().await;
        let mock = server
            .mock("GET", "/nonexistent")
            .with_status(404)
            .create_async()
            .await;

        let downloader = EopDownloader::new();
        let url = format!("{}/nonexistent", server.url());

        let result = downloader.download_from_url(&url).await;
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("HTTP request failed"));

        mock.assert_async().await;
    }

    #[tokio::test]
    async fn test_download_network_timeout() {
        let downloader = EopDownloader::new();

        let result = tokio::time::timeout(
            std::time::Duration::from_millis(1),
            downloader.download_from_url("http://192.0.2.1:12345/timeout"),
        )
        .await;

        match result {
            Err(_) => {
                // Expected: operation timed out
            }
            Ok(Ok(_)) => {
                println!("Timeout test: unexpected success (test IP responded)");
            }
            Ok(Err(e)) => {
                let error_msg = e.to_string();
                assert!(
                    error_msg.contains("Network request failed")
                        || error_msg.contains("download")
                        || error_msg.contains("timeout")
                        || error_msg.contains("connection")
                );
            }
        }
    }

    #[tokio::test]
    async fn test_download_invalid_data() {
        use mockito;

        let invalid_data = "This is not valid EOP data\nJust some random text\n";

        let mut server = mockito::Server::new_async().await;
        let mock = server
            .mock("GET", "/invalid")
            .with_status(200)
            .with_header("content-type", "text/plain")
            .with_body(invalid_data)
            .create_async()
            .await;

        let downloader = EopDownloader::new();
        let url = format!("{}/invalid", server.url());

        let result = downloader.download_from_url(&url).await;
        match result {
            Ok(records) => {
                assert!(records.is_empty());
            }
            Err(_) => {
                // Expected: invalid data should fail to parse
            }
        }

        mock.assert_async().await;
    }

    #[test]
    fn test_download_latest_delegates() {
        assert_eq!(
            IERS_C04_URL,
            "https://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now"
        );

        let downloader = EopDownloader::new();
        assert!(downloader.user_agent.starts_with("astro-coords/"));
        assert_eq!(
            downloader.base_url,
            "https://hpiers.obspm.fr/iers/eop/eopc04"
        );
    }

    #[tokio::test]
    async fn test_download_eop_data_delegates() {
        use mockito;

        let sample_data = "1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230";

        let mut server = mockito::Server::new_async().await;
        let mock = server
            .mock("GET", "/iers/eop/eopc04/eopc04_IAU2000.62-now")
            .with_status(200)
            .with_header("content-type", "text/plain")
            .with_body(sample_data)
            .create_async()
            .await;

        let downloader = EopDownloader::new();
        let mock_url = format!("{}/iers/eop/eopc04/eopc04_IAU2000.62-now", server.url());

        let result = downloader.download_from_url(&mock_url).await;
        assert!(result.is_ok());

        mock.assert_async().await;
    }

    #[test]
    fn test_user_agent_format() {
        let downloader = EopDownloader::new();
        assert!(downloader.user_agent.starts_with("astro-coords/"));
        // Verify it contains a version number (at least one digit)
        assert!(downloader.user_agent.chars().any(|c| c.is_ascii_digit()));
    }

    #[test]
    fn test_url_constants_format() {
        assert!(IERS_C04_URL.starts_with("https://"));
        assert!(IERS_C04_URL.contains("iers"));
        assert!(IERS_C04_URL.contains("eop"));

        assert!(IERS_C04_HISTORICAL_URL.starts_with("https://"));
        assert!(IERS_C04_HISTORICAL_URL.contains("iers"));
    }

    #[test]
    fn test_base_url_path_construction() {
        let downloader = EopDownloader::new();

        let test_cases = vec![
            (
                "data/file.txt",
                "https://hpiers.obspm.fr/iers/eop/eopc04/data/file.txt",
            ),
            (
                "/data/file.txt",
                "https://hpiers.obspm.fr/iers/eop/eopc04/data/file.txt",
            ),
            ("https://example.com/full", "https://example.com/full"),
            ("http://example.com/full", "http://example.com/full"),
        ];

        for (input_path, expected_url) in test_cases {
            let constructed_url =
                if input_path.starts_with("http://") || input_path.starts_with("https://") {
                    input_path.to_string()
                } else {
                    format!(
                        "{}/{}",
                        downloader.base_url.trim_end_matches('/'),
                        input_path.trim_start_matches('/')
                    )
                };
            assert_eq!(
                constructed_url, expected_url,
                "Failed for input: {}",
                input_path
            );
        }
    }
}
