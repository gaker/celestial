use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

#[cfg(feature = "cli")]
use reqwest::blocking::Client;

const BASE_URL: &str = "http://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02";

pub const MAIN_FILES: &[&str] = &["ELP_MAIN.S1", "ELP_MAIN.S2", "ELP_MAIN.S3"];

pub const PERT_FILES: &[&str] = &["ELP_PERT.S1", "ELP_PERT.S2", "ELP_PERT.S3"];

#[allow(dead_code)]
pub const FORTRAN_FILE: &str = "ELPMPP02.for";

pub fn file_url(filename: &str) -> String {
    format!("{}/{}", BASE_URL, filename)
}

#[cfg(feature = "cli")]
pub fn default_client() -> Result<Client, String> {
    Client::builder()
        .timeout(std::time::Duration::from_secs(120))
        .danger_accept_invalid_certs(true)
        .build()
        .map_err(|e| format!("Failed to create HTTP client: {}", e))
}

#[cfg(feature = "cli")]
#[allow(dead_code)]
pub fn download_file_from_url(url: &str, filename: &str, output_dir: &Path) -> Result<(), String> {
    let client = default_client()?;
    download_file_with_client(&client, url, filename, output_dir)
}

#[cfg(feature = "cli")]
pub fn download_file_with_client(
    client: &Client,
    url: &str,
    filename: &str,
    output_dir: &Path,
) -> Result<(), String> {
    let output_path = output_dir.join(filename);

    if output_path.exists() {
        println!("  {} already exists, skipping", filename);
        return Ok(());
    }

    println!("  Downloading {} ...", url);

    let response = client
        .get(url)
        .send()
        .map_err(|e| format!("Failed to fetch {}: {}", url, e))?;

    if !response.status().is_success() {
        return Err(format!("HTTP error {} for {}", response.status(), url));
    }

    let bytes = response
        .bytes()
        .map_err(|e| format!("Failed to read response: {}", e))?;

    let mut file = File::create(&output_path)
        .map_err(|e| format!("Failed to create {}: {}", output_path.display(), e))?;

    file.write_all(&bytes)
        .map_err(|e| format!("Failed to write {}: {}", output_path.display(), e))?;

    println!("  Saved {} ({} bytes)", filename, bytes.len());
    Ok(())
}

#[cfg(feature = "cli")]
pub fn download_all(client: &Client, output_dir: &Path) -> Result<(), String> {
    fs::create_dir_all(output_dir)
        .map_err(|e| format!("Failed to create output directory: {}", e))?;

    println!("Downloading ELP/MPP02 files to {}", output_dir.display());

    for filename in MAIN_FILES.iter().chain(PERT_FILES.iter()) {
        download_file_with_client(client, &file_url(filename), filename, output_dir)?;
    }

    println!("Download complete!");
    Ok(())
}

#[cfg(not(feature = "cli"))]
pub fn download_all(_client: &(), _output_dir: &Path) -> Result<(), String> {
    Err("Download requires the 'cli' feature".to_string())
}

pub fn find_elp_files(input_dir: &Path) -> Option<ElpFilePaths> {
    let main_files: Vec<_> = MAIN_FILES.iter().map(|f| input_dir.join(f)).collect();
    let pert_files: Vec<_> = PERT_FILES.iter().map(|f| input_dir.join(f)).collect();

    for path in main_files.iter().chain(pert_files.iter()) {
        if !path.exists() {
            return None;
        }
    }

    Some(ElpFilePaths {
        main_longitude: main_files[0].clone(),
        main_latitude: main_files[1].clone(),
        main_distance: main_files[2].clone(),
        pert_longitude: pert_files[0].clone(),
        pert_latitude: pert_files[1].clone(),
        pert_distance: pert_files[2].clone(),
    })
}

#[derive(Debug, Clone)]
pub struct ElpFilePaths {
    pub main_longitude: std::path::PathBuf,
    pub main_latitude: std::path::PathBuf,
    pub main_distance: std::path::PathBuf,
    pub pert_longitude: std::path::PathBuf,
    pub pert_latitude: std::path::PathBuf,
    pub pert_distance: std::path::PathBuf,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use tempfile::TempDir;

    #[test]
    fn test_file_url() {
        let url = file_url("ELP_MAIN.S1");
        assert!(url.contains("ELP_MAIN.S1"));
        assert!(url.starts_with("http://"));
        assert!(url.contains("cyrano-se.obspm.fr"));
    }

    #[test]
    fn test_file_url_for_all_main_files() {
        for filename in MAIN_FILES {
            let url = file_url(filename);
            assert_eq!(url, format!("{}/{}", BASE_URL, filename));
        }
    }

    #[test]
    fn test_file_url_for_all_pert_files() {
        for filename in PERT_FILES {
            let url = file_url(filename);
            assert_eq!(url, format!("{}/{}", BASE_URL, filename));
        }
    }

    #[test]
    fn test_file_url_for_fortran_file() {
        let url = file_url(FORTRAN_FILE);
        assert_eq!(
            url,
            "http://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/ELPMPP02.for"
        );
    }

    #[test]
    fn test_file_url_empty_filename() {
        let url = file_url("");
        assert_eq!(url, format!("{}/", BASE_URL));
    }

    #[test]
    fn test_file_url_with_special_characters() {
        let url = file_url("file with spaces.txt");
        assert!(url.ends_with("file with spaces.txt"));
    }

    #[test]
    fn test_file_lists() {
        assert_eq!(MAIN_FILES.len(), 3);
        assert_eq!(PERT_FILES.len(), 3);
    }

    #[test]
    fn test_main_files_naming_convention() {
        assert_eq!(MAIN_FILES[0], "ELP_MAIN.S1");
        assert_eq!(MAIN_FILES[1], "ELP_MAIN.S2");
        assert_eq!(MAIN_FILES[2], "ELP_MAIN.S3");
    }

    #[test]
    fn test_pert_files_naming_convention() {
        assert_eq!(PERT_FILES[0], "ELP_PERT.S1");
        assert_eq!(PERT_FILES[1], "ELP_PERT.S2");
        assert_eq!(PERT_FILES[2], "ELP_PERT.S3");
    }

    #[test]
    fn test_fortran_file_constant() {
        assert_eq!(FORTRAN_FILE, "ELPMPP02.for");
    }

    #[test]
    fn test_find_elp_files_returns_none_when_dir_empty() {
        let temp_dir = TempDir::new().unwrap();
        let result = find_elp_files(temp_dir.path());
        assert!(result.is_none());
    }

    #[test]
    fn test_find_elp_files_returns_none_when_main_files_missing() {
        let temp_dir = TempDir::new().unwrap();
        for filename in PERT_FILES {
            File::create(temp_dir.path().join(filename)).unwrap();
        }
        let result = find_elp_files(temp_dir.path());
        assert!(result.is_none());
    }

    #[test]
    fn test_find_elp_files_returns_none_when_pert_files_missing() {
        let temp_dir = TempDir::new().unwrap();
        for filename in MAIN_FILES {
            File::create(temp_dir.path().join(filename)).unwrap();
        }
        let result = find_elp_files(temp_dir.path());
        assert!(result.is_none());
    }

    #[test]
    fn test_find_elp_files_returns_none_when_partial_main_files() {
        let temp_dir = TempDir::new().unwrap();
        File::create(temp_dir.path().join("ELP_MAIN.S1")).unwrap();
        File::create(temp_dir.path().join("ELP_MAIN.S2")).unwrap();
        // Missing ELP_MAIN.S3
        for filename in PERT_FILES {
            File::create(temp_dir.path().join(filename)).unwrap();
        }
        let result = find_elp_files(temp_dir.path());
        assert!(result.is_none());
    }

    #[test]
    fn test_find_elp_files_returns_none_when_partial_pert_files() {
        let temp_dir = TempDir::new().unwrap();
        for filename in MAIN_FILES {
            File::create(temp_dir.path().join(filename)).unwrap();
        }
        File::create(temp_dir.path().join("ELP_PERT.S1")).unwrap();
        File::create(temp_dir.path().join("ELP_PERT.S2")).unwrap();
        // Missing ELP_PERT.S3
        let result = find_elp_files(temp_dir.path());
        assert!(result.is_none());
    }

    #[test]
    fn test_find_elp_files_returns_none_when_first_file_missing() {
        let temp_dir = TempDir::new().unwrap();
        // Missing ELP_MAIN.S1 (first file checked)
        File::create(temp_dir.path().join("ELP_MAIN.S2")).unwrap();
        File::create(temp_dir.path().join("ELP_MAIN.S3")).unwrap();
        for filename in PERT_FILES {
            File::create(temp_dir.path().join(filename)).unwrap();
        }
        let result = find_elp_files(temp_dir.path());
        assert!(result.is_none());
    }

    #[test]
    fn test_find_elp_files_returns_some_when_all_files_present() {
        let temp_dir = TempDir::new().unwrap();
        for filename in MAIN_FILES.iter().chain(PERT_FILES.iter()) {
            File::create(temp_dir.path().join(filename)).unwrap();
        }
        let result = find_elp_files(temp_dir.path());
        assert!(result.is_some());
    }

    #[test]
    fn test_find_elp_files_returns_correct_paths() {
        let temp_dir = TempDir::new().unwrap();
        for filename in MAIN_FILES.iter().chain(PERT_FILES.iter()) {
            File::create(temp_dir.path().join(filename)).unwrap();
        }
        let paths = find_elp_files(temp_dir.path()).unwrap();

        assert_eq!(paths.main_longitude, temp_dir.path().join("ELP_MAIN.S1"));
        assert_eq!(paths.main_latitude, temp_dir.path().join("ELP_MAIN.S2"));
        assert_eq!(paths.main_distance, temp_dir.path().join("ELP_MAIN.S3"));
        assert_eq!(paths.pert_longitude, temp_dir.path().join("ELP_PERT.S1"));
        assert_eq!(paths.pert_latitude, temp_dir.path().join("ELP_PERT.S2"));
        assert_eq!(paths.pert_distance, temp_dir.path().join("ELP_PERT.S3"));
    }

    #[test]
    fn test_find_elp_files_with_nonexistent_directory() {
        let nonexistent = Path::new("/nonexistent/path/that/does/not/exist");
        let result = find_elp_files(nonexistent);
        assert!(result.is_none());
    }

    #[test]
    fn test_elp_file_paths_clone() {
        let temp_dir = TempDir::new().unwrap();
        for filename in MAIN_FILES.iter().chain(PERT_FILES.iter()) {
            File::create(temp_dir.path().join(filename)).unwrap();
        }
        let paths = find_elp_files(temp_dir.path()).unwrap();
        let cloned = paths.clone();

        assert_eq!(paths.main_longitude, cloned.main_longitude);
        assert_eq!(paths.main_latitude, cloned.main_latitude);
        assert_eq!(paths.main_distance, cloned.main_distance);
        assert_eq!(paths.pert_longitude, cloned.pert_longitude);
        assert_eq!(paths.pert_latitude, cloned.pert_latitude);
        assert_eq!(paths.pert_distance, cloned.pert_distance);
    }

    #[test]
    fn test_elp_file_paths_debug() {
        let temp_dir = TempDir::new().unwrap();
        for filename in MAIN_FILES.iter().chain(PERT_FILES.iter()) {
            File::create(temp_dir.path().join(filename)).unwrap();
        }
        let paths = find_elp_files(temp_dir.path()).unwrap();
        let debug_str = format!("{:?}", paths);

        assert!(debug_str.contains("ElpFilePaths"));
        assert!(debug_str.contains("main_longitude"));
        assert!(debug_str.contains("main_latitude"));
        assert!(debug_str.contains("main_distance"));
        assert!(debug_str.contains("pert_longitude"));
        assert!(debug_str.contains("pert_latitude"));
        assert!(debug_str.contains("pert_distance"));
    }

    #[test]
    fn test_base_url_constant() {
        assert_eq!(
            BASE_URL,
            "http://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02"
        );
    }

    #[cfg(not(feature = "cli"))]
    #[test]
    fn test_download_all_without_cli_feature() {
        let temp_dir = TempDir::new().unwrap();
        let result = download_all(&(), temp_dir.path());
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Download requires the 'cli' feature".to_string()
        );
    }

    #[cfg(feature = "cli")]
    mod mock_http_tests {
        use super::*;
        use std::io::Write;
        use wiremock::matchers::{method, path};
        use wiremock::{Mock, MockServer, ResponseTemplate};

        #[tokio::test]
        async fn test_download_file_success() {
            let mock_server = MockServer::start().await;
            let test_content = b"test file content for ELP data";

            Mock::given(method("GET"))
                .and(path("/test.txt"))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(test_content.to_vec()))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/test.txt", mock_server.uri());

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "test.txt", temp_dir.path()).map(|_| temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            let downloaded = std::fs::read(temp_dir.path().join("test.txt")).unwrap();
            assert_eq!(downloaded, test_content);
        }

        #[tokio::test]
        async fn test_download_file_skips_existing() {
            let mock_server = MockServer::start().await;

            Mock::given(method("GET"))
                .and(path("/existing.txt"))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(b"new content".to_vec()))
                .expect(0) // Should NOT be called
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let existing_path = temp_dir.path().join("existing.txt");
            {
                let mut file = File::create(&existing_path).unwrap();
                file.write_all(b"original content").unwrap();
            }

            let url = format!("{}/existing.txt", mock_server.uri());
            let temp_path = temp_dir.path().to_path_buf();

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "existing.txt", &temp_path)
            })
            .await
            .unwrap();

            assert!(result.is_ok());
            let content = std::fs::read_to_string(&existing_path).unwrap();
            assert_eq!(content, "original content");
        }

        #[tokio::test]
        async fn test_download_file_http_404() {
            let mock_server = MockServer::start().await;

            Mock::given(method("GET"))
                .and(path("/missing.txt"))
                .respond_with(ResponseTemplate::new(404))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/missing.txt", mock_server.uri());
            let temp_path = temp_dir.path().to_path_buf();

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "missing.txt", &temp_path)
            })
            .await
            .unwrap();

            assert!(result.is_err());
            let err = result.unwrap_err();
            assert!(
                err.contains("HTTP error"),
                "Expected HTTP error, got: {}",
                err
            );
            assert!(err.contains("404"), "Expected 404 in error, got: {}", err);
        }

        #[tokio::test]
        async fn test_download_file_http_500() {
            let mock_server = MockServer::start().await;

            Mock::given(method("GET"))
                .and(path("/error.txt"))
                .respond_with(ResponseTemplate::new(500))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/error.txt", mock_server.uri());
            let temp_path = temp_dir.path().to_path_buf();

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "error.txt", &temp_path)
            })
            .await
            .unwrap();

            assert!(result.is_err());
            let err = result.unwrap_err();
            assert!(
                err.contains("HTTP error"),
                "Expected HTTP error, got: {}",
                err
            );
            assert!(err.contains("500"), "Expected 500 in error, got: {}", err);
        }

        #[tokio::test]
        async fn test_download_file_writes_correct_bytes() {
            let mock_server = MockServer::start().await;
            let binary_content: Vec<u8> = (0..256).map(|i| i as u8).collect();

            Mock::given(method("GET"))
                .and(path("/binary.dat"))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(binary_content.clone()))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/binary.dat", mock_server.uri());

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "binary.dat", temp_dir.path()).map(|_| temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            let downloaded = std::fs::read(temp_dir.path().join("binary.dat")).unwrap();
            assert_eq!(downloaded, binary_content);
        }

        #[tokio::test]
        async fn test_download_file_empty_response() {
            let mock_server = MockServer::start().await;

            Mock::given(method("GET"))
                .and(path("/empty.txt"))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(vec![]))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/empty.txt", mock_server.uri());

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "empty.txt", temp_dir.path()).map(|_| temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            let downloaded = std::fs::read(temp_dir.path().join("empty.txt")).unwrap();
            assert!(downloaded.is_empty());
        }

        #[tokio::test]
        async fn test_download_file_creates_file_in_output_dir() {
            let mock_server = MockServer::start().await;

            Mock::given(method("GET"))
                .and(path("/data.bin"))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(b"data".to_vec()))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/data.bin", mock_server.uri());

            assert!(!temp_dir.path().join("data.bin").exists());

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "data.bin", temp_dir.path()).map(|_| temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            assert!(temp_dir.path().join("data.bin").exists());
        }

        #[tokio::test]
        async fn test_download_all_success() {
            let mock_server = MockServer::start().await;

            for filename in MAIN_FILES.iter().chain(PERT_FILES.iter()) {
                let content = format!("mock content for {}", filename);
                Mock::given(method("GET"))
                    .and(path(format!("/{}", filename)))
                    .respond_with(ResponseTemplate::new(200).set_body_bytes(content.into_bytes()))
                    .mount(&mock_server)
                    .await;
            }

            let temp_dir = TempDir::new().unwrap();
            let base_url = mock_server.uri();

            let result = tokio::task::spawn_blocking(move || {
                let client = Client::builder()
                    .timeout(std::time::Duration::from_secs(10))
                    .build()
                    .unwrap();

                for filename in MAIN_FILES.iter().chain(PERT_FILES.iter()) {
                    let url = format!("{}/{}", base_url, filename);
                    download_file_with_client(&client, &url, filename, temp_dir.path())?;
                }
                Ok::<_, String>(temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            for filename in MAIN_FILES.iter().chain(PERT_FILES.iter()) {
                assert!(
                    temp_dir.path().join(filename).exists(),
                    "Missing file: {}",
                    filename
                );
            }
        }

        #[tokio::test]
        async fn test_download_all_fails_on_http_error() {
            let mock_server = MockServer::start().await;

            Mock::given(method("GET"))
                .and(path("/ELP_MAIN.S1"))
                .respond_with(ResponseTemplate::new(500))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let base_url = mock_server.uri();

            let result = tokio::task::spawn_blocking(move || {
                let client = Client::builder()
                    .timeout(std::time::Duration::from_secs(10))
                    .build()
                    .unwrap();

                let url = format!("{}/ELP_MAIN.S1", base_url);
                download_file_with_client(&client, &url, "ELP_MAIN.S1", temp_dir.path())
            })
            .await
            .unwrap();

            assert!(result.is_err());
            assert!(result.unwrap_err().contains("HTTP error"));
        }
    }
}
