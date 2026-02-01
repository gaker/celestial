use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

#[cfg(feature = "cli")]
use reqwest::blocking::Client;

const BASE_URL: &str = "https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution";

pub fn planet_filename(planet: u8) -> String {
    format!("VSOP2013p{}.dat", planet)
}

pub fn planet_url(planet: u8) -> String {
    format!("{}/{}", BASE_URL, planet_filename(planet))
}

#[cfg(feature = "cli")]
pub fn default_client() -> Result<Client, String> {
    Client::builder()
        .timeout(std::time::Duration::from_secs(300))
        .build()
        .map_err(|e| format!("Failed to create HTTP client: {}", e))
}

#[cfg(feature = "cli")]
pub fn download_planet(planet: u8, output_dir: &Path) -> Result<(), String> {
    let client = default_client()?;
    download_file_with_client(
        &client,
        &planet_url(planet),
        &planet_filename(planet),
        output_dir,
    )
}

#[cfg(feature = "cli")]
#[cfg(test)]
pub fn download_planet_with_base_url(
    client: &Client,
    planet: u8,
    base_url: &str,
    output_dir: &Path,
) -> Result<(), String> {
    let url = format!("{}/{}", base_url, planet_filename(planet));
    download_file_with_client(client, &url, &planet_filename(planet), output_dir)
}

#[cfg(feature = "cli")]
#[cfg(test)]
fn download_file_from_url(url: &str, filename: &str, output_dir: &Path) -> Result<(), String> {
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

    println!("Downloading VSOP2013 files to {}", output_dir.display());

    for planet in 1..=9 {
        download_file_with_client(
            client,
            &planet_url(planet),
            &planet_filename(planet),
            output_dir,
        )?;
    }

    println!("Download complete!");
    Ok(())
}

#[cfg(not(feature = "cli"))]
pub fn download_planet(_planet: u8, _output_dir: &Path) -> Result<(), String> {
    Err("Download requires the 'cli' feature".to_string())
}

#[cfg(not(feature = "cli"))]
pub fn download_all(_client: &(), _output_dir: &Path) -> Result<(), String> {
    Err("Download requires the 'cli' feature".to_string())
}

pub fn find_planet_files(input_dir: &Path) -> Vec<(u8, std::path::PathBuf)> {
    let mut files = Vec::new();
    for planet in 1..=9 {
        let path = input_dir.join(planet_filename(planet));
        if path.exists() {
            files.push((planet, path));
        }
    }
    files
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_planet_filename() {
        assert_eq!(planet_filename(1), "VSOP2013p1.dat");
        assert_eq!(planet_filename(9), "VSOP2013p9.dat");
    }

    #[test]
    fn test_planet_url() {
        let url = planet_url(3);
        assert!(url.contains("VSOP2013p3.dat"));
        assert!(url.starts_with("https://"));
    }

    #[test]
    fn test_base_url_constant() {
        assert!(BASE_URL.starts_with("https://"));
        assert!(BASE_URL.contains("imcce.fr"));
        assert!(BASE_URL.contains("vsop2013"));
    }

    #[test]
    fn test_planet_filename_all_planets() {
        for planet in 1..=9 {
            let filename = planet_filename(planet);
            assert!(filename.starts_with("VSOP2013p"));
            assert!(filename.ends_with(".dat"));
            assert!(filename.contains(&planet.to_string()));
        }
    }

    #[test]
    fn test_planet_url_all_planets() {
        for planet in 1..=9 {
            let url = planet_url(planet);
            assert!(url.starts_with(BASE_URL));
            assert!(url.contains(&planet_filename(planet)));
        }
    }

    #[test]
    fn test_find_planet_files_empty_dir() {
        let temp_dir = TempDir::new().unwrap();
        let files = find_planet_files(temp_dir.path());
        assert!(files.is_empty());
    }

    #[test]
    fn test_find_planet_files_with_some_files() {
        let temp_dir = TempDir::new().unwrap();

        // Create files for planets 1, 3, 5
        for planet in [1, 3, 5] {
            let path = temp_dir.path().join(planet_filename(planet));
            std::fs::write(&path, "test content").unwrap();
        }

        let files = find_planet_files(temp_dir.path());
        assert_eq!(files.len(), 3);

        let planets: Vec<u8> = files.iter().map(|(p, _)| *p).collect();
        assert!(planets.contains(&1));
        assert!(planets.contains(&3));
        assert!(planets.contains(&5));
    }

    #[test]
    fn test_find_planet_files_with_all_files() {
        let temp_dir = TempDir::new().unwrap();

        for planet in 1..=9 {
            let path = temp_dir.path().join(planet_filename(planet));
            std::fs::write(&path, "test content").unwrap();
        }

        let files = find_planet_files(temp_dir.path());
        assert_eq!(files.len(), 9);
    }

    #[test]
    fn test_find_planet_files_returns_correct_paths() {
        let temp_dir = TempDir::new().unwrap();
        let path = temp_dir.path().join(planet_filename(3));
        std::fs::write(&path, "test").unwrap();

        let files = find_planet_files(temp_dir.path());
        assert_eq!(files.len(), 1);
        assert_eq!(files[0].0, 3);
        assert_eq!(files[0].1, path);
    }

    #[cfg(not(feature = "cli"))]
    #[test]
    fn test_download_planet_without_cli_feature() {
        let temp_dir = TempDir::new().unwrap();
        let result = download_planet(1, temp_dir.path());
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Download requires the 'cli' feature");
    }

    #[cfg(not(feature = "cli"))]
    #[test]
    fn test_download_all_without_cli_feature() {
        let temp_dir = TempDir::new().unwrap();
        let result = download_all(&(), temp_dir.path());
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Download requires the 'cli' feature");
    }

    #[cfg(feature = "cli")]
    mod mock_http_tests {
        use super::*;
        use std::io::Write as IoWrite;
        use wiremock::matchers::{method, path};
        use wiremock::{Mock, MockServer, ResponseTemplate};

        #[tokio::test]
        async fn test_download_file_success() {
            let mock_server = MockServer::start().await;
            let test_content = b"VSOP2013 test data content";

            Mock::given(method("GET"))
                .and(path("/test.dat"))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(test_content.to_vec()))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/test.dat", mock_server.uri());

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "test.dat", temp_dir.path()).map(|_| temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            let downloaded = std::fs::read(temp_dir.path().join("test.dat")).unwrap();
            assert_eq!(downloaded, test_content);
        }

        #[tokio::test]
        async fn test_download_file_skips_existing() {
            let mock_server = MockServer::start().await;

            Mock::given(method("GET"))
                .and(path("/existing.dat"))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(b"new content".to_vec()))
                .expect(0) // Should NOT be called
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let existing_path = temp_dir.path().join("existing.dat");
            {
                let mut file = File::create(&existing_path).unwrap();
                file.write_all(b"original content").unwrap();
            }

            let url = format!("{}/existing.dat", mock_server.uri());
            let temp_path = temp_dir.path().to_path_buf();

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "existing.dat", &temp_path)
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
                .and(path("/missing.dat"))
                .respond_with(ResponseTemplate::new(404))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/missing.dat", mock_server.uri());
            let temp_path = temp_dir.path().to_path_buf();

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "missing.dat", &temp_path)
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
                .and(path("/error.dat"))
                .respond_with(ResponseTemplate::new(500))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/error.dat", mock_server.uri());
            let temp_path = temp_dir.path().to_path_buf();

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "error.dat", &temp_path)
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
        async fn test_download_file_writes_binary_content() {
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
                .and(path("/empty.dat"))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(vec![]))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/empty.dat", mock_server.uri());

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "empty.dat", temp_dir.path()).map(|_| temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            let downloaded = std::fs::read(temp_dir.path().join("empty.dat")).unwrap();
            assert!(downloaded.is_empty());
        }

        #[tokio::test]
        async fn test_download_file_creates_file() {
            let mock_server = MockServer::start().await;

            Mock::given(method("GET"))
                .and(path("/new.dat"))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(b"data".to_vec()))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let url = format!("{}/new.dat", mock_server.uri());

            assert!(!temp_dir.path().join("new.dat").exists());

            let result = tokio::task::spawn_blocking(move || {
                download_file_from_url(&url, "new.dat", temp_dir.path()).map(|_| temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            assert!(temp_dir.path().join("new.dat").exists());
        }

        #[tokio::test]
        async fn test_download_all_success() {
            let mock_server = MockServer::start().await;

            for planet in 1..=9 {
                let filename = planet_filename(planet);
                let content = format!("mock content for planet {}", planet);
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

                for planet in 1..=9 {
                    let filename = planet_filename(planet);
                    let url = format!("{}/{}", base_url, filename);
                    download_file_with_client(&client, &url, &filename, temp_dir.path())?;
                }
                Ok::<_, String>(temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            for planet in 1..=9 {
                let filename = planet_filename(planet);
                assert!(
                    temp_dir.path().join(&filename).exists(),
                    "Missing file: {}",
                    filename
                );
            }
        }

        #[tokio::test]
        async fn test_download_all_fails_on_http_error() {
            let mock_server = MockServer::start().await;

            Mock::given(method("GET"))
                .and(path("/VSOP2013p1.dat"))
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

                let url = format!("{}/VSOP2013p1.dat", base_url);
                download_file_with_client(&client, &url, "VSOP2013p1.dat", temp_dir.path())
            })
            .await
            .unwrap();

            assert!(result.is_err());
            assert!(result.unwrap_err().contains("HTTP error"));
        }

        #[tokio::test]
        async fn test_download_planet_with_base_url_success() {
            let mock_server = MockServer::start().await;
            let planet = 3u8;
            let filename = planet_filename(planet);
            let content = b"mocked VSOP2013 Earth data";

            Mock::given(method("GET"))
                .and(path(format!("/{}", filename)))
                .respond_with(ResponseTemplate::new(200).set_body_bytes(content.to_vec()))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let base_url = mock_server.uri();

            let result = tokio::task::spawn_blocking(move || {
                let client = Client::builder()
                    .timeout(std::time::Duration::from_secs(10))
                    .build()
                    .unwrap();

                download_planet_with_base_url(&client, planet, &base_url, temp_dir.path())
                    .map(|_| temp_dir)
            })
            .await
            .unwrap();

            let temp_dir = result.unwrap();
            let downloaded = std::fs::read(temp_dir.path().join(&filename)).unwrap();
            assert_eq!(downloaded, content);
        }

        #[tokio::test]
        async fn test_download_planet_with_base_url_http_error() {
            let mock_server = MockServer::start().await;
            let planet = 5u8;
            let filename = planet_filename(planet);

            Mock::given(method("GET"))
                .and(path(format!("/{}", filename)))
                .respond_with(ResponseTemplate::new(404))
                .mount(&mock_server)
                .await;

            let temp_dir = TempDir::new().unwrap();
            let base_url = mock_server.uri();

            let result = tokio::task::spawn_blocking(move || {
                let client = Client::builder()
                    .timeout(std::time::Duration::from_secs(10))
                    .build()
                    .unwrap();

                download_planet_with_base_url(&client, planet, &base_url, temp_dir.path())
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
    }

    #[cfg(feature = "cli")]
    #[test]
    fn test_default_client() {
        let result = default_client();
        assert!(result.is_ok());
    }
}
