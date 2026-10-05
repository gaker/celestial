use std::fs;
use std::path::Path;

use reqwest::blocking::Client;
use reqwest::StatusCode;

pub(crate) fn download_file(
    client: &Client,
    url: &str,
    filename: &str,
    output_dir: &Path,
) -> Result<(), String> {
    save_unless_present(&output_dir.join(filename), || fetch(client, url))
}

fn fetch(client: &Client, url: &str) -> Result<Vec<u8>, String> {
    println!("  Downloading {} ...", url);
    let response = client
        .get(url)
        .send()
        .map_err(|e| format!("Failed to fetch {}: {}", url, e))?;
    check_status(response.status(), url)?;
    let bytes = response
        .bytes()
        .map_err(|e| format!("Failed to read response: {}", e))?;
    Ok(bytes.to_vec())
}

fn check_status(status: StatusCode, url: &str) -> Result<(), String> {
    if status.is_success() {
        Ok(())
    } else {
        Err(format!("HTTP error {} for {}", status, url))
    }
}

fn save_unless_present(
    path: &Path,
    fetch: impl FnOnce() -> Result<Vec<u8>, String>,
) -> Result<(), String> {
    if path.exists() {
        println!("  {} already exists, skipping", path.display());
        return Ok(());
    }
    let bytes = fetch()?;
    fs::write(path, &bytes).map_err(|e| format!("Failed to write {}: {}", path.display(), e))?;
    println!("  Saved {} ({} bytes)", path.display(), bytes.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    const URL: &str = "http://example.invalid/file.dat";

    #[test]
    fn success_status_passes() {
        assert_eq!(check_status(StatusCode::OK, URL), Ok(()));
    }

    #[test]
    fn error_status_names_code_and_url() {
        assert_eq!(
            check_status(StatusCode::NOT_FOUND, URL),
            Err(format!("HTTP error 404 Not Found for {URL}"))
        );
        assert_eq!(
            check_status(StatusCode::INTERNAL_SERVER_ERROR, URL),
            Err(format!("HTTP error 500 Internal Server Error for {URL}"))
        );
    }

    #[test]
    fn writes_fetched_bytes() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("binary.dat");
        let bytes: Vec<u8> = (0..=255).collect();
        assert_eq!(save_unless_present(&path, || Ok(bytes.clone())), Ok(()));
        assert_eq!(fs::read(&path).unwrap(), bytes);
    }

    #[test]
    fn writes_empty_response() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("empty.dat");
        assert_eq!(save_unless_present(&path, || Ok(Vec::new())), Ok(()));
        assert_eq!(fs::read(&path).unwrap(), Vec::<u8>::new());
    }

    #[test]
    fn skips_existing_file_without_fetching() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("existing.dat");
        fs::write(&path, b"original").unwrap();
        let fetch = || Err("fetched".to_string());
        assert_eq!(save_unless_present(&path, fetch), Ok(()));
        assert_eq!(fs::read(&path).unwrap(), b"original");
    }

    #[test]
    fn fetch_error_leaves_no_file() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("missing.dat");
        let fetch = || Err("HTTP error 404 Not Found".to_string());
        assert_eq!(
            save_unless_present(&path, fetch),
            Err("HTTP error 404 Not Found".to_string())
        );
        assert!(!path.exists());
    }

    #[test]
    fn write_error_names_path() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("no-such-dir").join("file.dat");
        let err = save_unless_present(&path, || Ok(vec![1])).unwrap_err();
        assert!(
            err.starts_with(&format!("Failed to write {}: ", path.display())),
            "{err}"
        );
    }
}
