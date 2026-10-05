use std::fs;
use std::path::Path;

use reqwest::blocking::Client;

use crate::http::download_file;

const BASE_URL: &str = "http://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02";

pub(crate) const MAIN_FILES: &[&str] = &["ELP_MAIN.S1", "ELP_MAIN.S2", "ELP_MAIN.S3"];

pub(crate) const PERT_FILES: &[&str] = &["ELP_PERT.S1", "ELP_PERT.S2", "ELP_PERT.S3"];

pub(crate) fn file_url(filename: &str) -> String {
    format!("{}/{}", BASE_URL, filename)
}

pub(crate) fn default_client() -> Result<Client, String> {
    Client::builder()
        .timeout(std::time::Duration::from_secs(120))
        .danger_accept_invalid_certs(true)
        .build()
        .map_err(|e| format!("Failed to create HTTP client: {}", e))
}

pub(crate) fn download_all(client: &Client, output_dir: &Path) -> Result<(), String> {
    fs::create_dir_all(output_dir)
        .map_err(|e| format!("Failed to create output directory: {}", e))?;

    println!("Downloading ELP/MPP02 files to {}", output_dir.display());

    for filename in MAIN_FILES.iter().chain(PERT_FILES.iter()) {
        download_file(client, &file_url(filename), filename, output_dir)?;
    }

    println!("Download complete!");
    Ok(())
}

pub(crate) fn find_elp_files(input_dir: &Path) -> Option<ElpFilePaths> {
    let present = |name: &str| Some(input_dir.join(name)).filter(|path| path.exists());
    Some(ElpFilePaths {
        main_longitude: present(MAIN_FILES[0])?,
        main_latitude: present(MAIN_FILES[1])?,
        main_distance: present(MAIN_FILES[2])?,
        pert_longitude: present(PERT_FILES[0])?,
        pert_latitude: present(PERT_FILES[1])?,
        pert_distance: present(PERT_FILES[2])?,
    })
}

#[derive(Debug, Clone)]
pub(crate) struct ElpFilePaths {
    pub(crate) main_longitude: std::path::PathBuf,
    pub(crate) main_latitude: std::path::PathBuf,
    pub(crate) main_distance: std::path::PathBuf,
    pub(crate) pert_longitude: std::path::PathBuf,
    pub(crate) pert_latitude: std::path::PathBuf,
    pub(crate) pert_distance: std::path::PathBuf,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use tempfile::TempDir;

    fn all_files() -> Vec<&'static str> {
        MAIN_FILES.iter().chain(PERT_FILES).copied().collect()
    }

    fn dir_with(names: &[&str]) -> TempDir {
        let dir = TempDir::new().unwrap();
        for name in names {
            fs::write(dir.path().join(name), "").unwrap();
        }
        dir
    }

    #[test]
    fn names_and_urls() {
        assert_eq!(MAIN_FILES, ["ELP_MAIN.S1", "ELP_MAIN.S2", "ELP_MAIN.S3"]);
        assert_eq!(PERT_FILES, ["ELP_PERT.S1", "ELP_PERT.S2", "ELP_PERT.S3"]);
        assert_eq!(
            file_url("ELP_PERT.S2"),
            "http://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/ELP_PERT.S2"
        );
    }

    #[test]
    fn finds_all_six_files() {
        let dir = dir_with(&all_files());
        let paths = find_elp_files(dir.path()).unwrap();
        let found = [
            paths.main_longitude,
            paths.main_latitude,
            paths.main_distance,
            paths.pert_longitude,
            paths.pert_latitude,
            paths.pert_distance,
        ];
        let expected: Vec<PathBuf> = all_files().iter().map(|n| dir.path().join(n)).collect();
        assert_eq!(found.to_vec(), expected);
    }

    #[test]
    fn any_missing_file_means_none() {
        let empty = dir_with(&[]);
        assert!(find_elp_files(empty.path()).is_none());
        assert!(find_elp_files(&empty.path().join("nonexistent")).is_none());
        for missing in all_files() {
            let names: Vec<&str> = all_files().into_iter().filter(|&n| n != missing).collect();
            let dir = dir_with(&names);
            assert!(find_elp_files(dir.path()).is_none(), "{missing}");
        }
    }

    #[test]
    fn client_builds() {
        assert!(default_client().is_ok());
    }
}
