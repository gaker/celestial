use std::fs;
use std::path::Path;

use reqwest::blocking::Client;

use crate::http::download_file;

const BASE_URL: &str = "https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution";

pub(crate) fn planet_filename(planet: u8) -> String {
    format!("VSOP2013p{}.dat", planet)
}

pub(crate) fn planet_url(planet: u8) -> String {
    format!("{}/{}", BASE_URL, planet_filename(planet))
}

pub(crate) fn default_client() -> Result<Client, String> {
    Client::builder()
        .timeout(std::time::Duration::from_secs(300))
        .build()
        .map_err(|e| format!("Failed to create HTTP client: {}", e))
}

pub(crate) fn download_planet(
    client: &Client,
    planet: u8,
    output_dir: &Path,
) -> Result<(), String> {
    download_file(
        client,
        &planet_url(planet),
        &planet_filename(planet),
        output_dir,
    )
}

pub(crate) fn download_all(client: &Client, output_dir: &Path) -> Result<(), String> {
    fs::create_dir_all(output_dir)
        .map_err(|e| format!("Failed to create output directory: {}", e))?;

    println!("Downloading VSOP2013 files to {}", output_dir.display());

    for planet in 1..=9 {
        download_planet(client, planet, output_dir)?;
    }

    println!("Download complete!");
    Ok(())
}

pub(crate) fn find_planet_files(input_dir: &Path) -> Vec<(u8, std::path::PathBuf)> {
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
    fn names_and_urls() {
        assert_eq!(planet_filename(1), "VSOP2013p1.dat");
        assert_eq!(planet_filename(9), "VSOP2013p9.dat");
        assert_eq!(
            planet_url(3),
            "https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/VSOP2013p3.dat"
        );
    }

    #[test]
    fn finds_the_planet_files_present() {
        let dir = TempDir::new().unwrap();
        assert_eq!(find_planet_files(dir.path()), []);
        for planet in [5, 1, 3] {
            std::fs::write(dir.path().join(planet_filename(planet)), "").unwrap();
        }
        std::fs::write(dir.path().join("VSOP2013p10.dat"), "").unwrap();
        let expected: Vec<(u8, std::path::PathBuf)> = [1, 3, 5]
            .into_iter()
            .map(|p| (p, dir.path().join(planet_filename(p))))
            .collect();
        assert_eq!(find_planet_files(dir.path()), expected);
    }

    #[test]
    fn client_builds() {
        assert!(default_client().is_ok());
    }
}
