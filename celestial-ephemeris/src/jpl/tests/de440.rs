use super::fixture::*;
use crate::jpl::bodies;
use crate::jpl::spk::SpkFile;
use celestial_core::constants::{AU_KM, J2000_JD};
use std::path::PathBuf;

fn de440_path() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("DE440_PATH") {
        return Some(PathBuf::from(path)).filter(|p| p.exists());
    }
    let home = std::env::var("HOME").ok()?;
    [
        format!("{}/.local/share/ephemeris/de440.bsp", home),
        format!("{}/ephemeris/de440.bsp", home),
        "/usr/local/share/ephemeris/de440.bsp".to_string(),
    ]
    .into_iter()
    .map(PathBuf::from)
    .find(|p| p.exists())
}

fn de440() -> SpkFile {
    let path = de440_path().expect("DE440 file not found - set DE440_PATH");
    SpkFile::open(&path).expect("Failed to open DE440")
}

fn distance_au(body: i32, jd: f64) -> f64 {
    let (pos, _) = state(&de440(), body, bodies::SOLAR_SYSTEM_BARYCENTER, jd, 0.0).unwrap();
    pos.magnitude() / AU_KM
}

#[test]
#[ignore]
fn de440_opens() {
    assert!(!de440().segments().is_empty());
}

#[test]
#[ignore]
fn de440_earth_is_about_1_au_from_the_barycenter() {
    let au = distance_au(bodies::EARTH_MOON_BARYCENTER, J2000_JD);
    assert!(au > 0.98 && au < 1.02, "{}", au);
}

#[test]
#[ignore]
fn de440_mars_is_about_1_5_au_from_the_barycenter() {
    let au = distance_au(bodies::MARS_BARYCENTER, 2460000.5);
    assert!(au > 1.3 && au < 1.7, "{}", au);
}
