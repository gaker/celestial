#![cfg(feature = "integration-tests")]

use celestial_catalog::query::catalog::Catalog;

const TEST_CATALOG: &str = "data/catalog.bin";

#[test]
fn test_catalog_open() {
    let catalog = Catalog::open(TEST_CATALOG).expect("Failed to open catalog");
    let header = catalog.header();

    assert_eq!(header.order, 8);
    assert_eq!(header.nside, 256);
    assert_eq!(header.npix, 786432);
    assert!(header.total_stars > 0);
    assert_eq!(header.epoch, 2016.0);
}

#[test]
fn test_stars_in_pixel() {
    let catalog = Catalog::open(TEST_CATALOG).expect("Failed to open catalog");

    let stars = catalog.stars_in_pixel(100000);
    assert!(!stars.is_empty(), "Expected non-empty pixel");

    for star in stars {
        assert!(star.ra >= 0.0 && star.ra < 360.0, "Invalid RA");
        assert!(star.dec >= -90.0 && star.dec <= 90.0, "Invalid Dec");
        assert!(star.mag >= 0.0 && star.mag < 30.0, "Invalid magnitude");
    }
}

#[test]
fn test_out_of_bounds_pixel() {
    let catalog = Catalog::open(TEST_CATALOG).expect("Failed to open catalog");
    let header = catalog.header();

    let stars = catalog.stars_in_pixel(header.npix + 1000);
    assert_eq!(
        stars.len(),
        0,
        "Out of bounds pixel should return empty slice"
    );
}

#[test]
fn test_magnitude_sorting() {
    let catalog = Catalog::open(TEST_CATALOG).expect("Failed to open catalog");

    let stars = catalog.stars_in_pixel(100000);
    assert!(stars.len() > 1, "Need multiple stars to test sorting");

    for i in 1..stars.len() {
        assert!(
            stars[i].mag >= stars[i - 1].mag,
            "Stars not sorted by magnitude: {} >= {}",
            stars[i].mag,
            stars[i - 1].mag
        );
    }
}
