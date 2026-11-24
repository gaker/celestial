// Metadata, images, a catalog and pairs shared by the solve tests.

use celestial_catalog::query::catalog::Catalog;
use celestial_coords::frames::icrs::ICRSPosition;
use celestial_images::formats::image::Image;
use celestial_images::formats::pixel_data::PixelData;
use celestial_time::julian::JulianDate;

use crate::match_field::StarPair;
use crate::metadata::image_metadata::ImageMetadata;

pub(super) fn meta_at(ra: f64, dec: f64, scale_arcsec: f64) -> ImageMetadata {
    ImageMetadata {
        hint: ICRSPosition::from_degrees(ra, dec).unwrap(),
        scale_arcsec,
        epoch: JulianDate::new(2451545.0, 0.0),
        focal_mm: None,
        pixel_um: None,
    }
}

pub(super) fn zero_image(w: usize, h: usize) -> Image {
    Image::new(PixelData::from(vec![0u16; w * h]), [w, h])
}

pub(super) fn build_catalog_keep_alive() -> (Catalog, tempfile::NamedTempFile) {
    let stars = vec![crate::sky_scene::SynthStar {
        source_id: 1,
        ra: 0.0,
        dec: 0.0,
        mag: 9.0,
    }];
    let file = crate::sky_scene::build(4, &stars);
    let cat = Catalog::open(file.path()).unwrap();
    (cat, file)
}

pub(super) fn pairs_on_tangent_plane(
    center_ra: f64,
    center_dec: f64,
    scale_arcsec: f64,
) -> Vec<StarPair> {
    // Build 8 pairs on a synthetic TAN tangent plane so fit_wcs has data
    // it can actually fit. crpix is at (512, 512); scale = arcsec/px;
    // each grid offset in pixels maps to (xi, eta) in degrees, then we
    // tan-deproject to (ra, dec) at the given center.
    use celestial_catalog::query::quad::tan_deproject_star;
    use celestial_core::constants::DEG_TO_RAD;
    let crpix = 512.0;
    let s_deg = scale_arcsec / 3600.0;
    let mut pairs = Vec::new();
    for &(dx, dy) in &[
        (-200.0, -200.0),
        (200.0, -200.0),
        (-200.0, 200.0),
        (200.0, 200.0),
        (-100.0, 50.0),
        (100.0, -50.0),
        (0.0, 0.0),
        (150.0, 150.0),
    ] {
        let px_x = crpix + dx;
        let px_y = crpix + dy;
        // Match the (-s, 0, 0, +s) CD matrix of WcsSolution::north_up.
        let xi_deg = -s_deg * dx;
        let eta_deg = s_deg * dy;
        let (ra, dec) = tan_deproject_star(
            xi_deg * DEG_TO_RAD,
            eta_deg * DEG_TO_RAD,
            center_ra,
            center_dec,
        );
        pairs.push(StarPair {
            px_x,
            px_y,
            ra_deg: ra,
            dec_deg: dec,
            votes: 4,
            snr: 20.0,
        });
    }
    pairs
}
