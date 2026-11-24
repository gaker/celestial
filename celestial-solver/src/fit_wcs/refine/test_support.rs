// Detections and catalogs shared by the refinement tests.

use celestial_catalog::query::catalog::Catalog;

use crate::detect::DetectedStar;
use crate::sky_scene::{build, SynthStar};

pub(super) fn det(x: f64, y: f64, snr: f64) -> DetectedStar {
    DetectedStar {
        x,
        y,
        flux: 1000.0,
        snr,
        saturated: false,
        saturated_count: 0,
        background: 0.0,
    }
}

pub(super) fn synth_catalog(stars: Vec<SynthStar>) -> (Catalog, tempfile::NamedTempFile) {
    let file = build(4, &stars);
    let c = Catalog::open(file.path()).unwrap();
    (c, file)
}
