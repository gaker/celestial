use anyhow::Result;
use std::fs::create_dir_all;
use std::path::Path;

use celestial_images::formats::image::Image;
use celestial_solver::annotate::Annotation;
use celestial_solver::solve::SolveResult;

pub(crate) fn write(dir: &Path, img: &Image, result: &SolveResult) -> Result<()> {
    create_dir_all(dir)?;

    let mut ann = Annotation::from_image(img)?;
    ann.draw_detections(&result.stars, &result.pairs, 5);
    ann.to_image().save(dir.join("detections.png"))?;

    let mut ann = Annotation::from_image(img)?;
    ann.draw_solve_overlay(&result.wcs, 20.0);
    ann.to_image().save(dir.join("solve.png"))?;

    eprintln!("debug images written to {}", dir.display());
    Ok(())
}
