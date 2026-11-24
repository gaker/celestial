//! Debug image overlays for visualizing the solve pipeline.
//!
//! [`Annotation`] wraps an RGB buffer produced by auto-stretching the input image.
//! Draw primitives render matched/unmatched detections and solve residuals directly
//! onto that buffer, which can then be saved as PNG for inspection.
//!
//! Used by the `image-solver` binary's `--debug` flag.

mod raster;

use celestial_images::core::errors::ImageError;
use celestial_images::formats::image::Image;
use celestial_images::formats::pixel_data::PixelData;
use celestial_images::stretch::{StfParams, Stretch};

use crate::detect::DetectedStar;
use crate::error::SolveError;
use crate::fit_wcs::solution::WcsSolution;
use crate::match_field::StarPair;

/// Bright green, used for detections that landed in matched pairs.
const COLOR_MATCHED: [u8; 3] = [0, 255, 0];
/// Dim amber, used for detections that didn't match anything in the catalog.
const COLOR_UNMATCHED: [u8; 3] = [160, 120, 40];
/// Cyan, used for catalog-star positions projected through the fitted WCS.
const COLOR_CATALOG: [u8; 3] = [0, 200, 255];
/// Magenta, used for residual vectors between catalog positions and detected centroids.
const COLOR_RESIDUAL: [u8; 3] = [255, 0, 255];

/// RGB buffer for drawing solve-pipeline overlays.
///
/// Constructed from an [`Image`] via [`Annotation::from_image`], which stretches the
/// image for visibility and converts to RGB. Drawing methods mutate the buffer in
/// place; [`Annotation::to_image`] returns an [`Image`] suitable for saving.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::annotate::Annotation;
///
/// # let img = celestial_images::formats::image::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
///
/// let mut ann = Annotation::from_image(&img).unwrap();
/// ann.draw_detections(&result.stars, &result.pairs, 5);
/// ann.to_image().save("detections.png")?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct Annotation {
    /// Interleaved RGB bytes, row-major, `3 * width * height` elements.
    rgb: Vec<u8>,
    /// Image width in pixels.
    width: usize,
    /// Image height in pixels.
    height: usize,
}

impl Annotation {
    /// Builds an annotation buffer from an image, auto-stretched for visibility.
    ///
    /// Applies an STF (Screen Transfer Function) stretch to each channel to bring faint
    /// structure into the 8-bit range. A monochrome result is expanded to RGB. An RGB
    /// image keeps one RGB triple per pixel, but its channels are stretched
    /// independently, so its colour balance changes. Returns [`SolveError::Image`] if
    /// the stretch or conversion fails, or if the image is not one plane of 1 or 3
    /// channels.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_solver::annotate::Annotation;
    ///
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// let ann = Annotation::from_image(&img).unwrap();
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn from_image(image: &Image) -> Result<Self, SolveError> {
        let mut normed = image.clone();
        normed.normalize_to_f32();
        let stretched = normed.stf(&StfParams::default())?;
        let display = stretched.to_u8()?;
        Ok(Self {
            rgb: rgb_bytes(&display)?,
            width: display.width(),
            height: display.height(),
        })
    }

    /// Converts the buffer to an [`Image`] that can be saved as PNG / TIFF / etc.
    ///
    /// Clones the RGB buffer; the annotation is not consumed.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// let ann = celestial_solver::annotate::Annotation::from_image(&img).unwrap();
    /// let out = ann.to_image();
    /// out.save("annotated.png")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn to_image(&self) -> Image {
        Image::new(
            PixelData::U8(self.rgb.clone()),
            [self.width, self.height, 3],
        )
    }

    /// Draws circles around every detected star, colored by whether it matched.
    ///
    /// Detections within 1 pixel of a [`StarPair`] are drawn in bright green and the rest
    /// in dim amber. `radius` is the circle radius in
    /// pixels. Scale it to your image size (a radius of 5 is invisible on a 4K frame
    /// at preview zoom).
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// let mut ann = celestial_solver::annotate::Annotation::from_image(&img).unwrap();
    /// ann.draw_detections(&result.stars, &result.pairs, 8);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn draw_detections(&mut self, stars: &[DetectedStar], pairs: &[StarPair], radius: i32) {
        let margin = f64::from(radius) + 1.0;
        for s in stars {
            let Some((cx, cy)) = self.nearest_pixel_within((s.x, s.y), margin) else {
                continue;
            };
            let color = if is_matched(s, pairs) {
                COLOR_MATCHED
            } else {
                COLOR_UNMATCHED
            };
            self.draw_circle(cx, cy, i64::from(radius), color);
        }
    }

    /// Draws catalog-star crosshairs and scaled residual vectors.
    ///
    /// For each residual in `wcs.residuals`, draws a cyan crosshair at the
    /// projected catalog position and a magenta line from the detected
    /// centroid toward a scaled extrapolation of the error vector. `residual_scale`
    /// multiplies the residual length; sub-pixel residuals are invisible at scale 1,
    /// so pass something like 20 or 50 to make the pattern readable.
    ///
    /// A consistent directional pattern across the frame indicates sensor tilt; a
    /// radial pattern indicates optical distortion the TAN+SIP model didn't absorb;
    /// random scatter indicates noise-limited centroiding.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// let mut ann = celestial_solver::annotate::Annotation::from_image(&img).unwrap();
    /// ann.draw_solve_overlay(&result.wcs, 20.0);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn draw_solve_overlay(&mut self, wcs: &WcsSolution, residual_scale: f64) {
        const ARM: i64 = 6;
        for r in &wcs.residuals {
            let catalog = (r.px_x - r.err_x, r.px_y - r.err_y);
            if let Some((x, y)) = self.nearest_pixel_within(catalog, (ARM + 1) as f64) {
                self.draw_crosshair(x, y, ARM, COLOR_CATALOG);
            }
            let end = (
                r.px_x + r.err_x * (residual_scale - 1.0),
                r.px_y + r.err_y * (residual_scale - 1.0),
            );
            self.draw_line((r.px_x, r.px_y), end, COLOR_RESIDUAL);
        }
    }
}

// The stretched 8-bit samples as interleaved RGB: a monochrome plane is repeated in
// each channel and an RGB plane is used as it is.
fn rgb_bytes(display: &Image) -> Result<Vec<u8>, ImageError> {
    let samples = display.pixels.as_u8().ok_or_else(|| {
        ImageError::PixelConversion("the 8-bit stretch returned other pixels".to_string())
    })?;
    let (w, h, channels) = (display.width(), display.height(), display.channels());
    if samples.len() != w * h * channels {
        return Err(ImageError::InvalidDimensions(format!(
            "an overlay needs one image plane, got {} samples for {w}x{h}x{channels}",
            samples.len()
        )));
    }
    match channels {
        1 => Ok(samples.iter().flat_map(|&v| [v, v, v]).collect()),
        3 => Ok(samples.to_vec()),
        n => Err(ImageError::InvalidDimensions(format!(
            "an overlay needs 1 or 3 channels, got {n}"
        ))),
    }
}

fn is_matched(star: &DetectedStar, pairs: &[StarPair]) -> bool {
    const MATCH_TOL_PX: f64 = 1.0;
    pairs.iter().any(|p| {
        let dx = p.px_x - star.x;
        let dy = p.px_y - star.y;
        dx * dx + dy * dy <= MATCH_TOL_PX * MATCH_TOL_PX
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fit_wcs::solution::StarResidual;
    use celestial_images::formats::pixel_data::PixelData;

    fn blank_annotation(w: usize, h: usize) -> Annotation {
        Annotation {
            rgb: vec![0u8; w * h * 3],
            width: w,
            height: h,
        }
    }

    fn detection(x: f64, y: f64) -> DetectedStar {
        DetectedStar {
            x,
            y,
            flux: 1000.0,
            snr: 20.0,
            saturated: false,
            saturated_count: 0,
            background: 0.0,
        }
    }

    fn pair_at(x: f64, y: f64) -> StarPair {
        StarPair {
            px_x: x,
            px_y: y,
            ra_deg: 0.0,
            dec_deg: 0.0,
            votes: 4,
            snr: 20.0,
        }
    }

    fn wcs_with_residuals(residuals: Vec<StarResidual>) -> WcsSolution {
        WcsSolution {
            n_stars: residuals.len(),
            rms_px: 0.5,
            weighted_rms_px: 0.4,
            residuals,
            ..WcsSolution::diagonal((180.0, 30.0), (0.001, 0.001), (100, 100))
        }
    }

    fn any_pixel_matches(ann: &Annotation, color: [u8; 3]) -> bool {
        ann.rgb.as_chunks::<3>().0.contains(&color)
    }

    // An image with no finite pixel has no stretch, and the error says so where it was
    // `None`.
    #[test]
    fn from_image_of_a_blank_image_is_an_image_error() {
        let img = Image::new(PixelData::F32(vec![f32::NAN; 100]), [10, 10]);
        let err = Annotation::from_image(&img).err().map(|e| e.describe());
        assert_eq!(
            err.as_deref(),
            Some("Image: Invalid input: channel has no finite samples")
        );
    }

    #[test]
    fn from_image_produces_rgb_buffer() {
        let img = Image::new(PixelData::U8(vec![0; 100]), [10, 10]);
        let ann = Annotation::from_image(&img).unwrap();
        assert_eq!(ann.width, 10);
        assert_eq!(ann.height, 10);
        assert_eq!(ann.rgb.len(), 10 * 10 * 3);
    }

    // An RGB image is stretched per channel and keeps one RGB triple per pixel. The
    // 2x1 image has pixels (0, 10, 20) and (40, 30, 0), and the expected bytes are the
    // 8-bit STF of each channel, interleaved.
    #[test]
    fn from_image_of_an_rgb_image_keeps_one_triple_per_pixel() {
        let img = Image::new(PixelData::U8(vec![0, 10, 20, 40, 30, 0]), [2, 1, 3]);
        let ann = Annotation::from_image(&img).unwrap();
        assert_eq!((ann.width, ann.height), (2, 1));
        assert_eq!(ann.rgb, vec![0, 128, 128, 255, 230, 0]);
    }

    #[test]
    fn from_image_of_a_two_channel_image_is_an_image_error() {
        let img = Image::new(PixelData::U8(vec![0, 10, 20, 40]), [2, 1, 2]);
        let err = Annotation::from_image(&img).err().map(|e| e.describe());
        assert_eq!(
            err.as_deref(),
            Some("Image: Invalid dimensions: an overlay needs 1 or 3 channels, got 2")
        );
    }

    // A cube of two 2x1 planes is not one image to draw on.
    #[test]
    fn from_image_of_a_cube_is_an_image_error() {
        let img = Image::new(PixelData::U8(vec![0, 10, 20, 40]), [2, 1, 1, 2]);
        let err = Annotation::from_image(&img).err().map(|e| e.describe());
        assert_eq!(
            err.as_deref(),
            Some("Image: Invalid dimensions: an overlay needs one image plane, got 4 samples for 2x1x1")
        );
    }

    #[test]
    fn to_image_round_trips_dimensions_and_channel_count() {
        let ann = blank_annotation(8, 8);
        let img = ann.to_image();
        assert_eq!(img.width(), 8);
        assert_eq!(img.height(), 8);
        assert_eq!(img.channels(), 3);
    }

    #[test]
    fn matched_detections_use_matched_color() {
        let mut ann = blank_annotation(100, 100);
        let stars = vec![detection(50.0, 50.0)];
        let pairs = vec![pair_at(50.0, 50.0)];
        ann.draw_detections(&stars, &pairs, 5);
        assert!(any_pixel_matches(&ann, COLOR_MATCHED));
        assert!(!any_pixel_matches(&ann, COLOR_UNMATCHED));
    }

    #[test]
    fn unmatched_detections_use_unmatched_color() {
        let mut ann = blank_annotation(100, 100);
        let stars = vec![detection(20.0, 20.0), detection(80.0, 80.0)];
        let pairs = vec![pair_at(80.0, 80.0)];
        ann.draw_detections(&stars, &pairs, 5);
        assert!(any_pixel_matches(&ann, COLOR_MATCHED));
        assert!(any_pixel_matches(&ann, COLOR_UNMATCHED));
    }

    #[test]
    fn match_tolerance_allows_subpixel_jitter() {
        let mut ann = blank_annotation(100, 100);
        let stars = vec![detection(50.3, 50.7)];
        let pairs = vec![pair_at(50.0, 50.0)];
        ann.draw_detections(&stars, &pairs, 5);
        assert!(any_pixel_matches(&ann, COLOR_MATCHED));
    }

    #[test]
    fn match_tolerance_rejects_far_pairs() {
        let mut ann = blank_annotation(100, 100);
        let stars = vec![detection(50.0, 50.0)];
        let pairs = vec![pair_at(55.0, 55.0)];
        ann.draw_detections(&stars, &pairs, 5);
        assert!(any_pixel_matches(&ann, COLOR_UNMATCHED));
        assert!(!any_pixel_matches(&ann, COLOR_MATCHED));
    }

    #[test]
    fn empty_pairs_renders_all_detections_as_unmatched() {
        let mut ann = blank_annotation(100, 100);
        let stars = vec![detection(30.0, 30.0), detection(70.0, 70.0)];
        ann.draw_detections(&stars, &[], 5);
        assert!(any_pixel_matches(&ann, COLOR_UNMATCHED));
        assert!(!any_pixel_matches(&ann, COLOR_MATCHED));
    }

    #[test]
    fn solve_overlay_draws_catalog_crosshair_and_residual() {
        let mut ann = blank_annotation(200, 200);
        let residual = StarResidual {
            px_x: 100.0,
            px_y: 100.0,
            err_x: 0.5,
            err_y: 0.0,
            err_px: 0.5,
            ra_deg: 0.0,
            dec_deg: 0.0,
        };
        ann.draw_solve_overlay(&wcs_with_residuals(vec![residual]), 20.0);
        assert!(any_pixel_matches(&ann, COLOR_CATALOG));
        assert!(any_pixel_matches(&ann, COLOR_RESIDUAL));
    }

    #[test]
    fn residual_scale_controls_vector_length() {
        let mut short = blank_annotation(200, 200);
        let mut long = blank_annotation(200, 200);
        let residual = StarResidual {
            px_x: 100.0,
            px_y: 100.0,
            err_x: 0.5,
            err_y: 0.0,
            err_px: 0.5,
            ra_deg: 0.0,
            dec_deg: 0.0,
        };
        short.draw_solve_overlay(&wcs_with_residuals(vec![residual.clone()]), 2.0);
        long.draw_solve_overlay(&wcs_with_residuals(vec![residual]), 40.0);

        let count_magenta = |ann: &Annotation| {
            ann.rgb
                .as_chunks::<3>()
                .0
                .iter()
                .filter(|&&c| c == COLOR_RESIDUAL)
                .count()
        };
        assert!(count_magenta(&long) > count_magenta(&short));
    }

    fn residual(err_x: f64, err_y: f64) -> StarResidual {
        StarResidual {
            px_x: 100.0,
            px_y: 100.0,
            err_x,
            err_y,
            err_px: libm::hypot(err_x, err_y),
            ra_deg: 0.0,
            dec_deg: 0.0,
        }
    }

    fn pixels_of(ann: &Annotation, color: [u8; 3]) -> Vec<(usize, usize)> {
        let pixels = ann.rgb.as_chunks::<3>().0;
        (0..pixels.len())
            .filter(|&i| pixels[i] == color)
            .map(|i| (i % ann.width, i / ann.width))
            .collect()
    }

    // A residual of 1e9 px scaled 20 times ends about 2e10 px away. The vector stops at
    // the image edge, and the catalog crosshair 1e9 px to the left lands nowhere.
    #[test]
    fn huge_residual_is_drawn_to_the_image_edge() {
        let mut ann = blank_annotation(200, 200);
        ann.draw_solve_overlay(&wcs_with_residuals(vec![residual(1e9, 0.0)]), 20.0);
        let row: Vec<(usize, usize)> = (100..200).map(|x| (x, 100)).collect();
        assert_eq!(pixels_of(&ann, COLOR_RESIDUAL), row);
        assert_eq!(pixels_of(&ann, COLOR_CATALOG), Vec::new());
    }

    // The catalog position 1e10 px to the right is past the i32 range.
    #[test]
    fn crosshair_far_outside_the_image_draws_nothing() {
        let mut ann = blank_annotation(200, 200);
        ann.draw_solve_overlay(&wcs_with_residuals(vec![residual(-1e10, 0.0)]), 20.0);
        let row: Vec<(usize, usize)> = (0..=100).map(|x| (x, 100)).collect();
        assert_eq!(pixels_of(&ann, COLOR_RESIDUAL), row);
        assert_eq!(pixels_of(&ann, COLOR_CATALOG), Vec::new());
    }

    #[test]
    fn non_finite_residual_draws_nothing() {
        for err in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut ann = blank_annotation(200, 200);
            ann.draw_solve_overlay(&wcs_with_residuals(vec![residual(err, 1.0)]), 20.0);
            assert!(ann.rgb.iter().all(|&v| v == 0), "{err}");
        }
    }

    #[test]
    fn detection_far_outside_the_image_or_non_finite_draws_nothing() {
        for x in [1e12, -1e12, f64::NAN, f64::INFINITY] {
            let mut ann = blank_annotation(100, 100);
            ann.draw_detections(&[detection(x, 50.0)], &[], 5);
            assert!(ann.rgb.iter().all(|&v| v == 0), "{x}");
        }
    }

    // The ring of radius i32::MAX lies outside a 100x100 image centred on it.
    #[test]
    fn circle_enclosing_the_image_draws_nothing() {
        let mut ann = blank_annotation(100, 100);
        ann.draw_detections(&[detection(50.0, 50.0)], &[], i32::MAX);
        assert!(ann.rgb.iter().all(|&v| v == 0));
    }
}
