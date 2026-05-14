//! Debug image overlays for visualizing the solve pipeline.
//!
//! [`Annotation`] wraps an RGB buffer produced by auto-stretching the input image.
//! Draw primitives render matched/unmatched detections and solve residuals directly
//! onto that buffer, which can then be saved as PNG for inspection.
//!
//! Used by the `image-solver` binary's `--debug` flag.

use celestial_images::formats::Image;
use celestial_images::PixelData;
use celestial_images::stretch::{Stretch, StfParams};

use crate::detect::DetectedStar;
use crate::fit_wcs::WcsSolution;
use crate::match_field::StarPair;

/// Bright green, used for detections that landed in matched pairs.
pub const COLOR_MATCHED: [u8; 3] = [0, 255, 0];
/// Dim amber, used for detections that didn't match anything in the catalog.
pub const COLOR_UNMATCHED: [u8; 3] = [160, 120, 40];
/// Cyan, used for catalog-star positions projected through the fitted WCS.
pub const COLOR_CATALOG: [u8; 3] = [0, 200, 255];
/// Magenta, used for residual vectors between catalog positions and detected centroids.
pub const COLOR_RESIDUAL: [u8; 3] = [255, 0, 255];

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
/// # let img = celestial_images::formats::Image::open("f.fits")?;
/// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
/// let result = celestial_solver::solve(&img, &catalog).run()?;
///
/// let mut ann = Annotation::from_image(&img).unwrap();
/// ann.draw_detections(&result.stars, &result.pairs, 5);
/// ann.to_image().save("detections.png")?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct Annotation {
    /// Interleaved RGB bytes, row-major, `3 * width * height` elements.
    pub rgb: Vec<u8>,
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
}

impl Annotation {
    /// Builds an annotation buffer from an image, auto-stretched for visibility.
    ///
    /// Applies an STF (Screen Transfer Function) stretch to bring faint structure into
    /// the 8-bit range, then expands the monochrome result to RGB. Returns `None` if
    /// the stretch or conversion fails.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_solver::annotate::Annotation;
    ///
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// let ann = Annotation::from_image(&img).unwrap();
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn from_image(image: &Image) -> Option<Self> {
        let mut normed = image.clone();
        normed.normalize_to_f32();
        let stretched = normed.stf(&StfParams::default()).ok()?;
        let mono = stretched.to_u8().ok()?;
        let w = mono.width();
        let h = mono.height();
        let gray = mono.pixels.as_u8()?;
        let mut rgb = Vec::with_capacity(gray.len() * 3);
        for &v in gray {
            rgb.push(v);
            rgb.push(v);
            rgb.push(v);
        }
        Some(Self { rgb, width: w, height: h })
    }

    /// Converts the buffer to an [`Image`] that can be saved as PNG / TIFF / etc.
    ///
    /// Clones the RGB buffer; the annotation is not consumed.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
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
    /// Detections within 1 pixel of a [`StarPair`] use [`COLOR_MATCHED`] (bright green);
    /// the rest use [`COLOR_UNMATCHED`] (dim amber). `radius` is the circle radius in
    /// pixels — scale it to your image size (a radius of 5 is invisible on a 4K frame
    /// at preview zoom).
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// let mut ann = celestial_solver::annotate::Annotation::from_image(&img).unwrap();
    /// ann.draw_detections(&result.stars, &result.pairs, 8);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn draw_detections(
        &mut self,
        stars: &[DetectedStar],
        pairs: &[StarPair],
        radius: i32,
    ) {
        for s in stars {
            let color = if is_matched(s, pairs) {
                COLOR_MATCHED
            } else {
                COLOR_UNMATCHED
            };
            let cx = s.x.round() as i32;
            let cy = s.y.round() as i32;
            draw_circle(&mut self.rgb, self.width, self.height, cx, cy, radius, color);
        }
    }

    /// Draws catalog-star crosshairs and scaled residual vectors.
    ///
    /// For each residual in `wcs.residuals`, draws a [`COLOR_CATALOG`] crosshair at the
    /// projected catalog position and a [`COLOR_RESIDUAL`] line from the detected
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
    /// # let img = celestial_images::formats::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve(&img, &catalog).run()?;
    /// let mut ann = celestial_solver::annotate::Annotation::from_image(&img).unwrap();
    /// ann.draw_solve_overlay(&result.wcs, 20.0);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn draw_solve_overlay(&mut self, wcs: &WcsSolution, residual_scale: f64) {
        for r in &wcs.residuals {
            let cat_x = (r.px_x - r.err_x).round() as i32;
            let cat_y = (r.px_y - r.err_y).round() as i32;
            draw_crosshair(&mut self.rgb, self.width, self.height, cat_x, cat_y, 6, COLOR_CATALOG);

            let end_x = (r.px_x + r.err_x * (residual_scale - 1.0)).round() as i32;
            let end_y = (r.px_y + r.err_y * (residual_scale - 1.0)).round() as i32;
            let cx = r.px_x.round() as i32;
            let cy = r.px_y.round() as i32;
            draw_line(&mut self.rgb, self.width, self.height, (cx, cy), (end_x, end_y), COLOR_RESIDUAL);
        }
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

fn set_pixel(rgb: &mut [u8], w: usize, h: usize, x: i32, y: i32, color: [u8; 3]) {
    if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 {
        return;
    }
    let idx = (y as usize * w + x as usize) * 3;
    rgb[idx] = color[0];
    rgb[idx + 1] = color[1];
    rgb[idx + 2] = color[2];
}

fn draw_circle(rgb: &mut [u8], w: usize, h: usize, cx: i32, cy: i32, r: i32, color: [u8; 3]) {
    let mut x = r;
    let mut y = 0i32;
    let mut d = 1 - r;
    while x >= y {
        set_pixel(rgb, w, h, cx + x, cy + y, color);
        set_pixel(rgb, w, h, cx - x, cy + y, color);
        set_pixel(rgb, w, h, cx + x, cy - y, color);
        set_pixel(rgb, w, h, cx - x, cy - y, color);
        set_pixel(rgb, w, h, cx + y, cy + x, color);
        set_pixel(rgb, w, h, cx - y, cy + x, color);
        set_pixel(rgb, w, h, cx + y, cy - x, color);
        set_pixel(rgb, w, h, cx - y, cy - x, color);
        y += 1;
        if d <= 0 {
            d += 2 * y + 1;
        } else {
            x -= 1;
            d += 2 * (y - x) + 1;
        }
    }
}

fn draw_crosshair(rgb: &mut [u8], w: usize, h: usize, cx: i32, cy: i32, arm: i32, color: [u8; 3]) {
    for d in -arm..=arm {
        if d.abs() <= 1 {
            continue;
        }
        set_pixel(rgb, w, h, cx + d, cy, color);
        set_pixel(rgb, w, h, cx, cy + d, color);
    }
}

fn draw_line(
    rgb: &mut [u8], w: usize, h: usize,
    p0: (i32, i32), p1: (i32, i32),
    color: [u8; 3],
) {
    let (mut x, mut y) = p0;
    let (x1, y1) = p1;
    let dx = (x1 - x).abs();
    let dy = -(y1 - y).abs();
    let sx = if x < x1 { 1 } else { -1 };
    let sy = if y < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    loop {
        set_pixel(rgb, w, h, x, y, color);
        if x == x1 && y == y1 { break; }
        let e2 = 2 * err;
        if e2 >= dy { err += dy; x += sx; }
        if e2 <= dx { err += dx; y += sy; }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fit_wcs::solution::StarResidual;
    use celestial_images::formats::PixelData;

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
            crpix1: 50.0,
            crpix2: 50.0,
            crval1: 180.0,
            crval2: 30.0,
            cd1_1: 0.001,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: 0.001,
            width: 100,
            height: 100,
            focal_mm: None,
            pixel_um: None,
            n_stars: residuals.len(),
            rms_px: 0.5,
            weighted_rms_px: 0.4,
            residuals,
        }
    }

    fn pixel(ann: &Annotation, x: i32, y: i32) -> [u8; 3] {
        let idx = (y as usize * ann.width + x as usize) * 3;
        [ann.rgb[idx], ann.rgb[idx + 1], ann.rgb[idx + 2]]
    }

    fn any_pixel_matches(ann: &Annotation, color: [u8; 3]) -> bool {
        ann.rgb.chunks_exact(3).any(|c| c == color)
    }

    #[test]
    fn from_image_produces_rgb_buffer() {
        let img = Image::new(PixelData::U8(vec![0; 100]), [10, 10]);
        let ann = Annotation::from_image(&img).unwrap();
        assert_eq!(ann.width, 10);
        assert_eq!(ann.height, 10);
        assert_eq!(ann.rgb.len(), 10 * 10 * 3);
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
        };
        short.draw_solve_overlay(&wcs_with_residuals(vec![residual.clone()]), 2.0);
        long.draw_solve_overlay(&wcs_with_residuals(vec![residual]), 40.0);

        let count_magenta = |ann: &Annotation| {
            ann.rgb.chunks_exact(3).filter(|c| *c == COLOR_RESIDUAL).count()
        };
        assert!(count_magenta(&long) > count_magenta(&short));
    }

    #[test]
    fn set_pixel_ignores_out_of_bounds() {
        let mut ann = blank_annotation(10, 10);
        set_pixel(&mut ann.rgb, ann.width, ann.height, -1, 5, [255, 0, 0]);
        set_pixel(&mut ann.rgb, ann.width, ann.height, 5, 10, [255, 0, 0]);
        set_pixel(&mut ann.rgb, ann.width, ann.height, 100, 100, [255, 0, 0]);
        assert!(ann.rgb.iter().all(|&v| v == 0));
    }

    #[test]
    fn draw_circle_clips_near_edge() {
        let mut ann = blank_annotation(10, 10);
        draw_circle(&mut ann.rgb, ann.width, ann.height, 0, 0, 5, [255, 255, 255]);
        assert_eq!(pixel(&ann, 5, 0), [255, 255, 255]);
        assert_eq!(pixel(&ann, 0, 5), [255, 255, 255]);
    }
}
