// Pixel drawing on the RGB buffer of an `Annotation`. Every primitive clips to the
// image, and its number of steps is bounded by the image size or the circle radius.

use super::Annotation;

type Point = (f64, f64);

impl Annotation {
    // The pixel nearest to `p`, or None when `p` is non-finite or more than `margin`
    // pixels outside the image, where nothing drawn within `margin` of it lands on the
    // image. The bound also keeps the returned coordinates far inside the i64 range.
    pub(super) fn nearest_pixel_within(&self, p: Point, margin: f64) -> Option<(i64, i64)> {
        let near = |v: f64, len: usize| v >= -0.5 - margin && v <= len as f64 - 0.5 + margin;
        (near(p.0, self.width) && near(p.1, self.height)).then(|| nearest_pixel(p))
    }

    fn set_pixel(&mut self, x: i64, y: i64, color: [u8; 3]) {
        let (Ok(x), Ok(y)) = (usize::try_from(x), usize::try_from(y)) else {
            return;
        };
        if x >= self.width || y >= self.height {
            return;
        }
        let idx = (y * self.width + x) * 3;
        self.rgb[idx..idx + 3].copy_from_slice(&color);
    }

    pub(super) fn draw_circle(&mut self, cx: i64, cy: i64, r: i64, color: [u8; 3]) {
        if self.ring_encloses_image(cx, cy, r) {
            return;
        }
        let mut x = r;
        let mut y = 0i64;
        let mut d = 1 - r;
        while x >= y {
            for (a, b) in [(x, y), (y, x)] {
                for (sa, sb) in [(1, 1), (-1, 1), (1, -1), (-1, -1)] {
                    self.set_pixel(cx + sa * a, cy + sb * b, color);
                }
            }
            y += 1;
            if d <= 0 {
                d += 2 * y + 1;
            } else {
                x -= 1;
                d += 2 * (y - x) + 1;
            }
        }
    }

    // True when every pixel of the image is more than a pixel inside the ring, so the
    // ring draws nothing and stepping round it would only cost time.
    fn ring_encloses_image(&self, cx: i64, cy: i64, r: i64) -> bool {
        let far = |c: i64, len: usize| c.abs().max((len as i64 - 1 - c).abs()) as f64;
        libm::hypot(far(cx, self.width), far(cy, self.height)) + 1.0 < r as f64
    }

    pub(super) fn draw_crosshair(&mut self, cx: i64, cy: i64, arm: i64, color: [u8; 3]) {
        for d in -arm..=arm {
            if d.abs() <= 1 {
                continue;
            }
            self.set_pixel(cx + d, cy, color);
            self.set_pixel(cx, cy + d, color);
        }
    }

    // Bresenham's line between the pixels nearest to `p0` and `p1`, after clipping the
    // segment to the image so the number of steps is bounded by the image size.
    pub(super) fn draw_line(&mut self, p0: Point, p1: Point, color: [u8; 3]) {
        let Some((p0, p1)) = clip_to_image(p0, p1, self.width, self.height) else {
            return;
        };
        let ((mut x, mut y), (x1, y1)) = (nearest_pixel(p0), nearest_pixel(p1));
        let dx = (x1 - x).abs();
        let dy = -(y1 - y).abs();
        let sx = if x < x1 { 1 } else { -1 };
        let sy = if y < y1 { 1 } else { -1 };
        let mut err = dx + dy;
        loop {
            self.set_pixel(x, y, color);
            if x == x1 && y == y1 {
                break;
            }
            let e2 = 2 * err;
            if e2 >= dy {
                err += dy;
                x += sx;
            }
            if e2 <= dx {
                err += dx;
                y += sy;
            }
        }
    }
}

fn nearest_pixel(p: Point) -> (i64, i64) {
    (libm::round(p.0) as i64, libm::round(p.1) as i64)
}

// Liang-Barsky clipping of the segment to the pixel area of the image, the rectangle
// from -0.5 to w - 0.5 and h - 0.5. An endpoint inside is returned unchanged, so a
// segment inside the image draws the same pixels as without clipping. None when an
// endpoint or the difference between them is non-finite, or the segment misses the image.
fn clip_to_image(p0: Point, p1: Point, w: usize, h: usize) -> Option<(Point, Point)> {
    let (dx, dy) = (p1.0 - p0.0, p1.1 - p0.1);
    if ![p0.0, p0.1, dx, dy].iter().all(|v| v.is_finite()) {
        return None;
    }
    let (x_max, y_max) = (w as f64 - 0.5, h as f64 - 0.5);
    let edges = [
        (-dx, p0.0 + 0.5),
        (dx, x_max - p0.0),
        (-dy, p0.1 + 0.5),
        (dy, y_max - p0.1),
    ];
    let (t0, t1) = edges
        .into_iter()
        .try_fold((0.0, 1.0), |t, (p, q)| clip_parameter(t, p, q))?;
    // Rounding in t, or in p0 + t * d at large coordinates, can leave a clipped end
    // outside the area, so each end is clamped onto it. That bounds the steps drawn.
    let onto_area = |p: Point| (p.0.clamp(-0.5, x_max), p.1.clamp(-0.5, y_max));
    let at = |t: f64| (p0.0 + t * dx, p0.1 + t * dy);
    let start = if t0 > 0.0 { at(t0) } else { p0 };
    let end = if t1 < 1.0 { at(t1) } else { p1 };
    Some((onto_area(start), onto_area(end)))
}

// Narrows the parameter range `(t0, t1)` of the segment to the inside of one edge. `q`
// is the distance of the start inside the edge and `-p` the rate the segment moves
// further inside as the parameter grows. None when nothing of the range is inside.
fn clip_parameter((t0, t1): (f64, f64), p: f64, q: f64) -> Option<(f64, f64)> {
    if p == 0.0 {
        return (q >= 0.0).then_some((t0, t1));
    }
    let range = if p < 0.0 {
        (f64::max(t0, q / p), t1)
    } else {
        (t0, f64::min(t1, q / p))
    };
    (range.0 <= range.1).then_some(range)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blank(w: usize, h: usize) -> Annotation {
        Annotation {
            rgb: vec![0u8; w * h * 3],
            width: w,
            height: h,
        }
    }

    fn pixel(ann: &Annotation, x: usize, y: usize) -> [u8; 3] {
        let idx = (y * ann.width + x) * 3;
        [ann.rgb[idx], ann.rgb[idx + 1], ann.rgb[idx + 2]]
    }

    #[test]
    fn set_pixel_ignores_out_of_bounds() {
        let mut ann = blank(10, 10);
        ann.set_pixel(-1, 5, [255, 0, 0]);
        ann.set_pixel(5, 10, [255, 0, 0]);
        ann.set_pixel(100, 100, [255, 0, 0]);
        assert!(ann.rgb.iter().all(|&v| v == 0));
    }

    #[test]
    fn draw_circle_clips_near_edge() {
        let mut ann = blank(10, 10);
        ann.draw_circle(0, 0, 5, [255, 255, 255]);
        assert_eq!(pixel(&ann, 5, 0), [255, 255, 255]);
        assert_eq!(pixel(&ann, 0, 5), [255, 255, 255]);
    }

    #[test]
    fn clip_leaves_a_segment_inside_the_image_unchanged() {
        let segment = ((1.3, 2.7), (8.1, 0.2));
        assert_eq!(clip_to_image(segment.0, segment.1, 10, 10), Some(segment));
    }

    // From x = -8.5 to 23.5 the segment enters the pixel area at x = -0.5, a quarter of
    // the way along, and leaves it at x = 9.5.
    #[test]
    fn clip_cuts_a_crossing_segment_at_the_edges_of_the_pixel_area() {
        assert_eq!(
            clip_to_image((-8.5, 5.0), (23.5, 5.0), 10, 10),
            Some(((-0.5, 5.0), (9.5, 5.0)))
        );
    }

    #[test]
    fn clip_of_a_segment_that_misses_the_image_is_none() {
        assert_eq!(clip_to_image((-5.0, -5.0), (20.0, -1.0), 10, 10), None);
        assert_eq!(clip_to_image((-5.0, 3.0), (3.0, -5.0), 10, 10), None);
        assert_eq!(clip_to_image((12.0, 3.0), (12.0, 3.0), 10, 10), None);
    }

    // With p0.x = 1e300 the parameter where the segment leaves the image rounds to 1.0,
    // so the far endpoint, 1e9 px outside, is not moved onto the edge by the parameter.
    // Precision at 1e300 also loses where the segment enters, so only the bound holds.
    #[test]
    fn clip_keeps_both_ends_on_the_pixel_area_when_the_parameter_rounds() {
        let (start, end) = clip_to_image((1e300, -0.5), (-1e9, -0.5), 10, 10)
            .expect("the segment runs along the top edge of the pixel area");
        let on_area = |p: Point| (-0.5..=9.5).contains(&p.0) && (-0.5..=9.5).contains(&p.1);
        assert!(on_area(start), "start {start:?}");
        assert!(on_area(end), "end {end:?}");
    }

    #[test]
    fn clip_of_a_non_finite_endpoint_is_none() {
        assert_eq!(clip_to_image((1.0, 1.0), (f64::NAN, 2.0), 10, 10), None);
        assert_eq!(
            clip_to_image((f64::INFINITY, 1.0), (2.0, 2.0), 10, 10),
            None
        );
    }
}
