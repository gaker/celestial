use super::Pixel;

// 90% of the image maximum when two or more finite pixels share it and at least one of
// them is above the background there (a clipped sensor). Otherwise, for example on an
// all-zero or flat image, there is no saturation limit.
pub(super) fn estimate_saturation<T: Pixel>(image: &[T], background: &[f64]) -> f64 {
    let mut max_val = f64::NEG_INFINITY;
    let mut count_at_max: u32 = 0;
    let mut above_background = false;
    for (&px, &bg) in image.iter().zip(background) {
        let v = px.to_f64();
        if !v.is_finite() || v < max_val {
            continue;
        }
        if v > max_val {
            max_val = v;
            count_at_max = 0;
            above_background = false;
        }
        count_at_max = count_at_max.saturating_add(1);
        above_background |= v > bg;
    }
    if count_at_max >= 2 && above_background {
        max_val * 0.9
    } else {
        f64::INFINITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn estimate_saturation_unique_max_returns_infinity() {
        let image: Vec<f32> = (0..100).map(|i| i as f32).collect();
        assert_eq!(
            estimate_saturation(&image, &vec![0.0; image.len()]),
            f64::INFINITY
        );
    }

    #[test]
    fn estimate_saturation_tied_max_returns_fraction() {
        let image = vec![10.0_f32, 50.0, 30.0, 50.0, 20.0];
        let limit = estimate_saturation(&image, &vec![0.0; image.len()]);
        assert_eq!(limit, 45.0);
    }

    #[test]
    fn estimate_saturation_many_at_max() {
        let mut image = vec![0.0_f32; 1000];
        for v in image.iter_mut().take(50) {
            *v = 65535.0;
        }
        let limit = estimate_saturation(&image, &vec![0.0; image.len()]);
        assert_eq!(limit, 65535.0 * 0.9);
    }

    #[test]
    fn estimate_saturation_empty_image() {
        let image: Vec<f32> = vec![];
        assert_eq!(
            estimate_saturation(&image, &vec![0.0; image.len()]),
            f64::INFINITY
        );
    }

    #[test]
    fn estimate_saturation_single_pixel() {
        let image = vec![1000.0_f32];
        assert_eq!(
            estimate_saturation(&image, &vec![0.0; image.len()]),
            f64::INFINITY
        );
    }

    #[test]
    fn estimate_saturation_needs_the_maximum_above_the_background() {
        let image = vec![10.0_f32, 50.0, 30.0, 50.0, 20.0];
        assert_eq!(estimate_saturation(&image, &[50.0; 5]), f64::INFINITY);
        assert_eq!(
            estimate_saturation(&image, &[49.0, 50.0, 50.0, 50.0, 50.0]),
            f64::INFINITY
        );
        assert_eq!(
            estimate_saturation(&image, &[0.0, 0.0, 0.0, 49.0, 0.0]),
            50.0 * 0.9
        );
    }

    #[test]
    fn estimate_saturation_skips_non_finite_pixels() {
        let image = vec![f32::INFINITY, 50.0, f32::NAN, 50.0, f32::INFINITY];
        assert_eq!(estimate_saturation(&image, &[0.0; 5]), 50.0 * 0.9);
    }

    #[test]
    fn estimate_saturation_all_zero() {
        let image = vec![0.0_f32; 100];
        let limit = estimate_saturation(&image, &vec![0.0; image.len()]);
        assert_eq!(limit, f64::INFINITY);
    }
}
