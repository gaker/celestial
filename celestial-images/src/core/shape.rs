// Width, height and channels of an image with these axis lengths. A missing
// height or channel axis has length 1. An image with no axes holds no
// samples (NAXIS = 0 in FITS), so its width is 0.
pub(crate) fn plane_shape(dimensions: &[usize]) -> (usize, usize, usize) {
    let axis = |index: usize, missing: usize| dimensions.get(index).copied().unwrap_or(missing);
    (axis(0, 0), axis(1, 1), axis(2, 1))
}

// PNG, TIFF and XISF store a width and a height, so they cannot hold an
// image with no axes. `error` wraps the message in the caller's error type.
pub(crate) fn plane_shape_with_axes<E>(
    dimensions: &[usize],
    error: impl FnOnce(String) -> E,
) -> Result<(usize, usize, usize), E> {
    if dimensions.is_empty() {
        return Err(error("the image has no axes".to_string()));
    }
    Ok(plane_shape(dimensions))
}

// The axis lengths of a width x height image with `channels` samples per
// pixel. A single channel needs no third axis.
pub(crate) fn plane_dimensions(width: usize, height: usize, channels: usize) -> Vec<usize> {
    if channels > 1 {
        vec![width, height, channels]
    } else {
        vec![width, height]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_axes_take_their_defaults() {
        assert_eq!(plane_shape(&[]), (0, 1, 1));
        assert_eq!(plane_shape(&[5]), (5, 1, 1));
        assert_eq!(plane_shape(&[5, 4]), (5, 4, 1));
        assert_eq!(plane_shape(&[5, 4, 3]), (5, 4, 3));
        assert_eq!(plane_shape(&[5, 4, 3, 2]), (5, 4, 3));
    }

    #[test]
    fn one_channel_has_two_axes() {
        assert_eq!(plane_dimensions(1920, 1080, 1), vec![1920, 1080]);
        assert_eq!(plane_dimensions(100, 50, 0), vec![100, 50]);
        assert_eq!(plane_dimensions(1920, 1080, 3), vec![1920, 1080, 3]);
        assert_eq!(plane_dimensions(64, 64, 4), vec![64, 64, 4]);
    }
}
