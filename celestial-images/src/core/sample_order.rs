// `channels` samples per pixel. Trailing samples that do not fill a pixel
// keep their positions.
pub(crate) fn to_planar<T: Copy>(data: &[T], channels: usize) -> Vec<T> {
    let mut planar = data.to_vec();
    let Some(pixel_count) = data.len().checked_div(channels) else {
        return planar;
    };
    for (i, pixel) in data.chunks_exact(channels).enumerate() {
        for (c, &v) in pixel.iter().enumerate() {
            planar[c * pixel_count + i] = v;
        }
    }
    planar
}

pub(crate) fn to_interleaved<T: Copy>(data: &[T], channels: usize) -> Vec<T> {
    let mut interleaved = data.to_vec();
    let Some(pixel_count) = data.len().checked_div(channels) else {
        return interleaved;
    };
    for (i, pixel) in interleaved.chunks_exact_mut(channels).enumerate() {
        for (c, v) in pixel.iter_mut().enumerate() {
            *v = data[c * pixel_count + i];
        }
    }
    interleaved
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn planar_conversion_handles_any_channel_count() {
        let interleaved = [1, 10, 2, 20, 3, 30, 4];
        let planar = to_planar(&interleaved, 2);
        assert_eq!(planar, [1, 2, 3, 10, 20, 30, 4]);
        assert_eq!(to_interleaved(&planar, 2), interleaved);
        assert_eq!(to_planar(&interleaved, 0), interleaved);
    }
}
