#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImageKind {
    Mono,
    Rgb,
    Cube,
}

impl ImageKind {
    pub fn from_dimensions(dims: &[usize]) -> Self {
        match dims.len() {
            1 | 2 => Self::Mono,
            3 if dims[2] == 3 => Self::Rgb,
            _ => Self::Cube,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_and_two_dim_are_mono() {
        assert_eq!(ImageKind::from_dimensions(&[100]), ImageKind::Mono);
        assert_eq!(ImageKind::from_dimensions(&[100, 200]), ImageKind::Mono);
    }

    #[test]
    fn three_channel_three_dim_is_rgb() {
        assert_eq!(ImageKind::from_dimensions(&[100, 200, 3]), ImageKind::Rgb);
        assert_eq!(ImageKind::from_dimensions(&[1, 1, 3]), ImageKind::Rgb);
    }

    #[test]
    fn other_three_dim_is_cube() {
        assert_eq!(ImageKind::from_dimensions(&[100, 200, 1]), ImageKind::Cube);
        assert_eq!(ImageKind::from_dimensions(&[100, 200, 4]), ImageKind::Cube);
        assert_eq!(ImageKind::from_dimensions(&[100, 200, 7]), ImageKind::Cube);
    }

    #[test]
    fn four_plus_dim_is_cube() {
        assert_eq!(
            ImageKind::from_dimensions(&[10, 10, 3, 5]),
            ImageKind::Cube
        );
        assert_eq!(
            ImageKind::from_dimensions(&[10, 10, 3, 5, 2]),
            ImageKind::Cube
        );
    }

    #[test]
    fn empty_dimensions_is_cube() {
        assert_eq!(ImageKind::from_dimensions(&[]), ImageKind::Cube);
    }
}
