use crate::core::errors::{ImageError, Result};
use crate::ser::types::ColorId;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BayerPattern {
    Rggb,
    Grbg,
    Gbrg,
    Bggr,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Cfa {
    Red,
    Green,
    Blue,
}

// A CFA sample position, with a green site split by the colour of its left
// and right neighbours.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Site {
    Red,
    Blue,
    GreenInRedRow,
    GreenInBlueRow,
}

impl BayerPattern {
    pub fn from_color_id(color_id: ColorId) -> Option<Self> {
        match color_id {
            ColorId::BayerRggb => Some(Self::Rggb),
            ColorId::BayerGrbg => Some(Self::Grbg),
            ColorId::BayerGbrg => Some(Self::Gbrg),
            ColorId::BayerBggr => Some(Self::Bggr),
            _ => None,
        }
    }

    // The pattern names the 2x2 tile in row-major order: (0,0), (1,0), (0,1), (1,1).
    fn colour_at(self, x: usize, y: usize) -> Cfa {
        use Cfa::{Blue, Green, Red};
        let tile = match self {
            Self::Rggb => [Red, Green, Green, Blue],
            Self::Grbg => [Green, Red, Blue, Green],
            Self::Gbrg => [Green, Blue, Red, Green],
            Self::Bggr => [Blue, Green, Green, Red],
        };
        tile[(y % 2) * 2 + x % 2]
    }

    // The sites at even and odd x on row y.
    fn row_sites(self, y: usize) -> [Site; 2] {
        let site = |x: usize| match self.colour_at(x, y) {
            Cfa::Red => Site::Red,
            Cfa::Blue => Site::Blue,
            Cfa::Green if self.colour_at(x + 1, y) == Cfa::Red => Site::GreenInRedRow,
            Cfa::Green => Site::GreenInBlueRow,
        };
        [site(0), site(1)]
    }
}

pub fn debayer_bilinear_u8(
    raw: &[u8],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Result<Vec<u8>> {
    debayer_bilinear(raw, width, height, pattern)
}

pub fn debayer_bilinear_u16(
    raw: &[u16],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Result<Vec<u16>> {
    debayer_bilinear(raw, width, height, pattern)
}

// Mosaic samples are summed in u32: four u16 values cannot overflow it.
trait MosaicSample: Copy + Default + Send + Sync + Into<u32> {
    fn from_mean(mean: u32) -> Self;
}

// A mean of samples of one type never exceeds that type's maximum, so the
// casts below do not truncate.
impl MosaicSample for u8 {
    fn from_mean(mean: u32) -> Self {
        mean as u8
    }
}

impl MosaicSample for u16 {
    fn from_mean(mean: u32) -> Self {
        mean as u16
    }
}

struct Mosaic<'a, T> {
    raw: &'a [T],
    width: usize,
    height: usize,
    pattern: BayerPattern,
}

fn debayer_bilinear<T: MosaicSample>(
    raw: &[T],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Result<Vec<T>> {
    let row_len = checked_rgb_row_len(raw.len(), width, height)?;
    let mosaic = Mosaic {
        raw,
        width,
        height,
        pattern,
    };
    let mut rgb = vec![T::default(); raw.len() * 3];
    #[cfg(feature = "parallel")]
    let rows = rgb.par_chunks_mut(row_len).enumerate();
    #[cfg(not(feature = "parallel"))]
    let rows = rgb.chunks_mut(row_len).enumerate();
    rows.for_each(|(y, row)| mosaic.fill_row(y, row));
    Ok(rgb)
}

// Returns the length of one RGB output row after checking that the mosaic
// holds exactly width x height samples.
fn checked_rgb_row_len(len: usize, width: usize, height: usize) -> Result<usize> {
    let invalid = |reason: &str| {
        ImageError::InvalidDimensions(format!("Bayer mosaic {width}x{height}: {reason}"))
    };
    if width == 0 || height == 0 {
        return Err(invalid("width and height must be nonzero"));
    }
    let samples = width
        .checked_mul(height)
        .ok_or_else(|| invalid("sample count overflows usize"))?;
    if len != samples {
        return Err(invalid(&format!("expected {samples} samples, got {len}")));
    }
    samples
        .checked_mul(3)
        .map(|_| width * 3)
        .ok_or_else(|| invalid("RGB sample count overflows usize"))
}

// The 3x3 neighbourhood of one sample: rows above, at and below it, and
// columns left of, at and right of it. Neighbours outside the mosaic take the
// nearest edge sample.
struct Neighbours<'a, T> {
    rows: [&'a [T]; 3],
    cols: [usize; 3],
}

impl<T: MosaicSample> Neighbours<'_, T> {
    fn at(&self, row: usize, col: usize) -> u32 {
        self.rows[row][self.cols[col]].into()
    }

    fn horizontal(&self) -> u32 {
        (self.at(1, 0) + self.at(1, 2)) / 2
    }

    fn vertical(&self) -> u32 {
        (self.at(0, 1) + self.at(2, 1)) / 2
    }

    fn cross(&self) -> u32 {
        (self.at(1, 0) + self.at(1, 2) + self.at(0, 1) + self.at(2, 1)) / 4
    }

    fn diagonal(&self) -> u32 {
        (self.at(0, 0) + self.at(0, 2) + self.at(2, 0) + self.at(2, 2)) / 4
    }

    // Each site computes only the means its missing colours need.
    fn interpolate(&self, site: Site) -> [T; 3] {
        let centre = self.at(1, 1);
        let [r, g, b] = match site {
            Site::Red => [centre, self.cross(), self.diagonal()],
            Site::Blue => [self.diagonal(), self.cross(), centre],
            Site::GreenInRedRow => [self.horizontal(), centre, self.vertical()],
            Site::GreenInBlueRow => [self.vertical(), centre, self.horizontal()],
        };
        [T::from_mean(r), T::from_mean(g), T::from_mean(b)]
    }
}

impl<T: MosaicSample> Mosaic<'_, T> {
    fn row(&self, y: usize) -> &[T] {
        &self.raw[y * self.width..(y + 1) * self.width]
    }

    fn fill_row(&self, y: usize, out: &mut [T]) {
        let rows = [
            self.row(y.saturating_sub(1)),
            self.row(y),
            self.row((y + 1).min(self.height - 1)),
        ];
        let sites = self.pattern.row_sites(y);
        for (x, pixel) in out.as_chunks_mut::<3>().0.iter_mut().enumerate() {
            let cols = [x.saturating_sub(1), x, (x + 1).min(self.width - 1)];
            *pixel = Neighbours { rows, cols }.interpolate(sites[x % 2]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PATTERNS: [BayerPattern; 4] = [
        BayerPattern::Rggb,
        BayerPattern::Bggr,
        BayerPattern::Grbg,
        BayerPattern::Gbrg,
    ];

    fn mosaic_u8(n: usize) -> Vec<u8> {
        (0..n).map(|i| ((i * 37 + 11) % 251) as u8).collect()
    }

    #[test]
    fn bayer_pattern_from_color_id() {
        assert_eq!(
            BayerPattern::from_color_id(ColorId::BayerRggb),
            Some(BayerPattern::Rggb)
        );
        assert_eq!(
            BayerPattern::from_color_id(ColorId::BayerGrbg),
            Some(BayerPattern::Grbg)
        );
        assert_eq!(
            BayerPattern::from_color_id(ColorId::BayerGbrg),
            Some(BayerPattern::Gbrg)
        );
        assert_eq!(
            BayerPattern::from_color_id(ColorId::BayerBggr),
            Some(BayerPattern::Bggr)
        );
        assert_eq!(BayerPattern::from_color_id(ColorId::Mono), None);
        assert_eq!(BayerPattern::from_color_id(ColorId::Rgb), None);
    }

    #[test]
    fn debayer_2x2_rggb() {
        let raw: Vec<u8> = vec![
            100, 50, // R, G
            60, 200, // G, B
        ];

        let rgb = debayer_bilinear_u8(&raw, 2, 2, BayerPattern::Rggb).unwrap();

        assert_eq!(rgb.len(), 12);
        // Index layout: [R0,G0,B0, R1,G1,B1, R2,G2,B2, R3,G3,B3] for pixels (0,0), (1,0), (0,1), (1,1)
        assert_eq!(rgb[0], 100); // R at (0,0) - original R
        assert_eq!(rgb[11], 200); // B at (1,1) - original B (index 3*3+2=11)
    }

    #[test]
    fn debayer_4x4_output_size() {
        let raw: Vec<u8> = vec![0u8; 16];
        let rgb = debayer_bilinear_u8(&raw, 4, 4, BayerPattern::Rggb).unwrap();
        assert_eq!(rgb.len(), 48); // 4*4*3
    }

    #[test]
    fn debayer_u16() {
        let raw: Vec<u16> = vec![1000, 500, 600, 2000];

        let rgb = debayer_bilinear_u16(&raw, 2, 2, BayerPattern::Rggb).unwrap();

        assert_eq!(rgb.len(), 12);
        assert_eq!(rgb[0], 1000); // R at (0,0)
        assert_eq!(rgb[11], 2000); // B at (1,1) - index 3*3+2=11
    }

    #[test]
    fn debayer_2x2_rggb_exact_output() {
        let rgb = debayer_bilinear_u8(&[0, 100, 100, 0], 2, 2, BayerPattern::Rggb).unwrap();
        assert_eq!(rgb, [0, 50, 50, 50, 100, 50, 50, 100, 50, 50, 50, 0]);
    }

    #[test]
    fn debayer_3x3_every_pattern_exact_output() {
        let expected: [[u8; 27]; 4] = [
            [
                11, 48, 85, 48, 48, 103, 85, 103, 122, 122, 122, 140, 96, 96, 159, 70, 196, 177,
                233, 151, 133, 144, 19, 89, 56, 81, 107,
            ],
            [
                85, 48, 11, 103, 48, 48, 122, 103, 85, 140, 122, 122, 159, 96, 96, 177, 196, 70,
                133, 151, 233, 89, 19, 144, 107, 81, 56,
            ],
            [
                29, 11, 66, 48, 75, 103, 66, 85, 140, 77, 131, 122, 33, 159, 159, 52, 124, 196,
                126, 233, 177, 19, 116, 151, 37, 56, 126,
            ],
            [
                66, 11, 29, 103, 75, 48, 140, 85, 66, 122, 131, 77, 159, 159, 33, 196, 124, 52,
                177, 233, 126, 151, 116, 19, 126, 56, 37,
            ],
        ];
        for (pattern, want) in PATTERNS.into_iter().zip(expected) {
            let rgb = debayer_bilinear_u8(&mosaic_u8(9), 3, 3, pattern).unwrap();
            assert_eq!(rgb, want, "{pattern:?}");
        }
    }

    #[test]
    fn debayer_one_pixel_wide_and_tall() {
        let column = debayer_bilinear_u8(&mosaic_u8(3), 1, 3, BayerPattern::Rggb).unwrap();
        assert_eq!(column, [11, 20, 29, 48, 48, 48, 85, 75, 66]);
        let row = debayer_bilinear_u8(&mosaic_u8(3), 3, 1, BayerPattern::Gbrg).unwrap();
        assert_eq!(row, [11, 11, 29, 48, 48, 48, 85, 85, 66]);
        let single = debayer_bilinear_u16(&[0], 1, 1, BayerPattern::Bggr).unwrap();
        assert_eq!(single, [0, 0, 0]);
    }

    #[test]
    fn debayer_bggr_pattern() {
        let raw: Vec<u8> = vec![
            200, 60, // B, G
            50, 100, // G, R
        ];

        let rgb = debayer_bilinear_u8(&raw, 2, 2, BayerPattern::Bggr).unwrap();

        assert_eq!(rgb[2], 200); // B at (0,0) - original B
        assert_eq!(rgb[9], 100); // R at (1,1) - original R
    }

    #[test]
    fn debayer_u16_edge_rows_match_the_bilinear_reference() {
        let raw: Vec<u16> = (0..80).map(|i| ((i * 7919) % 65521) as u16).collect();
        let rgb = debayer_bilinear_u16(&raw, 20, 4, BayerPattern::Rggb).unwrap();
        assert_eq!(
            &rgb[0..12],
            &[0, 8814, 17628, 7919, 7919, 21588, 15838, 22672, 29507, 23757, 23757, 37426]
        );
        assert_eq!(
            &rgb[180..192],
            &[35584, 16493, 20452, 27123, 33957, 24412, 18662, 32331, 32331, 26581, 33415, 40250]
        );
    }

    #[test]
    fn debayer_short_buffer_is_rejected() {
        let result = debayer_bilinear_u8(&[1, 2, 3], 4, 4, BayerPattern::Rggb);
        assert!(matches!(result, Err(ImageError::InvalidDimensions(_))));
        let result = debayer_bilinear_u8(&[0; 17], 4, 4, BayerPattern::Rggb);
        assert!(matches!(result, Err(ImageError::InvalidDimensions(_))));
    }

    #[test]
    fn debayer_zero_width_is_rejected() {
        let result = debayer_bilinear_u16(&[], 0, 3, BayerPattern::Rggb);
        assert!(matches!(result, Err(ImageError::InvalidDimensions(_))));
        let result = debayer_bilinear_u16(&[], 3, 0, BayerPattern::Rggb);
        assert!(matches!(result, Err(ImageError::InvalidDimensions(_))));
    }

    #[test]
    fn debayer_overflowing_dimensions_are_rejected() {
        let result = debayer_bilinear_u8(&[0; 4], usize::MAX, 2, BayerPattern::Rggb);
        assert!(matches!(result, Err(ImageError::InvalidDimensions(_))));
    }

    #[test]
    fn colour_at_follows_the_pattern_name() {
        use Cfa::{Blue, Green, Red};
        let tile = |p: BayerPattern| {
            [
                p.colour_at(0, 0),
                p.colour_at(1, 0),
                p.colour_at(0, 1),
                p.colour_at(1, 1),
            ]
        };
        assert_eq!(tile(BayerPattern::Rggb), [Red, Green, Green, Blue]);
        assert_eq!(tile(BayerPattern::Grbg), [Green, Red, Blue, Green]);
        assert_eq!(tile(BayerPattern::Gbrg), [Green, Blue, Red, Green]);
        assert_eq!(tile(BayerPattern::Bggr), [Blue, Green, Green, Red]);
        assert_eq!(BayerPattern::Rggb.colour_at(3, 5), Blue);
    }

    #[test]
    fn row_sites_split_green_by_its_row() {
        use Site::{Blue, GreenInBlueRow, GreenInRedRow, Red};
        assert_eq!(BayerPattern::Rggb.row_sites(0), [Red, GreenInRedRow]);
        assert_eq!(BayerPattern::Rggb.row_sites(1), [GreenInBlueRow, Blue]);
        assert_eq!(BayerPattern::Grbg.row_sites(2), [GreenInRedRow, Red]);
        assert_eq!(BayerPattern::Gbrg.row_sites(3), [Red, GreenInRedRow]);
        assert_eq!(BayerPattern::Bggr.row_sites(0), [Blue, GreenInBlueRow]);
    }

    #[test]
    fn debayer_larger_image() {
        let raw: Vec<u8> = (0..64).collect();
        let rgb = debayer_bilinear_u8(&raw, 8, 8, BayerPattern::Rggb).unwrap();
        assert_eq!(rgb.len(), 192); // 8*8*3
    }
}
