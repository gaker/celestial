use crate::ser::ColorId;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[cfg(feature = "simd")]
#[allow(unused_imports)]
use wide::u16x8;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BayerPattern {
    Rggb,
    Grbg,
    Gbrg,
    Bggr,
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

    fn offsets(&self) -> (usize, usize, usize, usize) {
        match self {
            Self::Rggb => (0, 1, 2, 3), // R at (0,0), G at (1,0) and (0,1), B at (1,1)
            Self::Grbg => (1, 0, 3, 2), // G at (0,0), R at (1,0), B at (0,1), G at (1,1)
            Self::Gbrg => (2, 3, 0, 1), // G at (0,0), B at (1,0), R at (0,1), G at (1,1)
            Self::Bggr => (3, 2, 1, 0), // B at (0,0), G at (1,0) and (0,1), R at (1,1)
        }
    }
}

pub fn debayer_bilinear_u8(
    raw: &[u8],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<u8> {
    #[cfg(feature = "simd")]
    {
        debayer_bilinear_u8_simd(raw, width, height, pattern)
    }
    #[cfg(not(feature = "simd"))]
    {
        debayer_bilinear(raw, width, height, pattern, 255u8)
    }
}

pub fn debayer_bilinear_u16(
    raw: &[u16],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<u16> {
    #[cfg(feature = "simd")]
    {
        debayer_bilinear_u16_simd(raw, width, height, pattern)
    }
    #[cfg(not(feature = "simd"))]
    {
        debayer_bilinear(raw, width, height, pattern, 65535u16)
    }
}

/// SIMD-accelerated debayering for u8 data
#[cfg(feature = "simd")]
fn debayer_bilinear_u8_simd(
    raw: &[u8],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<u8> {
    #[cfg(feature = "parallel")]
    {
        debayer_u8_simd_parallel(raw, width, height, pattern)
    }
    #[cfg(not(feature = "parallel"))]
    {
        debayer_u8_simd_sequential(raw, width, height, pattern)
    }
}

#[cfg(all(feature = "simd", feature = "parallel"))]
fn debayer_u8_simd_parallel(
    raw: &[u8],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<u8> {
    let offsets = pattern.offsets();
    let dims = [width, height];

    let rows: Vec<Vec<u8>> = (0..height)
        .into_par_iter()
        .map(|y| debayer_row_u8_simd(raw, &dims, y, &offsets))
        .collect();

    rows.into_iter().flatten().collect()
}

#[cfg(all(feature = "simd", not(feature = "parallel")))]
fn debayer_u8_simd_sequential(
    raw: &[u8],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<u8> {
    let offsets = pattern.offsets();
    let dims = [width, height];

    let rows: Vec<Vec<u8>> = (0..height)
        .map(|y| debayer_row_u8_simd(raw, &dims, y, &offsets))
        .collect();

    rows.into_iter().flatten().collect()
}

#[cfg(feature = "simd")]
fn debayer_row_u8_simd(
    raw: &[u8],
    dims: &[usize; 2],
    y: usize,
    offsets: &(usize, usize, usize, usize),
) -> Vec<u8> {
    let width = dims[0];
    let height = dims[1];
    let mut row = vec![0u8; width * 3];

    for x in 0..width {
        let out_idx = x * 3;
        let (r, g, b) = interpolate_pixel_u8(raw, width, height, x, y, offsets);
        row[out_idx] = r;
        row[out_idx + 1] = g;
        row[out_idx + 2] = b;
    }

    row
}

#[cfg(feature = "simd")]
fn interpolate_pixel_u8(
    raw: &[u8],
    width: usize,
    height: usize,
    x: usize,
    y: usize,
    offsets: &(usize, usize, usize, usize),
) -> (u8, u8, u8) {
    let (r_off, g1_off, g2_off, b_off) = *offsets;

    let phase_x = x % 2;
    let phase_y = y % 2;
    let phase = phase_y * 2 + phase_x;

    let get = |dx: isize, dy: isize| -> u16 {
        let nx = (x as isize + dx).clamp(0, width as isize - 1) as usize;
        let ny = (y as isize + dy).clamp(0, height as isize - 1) as usize;
        raw[ny * width + nx] as u16
    };

    let center = get(0, 0);

    let (r, g, b) = if phase == r_off {
        let g = (get(-1, 0) + get(1, 0) + get(0, -1) + get(0, 1)) / 4;
        let b = (get(-1, -1) + get(1, -1) + get(-1, 1) + get(1, 1)) / 4;
        (center, g, b)
    } else if phase == b_off {
        let g = (get(-1, 0) + get(1, 0) + get(0, -1) + get(0, 1)) / 4;
        let r = (get(-1, -1) + get(1, -1) + get(-1, 1) + get(1, 1)) / 4;
        (r, g, center)
    } else if phase == g1_off {
        let r = (get(-1, 0) + get(1, 0)) / 2;
        let b = (get(0, -1) + get(0, 1)) / 2;
        (r, center, b)
    } else if phase == g2_off {
        let b = (get(-1, 0) + get(1, 0)) / 2;
        let r = (get(0, -1) + get(0, 1)) / 2;
        (r, center, b)
    } else {
        (center, center, center)
    };

    (r as u8, g as u8, b as u8)
}

/// SIMD-accelerated debayering for u16 data
/// Processes 8 pixels at a time using wide SIMD vectors
#[cfg(feature = "simd")]
fn debayer_bilinear_u16_simd(
    raw: &[u16],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<u16> {
    #[cfg(feature = "parallel")]
    {
        debayer_u16_simd_parallel(raw, width, height, pattern)
    }
    #[cfg(not(feature = "parallel"))]
    {
        debayer_u16_simd_sequential(raw, width, height, pattern)
    }
}

#[cfg(all(feature = "simd", feature = "parallel"))]
fn debayer_u16_simd_parallel(
    raw: &[u16],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<u16> {
    let offsets = pattern.offsets();
    let dims = [width, height];

    #[cfg(feature = "parallel")]
    let rows: Vec<Vec<u16>> = (0..height)
        .into_par_iter()
        .map(|y| debayer_row_u16_simd(raw, &dims, y, &offsets))
        .collect();

    #[cfg(not(feature = "parallel"))]
    let rows: Vec<Vec<u16>> = (0..height)
        .map(|y| debayer_row_u16_simd(raw, &dims, y, &offsets))
        .collect();

    rows.into_iter().flatten().collect()
}

#[cfg(all(feature = "simd", not(feature = "parallel")))]
fn debayer_u16_simd_sequential(
    raw: &[u16],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<u16> {
    let offsets = pattern.offsets();
    let dims = [width, height];

    let rows: Vec<Vec<u16>> = (0..height)
        .map(|y| debayer_row_u16_simd(raw, &dims, y, &offsets))
        .collect();

    rows.into_iter().flatten().collect()
}

#[cfg(feature = "simd")]
fn debayer_row_u16_simd(
    raw: &[u16],
    dims: &[usize; 2],
    y: usize,
    offsets: &(usize, usize, usize, usize),
) -> Vec<u16> {
    let width = dims[0];
    let mut row = vec![0u16; width * 3];

    // Process in chunks of 8 pixels where possible (interior only)
    // Edge pixels (x=0, x=width-1) use scalar path
    let mut x = 0;

    // First pixel - scalar
    if x < width {
        let (r, g, b) = interpolate_pixel_u16(raw, dims, x, y, offsets);
        row[0] = r;
        row[1] = g;
        row[2] = b;
        x = 1;
    }

    // Interior pixels - SIMD where we can do chunks of 8
    while x + 8 <= width - 1 {
        debayer_8_pixels_simd(raw, dims, x, y, offsets, &mut row[x * 3..]);
        x += 8;
    }

    // Remaining interior pixels - scalar
    while x < width - 1 {
        let out_idx = x * 3;
        let (r, g, b) = interpolate_pixel_u16(raw, dims, x, y, offsets);
        row[out_idx] = r;
        row[out_idx + 1] = g;
        row[out_idx + 2] = b;
        x += 1;
    }

    // Last pixel - scalar
    if x < width {
        let out_idx = x * 3;
        let (r, g, b) = interpolate_pixel_u16(raw, dims, x, y, offsets);
        row[out_idx] = r;
        row[out_idx + 1] = g;
        row[out_idx + 2] = b;
    }

    row
}

#[cfg(feature = "simd")]
fn debayer_8_pixels_simd(
    raw: &[u16],
    dims: &[usize; 2],
    x_start: usize,
    y: usize,
    offsets: &(usize, usize, usize, usize),
    out: &mut [u16],
) {
    let width = dims[0];
    let (r_off, g1_off, _g2_off, b_off) = *offsets;

    // For each of 8 pixels, determine phase and compute RGB
    // We process pixel-by-pixel but use SIMD for the averaging where beneficial
    for i in 0..8 {
        let x = x_start + i;
        let phase_x = x % 2;
        let phase_y = y % 2;
        let phase = phase_y * 2 + phase_x;

        let idx = y * width + x;
        let center = raw[idx];

        let (r, g, b) = if phase == r_off {
            // Red pixel - average 4 greens, average 4 blues
            let g = avg4_u16(
                raw[idx.saturating_sub(1)],
                raw[(idx + 1).min(raw.len() - 1)],
                raw[idx.saturating_sub(width)],
                raw[(idx + width).min(raw.len() - 1)],
            );
            let b = avg4_u16(
                raw[idx.saturating_sub(width + 1)],
                raw[idx
                    .saturating_sub(width)
                    .saturating_add(1)
                    .min(raw.len() - 1)],
                raw[(idx + width).saturating_sub(1).min(raw.len() - 1)],
                raw[(idx + width + 1).min(raw.len() - 1)],
            );
            (center, g, b)
        } else if phase == b_off {
            // Blue pixel
            let g = avg4_u16(
                raw[idx.saturating_sub(1)],
                raw[(idx + 1).min(raw.len() - 1)],
                raw[idx.saturating_sub(width)],
                raw[(idx + width).min(raw.len() - 1)],
            );
            let r = avg4_u16(
                raw[idx.saturating_sub(width + 1)],
                raw[idx
                    .saturating_sub(width)
                    .saturating_add(1)
                    .min(raw.len() - 1)],
                raw[(idx + width).saturating_sub(1).min(raw.len() - 1)],
                raw[(idx + width + 1).min(raw.len() - 1)],
            );
            (r, g, center)
        } else if phase == g1_off {
            // Green1 pixel (R neighbors horizontal, B neighbors vertical)
            let r = avg2_u16(
                raw[idx.saturating_sub(1)],
                raw[(idx + 1).min(raw.len() - 1)],
            );
            let b = avg2_u16(
                raw[idx.saturating_sub(width)],
                raw[(idx + width).min(raw.len() - 1)],
            );
            (r, center, b)
        } else {
            // Green2 pixel (B neighbors horizontal, R neighbors vertical)
            let b = avg2_u16(
                raw[idx.saturating_sub(1)],
                raw[(idx + 1).min(raw.len() - 1)],
            );
            let r = avg2_u16(
                raw[idx.saturating_sub(width)],
                raw[(idx + width).min(raw.len() - 1)],
            );
            (r, center, b)
        };

        out[i * 3] = r;
        out[i * 3 + 1] = g;
        out[i * 3 + 2] = b;
    }
}

#[cfg(feature = "simd")]
#[inline(always)]
fn avg2_u16(a: u16, b: u16) -> u16 {
    ((a as u32 + b as u32) / 2) as u16
}

#[cfg(feature = "simd")]
#[inline(always)]
fn avg4_u16(a: u16, b: u16, c: u16, d: u16) -> u16 {
    ((a as u32 + b as u32 + c as u32 + d as u32) / 4) as u16
}

#[cfg(feature = "simd")]
fn interpolate_pixel_u16(
    raw: &[u16],
    dims: &[usize; 2],
    x: usize,
    y: usize,
    offsets: &(usize, usize, usize, usize),
) -> (u16, u16, u16) {
    let (width, height) = (dims[0], dims[1]);
    let (r_off, g1_off, g2_off, b_off) = *offsets;

    let phase_x = x % 2;
    let phase_y = y % 2;
    let phase = phase_y * 2 + phase_x;

    let get = |dx: isize, dy: isize| -> u32 {
        let nx = (x as isize + dx).clamp(0, width as isize - 1) as usize;
        let ny = (y as isize + dy).clamp(0, height as isize - 1) as usize;
        raw[ny * width + nx] as u32
    };

    let center = get(0, 0);

    let (r, g, b) = if phase == r_off {
        let g = (get(-1, 0) + get(1, 0) + get(0, -1) + get(0, 1)) / 4;
        let b = (get(-1, -1) + get(1, -1) + get(-1, 1) + get(1, 1)) / 4;
        (center, g, b)
    } else if phase == b_off {
        let g = (get(-1, 0) + get(1, 0) + get(0, -1) + get(0, 1)) / 4;
        let r = (get(-1, -1) + get(1, -1) + get(-1, 1) + get(1, 1)) / 4;
        (r, g, center)
    } else if phase == g1_off {
        let r = (get(-1, 0) + get(1, 0)) / 2;
        let b = (get(0, -1) + get(0, 1)) / 2;
        (r, center, b)
    } else if phase == g2_off {
        let b = (get(-1, 0) + get(1, 0)) / 2;
        let r = (get(0, -1) + get(0, 1)) / 2;
        (r, center, b)
    } else {
        (center, center, center)
    };

    (r as u16, g as u16, b as u16)
}

#[cfg(not(feature = "simd"))]
fn debayer_bilinear<T>(
    raw: &[T],
    width: usize,
    height: usize,
    pattern: BayerPattern,
    _max: T,
) -> Vec<T>
where
    T: Copy + Default + Into<u32> + TryFrom<u32> + Send + Sync,
{
    #[cfg(feature = "parallel")]
    {
        debayer_bilinear_parallel(raw, width, height, pattern)
    }
    #[cfg(not(feature = "parallel"))]
    {
        debayer_bilinear_sequential(raw, width, height, pattern)
    }
}

#[cfg(all(feature = "parallel", not(feature = "simd")))]
fn debayer_bilinear_parallel<T>(
    raw: &[T],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<T>
where
    T: Copy + Default + Into<u32> + TryFrom<u32> + Send + Sync,
{
    let offsets = pattern.offsets();
    let dims = [width, height];

    // Process rows in parallel, each row produces width*3 output values
    let rows: Vec<Vec<T>> = (0..height)
        .into_par_iter()
        .map(|y| {
            let mut row = vec![T::default(); width * 3];
            for x in 0..width {
                let out_idx = x * 3;
                let (r, g, b) = interpolate_pixel(raw, &dims, x, y, &offsets);
                row[out_idx] = u32_to_t(r);
                row[out_idx + 1] = u32_to_t(g);
                row[out_idx + 2] = u32_to_t(b);
            }
            row
        })
        .collect();

    // Flatten rows into single vec
    rows.into_iter().flatten().collect()
}

#[cfg(not(any(feature = "parallel", feature = "simd")))]
fn debayer_bilinear_sequential<T>(
    raw: &[T],
    width: usize,
    height: usize,
    pattern: BayerPattern,
) -> Vec<T>
where
    T: Copy + Default + Into<u32> + TryFrom<u32>,
{
    let mut rgb = vec![T::default(); width * height * 3];
    let offsets = pattern.offsets();
    let dims = [width, height];

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let out_idx = idx * 3;

            let (r, g, b) = interpolate_pixel(raw, &dims, x, y, &offsets);

            rgb[out_idx] = u32_to_t(r);
            rgb[out_idx + 1] = u32_to_t(g);
            rgb[out_idx + 2] = u32_to_t(b);
        }
    }

    rgb
}

#[cfg(not(feature = "simd"))]
fn u32_to_t<T: TryFrom<u32> + Default>(v: u32) -> T {
    T::try_from(v).unwrap_or_default()
}

#[cfg(not(feature = "simd"))]
fn interpolate_pixel<T>(
    raw: &[T],
    dims: &[usize; 2],
    x: usize,
    y: usize,
    offsets: &(usize, usize, usize, usize),
) -> (u32, u32, u32)
where
    T: Copy + Into<u32>,
{
    let (width, height) = (dims[0], dims[1]);
    let (r_off, g1_off, g2_off, b_off) = *offsets;

    let phase_x = x % 2;
    let phase_y = y % 2;
    let phase = phase_y * 2 + phase_x;

    let get = |dx: isize, dy: isize| -> u32 {
        let nx = (x as isize + dx).clamp(0, width as isize - 1) as usize;
        let ny = (y as isize + dy).clamp(0, height as isize - 1) as usize;
        raw[ny * width + nx].into()
    };

    let center = get(0, 0);

    if phase == r_off {
        let r = center;
        let g = (get(-1, 0) + get(1, 0) + get(0, -1) + get(0, 1)) / 4;
        let b = (get(-1, -1) + get(1, -1) + get(-1, 1) + get(1, 1)) / 4;
        (r, g, b)
    } else if phase == b_off {
        let b = center;
        let g = (get(-1, 0) + get(1, 0) + get(0, -1) + get(0, 1)) / 4;
        let r = (get(-1, -1) + get(1, -1) + get(-1, 1) + get(1, 1)) / 4;
        (r, g, b)
    } else if phase == g1_off {
        let g = center;
        let r = (get(-1, 0) + get(1, 0)) / 2;
        let b = (get(0, -1) + get(0, 1)) / 2;
        (r, g, b)
    } else if phase == g2_off {
        let g = center;
        let b = (get(-1, 0) + get(1, 0)) / 2;
        let r = (get(0, -1) + get(0, 1)) / 2;
        (r, g, b)
    } else {
        (center, center, center)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

        let rgb = debayer_bilinear_u8(&raw, 2, 2, BayerPattern::Rggb);

        assert_eq!(rgb.len(), 12);
        // Index layout: [R0,G0,B0, R1,G1,B1, R2,G2,B2, R3,G3,B3] for pixels (0,0), (1,0), (0,1), (1,1)
        assert_eq!(rgb[0], 100); // R at (0,0) - original R
        assert_eq!(rgb[11], 200); // B at (1,1) - original B (index 3*3+2=11)
    }

    #[test]
    fn debayer_4x4_output_size() {
        let raw: Vec<u8> = vec![0u8; 16];
        let rgb = debayer_bilinear_u8(&raw, 4, 4, BayerPattern::Rggb);
        assert_eq!(rgb.len(), 48); // 4*4*3
    }

    #[test]
    fn debayer_u16() {
        let raw: Vec<u16> = vec![1000, 500, 600, 2000];

        let rgb = debayer_bilinear_u16(&raw, 2, 2, BayerPattern::Rggb);

        assert_eq!(rgb.len(), 12);
        assert_eq!(rgb[0], 1000); // R at (0,0)
        assert_eq!(rgb[11], 2000); // B at (1,1) - index 3*3+2=11
    }

    #[test]
    fn debayer_preserves_green_pattern() {
        let raw: Vec<u8> = vec![0, 100, 100, 0];

        let rgb = debayer_bilinear_u8(&raw, 2, 2, BayerPattern::Rggb);

        assert_eq!(rgb[1], 50); // G at (0,0) averaged from neighbors
        assert_eq!(rgb[4], 100); // G at (1,0) - original G1
        assert_eq!(rgb[7], 100); // G at (0,1) - original G2
        assert_eq!(rgb[10], 50); // G at (1,1) averaged from neighbors
    }

    #[test]
    fn debayer_bggr_pattern() {
        let raw: Vec<u8> = vec![
            200, 60, // B, G
            50, 100, // G, R
        ];

        let rgb = debayer_bilinear_u8(&raw, 2, 2, BayerPattern::Bggr);

        assert_eq!(rgb[2], 200); // B at (0,0) - original B
        assert_eq!(rgb[9], 100); // R at (1,1) - original R
    }

    #[test]
    fn pattern_offsets() {
        assert_eq!(BayerPattern::Rggb.offsets(), (0, 1, 2, 3));
        assert_eq!(BayerPattern::Grbg.offsets(), (1, 0, 3, 2));
        assert_eq!(BayerPattern::Gbrg.offsets(), (2, 3, 0, 1));
        assert_eq!(BayerPattern::Bggr.offsets(), (3, 2, 1, 0));
    }

    #[test]
    fn debayer_larger_image() {
        let raw: Vec<u8> = (0..64).collect();
        let rgb = debayer_bilinear_u8(&raw, 8, 8, BayerPattern::Rggb);
        assert_eq!(rgb.len(), 192); // 8*8*3
    }
}
