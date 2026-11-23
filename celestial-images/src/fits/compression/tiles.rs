use crate::fits::errors::{FitsError, Result};
use crate::fits::util::image::shape_text;

// The tiles of FITS 4.0 §10.1.2, numbered in the order of their first pixel,
// which is the order of the table rows. Axis 1 varies fastest.
pub(super) struct TileGrid {
    axes: Vec<usize>,
    tile: Vec<usize>,
    counts: Vec<usize>,
    tile_count: usize,
    pixel_count: usize,
}

// One tile: its first pixel and its size along each axis. Tiles at the end of
// an axis are cut short when ZTILEn does not divide ZNAXISn.
pub(super) struct Region {
    origin: Vec<usize>,
    shape: Vec<usize>,
}

impl TileGrid {
    pub(super) fn new(axes: Vec<usize>, tile: Vec<usize>) -> Result<Self> {
        if axes.len() != tile.len() || axes.iter().chain(&tile).any(|&n| n == 0) {
            return Err(FitsError::InvalidFormat(format!(
                "tile shape {} does not fit image shape {}",
                shape_text(&tile),
                shape_text(&axes)
            )));
        }
        let counts: Vec<usize> = axes
            .iter()
            .zip(&tile)
            .map(|(a, t)| a.div_ceil(*t))
            .collect();
        let too_large =
            || FitsError::InvalidFormat(format!("image shape {} is too large", shape_text(&axes)));
        let tile_count = checked_product(&counts).ok_or_else(too_large)?;
        let pixel_count = checked_product(&axes).ok_or_else(too_large)?;
        Ok(Self {
            axes,
            tile,
            counts,
            tile_count,
            pixel_count,
        })
    }

    pub(super) fn tile_count(&self) -> usize {
        self.tile_count
    }

    pub(super) fn pixel_count(&self) -> usize {
        self.pixel_count
    }

    pub(super) fn region(&self, index: usize) -> Region {
        let mut rest = index;
        let mut origin = Vec::with_capacity(self.axes.len());
        let mut shape = Vec::with_capacity(self.axes.len());
        for ((&axis, &tile), &count) in self.axes.iter().zip(&self.tile).zip(&self.counts) {
            let start = (rest % count) * tile;
            rest /= count;
            origin.push(start);
            shape.push(tile.min(axis - start));
        }
        Region { origin, shape }
    }

    // Copies a tile, stored axis 1 fastest, into the image at its place.
    pub(super) fn place(
        &self,
        region: &Region,
        width: usize,
        tile: &[u8],
        image: &mut [u8],
    ) -> Result<()> {
        self.check_sizes(region, width, tile.len(), image.len())?;
        let run = region.shape[0] * width;
        for (start, src) in self.line_starts(region).zip(tile.chunks_exact(run)) {
            image[start * width..start * width + run].copy_from_slice(src);
        }
        Ok(())
    }

    // Copies the pixels of a tile out of the image, axis 1 fastest: the
    // inverse of `place`.
    pub(super) fn extract(&self, region: &Region, width: usize, image: &[u8]) -> Result<Vec<u8>> {
        let len = region.pixel_count() * width;
        self.check_sizes(region, width, len, image.len())?;
        let run = region.shape[0] * width;
        let mut tile = Vec::with_capacity(len);
        for start in self.line_starts(region) {
            tile.extend_from_slice(&image[start * width..start * width + run]);
        }
        Ok(tile)
    }

    fn check_sizes(&self, region: &Region, width: usize, tile: usize, image: usize) -> Result<()> {
        if tile != region.pixel_count() * width || image != self.pixel_count * width {
            return Err(FitsError::InvalidFormat(format!(
                "a tile of {} bytes does not fit its {} region",
                tile,
                shape_text(&region.shape)
            )));
        }
        Ok(())
    }

    // The image offset, in pixels, of the first pixel of each line of the
    // tile, in the order the tile stores its lines.
    fn line_starts<'r>(&'r self, region: &'r Region) -> impl Iterator<Item = usize> + 'r {
        let lines = region.pixel_count() / region.shape[0];
        let mut position = vec![0usize; self.axes.len()];
        (0..lines).map(move |line| {
            if line > 0 {
                step(&mut position, &region.shape);
            }
            self.offset(region, &position)
        })
    }

    fn offset(&self, region: &Region, position: &[usize]) -> usize {
        let mut offset = 0;
        let mut stride = 1;
        for ((&axis, &origin), &pos) in self.axes.iter().zip(&region.origin).zip(position) {
            offset += (origin + pos) * stride;
            stride *= axis;
        }
        offset
    }
}

impl Region {
    pub(super) fn pixel_count(&self) -> usize {
        self.shape.iter().product()
    }
}

// Advances `position` to the next line of the tile: axis 2 first, carrying
// into the higher axes. Axis 1 stays 0 because each line is copied whole.
fn step(position: &mut [usize], shape: &[usize]) {
    for (pos, &len) in position.iter_mut().zip(shape).skip(1) {
        *pos += 1;
        if *pos < len {
            return;
        }
        *pos = 0;
    }
}

fn checked_product(values: &[usize]) -> Option<usize> {
    values.iter().try_fold(1usize, |acc, &v| acc.checked_mul(v))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(axes: &[usize], tile: &[usize]) -> TileGrid {
        TileGrid::new(axes.to_vec(), tile.to_vec()).unwrap()
    }

    #[test]
    fn regions_cover_an_uneven_grid() {
        let g = grid(&[5, 3], &[2, 2]);
        assert_eq!(g.tile_count(), 6);
        let shapes: Vec<(Vec<usize>, Vec<usize>)> = (0..6)
            .map(|i| {
                let r = g.region(i);
                (r.origin, r.shape)
            })
            .collect();
        assert_eq!(
            shapes,
            vec![
                (vec![0, 0], vec![2, 2]),
                (vec![2, 0], vec![2, 2]),
                (vec![4, 0], vec![1, 2]),
                (vec![0, 2], vec![2, 1]),
                (vec![2, 2], vec![2, 1]),
                (vec![4, 2], vec![1, 1]),
            ]
        );
    }

    #[test]
    fn place_assembles_tiles_into_the_image() {
        let g = grid(&[3, 2], &[2, 2]);
        let mut image = vec![0u8; 6];
        g.place(&g.region(0), 1, &[1, 2, 4, 5], &mut image).unwrap();
        g.place(&g.region(1), 1, &[3, 6], &mut image).unwrap();
        assert_eq!(image, vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn place_handles_three_axes_and_wide_pixels() {
        let g = grid(&[2, 2, 2], &[2, 1, 2]);
        let mut image = vec![0u8; 16];
        let tile0: Vec<u8> = [0u8, 1, 4, 5].iter().flat_map(|&v| [0, v]).collect();
        let tile1: Vec<u8> = [2u8, 3, 6, 7].iter().flat_map(|&v| [0, v]).collect();
        g.place(&g.region(0), 2, &tile0, &mut image).unwrap();
        g.place(&g.region(1), 2, &tile1, &mut image).unwrap();
        let values: Vec<u8> = image.chunks(2).map(|c| c[1]).collect();
        assert_eq!(values, vec![0, 1, 2, 3, 4, 5, 6, 7]);
    }

    #[test]
    fn extract_reverses_place_on_uneven_tiles() {
        let g = grid(&[5, 3, 2], &[2, 2, 1]);
        let image: Vec<u8> = (0..60).collect();
        let mut rebuilt = vec![0u8; 60];
        for i in 0..g.tile_count() {
            let region = g.region(i);
            let tile = g.extract(&region, 2, &image).unwrap();
            g.place(&region, 2, &tile, &mut rebuilt).unwrap();
        }
        assert_eq!(rebuilt, image);
        assert_eq!(
            g.extract(&g.region(2), 1, &image[..30]).unwrap(),
            vec![4, 9]
        );
        assert!(g.extract(&g.region(0), 1, &image).is_err());
    }

    #[test]
    fn tiles_larger_than_the_image_cover_it_once() {
        let g = grid(&[3, 2], &[10, 10]);
        assert_eq!(g.tile_count(), 1);
        assert_eq!(g.region(0).shape, vec![3, 2]);
    }

    #[test]
    fn new_rejects_zero_sizes_and_overflow() {
        assert!(TileGrid::new(vec![0, 2], vec![1, 1]).is_err());
        assert!(TileGrid::new(vec![2, 2], vec![0, 1]).is_err());
        assert!(TileGrid::new(vec![2], vec![1, 1]).is_err());
        assert!(TileGrid::new(vec![1 << 40, 1 << 40], vec![1, 1]).is_err());
    }

    fn message(result: Result<TileGrid>) -> String {
        match result {
            Err(FitsError::InvalidFormat(m)) => m,
            Err(other) => panic!("{other}"),
            Ok(_) => panic!("grid built"),
        }
    }

    #[test]
    fn errors_write_shapes_as_axis_lengths() {
        assert_eq!(
            message(TileGrid::new(vec![2], vec![1, 1])),
            "tile shape 1x1 does not fit image shape 2"
        );
        assert_eq!(
            message(TileGrid::new(vec![1 << 40, 1 << 40], vec![1, 1])),
            "image shape 1099511627776x1099511627776 is too large"
        );
    }
}
