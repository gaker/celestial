// The passes over the merged records: count per pixel, scatter into pixel order, sort each
// pixel by magnitude and optionally keep only the brightest stars of each pixel.

use super::layout::{create_catalog, map_file, Layout};
use super::Grid;
use crate::intermediate::{Header, Kind, HEADER_SIZE, RECORD_SIZE};
use anyhow::Context;
use celestial_catalog::query::catalog::{CatalogHeader, StarRecord};
use celestial_catalog::query::healpix::{ang2pix_nest, reduce_ra_deg};
use memmap2::Mmap;
use std::path::Path;

pub(super) struct MergedInput {
    pub(super) header: Header,
    map: Mmap,
}

impl MergedInput {
    pub(super) fn open(dir: &Path) -> anyhow::Result<Self> {
        let path = dir.join(Kind::Merged.file_name());
        let map = map_file(&path)?;
        let context = || format!("Failed to read merged catalog {path:?}");
        let bytes = map
            .first_chunk::<HEADER_SIZE>()
            .with_context(|| format!("{path:?} is shorter than the {HEADER_SIZE}-byte header"))?;
        let header = Header::parse(bytes, Kind::Merged, map.len() as u64).with_context(context)?;
        Ok(Self { header, map })
    }

    // Header::parse checked that the records fill the rest of the file exactly.
    fn records(&self) -> &[[u8; RECORD_SIZE]] {
        self.map[HEADER_SIZE..].as_chunks().0
    }
}

// Returns the record's pixel and the record with its RA reduced into [0, 360). A record whose
// position or magnitude is not finite, or whose Dec is outside [-90, 90], is an error. A
// non-finite magnitude would sort ahead of or behind every real star in its pixel.
fn place(
    order: u32,
    index: usize,
    bytes: &[u8; RECORD_SIZE],
) -> anyhow::Result<(usize, StarRecord)> {
    let mut star = StarRecord::from_le_bytes(bytes);
    let record = || format!("Merged record {index} (source_id {})", star.source_id);
    if !star.mag.is_finite() {
        anyhow::bail!("{} has magnitude {}", record(), star.mag);
    }
    let pixel = ang2pix_nest(order, star.ra, star.dec)
        .with_context(|| format!("{} has no HEALPix pixel", record()))?;
    star.ra = reduce_ra_deg(star.ra);
    Ok((usize::try_from(pixel)?, star))
}

fn count_stars_per_pixel(input: &MergedInput, grid: &Grid) -> anyhow::Result<Vec<u32>> {
    let mut counts = vec![0u32; grid.npix];
    for (index, bytes) in input.records().iter().enumerate() {
        let (pixel, _) = place(grid.order, index, bytes)?;
        let count = counts
            .get_mut(pixel)
            .with_context(|| format!("Pixel {pixel} is outside the {} pixels", grid.npix))?;
        *count = count
            .checked_add(1)
            .with_context(|| format!("Pixel {pixel} holds more than {} stars", u32::MAX))?;
    }
    Ok(counts)
}

// Writes a catalog of every merged record to `path`, each pixel's records sorted brightest
// first. Returns the number of records in each pixel.
pub(super) fn write_sorted(
    input: &MergedInput,
    grid: &Grid,
    header: &CatalogHeader,
    path: &Path,
) -> anyhow::Result<Vec<u32>> {
    println!("Pass 1: Counting stars per pixel...");
    let counts = count_stars_per_pixel(input, grid)?;
    let (mut map, layout) = create_catalog(path, header, &counts)?;
    let records = layout.records_mut(&mut map)?;
    println!("Pass 2: Scattering stars to output positions...");
    scatter(input, grid.order, &counts, records)?;
    println!("Pass 3: Sorting each pixel by magnitude...");
    for run in pixel_runs(&counts, records) {
        run.sort_by(|a, b| magnitude(a).total_cmp(&magnitude(b)));
    }
    map.flush()
        .with_context(|| format!("Failed to write {path:?}"))?;
    Ok(counts)
}

fn scatter(
    input: &MergedInput,
    order: u32,
    counts: &[u32],
    records: &mut [[u8; RECORD_SIZE]],
) -> anyhow::Result<()> {
    let mut next_slot: Vec<usize> = counts
        .iter()
        .scan(0usize, |start, &count| {
            let slot = *start;
            *start += count as usize;
            Some(slot)
        })
        .collect();
    for (index, bytes) in input.records().iter().enumerate() {
        let (pixel, star) = place(order, index, bytes)?;
        // The counting pass placed every record, so each pixel has a slot left for it.
        let slot = next_slot.get_mut(pixel).and_then(|slot| {
            let record = records.get_mut(*slot)?;
            *slot += 1;
            Some(record)
        });
        *slot.context("Merged records changed between the counting and scatter passes")? =
            star.to_le_bytes();
    }
    Ok(())
}

fn magnitude(record: &[u8; RECORD_SIZE]) -> f32 {
    StarRecord::from_le_bytes(record).mag
}

// Splits the records into each pixel's run. The counts sum to the number of records.
fn pixel_runs<'a>(
    counts: &'a [u32],
    mut records: &'a mut [[u8; RECORD_SIZE]],
) -> impl Iterator<Item = &'a mut [[u8; RECORD_SIZE]]> + 'a {
    counts.iter().map(move |&count| {
        let split = (count as usize).min(records.len());
        let (run, rest) = std::mem::take(&mut records).split_at_mut(split);
        records = rest;
        run
    })
}

// Writes to `output` a copy of the sorted catalog at `sorted` that keeps the first `cap`
// records of each pixel.
pub(super) fn compact_with_cap(
    sorted: &Path,
    output: &Path,
    header: &CatalogHeader,
    counts: &[u32],
    cap: u32,
) -> anyhow::Result<()> {
    let capped: Vec<u32> = counts.iter().map(|&count| count.min(cap)).collect();
    let capped_header = CatalogHeader {
        total_stars: capped.iter().map(|&count| u64::from(count)).sum(),
        ..header.clone()
    };
    let source_map = map_file(sorted)?;
    let source = Layout::new(header)?.records(&source_map)?;
    let (mut map, layout) = create_catalog(output, &capped_header, &capped)?;
    let target = layout.records_mut(&mut map)?;
    if copy_brightest(source, counts, cap, target) != target.len() {
        anyhow::bail!("{sorted:?} holds fewer records than its pixel counts give");
    }
    map.flush()
        .with_context(|| format!("Failed to write {output:?}"))
}

// Copies the first `cap` records of each pixel's run in `source` into `target`, in order.
// Returns the number of records copied.
fn copy_brightest(
    source: &[[u8; RECORD_SIZE]],
    counts: &[u32],
    cap: u32,
    target: &mut [[u8; RECORD_SIZE]],
) -> usize {
    let mut rest = source;
    let kept = counts.iter().flat_map(move |&count| {
        let (run, tail) = rest.split_at(rest.len().min(count as usize));
        rest = tail;
        run.iter().take(cap as usize)
    });
    target
        .iter_mut()
        .zip(kept)
        .map(|(slot, record)| *slot = *record)
        .count()
}
