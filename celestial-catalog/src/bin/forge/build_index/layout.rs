// Writes the catalog header and pixel table into a memory-mapped output file, using the
// library's byte writers for both.

use crate::intermediate::RECORD_SIZE;
use anyhow::Context;
use celestial_catalog::query::catalog::{pixel_entry_bytes, CatalogHeader};
use memmap2::{Mmap, MmapMut};
use std::fs::{File, OpenOptions};
use std::path::Path;

// Where the sections of a catalog file start, in bytes from the start of the file.
#[derive(Debug, Clone, Copy)]
pub(super) struct Layout {
    table_start: usize,
    records_start: usize,
    file_len: usize,
}

impl Layout {
    pub(super) fn new(header: &CatalogHeader) -> anyhow::Result<Self> {
        let table_start = header.to_bytes().len();
        let entry_len = pixel_entry_bytes(0, 0).len();
        let too_big = || {
            format!(
                "A catalog of {} pixels and {} stars does not fit in memory",
                header.npix, header.total_stars
            )
        };
        let table_len = usize::try_from(header.npix)?.checked_mul(entry_len);
        let records_len = usize::try_from(header.total_stars)?.checked_mul(RECORD_SIZE);
        let records_start = table_len.and_then(|len| len.checked_add(table_start));
        let file_len = records_start
            .zip(records_len)
            .and_then(|(start, len)| start.checked_add(len));
        Ok(Self {
            table_start,
            records_start: records_start.with_context(too_big)?,
            file_len: file_len.with_context(too_big)?,
        })
    }

    pub(super) fn records<'a>(&self, map: &'a Mmap) -> anyhow::Result<&'a [[u8; RECORD_SIZE]]> {
        self.check_len(map.len())?;
        Ok(map[self.records_start..].as_chunks().0)
    }

    pub(super) fn records_mut<'a>(
        &self,
        map: &'a mut MmapMut,
    ) -> anyhow::Result<&'a mut [[u8; RECORD_SIZE]]> {
        self.check_len(map.len())?;
        Ok(map[self.records_start..].as_chunks_mut().0)
    }

    fn check_len(&self, map_len: usize) -> anyhow::Result<()> {
        if map_len != self.file_len {
            anyhow::bail!(
                "Catalog file is {map_len} bytes, but its header implies {} bytes",
                self.file_len
            );
        }
        Ok(())
    }
}

// Creates the output file at its full size and writes the header and the table for `counts`,
// whose sum must be the header's total_stars. Returns the map, with the records still zero.
pub(super) fn create_catalog(
    path: &Path,
    header: &CatalogHeader,
    counts: &[u32],
) -> anyhow::Result<(MmapMut, Layout)> {
    let layout = Layout::new(header)?;
    let file = create_file(path, layout.file_len)?;
    // SAFETY: forge created the file above and is its only writer while the map exists.
    let mut map = unsafe { MmapMut::map_mut(&file) }
        .with_context(|| format!("Failed to memory-map {path:?}"))?;
    map[..layout.table_start].copy_from_slice(&header.to_bytes());
    let table = &mut map[layout.table_start..layout.records_start];
    let mut offset = 0u64;
    for (entry, &count) in table
        .chunks_exact_mut(pixel_entry_bytes(0, 0).len())
        .zip(counts)
    {
        entry.copy_from_slice(&pixel_entry_bytes(offset, count));
        offset += u64::from(count) * RECORD_SIZE as u64;
    }
    Ok((map, layout))
}

fn create_file(path: &Path, size: usize) -> anyhow::Result<File> {
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(path)
        .with_context(|| format!("Failed to create {path:?}"))?;
    file.set_len(size as u64)
        .with_context(|| format!("Failed to size {path:?} to {size} bytes"))?;
    Ok(file)
}

pub(super) fn map_file(path: &Path) -> anyhow::Result<Mmap> {
    let file = File::open(path).with_context(|| format!("Failed to open {path:?}"))?;
    // SAFETY: the map is read-only, and forge does not modify the file while it is mapped.
    unsafe { Mmap::map(&file) }.with_context(|| format!("Failed to memory-map {path:?}"))
}
