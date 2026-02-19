//! Build HEALPix-indexed binary catalog from merged catalog.
//!
//! Three-pass memory-mapped algorithm for handling arbitrarily large catalogs:
//! 1. Count stars per HEALPix pixel
//! 2. Scatter stars to their final positions in output file
//! 3. Sort each pixel's stars by magnitude in-place

use crate::cli::{BuildIndexArgs, Cli};
use celestial_catalog::query::healpix::ang2pix_nest;
use memmap2::{Mmap, MmapMut};
use std::fs::{self, File, OpenOptions};
use std::io::Read;
use std::path::Path;
use std::time::Instant;

const RECORD_SIZE: usize = 56;
const MERGED_HEADER_SIZE: usize = 24;
const OUTPUT_HEADER_SIZE: usize = 64;
const PIXEL_ENTRY_SIZE: usize = 16;
const CATALOG_MAGIC: &[u8; 4] = b"CCAT";
const CATALOG_VERSION: u32 = 1;
const EPOCH_J2016: f64 = 2016.0;

struct IndexStats {
    healpix_order: u32,
    nside: u32,
    npix: u64,
    total_stars: u64,
    stars_after_cap: Option<u64>,
    cells_capped: Option<u64>,
    min_stars: u32,
    max_stars: u32,
    mean_stars: f64,
    median_stars: u32,
    empty_pixels: u64,
    file_size: u64,
    elapsed_secs: f64,
}

pub fn run(args: &BuildIndexArgs, cli: &Cli) -> anyhow::Result<()> {
    validate_paths(args)?;
    print_plan(args, cli);
    let start = Instant::now();
    let merged_path = args.workdir.join("merged.bin");
    let (total_stars, mag_limit) = read_merged_header(&merged_path)?;
    let nside = 1u32 << args.healpix_order;
    let npix = 12u64 * (nside as u64) * (nside as u64);

    println!("Memory-mapping input file ({} stars)...", total_stars);
    let input_mmap = mmap_input(&merged_path)?;

    println!("Pass 1: Counting stars per pixel...");
    let counts = count_stars_per_pixel_mmap(&input_mmap, total_stars, args.healpix_order)?;

    println!("Pass 2: Scattering stars to output positions...");
    let temp_path = args.output.with_extension("bin.tmp");
    scatter_stars_to_output(
        &input_mmap,
        &temp_path,
        total_stars,
        args.healpix_order,
        mag_limit,
        &counts,
    )?;

    drop(input_mmap);

    println!("Pass 3: Sorting each pixel by magnitude...");
    sort_pixels_in_place(&temp_path, &counts)?;

    let final_counts = match args.max_per_cell {
        Some(cap) => {
            println!("Pass 4: Compacting to {} max stars per cell...", cap);
            compact_with_cap(
                &temp_path,
                &args.output,
                &counts,
                cap,
                mag_limit,
                args.healpix_order,
            )?;
            fs::remove_file(&temp_path)?;
            counts.iter().map(|&c| c.min(cap)).collect::<Vec<u32>>()
        }
        None => {
            fs::rename(&temp_path, &args.output)?;
            counts.clone()
        }
    };

    let elapsed = start.elapsed().as_secs_f64();
    let stats = compute_stats(
        args.healpix_order,
        nside,
        npix,
        total_stars,
        args.max_per_cell,
        &final_counts,
        &args.output,
        elapsed,
    )?;
    print_stats(&stats);

    println!("\nValidating output...");
    validate_output(&args.output, args.healpix_order)?;
    Ok(())
}

fn validate_paths(args: &BuildIndexArgs) -> anyhow::Result<()> {
    let merged = args.workdir.join("merged.bin");
    if !merged.exists() {
        anyhow::bail!("Merged catalog not found: {:?}", merged);
    }
    if let Some(parent) = args.output.parent() {
        if !parent.exists() {
            fs::create_dir_all(parent)?;
        }
    }
    Ok(())
}

fn print_plan(args: &BuildIndexArgs, cli: &Cli) {
    let nside = 1u32 << args.healpix_order;
    let npix = 12u64 * (nside as u64) * (nside as u64);
    println!("=== Build HEALPix Index ===");
    println!("Working directory: {:?}", args.workdir);
    println!("HEALPix order: {}", args.healpix_order);
    println!("nside: {}", nside);
    println!("npix: {}", npix);
    println!("Output: {:?}", args.output);
    match args.max_per_cell {
        Some(cap) => println!("Max per cell: {}", cap),
        None => println!("Max per cell: unlimited"),
    }
    println!("Verbose: {}", cli.verbose);
    println!();
}

fn read_merged_header(path: &Path) -> anyhow::Result<(u64, f32)> {
    let mut file = File::open(path)?;
    let mut header = [0u8; MERGED_HEADER_SIZE];
    file.read_exact(&mut header)?;
    let magic = &header[0..4];
    if magic != b"MERG" {
        anyhow::bail!("Invalid merged catalog magic: {:?}", magic);
    }
    let total_stars = u64::from_le_bytes(header[8..16].try_into().unwrap());
    let mag_limit = f32::from_le_bytes(header[16..20].try_into().unwrap());
    Ok((total_stars, mag_limit))
}

fn mmap_input(path: &Path) -> anyhow::Result<Mmap> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    Ok(mmap)
}

fn count_stars_per_pixel_mmap(mmap: &Mmap, total: u64, order: u32) -> anyhow::Result<Vec<u32>> {
    let nside = 1u32 << order;
    let npix = 12u64 * (nside as u64) * (nside as u64);
    let mut counts = vec![0u32; npix as usize];
    let data = &mmap[MERGED_HEADER_SIZE..];
    for i in 0..total as usize {
        let offset = i * RECORD_SIZE;
        let (ra_deg, dec_deg) = extract_ra_dec_from_slice(&data[offset..offset + RECORD_SIZE]);
        let pixel = ang2pix_nest(order, ra_deg, dec_deg);
        counts[pixel as usize] += 1;
    }
    Ok(counts)
}

fn extract_ra_dec_from_slice(buf: &[u8]) -> (f64, f64) {
    let ra = f64::from_le_bytes(buf[8..16].try_into().unwrap());
    let dec = f64::from_le_bytes(buf[16..24].try_into().unwrap());
    (ra, dec)
}

fn scatter_stars_to_output(
    input_mmap: &Mmap,
    output_path: &Path,
    total_stars: u64,
    order: u32,
    mag_limit: f32,
    counts: &[u32],
) -> anyhow::Result<()> {
    let nside = 1u32 << order;
    let npix = counts.len() as u64;
    let star_data_offset = OUTPUT_HEADER_SIZE + (npix as usize) * PIXEL_ENTRY_SIZE;
    let total_file_size = star_data_offset + (total_stars as usize) * RECORD_SIZE;

    let file = create_output_file(output_path, total_file_size)?;
    let mut output_mmap = unsafe { MmapMut::map_mut(&file)? };

    write_header_to_mmap(&mut output_mmap, order, nside, npix, total_stars, mag_limit);
    let byte_offsets = write_offset_table_to_mmap(&mut output_mmap, counts);

    let mut cursors: Vec<u64> = byte_offsets.clone();
    let input_data = &input_mmap[MERGED_HEADER_SIZE..];

    for i in 0..total_stars as usize {
        let src_offset = i * RECORD_SIZE;
        let record_bytes = &input_data[src_offset..src_offset + RECORD_SIZE];
        let (ra_deg, dec_deg) = extract_ra_dec_from_slice(record_bytes);
        let pixel = ang2pix_nest(order, ra_deg, dec_deg) as usize;
        let dst_offset = star_data_offset + cursors[pixel] as usize;
        output_mmap[dst_offset..dst_offset + RECORD_SIZE].copy_from_slice(record_bytes);
        cursors[pixel] += RECORD_SIZE as u64;
    }

    output_mmap.flush()?;
    Ok(())
}

fn create_output_file(path: &Path, size: usize) -> anyhow::Result<File> {
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(path)?;
    file.set_len(size as u64)?;
    Ok(file)
}

fn write_header_to_mmap(
    mmap: &mut MmapMut,
    order: u32,
    nside: u32,
    npix: u64,
    total_stars: u64,
    mag_limit: f32,
) {
    let mut offset = 0;
    mmap[offset..offset + 4].copy_from_slice(CATALOG_MAGIC);
    offset += 4;
    mmap[offset..offset + 4].copy_from_slice(&CATALOG_VERSION.to_le_bytes());
    offset += 4;
    mmap[offset..offset + 4].copy_from_slice(&order.to_le_bytes());
    offset += 4;
    mmap[offset..offset + 4].copy_from_slice(&nside.to_le_bytes());
    offset += 4;
    mmap[offset..offset + 8].copy_from_slice(&npix.to_le_bytes());
    offset += 8;
    mmap[offset..offset + 8].copy_from_slice(&total_stars.to_le_bytes());
    offset += 8;
    mmap[offset..offset + 8].copy_from_slice(&EPOCH_J2016.to_le_bytes());
    offset += 8;
    mmap[offset..offset + 4].copy_from_slice(&mag_limit.to_le_bytes());
}

fn write_offset_table_to_mmap(mmap: &mut MmapMut, counts: &[u32]) -> Vec<u64> {
    let mut byte_offsets = Vec::with_capacity(counts.len());
    let mut current_byte_offset: u64 = 0;

    for (i, &count) in counts.iter().enumerate() {
        let table_offset = OUTPUT_HEADER_SIZE + i * PIXEL_ENTRY_SIZE;
        mmap[table_offset..table_offset + 8].copy_from_slice(&current_byte_offset.to_le_bytes());
        mmap[table_offset + 8..table_offset + 12].copy_from_slice(&count.to_le_bytes());
        mmap[table_offset + 12..table_offset + 16].copy_from_slice(&0u32.to_le_bytes());
        byte_offsets.push(current_byte_offset);
        current_byte_offset += (count as u64) * (RECORD_SIZE as u64);
    }

    byte_offsets
}

fn sort_pixels_in_place(path: &Path, counts: &[u32]) -> anyhow::Result<()> {
    let file = OpenOptions::new().read(true).write(true).open(path)?;
    let mut mmap = unsafe { MmapMut::map_mut(&file)? };
    let npix = counts.len();
    let star_data_offset = OUTPUT_HEADER_SIZE + npix * PIXEL_ENTRY_SIZE;

    let mut current_offset = star_data_offset;
    for &count in counts {
        if count > 1 {
            sort_pixel_region(&mut mmap, current_offset, count as usize);
        }
        current_offset += (count as usize) * RECORD_SIZE;
    }

    mmap.flush()?;
    Ok(())
}

fn sort_pixel_region(mmap: &mut MmapMut, offset: usize, count: usize) {
    let region = &mut mmap[offset..offset + count * RECORD_SIZE];
    let records: &mut [[u8; RECORD_SIZE]] = unsafe {
        std::slice::from_raw_parts_mut(region.as_mut_ptr() as *mut [u8; RECORD_SIZE], count)
    };
    records.sort_by(|a, b| {
        let mag_a = f32::from_le_bytes(a[48..52].try_into().unwrap());
        let mag_b = f32::from_le_bytes(b[48..52].try_into().unwrap());
        mag_a
            .partial_cmp(&mag_b)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
}

fn compact_with_cap(
    sorted_path: &Path,
    output_path: &Path,
    counts: &[u32],
    cap: u32,
    mag_limit: f32,
    order: u32,
) -> anyhow::Result<()> {
    let capped: Vec<u32> = counts.iter().map(|&c| c.min(cap)).collect();
    let capped_total: u64 = capped.iter().map(|&c| c as u64).sum();
    let npix = counts.len() as u64;
    let nside = 1u32 << order;
    let star_data_offset = OUTPUT_HEADER_SIZE + (npix as usize) * PIXEL_ENTRY_SIZE;
    let out_size = star_data_offset + (capped_total as usize) * RECORD_SIZE;

    let src_file = File::open(sorted_path)?;
    let src_mmap = unsafe { Mmap::map(&src_file)? };
    let dst_file = create_output_file(output_path, out_size)?;
    let mut dst_mmap = unsafe { MmapMut::map_mut(&dst_file)? };

    write_header_to_mmap(&mut dst_mmap, order, nside, npix, capped_total, mag_limit);
    write_offset_table_to_mmap(&mut dst_mmap, &capped);
    copy_capped_records(&src_mmap, &mut dst_mmap, counts, &capped, star_data_offset);
    dst_mmap.flush()?;
    Ok(())
}

fn copy_capped_records(
    src: &Mmap,
    dst: &mut MmapMut,
    counts: &[u32],
    capped: &[u32],
    star_data_offset: usize,
) {
    let src_star_offset = OUTPUT_HEADER_SIZE + counts.len() * PIXEL_ENTRY_SIZE;
    let mut src_pos = src_star_offset;
    let mut dst_pos = star_data_offset;

    for (i, &orig) in counts.iter().enumerate() {
        let keep = capped[i] as usize;
        let copy_bytes = keep * RECORD_SIZE;
        dst[dst_pos..dst_pos + copy_bytes].copy_from_slice(&src[src_pos..src_pos + copy_bytes]);
        src_pos += (orig as usize) * RECORD_SIZE;
        dst_pos += copy_bytes;
    }
}

fn compute_stats(
    order: u32,
    nside: u32,
    npix: u64,
    total_stars_before: u64,
    max_per_cell: Option<u32>,
    final_counts: &[u32],
    output: &Path,
    elapsed: f64,
) -> anyhow::Result<IndexStats> {
    let total_after: u64 = final_counts.iter().map(|&c| c as u64).sum();
    let non_empty: Vec<u32> = final_counts.iter().copied().filter(|&c| c > 0).collect();
    let empty_pixels = npix - non_empty.len() as u64;
    let (min_stars, max_stars) = if non_empty.is_empty() {
        (0, 0)
    } else {
        (
            *non_empty.iter().min().unwrap(),
            *non_empty.iter().max().unwrap(),
        )
    };
    let mean_stars = if non_empty.is_empty() {
        0.0
    } else {
        non_empty.iter().map(|&c| c as f64).sum::<f64>() / non_empty.len() as f64
    };
    let median_stars = compute_median(&non_empty);
    let file_size = fs::metadata(output)?.len();

    let (stars_after_cap, cells_capped) = match max_per_cell {
        Some(cap) => {
            let capped = final_counts.iter().filter(|&&c| c >= cap).count() as u64;
            (Some(total_after), Some(capped))
        }
        None => (None, None),
    };

    Ok(IndexStats {
        healpix_order: order,
        nside,
        npix,
        total_stars: total_stars_before,
        stars_after_cap,
        cells_capped,
        min_stars,
        max_stars,
        mean_stars,
        median_stars,
        empty_pixels,
        file_size,
        elapsed_secs: elapsed,
    })
}

fn compute_median(values: &[u32]) -> u32 {
    if values.is_empty() {
        return 0;
    }
    let mut sorted = values.to_vec();
    sorted.sort();
    let mid = sorted.len() / 2;
    if sorted.len().is_multiple_of(2) {
        (sorted[mid - 1] + sorted[mid]) / 2
    } else {
        sorted[mid]
    }
}

fn print_stats(stats: &IndexStats) {
    println!();
    println!("=== Index Statistics ===");
    println!("HEALPix order: {}", stats.healpix_order);
    println!("nside: {}", stats.nside);
    println!("npix: {}", stats.npix);
    println!("Total stars before cap: {}", stats.total_stars);
    if let (Some(after), Some(capped)) = (stats.stars_after_cap, stats.cells_capped) {
        println!("Total stars after cap: {}", after);
        println!(
            "Stars removed by cap: {}",
            stats.total_stars.saturating_sub(after)
        );
        let non_empty = stats.npix - stats.empty_pixels;
        println!(
            "Cells at cap: {} / {} non-empty ({:.1}%)",
            capped,
            non_empty,
            if non_empty > 0 {
                capped as f64 / non_empty as f64 * 100.0
            } else {
                0.0
            }
        );
    }
    println!("Empty pixels: {}", stats.empty_pixels);
    println!(
        "Stars per pixel (non-empty): min={}, max={}, mean={:.1}, median={}",
        stats.min_stars, stats.max_stars, stats.mean_stars, stats.median_stars
    );
    println!(
        "Output file size: {} bytes ({:.2} GB)",
        stats.file_size,
        stats.file_size as f64 / 1_073_741_824.0
    );
    println!("Elapsed time: {:.2}s", stats.elapsed_secs);
}

fn validate_output(path: &Path, order: u32) -> anyhow::Result<()> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let header = read_output_header(&mmap)?;
    validate_header_fields(&header, order)?;
    let offsets = read_offset_table(&mmap, header.npix)?;
    validate_sample_pixels(&mmap, order, &header, &offsets)?;
    println!("Validation passed.");
    Ok(())
}

struct OutputHeader {
    npix: u64,
    _total_stars: u64,
}

fn read_output_header(mmap: &Mmap) -> anyhow::Result<OutputHeader> {
    let header = &mmap[0..OUTPUT_HEADER_SIZE];
    let magic = &header[0..4];
    if magic != CATALOG_MAGIC {
        anyhow::bail!("Invalid catalog magic");
    }
    let npix = u64::from_le_bytes(header[16..24].try_into().unwrap());
    let total_stars = u64::from_le_bytes(header[24..32].try_into().unwrap());
    Ok(OutputHeader {
        npix,
        _total_stars: total_stars,
    })
}

fn validate_header_fields(header: &OutputHeader, order: u32) -> anyhow::Result<()> {
    let expected_npix = 12u64 * (1u64 << order) * (1u64 << order);
    if header.npix != expected_npix {
        anyhow::bail!(
            "npix mismatch: expected {}, got {}",
            expected_npix,
            header.npix
        );
    }
    Ok(())
}

struct PixelEntry {
    offset: u64,
    count: u32,
}

fn read_offset_table(mmap: &Mmap, npix: u64) -> anyhow::Result<Vec<PixelEntry>> {
    let mut entries = Vec::with_capacity(npix as usize);
    for i in 0..npix as usize {
        let table_offset = OUTPUT_HEADER_SIZE + i * PIXEL_ENTRY_SIZE;
        let offset = u64::from_le_bytes(mmap[table_offset..table_offset + 8].try_into().unwrap());
        let count = u32::from_le_bytes(
            mmap[table_offset + 8..table_offset + 12]
                .try_into()
                .unwrap(),
        );
        entries.push(PixelEntry { offset, count });
    }
    Ok(entries)
}

fn validate_sample_pixels(
    mmap: &Mmap,
    order: u32,
    header: &OutputHeader,
    offsets: &[PixelEntry],
) -> anyhow::Result<()> {
    let star_data_offset = OUTPUT_HEADER_SIZE + (header.npix as usize) * PIXEL_ENTRY_SIZE;
    let samples = pick_sample_pixels(offsets);
    for (pixel_idx, entry) in samples {
        validate_single_pixel(mmap, order, pixel_idx, entry, star_data_offset)?;
    }
    Ok(())
}

fn pick_sample_pixels(offsets: &[PixelEntry]) -> Vec<(usize, &PixelEntry)> {
    let non_empty: Vec<(usize, &PixelEntry)> = offsets
        .iter()
        .enumerate()
        .filter(|(_, e)| e.count > 0)
        .collect();
    if non_empty.is_empty() {
        return vec![];
    }
    let mut samples = Vec::new();
    if !non_empty.is_empty() {
        samples.push(non_empty[0]);
    }
    if non_empty.len() > 1 {
        samples.push(non_empty[non_empty.len() / 2]);
    }
    if non_empty.len() > 2 {
        samples.push(non_empty[non_empty.len() - 1]);
    }
    samples
}

fn validate_single_pixel(
    mmap: &Mmap,
    order: u32,
    pixel_idx: usize,
    entry: &PixelEntry,
    star_data_offset: usize,
) -> anyhow::Result<()> {
    let pixel_start = star_data_offset + entry.offset as usize;
    let mut prev_mag: Option<f32> = None;
    for i in 0..entry.count {
        let record_offset = pixel_start + (i as usize) * RECORD_SIZE;
        let record_bytes = &mmap[record_offset..record_offset + RECORD_SIZE];
        let ra = f64::from_le_bytes(record_bytes[8..16].try_into().unwrap());
        let dec = f64::from_le_bytes(record_bytes[16..24].try_into().unwrap());
        let mag = f32::from_le_bytes(record_bytes[48..52].try_into().unwrap());
        let computed_pixel = ang2pix_nest(order, ra, dec);
        if computed_pixel != pixel_idx as u64 {
            anyhow::bail!(
                "Pixel mismatch: star at ({:.6}, {:.6}) expected pixel {}, got {}",
                ra,
                dec,
                pixel_idx,
                computed_pixel
            );
        }
        if let Some(prev) = prev_mag {
            if mag < prev {
                anyhow::bail!(
                    "Magnitude not sorted in pixel {}: star {} has mag {:.2}, prev was {:.2}",
                    pixel_idx,
                    i,
                    mag,
                    prev
                );
            }
        }
        prev_mag = Some(mag);
    }
    println!("  Pixel {}: {} stars verified", pixel_idx, entry.count);
    Ok(())
}
