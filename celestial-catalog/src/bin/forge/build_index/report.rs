// Opens the finished catalog with the library reader, prints its statistics and checks a few
// pixels' records against their pixel and magnitude order.

use anyhow::Context;
use celestial_catalog::query::catalog::Catalog;
use celestial_catalog::query::healpix::ang2pix_nest;
use std::cmp::Ordering;
use std::path::Path;
use std::time::Instant;

// Stars per non-empty pixel.
#[derive(Debug, PartialEq)]
struct Distribution {
    non_empty: usize,
    min: usize,
    max: usize,
    mean: f64,
    median: usize,
}

pub(super) fn check_and_report(
    path: &Path,
    merged_stars: u64,
    max_per_cell: Option<u32>,
    start: Instant,
) -> anyhow::Result<()> {
    let catalog =
        Catalog::open(path).with_context(|| format!("Built catalog {path:?} is invalid"))?;
    let counts: Vec<usize> = (0..catalog.header().npix)
        .map(|pixel| catalog.stars_in_pixel(pixel).len())
        .collect();
    print_stats(&catalog, &counts, merged_stars, max_per_cell);
    println!("Elapsed time: {:.2}s", start.elapsed().as_secs_f64());
    println!("\nValidating output...");
    for pixel in sample_pixels(&counts) {
        check_pixel(&catalog, pixel)?;
    }
    println!("Validation passed.");
    Ok(())
}

fn distribution(counts: &[usize]) -> Distribution {
    let mut non_empty: Vec<usize> = counts.iter().copied().filter(|&c| c > 0).collect();
    non_empty.sort_unstable();
    let n = non_empty.len();
    let at = |i: usize| non_empty.get(i).copied().unwrap_or(0);
    let median = if n.is_multiple_of(2) && n > 0 {
        (at(n / 2 - 1) + at(n / 2)) / 2
    } else {
        at(n / 2)
    };
    let sum: usize = non_empty.iter().sum();
    Distribution {
        non_empty: n,
        min: at(0),
        max: at(n.saturating_sub(1)),
        mean: if n == 0 { 0.0 } else { sum as f64 / n as f64 },
        median,
    }
}

fn print_stats(catalog: &Catalog, counts: &[usize], merged_stars: u64, cap: Option<u32>) {
    let header = catalog.header();
    let spread = distribution(counts);
    println!();
    println!("=== Index Statistics ===");
    println!("HEALPix order: {}", header.order);
    println!("nside: {}", header.nside);
    println!("npix: {}", header.npix);
    println!("Epoch: J{}", header.epoch);
    println!("Magnitude limit: {:.2}", header.mag_limit);
    println!("Total stars before cap: {}", merged_stars);
    if let Some(cap) = cap {
        print_cap_stats(counts, &spread, merged_stars, header.total_stars, cap);
    }
    println!("Empty pixels: {}", counts.len() - spread.non_empty);
    println!(
        "Stars per pixel (non-empty): min={}, max={}, mean={:.1}, median={}",
        spread.min, spread.max, spread.mean, spread.median
    );
    let size = catalog.file_size();
    println!(
        "Output file size: {} bytes ({:.2} GB)",
        size,
        size as f64 / 1_073_741_824.0
    );
}

fn print_cap_stats(counts: &[usize], spread: &Distribution, before: u64, after: u64, cap: u32) {
    let at_cap = counts.iter().filter(|&&c| c >= cap as usize).count();
    let share = if spread.non_empty > 0 {
        at_cap as f64 / spread.non_empty as f64 * 100.0
    } else {
        0.0
    };
    println!("Total stars after cap: {}", after);
    println!("Stars removed by cap: {}", before.saturating_sub(after));
    println!(
        "Cells at cap: {} / {} non-empty ({:.1}%)",
        at_cap, spread.non_empty, share
    );
}

// The first, middle and last non-empty pixels.
fn sample_pixels(counts: &[usize]) -> Vec<u64> {
    let non_empty: Vec<u64> = (0u64..)
        .zip(counts)
        .filter(|(_, &count)| count > 0)
        .map(|(pixel, _)| pixel)
        .collect();
    let mut samples: Vec<u64> = [0, non_empty.len() / 2, non_empty.len().saturating_sub(1)]
        .iter()
        .filter_map(|&i| non_empty.get(i).copied())
        .collect();
    samples.dedup();
    samples
}

fn check_pixel(catalog: &Catalog, pixel: u64) -> anyhow::Result<()> {
    let order = catalog.header().order;
    let stars = catalog.stars_in_pixel(pixel);
    for (i, star) in stars.iter().enumerate() {
        let found = ang2pix_nest(order, star.ra, star.dec)?;
        if found != pixel {
            anyhow::bail!(
                "Pixel mismatch: star at ({:.6}, {:.6}) is in pixel {pixel}, but belongs in {found}",
                star.ra,
                star.dec
            );
        }
        if i > 0 && stars[i - 1].mag.total_cmp(&star.mag) == Ordering::Greater {
            anyhow::bail!(
                "Magnitude not sorted in pixel {pixel}: star {i} has mag {:.2}, prev was {:.2}",
                star.mag,
                stars[i - 1].mag
            );
        }
    }
    println!("  Pixel {}: {} stars verified", pixel, stars.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distribution_of_counts() {
        let spread = distribution(&[0, 4, 0, 1, 2, 7]);
        let expected = Distribution {
            non_empty: 4,
            min: 1,
            max: 7,
            mean: 3.5,
            median: 3,
        };
        assert_eq!(spread, expected);
    }

    #[test]
    fn test_distribution_of_empty_pixels() {
        let spread = distribution(&[0, 0]);
        let expected = Distribution {
            non_empty: 0,
            min: 0,
            max: 0,
            mean: 0.0,
            median: 0,
        };
        assert_eq!(spread, expected);
    }

    #[test]
    fn test_sample_pixels() {
        assert_eq!(sample_pixels(&[0, 3, 0, 1, 5, 0, 2]), vec![1, 4, 6]);
        assert_eq!(sample_pixels(&[0, 3, 0]), vec![1]);
        assert_eq!(sample_pixels(&[0, 0]), Vec::<u64>::new());
    }
}
