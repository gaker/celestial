//! Hipparcos New Reduction (van Leeuwen 2007) catalog ingestion
//!
//! Parses hip2.dat fixed-width format from I/311.
//! Auto-downloads from CDS if not present locally.
//! Epoch-propagates positions from J1991.25 to J2016.0.

use crate::cli::{Cli, IngestHipparcosArgs};
use crate::gaia::{FLAG_HAS_PARALLAX, FLAG_HAS_PROPER_MOTION, FLAG_SOURCE_HIPPARCOS};
use anyhow::Context;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

const RECORD_SIZE: usize = 56;
const FINAL_MAGIC: &[u8; 4] = b"HIPP";
const FINAL_VERSION: u32 = 2;
const HIPPARCOS_EPOCH: f64 = 1991.25;
const GAIA_EPOCH: f64 = 2016.0;
const DELTA_T_YEARS: f64 = GAIA_EPOCH - HIPPARCOS_EPOCH;
const HIP_SYNTHETIC_ID_BASE: i64 = 0x4000_0000_0000_0000;
const HIP2_URL: &str = "https://cdsarc.cds.unistra.fr/ftp/I/311/hip2.dat.gz";
const CROSSMATCH_URL: &str =
    "https://cdn.gea.esac.esa.int/Gaia/gedr3/cross_match/hipparcos2_best_neighbour/Hipparcos2BestNeighbour.csv.gz";
const HIP2_FILENAME: &str = "hip2.dat";
const CROSSMATCH_FILENAME: &str = "Hipparcos2BestNeighbour.csv";

#[repr(C)]
struct StarRecord {
    source_id: i64,
    ra: f64,
    dec: f64,
    pmra: f64,
    pmdec: f64,
    parallax: f64,
    mag: f32,
    flags: u16,
    _padding: u16,
}

#[allow(dead_code)]
struct Hip2Star {
    hip: u32,
    ra_rad: f64,
    dec_rad: f64,
    parallax: f64,
    pmra: f64,
    pmdec: f64,
    hpmag: f32,
    b_v: Option<f64>,
    v_i: Option<f64>,
    solution_type: u16,
    num_components: u8,
}

struct IngestStats {
    total_parsed: u64,
    kept_after_mag: u64,
    with_gaia_match: u64,
    without_match: u64,
}

pub fn run(args: &IngestHipparcosArgs, cli: &Cli) -> anyhow::Result<()> {
    fs::create_dir_all(&args.workdir)?;
    fs::create_dir_all(&args.output)?;
    let hip_path = ensure_file(
        &args.workdir.join(HIP2_FILENAME),
        HIP2_URL,
        "hip2.dat (van Leeuwen 2007)",
    )?;
    let crossmatch_path = ensure_file(
        &args.workdir.join(CROSSMATCH_FILENAME),
        CROSSMATCH_URL,
        "Hipparcos2BestNeighbour.csv (Gaia eDR3)",
    )?;
    print_plan(args, cli);
    let crossmatch = load_crossmatch(&crossmatch_path)?;
    println!("Loaded {} cross-match entries", crossmatch.len());
    let stars = parse_hip2(&hip_path, args.mag_limit)?;
    println!(
        "Parsed {} stars (mag <= {:.1})",
        stars.len(),
        args.mag_limit
    );
    let stats = write_output(&stars, &crossmatch, &args.output, args.mag_limit)?;
    print_validation(&stars, &crossmatch);
    print_summary(&stats);
    Ok(())
}

fn ensure_file(path: &Path, url: &str, label: &str) -> anyhow::Result<std::path::PathBuf> {
    if path.exists() {
        println!("Found {}: {:?}", label, path);
        return Ok(path.to_path_buf());
    }
    let gz_path = path.with_extension(gz_extension(path));
    if gz_path.exists() {
        println!("Decompressing {:?}...", gz_path);
        decompress_gz(&gz_path, path)?;
        return Ok(path.to_path_buf());
    }
    println!("{} not found at {:?}", label, path);
    println!("Downloading: {}", url);
    download_and_decompress(url, path)?;
    Ok(path.to_path_buf())
}

fn gz_extension(path: &Path) -> String {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    format!("{}.gz", ext)
}

fn download_and_decompress(url: &str, dest: &Path) -> anyhow::Result<()> {
    let response =
        reqwest::blocking::get(url).with_context(|| format!("Failed to download {}", url))?;
    if !response.status().is_success() {
        anyhow::bail!("Download failed with status {}", response.status());
    }
    let compressed = response.bytes().context("Failed to read response")?;
    println!("Downloaded {:.1} MB", compressed.len() as f64 / 1_048_576.0);
    let decoder = flate2::read::GzDecoder::new(&compressed[..]);
    let mut reader = BufReader::new(decoder);
    let mut file = File::create(dest)?;
    let mut buf = vec![0u8; 256 * 1024];
    let mut written = 0u64;
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        file.write_all(&buf[..n])?;
        written += n as u64;
    }
    println!(
        "Decompressed to {:?} ({:.1} MB)",
        dest,
        written as f64 / 1_048_576.0
    );
    Ok(())
}

fn decompress_gz(gz_path: &Path, dest: &Path) -> anyhow::Result<()> {
    let gz_file = File::open(gz_path)?;
    let decoder = flate2::read::GzDecoder::new(BufReader::new(gz_file));
    let mut reader = BufReader::new(decoder);
    let mut file = File::create(dest)?;
    std::io::copy(&mut reader, &mut file)?;
    Ok(())
}

fn print_plan(args: &IngestHipparcosArgs, cli: &Cli) {
    println!("\n=== Hipparcos New Reduction Ingestion ===");
    println!("Workdir: {:?}", args.workdir);
    println!("Magnitude limit: {:.1}", args.mag_limit);
    println!("Output directory: {:?}", args.output);
    println!("Verbose: {}", cli.verbose);
    println!();
}

fn load_crossmatch(path: &Path) -> anyhow::Result<HashMap<u32, i64>> {
    let file = File::open(path).context("Failed to open cross-match file")?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        if line_num == 0 && line.contains("source_id") {
            continue;
        }
        if let Some((hip, gaia_id)) = parse_crossmatch_line(&line) {
            map.insert(hip, gaia_id);
        }
    }
    Ok(map)
}

fn parse_crossmatch_line(line: &str) -> Option<(u32, i64)> {
    let fields: Vec<&str> = line.split(',').collect();
    if fields.len() < 2 {
        return None;
    }
    let gaia_id: i64 = fields[0].trim().parse().ok()?;
    let hip: u32 = fields[1].trim().parse().ok()?;
    Some((hip, gaia_id))
}

fn parse_hip2(path: &Path, mag_limit: f32) -> anyhow::Result<Vec<Hip2Star>> {
    let file = File::open(path).context("Failed to open hip2.dat")?;
    let reader = BufReader::new(file);
    let mut stars = Vec::with_capacity(120_000);
    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        match parse_hip2_line(&line) {
            Some(star) if star.hpmag <= mag_limit => stars.push(star),
            Some(_) => {}
            None => {
                if line.len() > 10 {
                    eprintln!("Warning: failed to parse line {}", line_num + 1);
                }
            }
        }
    }
    Ok(stars)
}

fn parse_hip2_line(line: &str) -> Option<Hip2Star> {
    if line.len() < 171 {
        return None;
    }
    let bytes = line.as_bytes();
    let hip: u32 = col(bytes, 0, 6)?.trim().parse().ok()?;
    let sn: u16 = col(bytes, 7, 10)?.trim().parse().ok()?;
    let nc: u8 = col(bytes, 13, 14)?.trim().parse().ok()?;
    let ra_rad: f64 = col(bytes, 15, 28)?.trim().parse().ok()?;
    let dec_rad: f64 = col(bytes, 29, 42)?.trim().parse().ok()?;
    let parallax: f64 = col(bytes, 43, 50)?.trim().parse().unwrap_or(0.0);
    let pmra: f64 = col(bytes, 51, 59)?.trim().parse().unwrap_or(0.0);
    let pmdec: f64 = col(bytes, 60, 68)?.trim().parse().unwrap_or(0.0);
    let hpmag_str = col(bytes, 129, 136)?.trim();
    if hpmag_str.is_empty() {
        return None;
    }
    let hpmag: f32 = hpmag_str.parse().ok()?;
    let b_v = col(bytes, 152, 158).and_then(|s| s.trim().parse::<f64>().ok());
    let v_i = col(bytes, 165, 171).and_then(|s| s.trim().parse::<f64>().ok());

    Some(Hip2Star {
        hip,
        ra_rad,
        dec_rad,
        parallax,
        pmra,
        pmdec,
        hpmag,
        b_v,
        v_i,
        solution_type: sn,
        num_components: nc,
    })
}

fn col(bytes: &[u8], start: usize, end: usize) -> Option<&str> {
    if end > bytes.len() {
        return None;
    }
    std::str::from_utf8(&bytes[start..end]).ok()
}

fn propagate_position(star: &Hip2Star) -> (f64, f64) {
    let pmdec_rad_per_yr = star.pmdec / 3_600_000.0 * (celestial_core::constants::PI / 180.0);
    let pmra_rad_per_yr = star.pmra / 3_600_000.0 * (celestial_core::constants::PI / 180.0);
    let dec_2016 = star.dec_rad + pmdec_rad_per_yr * DELTA_T_YEARS;
    let cos_dec = libm::cos(star.dec_rad);
    let ra_2016 = if libm::fabs(cos_dec) > 1e-10 {
        star.ra_rad + pmra_rad_per_yr * DELTA_T_YEARS / cos_dec
    } else {
        star.ra_rad
    };
    (ra_2016, dec_2016)
}

fn compute_source_id(hip: u32, crossmatch: &HashMap<u32, i64>) -> i64 {
    crossmatch
        .get(&hip)
        .copied()
        .unwrap_or(HIP_SYNTHETIC_ID_BASE | hip as i64)
}

fn compute_flags(star: &Hip2Star) -> u16 {
    let mut flags = FLAG_SOURCE_HIPPARCOS;
    if star.pmra != 0.0 || star.pmdec != 0.0 {
        flags |= FLAG_HAS_PROPER_MOTION;
    }
    if star.parallax != 0.0 {
        flags |= FLAG_HAS_PARALLAX;
    }
    flags
}

fn write_output(
    stars: &[Hip2Star],
    crossmatch: &HashMap<u32, i64>,
    output_dir: &Path,
    mag_limit: f32,
) -> anyhow::Result<IngestStats> {
    let output_path = output_dir.join("hipparcos_ingest.bin");
    let file = File::create(&output_path)?;
    let mut writer = BufWriter::new(file);
    write_header(&mut writer, stars.len() as u64, mag_limit)?;
    let stats = write_records(&mut writer, stars, crossmatch)?;
    writer.flush()?;
    println!("Written {} stars to {:?}", stars.len(), output_path);
    Ok(stats)
}

fn write_header<W: Write>(w: &mut W, count: u64, mag_limit: f32) -> anyhow::Result<()> {
    w.write_all(FINAL_MAGIC)?;
    w.write_all(&FINAL_VERSION.to_le_bytes())?;
    w.write_all(&count.to_le_bytes())?;
    w.write_all(&mag_limit.to_le_bytes())?;
    w.write_all(&[0u8; 4])?;
    Ok(())
}

fn write_records<W: Write>(
    writer: &mut W,
    stars: &[Hip2Star],
    crossmatch: &HashMap<u32, i64>,
) -> anyhow::Result<IngestStats> {
    let mut stats = IngestStats {
        total_parsed: stars.len() as u64,
        kept_after_mag: stars.len() as u64,
        with_gaia_match: 0,
        without_match: 0,
    };
    for star in stars {
        let (ra_2016, dec_2016) = propagate_position(star);
        let source_id = compute_source_id(star.hip, crossmatch);
        if crossmatch.contains_key(&star.hip) {
            stats.with_gaia_match += 1;
        } else {
            stats.without_match += 1;
        }
        write_star_record(writer, star, ra_2016, dec_2016, source_id)?;
    }
    Ok(stats)
}

fn write_star_record<W: Write>(
    writer: &mut W,
    star: &Hip2Star,
    ra_rad: f64,
    dec_rad: f64,
    source_id: i64,
) -> anyhow::Result<()> {
    let record = StarRecord {
        source_id,
        ra: ra_rad * 180.0 / celestial_core::constants::PI,
        dec: dec_rad * 180.0 / celestial_core::constants::PI,
        pmra: star.pmra,
        pmdec: star.pmdec,
        parallax: star.parallax,
        mag: star.hpmag,
        flags: compute_flags(star),
        _padding: 0,
    };
    let bytes = unsafe {
        std::slice::from_raw_parts(&record as *const StarRecord as *const u8, RECORD_SIZE)
    };
    writer.write_all(bytes)?;
    Ok(())
}

fn print_validation(stars: &[Hip2Star], crossmatch: &HashMap<u32, i64>) {
    println!("\n=== Validation (known stars) ===");
    validate_star(stars, crossmatch, 32349, "Sirius");
    validate_star(stars, crossmatch, 91262, "Vega");
}

fn validate_star(stars: &[Hip2Star], crossmatch: &HashMap<u32, i64>, hip: u32, name: &str) {
    let Some(star) = stars.iter().find(|s| s.hip == hip) else {
        println!("HIP {} ({}): not found in filtered catalog", hip, name);
        return;
    };
    let (ra_2016, dec_2016) = propagate_position(star);
    let source_id = compute_source_id(hip, crossmatch);
    let ra_deg = ra_2016 * 180.0 / celestial_core::constants::PI;
    let dec_deg = dec_2016 * 180.0 / celestial_core::constants::PI;
    let match_status = if crossmatch.contains_key(&hip) {
        "Gaia match"
    } else {
        "synthetic ID"
    };
    println!(
        "HIP {} ({}): RA={:.6} Dec={:.6} deg (J2016.0), Hp={:.3}, {}",
        hip, name, ra_deg, dec_deg, star.hpmag, match_status
    );
    let orig_ra_deg = star.ra_rad * 180.0 / celestial_core::constants::PI;
    let orig_dec_deg = star.dec_rad * 180.0 / celestial_core::constants::PI;
    println!(
        "  Original (J1991.25): RA={:.6} Dec={:.6} deg",
        orig_ra_deg, orig_dec_deg
    );
    println!(
        "  PM: pmRA={:.2} mas/yr, pmDec={:.2} mas/yr",
        star.pmra, star.pmdec
    );
    if let Some(bv) = star.b_v {
        println!("  B-V={:.3}", bv);
    }
    println!("  Source ID: {}", source_id);
}

fn print_summary(stats: &IngestStats) {
    println!("\n=== Summary ===");
    println!("Total parsed: {}", stats.total_parsed);
    println!("Kept after mag filter: {}", stats.kept_after_mag);
    println!("With Gaia match: {}", stats.with_gaia_match);
    println!("Without match (synthetic ID): {}", stats.without_match);
}
