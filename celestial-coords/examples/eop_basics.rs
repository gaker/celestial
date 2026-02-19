use celestial_coords::eop::record::EopRecord;
use celestial_coords::EopProvider;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // --- Bundled IERS data ---
    // The library ships with real IERS C04 + finals2000A data (updated weekly).
    // C04 covers 1962-present observed values; finals extends with ~1yr predictions.

    let provider = EopProvider::bundled()?;

    let mjd = 59945.0; // 2023-01-01
    let params = provider.get(mjd)?;

    println!("Bundled IERS data lookup for MJD {mjd}:");
    println!("  {params}");
    println!("  LOD      = {:.7} s", params.lod);
    println!("  s'       = {:.4e} rad", params.s_prime);
    println!("  source   = {:?}", params.flags.source);
    println!("  quality  = {:?}", params.flags.quality);

    if let Some((start, end)) = provider.time_span() {
        println!(
            "  coverage = MJD {start:.0} to {end:.0} ({:.0} days)",
            end - start
        );
    }
    println!("  records  = {}", provider.record_count());
    println!();

    // --- Manual EOP records ---
    // Build records by hand for testing or when you have your own data source.

    let r1 = EopRecord::new(60000.0, 0.100, 0.250, -0.050, 0.0015)?
        .with_cip_offsets(0.120, -0.080)?
        .with_pole_rates(0.00012, -0.00008)?;

    let r2 = EopRecord::new(60001.0, 0.102, 0.248, -0.052, 0.0016)?
        .with_cip_offsets(0.125, -0.082)?
        .with_pole_rates(0.00013, -0.00009)?;

    let r3 = EopRecord::new(60002.0, 0.104, 0.246, -0.054, 0.0014)?
        .with_cip_offsets(0.130, -0.084)?
        .with_pole_rates(0.00014, -0.00010)?;

    let provider = EopProvider::from_records(vec![r1, r2, r3])?;

    let interp = provider.get(60000.5)?;
    println!("Interpolated at MJD 60000.5:");
    println!("  {interp}");
    println!("  dX       = {:.3} mas", interp.dx.unwrap_or(0.0));
    println!("  dY       = {:.3} mas", interp.dy.unwrap_or(0.0));
    println!("  xrt      = {:.6} arcsec/day", interp.xrt.unwrap_or(0.0));
    println!("  yrt      = {:.6} arcsec/day", interp.yrt.unwrap_or(0.0));
    println!("  has_cip   = {}", interp.flags.has_cip_offsets);
    println!("  has_rates = {}", interp.flags.has_pole_rates);
    println!();

    // --- Data span ---
    if let Some((start, end)) = provider.time_span() {
        println!(
            "Loaded data covers MJD {start:.0} to {end:.0} ({:.0} days)",
            end - start
        );
    }

    Ok(())
}
