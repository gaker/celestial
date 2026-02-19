use celestial_coords::EopProvider;
use std::path::PathBuf;

const FINALS_URL: &str = "https://datacenter.iers.org/data/9/finals2000A.all";

#[tokio::main(flavor = "current_thread")]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cache_dir = cache_dir();
    let finals_path = cache_dir.join("finals2000A.all");

    // Show bundled data coverage
    let bundled = EopProvider::bundled()?;
    let (bstart, bend) = bundled.time_span().unwrap();
    println!(
        "Bundled EOP data: MJD {bstart:.0} to {bend:.0} ({} records)",
        bundled.record_count()
    );

    // Download fresh finals2000A if not cached (or always, in production)
    if !finals_path.exists() {
        println!("\nDownloading finals2000A from IERS...");
        let body = reqwest::get(FINALS_URL).await?.text().await?;
        std::fs::create_dir_all(&cache_dir)?;
        std::fs::write(&finals_path, &body)?;
        println!("Saved to {}", finals_path.display());
    } else {
        println!("\nUsing cached {}", finals_path.display());
    }

    // Load from file only
    let from_file = EopProvider::from_finals_file(&finals_path)?;
    let (fstart, fend) = from_file.time_span().unwrap();
    println!(
        "\nFinals-only: MJD {fstart:.0} to {fend:.0} ({} records)",
        from_file.record_count()
    );

    // Bundled + update overlay (the telescope control pattern)
    let merged = EopProvider::bundled_with_update(&finals_path)?;
    let (mstart, mend) = merged.time_span().unwrap();
    println!(
        "Merged:      MJD {mstart:.0} to {mend:.0} ({} records)",
        merged.record_count()
    );

    // Compare a lookup across providers
    let test_mjd = bend - 10.0; // 10 days before end of bundled data
    let bundled_params = bundled.get(test_mjd)?;
    let merged_params = merged.get(test_mjd)?;
    println!("\nLookup at MJD {test_mjd:.1}:");
    println!("  Bundled: {bundled_params}");
    println!("  Merged:  {merged_params}");

    // Try a date beyond bundled range (if the finals data extends further)
    if mend > bend {
        let future_mjd = bend + 30.0;
        match merged.get(future_mjd) {
            Ok(params) => {
                println!("\nFuture lookup at MJD {future_mjd:.1} (beyond bundled):");
                println!("  {params}");
            }
            Err(e) => println!("\nMJD {future_mjd:.1} not available: {e}"),
        }
    }

    Ok(())
}

fn cache_dir() -> PathBuf {
    // XDG on Linux, ~/Library/Caches on macOS, AppData\Local on Windows
    #[cfg(target_os = "macos")]
    let base = std::env::var("HOME").map(|h| PathBuf::from(h).join("Library/Caches"));
    #[cfg(target_os = "linux")]
    let base = std::env::var("XDG_CACHE_HOME")
        .or_else(|_| std::env::var("HOME").map(|h| format!("{h}/.cache")))
        .map(PathBuf::from);
    #[cfg(target_os = "windows")]
    let base = std::env::var("LOCALAPPDATA").map(PathBuf::from);

    base.unwrap_or_else(|_| PathBuf::from("."))
        .join("celestial")
}
