use std::path::{Path, PathBuf};

use crate::download::{default_client, download_all, find_elp_files};
use crate::generate::{format_analysis, generate_moon_module, GenerateConfig};
use crate::parser::{parse_files, ElpData};

fn read_elp(input: &Path) -> Result<ElpData, String> {
    let paths = find_elp_files(input)
        .ok_or_else(|| format!("ELP/MPP02 files not found in {}", input.display()))?;
    println!("Parsing ELP/MPP02 files from {}...", input.display());
    parse_files(&paths).map_err(|e| format!("Parse error: {}", e))
}

pub(crate) fn cmd_download(output: PathBuf) -> Result<(), String> {
    std::fs::create_dir_all(&output).map_err(|e| format!("Failed to create output dir: {}", e))?;
    let client = default_client()?;
    download_all(&client, &output)
}

pub(crate) fn cmd_analyze(input: PathBuf, threshold: f64) -> Result<(), String> {
    let elp = read_elp(&input)?;
    print!("{}", format_analysis(&elp, threshold));
    Ok(())
}

pub(crate) fn cmd_generate(input: PathBuf, output: PathBuf, threshold: f64) -> Result<(), String> {
    let elp = read_elp(&input)?;
    let config = GenerateConfig {
        threshold,
        output_dir: output,
    };
    println!("\nGenerating Rust code...");
    println!("{}", generate_moon_module(&elp, &config)?);
    println!("\nGeneration complete!");
    println!("Output directory: {}", config.output_dir.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::download::{MAIN_FILES, PERT_FILES};
    use tempfile::TempDir;

    const MAIN_LINE: &str = "  0  2  0  0     -411.60287      168.48   -18433.81     -121.62        0.40       -0.18        0.00";
    const PERT_LINE: &str = "    1-0.1274921554086D+02 0.6368794709728D+01  0  0  1  0  0-18 16  0  0  0  0  0  0  0  0  0";

    fn elp_dir(main: &str, pert: &str) -> TempDir {
        let dir = TempDir::new().unwrap();
        for name in MAIN_FILES {
            std::fs::write(dir.path().join(name), main).unwrap();
        }
        for name in PERT_FILES {
            std::fs::write(dir.path().join(name), pert).unwrap();
        }
        dir
    }

    #[test]
    fn missing_or_bad_inputs() {
        let empty = TempDir::new().unwrap();
        let not_found = format!("ELP/MPP02 files not found in {}", empty.path().display());
        let damaged = elp_dir("", "");
        let cases = [
            (&empty, not_found),
            (
                &damaged,
                "Parse error: Invalid format: Empty file".to_string(),
            ),
        ];
        for (dir, message) in cases {
            let input = dir.path().to_path_buf();
            let output = dir.path().join("output");
            assert_eq!(cmd_analyze(input.clone(), 1e-3), Err(message.clone()));
            assert_eq!(cmd_generate(input, output.clone(), 1e-3), Err(message));
            assert!(!output.exists());
        }
    }

    #[test]
    fn analyzes_and_generates() {
        let main = format!("MAIN PROBLEM. LONGITUDE. 1\n{MAIN_LINE}\n");
        let pert = format!("PERTURBATIONS. LONGITUDE. 1 0\n{PERT_LINE}\n");
        let dir = elp_dir(&main, &pert);
        let input = dir.path().to_path_buf();
        assert_eq!(cmd_analyze(input.clone(), 1e-3), Ok(()));
        let output = dir.path().join("output");
        assert_eq!(cmd_generate(input, output.clone(), 1e-3), Ok(()));
        let source = std::fs::read_to_string(output.join("moon.rs")).unwrap();
        assert!(source.contains("--threshold 1e-3`"));
        assert!(source.contains("Terms retained: 6 of 6 (100.0%)"));
    }
}
