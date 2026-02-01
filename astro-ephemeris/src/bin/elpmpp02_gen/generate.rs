use std::fs::{self, File};
use std::io::Write;

#[cfg(feature = "cli")]
use chrono::Utc;

use super::parser::{Coordinate, ElpData, MainSeries, MainTerm, PertBlock, PertTerm};

pub struct GenerateConfig {
    pub threshold: f64,
    pub output_dir: std::path::PathBuf,
}

fn filter_main_terms(series: &MainSeries, threshold: f64) -> Vec<&MainTerm> {
    let mut filtered: Vec<_> = series
        .terms
        .iter()
        .filter(|t| t.amplitude() >= threshold)
        .collect();
    filtered.sort_by(|a, b| b.amplitude().partial_cmp(&a.amplitude()).unwrap());
    filtered
}

fn filter_pert_terms(block: &PertBlock, threshold: f64) -> Vec<&PertTerm> {
    let mut filtered: Vec<_> = block
        .terms
        .iter()
        .filter(|t| t.amplitude() >= threshold)
        .collect();
    filtered.sort_by(|a, b| b.amplitude().partial_cmp(&a.amplitude()).unwrap());
    filtered
}

fn coord_name(coord: Coordinate) -> &'static str {
    match coord {
        Coordinate::Longitude => "LONGITUDE",
        Coordinate::Latitude => "LATITUDE",
        Coordinate::Distance => "DISTANCE",
    }
}

#[allow(dead_code)]
fn coord_var(coord: Coordinate) -> &'static str {
    match coord {
        Coordinate::Longitude => "longitude",
        Coordinate::Latitude => "latitude",
        Coordinate::Distance => "distance",
    }
}

#[cfg(feature = "cli")]
pub fn generate_moon_module(
    elp: &ElpData,
    config: &GenerateConfig,
) -> Result<(usize, usize), String> {
    let threshold = config.threshold;
    let date = Utc::now().format("%Y-%m-%d").to_string();

    let mut out = String::new();

    out.push_str("#![allow(clippy::excessive_precision)]\n");
    out.push_str("#![allow(clippy::unreadable_literal)]\n");
    out.push_str("#![allow(clippy::approx_constant)]\n");
    out.push('\n');
    out.push_str("//! ELP/MPP02 coefficients for the Moon\n");
    out.push_str("//!\n");
    out.push_str("//! Generated from official ELP/MPP02 data files\n");
    out.push_str(&format!("//! Threshold: {:.0e}\n", threshold));
    out.push_str(&format!("//! Generated: {}\n", date));
    out.push_str("//!\n");
    out.push_str("//! Reference: Chapront & Francou (2003)\n");
    out.push_str(
        "//! \"The lunar theory ELP revisited. Introduction of new planetary perturbations\"\n",
    );
    out.push_str("//! Astronomy & Astrophysics, 404, 735-742\n");
    out.push('\n');

    out.push_str("/// Main problem term with Delaunay argument multipliers\n");
    out.push_str("#[derive(Debug, Clone, Copy)]\n");
    out.push_str("pub struct MainTerm {\n");
    out.push_str("    /// Multipliers for D, F, l, l' (Delaunay arguments)\n");
    out.push_str("    pub delaunay: [i8; 4],\n");
    out.push_str("    /// Polynomial coefficients A0 through A6\n");
    out.push_str("    pub coeffs: [f64; 7],\n");
    out.push_str("}\n\n");

    out.push_str("/// Perturbation term with full argument multipliers\n");
    out.push_str("#[derive(Debug, Clone, Copy)]\n");
    out.push_str("pub struct PertTerm {\n");
    out.push_str("    /// Amplitude (sqrt(S^2 + C^2))\n");
    out.push_str("    pub amplitude: f64,\n");
    out.push_str("    /// Phase angle (atan2(C, S))\n");
    out.push_str("    pub phase: f64,\n");
    out.push_str(
        "    /// Multipliers: [D, F, l, l', Me, Ve, Te, Ma, Ju, Sa, Ur, Ne, zeta, ?, ?, ?]\n",
    );
    out.push_str("    pub multipliers: [i8; 16],\n");
    out.push_str("}\n\n");

    out.push_str("/// Perturbation block for a specific time power\n");
    out.push_str("#[derive(Debug, Clone, Copy)]\n");
    out.push_str("pub struct PertBlock {\n");
    out.push_str("    /// Power of T (0, 1, 2, 3)\n");
    out.push_str("    pub power: u8,\n");
    out.push_str("    /// Terms for this power\n");
    out.push_str("    pub terms: &'static [PertTerm],\n");
    out.push_str("}\n\n");

    let mut total_original = 0usize;
    let mut total_retained = 0usize;

    for main_series in &elp.main {
        let coord = main_series.coordinate;
        let filtered = filter_main_terms(main_series, threshold);
        total_original += main_series.terms.len();
        total_retained += filtered.len();

        out.push_str(&format!(
            "/// Main problem terms for {} ({} of {} terms)\n",
            coord_name(coord),
            filtered.len(),
            main_series.terms.len()
        ));
        out.push_str(&format!(
            "pub const MAIN_{}: &[MainTerm] = &[\n",
            coord_name(coord)
        ));

        for term in &filtered {
            out.push_str(&format!(
                "    MainTerm {{ delaunay: [{}, {}, {}, {}], coeffs: [{:.11}, {:.5}, {:.5}, {:.5}, {:.5}, {:.5}, {:.5}] }},\n",
                term.delaunay[0], term.delaunay[1], term.delaunay[2], term.delaunay[3],
                term.coeffs[0], term.coeffs[1], term.coeffs[2], term.coeffs[3],
                term.coeffs[4], term.coeffs[5], term.coeffs[6]
            ));
        }
        out.push_str("];\n\n");
    }

    for pert_series in &elp.pert {
        let coord = pert_series.coordinate;

        for block in &pert_series.blocks {
            let filtered = filter_pert_terms(block, threshold);
            total_original += block.terms.len();
            total_retained += filtered.len();

            out.push_str(&format!(
                "/// Perturbation terms for {} T^{} ({} of {} terms)\n",
                coord_name(coord),
                block.time_power,
                filtered.len(),
                block.terms.len()
            ));
            out.push_str(&format!(
                "const PERT_{}_T{}: &[PertTerm] = &[\n",
                coord_name(coord),
                block.time_power
            ));

            for term in &filtered {
                let mults: Vec<String> = term.multipliers.iter().map(|m| m.to_string()).collect();
                out.push_str(&format!(
                    "    PertTerm {{ amplitude: {:.13}, phase: {:.13}, multipliers: [{}] }},\n",
                    term.amplitude(),
                    term.phase(),
                    mults.join(", ")
                ));
            }
            out.push_str("];\n\n");
        }

        out.push_str(&format!(
            "/// All perturbation blocks for {}\n",
            coord_name(coord)
        ));
        out.push_str(&format!(
            "pub const PERT_{}: &[PertBlock] = &[\n",
            coord_name(coord)
        ));
        for block in &pert_series.blocks {
            out.push_str(&format!(
                "    PertBlock {{ power: {}, terms: PERT_{}_T{} }},\n",
                block.time_power,
                coord_name(coord),
                block.time_power
            ));
        }
        out.push_str("];\n\n");
    }

    fs::create_dir_all(&config.output_dir)
        .map_err(|e| format!("Failed to create output directory: {}", e))?;

    let output_path = config.output_dir.join("moon.rs");
    let mut file = File::create(&output_path)
        .map_err(|e| format!("Failed to create {}: {}", output_path.display(), e))?;
    file.write_all(out.as_bytes())
        .map_err(|e| format!("Failed to write {}: {}", output_path.display(), e))?;

    println!(
        "Generated {} ({} of {} terms, {:.1}%)",
        output_path.display(),
        total_retained,
        total_original,
        (total_retained as f64 / total_original as f64) * 100.0
    );

    Ok((total_retained, total_original))
}

pub fn print_analysis(elp: &ElpData, threshold: f64) {
    println!("\nELP/MPP02 Analysis (threshold: {:.0e}):", threshold);
    println!("{:-<70}", "");

    println!("\nMain Problem Series:");
    for main_series in &elp.main {
        let filtered = filter_main_terms(main_series, threshold);
        let max_amp = main_series
            .terms
            .iter()
            .map(|t| t.amplitude())
            .fold(0.0f64, f64::max);
        let min_amp = main_series
            .terms
            .iter()
            .map(|t| t.amplitude())
            .fold(f64::MAX, f64::min);

        println!(
            "  {}: {} -> {} terms ({:.1}%), amp range: {:.2e} to {:.2e}",
            coord_name(main_series.coordinate),
            main_series.terms.len(),
            filtered.len(),
            (filtered.len() as f64 / main_series.terms.len() as f64) * 100.0,
            min_amp,
            max_amp
        );
    }

    println!("\nPerturbation Series:");
    for pert_series in &elp.pert {
        println!("  {}:", coord_name(pert_series.coordinate));
        for block in &pert_series.blocks {
            let filtered = filter_pert_terms(block, threshold);
            if block.terms.is_empty() {
                continue;
            }
            let max_amp = block
                .terms
                .iter()
                .map(|t| t.amplitude())
                .fold(0.0f64, f64::max);

            println!(
                "    T^{}: {} -> {} terms ({:.1}%), max amp: {:.2e}",
                block.time_power,
                block.terms.len(),
                filtered.len(),
                (filtered.len() as f64 / block.terms.len() as f64) * 100.0,
                max_amp
            );
        }
    }

    println!("\nTotals:");
    println!("  Main problem: {} terms", elp.total_main_terms());
    println!("  Perturbations: {} terms", elp.total_pert_terms());
    println!("  Total: {} terms", elp.total_terms());
}

#[cfg(not(feature = "cli"))]
pub fn generate_moon_module(
    _elp: &ElpData,
    _config: &GenerateConfig,
) -> Result<(usize, usize), String> {
    Err("Generate requires the 'cli' feature".to_string())
}

#[cfg(test)]
mod tests {
    use super::super::parser::PertSeries;
    use super::*;
    use tempfile::TempDir;

    fn make_main_term(delaunay: [i32; 4], a0: f64) -> MainTerm {
        MainTerm {
            delaunay,
            coeffs: [a0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        }
    }

    fn make_pert_term(sin_c: f64, cos_c: f64, mults: [i32; 16]) -> PertTerm {
        PertTerm {
            sin_coeff: sin_c,
            cos_coeff: cos_c,
            multipliers: mults,
        }
    }

    fn make_test_elp_data() -> ElpData {
        ElpData {
            main: [
                MainSeries {
                    coordinate: Coordinate::Longitude,
                    terms: vec![
                        make_main_term([0, 0, 1, 0], 100.0),
                        make_main_term([2, 0, -1, 0], 50.0),
                        make_main_term([2, 0, 0, 0], 10.0),
                        make_main_term([0, 0, 2, 0], 1.0),
                    ],
                },
                MainSeries {
                    coordinate: Coordinate::Latitude,
                    terms: vec![
                        make_main_term([0, 1, 0, 0], 80.0),
                        make_main_term([0, 1, 1, 0], 5.0),
                    ],
                },
                MainSeries {
                    coordinate: Coordinate::Distance,
                    terms: vec![
                        make_main_term([0, 0, 0, 0], 200.0),
                        make_main_term([2, 0, 0, 0], 20.0),
                    ],
                },
            ],
            pert: [
                PertSeries {
                    coordinate: Coordinate::Longitude,
                    blocks: vec![
                        PertBlock {
                            time_power: 0,
                            terms: vec![
                                make_pert_term(
                                    30.0,
                                    40.0,
                                    [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                ),
                                make_pert_term(
                                    3.0,
                                    4.0,
                                    [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                ),
                            ],
                        },
                        PertBlock {
                            time_power: 1,
                            terms: vec![make_pert_term(
                                6.0,
                                8.0,
                                [1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                            )],
                        },
                    ],
                },
                PertSeries {
                    coordinate: Coordinate::Latitude,
                    blocks: vec![PertBlock {
                        time_power: 0,
                        terms: vec![make_pert_term(
                            12.0,
                            16.0,
                            [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        )],
                    }],
                },
                PertSeries {
                    coordinate: Coordinate::Distance,
                    blocks: vec![PertBlock {
                        time_power: 0,
                        terms: vec![make_pert_term(
                            9.0,
                            12.0,
                            [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        )],
                    }],
                },
            ],
        }
    }

    #[test]
    fn test_coord_name() {
        assert_eq!(coord_name(Coordinate::Longitude), "LONGITUDE");
        assert_eq!(coord_name(Coordinate::Latitude), "LATITUDE");
        assert_eq!(coord_name(Coordinate::Distance), "DISTANCE");
    }

    #[test]
    fn test_coord_var() {
        assert_eq!(coord_var(Coordinate::Longitude), "longitude");
        assert_eq!(coord_var(Coordinate::Latitude), "latitude");
        assert_eq!(coord_var(Coordinate::Distance), "distance");
    }

    #[test]
    fn test_filter_main_terms_filters_by_threshold() {
        let series = MainSeries {
            coordinate: Coordinate::Longitude,
            terms: vec![
                make_main_term([0, 0, 1, 0], 100.0),
                make_main_term([2, 0, -1, 0], 50.0),
                make_main_term([2, 0, 0, 0], 10.0),
                make_main_term([0, 0, 2, 0], 1.0),
            ],
        };

        let filtered = filter_main_terms(&series, 20.0);
        assert_eq!(filtered.len(), 2);
        assert_eq!(filtered[0].amplitude(), 100.0);
        assert_eq!(filtered[1].amplitude(), 50.0);
    }

    #[test]
    fn test_filter_main_terms_sorts_by_amplitude_descending() {
        let series = MainSeries {
            coordinate: Coordinate::Latitude,
            terms: vec![
                make_main_term([0, 0, 1, 0], 10.0),
                make_main_term([2, 0, -1, 0], 100.0),
                make_main_term([2, 0, 0, 0], 50.0),
            ],
        };

        let filtered = filter_main_terms(&series, 0.0);
        assert_eq!(filtered.len(), 3);
        assert_eq!(filtered[0].amplitude(), 100.0);
        assert_eq!(filtered[1].amplitude(), 50.0);
        assert_eq!(filtered[2].amplitude(), 10.0);
    }

    #[test]
    fn test_filter_main_terms_empty_when_all_below_threshold() {
        let series = MainSeries {
            coordinate: Coordinate::Distance,
            terms: vec![
                make_main_term([0, 0, 1, 0], 1.0),
                make_main_term([2, 0, -1, 0], 2.0),
            ],
        };

        let filtered = filter_main_terms(&series, 100.0);
        assert!(filtered.is_empty());
    }

    #[test]
    fn test_filter_main_terms_includes_exact_threshold() {
        let series = MainSeries {
            coordinate: Coordinate::Longitude,
            terms: vec![make_main_term([0, 0, 1, 0], 50.0)],
        };

        let filtered = filter_main_terms(&series, 50.0);
        assert_eq!(filtered.len(), 1);
    }

    #[test]
    fn test_filter_pert_terms_filters_by_threshold() {
        let block = PertBlock {
            time_power: 0,
            terms: vec![
                make_pert_term(30.0, 40.0, [0; 16]), // amplitude = 50
                make_pert_term(3.0, 4.0, [0; 16]),   // amplitude = 5
                make_pert_term(0.6, 0.8, [0; 16]),   // amplitude = 1
            ],
        };

        let filtered = filter_pert_terms(&block, 10.0);
        assert_eq!(filtered.len(), 1);
        assert!((filtered[0].amplitude() - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_filter_pert_terms_sorts_by_amplitude_descending() {
        let block = PertBlock {
            time_power: 1,
            terms: vec![
                make_pert_term(3.0, 4.0, [0; 16]),   // amplitude = 5
                make_pert_term(30.0, 40.0, [0; 16]), // amplitude = 50
                make_pert_term(6.0, 8.0, [0; 16]),   // amplitude = 10
            ],
        };

        let filtered = filter_pert_terms(&block, 0.0);
        assert_eq!(filtered.len(), 3);
        assert!((filtered[0].amplitude() - 50.0).abs() < 1e-10);
        assert!((filtered[1].amplitude() - 10.0).abs() < 1e-10);
        assert!((filtered[2].amplitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_filter_pert_terms_empty_block() {
        let block = PertBlock {
            time_power: 0,
            terms: vec![],
        };

        let filtered = filter_pert_terms(&block, 0.0);
        assert!(filtered.is_empty());
    }

    #[test]
    fn test_generate_config_struct() {
        let config = GenerateConfig {
            threshold: 1e-5,
            output_dir: std::env::temp_dir().join("test"),
        };
        assert_eq!(config.threshold, 1e-5);
        assert!(config.output_dir.ends_with("test"));
    }

    #[test]
    fn test_print_analysis_runs_without_panic() {
        let elp = make_test_elp_data();
        print_analysis(&elp, 1e-5);
    }

    #[test]
    fn test_print_analysis_with_empty_pert_block() {
        let elp = ElpData {
            main: [
                MainSeries {
                    coordinate: Coordinate::Longitude,
                    terms: vec![make_main_term([0, 0, 1, 0], 100.0)],
                },
                MainSeries {
                    coordinate: Coordinate::Latitude,
                    terms: vec![make_main_term([0, 1, 0, 0], 80.0)],
                },
                MainSeries {
                    coordinate: Coordinate::Distance,
                    terms: vec![make_main_term([0, 0, 0, 0], 200.0)],
                },
            ],
            pert: [
                PertSeries {
                    coordinate: Coordinate::Longitude,
                    blocks: vec![PertBlock {
                        time_power: 0,
                        terms: vec![],
                    }],
                },
                PertSeries {
                    coordinate: Coordinate::Latitude,
                    blocks: vec![],
                },
                PertSeries {
                    coordinate: Coordinate::Distance,
                    blocks: vec![],
                },
            ],
        };
        print_analysis(&elp, 1e-5);
    }

    #[test]
    fn test_print_analysis_with_high_threshold() {
        let elp = make_test_elp_data();
        print_analysis(&elp, 1e10);
    }

    #[cfg(feature = "cli")]
    mod cli_tests {
        use super::*;

        #[test]
        fn test_generate_moon_module_creates_file() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 1.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            let result = generate_moon_module(&elp, &config);
            assert!(result.is_ok());

            let output_file = temp_dir.path().join("moon.rs");
            assert!(output_file.exists());
        }

        #[test]
        fn test_generate_moon_module_returns_term_counts() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            let (retained, original) = generate_moon_module(&elp, &config).unwrap();
            assert_eq!(original, elp.total_terms());
            assert_eq!(retained, original);
        }

        #[test]
        fn test_generate_moon_module_filters_terms() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 20.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            let (retained, original) = generate_moon_module(&elp, &config).unwrap();
            assert!(retained < original);
        }

        #[test]
        fn test_generate_moon_module_creates_nested_directory() {
            let temp_dir = TempDir::new().unwrap();
            let nested_path = temp_dir.path().join("nested").join("deep").join("path");
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 1.0,
                output_dir: nested_path.clone(),
            };

            let result = generate_moon_module(&elp, &config);
            assert!(result.is_ok());
            assert!(nested_path.join("moon.rs").exists());
        }

        #[test]
        fn test_generate_moon_module_output_contains_structs() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("pub struct MainTerm"));
            assert!(content.contains("pub struct PertTerm"));
            assert!(content.contains("pub struct PertBlock"));
        }

        #[test]
        fn test_generate_moon_module_output_contains_main_constants() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("pub const MAIN_LONGITUDE:"));
            assert!(content.contains("pub const MAIN_LATITUDE:"));
            assert!(content.contains("pub const MAIN_DISTANCE:"));
        }

        #[test]
        fn test_generate_moon_module_output_contains_pert_constants() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("pub const PERT_LONGITUDE:"));
            assert!(content.contains("pub const PERT_LATITUDE:"));
            assert!(content.contains("pub const PERT_DISTANCE:"));
        }

        #[test]
        fn test_generate_moon_module_output_contains_clippy_allows() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("#![allow(clippy::excessive_precision)]"));
            assert!(content.contains("#![allow(clippy::unreadable_literal)]"));
            assert!(content.contains("#![allow(clippy::approx_constant)]"));
        }

        #[test]
        fn test_generate_moon_module_output_contains_reference_comment() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("Chapront & Francou (2003)"));
            assert!(content.contains("ELP/MPP02 coefficients for the Moon"));
        }

        #[test]
        fn test_generate_moon_module_output_contains_threshold() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 1e-5,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("Threshold: 1e-5"));
        }

        #[test]
        fn test_generate_moon_module_output_contains_term_data() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("MainTerm { delaunay:"));
            assert!(content.contains("PertTerm { amplitude:"));
            assert!(content.contains("PertBlock { power:"));
        }

        #[test]
        fn test_generate_moon_module_term_counts_in_comments() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("(4 of 4 terms)"));
            assert!(content.contains("(2 of 2 terms)"));
        }

        #[test]
        fn test_generate_moon_module_pert_time_power_in_output() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("PERT_LONGITUDE_T0"));
            assert!(content.contains("PERT_LONGITUDE_T1"));
            assert!(content.contains("T^0"));
            assert!(content.contains("T^1"));
        }

        #[test]
        fn test_generate_moon_module_overwrites_existing_file() {
            let temp_dir = TempDir::new().unwrap();
            let output_file = temp_dir.path().join("moon.rs");
            std::fs::write(&output_file, "old content").unwrap();

            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(&output_file).unwrap();
            assert!(!content.contains("old content"));
            assert!(content.contains("ELP/MPP02 coefficients"));
        }

        #[test]
        fn test_generate_moon_module_with_high_threshold_filters_all() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 1e10,
                output_dir: temp_dir.path().to_path_buf(),
            };

            let (retained, original) = generate_moon_module(&elp, &config).unwrap();
            assert_eq!(retained, 0);
            assert!(original > 0);
        }

        #[test]
        fn test_generate_moon_module_output_is_valid_rust_syntax() {
            let temp_dir = TempDir::new().unwrap();
            let elp = make_test_elp_data();
            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();

            let open_braces = content.matches('{').count();
            let close_braces = content.matches('}').count();
            assert_eq!(open_braces, close_braces, "Mismatched braces");

            let open_brackets = content.matches('[').count();
            let close_brackets = content.matches(']').count();
            assert_eq!(open_brackets, close_brackets, "Mismatched brackets");

            assert!(content.contains("&[MainTerm]"));
            assert!(content.contains("&[PertTerm]"));
            assert!(content.contains("&[PertBlock]"));
        }

        #[test]
        fn test_generate_moon_module_delaunay_values_in_output() {
            let temp_dir = TempDir::new().unwrap();
            let mut elp = make_test_elp_data();
            elp.main[0].terms = vec![make_main_term([1, -2, 3, -4], 100.0)];

            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("[1, -2, 3, -4]"));
        }

        #[test]
        fn test_generate_moon_module_pert_multipliers_in_output() {
            let temp_dir = TempDir::new().unwrap();
            let mut elp = make_test_elp_data();
            let mults = [1, 2, 3, 4, 5, 6, 7, 8, -1, -2, -3, -4, 0, 0, 0, 0];
            elp.pert[0].blocks[0].terms = vec![make_pert_term(30.0, 40.0, mults)];

            let config = GenerateConfig {
                threshold: 0.0,
                output_dir: temp_dir.path().to_path_buf(),
            };

            generate_moon_module(&elp, &config).unwrap();

            let content = std::fs::read_to_string(temp_dir.path().join("moon.rs")).unwrap();
            assert!(content.contains("1, 2, 3, 4, 5, 6, 7, 8, -1, -2, -3, -4, 0, 0, 0, 0"));
        }
    }

    #[cfg(not(feature = "cli"))]
    #[test]
    fn test_generate_moon_module_without_cli_feature() {
        let elp = make_test_elp_data();
        let config = GenerateConfig {
            threshold: 1.0,
            output_dir: std::env::temp_dir().join("test"),
        };

        let result = generate_moon_module(&elp, &config);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Generate requires the 'cli' feature");
    }
}
