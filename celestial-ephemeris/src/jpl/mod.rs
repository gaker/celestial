mod chebyshev;
mod daf;
mod spk;

pub use spk::{SpkFile, SpkSegment};

#[derive(Debug, Clone, PartialEq)]
pub enum SpkError {
    Io(String),
    InvalidFormat(String),
    InvalidData(String),
    SegmentNotFound { body: i32, center: i32, epoch: f64 },
    UnsupportedType(i32),
}

impl std::fmt::Display for SpkError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpkError::Io(msg) => write!(f, "IO error: {}", msg),
            SpkError::InvalidFormat(msg) => write!(f, "Invalid SPK format: {}", msg),
            SpkError::InvalidData(msg) => write!(f, "Invalid SPK data: {}", msg),
            SpkError::SegmentNotFound {
                body,
                center,
                epoch,
            } => {
                write!(
                    f,
                    "No segment found for body {} relative to {} at JD {}",
                    body, center, epoch
                )
            }
            SpkError::UnsupportedType(t) => write!(f, "Unsupported SPK type: {}", t),
        }
    }
}

impl std::error::Error for SpkError {}

pub mod bodies {
    pub const SOLAR_SYSTEM_BARYCENTER: i32 = 0;
    pub const MERCURY_BARYCENTER: i32 = 1;
    pub const VENUS_BARYCENTER: i32 = 2;
    pub const EARTH_MOON_BARYCENTER: i32 = 3;
    pub const MARS_BARYCENTER: i32 = 4;
    pub const JUPITER_BARYCENTER: i32 = 5;
    pub const SATURN_BARYCENTER: i32 = 6;
    pub const URANUS_BARYCENTER: i32 = 7;
    pub const NEPTUNE_BARYCENTER: i32 = 8;
    pub const PLUTO_BARYCENTER: i32 = 9;
    pub const SUN: i32 = 10;
    pub const MOON: i32 = 301;
    pub const EARTH: i32 = 399;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_spk_error_display_io() {
        let err = SpkError::Io("file not found".to_string());
        let display = format!("{}", err);
        assert!(display.contains("IO error"));
        assert!(display.contains("file not found"));
    }

    #[test]
    fn test_spk_error_display_invalid_format() {
        let err = SpkError::InvalidFormat("bad format".to_string());
        let display = format!("{}", err);
        assert!(display.contains("Invalid SPK format"));
        assert!(display.contains("bad format"));
    }

    #[test]
    fn test_spk_error_display_invalid_data() {
        let err = SpkError::InvalidData("corrupted data".to_string());
        let display = format!("{}", err);
        assert!(display.contains("Invalid SPK data"));
        assert!(display.contains("corrupted data"));
    }

    #[test]
    fn test_spk_error_display_segment_not_found() {
        let err = SpkError::SegmentNotFound {
            body: 399,
            center: 0,
            epoch: 2451545.0,
        };
        let display = format!("{}", err);
        assert!(display.contains("No segment found"));
        assert!(display.contains("body 399"));
        assert!(display.contains("relative to 0"));
        assert!(display.contains("2451545"));
    }

    #[test]
    fn test_spk_error_display_unsupported_type() {
        let err = SpkError::UnsupportedType(99);
        let display = format!("{}", err);
        assert!(display.contains("Unsupported SPK type"));
        assert!(display.contains("99"));
    }

    #[test]
    fn test_spk_error_debug() {
        let err = SpkError::Io("test".to_string());
        let debug = format!("{:?}", err);
        assert!(debug.contains("Io"));
    }

    #[test]
    fn test_spk_error_clone() {
        let err = SpkError::InvalidFormat("test".to_string());
        let cloned = err.clone();
        assert_eq!(err, cloned);
    }

    #[test]
    fn test_spk_error_partial_eq() {
        let err1 = SpkError::Io("test".to_string());
        let err2 = SpkError::Io("test".to_string());
        let err3 = SpkError::Io("different".to_string());

        assert_eq!(err1, err2);
        assert_ne!(err1, err3);
    }

    #[test]
    fn test_spk_error_is_error() {
        let err: Box<dyn std::error::Error> = Box::new(SpkError::Io("test".to_string()));
        assert!(!err.to_string().is_empty());
    }

    #[test]
    fn test_bodies_constants() {
        assert_eq!(bodies::SOLAR_SYSTEM_BARYCENTER, 0);
        assert_eq!(bodies::SUN, 10);
        assert_eq!(bodies::EARTH, 399);
        assert_eq!(bodies::MOON, 301);
        assert_eq!(bodies::EARTH_MOON_BARYCENTER, 3);
    }

    fn get_de432s_path() -> Option<PathBuf> {
        let test_file =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/de432s.bsp");
        if test_file.exists() {
            Some(test_file)
        } else {
            None
        }
    }

    #[test]
    fn test_de432s_segment_count() {
        let path = match get_de432s_path() {
            Some(p) => p,
            None => {
                eprintln!("Skipping: de432s.bsp not found");
                return;
            }
        };
        let spk = SpkFile::open(&path).expect("Failed to open de432s.bsp");
        let segments = spk.segments();
        println!("Found {} segments in de432s.bsp:", segments.len());
        for seg in segments {
            println!(
                "  Body {:3} -> Center {:3}, type {}, JD {:.1} to {:.1}",
                seg.body_id,
                seg.center_id,
                seg.data_type,
                seg.start_jd(),
                seg.end_jd()
            );
        }
        assert!(
            segments.len() >= 14,
            "de432s.bsp should have at least 14 segments, found {}",
            segments.len()
        );
    }

    fn get_de440_path() -> Option<PathBuf> {
        if let Ok(path) = std::env::var("DE440_PATH") {
            let p = PathBuf::from(path);
            if p.exists() {
                return Some(p);
            }
        }
        let home = std::env::var("HOME").ok()?;
        let candidates = [
            format!("{}/.local/share/ephemeris/de440.bsp", home),
            format!("{}/ephemeris/de440.bsp", home),
            "/usr/local/share/ephemeris/de440.bsp".to_string(),
        ];
        for path in candidates {
            let p = PathBuf::from(&path);
            if p.exists() {
                return Some(p);
            }
        }
        None
    }

    #[test]
    #[ignore]
    fn test_open_de440() {
        let path = get_de440_path().expect("DE440 file not found - set DE440_PATH env var");
        let spk = SpkFile::open(&path).expect("Failed to open DE440");
        let segments = spk.segments();
        assert!(!segments.is_empty(), "No segments found in DE440");
        println!("Found {} segments", segments.len());
        for seg in segments.iter().take(10) {
            println!(
                "Body {} -> Center {}, JD {:.1} to {:.1}",
                seg.body_id,
                seg.center_id,
                seg.start_jd(),
                seg.end_jd()
            );
        }
    }

    #[test]
    #[ignore]
    fn test_compute_earth_position() {
        let path = get_de440_path().expect("DE440 file not found");
        let spk = SpkFile::open(&path).expect("Failed to open DE440");
        let jd_j2000 = 2451545.0;
        let (pos, vel) = spk
            .compute_state(
                bodies::EARTH_MOON_BARYCENTER,
                bodies::SOLAR_SYSTEM_BARYCENTER,
                jd_j2000,
            )
            .expect("Failed to compute Earth state");
        let distance_au =
            libm::sqrt(pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2)) / 149597870.7;
        println!("Earth-Moon barycenter at J2000.0:");
        println!(
            "  Position: [{:.3}, {:.3}, {:.3}] km",
            pos[0], pos[1], pos[2]
        );
        println!(
            "  Velocity: [{:.6}, {:.6}, {:.6}] km/s",
            vel[0], vel[1], vel[2]
        );
        println!("  Distance from SSB: {:.6} AU", distance_au);
        assert!(
            distance_au > 0.98 && distance_au < 1.02,
            "Earth should be ~1 AU from SSB"
        );
    }

    #[test]
    #[ignore]
    fn test_compute_mars_position() {
        let path = get_de440_path().expect("DE440 file not found");
        let spk = SpkFile::open(&path).expect("Failed to open DE440");
        let jd = 2460000.5;
        let (pos, vel) = spk
            .compute_state(bodies::MARS_BARYCENTER, bodies::SOLAR_SYSTEM_BARYCENTER, jd)
            .expect("Failed to compute Mars state");
        let distance_au =
            libm::sqrt(pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2)) / 149597870.7;
        println!("Mars barycenter at JD {}:", jd);
        println!(
            "  Position: [{:.3}, {:.3}, {:.3}] km",
            pos[0], pos[1], pos[2]
        );
        println!(
            "  Velocity: [{:.6}, {:.6}, {:.6}] km/s",
            vel[0], vel[1], vel[2]
        );
        println!("  Distance from SSB: {:.6} AU", distance_au);
        assert!(
            distance_au > 1.3 && distance_au < 1.7,
            "Mars should be 1.4-1.7 AU from SSB"
        );
    }
}
