use super::chebyshev::evaluate_position_velocity;
use super::daf::{DafFile, DafSummary};
use super::SpkError;
use celestial_core::constants::{J2000_JD, SECONDS_PER_DAY_F64};
use std::path::Path;

fn jd_to_seconds_from_j2000(jd: f64) -> f64 {
    (jd - J2000_JD) * SECONDS_PER_DAY_F64
}

fn seconds_from_j2000_to_jd(seconds: f64) -> f64 {
    J2000_JD + seconds / SECONDS_PER_DAY_F64
}

#[derive(Debug, Clone)]
pub struct SpkSegment {
    pub body_id: i32,
    pub center_id: i32,
    pub frame_id: i32,
    pub data_type: i32,
    pub start_epoch: f64,
    pub end_epoch: f64,
    pub start_index: usize,
    pub end_index: usize,
}

impl SpkSegment {
    pub fn from_summary(summary: &DafSummary) -> Result<Self, SpkError> {
        if summary.doubles.len() < 2 || summary.ints.len() < 6 {
            return Err(SpkError::InvalidData("Incomplete SPK summary".into()));
        }
        Ok(Self {
            start_epoch: summary.doubles[0],
            end_epoch: summary.doubles[1],
            body_id: summary.ints[0],
            center_id: summary.ints[1],
            frame_id: summary.ints[2],
            data_type: summary.ints[3],
            start_index: summary.ints[4] as usize,
            end_index: summary.ints[5] as usize,
        })
    }

    pub fn contains_epoch(&self, jd_tdb: f64) -> bool {
        let epoch_seconds = jd_to_seconds_from_j2000(jd_tdb);
        epoch_seconds >= self.start_epoch && epoch_seconds <= self.end_epoch
    }

    pub fn start_jd(&self) -> f64 {
        seconds_from_j2000_to_jd(self.start_epoch)
    }

    pub fn end_jd(&self) -> f64 {
        seconds_from_j2000_to_jd(self.end_epoch)
    }
}

struct Type2Metadata {
    init: f64,
    intlen: f64,
    rsize: usize,
    n_records: usize,
    n_coeffs: usize,
}

pub struct SpkFile {
    daf: DafFile,
    segments: Vec<SpkSegment>,
}

impl SpkFile {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, SpkError> {
        let daf = DafFile::open(path)?;
        let segments = Self::load_segments(&daf)?;
        Ok(Self { daf, segments })
    }

    fn load_segments(daf: &DafFile) -> Result<Vec<SpkSegment>, SpkError> {
        let mut segments = Vec::new();
        for summary_result in daf.iter_summaries() {
            let summary = summary_result?;
            let segment = SpkSegment::from_summary(&summary)?;
            if segment.data_type != 2 {
                continue;
            }
            segments.push(segment);
        }
        Ok(segments)
    }

    pub fn segments(&self) -> &[SpkSegment] {
        &self.segments
    }

    pub fn find_segment(&self, body: i32, center: i32, jd_tdb: f64) -> Option<&SpkSegment> {
        self.segments
            .iter()
            .find(|s| s.body_id == body && s.center_id == center && s.contains_epoch(jd_tdb))
    }

    fn read_type2_metadata(&self, segment: &SpkSegment) -> Result<Type2Metadata, SpkError> {
        let meta = self.daf.read_f64_array(segment.end_index - 3, 4)?;
        let init = meta[0];
        let intlen = meta[1];
        let rsize = meta[2] as usize;
        let n_records = meta[3] as usize;
        let n_coeffs = (rsize - 2) / 3;
        Ok(Type2Metadata {
            init,
            intlen,
            rsize,
            n_records,
            n_coeffs,
        })
    }

    fn find_record_index(&self, meta: &Type2Metadata, jd_tdb: f64) -> usize {
        let epoch_seconds = jd_to_seconds_from_j2000(jd_tdb);
        let seconds_from_init = epoch_seconds - meta.init;
        let record_index = libm::floor(seconds_from_init / meta.intlen) as usize;
        record_index.min(meta.n_records - 1)
    }

    fn read_type2_record(
        &self,
        segment: &SpkSegment,
        meta: &Type2Metadata,
        record_index: usize,
    ) -> Result<(Vec<f64>, f64, f64), SpkError> {
        let record_start = segment.start_index + record_index * meta.rsize;
        let record = self.daf.read_f64_array(record_start, meta.rsize)?;
        let mid = record[0];
        let radius = record[1];
        Ok((record, mid, radius))
    }

    pub fn compute_position(
        &self,
        body: i32,
        center: i32,
        jd_tdb: f64,
    ) -> Result<[f64; 3], SpkError> {
        let (pos, _) = self.compute_state(body, center, jd_tdb)?;
        Ok(pos)
    }

    pub fn compute_state(
        &self,
        body: i32,
        center: i32,
        jd_tdb: f64,
    ) -> Result<([f64; 3], [f64; 3]), SpkError> {
        let segment = self.find_segment(body, center, jd_tdb).ok_or({
            SpkError::SegmentNotFound {
                body,
                center,
                epoch: jd_tdb,
            }
        })?;
        let meta = self.read_type2_metadata(segment)?;
        let record_index = self.find_record_index(&meta, jd_tdb);
        let (record, mid, radius) = self.read_type2_record(segment, &meta, record_index)?;
        let epoch_seconds = jd_to_seconds_from_j2000(jd_tdb);
        let t_normalized = (epoch_seconds - mid) / radius;
        let coeffs = &record[2..];
        let n = meta.n_coeffs;
        let coeffs_x = &coeffs[0..n];
        let coeffs_y = &coeffs[n..2 * n];
        let coeffs_z = &coeffs[2 * n..3 * n];
        evaluate_position_velocity(coeffs_x, coeffs_y, coeffs_z, n, t_normalized, radius)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segment_contains_epoch() {
        let start_jd = celestial_core::constants::J2000_JD;
        let end_jd = 2451645.0;
        let segment = SpkSegment {
            body_id: 3,
            center_id: 0,
            frame_id: 1,
            data_type: 2,
            start_epoch: jd_to_seconds_from_j2000(start_jd),
            end_epoch: jd_to_seconds_from_j2000(end_jd),
            start_index: 0,
            end_index: 0,
        };
        assert!(segment.contains_epoch(2451545.0));
        assert!(segment.contains_epoch(2451600.0));
        assert!(segment.contains_epoch(2451645.0));
        assert!(!segment.contains_epoch(2451544.9));
        assert!(!segment.contains_epoch(2451645.1));
        assert_eq!(segment.start_jd(), start_jd);
        assert_eq!(segment.end_jd(), end_jd);
    }

    #[test]
    fn test_epoch_conversions() {
        let jd = 2451545.0;
        let seconds = jd_to_seconds_from_j2000(jd);
        assert_eq!(seconds, 0.0);
        assert_eq!(seconds_from_j2000_to_jd(seconds), jd);
        let jd2 = 2460000.5;
        let seconds2 = jd_to_seconds_from_j2000(jd2);
        assert_eq!(seconds_from_j2000_to_jd(seconds2), jd2);
        let expected_seconds = (jd2 - J2000_JD) * SECONDS_PER_DAY_F64;
        assert_eq!(seconds2, expected_seconds);
    }

    #[test]
    fn test_spk_segment_from_summary() {
        let summary = DafSummary {
            doubles: vec![0.0, celestial_core::constants::SECONDS_PER_DAY_F64], // start and end epoch
            ints: vec![399, 0, 1, 2, 100, 200], // body, center, frame, type, start_idx, end_idx
        };

        let segment = SpkSegment::from_summary(&summary).unwrap();
        assert_eq!(segment.body_id, 399);
        assert_eq!(segment.center_id, 0);
        assert_eq!(segment.frame_id, 1);
        assert_eq!(segment.data_type, 2);
        assert_eq!(segment.start_index, 100);
        assert_eq!(segment.end_index, 200);
        assert!((segment.start_epoch - 0.0).abs() < 1e-10);
        assert!((segment.end_epoch - celestial_core::constants::SECONDS_PER_DAY_F64).abs() < 1e-10);
    }

    #[test]
    fn test_spk_segment_from_summary_insufficient_doubles() {
        let summary = DafSummary {
            doubles: vec![0.0], // only 1 double, need 2
            ints: vec![399, 0, 1, 2, 100, 200],
        };

        let result = SpkSegment::from_summary(&summary);
        assert!(result.is_err());
        match result {
            Err(SpkError::InvalidData(msg)) => assert!(msg.contains("Incomplete")),
            _ => panic!("Expected InvalidData error"),
        }
    }

    #[test]
    fn test_spk_segment_from_summary_insufficient_ints() {
        let summary = DafSummary {
            doubles: vec![0.0, celestial_core::constants::SECONDS_PER_DAY_F64],
            ints: vec![399, 0, 1, 2, 100], // only 5 ints, need 6
        };

        let result = SpkSegment::from_summary(&summary);
        assert!(result.is_err());
        match result {
            Err(SpkError::InvalidData(msg)) => assert!(msg.contains("Incomplete")),
            _ => panic!("Expected InvalidData error"),
        }
    }

    #[test]
    fn test_spk_segment_clone() {
        let segment = SpkSegment {
            body_id: 399,
            center_id: 0,
            frame_id: 1,
            data_type: 2,
            start_epoch: 0.0,
            end_epoch: celestial_core::constants::SECONDS_PER_DAY_F64,
            start_index: 100,
            end_index: 200,
        };

        let cloned = segment.clone();
        assert_eq!(cloned.body_id, segment.body_id);
        assert_eq!(cloned.center_id, segment.center_id);
        assert_eq!(cloned.frame_id, segment.frame_id);
        assert_eq!(cloned.data_type, segment.data_type);
    }

    #[test]
    fn test_spk_segment_debug() {
        let segment = SpkSegment {
            body_id: 399,
            center_id: 0,
            frame_id: 1,
            data_type: 2,
            start_epoch: 0.0,
            end_epoch: celestial_core::constants::SECONDS_PER_DAY_F64,
            start_index: 100,
            end_index: 200,
        };

        let debug = format!("{:?}", segment);
        assert!(debug.contains("SpkSegment"));
        assert!(debug.contains("body_id: 399"));
    }

    #[test]
    fn test_type2_metadata_struct() {
        let meta = Type2Metadata {
            init: 0.0,
            intlen: celestial_core::constants::SECONDS_PER_DAY_F64,
            rsize: 50,
            n_records: 100,
            n_coeffs: 16,
        };

        assert_eq!(meta.init, 0.0);
        assert_eq!(meta.intlen, celestial_core::constants::SECONDS_PER_DAY_F64);
        assert_eq!(meta.rsize, 50);
        assert_eq!(meta.n_records, 100);
        assert_eq!(meta.n_coeffs, 16);
    }

    #[test]
    fn test_jd_to_seconds_roundtrip() {
        let test_jds = [
            2451545.0, // J2000.0
            2451545.5, // J2000.0 + 12 hours
            2460000.0, // future date
            2440000.0, // past date
        ];

        for &jd in &test_jds {
            let seconds = jd_to_seconds_from_j2000(jd);
            let reconstructed = seconds_from_j2000_to_jd(seconds);
            assert!(
                (jd - reconstructed).abs() < 1e-10,
                "Failed for JD {}: got {}",
                jd,
                reconstructed
            );
        }
    }

    #[test]
    fn test_spk_file_open_nonexistent() {
        let result = SpkFile::open("/nonexistent/path/file.bsp");
        assert!(result.is_err());
    }

    fn get_de432s_path() -> Option<std::path::PathBuf> {
        let test_file =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/de432s.bsp");
        if test_file.exists() {
            Some(test_file)
        } else {
            None
        }
    }

    #[test]
    fn test_spk_segments_accessor() {
        let path = match get_de432s_path() {
            Some(p) => p,
            None => {
                eprintln!("Skipping: de432s.bsp not found");
                return;
            }
        };
        let spk = SpkFile::open(&path).unwrap();
        let segments = spk.segments();
        assert!(!segments.is_empty());
    }

    #[test]
    fn test_spk_find_segment() {
        let path = match get_de432s_path() {
            Some(p) => p,
            None => {
                eprintln!("Skipping: de432s.bsp not found");
                return;
            }
        };
        let spk = SpkFile::open(&path).unwrap();

        // Earth-Moon barycenter (3) relative to SSB (0) at J2000
        let seg = spk.find_segment(3, 0, J2000_JD);
        assert!(
            seg.is_some(),
            "Should find segment for Earth-Moon barycenter"
        );

        // Non-existent body
        let seg = spk.find_segment(999, 0, J2000_JD);
        assert!(seg.is_none(), "Should not find segment for body 999");
    }

    #[test]
    fn test_spk_compute_position() {
        let path = match get_de432s_path() {
            Some(p) => p,
            None => {
                eprintln!("Skipping: de432s.bsp not found");
                return;
            }
        };
        let spk = SpkFile::open(&path).unwrap();

        // Compute Earth-Moon barycenter position at J2000
        let result = spk.compute_position(3, 0, J2000_JD);
        assert!(result.is_ok(), "Should compute position");

        let pos = result.unwrap();
        let dist = libm::sqrt(pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2));
        // Earth-Moon barycenter should be ~1 AU from SSB (in km)
        let dist_au = dist / 149597870.7;
        assert!(
            dist_au > 0.98 && dist_au < 1.02,
            "Earth should be ~1 AU from SSB, got {} AU",
            dist_au
        );
    }

    #[test]
    fn test_spk_compute_state() {
        let path = match get_de432s_path() {
            Some(p) => p,
            None => {
                eprintln!("Skipping: de432s.bsp not found");
                return;
            }
        };
        let spk = SpkFile::open(&path).unwrap();

        // Compute Earth-Moon barycenter state at J2000
        let result = spk.compute_state(3, 0, J2000_JD);
        assert!(result.is_ok(), "Should compute state");

        let (pos, vel) = result.unwrap();

        // Check position is reasonable
        let dist = libm::sqrt(pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2));
        assert!(dist > 1e8, "Position magnitude should be > 100 million km");

        // Check velocity is reasonable (Earth orbital velocity ~30 km/s)
        let vel_mag = libm::sqrt(vel[0].powi(2) + vel[1].powi(2) + vel[2].powi(2));
        assert!(
            vel_mag > 20.0 && vel_mag < 40.0,
            "Velocity should be ~30 km/s, got {} km/s",
            vel_mag
        );
    }

    #[test]
    fn test_spk_segment_not_found_error() {
        let path = match get_de432s_path() {
            Some(p) => p,
            None => {
                eprintln!("Skipping: de432s.bsp not found");
                return;
            }
        };
        let spk = SpkFile::open(&path).unwrap();

        // Request a non-existent body
        let result = spk.compute_position(999, 0, J2000_JD);
        assert!(result.is_err());
        match result {
            Err(SpkError::SegmentNotFound {
                body,
                center,
                epoch,
            }) => {
                assert_eq!(body, 999);
                assert_eq!(center, 0);
                assert!((epoch - J2000_JD).abs() < 1e-10);
            }
            _ => panic!("Expected SegmentNotFound error"),
        }
    }

    #[test]
    fn test_spk_compute_position_different_epochs() {
        let path = match get_de432s_path() {
            Some(p) => p,
            None => {
                eprintln!("Skipping: de432s.bsp not found");
                return;
            }
        };
        let spk = SpkFile::open(&path).unwrap();

        // Compute at different epochs
        let epochs = [J2000_JD, J2000_JD + 100.0, J2000_JD + 365.25];
        let mut prev_pos: Option<[f64; 3]> = None;

        for epoch in epochs {
            let result = spk.compute_position(3, 0, epoch);
            assert!(result.is_ok(), "Should compute position at epoch {}", epoch);

            let pos = result.unwrap();
            if let Some(prev) = prev_pos {
                // Positions should be different at different epochs
                let diff = libm::sqrt(
                    (pos[0] - prev[0]).powi(2)
                        + (pos[1] - prev[1]).powi(2)
                        + (pos[2] - prev[2]).powi(2),
                );
                assert!(
                    diff > 1e6,
                    "Positions should differ significantly between epochs"
                );
            }
            prev_pos = Some(pos);
        }
    }
}
