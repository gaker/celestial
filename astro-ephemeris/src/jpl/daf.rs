use super::SpkError;
use memmap2::Mmap;
use std::fs::File;
use std::path::Path;

const DAF_RECORD_SIZE: usize = 1024;
const FTPSTR: &[u8] = b"FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP";

pub struct DafFile {
    pub mmap: Mmap,
    pub endian: Endian,
    pub nd: usize,
    pub ni: usize,
    pub summary_size: usize,
    pub fward: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Endian {
    Little,
    Big,
}

impl Endian {
    pub fn read_f64(&self, bytes: &[u8]) -> f64 {
        let arr: [u8; 8] = bytes[..8].try_into().unwrap();
        match self {
            Endian::Little => f64::from_le_bytes(arr),
            Endian::Big => f64::from_be_bytes(arr),
        }
    }

    pub fn read_i32(&self, bytes: &[u8]) -> i32 {
        let arr: [u8; 4] = bytes[..4].try_into().unwrap();
        match self {
            Endian::Little => i32::from_le_bytes(arr),
            Endian::Big => i32::from_be_bytes(arr),
        }
    }
}

impl DafFile {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, SpkError> {
        let file = File::open(path.as_ref()).map_err(|e| SpkError::Io(e.to_string()))?;
        let mmap = unsafe { Mmap::map(&file).map_err(|e| SpkError::Io(e.to_string()))? };
        Self::from_mmap(mmap)
    }

    fn from_mmap(mmap: Mmap) -> Result<Self, SpkError> {
        if mmap.len() < DAF_RECORD_SIZE {
            return Err(SpkError::InvalidFormat("File too small for DAF".into()));
        }
        let locidw = &mmap[0..8];
        if !locidw.starts_with(b"DAF/") {
            return Err(SpkError::InvalidFormat(format!(
                "Invalid DAF signature: {:?}",
                String::from_utf8_lossy(locidw)
            )));
        }
        let endian = Self::detect_endian(&mmap)?;
        let nd = endian.read_i32(&mmap[8..12]) as usize;
        let ni = endian.read_i32(&mmap[12..16]) as usize;
        let fward = endian.read_i32(&mmap[76..80]) as usize;
        let summary_size = nd + ni.div_ceil(2);
        Self::verify_ftp(&mmap)?;
        Ok(Self {
            mmap,
            endian,
            nd,
            ni,
            summary_size,
            fward,
        })
    }

    fn detect_endian(mmap: &Mmap) -> Result<Endian, SpkError> {
        let nd_le = i32::from_le_bytes(mmap[8..12].try_into().unwrap());
        let nd_be = i32::from_be_bytes(mmap[8..12].try_into().unwrap());
        if (1..=100).contains(&nd_le) {
            Ok(Endian::Little)
        } else if (1..=100).contains(&nd_be) {
            Ok(Endian::Big)
        } else {
            Err(SpkError::InvalidFormat(
                "Cannot determine endianness".into(),
            ))
        }
    }

    fn verify_ftp(mmap: &Mmap) -> Result<(), SpkError> {
        if mmap.len() >= 1000 && &mmap[699..727] != FTPSTR {
            return Err(SpkError::InvalidFormat("FTP corruption detected".into()));
        }
        Ok(())
    }

    pub fn iter_summaries(&self) -> SummaryIterator<'_> {
        SummaryIterator {
            daf: self,
            record_num: self.fward,
            record_offset: 0,
            next_record: 0,
            index: 0,
            count: 0,
        }
    }

    pub fn read_array(&self, start: usize, end: usize) -> Result<&[u8], SpkError> {
        let byte_start = (start - 1) * 8;
        let byte_end = end * 8;
        if byte_end > self.mmap.len() {
            return Err(SpkError::InvalidData("Array range out of bounds".into()));
        }
        Ok(&self.mmap[byte_start..byte_end])
    }

    pub fn read_f64_array(&self, start: usize, count: usize) -> Result<Vec<f64>, SpkError> {
        let bytes = self.read_array(start, start + count - 1)?;
        let mut result = Vec::with_capacity(count);
        for i in 0..count {
            result.push(self.endian.read_f64(&bytes[i * 8..]));
        }
        Ok(result)
    }
}

pub struct SummaryIterator<'a> {
    daf: &'a DafFile,
    record_num: usize,
    record_offset: usize,
    next_record: usize,
    index: usize,
    count: usize,
}

impl<'a> Iterator for SummaryIterator<'a> {
    type Item = Result<DafSummary, SpkError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.index < self.count {
                let summary_start = 24 + self.index * self.daf.summary_size * 8;
                let offset = self.record_offset + summary_start;
                self.index += 1;
                if offset + self.daf.summary_size * 8 > self.daf.mmap.len() {
                    return Some(Err(SpkError::InvalidData("Summary out of bounds".into())));
                }
                let summary_bytes = &self.daf.mmap[offset..offset + self.daf.summary_size * 8];
                return Some(self.parse_summary(summary_bytes));
            }
            if self.next_record != 0 {
                self.record_num = self.next_record;
            } else if self.count > 0 {
                return None;
            }
            if self.record_num == 0 {
                return None;
            }
            self.record_offset = (self.record_num - 1) * DAF_RECORD_SIZE;
            if self.record_offset + DAF_RECORD_SIZE > self.daf.mmap.len() {
                return Some(Err(SpkError::InvalidData(
                    "Summary record out of bounds".into(),
                )));
            }
            let record = &self.daf.mmap[self.record_offset..self.record_offset + DAF_RECORD_SIZE];
            self.next_record = self.daf.endian.read_f64(&record[0..8]) as usize;
            self.count = self.daf.endian.read_f64(&record[16..24]) as usize;
            self.index = 0;
            if self.count == 0 && self.next_record == 0 {
                return None;
            }
        }
    }
}

impl<'a> SummaryIterator<'a> {
    fn parse_summary(&self, bytes: &[u8]) -> Result<DafSummary, SpkError> {
        let mut doubles = Vec::with_capacity(self.daf.nd);
        for i in 0..self.daf.nd {
            doubles.push(self.daf.endian.read_f64(&bytes[i * 8..]));
        }
        let int_offset = self.daf.nd * 8;
        let mut ints = Vec::with_capacity(self.daf.ni);
        for i in 0..self.daf.ni {
            ints.push(self.daf.endian.read_i32(&bytes[int_offset + i * 4..]));
        }
        Ok(DafSummary { doubles, ints })
    }
}

#[derive(Debug, Clone)]
pub struct DafSummary {
    pub doubles: Vec<f64>,
    pub ints: Vec<i32>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn create_minimal_daf_header(nd: i32, ni: i32, fward: i32) -> Vec<u8> {
        let mut data = vec![0u8; DAF_RECORD_SIZE];

        // LOCIDW: "DAF/SPK " (8 bytes)
        data[0..8].copy_from_slice(b"DAF/SPK ");

        // ND (number of double components in summary) at offset 8
        data[8..12].copy_from_slice(&nd.to_le_bytes());

        // NI (number of integer components in summary) at offset 12
        data[12..16].copy_from_slice(&ni.to_le_bytes());

        // Internal file name at offset 16 (60 bytes)
        data[16..76].copy_from_slice(&[b' '; 60]);

        // FWARD (first summary record) at offset 76
        data[76..80].copy_from_slice(&fward.to_le_bytes());

        // BWARD (backward pointer) at offset 80
        data[80..84].copy_from_slice(&0i32.to_le_bytes());

        // FREE (first free address) at offset 84
        data[84..88].copy_from_slice(&0i32.to_le_bytes());

        // FTP string at offset 699-727
        if data.len() >= 727 {
            data[699..727].copy_from_slice(FTPSTR);
        }

        data
    }

    fn create_summary_record(next_record: f64, prev_record: f64, count: f64) -> Vec<u8> {
        let mut record = vec![0u8; DAF_RECORD_SIZE];
        record[0..8].copy_from_slice(&next_record.to_le_bytes());
        record[8..16].copy_from_slice(&prev_record.to_le_bytes());
        record[16..24].copy_from_slice(&count.to_le_bytes());
        record
    }

    #[test]
    fn test_endian_read_f64_little() {
        let val: f64 = 123.456789;
        let bytes = val.to_le_bytes();
        let result = Endian::Little.read_f64(&bytes);
        assert!((result - val).abs() < 1e-15);
    }

    #[test]
    fn test_endian_read_f64_big() {
        let val: f64 = 987.654321;
        let bytes = val.to_be_bytes();
        let result = Endian::Big.read_f64(&bytes);
        assert!((result - val).abs() < 1e-15);
    }

    #[test]
    fn test_endian_read_i32_little() {
        let val: i32 = 12345;
        let bytes = val.to_le_bytes();
        let result = Endian::Little.read_i32(&bytes);
        assert_eq!(result, val);
    }

    #[test]
    fn test_endian_read_i32_big() {
        let val: i32 = -54321;
        let bytes = val.to_be_bytes();
        let result = Endian::Big.read_i32(&bytes);
        assert_eq!(result, val);
    }

    #[test]
    fn test_daf_file_too_small() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("small.daf");
        std::fs::write(&file_path, b"small").unwrap();

        let result = DafFile::open(&file_path);
        assert!(result.is_err());
        match result {
            Err(SpkError::InvalidFormat(msg)) => assert!(msg.contains("too small")),
            _ => panic!("Expected InvalidFormat error"),
        }
    }

    #[test]
    fn test_daf_invalid_signature() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("bad.daf");

        let mut data = vec![0u8; DAF_RECORD_SIZE];
        data[0..8].copy_from_slice(b"NOTADAF!");
        std::fs::write(&file_path, &data).unwrap();

        let result = DafFile::open(&file_path);
        assert!(result.is_err());
        match result {
            Err(SpkError::InvalidFormat(msg)) => assert!(msg.contains("Invalid DAF signature")),
            _ => panic!("Expected InvalidFormat error"),
        }
    }

    #[test]
    fn test_daf_cannot_determine_endianness() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("badendian.daf");

        let mut data = vec![0u8; DAF_RECORD_SIZE];
        data[0..8].copy_from_slice(b"DAF/SPK ");
        // Write an invalid ND value (neither LE nor BE gives 1-100)
        data[8..12].copy_from_slice(&[0xFF, 0xFF, 0xFF, 0xFF]);
        // FTP string
        data[699..727].copy_from_slice(FTPSTR);

        std::fs::write(&file_path, &data).unwrap();

        let result = DafFile::open(&file_path);
        assert!(result.is_err());
        match result {
            Err(SpkError::InvalidFormat(msg)) => assert!(msg.contains("endianness")),
            _ => panic!("Expected InvalidFormat error"),
        }
    }

    #[test]
    fn test_daf_ftp_corruption() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("corrupt.daf");

        let mut data = create_minimal_daf_header(2, 6, 2);
        // Corrupt the FTP string
        data[699..727].copy_from_slice(&[0x00; 28]);

        std::fs::write(&file_path, &data).unwrap();

        let result = DafFile::open(&file_path);
        assert!(result.is_err());
        match result {
            Err(SpkError::InvalidFormat(msg)) => assert!(msg.contains("FTP corruption")),
            _ => panic!("Expected InvalidFormat error"),
        }
    }

    #[test]
    fn test_daf_open_valid() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("valid.daf");

        let data = create_minimal_daf_header(2, 6, 0);
        std::fs::write(&file_path, &data).unwrap();

        let result = DafFile::open(&file_path);
        assert!(result.is_ok());

        let daf = result.unwrap();
        assert_eq!(daf.nd, 2);
        assert_eq!(daf.ni, 6);
        assert_eq!(daf.endian, Endian::Little);
    }

    #[test]
    fn test_daf_open_nonexistent_file() {
        let result = DafFile::open("/nonexistent/path/file.bsp");
        assert!(result.is_err());
        match result {
            Err(SpkError::Io(_)) => {}
            _ => panic!("Expected Io error"),
        }
    }

    #[test]
    fn test_daf_read_array_out_of_bounds() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("valid.daf");

        let data = create_minimal_daf_header(2, 6, 0);
        std::fs::write(&file_path, &data).unwrap();

        let daf = DafFile::open(&file_path).unwrap();

        // Try to read beyond the file
        let result = daf.read_array(1, 1000);
        assert!(result.is_err());
        match result {
            Err(SpkError::InvalidData(msg)) => assert!(msg.contains("out of bounds")),
            _ => panic!("Expected InvalidData error"),
        }
    }

    #[test]
    fn test_daf_read_array_valid() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("valid.daf");

        let data = create_minimal_daf_header(2, 6, 0);
        std::fs::write(&file_path, &data).unwrap();

        let daf = DafFile::open(&file_path).unwrap();

        // Read first 10 doubles (80 bytes)
        let result = daf.read_array(1, 10);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 80);
    }

    #[test]
    fn test_daf_read_f64_array() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("valid.daf");

        let mut data = create_minimal_daf_header(2, 6, 0);

        // Add some known f64 values at known positions
        // Start after the header (byte 88) to not overwrite required fields
        let test_vals = [1.0f64, 2.0, 3.0, 4.0, 5.0];
        let start_offset = 88;
        for (i, &val) in test_vals.iter().enumerate() {
            let offset = start_offset + i * 8;
            data[offset..offset + 8].copy_from_slice(&val.to_le_bytes());
        }

        std::fs::write(&file_path, &data).unwrap();

        let daf = DafFile::open(&file_path).unwrap();

        // Read 5 doubles starting at position 12 (byte 88)
        let result = daf.read_f64_array(12, 5).unwrap();
        assert_eq!(result.len(), 5);
        for (i, &val) in test_vals.iter().enumerate() {
            assert!((result[i] - val).abs() < 1e-14);
        }
    }

    #[test]
    fn test_iter_summaries_empty() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("empty_summaries.daf");

        // fward = 0 means no summary records
        let data = create_minimal_daf_header(2, 6, 0);
        std::fs::write(&file_path, &data).unwrap();

        let daf = DafFile::open(&file_path).unwrap();
        let mut iter = daf.iter_summaries();

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_iter_summaries_with_record() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("with_summaries.daf");

        // Create header pointing to record 2
        let header = create_minimal_daf_header(2, 6, 2);

        // Create summary record with 1 summary
        let mut summary_record = create_summary_record(0.0, 0.0, 1.0);

        // Add a summary starting at offset 24
        // Summary has 2 doubles + 3 integers (6 ints / 2 = 3 packed)
        // Total summary size = 2 + 3 = 5 words
        let d1 = 100.0f64;
        let d2 = 200.0f64;
        summary_record[24..32].copy_from_slice(&d1.to_le_bytes());
        summary_record[32..40].copy_from_slice(&d2.to_le_bytes());

        // 6 integers packed
        let i1: i32 = 1;
        let i2: i32 = 2;
        let i3: i32 = 3;
        let i4: i32 = 4;
        let i5: i32 = 5;
        let i6: i32 = 6;
        summary_record[40..44].copy_from_slice(&i1.to_le_bytes());
        summary_record[44..48].copy_from_slice(&i2.to_le_bytes());
        summary_record[48..52].copy_from_slice(&i3.to_le_bytes());
        summary_record[52..56].copy_from_slice(&i4.to_le_bytes());
        summary_record[56..60].copy_from_slice(&i5.to_le_bytes());
        summary_record[60..64].copy_from_slice(&i6.to_le_bytes());

        // Combine header and summary record
        let mut data = header;
        data.extend(summary_record);

        std::fs::write(&file_path, &data).unwrap();

        let daf = DafFile::open(&file_path).unwrap();
        let summaries: Vec<_> = daf.iter_summaries().collect();

        assert_eq!(summaries.len(), 1);
        let summary = summaries[0].as_ref().unwrap();
        assert_eq!(summary.doubles.len(), 2);
        assert!((summary.doubles[0] - 100.0).abs() < 1e-14);
        assert!((summary.doubles[1] - 200.0).abs() < 1e-14);
        assert_eq!(summary.ints.len(), 6);
        assert_eq!(summary.ints[0], 1);
        assert_eq!(summary.ints[5], 6);
    }

    #[test]
    fn test_endian_equality() {
        assert_eq!(Endian::Little, Endian::Little);
        assert_eq!(Endian::Big, Endian::Big);
        assert_ne!(Endian::Little, Endian::Big);
    }

    #[test]
    fn test_daf_summary_struct() {
        let summary = DafSummary {
            doubles: vec![1.0, 2.0],
            ints: vec![10, 20, 30],
        };

        assert_eq!(summary.doubles.len(), 2);
        assert_eq!(summary.ints.len(), 3);

        // Test Clone
        let cloned = summary.clone();
        assert_eq!(cloned.doubles, summary.doubles);
        assert_eq!(cloned.ints, summary.ints);
    }

    #[test]
    fn test_summary_iterator_record_out_of_bounds() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("bad_record.daf");

        // Create header pointing to a record that doesn't exist
        let data = create_minimal_daf_header(2, 6, 10);
        std::fs::write(&file_path, &data).unwrap();

        let daf = DafFile::open(&file_path).unwrap();
        let mut iter = daf.iter_summaries();

        let result = iter.next();
        assert!(result.is_some());
        assert!(result.unwrap().is_err());
    }

    #[test]
    fn test_big_endian_detection() {
        // Line 82: Big-endian DAF detection
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("big_endian.daf");

        let mut data = vec![0u8; DAF_RECORD_SIZE];

        // DAF signature
        data[0..8].copy_from_slice(b"DAF/SPK ");

        // ND = 2 in big-endian (value that's in 1..=100 when read as BE but not LE)
        // Write 2 as big-endian
        let nd_be: i32 = 2;
        data[8..12].copy_from_slice(&nd_be.to_be_bytes());

        // NI = 6 in big-endian
        let ni_be: i32 = 6;
        data[12..16].copy_from_slice(&ni_be.to_be_bytes());

        // fward = 0 in big-endian
        let fward_be: i32 = 0;
        data[76..80].copy_from_slice(&fward_be.to_be_bytes());

        // FTP string
        data[699..727].copy_from_slice(FTPSTR);

        std::fs::write(&file_path, &data).unwrap();

        let result = DafFile::open(&file_path);
        assert!(result.is_ok());
        let daf = result.unwrap();
        assert_eq!(daf.endian, Endian::Big);
    }

    #[test]
    fn test_summary_iterator_with_multiple_records() {
        // This tests line 152: when next_record != 0, we follow to the next record
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("multi_record.daf");

        // Create header pointing to record 2
        let header = create_minimal_daf_header(2, 6, 2);

        // Create first summary record that points to record 3
        let mut record2 = create_summary_record(3.0, 0.0, 1.0); // next_record=3, count=1

        // Add a summary at offset 24
        let d1 = 100.0f64;
        let d2 = 200.0f64;
        record2[24..32].copy_from_slice(&d1.to_le_bytes());
        record2[32..40].copy_from_slice(&d2.to_le_bytes());
        // 6 integers
        for i in 0..6 {
            let val: i32 = (i + 1) as i32;
            record2[40 + i * 4..44 + i * 4].copy_from_slice(&val.to_le_bytes());
        }

        // Create second summary record (record 3) with no next
        let mut record3 = create_summary_record(0.0, 0.0, 1.0); // next_record=0, count=1

        // Add a summary at offset 24
        let d3 = 300.0f64;
        let d4 = 400.0f64;
        record3[24..32].copy_from_slice(&d3.to_le_bytes());
        record3[32..40].copy_from_slice(&d4.to_le_bytes());
        for i in 0..6 {
            let val: i32 = (i + 10) as i32;
            record3[40 + i * 4..44 + i * 4].copy_from_slice(&val.to_le_bytes());
        }

        // Combine all records
        let mut data = header;
        data.extend(record2);
        data.extend(record3);

        std::fs::write(&file_path, &data).unwrap();

        let daf = DafFile::open(&file_path).unwrap();
        let summaries: Vec<_> = daf.iter_summaries().collect();

        assert_eq!(summaries.len(), 2);
        let s1 = summaries[0].as_ref().unwrap();
        let s2 = summaries[1].as_ref().unwrap();
        assert!((s1.doubles[0] - 100.0).abs() < 1e-10);
        assert!((s2.doubles[0] - 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_summary_iterator_empty_record_with_next() {
        // This tests line 170: count == 0 && next_record == 0 after reading record
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("empty_record.daf");

        // Create header pointing to record 2
        let header = create_minimal_daf_header(2, 6, 2);

        // Create summary record with 0 summaries and no next record
        let record = create_summary_record(0.0, 0.0, 0.0); // next_record=0, count=0

        let mut data = header;
        data.extend(record);

        std::fs::write(&file_path, &data).unwrap();

        let daf = DafFile::open(&file_path).unwrap();
        let summaries: Vec<_> = daf.iter_summaries().collect();

        assert!(summaries.is_empty());
    }

    #[test]
    fn test_summary_out_of_bounds() {
        // This tests line 146: summary extends beyond file length
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("truncated.daf");

        // Create header pointing to record 2
        let header = create_minimal_daf_header(2, 6, 2);

        // Create summary record claiming 100 summaries (but we won't have space)
        let record = create_summary_record(0.0, 0.0, 100.0); // count=100

        let mut data = header;
        data.extend(record);
        // Don't add enough data for all summaries

        std::fs::write(&file_path, &data).unwrap();

        let daf = DafFile::open(&file_path).unwrap();
        let mut iter = daf.iter_summaries();

        // First summary should succeed (offset 24 of record)
        let first = iter.next();
        assert!(first.is_some());
        let first_result = first.unwrap();
        assert!(first_result.is_ok());

        // Eventually we'll hit out of bounds
        let mut found_error = false;
        for result in iter {
            if result.is_err() {
                found_error = true;
                match result {
                    Err(SpkError::InvalidData(msg)) => assert!(msg.contains("out of bounds")),
                    _ => panic!("Expected InvalidData error"),
                }
                break;
            }
        }
        assert!(found_error, "Should have found an out of bounds error");
    }
}
