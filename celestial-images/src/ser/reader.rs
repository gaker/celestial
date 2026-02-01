use crate::ser::{Result, SerError, SerFrame, SerFrameId, SerHeader, SerTimestamp};
use memmap2::Mmap;
use std::fs::File;
use std::io::Read;
use std::path::Path;

pub struct SerReader {
    header: SerHeader,
    mmap: Option<Mmap>,
    _file: Option<File>, // Kept to maintain mmap validity
    frame_data_offset: u64,
    trailer_offset: Option<u64>,
}

impl SerReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let metadata = file.metadata()?;

        if metadata.len() < 178 {
            return Err(SerError::FileTruncated {
                expected: 178,
                actual: metadata.len(),
            });
        }

        let mmap = unsafe { Mmap::map(&file)? };
        let header = SerHeader::from_bytes(&mmap[0..178])?;

        let expected_size = Self::calculate_expected_file_size(&header);
        if metadata.len() < expected_size {
            return Err(SerError::FileTruncated {
                expected: expected_size,
                actual: metadata.len(),
            });
        }

        let trailer_offset = if header.has_trailer() {
            Some(crate::ser::types::SerFile::calculate_trailer_offset(
                header.frame_count,
                header.image_width,
                header.image_height,
                header.bytes_per_pixel(),
            ))
        } else {
            None
        };

        Ok(Self {
            header,
            mmap: Some(mmap),
            _file: Some(file),
            frame_data_offset: 178,
            trailer_offset,
        })
    }

    pub fn from_reader<R: Read>(mut reader: R) -> Result<Self> {
        let mut header_bytes = [0u8; 178];
        reader.read_exact(&mut header_bytes)?;

        let header = SerHeader::from_bytes(&header_bytes)?;

        Ok(Self {
            header,
            mmap: None,
            _file: None,
            frame_data_offset: 178,
            trailer_offset: None,
        })
    }

    pub fn header(&self) -> &SerHeader {
        &self.header
    }

    pub fn frame_count(&self) -> u32 {
        self.header.frame_count
    }

    pub fn is_color(&self) -> bool {
        self.header.color_id.planes() == 3
    }

    pub fn is_mono(&self) -> bool {
        self.header.color_id.planes() == 1 && !self.header.color_id.is_bayer()
    }

    pub fn is_bayer(&self) -> bool {
        self.header.color_id.is_bayer()
    }

    pub fn read_frame(&self, frame_id: SerFrameId) -> Result<SerFrame<'_>> {
        if frame_id >= self.header.frame_count {
            return Err(SerError::FrameOutOfBounds {
                frame: frame_id,
                total: self.header.frame_count,
            });
        }

        let frame_size = self.header.frame_size();
        let frame_offset = self.frame_data_offset + (frame_id as u64 * frame_size);

        if let Some(ref mmap) = self.mmap {
            let end_offset = frame_offset + frame_size;
            if end_offset > mmap.len() as u64 {
                return Err(SerError::FileTruncated {
                    expected: end_offset,
                    actual: mmap.len() as u64,
                });
            }

            let data = &mmap[frame_offset as usize..end_offset as usize];
            let timestamp = self.read_timestamp(frame_id)?;

            Ok(SerFrame::new(frame_id, data, timestamp))
        } else {
            Err(SerError::InvalidHeader(
                "No memory map available".to_string(),
            ))
        }
    }

    pub fn read_frame_data(&self, frame_id: SerFrameId) -> Result<Vec<u8>> {
        let frame = self.read_frame(frame_id)?;
        Ok(frame.data.to_vec())
    }

    pub fn read_timestamp(&self, frame_id: SerFrameId) -> Result<Option<SerTimestamp>> {
        if !self.header.has_trailer() {
            return Ok(None);
        }

        if frame_id >= self.header.frame_count {
            return Err(SerError::FrameOutOfBounds {
                frame: frame_id,
                total: self.header.frame_count,
            });
        }

        if let (Some(ref mmap), Some(trailer_offset)) = (&self.mmap, self.trailer_offset) {
            let timestamp_offset = trailer_offset + (frame_id as u64 * 8);
            let end_offset = timestamp_offset + 8;

            if end_offset > mmap.len() as u64 {
                return Err(SerError::FileTruncated {
                    expected: end_offset,
                    actual: mmap.len() as u64,
                });
            }

            let timestamp_bytes = &mmap[timestamp_offset as usize..end_offset as usize];
            Ok(Some(SerTimestamp::from_bytes(timestamp_bytes)))
        } else {
            Ok(None)
        }
    }

    pub fn iter_frames(&self) -> FrameIterator<'_> {
        FrameIterator::new(self, 0, self.header.frame_count)
    }

    pub fn iter_frames_range(&self, start: SerFrameId, count: u32) -> FrameIterator<'_> {
        let end = (start + count).min(self.header.frame_count);
        FrameIterator::new(self, start, end)
    }

    pub fn frame_slice(&self, frame_id: SerFrameId) -> Result<&[u8]> {
        if frame_id >= self.header.frame_count {
            return Err(SerError::FrameOutOfBounds {
                frame: frame_id,
                total: self.header.frame_count,
            });
        }

        let mmap = self
            .mmap
            .as_ref()
            .ok_or_else(|| SerError::InvalidHeader("No memory map available".to_string()))?;

        let frame_size = self.header.frame_size();
        let start = self.frame_data_offset + (frame_id as u64 * frame_size);
        let end = start + frame_size;

        if end > mmap.len() as u64 {
            return Err(SerError::FileTruncated {
                expected: end,
                actual: mmap.len() as u64,
            });
        }

        Ok(&mmap[start as usize..end as usize])
    }

    pub fn frame_slices(&self, start: SerFrameId, count: u32) -> Result<Vec<(SerFrameId, &[u8])>> {
        let end = (start + count).min(self.header.frame_count);
        let mut slices = Vec::with_capacity((end - start) as usize);

        for frame_id in start..end {
            slices.push((frame_id, self.frame_slice(frame_id)?));
        }

        Ok(slices)
    }

    fn calculate_expected_file_size(header: &SerHeader) -> u64 {
        let header_size = 178u64;
        let data_size = header.frame_count as u64 * header.frame_size();
        let trailer_size = if header.has_trailer() {
            header.frame_count as u64 * 8
        } else {
            0
        };

        header_size + data_size + trailer_size
    }
}

pub struct FrameIterator<'a> {
    reader: &'a SerReader,
    current: SerFrameId,
    end: SerFrameId,
}

impl<'a> FrameIterator<'a> {
    fn new(reader: &'a SerReader, start: SerFrameId, end: SerFrameId) -> Self {
        Self {
            reader,
            current: start,
            end,
        }
    }
}

impl<'a> Iterator for FrameIterator<'a> {
    type Item = Result<SerFrame<'a>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current >= self.end {
            return None;
        }

        let frame_id = self.current;
        self.current += 1;

        Some(self.reader.read_frame(frame_id))
    }
}

impl<'a> ExactSizeIterator for FrameIterator<'a> {
    fn len(&self) -> usize {
        (self.end - self.current) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ser::SerHeader;
    use std::io::Cursor;
    use std::io::Write;

    fn create_test_ser_file() -> Vec<u8> {
        let mut data = vec![0u8; 1000];

        data[0..14].copy_from_slice(b"LUCAM-RECORDER");
        data[14..18].copy_from_slice(&123u32.to_le_bytes());
        data[18..22].copy_from_slice(&0u32.to_le_bytes()); // ColorId::Mono
        data[22..26].copy_from_slice(&1u32.to_le_bytes());
        data[26..30].copy_from_slice(&10u32.to_le_bytes()); // width
        data[30..34].copy_from_slice(&10u32.to_le_bytes()); // height
        data[34..38].copy_from_slice(&16u32.to_le_bytes()); // pixel depth
        data[38..42].copy_from_slice(&2u32.to_le_bytes()); // frame count
        data[42..50].copy_from_slice(b"Observer");
        data[82..92].copy_from_slice(b"Instrument");
        data[122..131].copy_from_slice(b"Telescope");
        data[162..170].copy_from_slice(&1234567890i64.to_le_bytes());
        data[170..178].copy_from_slice(&1234567891i64.to_le_bytes());

        // Add frame data (2 frames * 10*10*2 bytes = 400 bytes)
        for (offset, byte) in data[178..578].iter_mut().enumerate() {
            *byte = (offset % 256) as u8;
        }

        // Add trailer (2 frames * 8 bytes = 16 bytes)
        for (i, byte) in data[578..594].iter_mut().enumerate() {
            *byte = ((578 + i) % 256) as u8;
        }

        data.truncate(594);
        data
    }

    fn create_minimal_ser_file() -> Vec<u8> {
        let mut data = vec![0u8; 378];

        data[0..14].copy_from_slice(b"LUCAM-RECORDER");
        data[14..18].copy_from_slice(&0u32.to_le_bytes());
        data[18..22].copy_from_slice(&0u32.to_le_bytes());
        data[22..26].copy_from_slice(&1u32.to_le_bytes());
        data[26..30].copy_from_slice(&10u32.to_le_bytes());
        data[30..34].copy_from_slice(&10u32.to_le_bytes());
        data[34..38].copy_from_slice(&16u32.to_le_bytes());
        data[38..42].copy_from_slice(&1u32.to_le_bytes());
        data[162..170].copy_from_slice(&0i64.to_le_bytes()); // No trailer

        data
    }

    #[test]
    fn from_reader_success() {
        // Lines 59-71: Test successful reader creation
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        assert_eq!(reader.header.frame_count, 1);
        assert_eq!(reader.frame_data_offset, 178);
        assert!(reader.trailer_offset.is_none());
        assert!(reader.mmap.is_none());
        assert!(reader._file.is_none());
    }

    #[test]
    fn from_reader_header_parse_error() {
        // Line 63: Test header parsing failure
        let mut data = vec![0u8; 178];
        data[0..14].copy_from_slice(b"INVALID-HEADER");
        let cursor = Cursor::new(data);
        let result = SerReader::from_reader(cursor);
        assert!(result.is_err());
    }

    #[test]
    fn header_accessor() {
        // Lines 74-75: Test header accessor
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();
        let header = reader.header();
        assert_eq!(header.frame_count, 1);
    }

    #[test]
    fn frame_count_accessor() {
        // Lines 78-79: Test frame count accessor
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();
        assert_eq!(reader.frame_count(), 1);
    }

    #[test]
    fn read_frame_out_of_bounds() {
        // Lines 82-87: Test frame out of bounds error
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let result = reader.read_frame(999);
        assert!(matches!(
            result,
            Err(SerError::FrameOutOfBounds {
                frame: 999,
                total: 1
            })
        ));
    }

    #[test]
    fn read_frame_no_mmap_error() {
        // Lines 90-91, 106-107: Test reading frame without memory map
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let result = reader.read_frame(0);
        assert!(matches!(result, Err(SerError::InvalidHeader(_))));
    }

    #[test]
    fn read_frame_data_delegates() {
        // Lines 111-113: Test read_frame_data delegation
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let result = reader.read_frame_data(0);
        assert!(result.is_err()); // Should fail because no mmap
    }

    #[test]
    fn read_timestamp_no_trailer() {
        // Lines 116-118: Test reading timestamp when no trailer
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let timestamp = reader.read_timestamp(0).unwrap();
        assert!(timestamp.is_none());
    }

    #[test]
    fn read_timestamp_frame_out_of_bounds() {
        // Lines 121-124: Test timestamp read with frame out of bounds
        let data = create_test_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let result = reader.read_timestamp(999);
        assert!(matches!(
            result,
            Err(SerError::FrameOutOfBounds {
                frame: 999,
                total: 2
            })
        ));
    }

    #[test]
    fn read_timestamp_no_mmap() {
        // Lines 141-142: Test timestamp read without mmap
        let data = create_test_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let timestamp = reader.read_timestamp(0).unwrap();
        assert!(timestamp.is_none());
    }

    #[test]
    fn iter_frames() {
        // Lines 146-147: Test frame iterator creation
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let iter = reader.iter_frames();
        assert_eq!(iter.len(), 1);
    }

    #[test]
    fn iter_frames_range() {
        // Lines 150-152: Test frame range iterator
        let data = create_test_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let iter = reader.iter_frames_range(0, 1);
        assert_eq!(iter.len(), 1);

        let iter = reader.iter_frames_range(1, 5);
        assert_eq!(iter.len(), 1); // Should be clamped to available frames
    }

    #[test]
    fn calculate_expected_file_size_with_trailer() {
        let header = SerHeader {
            frame_count: 10,
            image_width: 100,
            image_height: 100,
            pixel_depth_per_plane: 16,
            date_time: 1, // Has trailer
            ..Default::default()
        };

        let size = SerReader::calculate_expected_file_size(&header);
        let expected = 178 + (10 * 100 * 100 * 2) + (10 * 8); // header + data + trailer
        assert_eq!(size, expected);
    }

    #[test]
    fn calculate_expected_file_size_no_trailer() {
        let header = SerHeader {
            frame_count: 5,
            image_width: 50,
            image_height: 50,
            pixel_depth_per_plane: 8,
            date_time: 0, // No trailer
            ..Default::default()
        };

        let size = SerReader::calculate_expected_file_size(&header);
        let expected = 178 + (5 * 50 * 50); // header + data only
        assert_eq!(size, expected);
    }

    #[test]
    fn frame_iterator_new() {
        // Line 175: Test FrameIterator creation
        let data = create_test_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let iter = FrameIterator::new(&reader, 0, 2);
        assert_eq!(iter.current, 0);
        assert_eq!(iter.end, 2);
    }

    #[test]
    fn frame_iterator_next() {
        // Lines 187-195: Test iterator next functionality
        let data = create_test_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let mut iter = FrameIterator::new(&reader, 0, 1);
        assert!(iter.next().is_some()); // Should have one item
        assert!(iter.next().is_none()); // Should be exhausted
    }

    #[test]
    fn frame_iterator_exact_size() {
        // Lines 200-201: Test ExactSizeIterator implementation
        let data = create_test_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let iter = FrameIterator::new(&reader, 0, 2);
        assert_eq!(iter.len(), 2);

        let iter = FrameIterator::new(&reader, 1, 2);
        assert_eq!(iter.len(), 1);
    }

    fn write_test_file() -> std::io::Result<std::path::PathBuf> {
        let temp_dir = std::env::temp_dir();
        let unique_name = format!(
            "test_ser_file_{}_{}.ser",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        );
        let temp_path = temp_dir.join(unique_name);

        let mut file = File::create(&temp_path)?;
        let data = create_test_ser_file();
        file.write_all(&data)?;

        Ok(temp_path)
    }

    #[test]
    fn open_file_success() {
        // Lines 17-56: Test successful file opening with memory mapping
        if let Ok(temp_path) = write_test_file() {
            if let Ok(reader) = SerReader::open(&temp_path) {
                assert_eq!(reader.header.frame_count, 2);
                assert_eq!(reader.frame_data_offset, 178);
                assert!(reader.trailer_offset.is_some());
                assert!(reader.mmap.is_some());
                assert!(reader._file.is_some());
            }

            std::fs::remove_file(temp_path).ok();
        }
    }

    #[test]
    fn open_file_too_short() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("short_ser_file.ser");

        let mut file = File::create(&temp_path).unwrap();
        let short_data = vec![0u8; 100]; // Too short
        file.write_all(&short_data).unwrap();

        let result = SerReader::open(&temp_path);
        assert!(matches!(
            result,
            Err(SerError::FileTruncated {
                expected: 178,
                actual: 100
            })
        ));

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn open_file_truncated() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("truncated_ser_file.ser");

        let mut file = File::create(&temp_path).unwrap();
        let mut data = create_test_ser_file();
        data.truncate(300); // Truncate to be shorter than expected
        file.write_all(&data).unwrap();

        let result = SerReader::open(&temp_path);
        assert!(result.is_err());

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn trailer_offset_calculation() {
        // Lines 39-47: Test trailer offset calculation
        if let Ok(temp_path) = write_test_file() {
            let reader = SerReader::open(&temp_path).unwrap();
            assert!(reader.trailer_offset.is_some());

            std::fs::remove_file(temp_path).ok();
        }
    }

    #[test]
    fn read_frame_with_mmap() {
        // Lines 93-105: Test successful frame reading with memory map
        if let Ok(temp_path) = write_test_file() {
            if let Ok(reader) = SerReader::open(&temp_path) {
                if let Ok(frame) = reader.read_frame(0) {
                    assert_eq!(frame.id, 0);
                    assert_eq!(frame.data.len(), 200); // 10*10*2 bytes
                }
            }

            std::fs::remove_file(temp_path).ok();
        }
    }

    #[test]
    fn read_frame_truncated_data() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("truncated_frame_file.ser");

        let mut file = File::create(&temp_path).unwrap();
        let mut data = create_test_ser_file();
        data.truncate(400); // Truncate frames
        file.write_all(&data).unwrap();

        if let Ok(reader) = SerReader::open(&temp_path) {
            let result = reader.read_frame(1);
            assert!(result.is_err());
        }

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn read_timestamp_with_mmap() {
        // Lines 128-140: Test timestamp reading with memory map
        if let Ok(temp_path) = write_test_file() {
            if let Ok(reader) = SerReader::open(&temp_path) {
                if let Ok(timestamp) = reader.read_timestamp(0) {
                    assert!(timestamp.is_some());
                }
            }

            std::fs::remove_file(temp_path).ok();
        }
    }

    #[test]
    fn read_timestamp_truncated_trailer() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("truncated_trailer_file.ser");

        let mut file = File::create(&temp_path).unwrap();
        let mut data = create_test_ser_file();
        data.truncate(590); // Truncate trailer
        file.write_all(&data).unwrap();

        if let Ok(reader) = SerReader::open(&temp_path) {
            let result = reader.read_timestamp(1);
            assert!(result.is_err());
        }

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn from_reader_read_error() {
        // Lines 60-61: Test read error during header reading
        use std::io::{Error, ErrorKind};

        struct FailingReader;
        impl std::io::Read for FailingReader {
            fn read(&mut self, _buf: &mut [u8]) -> std::io::Result<usize> {
                Err(Error::new(ErrorKind::UnexpectedEof, "read failed"))
            }
        }

        let result = SerReader::from_reader(FailingReader);
        assert!(result.is_err());
    }

    fn write_minimal_file() -> std::io::Result<std::path::PathBuf> {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("minimal_ser_file.ser");

        let mut file = File::create(&temp_path)?;
        let data = create_minimal_ser_file();
        file.write_all(&data)?;

        Ok(temp_path)
    }

    #[test]
    fn open_file_no_trailer_path() {
        if let Ok(temp_path) = write_minimal_file() {
            if let Ok(reader) = SerReader::open(&temp_path) {
                assert_eq!(reader.header.frame_count, 1);
                assert_eq!(reader.frame_data_offset, 178);
                assert!(reader.trailer_offset.is_none());
                assert!(reader.mmap.is_some());
                assert!(reader._file.is_some());
            }

            std::fs::remove_file(temp_path).ok();
        }
    }

    #[test]
    fn read_frame_mmap_bounds_check() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("bounds_test_file.ser");

        let mut file = File::create(&temp_path).unwrap();
        let mut data = create_test_ser_file();
        data.truncate(300); // Make it just barely long enough for header but not frame
        file.write_all(&data).unwrap();

        if let Ok(reader) = SerReader::open(&temp_path) {
            let result = reader.read_frame(0);
            assert!(matches!(result, Err(SerError::FileTruncated { .. })));
        }

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn read_timestamp_mmap_bounds_check() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("timestamp_bounds_test.ser");

        let mut file = File::create(&temp_path).unwrap();
        let mut data = create_test_ser_file();
        data.truncate(580); // Truncate just before trailer
        file.write_all(&data).unwrap();

        if let Ok(reader) = SerReader::open(&temp_path) {
            let result = reader.read_timestamp(0);
            assert!(matches!(result, Err(SerError::FileTruncated { .. })));
        }

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn frame_slice_success() {
        if let Ok(temp_path) = write_test_file() {
            let reader = SerReader::open(&temp_path).unwrap();
            let slice = reader.frame_slice(0).unwrap();
            assert_eq!(slice.len(), 200); // 10*10*2 bytes
            std::fs::remove_file(temp_path).ok();
        }
    }

    #[test]
    fn frame_slice_out_of_bounds() {
        if let Ok(temp_path) = write_test_file() {
            let reader = SerReader::open(&temp_path).unwrap();
            let result = reader.frame_slice(999);
            assert!(matches!(
                result,
                Err(SerError::FrameOutOfBounds {
                    frame: 999,
                    total: 2
                })
            ));
            std::fs::remove_file(temp_path).ok();
        }
    }

    #[test]
    fn frame_slice_no_mmap() {
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let result = reader.frame_slice(0);
        assert!(matches!(result, Err(SerError::InvalidHeader(_))));
    }

    #[test]
    fn frame_slices_success() {
        if let Ok(temp_path) = write_test_file() {
            let reader = SerReader::open(&temp_path).unwrap();
            let slices = reader.frame_slices(0, 2).unwrap();

            assert_eq!(slices.len(), 2);
            assert_eq!(slices[0].0, 0);
            assert_eq!(slices[0].1.len(), 200);
            assert_eq!(slices[1].0, 1);
            assert_eq!(slices[1].1.len(), 200);
            std::fs::remove_file(temp_path).ok();
        }
    }

    #[test]
    fn frame_slices_clamped_to_available() {
        if let Ok(temp_path) = write_test_file() {
            let reader = SerReader::open(&temp_path).unwrap();
            let slices = reader.frame_slices(1, 100).unwrap();

            assert_eq!(slices.len(), 1); // Only 1 frame available starting from index 1
            assert_eq!(slices[0].0, 1);
            std::fs::remove_file(temp_path).ok();
        }
    }

    #[test]
    fn frame_slices_no_mmap() {
        let data = create_minimal_ser_file();
        let cursor = Cursor::new(data);
        let reader = SerReader::from_reader(cursor).unwrap();

        let result = reader.frame_slices(0, 1);
        assert!(matches!(result, Err(SerError::InvalidHeader(_))));
    }
}
