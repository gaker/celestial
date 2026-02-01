use crate::ser::{FrameBuffer, Result, SerError, SerFrameId, SerHeader, SerTimestamp};
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

pub struct SerWriter {
    writer: BufWriter<File>,
    header: SerHeader,
    buffer: Option<FrameBuffer>,
    frames_written: u32,
    timestamps: Vec<SerTimestamp>,
}

impl SerWriter {
    pub fn create<P: AsRef<Path>>(path: P, header: SerHeader) -> Result<Self> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        writer.write_all(&header.to_bytes())?;
        writer.flush()?;

        Ok(Self {
            writer,
            header,
            buffer: None,
            frames_written: 0,
            timestamps: Vec::new(),
        })
    }

    pub fn with_buffer<P: AsRef<Path>>(
        path: P,
        header: SerHeader,
        buffer_mb: usize,
    ) -> Result<Self> {
        let mut writer = Self::create(path, header.clone())?;
        writer.buffer = Some(FrameBuffer::with_capacity(&header, buffer_mb));
        Ok(writer)
    }

    pub fn header(&self) -> &SerHeader {
        &self.header
    }

    pub fn write_frame(&mut self, data: &[u8]) -> Result<SerFrameId> {
        self.write_frame_with_timestamp(data, None)
    }

    pub fn write_frame_with_timestamp(
        &mut self,
        data: &[u8],
        timestamp: Option<SerTimestamp>,
    ) -> Result<SerFrameId> {
        let expected_size = self.header.frame_size() as usize;
        if data.len() != expected_size {
            return Err(SerError::BufferSizeMismatch {
                expected: expected_size,
                actual: data.len(),
            });
        }

        if let Some(ref mut buffer) = self.buffer {
            let frame_id = buffer.push_frame(data.to_vec())?;
            if let Some(ts) = timestamp {
                self.timestamps.push(ts);
            }
            self.frames_written += 1;
            Ok(frame_id)
        } else {
            self.writer.write_all(data)?;
            if let Some(ts) = timestamp {
                self.timestamps.push(ts);
            }
            let frame_id = self.frames_written;
            self.frames_written += 1;
            Ok(frame_id)
        }
    }

    pub fn flush_buffer(&mut self) -> Result<()> {
        if let Some(ref mut buffer) = self.buffer {
            for i in 0..buffer.available_frames() {
                if let Some(frame_data) = buffer.get_frame(i as u32) {
                    self.writer.write_all(&frame_data)?;
                }
            }
            buffer.clear();
        }
        self.writer.flush()?;
        Ok(())
    }

    pub fn finalize(mut self) -> Result<()> {
        self.flush_buffer()?;

        if !self.timestamps.is_empty() && self.header.has_trailer() {
            for timestamp in &self.timestamps {
                let mut buf = [0u8; 8];
                timestamp.to_bytes(&mut buf);
                self.writer.write_all(&buf)?;
            }
        }

        let current_pos = self.writer.stream_position()?;
        self.writer.seek(SeekFrom::Start(38))?;
        let frame_count_bytes = self.frames_written.to_le_bytes();
        self.writer.write_all(&frame_count_bytes)?;
        self.writer.seek(SeekFrom::Start(current_pos))?;

        self.writer.flush()?;
        Ok(())
    }

    pub fn frames_written(&self) -> u32 {
        self.frames_written
    }

    pub fn is_buffered(&self) -> bool {
        self.buffer.is_some()
    }

    pub fn buffer_usage(&self) -> Option<(usize, usize)> {
        self.buffer
            .as_ref()
            .map(|b| (b.available_frames(), b.available_frames()))
    }
}

impl Drop for SerWriter {
    fn drop(&mut self) {
        let _ = self.flush_buffer();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn create_test_header() -> SerHeader {
        SerHeader {
            image_width: 10,
            image_height: 10,
            pixel_depth_per_plane: 16,
            frame_count: 0,
            date_time: 1234567890,
            ..Default::default()
        }
    }

    fn create_test_header_no_trailer() -> SerHeader {
        SerHeader {
            image_width: 10,
            image_height: 10,
            pixel_depth_per_plane: 16,
            frame_count: 0,
            date_time: 0,
            ..Default::default()
        }
    }

    fn create_test_frame_data() -> Vec<u8> {
        vec![42u8; 200]
    }

    #[test]
    fn create_writer_success() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("test_writer.ser");

        let header = create_test_header();
        let writer = SerWriter::create(&temp_path, header.clone()).unwrap();

        assert_eq!(writer.header.image_width, 10);
        assert_eq!(writer.frames_written, 0);
        assert!(writer.buffer.is_none());
        assert!(writer.timestamps.is_empty());

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn create_writer_file_error() {
        let invalid_path = "/invalid/path/that/does/not/exist/test.ser";
        let header = create_test_header();
        let result = SerWriter::create(invalid_path, header);
        assert!(result.is_err());
    }

    #[test]
    fn with_buffer_success() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("test_buffered_writer.ser");

        let header = create_test_header();
        let writer = SerWriter::with_buffer(&temp_path, header, 1).unwrap();

        assert!(writer.buffer.is_some());
        assert!(writer.is_buffered());

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn header_accessor() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("header_test.ser");

        let header = create_test_header();
        let writer = SerWriter::create(&temp_path, header.clone()).unwrap();

        let retrieved_header = writer.header();
        assert_eq!(retrieved_header.image_width, 10);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn write_frame_delegates() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("frame_delegate_test.ser");

        let header = create_test_header_no_trailer();
        let mut writer = SerWriter::create(&temp_path, header).unwrap();

        let frame_data = create_test_frame_data();
        let frame_id = writer.write_frame(&frame_data).unwrap();
        assert_eq!(frame_id, 0);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn write_frame_wrong_size() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("wrong_size_test.ser");

        let header = create_test_header();
        let mut writer = SerWriter::create(&temp_path, header).unwrap();

        let wrong_size_data = vec![0u8; 100];
        let result = writer.write_frame_with_timestamp(&wrong_size_data, None);
        assert!(matches!(
            result,
            Err(SerError::BufferSizeMismatch {
                expected: 200,
                actual: 100
            })
        ));

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn write_frame_with_buffer() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("buffered_frame_test.ser");

        let header = create_test_header();
        let mut writer = SerWriter::with_buffer(&temp_path, header, 1).unwrap();

        let frame_data = create_test_frame_data();
        let timestamp = SerTimestamp::new(123456);
        let frame_id = writer
            .write_frame_with_timestamp(&frame_data, Some(timestamp))
            .unwrap();

        assert_eq!(frame_id, 0);
        assert_eq!(writer.frames_written(), 1);
        assert_eq!(writer.timestamps.len(), 1);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn write_frame_no_buffer() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("direct_frame_test.ser");

        let header = create_test_header_no_trailer();
        let mut writer = SerWriter::create(&temp_path, header).unwrap();

        let frame_data = create_test_frame_data();
        let timestamp = SerTimestamp::new(789012);
        let frame_id = writer
            .write_frame_with_timestamp(&frame_data, Some(timestamp))
            .unwrap();

        assert_eq!(frame_id, 0);
        assert_eq!(writer.frames_written(), 1);
        assert_eq!(writer.timestamps.len(), 1);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn flush_buffer_with_buffer() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("flush_buffer_test.ser");

        let header = create_test_header();
        let mut writer = SerWriter::with_buffer(&temp_path, header, 1).unwrap();

        let frame_data = create_test_frame_data();
        writer.write_frame(&frame_data).unwrap();

        writer.flush_buffer().unwrap();

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn flush_buffer_no_buffer() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("flush_no_buffer_test.ser");

        let header = create_test_header();
        let mut writer = SerWriter::create(&temp_path, header).unwrap();

        writer.flush_buffer().unwrap();

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn finalize_with_timestamps() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("finalize_timestamps_test.ser");

        let header = create_test_header();
        let mut writer = SerWriter::create(&temp_path, header).unwrap();

        let frame_data = create_test_frame_data();
        let timestamp = SerTimestamp::new(555555);
        writer
            .write_frame_with_timestamp(&frame_data, Some(timestamp))
            .unwrap();

        writer.finalize().unwrap();

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn finalize_no_timestamps() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("finalize_no_timestamps_test.ser");

        let header = create_test_header_no_trailer();
        let mut writer = SerWriter::create(&temp_path, header).unwrap();

        let frame_data = create_test_frame_data();
        writer.write_frame(&frame_data).unwrap();

        writer.finalize().unwrap();

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn frames_written_accessor() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("frames_count_test.ser");

        let header = create_test_header();
        let mut writer = SerWriter::create(&temp_path, header).unwrap();

        assert_eq!(writer.frames_written(), 0);

        let frame_data = create_test_frame_data();
        writer.write_frame(&frame_data).unwrap();

        assert_eq!(writer.frames_written(), 1);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn is_buffered_accessor() {
        let temp_dir = std::env::temp_dir();
        let temp_path1 = temp_dir.join("not_buffered_test.ser");
        let temp_path2 = temp_dir.join("buffered_test.ser");

        let header = create_test_header();

        let writer1 = SerWriter::create(&temp_path1, header.clone()).unwrap();
        assert!(!writer1.is_buffered());

        let writer2 = SerWriter::with_buffer(&temp_path2, header, 1).unwrap();
        assert!(writer2.is_buffered());

        std::fs::remove_file(temp_path1).ok();
        std::fs::remove_file(temp_path2).ok();
    }

    #[test]
    fn buffer_usage_no_buffer() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("no_buffer_usage_test.ser");

        let header = create_test_header();
        let writer = SerWriter::create(&temp_path, header).unwrap();

        let usage = writer.buffer_usage();
        assert!(usage.is_none());

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn buffer_usage_with_buffer() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("buffer_usage_test.ser");

        let header = create_test_header();
        let mut writer = SerWriter::with_buffer(&temp_path, header, 1).unwrap();

        let frame_data = create_test_frame_data();
        writer.write_frame(&frame_data).unwrap();

        let usage = writer.buffer_usage();
        assert!(usage.is_some());
        let (available, _) = usage.unwrap();
        assert_eq!(available, 1);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn drop_flushes_buffer() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("drop_test.ser");

        let header = create_test_header();
        {
            let mut writer = SerWriter::with_buffer(&temp_path, header, 1).unwrap();
            let frame_data = create_test_frame_data();
            writer.write_frame(&frame_data).unwrap();
        }

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn write_and_read_round_trip() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("round_trip_test.ser");

        let header = create_test_header();
        let mut writer = SerWriter::create(&temp_path, header).unwrap();

        let frame_data1 = vec![1u8; 200];
        let frame_data2 = vec![2u8; 200];

        let timestamp1 = SerTimestamp::new(111111);
        let timestamp2 = SerTimestamp::new(222222);

        writer
            .write_frame_with_timestamp(&frame_data1, Some(timestamp1))
            .unwrap();
        writer
            .write_frame_with_timestamp(&frame_data2, Some(timestamp2))
            .unwrap();

        writer.finalize().unwrap();

        let mut file = File::open(&temp_path).unwrap();
        let mut header_bytes = [0u8; 178];
        file.read_exact(&mut header_bytes).unwrap();

        let read_header = SerHeader::from_bytes(&header_bytes).unwrap();
        assert_eq!(read_header.frame_count, 2);

        std::fs::remove_file(temp_path).ok();
    }

    #[test]
    fn complete_workflow_test() {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join("complete_workflow_test.ser");

        let header = create_test_header();
        let mut writer = SerWriter::with_buffer(&temp_path, header, 2).unwrap();

        let frame1 = vec![1u8; 200];
        let frame2 = vec![2u8; 200];
        let frame3 = vec![3u8; 200];

        writer.write_frame(&frame1).unwrap();
        writer.write_frame(&frame2).unwrap();
        writer.flush_buffer().unwrap();
        writer.write_frame(&frame3).unwrap();

        writer.finalize().unwrap();

        std::fs::remove_file(temp_path).ok();
    }
}
