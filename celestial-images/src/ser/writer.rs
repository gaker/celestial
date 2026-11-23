use crate::ser::error::{Result, SerError};
use crate::ser::header::SerHeader;
use crate::ser::timestamp::SerTimestamp;
use crate::ser::types::SerFrameId;
use std::fs::File;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

const BYTES_PER_MEGABYTE: usize = 1 << 20;
// BufWriter allocates the whole buffer up front, and a failed allocation aborts the process.
const MAX_BUFFER_MEGABYTES: usize = 1024;

pub struct SerWriter {
    inner: FrameWriter<File>,
}

impl SerWriter {
    pub fn create<P: AsRef<Path>>(path: P, header: SerHeader) -> Result<Self> {
        let (header, frame_size) = empty_header(header)?;
        let writer = BufWriter::new(File::create(path)?);
        Ok(Self {
            inner: FrameWriter::start(writer, header, frame_size)?,
        })
    }

    // The buffer holds `megabytes` MiB of frame data before it writes to the file.
    pub fn with_buffer<P: AsRef<Path>>(
        path: P,
        header: SerHeader,
        megabytes: usize,
    ) -> Result<Self> {
        if megabytes > MAX_BUFFER_MEGABYTES {
            return Err(SerError::BufferTooLarge { megabytes });
        }
        let capacity = megabytes * BYTES_PER_MEGABYTE;
        let (header, frame_size) = empty_header(header)?;
        let writer = BufWriter::with_capacity(capacity, File::create(path)?);
        Ok(Self {
            inner: FrameWriter::start(writer, header, frame_size)?,
        })
    }

    pub fn header(&self) -> &SerHeader {
        &self.inner.header
    }

    pub fn write_frame(&mut self, data: &[u8]) -> Result<SerFrameId> {
        self.inner.write_frame(data, None)
    }

    // Either every frame has a timestamp or none does.
    pub fn write_frame_with_timestamp(
        &mut self,
        data: &[u8],
        timestamp: Option<SerTimestamp>,
    ) -> Result<SerFrameId> {
        self.inner.write_frame(data, timestamp)
    }

    // Also rewrites the header, so the file on disk holds every frame written
    // so far. The timestamp trailer is written only by `finalize`.
    pub fn flush_buffer(&mut self) -> Result<()> {
        self.inner.flush()
    }

    pub fn finalize(mut self) -> Result<()> {
        self.inner.finish()
    }

    pub fn frames_written(&self) -> u32 {
        self.inner.header.frame_count
    }
}

struct FrameWriter<W: Write + Seek> {
    // After an I/O error this holds the error's message. The file may then end
    // in a partial frame, so nothing more is written to it.
    writer: std::result::Result<BufWriter<W>, String>,
    header: SerHeader,
    frame_size: usize,
    timestamps: Vec<SerTimestamp>,
    frame0_has_timestamp: Option<bool>,
    finished: bool,
}

impl<W: Write + Seek> FrameWriter<W> {
    fn start(writer: BufWriter<W>, header: SerHeader, frame_size: usize) -> Result<Self> {
        let mut writer = Self {
            writer: Ok(writer),
            header,
            frame_size,
            timestamps: Vec::new(),
            frame0_has_timestamp: None,
            finished: false,
        };
        let header = writer.header.to_bytes();
        writer.io(|w| w.write_all(&header))?;
        Ok(writer)
    }

    fn write_frame(&mut self, data: &[u8], timestamp: Option<SerTimestamp>) -> Result<SerFrameId> {
        self.check_writable()?;
        self.check_frame(data, timestamp.is_some())?;
        let frame_id = self.header.frame_count;
        let frame_count = frame_id.checked_add(1).ok_or(SerError::TooManyFrames)?;

        self.io(|w| w.write_all(data))?;
        self.timestamps.extend(timestamp);
        self.frame0_has_timestamp = Some(timestamp.is_some());
        self.header.frame_count = frame_count;
        Ok(frame_id)
    }

    fn check_frame(&self, data: &[u8], has_timestamp: bool) -> Result<()> {
        if data.len() != self.frame_size {
            return Err(SerError::BufferSizeMismatch {
                expected: self.frame_size,
                actual: data.len(),
            });
        }
        match self.frame0_has_timestamp {
            Some(frame0_has_timestamp) if frame0_has_timestamp != has_timestamp => {
                Err(SerError::TimestampPresence {
                    frame: self.header.frame_count,
                    frame0_has_timestamp,
                })
            }
            _ => Ok(()),
        }
    }

    fn flush(&mut self) -> Result<()> {
        self.write_header()?;
        self.io(|w| w.flush())
    }

    // The header goes before the trailer, so a failure while writing the
    // trailer leaves a header that counts every frame.
    fn finish(&mut self) -> Result<()> {
        self.finished = true;
        self.check_writable()?;
        if let Some(first) = self.timestamps.first() {
            // A DateTime of 0 or less tells readers that there is no trailer.
            let start = first.ticks() as i64;
            if self.header.date_time <= 0 {
                self.header.date_time = start;
            }
            if self.header.date_time_utc <= 0 {
                self.header.date_time_utc = start;
            }
        }
        self.write_header()?;
        let trailer: Vec<u8> = self
            .timestamps
            .iter()
            .flat_map(|t| t.to_le_bytes())
            .collect();
        self.io(|w| w.write_all(&trailer))?;
        self.io(|w| w.flush())
    }

    // Seeking flushes the buffered frames before the header is written.
    fn write_header(&mut self) -> Result<()> {
        let header = self.header.to_bytes();
        self.io(|w| {
            w.seek(SeekFrom::Start(0))?;
            w.write_all(&header)?;
            w.seek(SeekFrom::End(0))
        })?;
        Ok(())
    }

    fn check_writable(&self) -> Result<()> {
        match &self.writer {
            Ok(_) => Ok(()),
            Err(message) => Err(SerError::WriterFailed(message.clone())),
        }
    }

    fn io<T>(&mut self, op: impl FnOnce(&mut BufWriter<W>) -> io::Result<T>) -> Result<T> {
        let writer = match &mut self.writer {
            Ok(writer) => writer,
            Err(message) => return Err(SerError::WriterFailed(message.clone())),
        };
        let error = match op(writer) {
            Ok(value) => return Ok(value),
            Err(error) => error,
        };
        if let Ok(writer) = std::mem::replace(&mut self.writer, Err(error.to_string())) {
            // into_parts discards the buffered bytes. Dropping the BufWriter
            // would write them after the partial frame.
            drop(writer.into_parts());
        }
        Err(SerError::Io(error))
    }
}

// The frame count starts at 0 and counts the frames written.
fn empty_header(header: SerHeader) -> Result<(SerHeader, usize)> {
    let header = SerHeader {
        frame_count: 0,
        ..header
    };
    header.validate()?;
    let frame_size = header.frame_size()?;
    let frame_size = usize::try_from(frame_size).map_err(|_| {
        SerError::InvalidHeader(format!("Frame size {} does not fit in memory", frame_size))
    })?;
    Ok((header, frame_size))
}

impl<W: Write + Seek> Drop for FrameWriter<W> {
    fn drop(&mut self) {
        if !self.finished && self.writer.is_ok() {
            // Drop cannot return an error; call finalize to see one.
            let _ = self.finish();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::types::ByteOrder;
    use crate::ser::header::HEADER_SIZE;
    use crate::ser::reader::SerReader;
    use crate::ser::types::ColorId;
    use std::cell::{Cell, RefCell};
    use std::path::PathBuf;
    use std::rc::Rc;
    use tempfile::TempDir;

    // 2024-06-01T12:34:56.1234567 in .NET ticks.
    const T1: u64 = 638_528_420_961_234_567;

    fn create_test_header(date_time: i64) -> SerHeader {
        SerHeader {
            image_width: 10,
            image_height: 10,
            pixel_depth_per_plane: 16,
            frame_count: 7,
            date_time,
            ..Default::default()
        }
    }

    fn temp_path() -> (TempDir, PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.ser");
        (dir, path)
    }

    fn frame(value: u8) -> Vec<u8> {
        vec![value; 200]
    }

    #[test]
    fn create_writes_an_empty_header() {
        let (_dir, path) = temp_path();
        let mut writer = SerWriter::create(&path, create_test_header(0)).unwrap();

        assert_eq!(writer.header().image_width, 10);
        assert_eq!(writer.header().frame_count, 0);
        assert_eq!(writer.frames_written(), 0);
        writer.flush_buffer().unwrap();
        let expected = SerHeader {
            frame_count: 0,
            ..create_test_header(0)
        };
        assert_eq!(std::fs::read(&path).unwrap(), expected.to_bytes());
    }

    #[test]
    fn create_writer_file_error() {
        let result = SerWriter::create(
            "/invalid/path/that/does/not/exist/test.ser",
            create_test_header(0),
        );
        assert!(matches!(result, Err(SerError::Io(e)) if e.kind() == std::io::ErrorKind::NotFound));
    }

    #[test]
    fn invalid_header_is_rejected_before_the_file_is_created() {
        let (_dir, path) = temp_path();
        let header = SerHeader {
            image_width: 0,
            ..create_test_header(0)
        };
        let result = SerWriter::create(&path, header.clone());
        assert!(matches!(
            result,
            Err(SerError::InvalidDimensions {
                width: 0,
                height: 10
            })
        ));
        let result = SerWriter::with_buffer(&path, header, 1);
        assert!(matches!(result, Err(SerError::InvalidDimensions { .. })));
        assert!(!path.exists());
    }

    #[test]
    fn with_buffer_rejects_sizes_over_1024_mib() {
        let (_dir, path) = temp_path();
        std::fs::write(&path, b"keep").unwrap();
        for megabytes in [1025, usize::MAX] {
            let result = SerWriter::with_buffer(&path, create_test_header(0), megabytes);
            assert!(
                matches!(result, Err(SerError::BufferTooLarge { megabytes: m }) if m == megabytes),
                "{megabytes}"
            );
        }
        assert_eq!(std::fs::read(&path).unwrap(), b"keep");
    }

    #[test]
    fn write_frame_wrong_size() {
        let (_dir, path) = temp_path();
        let mut writer = SerWriter::create(&path, create_test_header(0)).unwrap();

        let result = writer.write_frame_with_timestamp(&[0u8; 100], None);
        assert!(matches!(
            result,
            Err(SerError::BufferSizeMismatch {
                expected: 200,
                actual: 100
            })
        ));
        assert_eq!(writer.frames_written(), 0);
    }

    #[test]
    fn buffered_writer_keeps_every_frame() {
        let (_dir, path) = temp_path();
        let header = SerHeader {
            image_width: 640,
            image_height: 480,
            ..create_test_header(0)
        };
        // One frame is 614,400 bytes, so a 1 MiB buffer holds one frame.
        let mut writer = SerWriter::with_buffer(&path, header, 1).unwrap();
        let frames: Vec<Vec<u8>> = (0..3).map(|i| vec![i as u8 + 1; 640 * 480 * 2]).collect();
        for (i, data) in frames.iter().enumerate() {
            assert_eq!(writer.write_frame(data).unwrap(), i as u32);
        }
        writer.finalize().unwrap();

        assert_eq!(std::fs::metadata(&path).unwrap().len(), 178 + 3 * 614_400);
        let reader = SerReader::open(&path).unwrap();
        assert_eq!(reader.frame_count(), 3);
        for (i, data) in frames.iter().enumerate() {
            assert_eq!(reader.frame_slice(i as u32).unwrap(), data.as_slice());
        }
    }

    #[test]
    fn zero_megabyte_buffer_writes_through() {
        let (_dir, path) = temp_path();
        let mut writer = SerWriter::with_buffer(&path, create_test_header(0), 0).unwrap();
        writer.write_frame(&frame(5)).unwrap();
        assert_eq!(std::fs::metadata(&path).unwrap().len(), 378);
        writer.finalize().unwrap();
        assert_eq!(
            SerReader::open(&path).unwrap().frame_slice(0).unwrap(),
            frame(5)
        );
    }

    #[test]
    fn mixed_timestamps_are_rejected() {
        let (_dir, path) = temp_path();
        let mut writer = SerWriter::create(&path, create_test_header(0)).unwrap();
        let ts = SerTimestamp::new(T1);
        writer
            .write_frame_with_timestamp(&frame(1), Some(ts))
            .unwrap();
        let result = writer.write_frame_with_timestamp(&frame(2), None);
        assert!(matches!(
            result,
            Err(SerError::TimestampPresence {
                frame: 1,
                frame0_has_timestamp: true
            })
        ));
        writer
            .write_frame_with_timestamp(&frame(3), Some(ts))
            .unwrap();
        writer.finalize().unwrap();

        let reader = SerReader::open(&path).unwrap();
        assert_eq!(reader.frame_count(), 2);
        assert_eq!(reader.frame_slice(1).unwrap(), frame(3));
        assert_eq!(reader.read_timestamp(1).unwrap(), Some(ts));
    }

    #[test]
    fn timestamp_after_untimed_frame_is_rejected() {
        let (_dir, path) = temp_path();
        let mut writer = SerWriter::create(&path, create_test_header(0)).unwrap();
        writer.write_frame(&frame(1)).unwrap();
        let result = writer.write_frame_with_timestamp(&frame(2), Some(SerTimestamp::new(T1)));
        assert!(matches!(
            result,
            Err(SerError::TimestampPresence {
                frame: 1,
                frame0_has_timestamp: false
            })
        ));
        assert_eq!(writer.frames_written(), 1);
    }

    #[test]
    fn frame_count_limit_is_an_error() {
        let (_dir, path) = temp_path();
        let mut writer = SerWriter::create(&path, create_test_header(0)).unwrap();
        writer.inner.header.frame_count = u32::MAX;
        assert!(matches!(
            writer.write_frame(&frame(1)),
            Err(SerError::TooManyFrames)
        ));
        writer.inner.header.frame_count = 0;
    }

    #[test]
    fn trailer_is_written_when_date_time_is_zero() {
        let (_dir, path) = temp_path();
        let mut writer = SerWriter::create(&path, create_test_header(0)).unwrap();
        let stamps = [SerTimestamp::new(T1), SerTimestamp::new(T1 + 400_000)];
        for (i, ts) in stamps.iter().enumerate() {
            writer
                .write_frame_with_timestamp(&frame(i as u8), Some(*ts))
                .unwrap();
        }
        writer.finalize().unwrap();

        let reader = SerReader::open(&path).unwrap();
        assert!(reader.has_timestamps());
        assert_eq!(reader.header().date_time, T1 as i64);
        assert_eq!(reader.header().date_time_utc, T1 as i64);
        assert_eq!(reader.read_timestamp(0).unwrap(), Some(stamps[0]));
        assert_eq!(reader.read_timestamp(1).unwrap(), Some(stamps[1]));
    }

    #[test]
    fn header_date_time_is_kept_when_set() {
        let (_dir, path) = temp_path();
        let header = SerHeader {
            date_time_utc: 5,
            ..create_test_header(9)
        };
        let mut writer = SerWriter::create(&path, header).unwrap();
        writer
            .write_frame_with_timestamp(&frame(1), Some(SerTimestamp::new(T1)))
            .unwrap();
        writer.finalize().unwrap();

        let reader = SerReader::open(&path).unwrap();
        assert_eq!(reader.header().date_time, 9);
        assert_eq!(reader.header().date_time_utc, 5);
    }

    #[test]
    fn no_trailer_without_timestamps() {
        let (_dir, path) = temp_path();
        let mut writer = SerWriter::create(&path, create_test_header(T1 as i64)).unwrap();
        writer.write_frame(&frame(1)).unwrap();
        writer.finalize().unwrap();

        assert_eq!(std::fs::metadata(&path).unwrap().len(), 378);
        let reader = SerReader::open(&path).unwrap();
        assert!(!reader.has_timestamps());
        assert_eq!(reader.header().date_time, T1 as i64);
    }

    #[test]
    fn drop_finalizes() {
        let (_dir, path) = temp_path();
        {
            let mut writer = SerWriter::create(&path, create_test_header(0)).unwrap();
            writer
                .write_frame_with_timestamp(&frame(7), Some(SerTimestamp::new(T1)))
                .unwrap();
        }

        let reader = SerReader::open(&path).unwrap();
        assert_eq!(reader.frame_count(), 1);
        assert_eq!(reader.frame_slice(0).unwrap(), frame(7));
        assert_eq!(
            reader.read_timestamp(0).unwrap(),
            Some(SerTimestamp::new(T1))
        );
    }

    #[test]
    fn flush_buffer_leaves_a_readable_file() {
        let (_dir, path) = temp_path();
        let mut writer = SerWriter::create(&path, create_test_header(0)).unwrap();
        let ts = SerTimestamp::new(T1);
        writer
            .write_frame_with_timestamp(&frame(1), Some(ts))
            .unwrap();
        writer.flush_buffer().unwrap();

        let reader = SerReader::open(&path).unwrap();
        assert_eq!(reader.frame_count(), 1);
        assert!(!reader.has_timestamps());

        writer
            .write_frame_with_timestamp(&frame(2), Some(ts))
            .unwrap();
        writer.finalize().unwrap();
        let reader = SerReader::open(&path).unwrap();
        assert_eq!(reader.frame_count(), 2);
        assert_eq!(reader.frame_slice(1).unwrap(), frame(2));
        assert_eq!(reader.read_timestamp(1).unwrap(), Some(ts));
    }

    #[test]
    fn write_and_read_round_trip() {
        let (_dir, path) = temp_path();
        let header = SerHeader {
            color_id: ColorId::Rgb,
            byte_order: ByteOrder::BigEndian,
            image_width: 3,
            image_height: 2,
            pixel_depth_per_plane: 12,
            observer: "\u{C5}sa \u{2603}".to_string(),
            telescope: "x".repeat(45),
            ..create_test_header(0)
        };
        let mut writer = SerWriter::create(&path, header).unwrap();
        let frames: Vec<Vec<u8>> = (0..3u8)
            .map(|i| (0..36).map(|j| i * 36 + j).collect())
            .collect();
        let stamps: Vec<SerTimestamp> = (0..3).map(|i| SerTimestamp::new(T1 + i)).collect();
        for (data, ts) in frames.iter().zip(&stamps) {
            writer.write_frame_with_timestamp(data, Some(*ts)).unwrap();
        }
        writer.finalize().unwrap();

        let reader = SerReader::open(&path).unwrap();
        let read = reader.header();
        assert_eq!(read.color_id, ColorId::Rgb);
        assert_eq!(read.byte_order, ByteOrder::BigEndian);
        assert_eq!(read.pixel_depth_per_plane, 12);
        assert_eq!(read.frame_count, 3);
        assert_eq!(read.observer, "?sa ?");
        assert_eq!(read.telescope, "x".repeat(40));
        for (i, (data, ts)) in frames.iter().zip(&stamps).enumerate() {
            let frame = reader.read_frame(i as u32).unwrap();
            assert_eq!(frame.data, data.as_slice());
            assert_eq!(frame.timestamp, Some(*ts));
        }
        let samples = reader.read_frame_u16(0).unwrap();
        assert_eq!(samples[..3], [0x0001, 0x0203, 0x0405]);
    }

    // An in-memory file that accepts bytes up to `limit` and then fails, as a
    // full disk or a file size limit does. Tests raise `limit` to model a
    // failure that goes away.
    #[derive(Clone, Default)]
    struct LimitedFile {
        bytes: Rc<RefCell<Vec<u8>>>,
        limit: Rc<Cell<usize>>,
        position: usize,
    }

    impl Write for LimitedFile {
        fn write(&mut self, data: &[u8]) -> io::Result<usize> {
            let room = self.limit.get().saturating_sub(self.position);
            if room == 0 && !data.is_empty() {
                return Err(io::Error::other("file size limit exceeded"));
            }
            let count = data.len().min(room);
            let end = self.position + count;
            let mut bytes = self.bytes.borrow_mut();
            if bytes.len() < end {
                bytes.resize(end, 0);
            }
            bytes[self.position..end].copy_from_slice(&data[..count]);
            self.position = end;
            Ok(count)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    impl Seek for LimitedFile {
        fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
            let length = self.bytes.borrow().len() as i64;
            let position = match from {
                SeekFrom::Start(offset) => offset as i64,
                SeekFrom::End(delta) => length + delta,
                SeekFrom::Current(delta) => self.position as i64 + delta,
            };
            self.position = position as usize;
            Ok(position as u64)
        }
    }

    fn limited_writer(limit: usize, capacity: usize) -> (FrameWriter<LimitedFile>, LimitedFile) {
        let file = LimitedFile::default();
        file.limit.set(limit);
        let (header, frame_size) = empty_header(create_test_header(0)).unwrap();
        let writer = BufWriter::with_capacity(capacity, file.clone());
        let writer = FrameWriter::start(writer, header, frame_size).unwrap();
        (writer, file)
    }

    fn written_header(file: &LimitedFile) -> SerHeader {
        SerHeader::from_bytes(&file.bytes.borrow()[..HEADER_SIZE]).unwrap()
    }

    fn assert_stopped(result: Result<()>) {
        assert_eq!(
            result.map_err(|e| e.to_string()),
            Err("SER writer stopped after an earlier I/O error: file size limit exceeded".into())
        );
    }

    #[test]
    fn header_is_rewritten_before_the_trailer() {
        let frames_end = HEADER_SIZE + 2 * 200;
        let (mut writer, file) = limited_writer(frames_end + 4, 64);
        writer
            .write_frame(&frame(1), Some(SerTimestamp::new(T1)))
            .unwrap();
        writer
            .write_frame(&frame(2), Some(SerTimestamp::new(T1 + 1)))
            .unwrap();

        assert!(matches!(writer.finish(), Err(SerError::Io(_))));
        let header = written_header(&file);
        assert_eq!(header.frame_count, 2);
        assert_eq!(header.date_time, T1 as i64);
        assert_eq!(header.date_time_utc, T1 as i64);
        assert_eq!(file.bytes.borrow().len(), frames_end + 4);
    }

    #[test]
    fn a_frame_that_fails_partway_stops_the_writer() {
        let limit = HEADER_SIZE + 200 + 50;
        let (mut writer, file) = limited_writer(limit, 64);
        writer.write_frame(&frame(1), None).unwrap();
        assert!(matches!(
            writer.write_frame(&frame(2), None),
            Err(SerError::Io(_))
        ));

        file.limit.set(usize::MAX);
        assert_stopped(writer.write_frame(&frame(3), None).map(|_| ()));
        assert_stopped(writer.write_frame(&[0; 3], None).map(|_| ()));
        assert_stopped(writer.flush());
        assert_stopped(writer.finish());
        drop(writer);
        assert_eq!(file.bytes.borrow().len(), limit);
        assert_eq!(written_header(&file).frame_count, 0);
    }

    #[test]
    fn buffered_frames_report_the_error_when_flushed() {
        let limit = HEADER_SIZE + 300;
        let (mut writer, file) = limited_writer(limit, BYTES_PER_MEGABYTE);
        for value in 1..=3 {
            writer.write_frame(&frame(value), None).unwrap();
        }

        assert!(matches!(writer.finish(), Err(SerError::Io(_))));
        file.limit.set(usize::MAX);
        assert_stopped(writer.finish());
        drop(writer);
        assert_eq!(file.bytes.borrow().len(), limit);
        assert_eq!(written_header(&file).frame_count, 0);
    }

    #[test]
    fn drop_does_nothing_after_an_error() {
        let limit = HEADER_SIZE + 300;
        let (mut writer, file) = limited_writer(limit, 0);
        writer.write_frame(&frame(1), None).unwrap();
        assert!(matches!(
            writer.write_frame(&frame(2), None),
            Err(SerError::Io(_))
        ));
        file.limit.set(usize::MAX);
        drop(writer);
        assert_eq!(file.bytes.borrow().len(), limit);
        assert_eq!(written_header(&file).frame_count, 0);
    }
}
