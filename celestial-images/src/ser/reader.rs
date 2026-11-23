use crate::core::types::ByteOrder;
use crate::ser::error::{Result, SerError};
use crate::ser::header::{SerHeader, HEADER_SIZE, TIMESTAMP_SIZE};
use crate::ser::timestamp::SerTimestamp;
use crate::ser::types::{SerFrame, SerFrameId};
use memmap2::Mmap;
use std::fs::File;
use std::ops::Range;
use std::path::Path;

pub struct SerReader {
    header: SerHeader,
    mmap: Mmap,
    has_trailer: bool,
}

impl SerReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let length = file.metadata()?.len();

        if length < HEADER_SIZE as u64 {
            return Err(SerError::FileTruncated {
                expected: HEADER_SIZE as u64,
                actual: length,
            });
        }

        let mmap = unsafe { Mmap::map(&file)? };
        let header = SerHeader::from_bytes(&mmap)?;
        let has_trailer = trailer_present(&header, length)?;

        Ok(Self {
            header,
            mmap,
            has_trailer,
        })
    }

    pub fn header(&self) -> &SerHeader {
        &self.header
    }

    pub fn frame_count(&self) -> u32 {
        self.header.frame_count
    }

    pub fn has_timestamps(&self) -> bool {
        self.has_trailer
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
        Ok(SerFrame {
            id: frame_id,
            data: self.frame_slice(frame_id)?,
            timestamp: self.read_timestamp(frame_id)?,
        })
    }

    pub fn read_frame_data(&self, frame_id: SerFrameId) -> Result<Vec<u8>> {
        Ok(self.frame_slice(frame_id)?.to_vec())
    }

    // Samples of depth 9 to 15 are returned as stored, in the low bits.
    pub fn read_frame_u16(&self, frame_id: SerFrameId) -> Result<Vec<u16>> {
        let depth = self.header.pixel_depth_per_plane;
        if depth <= 8 {
            return Err(SerError::NotSixteenBit { depth });
        }
        let (samples, _) = self.frame_slice(frame_id)?.as_chunks::<2>();
        let decode = match self.header.byte_order {
            ByteOrder::LittleEndian => u16::from_le_bytes,
            ByteOrder::BigEndian => u16::from_be_bytes,
        };
        Ok(samples.iter().map(|&sample| decode(sample)).collect())
    }

    pub fn read_timestamp(&self, frame_id: SerFrameId) -> Result<Option<SerTimestamp>> {
        self.check_frame_id(frame_id)?;
        if !self.has_trailer {
            return Ok(None);
        }
        let start = self.header.timestamp_offset(frame_id)?;
        let bytes = self.mapped_bytes(start, TIMESTAMP_SIZE)?;
        let mut ticks = [0u8; 8];
        ticks.copy_from_slice(bytes);
        Ok(Some(SerTimestamp::from_le_bytes(ticks)))
    }

    pub fn iter_frames(&self) -> FrameIterator<'_> {
        FrameIterator {
            reader: self,
            ids: 0..self.header.frame_count,
        }
    }

    pub fn iter_frames_range(&self, start: SerFrameId, count: u32) -> Result<FrameIterator<'_>> {
        Ok(FrameIterator {
            reader: self,
            ids: self.frame_ids(start, count)?,
        })
    }

    pub fn frame_slice(&self, frame_id: SerFrameId) -> Result<&[u8]> {
        self.check_frame_id(frame_id)?;
        let start = self.header.frame_offset(frame_id)?;
        self.mapped_bytes(start, self.header.frame_size()?)
    }

    pub fn frame_slices(&self, start: SerFrameId, count: u32) -> Result<Vec<(SerFrameId, &[u8])>> {
        self.frame_ids(start, count)?
            .map(|frame_id| Ok((frame_id, self.frame_slice(frame_id)?)))
            .collect()
    }

    // `count` is cut to the frames that exist; a start past the last frame is an error.
    fn frame_ids(&self, start: SerFrameId, count: u32) -> Result<Range<SerFrameId>> {
        let total = self.header.frame_count;
        if start > total {
            return Err(SerError::FrameOutOfBounds {
                frame: start,
                total,
            });
        }
        Ok(start..start.saturating_add(count).min(total))
    }

    fn check_frame_id(&self, frame_id: SerFrameId) -> Result<()> {
        if frame_id >= self.header.frame_count {
            return Err(SerError::FrameOutOfBounds {
                frame: frame_id,
                total: self.header.frame_count,
            });
        }
        Ok(())
    }

    fn mapped_bytes(&self, start: u64, length: u64) -> Result<&[u8]> {
        let end = start.saturating_add(length);
        let truncated = SerError::FileTruncated {
            expected: end,
            actual: self.mmap.len() as u64,
        };
        match (usize::try_from(start), usize::try_from(end)) {
            (Ok(start), Ok(end)) => self.mmap.get(start..end).ok_or(truncated),
            _ => Err(truncated),
        }
    }
}

// The frames must all be present. A trailer is read only when it holds a
// timestamp for every frame, and bytes after the frames or the trailer are
// ignored. The header's DateTime does not decide this.
fn trailer_present(header: &SerHeader, length: u64) -> Result<bool> {
    let frames_end = header.frames_end()?;
    if length < frames_end {
        return Err(SerError::FileTruncated {
            expected: frames_end,
            actual: length,
        });
    }
    Ok(header.frame_count > 0 && length >= header.trailer_end()?)
}

pub struct FrameIterator<'a> {
    reader: &'a SerReader,
    ids: Range<SerFrameId>,
}

impl<'a> Iterator for FrameIterator<'a> {
    type Item = Result<SerFrame<'a>>;

    fn next(&mut self) -> Option<Self::Item> {
        self.ids
            .next()
            .map(|frame_id| self.reader.read_frame(frame_id))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.ids.size_hint()
    }
}

impl ExactSizeIterator for FrameIterator<'_> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ser::types::ColorId;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // 2024-06-01T12:34:56.1234567 in .NET ticks.
    const T1: u64 = 638_528_420_961_234_567;

    fn shown<T>(result: Result<T>) -> std::result::Result<T, String> {
        result.map_err(|e| e.to_string())
    }

    fn error_text<T>(result: Result<T>) -> Option<String> {
        result.err().map(|e| e.to_string())
    }

    fn test_header(frame_count: u32, date_time: i64) -> SerHeader {
        SerHeader {
            lu_id: 123,
            byte_order: ByteOrder::BigEndian,
            image_width: 10,
            image_height: 10,
            pixel_depth_per_plane: 16,
            frame_count,
            observer: "Observer".to_string(),
            instrument: "Instrument".to_string(),
            telescope: "Telescope".to_string(),
            date_time,
            date_time_utc: date_time,
            ..SerHeader::default()
        }
    }

    // Frame i holds bytes (i * 200 + offset) % 256.
    fn ser_bytes(header: &SerHeader, frames: u32, trailer: &[u64]) -> Vec<u8> {
        let mut data = header.to_bytes().to_vec();
        let frame_size = header.frame_size().unwrap() as usize;
        let pixels = frames as usize * frame_size;
        data.extend((0..pixels).map(|i| (i % 256) as u8));
        for ticks in trailer {
            data.extend_from_slice(&ticks.to_le_bytes());
        }
        data
    }

    // 2 frames of 10x10 16-bit mono with a timestamp trailer.
    fn create_test_ser_file() -> Vec<u8> {
        ser_bytes(&test_header(2, T1 as i64), 2, &[T1, T1 + 1])
    }

    // 1 frame and no trailer.
    fn create_minimal_ser_file() -> Vec<u8> {
        ser_bytes(&test_header(1, 0), 1, &[])
    }

    fn write_temp(data: &[u8]) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(data).unwrap();
        file.flush().unwrap();
        file
    }

    fn open_bytes(data: &[u8]) -> Result<SerReader> {
        SerReader::open(write_temp(data).path())
    }

    #[test]
    fn open_rejects_a_bad_file_id() {
        let mut data = create_minimal_ser_file();
        data[0..14].copy_from_slice(b"INVALID-HEADER");
        let result = open_bytes(&data);
        assert!(matches!(result, Err(SerError::InvalidFileId { .. })));
    }

    #[test]
    fn header_and_frame_count_accessors() {
        let reader = open_bytes(&create_minimal_ser_file()).unwrap();
        assert_eq!(reader.header().frame_count, 1);
        assert_eq!(reader.header().observer, "Observer");
        assert_eq!(reader.frame_count(), 1);
    }

    #[test]
    fn color_kind_accessors() {
        let cases = [
            (ColorId::Mono, (false, true, false)),
            (ColorId::BayerGbrg, (false, false, true)),
            (ColorId::Bgr, (true, false, false)),
        ];
        for (color_id, kinds) in cases {
            let header = SerHeader {
                color_id,
                ..test_header(0, 0)
            };
            let reader = open_bytes(&header.to_bytes()).unwrap();
            assert_eq!(
                (reader.is_color(), reader.is_mono(), reader.is_bayer()),
                kinds,
                "{color_id}"
            );
        }
    }

    #[test]
    fn frame_out_of_bounds() {
        let file = write_temp(&create_test_ser_file());
        let reader = SerReader::open(file.path()).unwrap();
        let out_of_bounds = |frame| {
            Some(format!(
                "Frame index {frame} out of bounds (total frames: 2)"
            ))
        };

        assert_eq!(error_text(reader.read_frame(2)), out_of_bounds(2));
        assert_eq!(error_text(reader.frame_slice(999)), out_of_bounds(999));
        assert_eq!(error_text(reader.read_timestamp(999)), out_of_bounds(999));
    }

    #[test]
    fn open_file_success() {
        let data = create_test_ser_file();
        let reader = open_bytes(&data).unwrap();

        assert_eq!(reader.frame_count(), 2);
        assert_eq!(reader.header().instrument, "Instrument");
        assert!(reader.has_timestamps());
        assert_eq!(reader.mmap.len(), data.len());
    }

    #[test]
    fn open_file_no_trailer() {
        let file = write_temp(&create_minimal_ser_file());
        let reader = SerReader::open(file.path()).unwrap();

        assert_eq!(reader.frame_count(), 1);
        assert!(!reader.has_timestamps());
        assert_eq!(shown(reader.read_timestamp(0)), Ok(None));
    }

    #[test]
    fn open_file_too_short() {
        let result = open_bytes(&[0u8; 100]);
        assert!(matches!(
            result,
            Err(SerError::FileTruncated {
                expected: 178,
                actual: 100
            })
        ));
    }

    #[test]
    fn frames_cut_short_are_truncated() {
        let mut data = create_test_ser_file();
        data.truncate(300);
        assert!(matches!(
            open_bytes(&data),
            Err(SerError::FileTruncated {
                expected: 578,
                actual: 300
            })
        ));
    }

    fn timestamps(reader: &SerReader) -> Vec<Option<SerTimestamp>> {
        (0..reader.frame_count())
            .map(|id| reader.read_timestamp(id).unwrap())
            .collect()
    }

    #[test]
    fn partial_trailer_opens_without_timestamps() {
        let full = create_test_ser_file();
        let reader = open_bytes(&full[..590]).unwrap();
        assert_eq!(reader.frame_count(), 2);
        assert!(!reader.has_timestamps());
        assert_eq!(timestamps(&reader), [None, None]);
        assert_eq!(reader.read_frame_data(1).unwrap(), &full[378..578]);
    }

    #[test]
    fn bytes_after_the_trailer_are_ignored() {
        let mut data = create_test_ser_file();
        data.extend([0xAB; 9]);
        let reader = open_bytes(&data).unwrap();
        assert_eq!(reader.frame_count(), 2);
        let expected = [Some(SerTimestamp::new(T1)), Some(SerTimestamp::new(T1 + 1))];
        assert_eq!(timestamps(&reader), expected);
    }

    #[test]
    fn bytes_after_the_frames_are_ignored() {
        let frames = ser_bytes(&test_header(2, 0), 2, &[]);
        for extra in [1, 9, 15] {
            let mut data = frames.clone();
            data.extend(vec![0xAB; extra]);
            let reader = open_bytes(&data).unwrap();
            assert_eq!(reader.frame_count(), 2, "{extra}");
            assert_eq!(timestamps(&reader), [None, None], "{extra}");
            assert_eq!(reader.read_frame_data(1).unwrap(), &frames[378..578]);
        }
    }

    #[test]
    fn zero_frame_count_with_frame_data_opens_empty() {
        let mut data = ser_bytes(&test_header(2, 0), 2, &[T1, T1 + 1]);
        data[38..42].copy_from_slice(&0u32.to_le_bytes());
        let reader = open_bytes(&data).unwrap();
        assert_eq!(reader.frame_count(), 0);
        assert!(!reader.has_timestamps());
        assert_eq!(reader.iter_frames().count(), 0);
        assert_eq!(
            error_text(reader.read_frame(0)),
            Some("Frame index 0 out of bounds (total frames: 0)".to_string())
        );
    }

    #[test]
    fn positive_date_time_without_a_trailer_opens() {
        let data = ser_bytes(&test_header(3, T1 as i64), 3, &[]);
        let reader = open_bytes(&data).unwrap();
        assert!(!reader.has_timestamps());
        assert_eq!(shown(reader.read_timestamp(2)), Ok(None));
    }

    #[test]
    fn zero_date_time_with_a_trailer_reads_it() {
        let data = ser_bytes(&test_header(2, 0), 2, &[T1, T1 + 7]);
        let reader = open_bytes(&data).unwrap();
        assert!(reader.has_timestamps());
        assert_eq!(
            shown(reader.read_timestamp(1)),
            Ok(Some(SerTimestamp::new(T1 + 7)))
        );
    }

    #[test]
    fn local_kind_date_time_still_reads_the_trailer() {
        let data = ser_bytes(&test_header(2, (T1 | 1 << 63) as i64), 2, &[T1, T1 + 1]);
        let reader = open_bytes(&data).unwrap();
        assert_eq!(reader.header().date_time, T1 as i64);
        assert_eq!(
            shown(reader.read_timestamp(1)),
            Ok(Some(SerTimestamp::new(T1 + 1)))
        );
    }

    #[test]
    fn trailer_ticks_have_the_kind_bits_masked() {
        let data = ser_bytes(&test_header(1, T1 as i64), 1, &[T1 | 1 << 62]);
        let reader = open_bytes(&data).unwrap();
        let ts = reader.read_timestamp(0).unwrap().unwrap();
        assert_eq!(ts.ticks(), T1);
    }

    #[test]
    fn zero_frames_has_no_trailer() {
        let reader = open_bytes(&test_header(0, T1 as i64).to_bytes()).unwrap();
        assert_eq!(reader.frame_count(), 0);
        assert!(!reader.has_timestamps());
        assert_eq!(reader.iter_frames().len(), 0);
    }

    #[test]
    fn frames_wider_than_65535_are_read() {
        let header = SerHeader {
            image_width: 70_000,
            image_height: 1,
            pixel_depth_per_plane: 8,
            frame_count: 1,
            ..SerHeader::default()
        };
        let data = ser_bytes(&header, 1, &[]);
        let reader = open_bytes(&data).unwrap();
        assert_eq!(shown(reader.frame_slice(0)), Ok(&data[HEADER_SIZE..]));
    }

    #[test]
    fn overflowing_header_is_rejected() {
        let header = SerHeader {
            image_width: 65535,
            image_height: 65535,
            color_id: ColorId::Rgb,
            frame_count: u32::MAX,
            ..SerHeader::default()
        };
        assert_eq!(
            shown(open_bytes(&header.to_bytes()).map(|_| ())),
            Err(
                "Invalid SER header: 4294967295 frames of 65535x65535 RGB at depth 16 \
                 overflow a 64-bit file size"
                    .to_string()
            )
        );
    }

    #[test]
    fn read_frame_returns_data_and_timestamp() {
        let data = create_test_ser_file();
        let reader = open_bytes(&data).unwrap();

        let frame = reader.read_frame(1).unwrap();
        assert_eq!(frame.id, 1);
        assert_eq!(frame.data, &data[378..578]);
        assert_eq!(frame.timestamp, Some(SerTimestamp::new(T1 + 1)));
        assert_eq!(
            shown(reader.read_frame_data(0)),
            Ok(data[178..378].to_vec())
        );
    }

    #[test]
    fn read_frame_u16_applies_the_byte_order() {
        let data = create_test_ser_file();
        let reader = open_bytes(&data).unwrap();
        let big: Vec<u16> = (0..100u16)
            .map(|i| u16::from_be_bytes([(2 * i % 256) as u8, (2 * i % 256 + 1) as u8]))
            .collect();
        assert_eq!(shown(reader.read_frame_u16(0)), Ok(big));

        let mut little_file = data;
        little_file[22..26].copy_from_slice(&0u32.to_le_bytes());
        let reader = open_bytes(&little_file).unwrap();
        let little: Vec<u16> = (0..100u16)
            .map(|i| u16::from_le_bytes([(2 * i % 256) as u8, (2 * i % 256 + 1) as u8]))
            .collect();
        assert_eq!(shown(reader.read_frame_u16(0)), Ok(little));
    }

    #[test]
    fn read_frame_u16_keeps_low_depth_samples_as_stored() {
        let header = SerHeader {
            pixel_depth_per_plane: 12,
            byte_order: ByteOrder::LittleEndian,
            image_width: 2,
            image_height: 1,
            frame_count: 1,
            ..SerHeader::default()
        };
        let mut data = header.to_bytes().to_vec();
        data.extend_from_slice(&[0xFF, 0x0F, 0x01, 0x00]);
        let reader = open_bytes(&data).unwrap();
        assert_eq!(shown(reader.read_frame_u16(0)), Ok(vec![4095, 1]));
        assert_eq!(reader.header().pixel_depth_per_plane, 12);
    }

    #[test]
    fn read_frame_u16_rejects_8_bit_frames() {
        let header = SerHeader {
            pixel_depth_per_plane: 8,
            frame_count: 1,
            ..SerHeader::default()
        };
        let mut data = header.to_bytes().to_vec();
        data.push(9);
        let reader = open_bytes(&data).unwrap();
        assert!(matches!(
            reader.read_frame_u16(0),
            Err(SerError::NotSixteenBit { depth: 8 })
        ));
        assert_eq!(shown(reader.frame_slice(0)), Ok(&[9u8][..]));
    }

    #[test]
    fn frame_slices_success() {
        let data = create_test_ser_file();
        let reader = open_bytes(&data).unwrap();
        let slices = reader.frame_slices(0, 2).unwrap();

        assert_eq!(slices.len(), 2);
        assert_eq!(slices[0], (0, &data[178..378]));
        assert_eq!(slices[1], (1, &data[378..578]));
    }

    #[test]
    fn frame_slices_clamped_to_available() {
        let data = create_test_ser_file();
        let reader = open_bytes(&data).unwrap();

        let slices = reader.frame_slices(1, 100).unwrap();
        assert_eq!(slices, vec![(1, &data[378..578])]);
        assert_eq!(
            shown(reader.frame_slices(1, u32::MAX).map(|s| s.len())),
            Ok(1)
        );
        assert_eq!(shown(reader.frame_slices(2, 5).map(|s| s.len())), Ok(0));
    }

    #[test]
    fn frame_slices_past_the_end_is_an_error() {
        let reader = open_bytes(&create_test_ser_file()).unwrap();
        assert_eq!(
            shown(reader.frame_slices(5, 1).map(|s| s.len())),
            Err("Frame index 5 out of bounds (total frames: 2)".to_string())
        );
    }

    #[test]
    fn iter_frames_reads_every_frame() {
        let data = create_test_ser_file();
        let reader = open_bytes(&data).unwrap();
        let mut iter = reader.iter_frames();
        assert_eq!(iter.len(), 2);

        let first = iter.next().unwrap().unwrap();
        assert_eq!((first.id, first.data), (0, &data[178..378]));
        assert_eq!(iter.len(), 1);
        let second = iter.next().unwrap().unwrap();
        assert_eq!(second.timestamp, Some(SerTimestamp::new(T1 + 1)));
        assert!(iter.next().is_none());
        assert_eq!(iter.len(), 0);
    }

    #[test]
    fn iter_frames_range() {
        let reader = open_bytes(&create_test_ser_file()).unwrap();
        let ids = |start, count| {
            shown(
                reader
                    .iter_frames_range(start, count)
                    .map(|iter| iter.map(|frame| frame.unwrap().id).collect::<Vec<_>>()),
            )
        };

        assert_eq!(ids(0, 1), Ok(vec![0]));
        assert_eq!(ids(1, 5), Ok(vec![1]));
        assert_eq!(ids(0, u32::MAX), Ok(vec![0, 1]));
        assert_eq!(ids(2, 1), Ok(vec![]));
        assert_eq!(
            shown(reader.iter_frames_range(1, 5).map(|i| i.len())),
            Ok(1)
        );
    }

    #[test]
    fn iter_frames_range_past_the_end_is_an_error() {
        let reader = open_bytes(&create_test_ser_file()).unwrap();
        for start in [3, 5, u32::MAX] {
            assert_eq!(
                shown(reader.iter_frames_range(start, 2).map(|i| i.len())),
                Err(format!(
                    "Frame index {start} out of bounds (total frames: 2)"
                ))
            );
        }
    }
}
