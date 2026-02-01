use crate::ser::{Result, SerError, SerFrameId, SerHeader};
use std::collections::VecDeque;
use std::sync::{Arc, Mutex};

pub struct FrameBuffer {
    frames: Arc<Mutex<VecDeque<Vec<u8>>>>,
    frame_size: u64,
    max_frames: usize,
    current_id: SerFrameId,
}

impl FrameBuffer {
    pub fn new(header: &SerHeader, buffer_frames: usize) -> Self {
        Self {
            frames: Arc::new(Mutex::new(VecDeque::with_capacity(buffer_frames))),
            frame_size: header.frame_size(),
            max_frames: buffer_frames,
            current_id: 0,
        }
    }

    pub fn push_frame(&mut self, data: Vec<u8>) -> Result<SerFrameId> {
        if data.len() != self.frame_size as usize {
            return Err(SerError::BufferSizeMismatch {
                expected: self.frame_size as usize,
                actual: data.len(),
            });
        }

        let mut frames = self.frames.lock().unwrap();
        if frames.len() >= self.max_frames {
            frames.pop_front();
        }

        frames.push_back(data);
        let id = self.current_id;
        self.current_id += 1;
        Ok(id)
    }

    pub fn get_frame(&self, id: SerFrameId) -> Option<Vec<u8>> {
        let frames = self.frames.lock().unwrap();
        let relative_id = if self.current_id >= frames.len() as u32 {
            self.current_id - frames.len() as u32
        } else {
            0
        };

        if id >= relative_id && id < self.current_id {
            let index = (id - relative_id) as usize;
            frames.get(index).cloned()
        } else {
            None
        }
    }

    pub fn available_frames(&self) -> usize {
        self.frames.lock().unwrap().len()
    }

    pub fn is_full(&self) -> bool {
        self.frames.lock().unwrap().len() >= self.max_frames
    }

    pub fn clear(&mut self) {
        self.frames.lock().unwrap().clear();
        self.current_id = 0;
    }

    pub fn with_capacity(header: &SerHeader, megabytes: usize) -> Self {
        let frame_size = header.frame_size() as usize;
        let max_frames = if frame_size > 0 {
            (megabytes * 1024 * 1024) / frame_size
        } else {
            1024
        };
        Self::new(header, max_frames.max(1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_header() -> SerHeader {
        SerHeader {
            image_width: 10,
            image_height: 10,
            pixel_depth_per_plane: 16,
            ..Default::default()
        }
    }

    #[test]
    fn new_creates_empty_buffer() {
        let header = test_header();
        let buffer = FrameBuffer::new(&header, 10);
        assert_eq!(buffer.available_frames(), 0);
        assert!(!buffer.is_full());
        assert_eq!(buffer.frame_size, 200);
        assert_eq!(buffer.max_frames, 10);
        assert_eq!(buffer.current_id, 0);
    }

    #[test]
    fn push_frame_adds_valid_frame() {
        let header = test_header();
        let mut buffer = FrameBuffer::new(&header, 2);
        let frame_data = vec![42u8; 200];

        let id = buffer.push_frame(frame_data.clone()).unwrap();
        assert_eq!(id, 0);
        assert_eq!(buffer.available_frames(), 1);
        assert!(!buffer.is_full());

        let retrieved = buffer.get_frame(0).unwrap();
        assert_eq!(retrieved, frame_data);
    }

    #[test]
    fn push_frame_rejects_wrong_size() {
        let header = test_header();
        let mut buffer = FrameBuffer::new(&header, 10);
        let wrong_size_data = vec![0u8; 100];

        let result = buffer.push_frame(wrong_size_data);
        assert!(matches!(
            result,
            Err(SerError::BufferSizeMismatch {
                expected: 200,
                actual: 100
            })
        ));
    }

    #[test]
    fn buffer_overflow_removes_oldest() {
        let header = test_header();
        let mut buffer = FrameBuffer::new(&header, 2);

        let frame1 = vec![1u8; 200];
        let frame2 = vec![2u8; 200];
        let frame3 = vec![3u8; 200];

        buffer.push_frame(frame1.clone()).unwrap();
        buffer.push_frame(frame2.clone()).unwrap();
        assert!(buffer.is_full());

        buffer.push_frame(frame3.clone()).unwrap();
        assert_eq!(buffer.available_frames(), 2);

        assert!(buffer.get_frame(0).is_none());
        assert_eq!(buffer.get_frame(1).unwrap(), frame2);
        assert_eq!(buffer.get_frame(2).unwrap(), frame3);
    }

    #[test]
    fn get_frame_returns_none_for_invalid_id() {
        let header = test_header();
        let mut buffer = FrameBuffer::new(&header, 5);
        buffer.push_frame(vec![42u8; 200]).unwrap();

        assert!(buffer.get_frame(999).is_none());
    }

    #[test]
    fn clear_empties_buffer() {
        let header = test_header();
        let mut buffer = FrameBuffer::new(&header, 5);
        buffer.push_frame(vec![0u8; 200]).unwrap();
        buffer.push_frame(vec![1u8; 200]).unwrap();

        assert_eq!(buffer.available_frames(), 2);
        buffer.clear();
        assert_eq!(buffer.available_frames(), 0);
        assert_eq!(buffer.current_id, 0);
    }

    #[test]
    fn with_capacity_calculates_frame_count() {
        let header = test_header();
        let buffer = FrameBuffer::with_capacity(&header, 1);

        let expected_frames = (1024 * 1024) / 200;
        assert_eq!(buffer.max_frames, expected_frames);
    }

    #[test]
    fn with_capacity_handles_zero_frame_size() {
        let header = SerHeader {
            image_width: 0,
            image_height: 0,
            ..Default::default()
        };

        let buffer = FrameBuffer::with_capacity(&header, 1);
        assert_eq!(buffer.max_frames, 1024);
    }

    #[test]
    fn frame_id_sequence() {
        let header = test_header();
        let mut buffer = FrameBuffer::new(&header, 10);

        let id1 = buffer.push_frame(vec![1u8; 200]).unwrap();
        let id2 = buffer.push_frame(vec![2u8; 200]).unwrap();
        let id3 = buffer.push_frame(vec![3u8; 200]).unwrap();

        assert_eq!(id1, 0);
        assert_eq!(id2, 1);
        assert_eq!(id3, 2);
    }

    #[test]
    fn concurrent_access_safety() {
        let header = test_header();
        let buffer = FrameBuffer::new(&header, 10);

        let buffer_clone = FrameBuffer {
            frames: buffer.frames.clone(),
            frame_size: buffer.frame_size,
            max_frames: buffer.max_frames,
            current_id: 0,
        };

        assert_eq!(buffer.available_frames(), buffer_clone.available_frames());
    }

    #[test]
    fn get_frame_early_access_pattern() {
        let header = test_header();
        let mut buffer = FrameBuffer::new(&header, 10);

        buffer.push_frame(vec![1u8; 200]).unwrap();
        buffer.push_frame(vec![2u8; 200]).unwrap();

        let frames_len = buffer.available_frames();
        assert_eq!(frames_len, 2);

        assert!(buffer.get_frame(0).is_some());
        assert!(buffer.get_frame(1).is_some());

        let mut special_buffer = FrameBuffer::new(&header, 10);
        special_buffer.push_frame(vec![42u8; 200]).unwrap();

        special_buffer.current_id = 0;

        let result = special_buffer.get_frame(0);
        assert!(result.is_none());
    }
}
