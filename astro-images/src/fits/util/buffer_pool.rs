use std::collections::VecDeque;

#[derive(Debug)]
pub struct BufferPool {
    pool: VecDeque<Vec<u8>>,
    max_capacity: usize,
    target_capacity: usize,
}

impl BufferPool {
    pub fn new(target_capacity: usize, max_capacity: usize) -> Self {
        Self {
            pool: VecDeque::new(),
            max_capacity,
            target_capacity,
        }
    }

    pub fn get(&mut self, size: usize) -> Vec<u8> {
        if let Some(mut buffer) = self.pool.pop_front() {
            if buffer.capacity() >= size {
                buffer.clear();
                buffer.resize(size, 0);
                return buffer;
            }
        }

        vec![0u8; size]
    }

    pub fn return_buffer(&mut self, mut buffer: Vec<u8>) {
        if self.pool.len() < self.max_capacity && buffer.capacity() >= self.target_capacity {
            buffer.clear();
            self.pool.push_back(buffer);
        }
    }

    pub fn clear(&mut self) {
        self.pool.clear();
    }
}

impl Default for BufferPool {
    fn default() -> Self {
        Self::new(128 * 1024, 50)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn buffer_pool_functionality() {
        let mut pool = BufferPool::new(1024, 5);

        let buffer1 = pool.get(1024);
        assert_eq!(buffer1.len(), 1024);
        assert!(buffer1.iter().all(|&x| x == 0));

        pool.return_buffer(buffer1);
        assert_eq!(pool.pool.len(), 1);

        let buffer2 = pool.get(512);
        assert_eq!(buffer2.len(), 512);
        assert_eq!(pool.pool.len(), 0);
    }

    #[test]
    fn buffer_pool_capacity_filtering() {
        let mut pool = BufferPool::new(1024, 2);

        let small_buffer = vec![0u8; 100];
        pool.return_buffer(small_buffer);
        assert_eq!(pool.pool.len(), 0);

        let large_buffer = vec![0u8; 2048];
        pool.return_buffer(large_buffer);
        assert_eq!(pool.pool.len(), 1);
    }

    #[test]
    fn buffer_pool_max_capacity() {
        let mut pool = BufferPool::new(512, 2);

        pool.return_buffer(vec![0u8; 1024]);
        pool.return_buffer(vec![0u8; 1024]);
        assert_eq!(pool.pool.len(), 2);

        pool.return_buffer(vec![0u8; 1024]);
        assert_eq!(pool.pool.len(), 2);
    }

    #[test]
    fn buffer_pool_clear() {
        let mut pool = BufferPool::new(512, 5);

        pool.return_buffer(vec![0u8; 1024]);
        pool.return_buffer(vec![0u8; 1024]);
        assert_eq!(pool.pool.len(), 2);

        pool.clear();
        assert_eq!(pool.pool.len(), 0);
    }

    #[test]
    fn buffer_pool_default_optimized_for_astronomy() {
        let mut pool = BufferPool::default();

        let large_buffer = pool.get(128 * 1024);
        assert_eq!(large_buffer.len(), 128 * 1024);
        assert!(large_buffer.iter().all(|&x| x == 0));

        pool.return_buffer(large_buffer);
        assert_eq!(pool.pool.len(), 1);

        let reused_buffer = pool.get(64 * 1024);
        assert_eq!(reused_buffer.len(), 64 * 1024);
        assert_eq!(pool.pool.len(), 0);
    }
}
