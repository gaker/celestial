use super::record::EopRecord;
use crate::{CoordError, CoordResult};
use std::collections::HashMap;
use std::marker::PhantomData;
use std::ptr::NonNull;

struct LruNode {
    key: i32,
    record: EopRecord,
    prev: Option<NonNull<LruNode>>,
    next: Option<NonNull<LruNode>>,
}

pub struct EopMemoryCache {
    cache: HashMap<i32, NonNull<LruNode>>,

    head: Option<NonNull<LruNode>>,

    tail: Option<NonNull<LruNode>>,

    max_size: usize,

    size: usize,

    _marker: PhantomData<Box<LruNode>>,
}

impl EopMemoryCache {
    pub fn new(max_size: usize) -> Self {
        Self {
            cache: HashMap::with_capacity(max_size),
            head: None,
            tail: None,
            max_size,
            size: 0,
            _marker: PhantomData,
        }
    }

    pub fn get(&mut self, mjd: f64) -> Option<EopRecord> {
        let key = mjd.floor() as i32;

        if let Some(&node_ptr) = self.cache.get(&key) {
            unsafe {
                let record = (*node_ptr.as_ptr()).record.clone();
                self.move_to_front(node_ptr);
                Some(record)
            }
        } else {
            None
        }
    }

    pub fn insert(&mut self, record: EopRecord) {
        let key = record.mjd.floor() as i32;

        if let Some(&existing_node) = self.cache.get(&key) {
            unsafe {
                (*existing_node.as_ptr()).record = record;
                self.move_to_front(existing_node);
            }
            return;
        }

        let new_node = Box::into_raw(Box::new(LruNode {
            key,
            record,
            prev: None,
            next: None,
        }));
        let new_node_ptr = unsafe { NonNull::new_unchecked(new_node) };

        self.cache.insert(key, new_node_ptr);

        unsafe {
            self.add_to_front(new_node_ptr);
        }
        self.size += 1;

        if self.size > self.max_size {
            self.evict_lru();
        }
    }

    unsafe fn move_to_front(&mut self, node_ptr: NonNull<LruNode>) {
        self.remove_from_list(node_ptr);
        self.add_to_front(node_ptr);
    }

    unsafe fn add_to_front(&mut self, node_ptr: NonNull<LruNode>) {
        let node = node_ptr.as_ptr();
        (*node).prev = None;
        (*node).next = self.head;

        if let Some(old_head) = self.head {
            (*old_head.as_ptr()).prev = Some(node_ptr);
        } else {
            self.tail = Some(node_ptr);
        }

        self.head = Some(node_ptr);
    }

    unsafe fn remove_from_list(&mut self, node_ptr: NonNull<LruNode>) {
        let node = node_ptr.as_ptr();

        match ((*node).prev, (*node).next) {
            (None, None) => {
                self.head = None;
                self.tail = None;
            }
            (None, Some(next)) => {
                self.head = Some(next);
                (*next.as_ptr()).prev = None;
            }
            (Some(prev), None) => {
                self.tail = Some(prev);
                (*prev.as_ptr()).next = None;
            }
            (Some(prev), Some(next)) => {
                (*prev.as_ptr()).next = Some(next);
                (*next.as_ptr()).prev = Some(prev);
            }
        }
    }

    fn evict_lru(&mut self) {
        if let Some(tail_ptr) = self.tail {
            unsafe {
                let key = (*tail_ptr.as_ptr()).key;
                self.cache.remove(&key);
                self.remove_from_list(tail_ptr);
                let _ = Box::from_raw(tail_ptr.as_ptr());
                self.size -= 1;
            }
        }
    }

    pub fn clear(&mut self) {
        while let Some(head_ptr) = self.head {
            unsafe {
                self.remove_from_list(head_ptr);
                let _ = Box::from_raw(head_ptr.as_ptr());
            }
        }

        self.cache.clear();
        self.head = None;
        self.tail = None;
        self.size = 0;
    }

    pub fn len(&self) -> usize {
        self.size
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    pub fn hit_rate(&self) -> f64 {
        if self.max_size == 0 {
            0.0
        } else {
            self.size as f64 / self.max_size as f64
        }
    }
}

impl Drop for EopMemoryCache {
    fn drop(&mut self) {
        self.clear();
    }
}

pub struct EopDiskCache {
    cache_dir: std::path::PathBuf,

    use_compression: bool,
}

impl EopDiskCache {
    pub fn new<P: AsRef<std::path::Path>>(cache_dir: P) -> CoordResult<Self> {
        let cache_dir = cache_dir.as_ref().to_path_buf();

        if !cache_dir.exists() {
            std::fs::create_dir_all(&cache_dir).map_err(|e| {
                CoordError::external_library("Failed to create cache directory", &e.to_string())
            })?;
        }

        Ok(Self {
            cache_dir,
            use_compression: true,
        })
    }

    /// Configure compression (currently unimplemented).
    ///
    /// **Note**: This only changes the file extension to `.eop.zst` vs `.eop`.
    /// Actual compression is not implemented. To add compression, integrate a compression
    /// library (e.g., `zstd` crate) in the `store()` and `load()` methods.
    pub fn with_compression(mut self, use_compression: bool) -> Self {
        self.use_compression = use_compression;
        self
    }

    fn cache_file_path(&self, data_type: &str, start_mjd: f64, end_mjd: f64) -> std::path::PathBuf {
        let filename = if self.use_compression {
            format!("{}_{:.0}_{:.0}.eop.zst", data_type, start_mjd, end_mjd)
        } else {
            format!("{}_{:.0}_{:.0}.eop", data_type, start_mjd, end_mjd)
        };

        self.cache_dir.join(filename)
    }

    pub fn store(&self, data_type: &str, records: &[EopRecord]) -> CoordResult<()> {
        if records.is_empty() {
            return Ok(());
        }

        let start_mjd = records[0].mjd;
        let end_mjd = records.last().unwrap().mjd;
        let cache_file = self.cache_file_path(data_type, start_mjd, end_mjd);

        let binary_data = self.encode_eop_records(records)?;

        // Note: Actual compression not implemented yet (see with_compression() docs)
        std::fs::write(&cache_file, binary_data).map_err(|e| {
            CoordError::external_library("Failed to write cache file", &e.to_string())
        })?;

        Ok(())
    }

    fn encode_eop_records(&self, records: &[EopRecord]) -> CoordResult<Vec<u8>> {
        let mut buffer = Vec::new();

        // EOP2: Version 2 uses i32 for MJD delta (was i16 in EOP1)
        // Note: Only stores dx/dy presence flags, not full EopFlags (source/quality metadata).
        // This is a deliberate tradeoff to keep the binary format compact.
        buffer.extend_from_slice(b"EOP2");
        buffer.extend_from_slice(&(records.len() as u32).to_le_bytes());

        // Reserve 4 bytes for CRC32 checksum (written at end)
        let crc_offset = buffer.len();
        buffer.extend_from_slice(&[0u8; 4]);

        if records.is_empty() {
            return Ok(buffer);
        }

        let first = &records[0];
        buffer.extend_from_slice(&first.mjd.to_le_bytes());
        buffer.extend_from_slice(&first.x_p_encoded.to_le_bytes());
        buffer.extend_from_slice(&first.y_p_encoded.to_le_bytes());
        buffer.extend_from_slice(&first.ut1_utc_encoded.to_le_bytes());
        buffer.extend_from_slice(&first.lod_encoded.to_le_bytes());

        let mut flags = 0u8;
        if first.dx_encoded.is_some() {
            flags |= 0x01;
        }
        if first.dy_encoded.is_some() {
            flags |= 0x02;
        }
        buffer.push(flags);

        if let Some(dx) = first.dx_encoded {
            buffer.extend_from_slice(&dx.to_le_bytes());
        }
        if let Some(dy) = first.dy_encoded {
            buffer.extend_from_slice(&dy.to_le_bytes());
        }

        for i in 1..records.len() {
            let current = &records[i];
            let previous = &records[i - 1];

            // Use i32 for MJD delta to handle gaps > 32 days (i16 max = 32,767 ms = 32.767 days)
            let mjd_delta = ((current.mjd - previous.mjd) * 1000.0).round() as i32;
            let x_p_delta = (current.x_p_encoded - previous.x_p_encoded) as i16;
            let y_p_delta = (current.y_p_encoded - previous.y_p_encoded) as i16;
            let ut1_utc_delta = (current.ut1_utc_encoded - previous.ut1_utc_encoded) as i16;
            let lod_delta = (current.lod_encoded - previous.lod_encoded) as i16;

            buffer.extend_from_slice(&mjd_delta.to_le_bytes());
            buffer.extend_from_slice(&x_p_delta.to_le_bytes());
            buffer.extend_from_slice(&y_p_delta.to_le_bytes());
            buffer.extend_from_slice(&ut1_utc_delta.to_le_bytes());
            buffer.extend_from_slice(&lod_delta.to_le_bytes());

            let mut flags = 0u8;
            if current.dx_encoded.is_some() {
                flags |= 0x01;
            }
            if current.dy_encoded.is_some() {
                flags |= 0x02;
            }
            buffer.push(flags);

            if let Some(dx) = current.dx_encoded {
                buffer.extend_from_slice(&dx.to_le_bytes());
            }
            if let Some(dy) = current.dy_encoded {
                buffer.extend_from_slice(&dy.to_le_bytes());
            }
        }

        // Compute CRC32 checksum of data after the checksum field
        let data_for_checksum = &buffer[crc_offset + 4..];
        let mut hasher = crc32fast::Hasher::new();
        hasher.update(data_for_checksum);
        let checksum = hasher.finalize();

        // Write checksum to reserved location
        buffer[crc_offset..crc_offset + 4].copy_from_slice(&checksum.to_le_bytes());

        Ok(buffer)
    }

    pub fn load(
        &self,
        data_type: &str,
        start_mjd: f64,
        end_mjd: f64,
    ) -> CoordResult<Option<Vec<EopRecord>>> {
        let cache_file = self.cache_file_path(data_type, start_mjd, end_mjd);

        if !cache_file.exists() {
            return Ok(None);
        }

        let file_data = std::fs::read(&cache_file).map_err(|e| {
            CoordError::external_library("Failed to read cache file", &e.to_string())
        })?;

        let records = self.decode_eop_records(&file_data)?;
        Ok(Some(records))
    }

    fn validate_header(data: &[u8]) -> CoordResult<(u8, usize)> {
        if data.len() < 8 {
            return Err(CoordError::external_library(
                "Invalid cache file",
                "File too small to contain valid header",
            ));
        }

        let version = if &data[0..4] == b"EOP2" {
            2
        } else if &data[0..4] == b"EOP1" {
            1
        } else {
            return Err(CoordError::external_library(
                "Invalid cache file",
                "Wrong magic bytes or unsupported version",
            ));
        };

        let record_count = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
        Ok((version, record_count))
    }

    fn validate_crc32(data: &[u8]) -> CoordResult<()> {
        if data.len() < 12 {
            return Err(CoordError::external_library(
                "Invalid cache file",
                "EOP2 file too small for checksum",
            ));
        }

        let stored_checksum = u32::from_le_bytes([data[8], data[9], data[10], data[11]]);
        let mut hasher = crc32fast::Hasher::new();
        hasher.update(&data[12..]);
        let computed_checksum = hasher.finalize();

        if stored_checksum != computed_checksum {
            return Err(CoordError::external_library(
                "Cache file corrupted",
                &format!(
                    "CRC32 mismatch: expected 0x{:08X}, got 0x{:08X}",
                    stored_checksum, computed_checksum
                ),
            ));
        }

        Ok(())
    }

    fn decode_eop_records(&self, data: &[u8]) -> CoordResult<Vec<EopRecord>> {
        let (version, record_count) = Self::validate_header(data)?;
        let mut records = Vec::with_capacity(record_count);

        let mut offset = 8;
        if version == 2 {
            Self::validate_crc32(data)?;
            offset = 12;
        }

        if record_count == 0 {
            return Ok(records);
        }

        if offset + 24 > data.len() {
            return Err(CoordError::external_library(
                "Invalid cache file",
                "Insufficient data for first record",
            ));
        }

        let mjd = f64::from_le_bytes(data[offset..offset + 8].try_into().unwrap());
        offset += 8;
        let x_p_encoded = i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
        offset += 4;
        let y_p_encoded = i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
        offset += 4;
        let ut1_utc_encoded = i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
        offset += 4;
        let lod_encoded = i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
        offset += 4;

        let flags = data[offset];
        offset += 1;

        let mut dx_encoded = None;
        let mut dy_encoded = None;

        if flags & 0x01 != 0 {
            if offset + 2 > data.len() {
                return Err(CoordError::external_library(
                    "Invalid cache file",
                    "Missing dX data",
                ));
            }
            dx_encoded = Some(i16::from_le_bytes(
                data[offset..offset + 2].try_into().unwrap(),
            ));
            offset += 2;
        }

        if flags & 0x02 != 0 {
            if offset + 2 > data.len() {
                return Err(CoordError::external_library(
                    "Invalid cache file",
                    "Missing dY data",
                ));
            }
            dy_encoded = Some(i16::from_le_bytes(
                data[offset..offset + 2].try_into().unwrap(),
            ));
            offset += 2;
        }

        let first_record = EopRecord {
            mjd,
            x_p_encoded,
            y_p_encoded,
            ut1_utc_encoded,
            lod_encoded,
            dx_encoded,
            dy_encoded,
            flags: super::record::EopFlags::default(),
        };
        records.push(first_record);

        for _ in 1..record_count {
            // Version 2: 13 bytes (i32 MJD + 4×i16 + flags)
            // Version 1: 11 bytes (i16 MJD + 4×i16 + flags)
            let record_size = if version == 2 { 13 } else { 11 };
            if offset + record_size > data.len() {
                return Err(CoordError::external_library(
                    "Invalid cache file",
                    "Insufficient data for delta record",
                ));
            }

            // Read MJD delta based on version
            let mjd_delta = if version == 2 {
                // EOP2: i32 for MJD delta (handles gaps > 32 days)
                let delta = i32::from_le_bytes(data[offset..offset + 4].try_into().unwrap());
                offset += 4;
                delta
            } else {
                // EOP1: i16 for MJD delta (legacy, max 32.767 days)
                let delta = i16::from_le_bytes(data[offset..offset + 2].try_into().unwrap());
                offset += 2;
                delta as i32
            };
            let x_p_delta = i16::from_le_bytes(data[offset..offset + 2].try_into().unwrap());
            offset += 2;
            let y_p_delta = i16::from_le_bytes(data[offset..offset + 2].try_into().unwrap());
            offset += 2;
            let ut1_utc_delta = i16::from_le_bytes(data[offset..offset + 2].try_into().unwrap());
            offset += 2;
            let lod_delta = i16::from_le_bytes(data[offset..offset + 2].try_into().unwrap());
            offset += 2;

            let flags = data[offset];
            offset += 1;

            let previous = &records.last().unwrap();

            let mjd = previous.mjd + (mjd_delta as f64) / 1000.0;
            let x_p_encoded = previous.x_p_encoded + x_p_delta as i32;
            let y_p_encoded = previous.y_p_encoded + y_p_delta as i32;
            let ut1_utc_encoded = previous.ut1_utc_encoded + ut1_utc_delta as i32;
            let lod_encoded = previous.lod_encoded + lod_delta as i32;

            let mut dx_encoded = None;
            let mut dy_encoded = None;

            if flags & 0x01 != 0 {
                if offset + 2 > data.len() {
                    return Err(CoordError::external_library(
                        "Invalid cache file",
                        "Missing dX data",
                    ));
                }
                dx_encoded = Some(i16::from_le_bytes(
                    data[offset..offset + 2].try_into().unwrap(),
                ));
                offset += 2;
            }

            if flags & 0x02 != 0 {
                if offset + 2 > data.len() {
                    return Err(CoordError::external_library(
                        "Invalid cache file",
                        "Missing dY data",
                    ));
                }
                dy_encoded = Some(i16::from_le_bytes(
                    data[offset..offset + 2].try_into().unwrap(),
                ));
                offset += 2;
            }

            let flags = super::record::EopFlags {
                has_polar_motion: true,
                has_ut1_utc: true,
                has_cip_offsets: dx_encoded.is_some() && dy_encoded.is_some(),
                ..super::record::EopFlags::default()
            };

            let record = EopRecord {
                mjd,
                x_p_encoded,
                y_p_encoded,
                ut1_utc_encoded,
                lod_encoded,
                dx_encoded,
                dy_encoded,
                flags,
            };
            records.push(record);
        }

        Ok(records)
    }

    pub fn has_data(&self, data_type: &str, start_mjd: f64, end_mjd: f64) -> bool {
        let cache_file = self.cache_file_path(data_type, start_mjd, end_mjd);
        cache_file.exists()
    }

    pub fn clear_all(&self) -> CoordResult<()> {
        let entries = std::fs::read_dir(&self.cache_dir).map_err(|e| {
            CoordError::external_library("Failed to read cache directory", &e.to_string())
        })?;

        for entry in entries {
            let entry = entry.map_err(|e| {
                CoordError::external_library("Failed to read directory entry", &e.to_string())
            })?;

            let path = entry.path();
            if path.is_file()
                && path
                    .extension()
                    .is_some_and(|ext| ext == "eop" || ext == "zst" || ext == "json")
            {
                std::fs::remove_file(&path).map_err(|e| {
                    CoordError::external_library("Failed to delete cache file", &e.to_string())
                })?;
            }
        }

        Ok(())
    }

    pub fn cache_size(&self) -> CoordResult<u64> {
        let mut total_size = 0u64;

        let entries = std::fs::read_dir(&self.cache_dir).map_err(|e| {
            CoordError::external_library("Failed to read cache directory", &e.to_string())
        })?;

        for entry in entries {
            let entry = entry.map_err(|e| {
                CoordError::external_library("Failed to read directory entry", &e.to_string())
            })?;

            let metadata = entry.metadata().map_err(|e| {
                CoordError::external_library("Failed to read file metadata", &e.to_string())
            })?;

            if metadata.is_file() {
                total_size += metadata.len();
            }
        }

        Ok(total_size)
    }

    pub fn get_for_mjd(&self, mjd: f64) -> Option<EopRecord> {
        let entries = std::fs::read_dir(&self.cache_dir).ok()?;

        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }

            if let Some(filename) = path.file_name().and_then(|n| n.to_str()) {
                if let Some((start, end)) = Self::parse_cache_filename(filename) {
                    if mjd >= start && mjd <= end {
                        if let Ok(Some(records)) = self.load("", start, end) {
                            for record in records {
                                if (record.mjd - mjd).abs() < 0.5 {
                                    return Some(record);
                                }
                            }
                        }
                    }
                }
            }
        }
        None
    }

    fn parse_cache_filename(filename: &str) -> Option<(f64, f64)> {
        let parts: Vec<&str> = filename.split('_').collect();
        if parts.len() >= 3 {
            let start = parts[1].parse::<f64>().ok()?;
            let end_part = parts[2].split('.').next()?;
            let end = end_part.parse::<f64>().ok()?;
            Some((start, end))
        } else {
            None
        }
    }
}

pub struct EopCacheSystem {
    memory_cache: EopMemoryCache,

    disk_cache: Option<EopDiskCache>,

    hits: u64,
    misses: u64,
}

impl EopCacheSystem {
    pub fn memory_only(memory_size: usize) -> Self {
        Self {
            memory_cache: EopMemoryCache::new(memory_size),
            disk_cache: None,
            hits: 0,
            misses: 0,
        }
    }

    pub fn with_disk_cache<P: AsRef<std::path::Path>>(
        memory_size: usize,
        cache_dir: P,
    ) -> CoordResult<Self> {
        Ok(Self {
            memory_cache: EopMemoryCache::new(memory_size),
            disk_cache: Some(EopDiskCache::new(cache_dir)?),
            hits: 0,
            misses: 0,
        })
    }

    pub fn get(&mut self, mjd: f64) -> Option<EopRecord> {
        if let Some(record) = self.memory_cache.get(mjd) {
            self.hits += 1;
            return Some(record);
        }

        if let Some(ref disk_cache) = self.disk_cache {
            if let Some(record) = disk_cache.get_for_mjd(mjd) {
                self.memory_cache.insert(record.clone());
                self.hits += 1;
                return Some(record);
            }
        }

        self.misses += 1;
        None
    }

    pub fn insert(&mut self, record: EopRecord) {
        self.memory_cache.insert(record);
    }

    pub fn insert_batch(&mut self, records: Vec<EopRecord>) {
        for record in records {
            self.memory_cache.insert(record);
        }
    }

    pub fn hit_ratio(&self) -> f64 {
        let total = self.hits + self.misses;
        if total == 0 {
            0.0
        } else {
            self.hits as f64 / total as f64
        }
    }

    pub fn clear(&mut self) -> CoordResult<()> {
        self.memory_cache.clear();

        if let Some(ref disk_cache) = self.disk_cache {
            disk_cache.clear_all()?;
        }

        self.hits = 0;
        self.misses = 0;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_record(mjd: f64) -> EopRecord {
        EopRecord::new(mjd, 0.1, 0.2, 0.01, 0.001).unwrap()
    }

    fn unique_cache_path(prefix: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        let unique = format!(
            "{}_{}_{}",
            prefix,
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        );
        path.push(unique);
        path
    }

    #[test]
    fn test_memory_cache() {
        let mut cache = EopMemoryCache::new(3);

        let record1 = create_test_record(59945.0);
        let record2 = create_test_record(59946.0);

        cache.insert(record1);
        cache.insert(record2);

        assert_eq!(cache.len(), 2);
        assert!(cache.get(59945.0).is_some());
        assert!(cache.get(59946.0).is_some());
        assert!(cache.get(59947.0).is_none());
    }

    #[test]
    fn test_memory_cache_lru_eviction() {
        let mut cache = EopMemoryCache::new(2);

        let record1 = create_test_record(59945.0);
        let record2 = create_test_record(59946.0);
        let record3 = create_test_record(59947.0);

        cache.insert(record1);
        cache.insert(record2);
        assert_eq!(cache.len(), 2);

        cache.insert(record3);
        assert_eq!(cache.len(), 2);
        assert!(cache.get(59945.0).is_none());
        assert!(cache.get(59946.0).is_some());
        assert!(cache.get(59947.0).is_some());
    }

    #[test]
    fn test_cache_system_memory_only() {
        let mut cache = EopCacheSystem::memory_only(10);

        let record = create_test_record(59945.0);
        cache.insert(record);

        assert!(cache.get(59945.0).is_some());
        assert!(cache.get(59946.0).is_none());

        assert!(cache.hit_ratio() > 0.0);
    }

    #[test]
    fn test_encode_decode_large_mjd_gaps() {
        // Test that we can handle gaps > 32 days (i16 overflow threshold)
        let cache_path = std::env::temp_dir().join("test_eop_large_gaps.cache");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let records = [
            create_test_record(59945.0), // Day 1
            create_test_record(59946.0), // Day 2 (+1 day)
            create_test_record(59980.0), // Day 36 (+34 days - would overflow i16!)
            create_test_record(60045.0), // Day 101 (+65 days - definitely overflow!)
            create_test_record(60046.0), // Day 102 (+1 day)
        ];

        // Encode
        let encoded = cache.encode_eop_records(&records).unwrap();

        // Verify EOP2 format
        assert_eq!(&encoded[0..4], b"EOP2", "Should use EOP2 format");

        // Verify CRC32 checksum is present
        let stored_checksum =
            u32::from_le_bytes([encoded[8], encoded[9], encoded[10], encoded[11]]);
        assert!(stored_checksum != 0, "CRC32 should be non-zero");

        // Decode
        let decoded = cache.decode_eop_records(&encoded).unwrap();

        // Verify all records decoded correctly
        assert_eq!(decoded.len(), records.len());
        for (i, (original, decoded)) in records.iter().zip(decoded.iter()).enumerate() {
            assert!(
                (original.mjd - decoded.mjd).abs() < 0.001,
                "MJD mismatch at record {}: expected {}, got {}",
                i,
                original.mjd,
                decoded.mjd
            );
        }

        // Verify large gaps specifically
        assert!(
            (decoded[2].mjd - decoded[1].mjd - 34.0).abs() < 0.001,
            "34-day gap failed"
        );
        assert!(
            (decoded[3].mjd - decoded[2].mjd - 65.0).abs() < 0.001,
            "65-day gap failed"
        );
    }

    #[test]
    fn test_crc32_corruption_detection() {
        let cache_path = std::env::temp_dir().join("test_eop_crc32.cache");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let records = [
            create_test_record(59945.0),
            create_test_record(59946.0),
            create_test_record(59947.0),
        ];

        let mut encoded = cache.encode_eop_records(&records).unwrap();

        // Corrupt a byte in the data section (after checksum)
        encoded[20] ^= 0xFF;

        // Decode should fail due to checksum mismatch
        let result = cache.decode_eop_records(&encoded);
        assert!(result.is_err(), "Should detect corruption");

        if let Err(e) = result {
            let error_msg = format!("{:?}", e);
            assert!(error_msg.contains("CRC32 mismatch") || error_msg.contains("corrupted"));
        }
    }

    #[test]
    fn test_decode_legacy_eop1_format() {
        // Test that we can still read old EOP1 format (with i16 MJD delta limitation)
        let cache_path = std::env::temp_dir().join("test_eop_legacy.cache");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        // Create records with gaps that fit in i16 (< 32.767 days)
        let records = [
            create_test_record(59945.0),
            create_test_record(59946.0), // +1 day
            create_test_record(59970.0), // +24 days (safe for i16)
        ];

        // Manually create EOP1 format
        let mut buffer = Vec::new();
        buffer.extend_from_slice(b"EOP1");
        buffer.extend_from_slice(&(records.len() as u32).to_le_bytes());

        // First record (full)
        buffer.extend_from_slice(&records[0].mjd.to_le_bytes());
        buffer.extend_from_slice(&records[0].x_p_encoded.to_le_bytes());
        buffer.extend_from_slice(&records[0].y_p_encoded.to_le_bytes());
        buffer.extend_from_slice(&records[0].ut1_utc_encoded.to_le_bytes());
        buffer.extend_from_slice(&records[0].lod_encoded.to_le_bytes());
        buffer.push(0u8); // flags

        // Delta records (using i16 for MJD)
        for i in 1..records.len() {
            let mjd_delta = ((records[i].mjd - records[i - 1].mjd) * 1000.0).round() as i16;
            buffer.extend_from_slice(&mjd_delta.to_le_bytes());
            buffer.extend_from_slice(&0i16.to_le_bytes()); // x_p_delta
            buffer.extend_from_slice(&0i16.to_le_bytes()); // y_p_delta
            buffer.extend_from_slice(&0i16.to_le_bytes()); // ut1_utc_delta
            buffer.extend_from_slice(&0i16.to_le_bytes()); // lod_delta
            buffer.push(0u8); // flags
        }

        // Decode EOP1 format
        let decoded = cache.decode_eop_records(&buffer).unwrap();

        // Verify decoding works
        assert_eq!(decoded.len(), records.len());
        for (i, (original, decoded)) in records.iter().zip(decoded.iter()).enumerate() {
            assert!(
                (original.mjd - decoded.mjd).abs() < 0.001,
                "MJD mismatch at record {}: expected {}, got {}",
                i,
                original.mjd,
                decoded.mjd
            );
        }
    }

    #[test]
    fn test_mjd_delta_range() {
        // Test the actual limits of i32 vs i16
        // i16 max: 32,767 ms = 32.767 days
        // i32 max: 2,147,483,647 ms = 2,147,483.647 days (~ 5,881 years!)

        let cache_path = std::env::temp_dir().join("test_eop_limits.cache");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        // Test near i16 boundary
        let records_near_limit = [
            create_test_record(59945.0),
            create_test_record(59945.0 + 32.5), // 32.5 days - just over i16 limit
        ];

        let encoded = cache.encode_eop_records(&records_near_limit).unwrap();
        let decoded = cache.decode_eop_records(&encoded).unwrap();

        assert_eq!(decoded.len(), 2);
        assert!((decoded[1].mjd - decoded[0].mjd - 32.5).abs() < 0.001);

        // Test very large gap (365 days)
        let records_year_gap = [
            create_test_record(59945.0),
            create_test_record(60310.0), // +365 days
        ];

        let encoded = cache.encode_eop_records(&records_year_gap).unwrap();
        let decoded = cache.decode_eop_records(&encoded).unwrap();

        assert_eq!(decoded.len(), 2);
        assert!((decoded[1].mjd - decoded[0].mjd - 365.0).abs() < 0.001);
    }

    #[test]
    fn test_memory_cache_is_empty() {
        let mut cache = EopMemoryCache::new(3);
        assert!(cache.is_empty());

        let record = create_test_record(59945.0);
        cache.insert(record);
        assert!(!cache.is_empty());

        cache.clear();
        assert!(cache.is_empty());
    }

    #[test]
    fn test_memory_cache_hit_rate() {
        let cache = EopMemoryCache::new(10);
        assert_eq!(cache.hit_rate(), 0.0);

        let mut cache = EopMemoryCache::new(10);
        cache.insert(create_test_record(59945.0));
        cache.insert(create_test_record(59946.0));
        cache.insert(create_test_record(59947.0));

        // 3 items in cache of max_size 10
        assert!((cache.hit_rate() - 0.3).abs() < 0.01);
    }

    #[test]
    fn test_memory_cache_middle_node_removal() {
        let mut cache = EopMemoryCache::new(5);

        // Insert 5 items to fill cache
        cache.insert(create_test_record(59945.0));
        cache.insert(create_test_record(59946.0));
        cache.insert(create_test_record(59947.0));
        cache.insert(create_test_record(59948.0));
        cache.insert(create_test_record(59949.0));

        // Access middle item to move it to front
        cache.get(59947.0);

        // Insert new item - should evict LRU (59945)
        cache.insert(create_test_record(59950.0));

        assert!(cache.get(59945.0).is_none());
        assert!(cache.get(59947.0).is_some());
        assert!(cache.get(59950.0).is_some());
    }

    #[test]
    fn test_disk_cache_with_compression() {
        let cache_path = std::env::temp_dir().join("test_eop_compression");
        let cache = EopDiskCache::new(&cache_path)
            .unwrap()
            .with_compression(true);

        let records = vec![create_test_record(59945.0), create_test_record(59946.0)];

        cache.store("test", &records).unwrap();
        let loaded = cache.load("test", 59945.0, 59946.0).unwrap();

        assert!(loaded.is_some());
        assert_eq!(loaded.unwrap().len(), 2);
    }

    #[test]
    fn test_disk_cache_without_compression() {
        let cache_path = std::env::temp_dir().join("test_eop_no_compression");
        let cache = EopDiskCache::new(&cache_path)
            .unwrap()
            .with_compression(false);

        let records = vec![create_test_record(59945.0), create_test_record(59946.0)];

        cache.store("test", &records).unwrap();
        let loaded = cache.load("test", 59945.0, 59946.0).unwrap();

        assert!(loaded.is_some());
        assert_eq!(loaded.unwrap().len(), 2);
    }

    #[test]
    fn test_disk_cache_store_empty_records() {
        let cache_path = unique_cache_path("eop_empty_store");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        cache.store("test", &[]).unwrap();
        let file_count = std::fs::read_dir(&cache_path)
            .map(|entries| entries.count())
            .unwrap_or_default();
        assert_eq!(file_count, 0);

        let _ = std::fs::remove_dir_all(&cache_path);
    }

    #[test]
    fn test_disk_cache_get_for_mjd_hits_disk() {
        let cache_path = unique_cache_path("eop_get_for_mjd");
        let mut system = EopCacheSystem::with_disk_cache(1, &cache_path).unwrap();

        let records = vec![create_test_record(59945.0), create_test_record(59946.0)];

        {
            let disk_cache = system
                .disk_cache
                .as_ref()
                .expect("disk cache not configured");
            disk_cache.store("", &records).unwrap();
            assert!(disk_cache.has_data("", 59945.0, 59946.0));
        }

        let size = system.disk_cache.as_ref().unwrap().cache_size().unwrap();
        assert!(size > 0);

        let retrieved = system.get(59945.0).unwrap();
        assert!((retrieved.mjd - 59945.0).abs() < 1e-9);
        assert_eq!(system.memory_cache.len(), 1);

        let disk_cache = system.disk_cache.as_ref().unwrap();
        disk_cache.clear_all().unwrap();
        assert_eq!(disk_cache.cache_size().unwrap(), 0);

        let _ = std::fs::remove_dir_all(&cache_path);
    }

    #[test]
    fn test_parse_cache_filename_variants() {
        let parsed =
            EopDiskCache::parse_cache_filename("test_59945_59946.eop.zst").expect("valid filename");
        assert_eq!(parsed.0, 59945.0);
        assert_eq!(parsed.1, 59946.0);

        assert!(EopDiskCache::parse_cache_filename("invalid").is_none());
        assert!(EopDiskCache::parse_cache_filename("missing_parts_123").is_none());
    }

    #[test]
    fn test_cache_system_hit_ratio() {
        let mut cache = EopCacheSystem::memory_only(10);

        // No hits/misses yet
        assert_eq!(cache.hit_ratio(), 0.0);

        let record = create_test_record(59945.0);
        cache.insert(record);

        // First access is a hit
        cache.get(59945.0);

        // Try to get non-existent record (miss)
        cache.get(59999.0);

        // Hit ratio should be 0.5 (1 hit, 1 miss)
        assert!((cache.hit_ratio() - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_memory_cache_zero_size() {
        let cache = EopMemoryCache::new(0);
        assert_eq!(cache.hit_rate(), 0.0);
    }

    #[test]
    fn test_disk_cache_create_directory() {
        // Test creating a new cache directory
        let cache_path = std::env::temp_dir().join("test_eop_new_dir/subdir");

        // Remove if exists
        let _ = std::fs::remove_dir_all(&cache_path);

        let cache = EopDiskCache::new(&cache_path);
        assert!(cache.is_ok());
        assert!(cache_path.exists());
    }

    #[test]
    fn test_encode_decode_with_dx_dy() {
        let cache_path = std::env::temp_dir().join("test_eop_dx_dy");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        // Create records with dx/dy fields
        let mut record1 = create_test_record(59945.0);
        record1.dx_encoded = Some(100);
        record1.dy_encoded = Some(200);

        let mut record2 = create_test_record(59946.0);
        record2.dx_encoded = Some(150);
        record2.dy_encoded = None; // Only dx

        let mut record3 = create_test_record(59947.0);
        record3.dx_encoded = None;
        record3.dy_encoded = Some(250); // Only dy

        let records = vec![record1, record2, record3];

        let encoded = cache.encode_eop_records(&records).unwrap();
        let decoded = cache.decode_eop_records(&encoded).unwrap();

        assert_eq!(decoded.len(), 3);

        // Check first record has both dx and dy
        assert_eq!(decoded[0].dx_encoded, Some(100));
        assert_eq!(decoded[0].dy_encoded, Some(200));

        // Check second has only dx
        assert_eq!(decoded[1].dx_encoded, Some(150));
        assert_eq!(decoded[1].dy_encoded, None);

        // Check third has only dy
        assert_eq!(decoded[2].dx_encoded, None);
        assert_eq!(decoded[2].dy_encoded, Some(250));
    }

    #[test]
    fn test_decode_empty_record_count() {
        let cache_path = std::env::temp_dir().join("test_eop_empty_count");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        // Create valid header with 0 records
        let mut data = Vec::from(b"EOP2" as &[u8]);
        data.extend_from_slice(&0u32.to_le_bytes()); // record_count = 0
        data.extend_from_slice(&0u32.to_le_bytes()); // crc32 = 0

        let result = cache.decode_eop_records(&data);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 0);
    }

    #[test]
    fn test_decode_insufficient_data_for_first_record() {
        let cache_path = std::env::temp_dir().join("test_eop_insufficient");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        // Valid header but not enough data for first record
        let mut data = Vec::from(b"EOP2" as &[u8]);
        data.extend_from_slice(&1u32.to_le_bytes()); // record_count = 1
        data.extend_from_slice(&0u32.to_le_bytes()); // crc32
                                                     // Only add 10 bytes, need 24+
        data.extend_from_slice(&[0u8; 10]);

        let result = cache.decode_eop_records(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_decode_missing_dx_data() {
        let cache_path = std::env::temp_dir().join("test_eop_missing_dx");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let mut data = Vec::from(b"EOP2" as &[u8]);
        data.extend_from_slice(&1u32.to_le_bytes()); // record_count = 1
        data.extend_from_slice(&0u32.to_le_bytes()); // crc32

        // First record
        data.extend_from_slice(&59945.0_f64.to_le_bytes()); // mjd
        data.extend_from_slice(&100i32.to_le_bytes()); // x_p
        data.extend_from_slice(&200i32.to_le_bytes()); // y_p
        data.extend_from_slice(&10i32.to_le_bytes()); // ut1_utc
        data.extend_from_slice(&1i32.to_le_bytes()); // lod
        data.push(0x01); // flags: dx present
                         // Don't add dx data - should error

        let result = cache.decode_eop_records(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_decode_missing_dy_data() {
        let cache_path = std::env::temp_dir().join("test_eop_missing_dy");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let mut data = Vec::from(b"EOP2" as &[u8]);
        data.extend_from_slice(&1u32.to_le_bytes());
        data.extend_from_slice(&0u32.to_le_bytes());

        data.extend_from_slice(&59945.0_f64.to_le_bytes());
        data.extend_from_slice(&100i32.to_le_bytes());
        data.extend_from_slice(&200i32.to_le_bytes());
        data.extend_from_slice(&10i32.to_le_bytes());
        data.extend_from_slice(&1i32.to_le_bytes());
        data.push(0x02); // flags: dy present
                         // Don't add dy data - should error

        let result = cache.decode_eop_records(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_decode_truncated_subsequent_record() {
        let cache_path = std::env::temp_dir().join("test_eop_truncated");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let mut data = Vec::from(b"EOP2" as &[u8]);
        data.extend_from_slice(&2u32.to_le_bytes()); // 2 records
        data.extend_from_slice(&0u32.to_le_bytes());

        // First record complete
        data.extend_from_slice(&59945.0_f64.to_le_bytes());
        data.extend_from_slice(&100i32.to_le_bytes());
        data.extend_from_slice(&200i32.to_le_bytes());
        data.extend_from_slice(&10i32.to_le_bytes());
        data.extend_from_slice(&1i32.to_le_bytes());
        data.push(0x00); // no dx/dy

        // Second record truncated - only add 3 bytes, need 5
        data.extend_from_slice(&[0u8; 3]);

        let result = cache.decode_eop_records(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_decode_subsequent_record_missing_dx() {
        let cache_path = std::env::temp_dir().join("test_eop_subsequent_dx");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let mut data = Vec::from(b"EOP2" as &[u8]);
        data.extend_from_slice(&2u32.to_le_bytes());
        data.extend_from_slice(&0u32.to_le_bytes());

        // First record
        data.extend_from_slice(&59945.0_f64.to_le_bytes());
        data.extend_from_slice(&100i32.to_le_bytes());
        data.extend_from_slice(&200i32.to_le_bytes());
        data.extend_from_slice(&10i32.to_le_bytes());
        data.extend_from_slice(&1i32.to_le_bytes());
        data.push(0x00);

        // Second record
        data.extend_from_slice(&1i16.to_le_bytes()); // mjd_delta
        data.extend_from_slice(&10i16.to_le_bytes()); // x_p_delta
        data.extend_from_slice(&20i16.to_le_bytes()); // y_p_delta
        data.push(0x01); // dx present
                         // Missing dx data

        let result = cache.decode_eop_records(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_decode_subsequent_record_missing_dy() {
        let cache_path = std::env::temp_dir().join("test_eop_subsequent_dy");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let mut data = Vec::from(b"EOP2" as &[u8]);
        data.extend_from_slice(&2u32.to_le_bytes());
        data.extend_from_slice(&0u32.to_le_bytes());

        // First record
        data.extend_from_slice(&59945.0_f64.to_le_bytes());
        data.extend_from_slice(&100i32.to_le_bytes());
        data.extend_from_slice(&200i32.to_le_bytes());
        data.extend_from_slice(&10i32.to_le_bytes());
        data.extend_from_slice(&1i32.to_le_bytes());
        data.push(0x00);

        // Second record
        data.extend_from_slice(&1i16.to_le_bytes());
        data.extend_from_slice(&10i16.to_le_bytes());
        data.extend_from_slice(&20i16.to_le_bytes());
        data.push(0x02); // dy present
                         // Missing dy data

        let result = cache.decode_eop_records(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_load_nonexistent_file() {
        let cache_path = std::env::temp_dir().join("test_eop_nonexist");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let result = cache.load("nonexistent", 59945.0, 59946.0);
        assert!(result.is_ok());
        assert!(result.unwrap().is_none());
    }

    #[test]
    fn test_decode_file_too_small() {
        let cache_path = std::env::temp_dir().join("test_eop_too_small");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let data = vec![1, 2, 3]; // Only 3 bytes
        let result = cache.decode_eop_records(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_decode_wrong_magic_bytes() {
        let cache_path = std::env::temp_dir().join("test_eop_wrong_magic");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let data = vec![b'X', b'X', b'X', b'X', 0, 0, 0, 0];
        let result = cache.decode_eop_records(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_decode_eop2_file_too_small_for_crc() {
        let cache_path = std::env::temp_dir().join("test_eop_no_crc");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        // EOP2 header but no CRC field
        let data = Vec::from(b"EOP2\x01\x00\x00\x00" as &[u8]);
        let result = cache.decode_eop_records(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_disk_cache_has_data() {
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let cache_path = std::env::temp_dir().join(format!("test_eop_has_data_{}", timestamp));

        let cache = EopDiskCache::new(&cache_path).unwrap();

        let records = vec![create_test_record(59945.0), create_test_record(59946.0)];

        // No data initially
        assert!(!cache.has_data("test", 59945.0, 59946.0));

        // Store data
        cache.store("test", &records).unwrap();

        // Now has data
        assert!(cache.has_data("test", 59945.0, 59946.0));

        // Clean up
        let _ = std::fs::remove_dir_all(&cache_path);
    }

    #[test]
    fn test_disk_cache_clear_all() {
        let cache_path = std::env::temp_dir().join("test_eop_clear_all");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        let records = vec![create_test_record(59945.0), create_test_record(59946.0)];

        cache.store("test1", &records).unwrap();
        cache.store("test2", &records).unwrap();

        assert!(cache.has_data("test1", 59945.0, 59946.0));
        assert!(cache.has_data("test2", 59945.0, 59946.0));

        // Clear all cache files
        cache.clear_all().unwrap();

        assert!(!cache.has_data("test1", 59945.0, 59946.0));
        assert!(!cache.has_data("test2", 59945.0, 59946.0));
    }

    #[test]
    fn test_disk_cache_size() {
        let cache_path = std::env::temp_dir().join("test_eop_cache_size");
        let _ = std::fs::remove_dir_all(&cache_path);
        let cache = EopDiskCache::new(&cache_path).unwrap();

        // Empty initially
        let size = cache.cache_size().unwrap();
        assert_eq!(size, 0);

        // Store some data
        let records = vec![create_test_record(59945.0), create_test_record(59946.0)];
        cache.store("test", &records).unwrap();

        // Should have size > 0
        let size = cache.cache_size().unwrap();
        assert!(size > 0);
    }

    #[test]
    fn test_cache_system_with_disk_cache() {
        let cache_path = std::env::temp_dir().join("test_eop_system_disk");
        let cache = EopCacheSystem::with_disk_cache(10, &cache_path).unwrap();

        assert!(cache.disk_cache.is_some());
    }

    #[test]
    fn test_cache_system_insert_batch() {
        let mut cache = EopCacheSystem::memory_only(10);

        let records = vec![
            create_test_record(59945.0),
            create_test_record(59946.0),
            create_test_record(59947.0),
        ];

        cache.insert_batch(records);

        assert!(cache.get(59945.0).is_some());
        assert!(cache.get(59946.0).is_some());
        assert!(cache.get(59947.0).is_some());
    }

    #[test]
    fn test_cache_system_clear() {
        let mut cache = EopCacheSystem::memory_only(10);

        cache.insert(create_test_record(59945.0));
        cache.insert(create_test_record(59946.0));

        assert!(cache.get(59945.0).is_some());

        cache.clear().unwrap();

        assert!(cache.get(59945.0).is_none());
    }

    #[test]
    fn test_disk_cache_get_for_mjd() {
        let cache_path = std::env::temp_dir().join("test_eop_get_mjd");
        let cache = EopDiskCache::new(&cache_path).unwrap();

        // get_for_mjd returns None when no cached data exists
        assert!(cache.get_for_mjd(59945.0).is_none());
    }
}
