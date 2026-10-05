use crate::jpl::spk::SpkFile;
use crate::jpl::SpkError;
use celestial_core::matrix::Vector3;
use celestial_time::julian::JulianDate;
use celestial_time::scales::tdb::TDB;

pub(super) type State = (Vector3, Vector3);

const SUMMARY_WORDS: usize = 5;

pub(super) const PAIRS: [(i32, i32); 14] = [
    (1, 0),
    (2, 0),
    (3, 0),
    (4, 0),
    (5, 0),
    (6, 0),
    (7, 0),
    (8, 0),
    (9, 0),
    (10, 0),
    (301, 3),
    (399, 3),
    (199, 1),
    (299, 2),
];

pub(super) fn de432s() -> Vec<u8> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/de432s.bsp");
    std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {}", path.display(), e))
}

pub(super) fn load(bytes: &[u8]) -> Result<SpkFile, SpkError> {
    SpkFile::from_reader(bytes, bytes.len() as u64)
}

pub(super) fn rejected(bytes: &[u8]) -> String {
    match load(bytes) {
        Err(SpkError::InvalidFormat(msg)) => msg,
        Err(other) => panic!("expected InvalidFormat, got {:?}", other),
        Ok(_) => panic!("a corrupt kernel opened"),
    }
}

pub(super) fn corrupted(edit: impl FnOnce(&mut Vec<u8>)) -> Vec<u8> {
    let mut bytes = de432s();
    edit(&mut bytes);
    bytes
}

pub(super) fn tdb(jd1: f64, jd2: f64) -> TDB {
    TDB::from_julian_date(JulianDate::new(jd1, jd2))
}

pub(super) fn state(
    spk: &SpkFile,
    body: i32,
    center: i32,
    jd1: f64,
    jd2: f64,
) -> Result<State, SpkError> {
    spk.compute_state(body, center, &tdb(jd1, jd2))
}

pub(super) fn position(
    spk: &SpkFile,
    body: i32,
    center: i32,
    jd1: f64,
    jd2: f64,
) -> Result<Vector3, SpkError> {
    spk.compute_position(body, center, &tdb(jd1, jd2))
}

pub(super) fn get_f64(bytes: &[u8], at: usize) -> f64 {
    f64::from_le_bytes(bytes[at..at + 8].try_into().unwrap())
}

pub(super) fn get_i32(bytes: &[u8], at: usize) -> i32 {
    i32::from_le_bytes(bytes[at..at + 4].try_into().unwrap())
}

pub(super) fn put_f64(bytes: &mut [u8], at: usize, value: f64) {
    bytes[at..at + 8].copy_from_slice(&value.to_le_bytes());
}

pub(super) fn put_i32(bytes: &mut [u8], at: usize, value: i32) {
    bytes[at..at + 4].copy_from_slice(&value.to_le_bytes());
}

pub(super) fn summary_record(bytes: &[u8]) -> usize {
    (get_i32(bytes, 76) as usize - 1) * 1024
}

pub(super) fn summary_count(bytes: &[u8]) -> usize {
    get_f64(bytes, summary_record(bytes) + 16) as usize
}

pub(super) fn summary_at(bytes: &[u8], index: usize) -> usize {
    summary_record(bytes) + 24 + index * SUMMARY_WORDS * 8
}

pub(super) fn summary(bytes: &[u8], body: i32) -> usize {
    (0..summary_count(bytes))
        .map(|k| summary_at(bytes, k))
        .find(|&at| get_i32(bytes, at + 16) == body)
        .unwrap_or_else(|| panic!("no summary for body {}", body))
}

pub(super) fn begin_address(bytes: &[u8], body: i32) -> usize {
    get_i32(bytes, summary(bytes, body) + 32) as usize
}

pub(super) fn end_address(bytes: &[u8], body: i32) -> usize {
    get_i32(bytes, summary(bytes, body) + 36) as usize
}

pub(super) fn word(address: usize) -> usize {
    (address - 1) * 8
}

pub(super) fn directory(bytes: &[u8], body: i32) -> usize {
    word(end_address(bytes, body) - 3)
}
