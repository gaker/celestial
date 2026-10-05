use super::fixture::*;
use crate::jpl::SpkError;
use celestial_core::constants::J2000_JD;

const RECORD_WORDS: usize = 128;
const PER_RECORD: usize = 25;
const UNSUPPORTED_TYPE: i32 = 3;

// A kernel with one summary record (plus its name record) per 25 links, and
// one data record that every segment points into. Segments use an unsupported
// data type, so only the DAF structure and the chain are exercised.
fn kernel(links: &[(i32, i32)], big: bool) -> Vec<u8> {
    let records = links.len().div_ceil(PER_RECORD).max(1);
    let data = 2 + 2 * records;
    let begin = ((data - 1) * RECORD_WORDS + 1) as i32;
    let mut words = vec![0u64; data * RECORD_WORDS];
    for r in 0..records {
        let chunk = links.chunks(PER_RECORD).nth(r).unwrap_or(&[]);
        let next = if r + 1 == records { 0 } else { 2 + 2 * (r + 1) };
        let at = (1 + 2 * r) * RECORD_WORDS;
        write_summaries(&mut words[at..at + RECORD_WORDS], next, chunk, begin, big);
    }
    let mut bytes = file_record(2 * records, data * RECORD_WORDS + 1, big);
    for word in &words[RECORD_WORDS..] {
        bytes.extend(if big {
            word.to_be_bytes()
        } else {
            word.to_le_bytes()
        });
    }
    bytes
}

fn write_summaries(record: &mut [u64], next: usize, links: &[(i32, i32)], begin: i32, big: bool) {
    record[0] = (next as f64).to_bits();
    record[2] = (links.len() as f64).to_bits();
    for (summary, &(body, center)) in record[3..].chunks_mut(5).zip(links) {
        summary[0] = (-1e9f64).to_bits();
        summary[1] = 1e9f64.to_bits();
        summary[2] = pack(body, center, big);
        summary[3] = pack(1, UNSUPPORTED_TYPE, big);
        summary[4] = pack(begin, begin + 3, big);
    }
}

fn pack(first: i32, second: i32, big: bool) -> u64 {
    let (low, high) = if big {
        (second, first)
    } else {
        (first, second)
    };
    (u64::from(high as u32) << 32) | u64::from(low as u32)
}

fn file_record(bward: usize, free: usize, big: bool) -> Vec<u8> {
    let int = |x: i32| {
        if big {
            x.to_be_bytes()
        } else {
            x.to_le_bytes()
        }
    };
    let mut record = vec![0u8; RECORD_WORDS * 8];
    record[..8].copy_from_slice(b"DAF/SPK ");
    record[8..12].copy_from_slice(&int(2));
    record[12..16].copy_from_slice(&int(6));
    record[76..80].copy_from_slice(&int(2));
    record[80..84].copy_from_slice(&int(bward as i32));
    record[84..88].copy_from_slice(&int(free as i32));
    record[88..96].copy_from_slice(if big { b"BIG-IEEE" } else { b"LTL-IEEE" });
    record[699..727].copy_from_slice(b"FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP");
    record
}

fn chain(len: i32) -> Vec<(i32, i32)> {
    (0..len).map(|k| (1000 + k, 1001 + k)).collect()
}

#[test]
fn summaries_continue_across_records_in_either_byte_order() {
    let links = chain(30);
    for big in [false, true] {
        let spk = load(&kernel(&links, big)).unwrap();
        let found: Vec<_> = spk
            .segments()
            .iter()
            .map(|s| (s.body(), s.center(), s.frame(), s.data_type()))
            .collect();
        let expected: Vec<_> = links.iter().map(|&(b, c)| (b, c, 1, 3)).collect();
        assert_eq!(found, expected);
    }
}

#[test]
fn blank_format_field_falls_back_to_the_byte_order_of_nd() {
    for fill in [b' ', 0] {
        for big in [false, true] {
            let mut bytes = kernel(&chain(2), big);
            bytes[88..96].fill(fill);
            let bodies: Vec<_> = load(&bytes)
                .unwrap()
                .segments()
                .iter()
                .map(|s| s.body())
                .collect();
            assert_eq!(bodies, [1000, 1001]);
        }
    }
}

#[test]
fn empty_kernel_opens_and_finds_nothing() {
    let spk = load(&kernel(&[], false)).unwrap();
    assert!(spk.segments().is_empty());
    let err = state(&spk, 399, 0, J2000_JD, 0.0).unwrap_err();
    assert!(matches!(err, SpkError::SegmentNotFound { .. }), "{:?}", err);
}

#[test]
fn segments_that_loop_are_invalid_data() {
    let spk = load(&kernel(&[(1000, 1001), (1001, 1000)], false)).unwrap();
    match state(&spk, 1000, 0, J2000_JD, 0.0) {
        Err(SpkError::InvalidData(msg)) => assert_eq!(
            msg,
            "segments covering TDB 0 s form a loop through body 1000"
        ),
        other => panic!("expected InvalidData, got {:?}", other),
    }
}

#[test]
fn chains_are_followed_up_to_32_links() {
    let mut links = chain(32);
    links[31].1 = 0;
    let spk = load(&kernel(&links, false)).unwrap();
    let err = state(&spk, 1000, 0, J2000_JD, 0.0).unwrap_err();
    assert!(matches!(err, SpkError::UnsupportedType(3)), "{:?}", err);
    let spk = load(&kernel(&chain(33), false)).unwrap();
    match state(&spk, 1000, 0, J2000_JD, 0.0) {
        Err(SpkError::InvalidData(msg)) => assert_eq!(
            msg,
            "segments covering TDB 0 s chain body 1000 through more than 32 centers"
        ),
        other => panic!("expected InvalidData, got {:?}", other),
    }
}
