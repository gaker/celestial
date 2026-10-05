use super::fixture::*;
use crate::jpl::bodies::{EARTH_MOON_BARYCENTER as EMB, MOON};
use crate::jpl::SpkError;
use celestial_core::constants::J2000_JD;

const START_JD: f64 = 2433264.5;

fn moon_error(bytes: &[u8], jd1: f64, jd2: f64) -> SpkError {
    let spk = load(bytes).unwrap();
    match state(&spk, MOON, EMB, jd1, jd2) {
        Err(e) => e,
        Ok(s) => panic!("expected an error, got {:?}", s),
    }
}

#[test]
fn segment_addresses_must_lie_inside_the_file() {
    let end = end_address(&de432s(), MOON) as i32;
    for (begin, end) in [(0, end), (end + 1, end), (612665, 1361921), (-5, end)] {
        let bytes = corrupted(|b| {
            let at = summary(b, MOON);
            put_i32(b, at + 32, begin);
            put_i32(b, at + 36, end);
        });
        assert_eq!(
            rejected(&bytes),
            format!(
                "segment 10 (body 301, center 3): data addresses {} to {} are not \
                 an ordered range inside 1 to 1361920",
                begin, end
            )
        );
    }
}

#[test]
fn type2_segment_must_hold_a_directory() {
    for length in [1, 3] {
        let bytes = corrupted(|b| {
            let at = summary(b, MOON);
            let begin = get_i32(b, at + 32);
            put_i32(b, at + 36, begin + length - 1);
        });
        assert_eq!(
            rejected(&bytes),
            format!(
                "segment 10 (body 301, center 3): {} words cannot hold a type 2 directory",
                length
            )
        );
    }
}

#[test]
fn type2_directory_must_describe_the_segment() {
    let cases: [(usize, f64, &str); 9] = [
        (16, 0.0, "RSIZE 0 is not 2 + 3n for a whole n >= 1"),
        (16, 6.0, "RSIZE 6 is not 2 + 3n for a whole n >= 1"),
        (16, 41.5, "RSIZE 41.5 is not 2 + 3n for a whole n >= 1"),
        (8, f64::NAN, "INTLEN NaN is not a positive finite number"),
        (8, 0.0, "INTLEN 0 is not a positive finite number"),
        (0, f64::INFINITY, "INIT inf is not finite"),
        (24, 0.0, "N 0 is not a whole number from 1 to 1361920"),
        (
            24,
            9135.0,
            "9135 records of 41 words plus 4 is not the segment's 374580 words",
        ),
        (
            0,
            -1579435199.0,
            "records from -1579435199 s for 9136 x 345600 s do not \
                            cover the segment's -1579435200 to 1577966400 s",
        ),
    ];
    for (offset, value, reason) in cases {
        let bytes = corrupted(|b| {
            let at = directory(b, MOON) + offset;
            put_f64(b, at, value)
        });
        assert_eq!(
            rejected(&bytes),
            format!("segment 10 (body 301, center 3): {}", reason)
        );
    }
}

#[test]
fn segment_coverage_must_be_finite_ordered_and_inside_the_records() {
    let cases: [(usize, f64, &str); 3] = [
        (
            8,
            1577966400.0 + 50.0 * 365.25 * 86400.0,
            "records from -1579435200 s for \
            9136 x 345600 s do not cover the segment's -1579435200 to 3155846400 s",
        ),
        (
            0,
            f64::NAN,
            "coverage NaN to 1577966400 s is not a finite, ordered interval",
        ),
        (
            0,
            1577966401.0,
            "coverage 1577966401 to 1577966400 s is not a finite, ordered interval",
        ),
    ];
    for (offset, value, reason) in cases {
        let bytes = corrupted(|b| {
            let at = summary(b, MOON) + offset;
            put_f64(b, at, value)
        });
        assert_eq!(
            rejected(&bytes),
            format!("segment 10 (body 301, center 3): {}", reason)
        );
    }
}

#[test]
fn non_finite_coefficient_is_an_error() {
    let bytes = corrupted(|b| {
        let at = word(begin_address(b, MOON) + 2);
        put_f64(b, at, f64::NAN)
    });
    match moon_error(&bytes, START_JD, 0.0) {
        SpkError::InvalidData(msg) => assert_eq!(
            msg,
            "record 0 of segment 10 (body 301, center 3) gives a non-finite state"
        ),
        other => panic!("expected InvalidData, got {:?}", other),
    }
}

#[test]
fn corrupt_record_midpoint_or_radius_is_an_error() {
    let radius = get_f64(&de432s(), word(begin_address(&de432s(), MOON) + 1));
    let cases = [
        (
            1,
            -radius,
            "record 0 of segment 10 (body 301, center 3) has radius -172800 s",
        ),
        (
            0,
            -1579262400.0 + 10.0 * radius,
            "TDB -1579435200 s is 11 radii from the \
            midpoint of record 0 of segment 10 (body 301, center 3)",
        ),
    ];
    for (offset, value, reason) in cases {
        let bytes = corrupted(|b| {
            let at = word(begin_address(b, MOON) + offset);
            put_f64(b, at, value)
        });
        match moon_error(&bytes, START_JD, 0.0) {
            SpkError::InvalidData(msg) => assert_eq!(msg, reason),
            other => panic!("expected InvalidData, got {:?}", other),
        }
    }
}

#[test]
fn unsupported_type_and_frame_are_reported() {
    let bytes = corrupted(|b| {
        let at = summary(b, MOON) + 28;
        put_i32(b, at, 3)
    });
    assert!(matches!(
        moon_error(&bytes, J2000_JD, 0.37),
        SpkError::UnsupportedType(3)
    ));
    let bytes = corrupted(|b| {
        let at = summary(b, MOON) + 24;
        put_i32(b, at, 17)
    });
    assert!(matches!(
        moon_error(&bytes, J2000_JD, 0.37),
        SpkError::UnsupportedFrame(17)
    ));
}

#[test]
fn later_segment_takes_priority() {
    let pristine = de432s();
    let bytes = corrupted(|b| {
        let n = summary_count(b);
        let (first, moon, slot) = (summary_at(b, 0), summary(b, MOON), summary_at(b, n));
        let count = summary_record(b) + 16;
        b.copy_within(first..first + 40, slot);
        b.copy_within(moon..moon + 40, first);
        put_i32(b, first + 16, 399);
        put_f64(b, count, (n + 1) as f64);
    });
    let real = state(&load(&pristine).unwrap(), 399, EMB, J2000_JD, 0.37).unwrap();
    assert_eq!(
        state(&load(&bytes).unwrap(), 399, EMB, J2000_JD, 0.37).unwrap(),
        real
    );
}
