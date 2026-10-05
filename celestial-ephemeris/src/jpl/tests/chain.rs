use super::fixture::*;
use crate::jpl::SpkError;
use celestial_core::constants::J2000_JD;
use celestial_core::matrix::Vector3;

const EPOCHS: [(f64, f64); 4] = [
    (J2000_JD, 0.37),
    (2433264.5, 0.0),
    (2469808.5, 0.0),
    (2460000.5, 0.123),
];

fn plus(a: State, b: State) -> State {
    (a.0 + b.0, a.1 + b.1)
}

fn minus(a: State, b: State) -> State {
    (a.0 - b.0, a.1 - b.1)
}

#[test]
fn bodies_without_a_direct_segment_chain_through_their_centers() {
    let spk = load(&de432s()).unwrap();
    for (jd1, jd2) in EPOCHS {
        let s = |body, center| state(&spk, body, center, jd1, jd2).unwrap();
        let earth_ssb = plus(s(399, 3), s(3, 0));
        let mercury_ssb = plus(s(199, 1), s(1, 0));
        let zero = (Vector3::zeros(), Vector3::zeros());
        assert_eq!(s(301, 399), minus(s(301, 3), s(399, 3)));
        assert_eq!(s(399, 301), minus(s(399, 3), s(301, 3)));
        assert_eq!(s(399, 0), earth_ssb);
        assert_eq!(s(399, 10), minus(earth_ssb, s(10, 0)));
        assert_eq!(s(10, 399), minus(s(10, 0), earth_ssb));
        assert_eq!(s(3, 10), minus(s(3, 0), s(10, 0)));
        assert_eq!(s(0, 3), minus(zero, s(3, 0)));
        assert_eq!(s(199, 399), minus(mercury_ssb, earth_ssb));
        assert_eq!(s(3, 3), zero);
    }
}

#[test]
fn position_is_the_position_of_the_state() {
    let spk = load(&de432s()).unwrap();
    for (jd1, jd2) in EPOCHS {
        for (body, center) in PAIRS.into_iter().chain([(301, 399), (0, 3), (3, 3)]) {
            assert_eq!(
                position(&spk, body, center, jd1, jd2).unwrap(),
                state(&spk, body, center, jd1, jd2).unwrap().0
            );
        }
    }
}

#[test]
fn coverage_edges_are_inclusive() {
    let spk = load(&de432s()).unwrap();
    for (jd1, jd2) in [(2433264.5, 0.0), (2469808.5, 0.0)] {
        for (body, center) in PAIRS {
            assert!(state(&spk, body, center, jd1, jd2).is_ok());
        }
    }
}

#[test]
fn epochs_outside_the_coverage_are_not_found() {
    let spk = load(&de432s()).unwrap();
    for (jd1, jd2) in [(2433264.5, -1e-6), (2469808.5, 1e-6), (J2000_JD, 1e6)] {
        let err = state(&spk, 301, 3, jd1, jd2).unwrap_err();
        assert!(
            matches!(
                err,
                SpkError::SegmentNotFound {
                    body: 301,
                    center: 3,
                    ..
                }
            ),
            "{:?}",
            err
        );
    }
}

#[test]
fn unknown_body_is_not_found() {
    let spk = load(&de432s()).unwrap();
    let err = state(&spk, 999, 0, J2000_JD, 0.25).unwrap_err();
    assert_eq!(
        err.to_string(),
        "no SPK segments connect body 999 to 0 at TDB JD 2451545.25"
    );
}

#[test]
fn non_finite_epochs_are_invalid() {
    let spk = load(&de432s()).unwrap();
    for (jd1, jd2) in [
        (f64::NAN, 0.0),
        (J2000_JD, f64::INFINITY),
        (f64::NEG_INFINITY, 0.0),
    ] {
        let err = state(&spk, 301, 3, jd1, jd2).unwrap_err();
        assert!(matches!(err, SpkError::InvalidEpoch { .. }), "{:?}", err);
        let err = position(&spk, 301, 3, jd1, jd2).unwrap_err();
        assert!(matches!(err, SpkError::InvalidEpoch { .. }), "{:?}", err);
    }
    let err = state(&spk, 301, 3, f64::NAN, 0.0).unwrap_err();
    assert_eq!(
        err.to_string(),
        "TDB epoch JD NaN is not a finite number of seconds from J2000"
    );
}

#[test]
fn epoch_keeps_both_parts() {
    let spk = load(&de432s()).unwrap();
    let later = state(&spk, 301, 3, 2460000.5, 1e-11).unwrap();
    assert_ne!(later, state(&spk, 301, 3, 2460000.5, 0.0).unwrap());
    assert_eq!(
        later,
        state(&spk, 301, 3, J2000_JD, 8455.5 + 1e-11).unwrap()
    );
}

#[test]
fn big_endian_file_gives_identical_states() {
    let little = de432s();
    let big = load(&big_endian(&little)).unwrap();
    let little = load(&little).unwrap();
    for (jd1, jd2) in EPOCHS {
        for (body, center) in PAIRS {
            assert_eq!(
                state(&big, body, center, jd1, jd2).unwrap(),
                state(&little, body, center, jd1, jd2).unwrap()
            );
        }
    }
}

fn big_endian(little: &[u8]) -> Vec<u8> {
    let mut big = little.to_vec();
    for at in [8, 12, 76, 80, 84] {
        big[at..at + 4].reverse();
    }
    big[88..96].copy_from_slice(b"BIG-IEEE");
    let record = summary_record(little);
    for k in 0..3 {
        big[record + 8 * k..record + 8 * k + 8].reverse();
    }
    for k in 0..summary_count(little) {
        swap_summary(&mut big, summary_at(little, k));
    }
    for address in 513..get_i32(little, 84) as usize {
        big[word(address)..word(address) + 8].reverse();
    }
    big
}

fn swap_summary(bytes: &mut [u8], at: usize) {
    bytes[at..at + 8].reverse();
    bytes[at + 8..at + 16].reverse();
    for i in 0..6 {
        bytes[at + 16 + 4 * i..at + 20 + 4 * i].reverse();
    }
}
