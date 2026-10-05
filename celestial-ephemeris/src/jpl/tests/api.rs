use super::fixture::*;
use crate::jpl::spk::SpkFile;
use crate::jpl::{bodies, SpkError};
use celestial_core::constants::J2000_JD;
use celestial_core::errors::{AstroError, MathErrorKind};
use celestial_core::matrix::Vector3;
use std::io::{self, ErrorKind, Read};

struct Trickle<'a> {
    bytes: &'a [u8],
    interrupt: bool,
}

impl Read for Trickle<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.interrupt = !self.interrupt;
        if self.interrupt {
            return Err(ErrorKind::Interrupted.into());
        }
        let n = buf.len().min(1001).min(self.bytes.len());
        buf[..n].copy_from_slice(&self.bytes[..n]);
        self.bytes = &self.bytes[n..];
        Ok(n)
    }
}

struct Failing<'a> {
    bytes: &'a [u8],
}

impl Read for Failing<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.bytes.is_empty() {
            return Err(io::Error::other("disk on fire"));
        }
        let n = buf.len().min(self.bytes.len());
        buf[..n].copy_from_slice(&self.bytes[..n]);
        self.bytes = &self.bytes[n..];
        Ok(n)
    }
}

#[test]
fn de432s_segments_are_listed_in_file_order() {
    let spk = load(&de432s()).unwrap();
    let pairs: Vec<_> = spk
        .segments()
        .iter()
        .map(|s| (s.body(), s.center()))
        .collect();
    assert_eq!(pairs, PAIRS);
    for segment in spk.segments() {
        assert_eq!((segment.frame(), segment.data_type()), (1, 2));
        let (start, end) = (
            segment.start().to_julian_date(),
            segment.end().to_julian_date(),
        );
        assert_eq!((start.jd1(), start.jd2()), (J2000_JD, -18280.5));
        assert_eq!((end.jd1(), end.jd2()), (J2000_JD, 18263.5));
    }
}

#[test]
fn states_match_the_kernel() {
    let spk = load(&de432s()).unwrap();
    let emb = state(
        &spk,
        bodies::EARTH_MOON_BARYCENTER,
        bodies::SOLAR_SYSTEM_BARYCENTER,
        J2000_JD,
        0.0,
    );
    assert_eq!(
        emb.unwrap(),
        (
            Vector3::new(-27570175.588424373, 132358187.63084124, 57417722.61070698),
            Vector3::new(-29.77712821627309, -5.037847183456792, -2.184306335195768)
        )
    );
    let moon = state(&spk, bodies::MOON, bodies::EARTH, 2460000.5, 0.123);
    assert_eq!(
        moon.unwrap(),
        (
            Vector3::new(293419.33570357383, 224376.11219133495, 100117.60279793796),
            Vector3::new(-0.5985623931708289, 0.722441550541203, 0.4123227368870718)
        )
    );
}

#[test]
fn every_body_constant_reaches_the_barycenter() {
    let spk = load(&de432s()).unwrap();
    for body in [
        bodies::MERCURY_BARYCENTER,
        bodies::VENUS_BARYCENTER,
        bodies::EARTH_MOON_BARYCENTER,
        bodies::MARS_BARYCENTER,
        bodies::JUPITER_BARYCENTER,
        bodies::SATURN_BARYCENTER,
        bodies::URANUS_BARYCENTER,
        bodies::NEPTUNE_BARYCENTER,
        bodies::PLUTO_BARYCENTER,
        bodies::SUN,
        bodies::MERCURY,
        bodies::VENUS,
        bodies::MOON,
        bodies::EARTH,
    ] {
        let center = bodies::SOLAR_SYSTEM_BARYCENTER;
        assert!(state(&spk, body, center, J2000_JD, 0.0).is_ok(), "{}", body);
    }
}

#[test]
fn open_and_any_reader_give_the_same_kernel() {
    let bytes = de432s();
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/de432s.bsp");
    let opened = SpkFile::open(&path).unwrap();
    let reader = Trickle {
        bytes: &bytes,
        interrupt: false,
    };
    let trickled = SpkFile::from_reader(reader, 0).unwrap();
    let read = load(&bytes).unwrap();
    for (body, center) in PAIRS {
        let expected = state(&read, body, center, J2000_JD, 0.37).unwrap();
        assert_eq!(
            state(&opened, body, center, J2000_JD, 0.37).unwrap(),
            expected
        );
        assert_eq!(
            state(&trickled, body, center, J2000_JD, 0.37).unwrap(),
            expected
        );
    }
}

#[test]
fn read_failures_are_io_errors() {
    let err = SpkFile::open("/nonexistent/path/file.bsp").unwrap_err();
    assert!(
        matches!(err, SpkError::Io(ref e) if e.kind() == ErrorKind::NotFound),
        "{:?}",
        err
    );
    let bytes = de432s();
    let err = SpkFile::from_reader(
        Failing {
            bytes: &bytes[..4096],
        },
        0,
    )
    .unwrap_err();
    assert_eq!(err.to_string(), "SPK read failed: disk on fire");
}

#[test]
fn errors_read_as_sentences() {
    let cases = [
        (
            SpkError::Io(io::Error::other("gone")),
            "SPK read failed: gone",
        ),
        (SpkError::InvalidFormat("x".into()), "invalid SPK file: x"),
        (SpkError::InvalidData("x".into()), "invalid SPK data: x"),
        (
            SpkError::InvalidEpoch { jd: f64::INFINITY },
            "TDB epoch JD inf is not a finite number of seconds from J2000",
        ),
        (
            SpkError::SegmentNotFound {
                body: 399,
                center: 0,
                jd: J2000_JD,
            },
            "no SPK segments connect body 399 to 0 at TDB JD 2451545",
        ),
        (
            SpkError::UnsupportedType(3),
            "SPK data type 3 is not supported; only type 2 is",
        ),
        (
            SpkError::UnsupportedFrame(17),
            "SPK frame 17 is not supported; only frame 1 (J2000) is",
        ),
    ];
    for (error, text) in cases {
        assert_eq!(error.to_string(), text);
    }
}

#[test]
fn spk_errors_convert_to_data_errors() {
    let operation = |error: SpkError| match AstroError::from(error) {
        AstroError::DataError { operation, .. } => operation,
        other => panic!("expected DataError, got {:?}", other),
    };
    let not_found = SpkError::SegmentNotFound {
        body: 399,
        center: 0,
        jd: J2000_JD,
    };
    let cases = [
        (SpkError::Io(io::Error::other("x")), "read"),
        (SpkError::InvalidFormat("x".into()), "open"),
        (SpkError::InvalidData("x".into()), "evaluate"),
        (SpkError::UnsupportedType(3), "evaluate"),
        (SpkError::UnsupportedFrame(17), "evaluate"),
        (not_found, "evaluate"),
    ];
    for (error, expected) in cases {
        assert_eq!(operation(error), expected);
    }
}

#[test]
fn spk_error_text_survives_conversion() {
    assert_eq!(
        AstroError::from(SpkError::InvalidFormat("x".into())).to_string(),
        "Data error (SPK - open): invalid SPK file: x"
    );
}

#[test]
fn invalid_epoch_converts_to_a_math_error() {
    let epoch = AstroError::from(SpkError::InvalidEpoch { jd: f64::NAN });
    assert!(
        matches!(
            epoch,
            AstroError::MathError {
                kind: MathErrorKind::NotFinite,
                ..
            }
        ),
        "{:?}",
        epoch
    );
}
