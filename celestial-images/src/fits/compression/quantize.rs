use crate::fits::errors::{FitsError, Result};
use std::sync::OnceLock;

const RANDOM_COUNT: usize = 10000;

// fpack writes zero pixels of SUBTRACTIVE_DITHER_2 images as -2147483646.
// FITS 4.0 §10.2.1 names -2147483647, the value fpack uses for ZBLANK.
const DITHER_2_ZERO: i64 = -2147483646;

// Dithered files written before ZDITHER0 was defined omit it and used seed 1.
const DEFAULT_DITHER_SEED: i64 = 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Method {
    NoDither,
    SubtractiveDither1,
    SubtractiveDither2,
}

// The ZSCALE, ZZERO and ZBLANK values of one tile, and its seed for the
// dither sequence.
pub(super) struct TileScaling {
    pub(super) scale: f64,
    pub(super) zero: f64,
    pub(super) blank: Option<i64>,
    pub(super) method: Method,
    pub(super) dither_seed: Option<i64>,
}

// FITS 4.0 Appendix I: the Park and Miller (1988) generator, run in double
// precision and stored in single precision.
fn random_values() -> &'static [f32; RANDOM_COUNT] {
    static VALUES: OnceLock<[f32; RANDOM_COUNT]> = OnceLock::new();
    VALUES.get_or_init(|| {
        let a = 16807.0;
        let m = 2147483647.0;
        let mut seed: f64 = 1.0;
        let mut values = [0.0f32; RANDOM_COUNT];
        for value in values.iter_mut() {
            let temp = a * seed;
            seed = temp - m * libm::trunc(temp / m);
            *value = (seed / m) as f32;
        }
        values
    })
}

// The dither sequence of the tiled-image convention (version 2.3, §4.1.2).
// `first` counts from 0 for ZDITHER0 = 1 on table row 1, and `next` restarts
// at a new offset after random value 9999, as fpack does. FITS 4.0 §10.2.1
// says the restart comes after value 499, which neither fpack nor astropy do.
struct DitherSequence {
    first: usize,
    next: usize,
}

impl DitherSequence {
    fn new(row: usize, zdither0: i64) -> Result<Self> {
        let seed = usize::try_from(zdither0)
            .ok()
            .filter(|s| (1..=RANDOM_COUNT).contains(s))
            .ok_or_else(|| FitsError::InvalidKeywordValue {
                keyword: "ZDITHER0".to_string(),
                value: zdither0.to_string(),
            })?;
        let first = (row % RANDOM_COUNT + seed - 1) % RANDOM_COUNT;
        Ok(Self {
            first,
            next: offset(first),
        })
    }

    fn advance(&mut self) -> f64 {
        let value = f64::from(random_values()[self.next]);
        self.next += 1;
        if self.next == RANDOM_COUNT {
            self.first = (self.first + 1) % RANDOM_COUNT;
            self.next = offset(self.first);
        }
        value
    }
}

fn offset(first: usize) -> usize {
    (random_values()[first] * 500.0) as usize
}

// FITS 4.0 eq. 14. fpack evaluates (I - R + 0.5) * ZSCALE + ZZERO, and
// I * ZSCALE + ZZERO without dither, with one rounding. `fma` matches it bit
// for bit where the two terms nearly cancel. astropy rounds twice without
// dither.
pub(super) fn unquantize(values: &[i64], tile: &TileScaling, row: usize) -> Result<Vec<f64>> {
    if !tile.scale.is_finite() || !tile.zero.is_finite() {
        return Err(FitsError::InvalidFormat(format!(
            "tile {} has ZSCALE {} and ZZERO {}",
            row + 1,
            tile.scale,
            tile.zero
        )));
    }
    let mut dither = dither_sequence(tile, row)?;
    let zero_marker = (tile.method == Method::SubtractiveDither2).then_some(DITHER_2_ZERO);
    let out = values
        .iter()
        .map(|&v| {
            // Undefined and zero pixels still use up one random value.
            let r = dither.as_mut().map(DitherSequence::advance);
            if Some(v) == tile.blank {
                f64::NAN
            } else if Some(v) == zero_marker {
                0.0
            } else {
                let level = r.map_or(v as f64, |r| v as f64 - r + 0.5);
                libm::fma(level, tile.scale, tile.zero)
            }
        })
        .collect();
    Ok(out)
}

fn dither_sequence(tile: &TileScaling, row: usize) -> Result<Option<DitherSequence>> {
    match tile.method {
        Method::NoDither => Ok(None),
        Method::SubtractiveDither1 | Method::SubtractiveDither2 => {
            let seed = tile.dither_seed.unwrap_or(DEFAULT_DITHER_SEED);
            DitherSequence::new(row, seed).map(Some)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tile(method: Method, blank: Option<i64>, seed: Option<i64>) -> TileScaling {
        TileScaling {
            scale: 2.0,
            zero: 10.0,
            blank,
            method,
            dither_seed: seed,
        }
    }

    // FITS 4.0 Appendix I: the 10000th seed is 1043618065.
    #[test]
    fn random_values_end_at_the_published_seed() {
        let last = random_values()[RANDOM_COUNT - 1];
        assert_eq!(last, (1043618065.0f64 / 2147483647.0) as f32);
        assert_eq!(random_values()[0], (16807.0f64 / 2147483647.0) as f32);
    }

    #[test]
    fn no_dither_is_linear() {
        let out = unquantize(&[0, 1, -3], &tile(Method::NoDither, None, None), 0).unwrap();
        assert_eq!(out, vec![10.0, 12.0, 4.0]);
    }

    #[test]
    fn blank_and_dither_2_zero_map_to_nan_and_zero() {
        let t = tile(Method::SubtractiveDither2, Some(-2147483647), Some(1));
        let out = unquantize(&[-2147483647, DITHER_2_ZERO, 0], &t, 0).unwrap();
        assert!(out[0].is_nan());
        assert_eq!(out[1], 0.0);
        let r = f64::from(random_values()[offset(0) + 2]);
        assert_eq!(out[2], libm::fma(0.5 - r, 2.0, 10.0));
    }

    #[test]
    fn dither_1_does_not_treat_the_dither_2_marker_as_zero() {
        let t = tile(Method::SubtractiveDither1, None, Some(1));
        let out = unquantize(&[DITHER_2_ZERO], &t, 0).unwrap();
        assert_ne!(out[0], 0.0);
    }

    #[test]
    fn dither_sequence_starts_from_row_and_seed() {
        let mut seq = DitherSequence::new(0, 1).unwrap();
        assert_eq!(seq.first, 0);
        assert_eq!(seq.next, offset(0));
        let first = seq.advance();
        assert_eq!(first, f64::from(random_values()[offset(0)]));

        let seq = DitherSequence::new(3, 10000).unwrap();
        assert_eq!(seq.first, 2);
        let seq = DitherSequence::new(20000, 1).unwrap();
        assert_eq!(seq.first, 0);
    }

    #[test]
    fn dither_sequence_restarts_after_value_9999() {
        let mut seq = DitherSequence::new(0, 10000).unwrap();
        assert_eq!(seq.first, 9999);
        while seq.next != RANDOM_COUNT - 1 {
            seq.advance();
        }
        seq.advance();
        assert_eq!(seq.first, 0);
        assert_eq!(seq.next, offset(0));
    }

    // The convention writes INT(RN(I0) * 500.) without a precision. Single
    // and double precision give the same offsets for all 10000 values.
    #[test]
    fn offset_does_not_depend_on_the_multiply_precision() {
        for (i, &r) in random_values().iter().enumerate() {
            assert_eq!(offset(i), (f64::from(r) * 500.0) as usize, "{i}");
        }
    }

    #[test]
    fn dither_without_a_seed_uses_seed_1() {
        let values = [3, -7, 12];
        let unseeded = tile(Method::SubtractiveDither1, None, None);
        let seeded = tile(Method::SubtractiveDither1, None, Some(1));
        assert_eq!(
            unquantize(&values, &unseeded, 5).unwrap(),
            unquantize(&values, &seeded, 5).unwrap()
        );
    }

    #[test]
    fn dither_requires_a_seed_from_1_to_10000() {
        for seed in [Some(0), Some(10001), Some(-5)] {
            let t = tile(Method::SubtractiveDither1, None, seed);
            assert!(unquantize(&[1], &t, 0).is_err(), "{seed:?}");
        }
    }

    #[test]
    fn non_finite_scaling_is_an_error() {
        let mut t = tile(Method::NoDither, None, None);
        t.scale = f64::INFINITY;
        assert!(unquantize(&[1], &t, 0).is_err());
        t.scale = 1.0;
        t.zero = f64::NAN;
        assert!(unquantize(&[1], &t, 0).is_err());
    }
}
