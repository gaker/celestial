use crate::ser::error::{Result, SerError};
use celestial_core::constants::NANOSECONDS_PER_SECOND;
use celestial_time::scales::utc::UTC;

// SER timestamps are .NET ticks: 100 ns since 0001-01-01T00:00:00, with no
// leap seconds. .NET keeps the DateTime Kind in bits 62 and 63, so the ticks
// are the low 62 bits.
pub(crate) const TICKS_MASK: u64 = 0x3FFF_FFFF_FFFF_FFFF;
const TICKS_PER_SECOND: i64 = 10_000_000;
const NANOS_PER_TICK: u32 = 100;
const UNIX_EPOCH_TICKS: i64 = 621_355_968_000_000_000;
// 9999-12-31T23:59:59.9999999, the last instant .NET DateTime holds.
const MAX_TICKS: i64 = 3_155_378_975_999_999_999;
// Unix seconds past 9999-12-31 or before 0001-01-01 are out of range anyway,
// and this bound keeps the f64 to i64 cast exact.
const MAX_UNIX_SECONDS_F64: f64 = 1.0e12;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SerTimestamp {
    ticks: u64,
}

impl SerTimestamp {
    pub fn new(ticks: u64) -> Self {
        Self {
            ticks: ticks & TICKS_MASK,
        }
    }

    pub fn ticks(&self) -> u64 {
        self.ticks
    }

    pub(crate) fn from_le_bytes(bytes: [u8; 8]) -> Self {
        Self::new(u64::from_le_bytes(bytes))
    }

    pub(crate) fn to_le_bytes(self) -> [u8; 8] {
        self.ticks.to_le_bytes()
    }

    // Nanoseconds below the 100 ns tick are dropped.
    pub fn from_unix(seconds: i64, nanos: u32) -> Result<Self> {
        if nanos >= NANOSECONDS_PER_SECOND {
            return Err(SerError::TimestampError(format!(
                "nanos must be less than {}, got {}",
                NANOSECONDS_PER_SECOND, nanos
            )));
        }
        let ticks = i128::from(seconds) * i128::from(TICKS_PER_SECOND)
            + i128::from(UNIX_EPOCH_TICKS)
            + i128::from(nanos / NANOS_PER_TICK);
        match u64::try_from(ticks) {
            Ok(ticks) if ticks <= MAX_TICKS as u64 => Ok(Self { ticks }),
            _ => Err(SerError::TimestampError(format!(
                "Unix time {}.{:09} s is outside 0001-01-01 to 9999-12-31",
                seconds, nanos
            ))),
        }
    }

    pub fn to_unix(&self) -> (i64, u32) {
        // The mask keeps ticks below 2^62, so the cast is exact.
        let since_epoch = self.ticks as i64 - UNIX_EPOCH_TICKS;
        let nanos = since_epoch.rem_euclid(TICKS_PER_SECOND) as u32 * NANOS_PER_TICK;
        (since_epoch.div_euclid(TICKS_PER_SECOND), nanos)
    }

    // Rounds to the nearest tick. An f64 near 1.7e9 s has a resolution of
    // about 0.24 µs, so a value read back from to_unix_seconds can differ from
    // the original by a few ticks.
    pub fn from_unix_seconds(seconds: f64) -> Result<Self> {
        if !(-MAX_UNIX_SECONDS_F64..=MAX_UNIX_SECONDS_F64).contains(&seconds) {
            return Err(SerError::TimestampError(format!(
                "Unix time {} s is not finite or outside 0001-01-01 to 9999-12-31",
                seconds
            )));
        }
        let whole = libm::floor(seconds);
        let ticks = libm::round((seconds - whole) * TICKS_PER_SECOND as f64) as u32;
        let carry = i64::from(ticks / TICKS_PER_SECOND as u32);
        let ticks = ticks % TICKS_PER_SECOND as u32;
        Self::from_unix(whole as i64 + carry, ticks * NANOS_PER_TICK)
    }

    pub fn to_unix_seconds(&self) -> f64 {
        let (seconds, nanos) = self.to_unix();
        seconds as f64 + f64::from(nanos) / f64::from(NANOSECONDS_PER_SECOND)
    }

    pub fn from_celestial_time(utc: UTC) -> Result<Self> {
        let (seconds, nanos) = utc.to_unix()?;
        Self::from_unix(seconds, nanos)
    }

    pub fn to_utc(&self) -> Result<UTC> {
        let (seconds, nanos) = self.to_unix();
        Ok(UTC::new(seconds, nanos)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_time::TimeError;
    use std::str::FromStr;

    // SerError holds io::Error, which has no PartialEq, so results are
    // compared with the error's message.
    fn shown<T>(result: Result<T>) -> std::result::Result<T, String> {
        result.map_err(|e| e.to_string())
    }

    fn timestamp_error(message: &str) -> std::result::Result<SerTimestamp, String> {
        Err(SerError::TimestampError(message.to_string()).to_string())
    }

    // 2024-06-01T12:34:56.1234567 UTC, from Python datetime (ts.py).
    const T1: u64 = 638_528_420_961_234_567;
    const T1_UNIX: (i64, u32) = (1_717_245_296, 123_456_700);

    #[test]
    fn new_masks_the_kind_bits() {
        assert_eq!(SerTimestamp::new(T1).ticks(), T1);
        assert_eq!(SerTimestamp::new(T1 | 1 << 62).ticks(), T1);
        assert_eq!(SerTimestamp::new(T1 | 1 << 63).ticks(), T1);
        assert_eq!(SerTimestamp::new(u64::MAX).ticks(), TICKS_MASK);
    }

    #[test]
    fn le_bytes_round_trip() {
        let bytes = [0xEF, 0xCD, 0xAB, 0x90, 0x78, 0x56, 0x34, 0x12];
        let ts = SerTimestamp::from_le_bytes(bytes);
        assert_eq!(ts.ticks(), 0x1234567890ABCDEF);
        assert_eq!(ts.to_le_bytes(), bytes);
        let utc_kind = (T1 | 1 << 62).to_le_bytes();
        assert_eq!(SerTimestamp::from_le_bytes(utc_kind).ticks(), T1);
    }

    #[test]
    fn unix_conversion_is_exact() {
        let ts = SerTimestamp::new(T1);
        assert_eq!(ts.to_unix(), T1_UNIX);
        assert_eq!(shown(SerTimestamp::from_unix(T1_UNIX.0, T1_UNIX.1)), Ok(ts));
        // Nanoseconds below a tick are dropped.
        assert_eq!(
            shown(SerTimestamp::from_unix(T1_UNIX.0, T1_UNIX.1 + 99)),
            Ok(ts)
        );
        assert_eq!(SerTimestamp::new(UNIX_EPOCH_TICKS as u64).to_unix(), (0, 0));
    }

    #[test]
    fn unix_conversion_before_1970() {
        let ts = SerTimestamp::new(UNIX_EPOCH_TICKS as u64 - 5_000_000);
        assert_eq!(ts.to_unix(), (-1, 500_000_000));
        assert_eq!(ts.to_unix_seconds(), -0.5);
        assert_eq!(shown(SerTimestamp::from_unix(-1, 500_000_000)), Ok(ts));
        assert_eq!(SerTimestamp::new(0).to_unix(), (-62_135_596_800, 0));
        assert_eq!(
            shown(SerTimestamp::from_unix(-62_135_596_800, 0)),
            Ok(SerTimestamp::new(0))
        );
    }

    #[test]
    fn from_unix_rejects_out_of_range_input() {
        let out_of_range = |text: &str| {
            timestamp_error(&format!(
                "Unix time {text} s is outside 0001-01-01 to 9999-12-31"
            ))
        };
        let cases = [
            ((-62_135_596_801, 999_999_999), "-62135596801.999999999"),
            ((253_402_300_800, 0), "253402300800.000000000"),
            ((i64::MAX, 0), "9223372036854775807.000000000"),
            ((i64::MIN, 0), "-9223372036854775808.000000000"),
        ];
        for ((seconds, nanos), text) in cases {
            assert_eq!(
                shown(SerTimestamp::from_unix(seconds, nanos)),
                out_of_range(text)
            );
        }
        assert_eq!(
            shown(SerTimestamp::from_unix(253_402_300_799, 999_999_999)),
            Ok(SerTimestamp::new(MAX_TICKS as u64))
        );
        assert_eq!(
            shown(SerTimestamp::from_unix(0, 1_000_000_000)),
            timestamp_error("nanos must be less than 1000000000, got 1000000000")
        );
    }

    #[test]
    fn from_unix_seconds_rounds_to_the_nearest_tick() {
        let cases = [
            (1_672_531_200.5, (1_672_531_200, 500_000_000)),
            // The f64 nearest 1672531200.1234567 is 1672531200.12345671654.
            (1_672_531_200.123_456_7, (1_672_531_200, 123_456_700)),
            (-0.000_000_01, (0, 0)),
            (0.999_999_99, (1, 0)),
            (-1.25, (-2, 750_000_000)),
        ];
        for (seconds, unix) in cases {
            let ts = SerTimestamp::from_unix_seconds(seconds).unwrap();
            assert_eq!(ts.to_unix(), unix, "{seconds}");
        }
        let half = SerTimestamp::from_unix_seconds(1_672_531_200.5).unwrap();
        assert_eq!(half.to_unix_seconds(), 1_672_531_200.5);
    }

    #[test]
    fn from_unix_seconds_rejects_non_finite_and_huge_input() {
        for seconds in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e13, -1e13] {
            assert_eq!(
                shown(SerTimestamp::from_unix_seconds(seconds)),
                timestamp_error(&format!(
                    "Unix time {seconds} s is not finite or outside 0001-01-01 to 9999-12-31"
                )),
                "{seconds}"
            );
        }
        assert_eq!(
            shown(SerTimestamp::from_unix_seconds(9e11)),
            timestamp_error(
                "Unix time 900000000000.000000000 s is outside 0001-01-01 to 9999-12-31"
            )
        );
    }

    #[test]
    fn celestial_time_conversion_is_exact_on_a_leap_second_day() {
        let cases = [
            ((1_483_185_600, 0), 636_187_824_000_000_000),
            ((1_483_228_799, 0), 636_188_255_990_000_000),
            ((1_483_228_799, 999_999_900), 636_188_255_999_999_999),
            ((1_483_228_800, 0), 636_188_256_000_000_000),
        ];
        for ((seconds, nanos), ticks) in cases {
            let utc = UTC::new(seconds, nanos).unwrap();
            let ts = SerTimestamp::from_celestial_time(utc).unwrap();
            assert_eq!(ts.ticks(), ticks, "{seconds} s {nanos} ns");
            assert_eq!(shown(ts.to_utc()), Ok(utc), "{seconds} s {nanos} ns");
        }
    }

    #[test]
    fn celestial_time_round_trip_is_exact() {
        let ts = SerTimestamp::new(T1);
        let utc = ts.to_utc().unwrap();
        assert_eq!(utc, UTC::from_str("2024-06-01T12:34:56.1234567").unwrap());
        assert_eq!(shown(SerTimestamp::from_celestial_time(utc)), Ok(ts));
        let pre_1970 = SerTimestamp::new(UNIX_EPOCH_TICKS as u64 - 1);
        let utc = pre_1970.to_utc().unwrap();
        assert_eq!(shown(SerTimestamp::from_celestial_time(utc)), Ok(pre_1970));
    }

    #[test]
    fn celestial_time_inside_a_leap_second_is_an_error() {
        let utc = UTC::from_str("2016-12-31T23:59:60.5").unwrap();
        let error = SerError::Time(TimeError::ConversionError(
            "2016-12-31 86400.5 s after midnight is in the day's extra time, \
             which Unix time cannot represent"
                .to_string(),
        ));
        assert_eq!(
            shown(SerTimestamp::from_celestial_time(utc)),
            Err(error.to_string())
        );
    }
}
