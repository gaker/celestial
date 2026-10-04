// Every scale is a split Julian Date read on that scale. `time_scale!` writes
// what all eight share. `uniform_day_scale!` adds the constructors and parsing
// for scales whose days are all 86,400 s; UTC writes its own, because a UTC day
// can end in a leap second.

macro_rules! time_scale {
    ($(#[$meta:meta])* $name:ident $(, j2000_note = $note:literal)?) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        pub struct $name($crate::julian::JulianDate);

        impl $name {
            #[doc = concat!("Creates ", stringify!($name), " from a Julian Date.")]
            pub fn from_julian_date(jd: $crate::julian::JulianDate) -> Self {
                Self(jd)
            }

            #[doc = concat!(
                "Returns ", stringify!($name), " at JD 2451545.0 on the ", stringify!($name),
                " scale (2000-01-01 12:00:00 ", stringify!($name), ")."
            )]
            $(#[doc = ""] #[doc = $note])?
            pub fn j2000() -> Self {
                Self($crate::julian::JulianDate::j2000())
            }

            /// Returns the underlying Julian Date.
            pub fn to_julian_date(&self) -> $crate::julian::JulianDate {
                self.0
            }

            #[doc = concat!("Returns a new ", stringify!($name), " offset by the given seconds.")]
            ///
            /// The fraction of a day is added to jd2 and whole days are carried into jd1.
            pub fn add_seconds(&self, seconds: f64) -> Self {
                Self(self.0.add_seconds(seconds))
            }

            #[doc = concat!("Returns a new ", stringify!($name), " offset by the given days.")]
            ///
            /// The fraction of a day is added to jd2 and whole days are carried into jd1.
            pub fn add_days(&self, days: f64) -> Self {
                Self(self.0.add_days(days))
            }
        }

        #[doc = concat!("Formats as `", stringify!($name), " JD {jd1} + {jd2}`.")]
        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, concat!(stringify!($name), " {}"), self.0)
            }
        }

        #[doc = concat!("Converts a Julian Date to ", stringify!($name), ".")]
        impl From<$crate::julian::JulianDate> for $name {
            fn from(jd: $crate::julian::JulianDate) -> Self {
                Self::from_julian_date(jd)
            }
        }
    };
}

macro_rules! uniform_day_scale {
    ($name:ident, $from_calendar:ident) => {
        impl $name {
            #[doc = concat!(
                "Creates ", stringify!($name), " from seconds and nanoseconds since ",
                "1970-01-01 00:00:00 ", stringify!($name), "."
            )]
            ///
            /// The count is read on this scale, not as UTC, and every day is 86,400 s.
            pub fn new(seconds: i64, nanos: u32) -> Self {
                Self($crate::julian::JulianDate::from_unix_time(seconds, nanos))
            }
        }

        #[doc = concat!(
            "Parses an ISO 8601 string such as \"2000-01-01T12:00:00.123\" as ",
            stringify!($name), ", with no scale conversion."
        )]
        impl std::str::FromStr for $name {
            type Err = $crate::TimeError;

            fn from_str(s: &str) -> $crate::TimeResult<Self> {
                let parsed = $crate::parsing::parse_iso8601(s)?;
                parsed.to_julian_date().map(Self::from_julian_date)
            }
        }

        #[doc = concat!(
            "Creates ", stringify!($name), " from Gregorian calendar components, read directly as ",
            stringify!($name), " with no leap-second or scale corrections. For a UTC date and ",
            "time, use `utc_from_calendar` and convert."
        )]
        ///
        /// # Arguments
        ///
        /// * `year` - Astronomical year (0 is 1 BCE, -1 is 2 BCE)
        /// * `month` - Month (1-12)
        /// * `day` - Day of month (1-31)
        /// * `hour` - Hour (0-23)
        /// * `minute` - Minute (0-59)
        /// * `second` - Second with fractional part (0.0 to <60.0)
        pub fn $from_calendar(
            year: i32,
            month: u8,
            day: u8,
            hour: u8,
            minute: u8,
            second: f64,
        ) -> $crate::TimeResult<$name> {
            $crate::julian::JulianDate::from_calendar(year, month, day, hour, minute, second)
                .map($name::from_julian_date)
        }
    };
}

#[cfg(test)]
mod tests {
    use crate::constants::UNIX_EPOCH_JD;
    use crate::julian::JulianDate;
    use crate::scales::gps::{gps_from_calendar, GPS};
    use crate::scales::tai::{tai_from_calendar, TAI};
    use crate::scales::tcb::{tcb_from_calendar, TCB};
    use crate::scales::tcg::{tcg_from_calendar, TCG};
    use crate::scales::tdb::{tdb_from_calendar, TDB};
    use crate::scales::tt::{tt_from_calendar, TT};
    use crate::scales::ut1::{ut1_from_calendar, UT1};
    use crate::scales::utc::UTC;
    use crate::TimeError;
    use celestial_core::constants::{J2000_JD, SECONDS_PER_DAY_F64};
    use std::str::FromStr;

    macro_rules! check_time_scale {
        ($name:ident) => {{
            let label = stringify!($name);
            let jd = JulianDate::new(J2000_JD, 0.123456789);
            let stored = $name::from_julian_date(jd).to_julian_date();
            assert_eq!(stored.parts(), jd.parts(), "{label}");
            assert_eq!($name::from(jd), $name::from_julian_date(jd), "{label}");
            let j2000 = $name::j2000();
            assert_eq!(j2000.to_julian_date().parts(), (J2000_JD, 0.0), "{label}");
            assert_eq!(
                j2000.add_days(1.0).to_julian_date().parts(),
                (J2000_JD + 1.0, 0.0),
                "{label}"
            );
            assert_eq!(
                j2000.add_seconds(3600.0).to_julian_date().parts(),
                (J2000_JD, 3600.0 / SECONDS_PER_DAY_F64),
                "{label}"
            );
            let half_day = $name::from_julian_date(JulianDate::new(J2000_JD, 0.5));
            assert_eq!(half_day.to_string(), format!("{label} JD 2451545 + 0.5"));
        }};
    }

    macro_rules! check_uniform_day_scale {
        ($name:ident, $from_calendar:ident) => {{
            let label = stringify!($name);
            let unix_epoch = $name::new(0, 0).to_julian_date().parts();
            assert_eq!(unix_epoch, (UNIX_EPOCH_JD, 0.0), "{label}");
            assert_eq!(
                $name::new(1_700_000_000, 1).to_julian_date().parts(),
                (2460262.5, 0.9259259259259376),
                "{label}"
            );
            assert_eq!(
                $from_calendar(2000, 1, 1, 12, 0, 0.0)
                    .unwrap()
                    .to_julian_date()
                    .parts(),
                (2451544.5, 0.5),
                "{label}"
            );
            assert_eq!(
                $name::from_str("2000-01-01T12:00:00.123"),
                $from_calendar(2000, 1, 1, 12, 0, 0.123),
                "{label}"
            );
            assert_eq!(
                $name::from_str("2000-01-01T12:00:00Z"),
                Err(TimeError::ParseError(
                    "Only UTC takes a Z suffix: '2000-01-01T12:00:00Z'".into()
                )),
                "{label}"
            );
            assert_eq!(
                $name::from_str("invalid-date"),
                Err(TimeError::ParseError(
                    "Invalid datetime format: 'invalid-date'. Expected YYYY-MM-DDTHH:MM:SS".into()
                )),
                "{label}"
            );
        }};
    }

    #[test]
    fn test_every_scale_has_the_shared_api() {
        check_time_scale!(GPS);
        check_time_scale!(TAI);
        check_time_scale!(TCB);
        check_time_scale!(TCG);
        check_time_scale!(TDB);
        check_time_scale!(TT);
        check_time_scale!(UT1);
        check_time_scale!(UTC);
    }

    #[test]
    fn test_scales_with_86400_s_days() {
        check_uniform_day_scale!(GPS, gps_from_calendar);
        check_uniform_day_scale!(TAI, tai_from_calendar);
        check_uniform_day_scale!(TCB, tcb_from_calendar);
        check_uniform_day_scale!(TCG, tcg_from_calendar);
        check_uniform_day_scale!(TDB, tdb_from_calendar);
        check_uniform_day_scale!(TT, tt_from_calendar);
        check_uniform_day_scale!(UT1, ut1_from_calendar);
    }

    #[cfg(feature = "serde")]
    macro_rules! check_serde_round_trip {
        ($name:ident, $from_calendar:ident) => {{
            let cases = [
                $name::j2000(),
                $name::from_julian_date(JulianDate::new(J2000_JD, 0.123456789)),
                $from_calendar(2024, 6, 15, 14, 30, 45.123).unwrap(),
                $from_calendar(1990, 12, 31, 23, 59, 59.999999999).unwrap(),
            ];
            for original in cases {
                let json = serde_json::to_string(&original).unwrap();
                let back: $name = serde_json::from_str(&json).unwrap();
                let back = back.to_julian_date().parts();
                assert_eq!(back, original.to_julian_date().parts(), "{json}");
            }
        }};
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_serde_round_trip_is_exact() {
        use crate::scales::utc::utc_from_calendar;

        check_serde_round_trip!(GPS, gps_from_calendar);
        check_serde_round_trip!(TAI, tai_from_calendar);
        check_serde_round_trip!(TCB, tcb_from_calendar);
        check_serde_round_trip!(TCG, tcg_from_calendar);
        check_serde_round_trip!(TDB, tdb_from_calendar);
        check_serde_round_trip!(TT, tt_from_calendar);
        check_serde_round_trip!(UT1, ut1_from_calendar);
        check_serde_round_trip!(UTC, utc_from_calendar);
    }
}
