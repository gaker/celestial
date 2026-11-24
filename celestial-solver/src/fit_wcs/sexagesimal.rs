//! Sexagesimal text for right ascension, declination and geodetic angles.
//!
//! `celestial_core::angle::{HmsFmt, DmsFmt}` truncate the units and minutes before
//! rounding the seconds, so they also print a seconds field of 60, and their text
//! (`12ʰ 30ᵐ 00ˢ`) differs from the space-separated fields these reports use.

use celestial_catalog::query::healpix::reduce_ra_deg;

// `value`, a non-negative finite angle, as whole units and "MM SS" text with `decimals`
// places of seconds. The value is rounded to the printed precision before it is split,
// so 59.9996 seconds printed to 3 places carries into the minutes instead of printing
// 60.000.
fn split_rounded(value: f64, decimals: u32) -> (u64, String) {
    let steps_per_second = 10_u64.pow(decimals);
    let steps = libm::round(value * 3600.0 * steps_per_second as f64) as u64;
    let seconds = steps / steps_per_second;
    let (minutes, second) = (seconds / 60 % 60, seconds % 60);
    let text = if decimals == 0 {
        format!("{minutes:02} {second:02}")
    } else {
        let fraction = steps % steps_per_second;
        let width = decimals as usize;
        format!("{minutes:02} {second:02}.{fraction:0width$}")
    };
    (seconds / 3600, text)
}

// Right ascension in degrees as "HH MM SS.sss", wrapped to [0h, 24h).
pub(crate) fn fmt_ra(ra_deg: f64) -> String {
    if !ra_deg.is_finite() {
        return ra_deg.to_string();
    }
    let (hours, rest) = split_rounded(reduce_ra_deg(ra_deg) / 15.0, 3);
    format!("{:02} {rest}", hours % 24)
}

// Declination in degrees as "+DD MM SS.ss".
pub(crate) fn fmt_dec(dec_deg: f64) -> String {
    if !dec_deg.is_finite() {
        return dec_deg.to_string();
    }
    let sign = if dec_deg < 0.0 { "-" } else { "+" };
    let (degrees, rest) = split_rounded(libm::fabs(dec_deg), 2);
    format!("{sign}{degrees:02} {rest}")
}

// Latitude in degrees as "DD MM SS N" or "DD MM SS S".
pub(crate) fn format_lat_deg(deg: f64) -> String {
    let hemisphere = if deg < 0.0 { "S" } else { "N" };
    format_geodetic(deg, hemisphere)
}

// Longitude in degrees as "DDD MM SS E" or "DDD MM SS W", east positive.
pub(crate) fn format_lon_deg(deg: f64) -> String {
    let hemisphere = if deg < 0.0 { "W" } else { "E" };
    format_geodetic(deg, hemisphere)
}

fn format_geodetic(deg: f64, hemisphere: &str) -> String {
    if !deg.is_finite() {
        return deg.to_string();
    }
    let (degrees, rest) = split_rounded(libm::fabs(deg), 0);
    format!("{degrees:02} {rest} {hemisphere}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fmt_ra_formats_decimal_degrees_as_hms() {
        // 180° = 12h
        assert_eq!(fmt_ra(180.0), "12 00 00.000");
        // 0° = 00h
        assert_eq!(fmt_ra(0.0), "00 00 00.000");
        // 15° = 1h
        assert_eq!(fmt_ra(15.0), "01 00 00.000");
    }

    #[test]
    fn fmt_dec_shows_sign_and_dms() {
        assert_eq!(fmt_dec(0.0), "+00 00 00.00");
        assert_eq!(fmt_dec(45.5), "+45 30 00.00");
        assert_eq!(fmt_dec(-30.25), "-30 15 00.00");
    }

    // The seconds round up to 60 and carry into the minutes and hours.
    #[test]
    fn fmt_ra_and_fmt_dec_carry_rounded_seconds() {
        assert_eq!(
            fmt_ra(15.0 * (1.0 + 59.0 / 60.0 + 59.9996 / 3600.0)),
            "02 00 00.000"
        );
        assert_eq!(fmt_ra(359.999_999_9), "00 00 00.000");
        assert_eq!(
            fmt_dec(10.0 + 59.0 / 60.0 + 59.996 / 3600.0),
            "+11 00 00.00"
        );
        assert_eq!(
            fmt_dec(-(10.0 + 59.0 / 60.0 + 59.996 / 3600.0)),
            "-11 00 00.00"
        );
        assert_eq!(fmt_dec(-0.000_000_1), "-00 00 00.00");
    }

    #[test]
    fn format_lat_deg_positive_is_north() {
        let s = format_lat_deg(33.5);
        assert!(s.ends_with(" N"), "got {s}");
        assert!(s.starts_with("33 "));
    }

    #[test]
    fn format_lat_deg_negative_is_south() {
        let s = format_lat_deg(-45.25);
        assert!(s.ends_with(" S"), "got {s}");
        assert!(s.starts_with("45 "));
    }

    #[test]
    fn format_lon_deg_positive_is_east() {
        let s = format_lon_deg(120.0);
        assert!(s.ends_with(" E"), "got {s}");
    }

    #[test]
    fn format_lon_deg_negative_is_west() {
        let s = format_lon_deg(-117.0);
        assert!(s.ends_with(" W"), "got {s}");
    }

    #[test]
    fn format_lon_and_lat_deg_carry_rounded_seconds() {
        assert_eq!(
            format_lon_deg(-(92.0 + 59.0 / 60.0 + 59.7 / 3600.0)),
            "93 00 00 W"
        );
        assert_eq!(
            format_lat_deg(38.0 + 59.0 / 60.0 + 59.5 / 3600.0),
            "39 00 00 N"
        );
    }
}
