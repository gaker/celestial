use super::xml::write_property;
use super::ImageData;
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::compression::XisfCompression;
use crate::xisf::header::property::XisfProperty;
use celestial_time::constants::UNIX_EPOCH_JD;
use celestial_time::scales::conversions::utc_tai::julian_to_calendar;
use quick_xml::events::{BytesEnd, BytesStart, Event};
use quick_xml::Writer;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

const SECONDS_PER_DAY: u64 = 86_400;
const CREATOR_APPLICATION: &str = concat!("celestial-images ", env!("CARGO_PKG_VERSION"));

// XISF 1.0 section 11.4: every header has one Metadata element, with the
// XISF:CreationTime and XISF:CreatorApplication properties (section
// 11.4.1). Of the optional properties of section 11.4.2,
// XISF:ChecksumAlgorithms lists the algorithm of the data block checksums
// and XISF:CompressionCodecs the codecs of the compressed image blocks.
// `created` is a TimePoint from `time_point`.
pub(super) fn write_metadata<Wr: Write>(
    writer: &mut Writer<Wr>,
    created: &str,
    sha1_checksums: bool,
    images: &[ImageData],
) -> Result<()> {
    writer.write_event(Event::Start(BytesStart::new("Metadata")))?;
    // Section 11.1.7 puts a TimePoint value in the value attribute.
    let mut creation_time = BytesStart::new("Property");
    creation_time.push_attribute(("id", "XISF:CreationTime"));
    creation_time.push_attribute(("type", "TimePoint"));
    creation_time.push_attribute(("value", created));
    writer.write_event(Event::Empty(creation_time))?;
    let creator = XisfProperty::string("XISF:CreatorApplication", CREATOR_APPLICATION);
    write_property(&creator, None, writer)?;
    if sha1_checksums {
        let algorithms = XisfProperty::string("XISF:ChecksumAlgorithms", "sha-1");
        write_property(&algorithms, None, writer)?;
    }
    if let Some(codecs) = compression_codecs(images) {
        let codecs = XisfProperty::string("XISF:CompressionCodecs", codecs);
        write_property(&codecs, None, writer)?;
    }
    writer.write_event(Event::End(BytesEnd::new("Metadata")))?;
    Ok(())
}

// The codec names of section 10.6 as a comma-separated list, in the order
// of the first image that uses each one. The writer does not shuffle bytes,
// so no name has the +sh suffix. None when no image block is compressed.
fn compression_codecs(images: &[ImageData]) -> Option<String> {
    let mut codecs: Vec<XisfCompression> = Vec::new();
    for image in images {
        let codec = image.info.compression;
        if codec.is_compressed() && !codecs.contains(&codec) {
            codecs.push(codec);
        }
    }
    let names: Vec<&str> = codecs.iter().map(XisfCompression::as_str).collect();
    (!names.is_empty()).then(|| names.join(","))
}

// XISF 1.0 section 8.4.4.4 takes an ISO 8601 extended UTC time point. It is
// written to the second with the Z zone designator, as in the examples of
// section 11.4.1. Unix time counts 86,400 s in every day, so the date and
// the time of day follow from it without the leap second table.
pub(super) fn time_point(time: SystemTime) -> Result<String> {
    let since_epoch = time.duration_since(UNIX_EPOCH).map_err(|_| {
        XisfError::CreationTime("the system clock is before 1970-01-01".to_string())
    })?;
    let seconds = since_epoch.as_secs();
    let (days, second_of_day) = (seconds / SECONDS_PER_DAY, seconds % SECONDS_PER_DAY);
    // A u64 count of days is below 2^53, so it converts to f64 exactly.
    let (year, month, day, _) = julian_to_calendar(UNIX_EPOCH_JD, days as f64)
        .map_err(|e| XisfError::CreationTime(e.to_string()))?;
    Ok(format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        second_of_day / 3600,
        second_of_day % 3600 / 60,
        second_of_day % 60
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    fn at(unix_seconds: u64) -> SystemTime {
        UNIX_EPOCH + Duration::from_secs(unix_seconds)
    }

    #[test]
    fn time_point_matches_the_example_of_section_11_4_1() {
        assert_eq!(
            time_point(at(1_418_128_695)).unwrap(),
            "2014-12-09T12:38:15Z"
        );
    }

    #[test]
    fn time_point_drops_the_fraction_of_a_second() {
        let time = at(1_418_128_695) + Duration::from_nanos(999_999_999);
        assert_eq!(time_point(time).unwrap(), "2014-12-09T12:38:15Z");
    }

    // 2016-12-31 ends in a leap second, which Unix time does not count.
    #[test]
    fn time_point_around_a_leap_second() {
        assert_eq!(
            time_point(at(1_483_228_799)).unwrap(),
            "2016-12-31T23:59:59Z"
        );
        assert_eq!(
            time_point(at(1_483_228_800)).unwrap(),
            "2017-01-01T00:00:00Z"
        );
        assert_eq!(time_point(at(0)).unwrap(), "1970-01-01T00:00:00Z");
    }

    // Before 1972 the UTC second was not the SI second, so a conversion
    // through the leap second table gives fractions of a second. Unix time
    // counts 86,400 s per day throughout.
    #[test]
    fn time_point_before_1972_counts_86400_seconds_per_day() {
        assert_eq!(time_point(at(63_071_998)).unwrap(), "1971-12-31T23:59:58Z");
        assert_eq!(time_point(at(62_986_001)).unwrap(), "1971-12-31T00:06:41Z");
    }

    #[test]
    fn time_point_before_1970_is_an_error() {
        let err = time_point(UNIX_EPOCH - Duration::from_secs(1)).unwrap_err();
        assert!(matches!(err, XisfError::CreationTime(_)), "{err}");
    }
}
