use super::record::EopRecord;
use crate::errors::{CoordError, CoordResult};
use celestial_core::constants::MJD_ZERO_POINT;
use celestial_time::scales::common::get_tai_utc_offset;
use celestial_time::scales::conversions::utc_tai::julian_to_calendar;

// TAI-UTC is cached per record because UT1-UTC steps by a whole second at every
// leap; interpolation runs on the continuous UT1-TAI instead.
#[derive(Clone)]
pub(super) struct EopTable {
    pub(super) records: Vec<EopRecord>,
    pub(super) tai_utc: Vec<f64>,
}

impl EopTable {
    pub(super) fn new(records: Vec<EopRecord>) -> CoordResult<Self> {
        if records.is_empty() {
            return Err(CoordError::data_unavailable("No EOP records supplied"));
        }
        let records = sorted_by_mjd(records)?;
        let tai_utc = records
            .iter()
            .map(|r| tai_minus_utc(r.mjd))
            .collect::<CoordResult<_>>()?;
        Ok(Self { records, tai_utc })
    }

    pub(super) fn contains_mjd(&self, mjd: f64) -> bool {
        self.records
            .binary_search_by(|r| r.mjd.total_cmp(&mjd))
            .is_ok()
    }
}

fn sorted_by_mjd(mut records: Vec<EopRecord>) -> CoordResult<Vec<EopRecord>> {
    if let Some(bad) = records.iter().find(|r| !r.mjd.is_finite()) {
        return Err(CoordError::invalid_coordinate(format!(
            "EOP record MJD {} is not finite",
            bad.mjd
        )));
    }
    records.sort_by(|a, b| a.mjd.total_cmp(&b.mjd));
    if let Some(pair) = records.windows(2).find(|w| w[0].mjd == w[1].mjd) {
        return Err(CoordError::invalid_coordinate(format!(
            "Duplicate EOP records at MJD {}",
            pair[0].mjd
        )));
    }
    Ok(records)
}

pub(super) fn tai_minus_utc(mjd: f64) -> CoordResult<f64> {
    let (year, month, day, fraction) = julian_to_calendar(MJD_ZERO_POINT, mjd)?;
    Ok(get_tai_utc_offset(year, month, day, fraction)?)
}
