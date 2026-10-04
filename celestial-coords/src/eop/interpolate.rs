use super::record::{EopParameters, EopRecord};
use super::table::{tai_minus_utc, EopTable};
use crate::errors::{CoordError, CoordResult};
use std::sync::Arc;

const LAGRANGE_POINTS: usize = 5;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterpolationMethod {
    Linear,

    Lagrange5,
}

pub struct EopInterpolator {
    table: Arc<EopTable>,

    method: InterpolationMethod,

    max_gap_days: f64,
}

impl EopInterpolator {
    pub fn new(records: Vec<EopRecord>) -> CoordResult<Self> {
        Ok(Self::from_table(Arc::new(EopTable::new(records)?)))
    }

    pub(super) fn from_table(table: Arc<EopTable>) -> Self {
        Self {
            table,
            method: InterpolationMethod::Linear,
            max_gap_days: 5.0,
        }
    }

    pub fn with_method(mut self, method: InterpolationMethod) -> Self {
        self.method = method;
        self
    }

    pub fn with_max_gap(mut self, max_gap_days: f64) -> CoordResult<Self> {
        if max_gap_days.is_nan() || max_gap_days <= 0.0 {
            return Err(CoordError::invalid_coordinate(format!(
                "Maximum EOP gap {} days must be positive",
                max_gap_days
            )));
        }
        self.max_gap_days = max_gap_days;
        Ok(self)
    }

    pub fn get(&self, mjd: f64) -> CoordResult<EopParameters> {
        if !mjd.is_finite() {
            return Err(CoordError::invalid_coordinate(format!(
                "EOP lookup MJD {} is not finite",
                mjd
            )));
        }
        let records = &self.table.records;
        if let Ok(idx) = records.binary_search_by(|r| r.mjd.total_cmp(&mjd)) {
            return Ok(records[idx].to_parameters());
        }
        let (before_idx, after_idx) = self.find_interpolation_interval(mjd)?;
        self.check_gap(before_idx, after_idx)?;
        let tai_utc = tai_minus_utc(mjd)?;
        match self.method {
            InterpolationMethod::Linear => {
                Ok(self.linear_interpolate(mjd, before_idx, after_idx, tai_utc))
            }
            InterpolationMethod::Lagrange5 => self.lagrange_interpolate(mjd, tai_utc),
        }
    }

    fn find_interpolation_interval(&self, mjd: f64) -> CoordResult<(usize, usize)> {
        let (first, last) = self.time_span();
        if mjd < first {
            return Err(CoordError::data_unavailable(format!(
                "MJD {:.1} is before first available record (MJD {:.1})",
                mjd, first
            )));
        }
        if mjd > last {
            return Err(CoordError::data_unavailable(format!(
                "MJD {:.1} is after last available record (MJD {:.1})",
                mjd, last
            )));
        }
        let after_idx = self.table.records.partition_point(|r| r.mjd <= mjd);
        Ok((after_idx - 1, after_idx))
    }

    fn check_gap(&self, before_idx: usize, after_idx: usize) -> CoordResult<()> {
        let records = &self.table.records;
        let gap = records[after_idx].mjd - records[before_idx].mjd;
        if gap > self.max_gap_days {
            return Err(CoordError::data_unavailable(format!(
                "Gap of {:.1} days exceeds maximum interpolation gap of {:.1} days",
                gap, self.max_gap_days
            )));
        }
        Ok(())
    }

    // UT1-UTC re-expressed against the TAI-UTC in force at the target, so the
    // interpolated quantity is UT1-TAI and a leap inside the span does not leak in.
    fn leap_adjusted(&self, idx: usize, tai_utc: f64) -> EopParameters {
        let mut params = self.table.records[idx].to_parameters();
        params.ut1_utc += tai_utc - self.table.tai_utc[idx];
        params
    }

    fn linear_interpolate(
        &self,
        mjd: f64,
        before_idx: usize,
        after_idx: usize,
        tai_utc: f64,
    ) -> EopParameters {
        let p1 = self.leap_adjusted(before_idx, tai_utc);
        let p2 = self.leap_adjusted(after_idx, tai_utc);
        let t = (mjd - p1.mjd) / (p2.mjd - p1.mjd);
        let lerp = |a: f64, b: f64| a + t * (b - a);
        let lerp_opt = |a: Option<f64>, b: Option<f64>| Some(lerp(a?, b?));
        flag_present_columns(EopParameters {
            mjd,
            x_p: lerp(p1.x_p, p2.x_p),
            y_p: lerp(p1.y_p, p2.y_p),
            ut1_utc: lerp(p1.ut1_utc, p2.ut1_utc),
            lod: lerp_opt(p1.lod, p2.lod),
            dx: lerp_opt(p1.dx, p2.dx),
            dy: lerp_opt(p1.dy, p2.dy),
            xrt: lerp_opt(p1.xrt, p2.xrt),
            yrt: lerp_opt(p1.yrt, p2.yrt),
            flags: p1.flags,
        })
    }

    fn lagrange_interpolate(&self, mjd: f64, tai_utc: f64) -> CoordResult<EopParameters> {
        let start = self.lagrange_start(mjd)?;
        let points: [EopParameters; LAGRANGE_POINTS] =
            std::array::from_fn(|k| self.leap_adjusted(start + k, tai_utc));
        let weights = lagrange_weights(mjd, &points);
        let sum = |value: fn(&EopParameters) -> f64| weighted_sum(&weights, &points, value);
        let sum_opt =
            |value: fn(&EopParameters) -> Option<f64>| weighted_sum_opt(&weights, &points, value);
        Ok(flag_present_columns(EopParameters {
            mjd,
            x_p: sum(|p| p.x_p),
            y_p: sum(|p| p.y_p),
            ut1_utc: sum(|p| p.ut1_utc),
            lod: sum_opt(|p| p.lod),
            dx: sum_opt(|p| p.dx),
            dy: sum_opt(|p| p.dy),
            xrt: sum_opt(|p| p.xrt),
            yrt: sum_opt(|p| p.yrt),
            flags: points[0].flags,
        }))
    }

    fn lagrange_start(&self, mjd: f64) -> CoordResult<usize> {
        let len = self.table.records.len();
        if len < LAGRANGE_POINTS {
            return Err(CoordError::data_unavailable(format!(
                "Not enough records for {}-point Lagrange interpolation",
                LAGRANGE_POINTS
            )));
        }
        let center = self.find_center_index(mjd);
        Ok(center
            .saturating_sub(LAGRANGE_POINTS / 2)
            .min(len - LAGRANGE_POINTS))
    }

    fn find_center_index(&self, mjd: f64) -> usize {
        let records = &self.table.records;
        let i = records.partition_point(|r| r.mjd < mjd);
        if i == 0 {
            return 0;
        }
        if i == records.len() {
            return i - 1;
        }
        let before = mjd - records[i - 1].mjd;
        let after = records[i].mjd - mjd;
        if before <= after {
            i - 1
        } else {
            i
        }
    }

    pub fn time_span(&self) -> (f64, f64) {
        let records = &self.table.records;
        (records[0].mjd, records[records.len() - 1].mjd)
    }

    pub fn record_count(&self) -> usize {
        self.table.records.len()
    }

    // A record whose MJD is already present replaces the existing one.
    pub(crate) fn extend(&mut self, records: Vec<EopRecord>) -> CoordResult<()> {
        let update = EopTable::new(records)?;
        let kept = self
            .table
            .records
            .iter()
            .filter(|r| !update.contains_mjd(r.mjd))
            .cloned();
        let merged = update.records.iter().cloned().chain(kept).collect();
        self.table = Arc::new(EopTable::new(merged)?);
        Ok(())
    }
}

fn lagrange_weights(mjd: f64, points: &[EopParameters; LAGRANGE_POINTS]) -> [f64; LAGRANGE_POINTS] {
    std::array::from_fn(|i| {
        let mut weight = 1.0;
        for (j, point) in points.iter().enumerate() {
            if i != j {
                weight *= (mjd - point.mjd) / (points[i].mjd - point.mjd);
            }
        }
        weight
    })
}

fn weighted_sum(
    weights: &[f64; LAGRANGE_POINTS],
    points: &[EopParameters; LAGRANGE_POINTS],
    value: impl Fn(&EopParameters) -> f64,
) -> f64 {
    let mut sum = 0.0;
    for (weight, point) in weights.iter().zip(points) {
        sum += value(point) * weight;
    }
    sum
}

// An interpolated optional column is present only when every record used has it.
fn flag_present_columns(mut params: EopParameters) -> EopParameters {
    params.flags.has_cip_offsets = params.dx.is_some() && params.dy.is_some();
    params.flags.has_pole_rates = params.xrt.is_some() && params.yrt.is_some();
    params
}

fn weighted_sum_opt(
    weights: &[f64; LAGRANGE_POINTS],
    points: &[EopParameters; LAGRANGE_POINTS],
    value: impl Fn(&EopParameters) -> Option<f64>,
) -> Option<f64> {
    let mut sum = 0.0;
    for (weight, point) in weights.iter().zip(points) {
        sum += value(point)? * weight;
    }
    Some(sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eop::EopProvider;

    fn create_test_records() -> Vec<EopRecord> {
        let mut records = Vec::new();

        for i in 0..5 {
            let mjd = 59945.0 + i as f64;
            let x_p = 0.1 + 0.001 * i as f64;
            let y_p = 0.2 + 0.002 * i as f64;
            let ut1_utc = 0.01 + 0.0001 * i as f64;
            let lod = 0.001 + 0.00001 * i as f64;

            let record = EopRecord::new(mjd, x_p, y_p, ut1_utc)
                .unwrap()
                .with_lod(lod)
                .unwrap();
            records.push(record);
        }

        records
    }

    fn interpolator(records: Vec<EopRecord>) -> EopInterpolator {
        EopInterpolator::new(records).unwrap()
    }

    #[test]
    fn test_linear_interpolation() {
        let params = interpolator(create_test_records()).get(59946.5).unwrap();
        assert_eq!(params.x_p, 0.1015);
        assert_eq!(params.y_p, 0.203);
        assert_eq!(params.mjd, 59946.5);
    }

    #[test]
    fn test_exact_match() {
        let params = interpolator(create_test_records()).get(59947.0).unwrap();
        assert_eq!(params.mjd, 59947.0);
        assert_eq!(params.x_p, 0.102);
        assert_eq!(params.y_p, 0.204);
    }

    #[test]
    fn test_linear_interpolation_cip_offsets() {
        let mut records = Vec::new();

        for i in 0..=1 {
            let mjd = 59945.0 + i as f64;
            let mut record = EopRecord::new(mjd, 0.1 + 0.001 * i as f64, 0.2, 0.01).unwrap();
            let dx = 1.0 + i as f64;
            let dy = -0.2 - 0.1 * i as f64;
            record = record.with_cip_offsets(dx, dy).unwrap();
            records.push(record);
        }

        let params = interpolator(records).get(59945.5).unwrap();
        assert_eq!(params.dx, Some(1.5));
        assert_eq!(params.dy, Some(-0.25));
        assert!(params.flags.has_cip_offsets);
    }

    // Every record but the last carries CIP offsets and pole rates.
    fn records_missing_columns_at_the_end(n: usize) -> Vec<EopRecord> {
        let record = |i: usize| EopRecord::new(59945.0 + i as f64, 0.1, 0.2, 0.01).unwrap();
        let full = |i: usize| {
            record(i)
                .with_cip_offsets(0.2, -0.1)
                .and_then(|r| r.with_pole_rates(0.1, 0.1))
                .unwrap()
        };
        (0..n)
            .map(|i| if i + 1 < n { full(i) } else { record(i) })
            .collect()
    }

    fn assert_columns_dropped(params: &EopParameters) {
        assert_eq!(
            (params.dx, params.dy, params.xrt, params.yrt),
            (None, None, None, None)
        );
        assert!(!params.flags.has_cip_offsets, "{:?}", params.flags);
        assert!(!params.flags.has_pole_rates, "{:?}", params.flags);
    }

    #[test]
    fn linear_flags_follow_the_dropped_columns() {
        let interpolator = interpolator(records_missing_columns_at_the_end(2));
        assert_columns_dropped(&interpolator.get(59945.5).unwrap());
    }

    #[test]
    fn lagrange_flags_follow_the_dropped_columns() {
        let interpolator = interpolator(records_missing_columns_at_the_end(5))
            .with_method(InterpolationMethod::Lagrange5);
        assert_columns_dropped(&interpolator.get(59946.5).unwrap());
    }

    #[test]
    fn test_lagrange_interpolation() {
        let interpolator =
            interpolator(create_test_records()).with_method(InterpolationMethod::Lagrange5);
        let params = interpolator.get(59947.0).unwrap();
        assert_eq!(params.x_p, 0.102);
        assert_eq!(params.y_p, 0.204);
    }

    #[test]
    fn test_lagrange_interpolation_cip_offsets() {
        let mut records = Vec::new();
        for i in 0..6 {
            let mjd = 59945.0 + i as f64;
            let mut record = EopRecord::new(mjd, 0.1 + 0.001 * i as f64, 0.2, 0.01).unwrap();
            let dx = 1.0 + 0.5 * i as f64;
            let dy = -0.2 + 0.05 * i as f64;
            record = record.with_cip_offsets(dx, dy).unwrap();
            records.push(record);
        }

        let interpolator = interpolator(records).with_method(InterpolationMethod::Lagrange5);
        let params = interpolator.get(59947.5).unwrap();
        assert_eq!(params.dx, Some(2.25));
        assert_eq!(params.dy, Some(-0.075));
        assert!(params.flags.has_cip_offsets);
    }

    // Bundled C04 at the 2016-12-31 leap: UT1-UTC is -0.4077697 s on MJD 57753 and
    // +0.591287 s on 57754. Expected values are the UT1-TAI interpolation, recomputed
    // outside the crate with the same operation order.
    #[test]
    fn interpolation_inside_leap_day_follows_ut1_minus_tai() {
        let provider = EopProvider::bundled().unwrap();
        assert_eq!(provider.get(57753.75).unwrap().ut1_utc, -0.408477175);
        let lagrange = provider.with_interpolation(InterpolationMethod::Lagrange5);
        assert_eq!(lagrange.get(57753.75).unwrap().ut1_utc, -0.4084669794433594);
    }

    #[test]
    fn lagrange_next_to_leap_does_not_ring() {
        let provider = EopProvider::bundled()
            .unwrap()
            .with_interpolation(InterpolationMethod::Lagrange5);
        assert_eq!(provider.get(57752.5).unwrap().ut1_utc, -0.40733358593749996);
    }

    #[test]
    fn non_finite_lookup_is_an_error() {
        let interpolator = interpolator(create_test_records());
        for mjd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(interpolator.get(mjd).is_err());
        }
    }

    #[test]
    fn non_finite_record_mjd_is_rejected() {
        let mut records = create_test_records();
        records[2].mjd = f64::NAN;
        assert!(EopInterpolator::new(records).is_err());
    }

    #[test]
    fn duplicate_record_mjd_is_rejected() {
        let mut records = create_test_records();
        records[3].mjd = records[1].mjd;
        let err = EopInterpolator::new(records).err().unwrap();
        assert!(err
            .to_string()
            .contains("Duplicate EOP records at MJD 59946"));
    }

    #[test]
    fn max_gap_must_be_positive() {
        for gap in [f64::NAN, 0.0, -1.0] {
            assert!(interpolator(create_test_records())
                .with_max_gap(gap)
                .is_err());
        }
    }

    #[test]
    fn test_out_of_range() {
        let interpolator = interpolator(create_test_records());
        assert!(interpolator.get(59944.0).is_err());
        assert!(interpolator.get(59950.0).is_err());
    }

    #[test]
    fn test_max_gap_enforcement() {
        let mut records = create_test_records();

        records[3].mjd = 59955.0;
        records[4].mjd = 59956.0;

        let interpolator = interpolator(records).with_max_gap(3.0).unwrap();
        assert!(interpolator.get(59950.0).is_err());
    }

    #[test]
    fn test_time_span() {
        assert_eq!(
            interpolator(create_test_records()).time_span(),
            (59945.0, 59949.0)
        );
    }

    #[test]
    fn test_empty_records() {
        let err = EopInterpolator::new(vec![]).err().unwrap();
        assert!(err.to_string().contains("No EOP records supplied"));
    }

    #[test]
    fn test_lagrange_insufficient_points() {
        let records = vec![
            EopRecord::new(59945.0, 0.1, 0.2, 0.01).unwrap(),
            EopRecord::new(59946.0, 0.101, 0.202, 0.0101).unwrap(),
        ];

        let interpolator = interpolator(records).with_method(InterpolationMethod::Lagrange5);
        assert!(interpolator.get(59945.5).is_err());
    }

    #[test]
    fn test_record_count() {
        assert_eq!(interpolator(create_test_records()).record_count(), 5);
    }

    #[test]
    fn test_find_center_index() {
        let interpolator = interpolator(create_test_records());
        assert_eq!(interpolator.find_center_index(59947.0), 2);
        assert_eq!(interpolator.find_center_index(59945.1), 0);
        assert_eq!(interpolator.find_center_index(59948.6), 4);
    }

    #[test]
    fn test_lagrange_edge_cases() {
        let mut records = Vec::new();
        for i in 0..10 {
            let mjd = 59945.0 + i as f64;
            records.push(EopRecord::new(mjd, 0.1 + 0.001 * i as f64, 0.2, 0.01).unwrap());
        }

        let interpolator = interpolator(records).with_method(InterpolationMethod::Lagrange5);
        assert_eq!(interpolator.lagrange_start(59945.5).unwrap(), 0);
        assert_eq!(interpolator.lagrange_start(59949.5).unwrap(), 2);
        assert_eq!(interpolator.lagrange_start(59953.5).unwrap(), 5);
    }

    #[test]
    fn extend_replaces_records_at_the_same_mjd() {
        let mut interpolator = interpolator(create_test_records());
        let update = vec![
            EopRecord::new(59949.0, 0.3, 0.3, -0.5).unwrap(),
            EopRecord::new(59950.0, 0.3, 0.3, -0.5).unwrap(),
        ];
        interpolator.extend(update).unwrap();
        assert_eq!(interpolator.record_count(), 6);
        assert_eq!(interpolator.get(59949.0).unwrap().ut1_utc, -0.5);
        assert_eq!(interpolator.get(59948.0).unwrap().ut1_utc, 0.0103);
    }
}
