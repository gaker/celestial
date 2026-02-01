use super::record::{EopParameters, EopRecord};
use crate::{CoordError, CoordResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterpolationMethod {
    Linear,

    Lagrange5,

    CubicSpline,
}

pub struct EopInterpolator {
    records: Vec<EopRecord>,

    method: InterpolationMethod,

    max_gap_days: f64,
}

impl EopInterpolator {
    pub fn new(mut records: Vec<EopRecord>) -> Self {
        records.sort_by(|a, b| a.mjd.partial_cmp(&b.mjd).unwrap());

        Self {
            records,
            method: InterpolationMethod::Linear,
            max_gap_days: 5.0,
        }
    }

    pub fn with_method(mut self, method: InterpolationMethod) -> Self {
        self.method = method;
        self
    }

    pub fn with_max_gap(mut self, max_gap_days: f64) -> Self {
        self.max_gap_days = max_gap_days;
        self
    }

    pub fn get(&self, mjd: f64) -> CoordResult<EopParameters> {
        if self.records.is_empty() {
            return Err(CoordError::data_unavailable(
                "No EOP records available for interpolation",
            ));
        }

        if let Some(record) = self.records.iter().find(|r| (r.mjd - mjd).abs() < 1e-10) {
            return Ok(record.to_parameters());
        }

        let (before_idx, after_idx) = self.find_interpolation_interval(mjd)?;

        let gap = self.records[after_idx].mjd - self.records[before_idx].mjd;
        if gap > self.max_gap_days {
            return Err(CoordError::data_unavailable(format!(
                "Gap of {:.1} days exceeds maximum interpolation gap of {:.1} days",
                gap, self.max_gap_days
            )));
        }

        match self.method {
            InterpolationMethod::Linear => self.linear_interpolate(mjd, before_idx, after_idx),
            InterpolationMethod::Lagrange5 => self.lagrange_interpolate(mjd, 5),
            InterpolationMethod::CubicSpline => self.linear_interpolate(mjd, before_idx, after_idx),
        }
    }

    fn find_interpolation_interval(&self, mjd: f64) -> CoordResult<(usize, usize)> {
        if mjd < self.records[0].mjd {
            return Err(CoordError::data_unavailable(format!(
                "MJD {:.1} is before first available record (MJD {:.1})",
                mjd, self.records[0].mjd
            )));
        }

        if mjd > self.records.last().unwrap().mjd {
            return Err(CoordError::data_unavailable(format!(
                "MJD {:.1} is after last available record (MJD {:.1})",
                mjd,
                self.records.last().unwrap().mjd
            )));
        }

        let mut left = 0;
        let mut right = self.records.len() - 1;

        while right - left > 1 {
            let mid = (left + right) / 2;
            if self.records[mid].mjd <= mjd {
                left = mid;
            } else {
                right = mid;
            }
        }

        Ok((left, right))
    }

    fn linear_interpolate(
        &self,
        mjd: f64,
        before_idx: usize,
        after_idx: usize,
    ) -> CoordResult<EopParameters> {
        let r1 = &self.records[before_idx];
        let r2 = &self.records[after_idx];

        let t = (mjd - r1.mjd) / (r2.mjd - r1.mjd);

        let p1 = r1.to_parameters();
        let p2 = r2.to_parameters();

        let x_p = p1.x_p + t * (p2.x_p - p1.x_p);
        let y_p = p1.y_p + t * (p2.y_p - p1.y_p);
        let ut1_utc = p1.ut1_utc + t * (p2.ut1_utc - p1.ut1_utc);
        let lod = p1.lod + t * (p2.lod - p1.lod);

        let dx = match (p1.dx, p2.dx) {
            (Some(dx1), Some(dx2)) => Some(dx1 + t * (dx2 - dx1)),
            _ => None,
        };

        let dy = match (p1.dy, p2.dy) {
            (Some(dy1), Some(dy2)) => Some(dy1 + t * (dy2 - dy1)),
            _ => None,
        };

        let mut params = EopParameters {
            mjd,
            x_p,
            y_p,
            ut1_utc,
            lod,
            dx,
            dy,
            s_prime: 0.0,
            flags: p1.flags,
        };

        params.compute_s_prime();

        Ok(params)
    }

    fn lagrange_interpolate(&self, mjd: f64, n: usize) -> CoordResult<EopParameters> {
        if n > self.records.len() {
            return Err(CoordError::invalid_coordinate(format!(
                "Not enough records for {}-point Lagrange interpolation",
                n
            )));
        }

        let center_idx = self.find_center_index(mjd)?;
        let half_n = n / 2;

        let start_idx = center_idx.saturating_sub(half_n);
        let end_idx = (start_idx + n).min(self.records.len());
        let actual_start = end_idx.saturating_sub(n);

        if end_idx - actual_start < n {
            return Err(CoordError::data_unavailable(
                "Insufficient records for Lagrange interpolation",
            ));
        }

        let points: Vec<_> = self.records[actual_start..actual_start + n]
            .iter()
            .map(|r| r.to_parameters())
            .collect();

        let x_p = self.lagrange_interpolate_value(mjd, &points, |p| p.x_p);
        let y_p = self.lagrange_interpolate_value(mjd, &points, |p| p.y_p);
        let ut1_utc = self.lagrange_interpolate_value(mjd, &points, |p| p.ut1_utc);
        let lod = self.lagrange_interpolate_value(mjd, &points, |p| p.lod);

        let has_dx = points.iter().all(|p| p.dx.is_some());
        let has_dy = points.iter().all(|p| p.dy.is_some());

        let dx = if has_dx {
            Some(self.lagrange_interpolate_value(mjd, &points, |p| p.dx.unwrap()))
        } else {
            None
        };

        let dy = if has_dy {
            Some(self.lagrange_interpolate_value(mjd, &points, |p| p.dy.unwrap()))
        } else {
            None
        };

        let mut params = EopParameters {
            mjd,
            x_p,
            y_p,
            ut1_utc,
            lod,
            dx,
            dy,
            s_prime: 0.0,
            flags: points[0].flags,
        };

        params.compute_s_prime();

        Ok(params)
    }

    fn find_center_index(&self, mjd: f64) -> CoordResult<usize> {
        let mut best_idx = 0;
        let mut best_diff = (self.records[0].mjd - mjd).abs();

        for (idx, record) in self.records.iter().enumerate() {
            let diff = (record.mjd - mjd).abs();
            if diff < best_diff {
                best_diff = diff;
                best_idx = idx;
            }
        }

        Ok(best_idx)
    }

    fn lagrange_interpolate_value<F>(&self, mjd: f64, points: &[EopParameters], extract: F) -> f64
    where
        F: Fn(&EopParameters) -> f64,
    {
        let n = points.len();
        let mut result = 0.0;

        for i in 0..n {
            let yi = extract(&points[i]);
            let xi = points[i].mjd;

            let mut li = 1.0;
            for (j, point) in points.iter().enumerate().take(n) {
                if i != j {
                    let xj = point.mjd;
                    li *= (mjd - xj) / (xi - xj);
                }
            }

            result += yi * li;
        }

        result
    }

    pub fn time_span(&self) -> Option<(f64, f64)> {
        if self.records.is_empty() {
            None
        } else {
            Some((self.records[0].mjd, self.records.last().unwrap().mjd))
        }
    }

    pub fn record_count(&self) -> usize {
        self.records.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_records() -> Vec<EopRecord> {
        let mut records = Vec::new();

        for i in 0..5 {
            let mjd = 59945.0 + i as f64;
            let x_p = 0.1 + 0.001 * i as f64;
            let y_p = 0.2 + 0.002 * i as f64;
            let ut1_utc = 0.01 + 0.0001 * i as f64;
            let lod = 0.001 + 0.00001 * i as f64;

            let record = EopRecord::new(mjd, x_p, y_p, ut1_utc, lod).unwrap();
            records.push(record);
        }

        records
    }

    #[test]
    fn test_linear_interpolation() {
        let records = create_test_records();
        let interpolator = EopInterpolator::new(records);

        let mjd = 59946.5;
        let params = interpolator.get(mjd).unwrap();

        let expected_x_p = (0.101 + 0.102) / 2.0;
        let expected_y_p = (0.202 + 0.204) / 2.0;

        assert!((params.x_p - expected_x_p).abs() < 1e-10);
        assert!((params.y_p - expected_y_p).abs() < 1e-10);
        assert_eq!(params.mjd, mjd);
    }

    #[test]
    fn test_exact_match() {
        let records = create_test_records();
        let interpolator = EopInterpolator::new(records);

        let mjd = 59947.0;
        let params = interpolator.get(mjd).unwrap();

        assert_eq!(params.mjd, mjd);
        assert!((params.x_p - 0.102).abs() < 1e-10);
        assert!((params.y_p - 0.204).abs() < 1e-10);
    }

    #[test]
    fn test_linear_interpolation_cip_offsets() {
        let mut records = Vec::new();

        for i in 0..=1 {
            let mjd = 59945.0 + i as f64;
            let mut record = EopRecord::new(mjd, 0.1 + 0.001 * i as f64, 0.2, 0.01, 0.001).unwrap();
            let dx = 1.0 + i as f64;
            let dy = -0.2 - 0.1 * i as f64;
            record = record.with_cip_offsets(dx, dy).unwrap();
            records.push(record);
        }

        let interpolator = EopInterpolator::new(records);
        let params = interpolator.get(59945.5).unwrap();

        assert!(params.dx.is_some());
        assert!(params.dy.is_some());
        assert!((params.dx.unwrap() - 1.5).abs() < 1e-10);
        assert!((params.dy.unwrap() - (-0.25)).abs() < 1e-10);
        assert!(params.flags.has_cip_offsets);
    }

    #[test]
    fn test_lagrange_interpolation() {
        let records = create_test_records();
        let interpolator =
            EopInterpolator::new(records).with_method(InterpolationMethod::Lagrange5);

        let mjd = 59947.0;
        let params = interpolator.get(mjd).unwrap();

        assert!((params.x_p - 0.102).abs() < 1e-10);
        assert!((params.y_p - 0.204).abs() < 1e-10);
    }

    #[test]
    fn test_lagrange_interpolation_cip_offsets() {
        let mut records = Vec::new();
        for i in 0..6 {
            let mjd = 59945.0 + i as f64;
            let mut record = EopRecord::new(mjd, 0.1 + 0.001 * i as f64, 0.2, 0.01, 0.001).unwrap();
            let dx = 1.0 + 0.5 * i as f64;
            let dy = -0.2 + 0.05 * i as f64;
            record = record.with_cip_offsets(dx, dy).unwrap();
            records.push(record);
        }

        let interpolator =
            EopInterpolator::new(records).with_method(InterpolationMethod::Lagrange5);

        let target_mjd = 59947.5;
        let params = interpolator.get(target_mjd).unwrap();

        let expected_dx = 1.0 + 0.5 * (target_mjd - 59945.0);
        let expected_dy = -0.2 + 0.05 * (target_mjd - 59945.0);

        assert!(params.dx.is_some());
        assert!(params.dy.is_some());
        assert!((params.dx.unwrap() - expected_dx).abs() < 1e-10);
        assert!((params.dy.unwrap() - expected_dy).abs() < 1e-10);
        assert!(params.flags.has_cip_offsets);
    }

    #[test]
    fn test_out_of_range() {
        let records = create_test_records();
        let interpolator = EopInterpolator::new(records);

        let result = interpolator.get(59944.0);
        assert!(result.is_err());

        let result = interpolator.get(59950.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_max_gap_enforcement() {
        let mut records = create_test_records();

        records[3].mjd = 59955.0;
        records[4].mjd = 59956.0;

        let interpolator = EopInterpolator::new(records).with_max_gap(3.0);

        let result = interpolator.get(59950.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_time_span() {
        let records = create_test_records();
        let interpolator = EopInterpolator::new(records);

        let (start, end) = interpolator.time_span().unwrap();
        assert_eq!(start, 59945.0);
        assert_eq!(end, 59949.0);
    }

    #[test]
    fn test_empty_records() {
        let interpolator = EopInterpolator::new(vec![]);
        let result = interpolator.get(59945.0);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("No EOP records available"));
    }

    #[test]
    fn test_empty_records_time_span() {
        let interpolator = EopInterpolator::new(vec![]);
        assert_eq!(interpolator.time_span(), None);
    }

    #[test]
    fn test_cubic_spline_fallback_to_linear() {
        let records = create_test_records();
        let interpolator =
            EopInterpolator::new(records).with_method(InterpolationMethod::CubicSpline);

        let mjd = 59946.5;
        let params = interpolator.get(mjd).unwrap();

        let expected_x_p = (0.101 + 0.102) / 2.0;
        assert!((params.x_p - expected_x_p).abs() < 1e-10);
    }

    #[test]
    fn test_lagrange_insufficient_points() {
        let records = vec![
            EopRecord::new(59945.0, 0.1, 0.2, 0.01, 0.001).unwrap(),
            EopRecord::new(59946.0, 0.101, 0.202, 0.0101, 0.0011).unwrap(),
        ];

        let interpolator =
            EopInterpolator::new(records).with_method(InterpolationMethod::Lagrange5);

        let result = interpolator.get(59945.5);
        assert!(result.is_err());
    }

    #[test]
    fn test_record_count() {
        let records = create_test_records();
        let interpolator = EopInterpolator::new(records);
        assert_eq!(interpolator.record_count(), 5);

        let empty = EopInterpolator::new(vec![]);
        assert_eq!(empty.record_count(), 0);
    }

    #[test]
    fn test_find_center_index() {
        let records = create_test_records();
        let interpolator = EopInterpolator::new(records);

        let center = interpolator.find_center_index(59947.0).unwrap();
        assert_eq!(center, 2);

        let edge = interpolator.find_center_index(59945.1).unwrap();
        assert_eq!(edge, 0);
    }

    #[test]
    fn test_lagrange_edge_cases() {
        let mut records = Vec::new();
        for i in 0..10 {
            let mjd = 59945.0 + i as f64;
            records.push(EopRecord::new(mjd, 0.1 + 0.001 * i as f64, 0.2, 0.01, 0.001).unwrap());
        }

        let interpolator =
            EopInterpolator::new(records).with_method(InterpolationMethod::Lagrange5);

        let near_start = interpolator.get(59946.0).unwrap();
        assert!(near_start.x_p > 0.1);

        let near_end = interpolator.get(59953.0).unwrap();
        assert!(near_end.x_p > 0.1);
    }
}
