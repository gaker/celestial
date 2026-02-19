pub mod bundled;
pub mod interpolate;
pub mod parse;
pub mod record;

pub use record::{EopParameters, EopRecord};

use interpolate::{EopInterpolator, InterpolationMethod};

use crate::{CoordError, CoordResult};
use std::path::Path;

pub struct EopProvider {
    interpolator: EopInterpolator,
}

impl EopProvider {
    pub fn bundled() -> CoordResult<Self> {
        let records = bundled::load_bundled_combined()?;
        Self::from_records(records)
    }

    pub fn bundled_c04() -> CoordResult<Self> {
        let records = bundled::load_bundled_c04()?;
        Self::from_records(records)
    }

    pub fn from_records(records: Vec<EopRecord>) -> CoordResult<Self> {
        if records.is_empty() {
            return Err(CoordError::data_unavailable(
                "Cannot create EopProvider with empty records",
            ));
        }
        Ok(Self {
            interpolator: EopInterpolator::new(records),
        })
    }

    pub fn with_interpolation(mut self, method: InterpolationMethod) -> Self {
        self.interpolator = self.interpolator.with_method(method);
        self
    }

    pub fn get(&self, mjd: f64) -> CoordResult<EopParameters> {
        self.interpolator.get(mjd)
    }

    pub fn time_span(&self) -> Option<(f64, f64)> {
        self.interpolator.time_span()
    }

    pub fn record_count(&self) -> usize {
        self.interpolator.record_count()
    }

    pub fn from_finals_str(content: &str) -> CoordResult<Self> {
        let records = parse::parse_finals(content)?;
        Self::from_records(records)
    }

    pub fn from_finals_file(path: impl AsRef<Path>) -> CoordResult<Self> {
        let content = std::fs::read_to_string(path.as_ref()).map_err(|e| {
            CoordError::external_library("reading finals2000A file", &e.to_string())
        })?;
        Self::from_finals_str(&content)
    }

    pub fn bundled_with_update(path: impl AsRef<Path>) -> CoordResult<Self> {
        let mut provider = Self::bundled()?;
        let update_content = std::fs::read_to_string(path.as_ref()).map_err(|e| {
            CoordError::external_library("reading finals2000A update file", &e.to_string())
        })?;
        let update_records = parse::parse_finals(&update_content)?;
        provider.interpolator.extend(update_records);
        Ok(provider)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bundled_provider() {
        let provider = EopProvider::bundled().unwrap();
        assert!(provider.record_count() > 0);
        assert!(provider.time_span().is_some());
    }

    #[test]
    fn test_bundled_lookup() {
        let provider = EopProvider::bundled().unwrap();
        let params = provider.get(59945.0).unwrap();
        assert_eq!(params.mjd, 59945.0);
        assert!(params.x_p.abs() < 1.0);
        assert!(params.y_p.abs() < 1.0);
        assert!(params.ut1_utc.abs() < 1.0);
    }

    #[test]
    fn test_from_records() {
        let records = vec![
            EopRecord::new(60000.0, 0.1, 0.2, 0.01, 0.001).unwrap(),
            EopRecord::new(60001.0, 0.101, 0.202, 0.011, 0.001).unwrap(),
        ];
        let provider = EopProvider::from_records(records).unwrap();
        let params = provider.get(60000.5).unwrap();
        assert!((params.x_p - 0.1005).abs() < 1e-7);
    }

    #[test]
    fn test_empty_records_rejected() {
        let result = EopProvider::from_records(vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_out_of_range() {
        let provider = EopProvider::bundled().unwrap();
        assert!(provider.get(70000.0).is_err());
    }

    #[test]
    fn test_immutable_get() {
        let provider = EopProvider::bundled().unwrap();
        let _p1 = provider.get(59945.0).unwrap();
        let _p2 = provider.get(59945.0).unwrap();
    }
}
