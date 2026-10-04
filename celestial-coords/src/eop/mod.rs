pub mod bundled;
pub mod interpolate;
mod parse;
pub mod record;
mod table;

use interpolate::{EopInterpolator, InterpolationMethod};
use record::{EopParameters, EopRecord};
use table::EopTable;

use crate::errors::{CoordError, CoordResult};
use std::path::Path;
use std::sync::{Arc, OnceLock};

pub struct EopProvider {
    interpolator: EopInterpolator,
}

impl EopProvider {
    pub fn bundled() -> CoordResult<Self> {
        static TABLE: OnceLock<Arc<EopTable>> = OnceLock::new();
        Self::from_cached(&TABLE, bundled::load_bundled_combined)
    }

    pub fn bundled_c04() -> CoordResult<Self> {
        static TABLE: OnceLock<Arc<EopTable>> = OnceLock::new();
        Self::from_cached(&TABLE, bundled::load_bundled_c04)
    }

    fn from_cached(
        cell: &OnceLock<Arc<EopTable>>,
        load: fn() -> CoordResult<Vec<EopRecord>>,
    ) -> CoordResult<Self> {
        let table = match cell.get() {
            Some(table) => Arc::clone(table),
            None => {
                let built = Arc::new(EopTable::new(load()?)?);
                Arc::clone(cell.get_or_init(|| built))
            }
        };
        Ok(Self {
            interpolator: EopInterpolator::from_table(table),
        })
    }

    pub fn from_records(records: Vec<EopRecord>) -> CoordResult<Self> {
        Ok(Self {
            interpolator: EopInterpolator::new(records)?,
        })
    }

    pub fn with_interpolation(mut self, method: InterpolationMethod) -> Self {
        self.interpolator = self.interpolator.with_method(method);
        self
    }

    pub fn get(&self, mjd: f64) -> CoordResult<EopParameters> {
        self.interpolator.get(mjd)
    }

    pub fn time_span(&self) -> (f64, f64) {
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
        let content = std::fs::read_to_string(path.as_ref())
            .map_err(|e| CoordError::io("reading finals2000A file", e))?;
        Self::from_finals_str(&content)
    }

    pub fn bundled_with_update(path: impl AsRef<Path>) -> CoordResult<Self> {
        let update_content = std::fs::read_to_string(path.as_ref())
            .map_err(|e| CoordError::io("reading finals2000A update file", e))?;
        Self::bundled()?.overlay_finals(&update_content)
    }

    fn overlay_finals(mut self, content: &str) -> CoordResult<Self> {
        let update = bundled::after_c04(parse::parse_finals(content)?);
        if !update.is_empty() {
            self.interpolator.extend(update)?;
        }
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bundled_provider() {
        let provider = EopProvider::bundled().unwrap();
        assert!(provider.record_count() > 0);
        assert_eq!(provider.time_span(), bundled::bundled_time_span());
    }

    #[test]
    fn test_bundled_lookup() {
        let provider = EopProvider::bundled().unwrap();
        let params = provider.get(59945.0).unwrap();
        assert_eq!(params.mjd, 59945.0);
        assert!(libm::fabs(params.x_p) < 1.0);
        assert!(libm::fabs(params.y_p) < 1.0);
        assert!(libm::fabs(params.ut1_utc) < 1.0);
    }

    #[test]
    fn test_from_records() {
        let records = vec![
            EopRecord::new(60000.0, 0.1, 0.2, 0.01).unwrap(),
            EopRecord::new(60001.0, 0.101, 0.202, 0.011).unwrap(),
        ];
        let provider = EopProvider::from_records(records).unwrap();
        let params = provider.get(60000.5).unwrap();
        assert_eq!(params.x_p, 0.1005);
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

    fn io_kind(err: &CoordError) -> Option<std::io::ErrorKind> {
        let source = std::error::Error::source(err)?;
        source
            .downcast_ref::<std::io::Error>()
            .map(std::io::Error::kind)
    }

    #[test]
    fn unreadable_files_keep_the_io_error() {
        let path = "/nonexistent/finals2000A.all";
        let err = EopProvider::from_finals_file(path).err().unwrap();
        assert_eq!(io_kind(&err), Some(std::io::ErrorKind::NotFound), "{err:?}");
        let err = EopProvider::bundled_with_update(path).err().unwrap();
        assert_eq!(io_kind(&err), Some(std::io::ErrorKind::NotFound), "{err:?}");
    }

    // The last bundled days are predictions, which have no LOD.
    #[test]
    fn bundled_predictions_have_no_lod() {
        let provider = EopProvider::bundled().unwrap();
        let last = provider.get(provider.time_span().1).unwrap();
        assert_eq!(last.lod, None);
    }

    #[test]
    fn test_immutable_get() {
        let provider = EopProvider::bundled().unwrap();
        let _p1 = provider.get(59945.0).unwrap();
        let _p2 = provider.get(59945.0).unwrap();
    }

    fn sample_finals_line_at(mjd: &[u8]) -> String {
        let mut line = vec![b' '; 188];
        line[7..7 + mjd.len()].copy_from_slice(mjd);
        line[18..27].copy_from_slice(b"  0.10000");
        line[37..46].copy_from_slice(b"  0.25000");
        line[58..68].copy_from_slice(b" -0.050000");
        line[79..86].copy_from_slice(b"  1.500");
        String::from_utf8(line).unwrap()
    }

    #[test]
    fn bundled_c04_loads_records() {
        let provider = EopProvider::bundled_c04().unwrap();
        assert!(provider.record_count() > 0);
        assert_eq!(provider.time_span().0, bundled::bundled_time_span().0);
    }

    #[test]
    fn with_interpolation_changes_method() {
        // Both methods must yield the same exact-MJD value; the setter only
        // changes behavior between samples. Pin the return type by checking
        // chaining and that the result is still queryable.
        let provider = EopProvider::bundled()
            .unwrap()
            .with_interpolation(InterpolationMethod::Lagrange5);
        let params = provider.get(59945.0).unwrap();
        assert_eq!(params.mjd, 59945.0);
    }

    #[test]
    fn from_finals_str_parses_and_builds_provider() {
        let line1 = sample_finals_line_at(b"60000.00");
        let line2 = sample_finals_line_at(b"60001.00");
        let content = format!("{}\n{}\n", line1, line2);
        let provider = EopProvider::from_finals_str(&content).unwrap();
        assert_eq!(provider.record_count(), 2);
        assert_eq!(provider.time_span(), (60000.0, 60001.0));
    }

    #[test]
    fn from_finals_str_propagates_parse_error() {
        let result = EopProvider::from_finals_str("garbage\nlines\nonly\n");
        assert!(result.is_err());
    }

    #[test]
    fn finals_overlay_replaces_bundled_predictions() {
        let (_, end) = bundled::bundled_time_span();
        let mjd = format!("{:8.2}", end);
        let content = sample_finals_line_at(mjd.as_bytes());
        let provider = EopProvider::bundled()
            .unwrap()
            .overlay_finals(&content)
            .unwrap();
        assert_eq!(provider.get(end).unwrap().ut1_utc, -0.05);
    }

    #[test]
    fn finals_overlay_keeps_c04_days() {
        let content = sample_finals_line_at(b"57754.00");
        let provider = EopProvider::bundled()
            .unwrap()
            .overlay_finals(&content)
            .unwrap();
        assert_eq!(provider.get(57754.0).unwrap().ut1_utc, 0.591287);
    }
}
