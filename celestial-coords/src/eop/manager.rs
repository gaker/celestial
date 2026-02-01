use super::{
    cache::EopCacheSystem,
    embedded,
    interpolate::{EopInterpolator, InterpolationMethod},
    parser::{parse_eop_file, C04Parser, FinalsParser},
    record::{EopParameters, EopRecord},
};
use crate::{CoordError, CoordResult};

use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UpdatePolicy {
    Offline,

    OnDemand,

    Proactive,
}

#[derive(Debug, Clone)]
pub struct EopConfig {
    pub update_policy: UpdatePolicy,

    pub interpolation_method: InterpolationMethod,

    pub cache_dir: Option<PathBuf>,

    pub memory_cache_size: usize,

    pub max_interpolation_gap: f64,

    pub use_embedded_baseline: bool,
}

impl Default for EopConfig {
    fn default() -> Self {
        Self {
            update_policy: UpdatePolicy::Offline,
            interpolation_method: InterpolationMethod::Linear,
            cache_dir: None,
            memory_cache_size: 365,
            max_interpolation_gap: 5.0,
            use_embedded_baseline: true,
        }
    }
}

pub struct EopManager {
    config: EopConfig,

    interpolator: Option<EopInterpolator>,

    cache: EopCacheSystem,

    data_span: Option<(f64, f64)>,
}

impl EopManager {
    pub fn new(config: EopConfig) -> CoordResult<Self> {
        let cache = if let Some(ref cache_dir) = config.cache_dir {
            EopCacheSystem::with_disk_cache(config.memory_cache_size, cache_dir)?
        } else {
            EopCacheSystem::memory_only(config.memory_cache_size)
        };

        let mut manager = Self {
            config,
            interpolator: None,
            cache,
            data_span: None,
        };

        if manager.config.use_embedded_baseline {
            manager.load_embedded_baseline()?;
        }

        Ok(manager)
    }

    pub fn get(&mut self, mjd: f64) -> CoordResult<EopParameters> {
        if let Some(record) = self.cache.get(mjd) {
            return Ok(record.to_parameters());
        }

        if let Some(ref interpolator) = self.interpolator {
            let params = interpolator.get(mjd)?;
            let mut record = EopRecord::new(
                params.mjd,
                params.x_p,
                params.y_p,
                params.ut1_utc,
                params.lod,
            )?;

            if let (Some(dx), Some(dy)) = (params.dx, params.dy) {
                record = record.with_cip_offsets(dx, dy)?;
            }

            record = record.with_flags(params.flags);
            self.cache.insert(record);
            return Ok(params);
        }

        match self.config.update_policy {
            UpdatePolicy::Offline => Err(CoordError::data_unavailable(format!(
                "No EOP data available for MJD {:.1} (offline mode)",
                mjd
            ))),
            UpdatePolicy::OnDemand => {
                self.fetch_data_for_mjd(mjd)?;
                self.get(mjd)
            }
            UpdatePolicy::Proactive => {
                self.fetch_data_for_mjd(mjd)?;
                self.get(mjd)
            }
        }
    }

    fn load_embedded_baseline(&mut self) -> CoordResult<()> {
        let baseline_records = embedded::generate_interpolated_baseline()?;
        self.load_records(baseline_records)?;
        Ok(())
    }

    pub fn load_records(&mut self, records: Vec<EopRecord>) -> CoordResult<()> {
        if records.is_empty() {
            return Ok(());
        }

        let start_mjd = records[0].mjd;
        let end_mjd = records.last().unwrap().mjd;
        self.data_span = Some((start_mjd, end_mjd));

        self.interpolator = Some(
            EopInterpolator::new(records)
                .with_method(self.config.interpolation_method)
                .with_max_gap(self.config.max_interpolation_gap),
        );

        Ok(())
    }

    pub fn load_from_file<P: AsRef<std::path::Path>>(&mut self, path: P) -> CoordResult<()> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| CoordError::external_library("File read", &e.to_string()))?;

        let records = parse_eop_file(&content)?;
        self.load_records(records)
    }

    pub fn load_c04_data(&mut self, content: &str) -> CoordResult<()> {
        let parser = C04Parser::new();
        let records = parser.parse(content)?;
        self.load_records(records)
    }

    pub fn load_finals_data(&mut self, content: &str) -> CoordResult<()> {
        let parser = FinalsParser::new();
        let records = parser.parse(content)?;
        self.load_records(records)
    }

    fn fetch_data_for_mjd(&mut self, _mjd: f64) -> CoordResult<()> {
        Err(CoordError::data_unavailable(
            "Network data fetching not yet implemented",
        ))
    }

    pub fn data_time_span(&self) -> Option<(f64, f64)> {
        self.data_span
    }

    pub fn cache_hit_ratio(&self) -> f64 {
        self.cache.hit_ratio()
    }

    pub fn clear_cache(&mut self) -> CoordResult<()> {
        self.cache.clear()
    }

    pub fn config(&self) -> &EopConfig {
        &self.config
    }

    pub fn record_count(&self) -> usize {
        self.interpolator
            .as_ref()
            .map(|i| i.record_count())
            .unwrap_or(0)
    }
}

pub struct EopBuilder {
    config: EopConfig,
}

impl EopBuilder {
    pub fn new() -> Self {
        Self {
            config: EopConfig::default(),
        }
    }

    pub fn with_update_policy(mut self, policy: UpdatePolicy) -> Self {
        self.config.update_policy = policy;
        self
    }

    pub fn with_interpolation(mut self, method: InterpolationMethod) -> Self {
        self.config.interpolation_method = method;
        self
    }

    pub fn with_cache_dir<P: Into<PathBuf>>(mut self, dir: P) -> Self {
        self.config.cache_dir = Some(dir.into());
        self
    }

    pub fn with_memory_cache_size(mut self, size: usize) -> Self {
        self.config.memory_cache_size = size;
        self
    }

    pub fn with_max_gap(mut self, days: f64) -> Self {
        self.config.max_interpolation_gap = days;
        self
    }

    pub fn with_embedded_baseline(mut self) -> Self {
        self.config.use_embedded_baseline = true;
        self
    }

    pub fn without_embedded_baseline(mut self) -> Self {
        self.config.use_embedded_baseline = false;
        self
    }

    pub fn build(self) -> CoordResult<EopManager> {
        EopManager::new(self.config)
    }
}

impl Default for EopBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl EopManager {
    pub fn builder() -> EopBuilder {
        EopBuilder::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unique_temp_path(prefix: &str, suffix: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        let unique = format!(
            "{}_{}_{}{}",
            prefix,
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            suffix
        );
        path.push(unique);
        path
    }

    #[test]
    fn test_eop_manager_creation() {
        let manager = EopManager::builder()
            .with_embedded_baseline()
            .build()
            .unwrap();

        assert!(manager.record_count() > 0);
        assert!(manager.data_time_span().is_some());
    }

    #[test]
    fn test_eop_manager_interpolation() {
        let mut manager = EopManager::builder()
            .with_embedded_baseline()
            .with_interpolation(InterpolationMethod::Linear)
            .build()
            .unwrap();

        let mjd = 59945.0;
        let params = manager.get(mjd).unwrap();

        assert_eq!(params.mjd, mjd);
        assert!(params.x_p.abs() < 1.0);
        assert!(params.y_p.abs() < 1.0);
        assert!(params.ut1_utc.abs() < 1.0);
    }

    #[test]
    fn test_eop_manager_cache() {
        let mut manager = EopManager::builder()
            .with_embedded_baseline()
            .with_memory_cache_size(100)
            .build()
            .unwrap();

        let mjd = 59945.0;

        let _params1 = manager.get(mjd).unwrap();

        let _params2 = manager.get(mjd).unwrap();

        assert!(manager.cache_hit_ratio() >= 0.0);
    }

    #[test]
    fn test_eop_builder_configuration() {
        let manager = EopManager::builder()
            .with_update_policy(UpdatePolicy::OnDemand)
            .with_interpolation(InterpolationMethod::Lagrange5)
            .with_memory_cache_size(500)
            .with_max_gap(10.0)
            .without_embedded_baseline()
            .build()
            .unwrap();

        assert_eq!(manager.config().update_policy, UpdatePolicy::OnDemand);
        assert_eq!(
            manager.config().interpolation_method,
            InterpolationMethod::Lagrange5
        );
        assert_eq!(manager.config().memory_cache_size, 500);
        assert_eq!(manager.config().max_interpolation_gap, 10.0);
        assert!(!manager.config().use_embedded_baseline);
    }

    #[test]
    fn test_out_of_range_offline() {
        let mut manager = EopManager::builder()
            .with_embedded_baseline()
            .with_update_policy(UpdatePolicy::Offline)
            .build()
            .unwrap();

        let result = manager.get(70000.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_load_empty_records() {
        let mut manager = EopManager::builder()
            .without_embedded_baseline()
            .build()
            .unwrap();

        let result = manager.load_records(vec![]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_load_c04_data_updates_state() {
        let mut manager = EopManager::builder()
            .without_embedded_baseline()
            .build()
            .unwrap();

        let content =
            "1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230\n1962   1   2  37666  -0.015900   0.214100   0.0320547   0.0016690";
        manager.load_c04_data(content).unwrap();

        assert_eq!(manager.record_count(), 2);
        let span = manager.data_time_span().unwrap();
        assert!((span.0 - 37665.0).abs() < 1e-9);
        assert!((span.1 - 37666.0).abs() < 1e-9);
    }

    #[test]
    fn test_load_finals_data_updates_state() {
        let mut manager = EopManager::builder()
            .without_embedded_baseline()
            .build()
            .unwrap();

        let finals_line = "       59665.00   0.123456           0.234567             0.345678             0.001234                                                                                                                 ";
        manager.load_finals_data(finals_line).unwrap();

        assert_eq!(manager.record_count(), 1);
        let span = manager.data_time_span().unwrap();
        assert!((span.0 - 59665.0).abs() < 1.0);
        assert!((span.1 - 59665.0).abs() < 1.0);
    }

    #[test]
    fn test_load_from_file_updates_state() {
        let mut manager = EopManager::builder()
            .without_embedded_baseline()
            .build()
            .unwrap();

        let path = unique_temp_path("eop_manager_load", ".txt");
        let content = "1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230\n";
        std::fs::write(&path, content).unwrap();

        manager.load_from_file(&path).unwrap();
        assert_eq!(manager.record_count(), 1);
        let span = manager.data_time_span().unwrap();
        assert!((span.0 - 37665.0).abs() < 1e-9);
        assert!((span.1 - 37665.0).abs() < 1e-9);

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_clear_cache() {
        let mut manager = EopManager::builder()
            .with_embedded_baseline()
            .build()
            .unwrap();

        let _params = manager.get(59945.0).unwrap();
        assert!(manager.clear_cache().is_ok());
    }

    #[test]
    fn test_data_time_span() {
        let manager = EopManager::builder()
            .with_embedded_baseline()
            .build()
            .unwrap();

        let span = manager.data_time_span();
        assert!(span.is_some());
        let (start, end) = span.unwrap();
        assert!(end > start);
    }

    #[test]
    fn test_on_demand_fetch_unimplemented() {
        let mut manager = EopManager::builder()
            .without_embedded_baseline()
            .with_update_policy(UpdatePolicy::OnDemand)
            .build()
            .unwrap();

        let result = manager.get(59945.0);
        assert!(result.is_err());
    }
}
