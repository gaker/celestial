use std::collections::HashMap;

use crate::error::{WcsError, WcsResult};

pub trait KeywordProvider {
    fn get_string(&self, key: &str) -> Option<String>;
    fn get_float(&self, key: &str) -> Option<f64>;
    fn get_int(&self, key: &str) -> Option<i64>;

    fn require_float(&self, key: &str) -> WcsResult<f64> {
        self.get_float(key)
            .ok_or_else(|| WcsError::missing_keyword(key))
    }

    fn require_string(&self, key: &str) -> WcsResult<String> {
        self.get_string(key)
            .ok_or_else(|| WcsError::missing_keyword(key))
    }
}

#[derive(Debug, Clone, Default)]
pub struct KeywordMap {
    strings: HashMap<String, String>,
    floats: HashMap<String, f64>,
    ints: HashMap<String, i64>,
}

impl KeywordMap {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_string(&mut self, key: impl Into<String>, value: impl Into<String>) -> &mut Self {
        self.strings.insert(key.into(), value.into());
        self
    }

    pub fn set_float(&mut self, key: impl Into<String>, value: f64) -> &mut Self {
        self.floats.insert(key.into(), value);
        self
    }

    pub fn set_int(&mut self, key: impl Into<String>, value: i64) -> &mut Self {
        self.ints.insert(key.into(), value);
        self
    }
}

impl KeywordProvider for KeywordMap {
    fn get_string(&self, key: &str) -> Option<String> {
        self.strings.get(key).cloned()
    }

    fn get_float(&self, key: &str) -> Option<f64> {
        self.floats.get(key).copied()
    }

    fn get_int(&self, key: &str) -> Option<i64> {
        self.ints.get(key).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_keyword_map_get_set_all_types() {
        // Chained setters and per-type get/miss in one map.
        let mut map = KeywordMap::new();
        map.set_string("CTYPE1", "RA---TAN")
            .set_float("CRPIX1", 512.0)
            .set_int("NAXIS", 2);

        // Hits.
        assert_eq!(map.get_string("CTYPE1"), Some("RA---TAN".to_string()));
        assert_eq!(map.get_float("CRPIX1"), Some(512.0));
        assert_eq!(map.get_int("NAXIS"), Some(2));

        // Misses for each type return None rather than a default.
        assert_eq!(map.get_string("CTYPE2"), None);
        assert_eq!(map.get_float("CRPIX2"), None);
        assert_eq!(map.get_int("NAXIS1"), None);
    }

    #[test]
    fn test_require_methods_present_and_missing() {
        let mut map = KeywordMap::new();
        map.set_float("CRVAL1", 180.0)
            .set_string("CTYPE1", "RA---TAN");

        // Present keys return the value.
        assert_eq!(map.require_float("CRVAL1").unwrap(), 180.0);
        assert_eq!(map.require_string("CTYPE1").unwrap(), "RA---TAN");

        // Missing keys surface the keyword name in the error.
        let err = map.require_float("CRVAL2").unwrap_err();
        assert!(err.to_string().contains("CRVAL2"));

        let err = map.require_string("CTYPE2").unwrap_err();
        assert!(err.to_string().contains("CTYPE2"));
    }
}
