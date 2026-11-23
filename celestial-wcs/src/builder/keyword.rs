#[derive(Debug, Clone, PartialEq)]
pub enum WcsKeywordValue {
    Real(f64),
    Integer(i64),
    String(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct WcsKeyword {
    pub name: String,
    pub value: WcsKeywordValue,
}

impl WcsKeyword {
    pub fn real(name: impl Into<String>, value: f64) -> Self {
        Self {
            name: name.into(),
            value: WcsKeywordValue::Real(value),
        }
    }

    pub fn integer(name: impl Into<String>, value: i64) -> Self {
        Self {
            name: name.into(),
            value: WcsKeywordValue::Integer(value),
        }
    }

    pub fn string(name: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            value: WcsKeywordValue::String(value.into()),
        }
    }
}

#[cfg(test)]
pub(super) fn keyword_map(keywords: Vec<WcsKeyword>) -> crate::header::KeywordMap {
    let mut header = crate::header::KeywordMap::new();
    for keyword in keywords {
        match keyword.value {
            WcsKeywordValue::Real(v) => header.set_float(keyword.name, v),
            WcsKeywordValue::Integer(v) => header.set_int(keyword.name, v),
            WcsKeywordValue::String(v) => header.set_string(keyword.name, v),
        };
    }
    header
}
