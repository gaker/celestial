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
