#[derive(Debug, Clone, PartialEq)]
pub enum DataValue<T> {
    Value(T),
    Null,
}

#[derive(Debug, Clone, PartialEq)]
pub enum TableValue {
    Logical(bool),
    Byte(u8),
    I16(i16),
    I32(i32),
    I64(i64),
    F32(f32),
    F64(f64),
    String(String),
    Complex32(f32, f32),
    Complex64(f64, f64),
    Null,
}

impl TableValue {
    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Self::Byte(v) => Some(*v as i64),
            Self::I16(v) => Some(*v as i64),
            Self::I32(v) => Some(*v as i64),
            Self::I64(v) => Some(*v),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Self::Byte(v) => Some(*v as f64),
            Self::I16(v) => Some(*v as f64),
            Self::I32(v) => Some(*v as f64),
            Self::I64(v) => Some(*v as f64),
            Self::F32(v) => Some(*v as f64),
            Self::F64(v) => Some(*v),
            _ => None,
        }
    }

    pub fn as_string(&self) -> Option<&str> {
        match self {
            Self::String(s) => Some(s),
            _ => None,
        }
    }

    pub fn is_null(&self) -> bool {
        matches!(self, Self::Null)
    }
}

impl<T> DataValue<T> {
    pub fn is_null(&self) -> bool {
        matches!(self, Self::Null)
    }

    pub fn value(&self) -> Option<&T> {
        match self {
            Self::Value(v) => Some(v),
            Self::Null => None,
        }
    }

    pub fn unwrap_or(self, default: T) -> T {
        match self {
            Self::Value(v) => v,
            Self::Null => default,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn data_value_creation() {
        let val = DataValue::Value(42i16);
        let null = DataValue::<i16>::Null;

        assert!(!val.is_null());
        assert!(null.is_null());
        assert_eq!(val.value(), Some(&42));
        assert_eq!(null.value(), None);
    }

    #[test]
    fn data_value_unwrap_or() {
        let val = DataValue::Value(42i16);
        let null = DataValue::<i16>::Null;

        assert_eq!(val.unwrap_or(0), 42);
        assert_eq!(null.unwrap_or(0), 0);
    }

    #[test]
    fn data_value_operations() {
        let value = DataValue::Value(42);
        let null = DataValue::<i32>::Null;

        assert!(!value.is_null());
        assert!(null.is_null());

        assert_eq!(value.value(), Some(&42));
        assert_eq!(null.value(), None);

        let cloned_value = value.clone();
        assert_eq!(cloned_value, DataValue::Value(42));

        let cloned_null = null.clone();
        assert_eq!(cloned_null, DataValue::<i32>::Null);

        assert_eq!(value.unwrap_or(-1), 42);
        assert_eq!(null.clone().unwrap_or(-1), -1);
        assert_eq!(null.unwrap_or(999), 999);
    }

    #[test]
    fn table_value_as_i64_integer_types() {
        assert_eq!(TableValue::Byte(42).as_i64(), Some(42));
        assert_eq!(TableValue::I16(1000).as_i64(), Some(1000));
        assert_eq!(TableValue::I32(-500).as_i64(), Some(-500));
        assert_eq!(TableValue::I64(123456789).as_i64(), Some(123456789));
    }

    #[test]
    fn table_value_as_i64_non_integer_types() {
        assert_eq!(TableValue::F32(1.5).as_i64(), None);
        assert_eq!(TableValue::F64(2.75).as_i64(), None);
        assert_eq!(TableValue::String("test".to_string()).as_i64(), None);
        assert_eq!(TableValue::Null.as_i64(), None);
        assert_eq!(TableValue::Logical(true).as_i64(), None);
    }

    #[test]
    fn table_value_as_f64_numeric_types() {
        assert_eq!(TableValue::Byte(42).as_f64(), Some(42.0));
        assert_eq!(TableValue::I16(1000).as_f64(), Some(1000.0));
        assert_eq!(TableValue::I32(-500).as_f64(), Some(-500.0));
        assert_eq!(TableValue::I64(123456789).as_f64(), Some(123456789.0));
        assert_eq!(TableValue::F32(1.5).as_f64(), Some(1.5));
        assert_eq!(TableValue::F64(2.75).as_f64(), Some(2.75));
    }

    #[test]
    fn table_value_as_f64_non_numeric_types() {
        assert_eq!(TableValue::String("test".to_string()).as_f64(), None);
        assert_eq!(TableValue::Null.as_f64(), None);
        assert_eq!(TableValue::Logical(true).as_f64(), None);
        assert_eq!(TableValue::Complex32(1.0, 2.0).as_f64(), None);
    }

    #[test]
    fn table_value_as_string() {
        assert_eq!(
            TableValue::String("hello".to_string()).as_string(),
            Some("hello")
        );
        assert_eq!(TableValue::I32(42).as_string(), None);
        assert_eq!(TableValue::Null.as_string(), None);
    }

    #[test]
    fn table_value_is_null() {
        assert!(TableValue::Null.is_null());
        assert!(!TableValue::I32(0).is_null());
        assert!(!TableValue::String("".to_string()).is_null());
    }

    #[test]
    fn table_value_equality() {
        assert_eq!(TableValue::I32(42), TableValue::I32(42));
        assert_ne!(TableValue::I32(42), TableValue::I32(0));
        assert_ne!(TableValue::I32(42), TableValue::I64(42));
        assert_eq!(TableValue::Null, TableValue::Null);
    }

    #[test]
    fn table_value_clone() {
        let original = TableValue::String("test".to_string());
        let cloned = original.clone();
        assert_eq!(original, cloned);
    }

    #[test]
    fn table_value_complex_types() {
        let c32 = TableValue::Complex32(1.0, 2.0);
        let c64 = TableValue::Complex64(3.0, 4.0);

        assert_eq!(c32, TableValue::Complex32(1.0, 2.0));
        assert_eq!(c64, TableValue::Complex64(3.0, 4.0));
        assert_ne!(c32, c64);
    }
}
