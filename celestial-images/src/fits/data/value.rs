#[derive(Debug, Clone, PartialEq)]
pub enum TableValue {
    Logical(bool),
    Byte(u8),
    I8(i8),
    I16(i16),
    U16(u16),
    I32(i32),
    U32(u32),
    I64(i64),
    U64(u64),
    F32(f32),
    F64(f64),
    String(String),
    Complex32(f32, f32),
    Complex64(f64, f64),
    Bits(Vec<bool>),
    Array(Vec<TableValue>),
    Null,
}

impl TableValue {
    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Self::Byte(v) => Some(i64::from(*v)),
            Self::I8(v) => Some(i64::from(*v)),
            Self::I16(v) => Some(i64::from(*v)),
            Self::U16(v) => Some(i64::from(*v)),
            Self::I32(v) => Some(i64::from(*v)),
            Self::U32(v) => Some(i64::from(*v)),
            Self::I64(v) => Some(*v),
            Self::U64(v) => i64::try_from(*v).ok(),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Self::Byte(v) => Some(f64::from(*v)),
            Self::I8(v) => Some(f64::from(*v)),
            Self::I16(v) => Some(f64::from(*v)),
            Self::U16(v) => Some(f64::from(*v)),
            Self::I32(v) => Some(f64::from(*v)),
            Self::U32(v) => Some(f64::from(*v)),
            Self::I64(v) => Some(*v as f64),
            Self::U64(v) => Some(*v as f64),
            Self::F32(v) => Some(f64::from(*v)),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_value_as_i64_integer_types() {
        assert_eq!(TableValue::Byte(42).as_i64(), Some(42));
        assert_eq!(TableValue::I16(1000).as_i64(), Some(1000));
        assert_eq!(TableValue::I32(-500).as_i64(), Some(-500));
        assert_eq!(TableValue::I64(123456789).as_i64(), Some(123456789));
        assert_eq!(TableValue::I8(-5).as_i64(), Some(-5));
        assert_eq!(TableValue::U16(65535).as_i64(), Some(65535));
        assert_eq!(TableValue::U32(u32::MAX).as_i64(), Some(4294967295));
        assert_eq!(TableValue::U64(1 << 63).as_i64(), None);
        assert_eq!(TableValue::U64(7).as_f64(), Some(7.0));
        assert_eq!(TableValue::Array(vec![TableValue::I16(1)]).as_i64(), None);
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
