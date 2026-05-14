#[derive(Debug, Clone, PartialEq)]
pub enum XisfPropertyValue {
    String(String),
    Boolean(bool),
    Float64(f64),
    Int32(i32),
    F64Vector(Vec<f64>),
    F64Matrix {
        rows: usize,
        cols: usize,
        data: Vec<f64>,
    },
}

impl XisfPropertyValue {
    pub fn type_name(&self) -> &'static str {
        match self {
            Self::String(_) => "String",
            Self::Boolean(_) => "Boolean",
            Self::Float64(_) => "Float64",
            Self::Int32(_) => "Int32",
            Self::F64Vector(_) => "F64Vector",
            Self::F64Matrix { .. } => "F64Matrix",
        }
    }

    pub fn format_value(&self) -> String {
        match self {
            Self::String(s) => s.clone(),
            Self::Boolean(b) => if *b { "1" } else { "0" }.to_string(),
            Self::Float64(v) => format!("{:?}", v),
            Self::Int32(v) => format!("{}", v),
            Self::F64Vector(v) => v
                .iter()
                .map(|x| format!("{:?}", x))
                .collect::<Vec<_>>()
                .join(" "),
            Self::F64Matrix { data, .. } => data
                .iter()
                .map(|x| format!("{:?}", x))
                .collect::<Vec<_>>()
                .join(" "),
        }
    }

    pub fn is_scalar(&self) -> bool {
        matches!(self, Self::Float64(_) | Self::Int32(_) | Self::Boolean(_))
    }

    pub fn needs_data_block(&self) -> bool {
        matches!(self, Self::F64Vector(_) | Self::F64Matrix { .. })
    }

    pub fn to_le_bytes(&self) -> Option<Vec<u8>> {
        match self {
            Self::F64Vector(v) => {
                let mut bytes = Vec::with_capacity(v.len() * 8);
                for &val in v {
                    bytes.extend_from_slice(&val.to_le_bytes());
                }
                Some(bytes)
            }
            Self::F64Matrix { data, .. } => {
                let mut bytes = Vec::with_capacity(data.len() * 8);
                for &val in data {
                    bytes.extend_from_slice(&val.to_le_bytes());
                }
                Some(bytes)
            }
            _ => None,
        }
    }

    pub fn data_size(&self) -> usize {
        match self {
            Self::F64Vector(v) => v.len() * 8,
            Self::F64Matrix { data, .. } => data.len() * 8,
            _ => 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct XisfProperty {
    pub id: String,
    pub value: XisfPropertyValue,
}

impl XisfProperty {
    pub fn new(id: impl Into<String>, value: XisfPropertyValue) -> Self {
        Self {
            id: id.into(),
            value,
        }
    }

    pub fn string(id: impl Into<String>, value: impl Into<String>) -> Self {
        Self::new(id, XisfPropertyValue::String(value.into()))
    }

    pub fn boolean(id: impl Into<String>, value: bool) -> Self {
        Self::new(id, XisfPropertyValue::Boolean(value))
    }

    pub fn float64(id: impl Into<String>, value: f64) -> Self {
        Self::new(id, XisfPropertyValue::Float64(value))
    }

    pub fn int32(id: impl Into<String>, value: i32) -> Self {
        Self::new(id, XisfPropertyValue::Int32(value))
    }

    pub fn f64_vector(id: impl Into<String>, values: Vec<f64>) -> Self {
        Self::new(id, XisfPropertyValue::F64Vector(values))
    }

    pub fn f64_matrix(id: impl Into<String>, rows: usize, cols: usize, data: Vec<f64>) -> Self {
        Self::new(id, XisfPropertyValue::F64Matrix { rows, cols, data })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::PI;

    #[test]
    fn type_name_covers_all_variants() {
        assert_eq!(XisfPropertyValue::String("x".into()).type_name(), "String");
        assert_eq!(XisfPropertyValue::Boolean(true).type_name(), "Boolean");
        assert_eq!(XisfPropertyValue::Float64(0.0).type_name(), "Float64");
        assert_eq!(XisfPropertyValue::Int32(0).type_name(), "Int32");
        assert_eq!(XisfPropertyValue::F64Vector(vec![]).type_name(), "F64Vector");
        assert_eq!(
            XisfPropertyValue::F64Matrix {
                rows: 0,
                cols: 0,
                data: vec![]
            }
            .type_name(),
            "F64Matrix"
        );
    }

    #[test]
    fn format_value_string_passthrough() {
        assert_eq!(XisfPropertyValue::String("hello".into()).format_value(), "hello");
        assert_eq!(XisfPropertyValue::String(String::new()).format_value(), "");
    }

    #[test]
    fn format_value_boolean_is_zero_or_one() {
        assert_eq!(XisfPropertyValue::Boolean(true).format_value(), "1");
        assert_eq!(XisfPropertyValue::Boolean(false).format_value(), "0");
    }

    #[test]
    fn format_value_int32_uses_display() {
        assert_eq!(XisfPropertyValue::Int32(42).format_value(), "42");
        assert_eq!(XisfPropertyValue::Int32(-1).format_value(), "-1");
        assert_eq!(XisfPropertyValue::Int32(i32::MAX).format_value(), "2147483647");
        assert_eq!(XisfPropertyValue::Int32(i32::MIN).format_value(), "-2147483648");
    }

    #[test]
    fn format_value_float_preserves_integer_decimals() {
        assert_eq!(XisfPropertyValue::Float64(1.0).format_value(), "1.0");
        assert_eq!(XisfPropertyValue::Float64(0.0).format_value(), "0.0");
        assert_eq!(XisfPropertyValue::Float64(-3.0).format_value(), "-3.0");
    }

    #[test]
    fn format_value_float_roundtrips_through_parse() {
        for value in [1.234567890123456e-9, -9.87654321e10, PI] {
            let formatted = XisfPropertyValue::Float64(value).format_value();
            let parsed: f64 = formatted.parse().expect("formatted value parses");
            assert_eq!(parsed, value);
        }
    }

    #[test]
    fn format_value_float_handles_special_values() {
        assert_eq!(XisfPropertyValue::Float64(f64::INFINITY).format_value(), "inf");
        assert_eq!(XisfPropertyValue::Float64(f64::NEG_INFINITY).format_value(), "-inf");
        assert_eq!(XisfPropertyValue::Float64(f64::NAN).format_value(), "NaN");
    }

    #[test]
    fn format_value_vector_space_separated() {
        assert_eq!(
            XisfPropertyValue::F64Vector(vec![1.5, 2.5, 3.5]).format_value(),
            "1.5 2.5 3.5"
        );
        assert_eq!(XisfPropertyValue::F64Vector(vec![]).format_value(), "");
        assert_eq!(XisfPropertyValue::F64Vector(vec![42.0]).format_value(), "42.0");
    }

    #[test]
    fn format_value_matrix_ignores_shape_in_text() {
        let m = XisfPropertyValue::F64Matrix {
            rows: 2,
            cols: 2,
            data: vec![1.0, 0.0, 0.0, 1.0],
        };
        assert_eq!(m.format_value(), "1.0 0.0 0.0 1.0");
    }

    #[test]
    fn is_scalar_true_only_for_primitive_scalars() {
        assert!(XisfPropertyValue::Float64(1.0).is_scalar());
        assert!(XisfPropertyValue::Int32(1).is_scalar());
        assert!(XisfPropertyValue::Boolean(true).is_scalar());
        assert!(!XisfPropertyValue::String("x".into()).is_scalar());
        assert!(!XisfPropertyValue::F64Vector(vec![]).is_scalar());
        assert!(!XisfPropertyValue::F64Matrix {
            rows: 0,
            cols: 0,
            data: vec![]
        }
        .is_scalar());
    }

    #[test]
    fn needs_data_block_true_only_for_bulk_types() {
        assert!(XisfPropertyValue::F64Vector(vec![]).needs_data_block());
        assert!(XisfPropertyValue::F64Matrix {
            rows: 0,
            cols: 0,
            data: vec![]
        }
        .needs_data_block());
        assert!(!XisfPropertyValue::String("x".into()).needs_data_block());
        assert!(!XisfPropertyValue::Float64(0.0).needs_data_block());
        assert!(!XisfPropertyValue::Int32(0).needs_data_block());
        assert!(!XisfPropertyValue::Boolean(false).needs_data_block());
    }

    #[test]
    fn is_scalar_and_needs_data_block_are_disjoint() {
        let cases = [
            XisfPropertyValue::String("x".into()),
            XisfPropertyValue::Boolean(true),
            XisfPropertyValue::Float64(1.0),
            XisfPropertyValue::Int32(1),
            XisfPropertyValue::F64Vector(vec![1.0]),
            XisfPropertyValue::F64Matrix {
                rows: 1,
                cols: 1,
                data: vec![1.0],
            },
        ];
        for v in cases {
            assert!(!(v.is_scalar() && v.needs_data_block()));
        }
    }

    #[test]
    fn to_le_bytes_vector_encodes_in_little_endian() {
        let v = XisfPropertyValue::F64Vector(vec![1.0, 2.0]);
        let bytes = v.to_le_bytes().expect("vector produces bytes");
        assert_eq!(bytes.len(), 16);
        assert_eq!(&bytes[0..8], &1.0f64.to_le_bytes());
        assert_eq!(&bytes[8..16], &2.0f64.to_le_bytes());
    }

    #[test]
    fn to_le_bytes_matrix_encodes_row_major() {
        let m = XisfPropertyValue::F64Matrix {
            rows: 2,
            cols: 2,
            data: vec![1.0, 2.0, 3.0, 4.0],
        };
        let bytes = m.to_le_bytes().expect("matrix produces bytes");
        assert_eq!(bytes.len(), 32);
        for (i, expected) in [1.0f64, 2.0, 3.0, 4.0].iter().enumerate() {
            let slice = &bytes[i * 8..(i + 1) * 8];
            assert_eq!(slice, &expected.to_le_bytes());
        }
    }

    #[test]
    fn to_le_bytes_none_for_scalar_types() {
        assert!(XisfPropertyValue::String("x".into()).to_le_bytes().is_none());
        assert!(XisfPropertyValue::Boolean(true).to_le_bytes().is_none());
        assert!(XisfPropertyValue::Float64(1.0).to_le_bytes().is_none());
        assert!(XisfPropertyValue::Int32(1).to_le_bytes().is_none());
    }

    #[test]
    fn to_le_bytes_empty_vector_returns_empty_buffer() {
        let bytes = XisfPropertyValue::F64Vector(vec![]).to_le_bytes().unwrap();
        assert!(bytes.is_empty());
    }

    #[test]
    fn data_size_matches_to_le_bytes_length() {
        let cases: [XisfPropertyValue; 6] = [
            XisfPropertyValue::String("x".into()),
            XisfPropertyValue::Boolean(true),
            XisfPropertyValue::Float64(1.0),
            XisfPropertyValue::Int32(1),
            XisfPropertyValue::F64Vector(vec![1.0, 2.0, 3.0]),
            XisfPropertyValue::F64Matrix {
                rows: 2,
                cols: 3,
                data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            },
        ];
        for v in cases {
            let expected = v.to_le_bytes().map(|b| b.len()).unwrap_or(0);
            assert_eq!(v.data_size(), expected);
        }
    }

    #[test]
    fn property_new_stores_id_and_value() {
        let prop = XisfProperty::new("Test:Id", XisfPropertyValue::Int32(7));
        assert_eq!(prop.id, "Test:Id");
        assert_eq!(prop.value, XisfPropertyValue::Int32(7));
    }

    #[test]
    fn property_string_constructor() {
        let prop = XisfProperty::string("Id", "hello");
        assert_eq!(prop.id, "Id");
        assert_eq!(prop.value, XisfPropertyValue::String("hello".into()));
    }

    #[test]
    fn property_boolean_constructor() {
        let prop = XisfProperty::boolean("Id", true);
        assert_eq!(prop.value, XisfPropertyValue::Boolean(true));
    }

    #[test]
    fn property_float64_constructor() {
        let prop = XisfProperty::float64("Id", PI);
        assert_eq!(prop.value, XisfPropertyValue::Float64(PI));
    }

    #[test]
    fn property_int32_constructor() {
        let prop = XisfProperty::int32("Id", -42);
        assert_eq!(prop.value, XisfPropertyValue::Int32(-42));
    }

    #[test]
    fn property_f64_vector_constructor() {
        let prop = XisfProperty::f64_vector("Id", vec![1.0, 2.0]);
        assert_eq!(prop.value, XisfPropertyValue::F64Vector(vec![1.0, 2.0]));
    }

    #[test]
    fn property_f64_matrix_constructor_preserves_shape() {
        let prop = XisfProperty::f64_matrix("Id", 2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        match prop.value {
            XisfPropertyValue::F64Matrix { rows, cols, data } => {
                assert_eq!(rows, 2);
                assert_eq!(cols, 3);
                assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
            }
            _ => panic!("expected F64Matrix"),
        }
    }

    #[test]
    fn property_constructors_accept_string_and_str_ids() {
        let from_str = XisfProperty::string("Id", "v");
        let from_string = XisfProperty::string(String::from("Id"), "v");
        assert_eq!(from_str.id, from_string.id);
    }
}
