use crate::xisf::header::property::XisfPropertyValue;

// The value of a scalar or String property: `value` is its value attribute
// and `text` its character content.
pub(super) fn scalar_value(
    type_name: &str,
    value: &str,
    text: Option<String>,
) -> Option<XisfPropertyValue> {
    match type_name {
        "Boolean" => parse_bool(value).map(XisfPropertyValue::Boolean),
        "Int32" => parse_int32(value).map(XisfPropertyValue::Int32),
        "Float64" => parse_f64(value).map(XisfPropertyValue::Float64),
        "String" => text.map(XisfPropertyValue::String),
        _ => None,
    }
}

// XISF 1.0 section 8.3.4.
fn parse_bool(text: &str) -> Option<bool> {
    match text.trim() {
        "true" | "1" => Some(true),
        "false" | "0" => Some(false),
        _ => None,
    }
}

// XISF 1.0 section 8.3.2: binary, octal and hexadecimal forms carry no sign
// and give the 32-bit two's complement pattern.
fn parse_int32(text: &str) -> Option<i32> {
    let text = text.trim();
    let radix = match text.get(..2) {
        Some("0b" | "0B") => 2,
        Some("0o" | "0O") => 8,
        Some("0x" | "0X") => 16,
        _ => return text.parse().ok(),
    };
    let digits = &text[2..];
    if digits.starts_with(['+', '-']) {
        return None;
    }
    u32::from_str_radix(digits, radix)
        .ok()
        .map(|bits| bits as i32)
}

// Rust's parser accepts every form of XISF 1.0 section 8.3.3, including NaN,
// +Inf, -Inf, nan and inf.
fn parse_f64(text: &str) -> Option<f64> {
    text.trim().parse().ok()
}
