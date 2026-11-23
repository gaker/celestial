// Writes the shortest digits that parse back to the same f64, with a
// decimal point so readers type the value as real (FITS 4.0 §4.2.4).
// Returns None for NaN and infinities, which FITS cannot represent.
pub(super) fn format_real(value: f64) -> Option<String> {
    if !value.is_finite() {
        return None;
    }
    let shortest = format!("{:e}", value);
    let (mantissa, exponent) = shortest.split_once('e')?;
    let exponent: i32 = exponent.parse().ok()?;
    let (sign, mantissa) = match mantissa.strip_prefix('-') {
        Some(rest) => ("-", rest),
        None => ("", mantissa),
    };
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();

    let fixed = fixed_notation(&digits, exponent);
    let scientific = exponent_notation(&digits, exponent);
    let body = if scientific.len() < fixed.len() {
        scientific
    } else {
        fixed
    };
    Some(format!("{}{}", sign, body))
}

fn fixed_notation(digits: &str, exponent: i32) -> String {
    if exponent < 0 {
        let zeros = "0".repeat(exponent.unsigned_abs() as usize - 1);
        return format!("0.{}{}", zeros, digits);
    }
    let int_len = exponent as usize + 1;
    if digits.len() <= int_len {
        let zeros = "0".repeat(int_len - digits.len());
        format!("{}{}.0", digits, zeros)
    } else {
        format!("{}.{}", &digits[..int_len], &digits[int_len..])
    }
}

fn exponent_notation(digits: &str, exponent: i32) -> String {
    let (lead, rest) = digits.split_at(1);
    let rest = if rest.is_empty() { "0" } else { rest };
    let sign = if exponent < 0 { '-' } else { '+' };
    format!("{}.{}E{}{:02}", lead, rest, sign, exponent.unsigned_abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_and_exponent_notation_place_the_decimal_point() {
        assert_eq!(fixed_notation("512", 2), "512.0");
        assert_eq!(fixed_notation("1", 4), "10000.0");
        assert_eq!(fixed_notation("314", 0), "3.14");
        assert_eq!(fixed_notation("5", -3), "0.005");
        assert_eq!(exponent_notation("1", 25), "1.0E+25");
        assert_eq!(exponent_notation("457", -19), "4.57E-19");
        assert_eq!(exponent_notation("5", -324), "5.0E-324");
    }

    #[test]
    fn format_real_picks_the_shorter_form() {
        assert_eq!(format_real(0.0).as_deref(), Some("0.0"));
        assert_eq!(format_real(-0.0).as_deref(), Some("-0.0"));
        assert_eq!(format_real(1e-5).as_deref(), Some("0.00001"));
        assert_eq!(format_real(1e-6).as_deref(), Some("1.0E-06"));
        assert_eq!(format_real(1e16).as_deref(), Some("1.0E+16"));
        assert_eq!(format_real(123456.0).as_deref(), Some("123456.0"));
    }
}
