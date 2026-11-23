pub(crate) fn power_term(x: f64, y: f64, p: u32, q: u32) -> f64 {
    power(x, p) * power(y, q)
}

// powi leaves its rounding sequence unspecified, so it may differ between
// platforms, or between a constant-folded call and a run-time one.
pub(crate) fn power(x: f64, n: u32) -> f64 {
    (0..n).fold(1.0, |acc, _| acc * x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_power_term() {
        // x^p * y^q for representative (p, q) covering both constant and
        // mixed-degree cases.
        let (x, y) = (2.0, 3.0);
        assert_eq!(power_term(x, y, 0, 0), 1.0);
        assert_eq!(power_term(x, y, 1, 0), 2.0);
        assert_eq!(power_term(x, y, 0, 1), 3.0);
        assert_eq!(power_term(x, y, 1, 1), 6.0);
        assert_eq!(power_term(x, y, 2, 0), 4.0);
        assert_eq!(power_term(x, y, 0, 2), 9.0);
        assert_eq!(power_term(x, y, 2, 3), 4.0 * 27.0);
        assert_eq!(power_term(x, y, 3, 2), 8.0 * 9.0);
    }
}
