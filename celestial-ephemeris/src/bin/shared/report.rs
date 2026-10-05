// An empty series keeps 0% rather than NaN.
pub(crate) fn percent(part: usize, whole: usize) -> f64 {
    if whole > 0 {
        (part as f64 / whole as f64) * 100.0
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percentages() {
        assert_eq!(percent(1, 4), 25.0);
        assert_eq!(percent(3, 3), 100.0);
        assert_eq!(percent(0, 0), 0.0);
    }
}
