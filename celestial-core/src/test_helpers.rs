#[inline]
fn f64_to_ordered_u64(x: f64) -> u64 {
    let bits = x.to_bits();
    if bits & 0x8000_0000_0000_0000 != 0 {
        !bits
    } else {
        bits | 0x8000_0000_0000_0000
    }
}

#[inline]
pub fn ulp_diff(a: f64, b: f64) -> u64 {
    let ua = f64_to_ordered_u64(a);
    let ub = f64_to_ordered_u64(b);
    ua.abs_diff(ub)
}

#[track_caller]
pub fn assert_ulp_le(a: f64, b: f64, max_ulp: u64, ctx: &str) {
    if a == 0.0 && b == 0.0 {
        return;
    }
    assert!(a.is_finite() && b.is_finite(), "non-finite value in {ctx}");
    let d = ulp_diff(a, b);
    let (a_bits, b_bits) = (a.to_bits(), b.to_bits());
    assert!(
        d <= max_ulp,
        "{ctx}: ULP={d} exceeds {max_ulp}, a={a} (0x{a_bits:016x}) b={b} (0x{b_bits:016x})"
    );
}

#[macro_export]
macro_rules! assert_ulp_le {
    ($a:expr, $b:expr, $max_ulp:expr) => {
        $crate::test_helpers::assert_ulp_le(
            $a,
            $b,
            $max_ulp,
            &format!(
                "ULP check failed: {} vs {} (max_ulp={})",
                stringify!($a),
                stringify!($b),
                $max_ulp
            ),
        )
    };
    ($a:expr, $b:expr, $max_ulp:expr, $($arg:tt)*) => {
        $crate::test_helpers::assert_ulp_le($a, $b, $max_ulp, &format!($($arg)*))
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    fn next_up(x: f64) -> f64 {
        f64::from_bits(x.to_bits() + 1)
    }

    #[test]
    fn test_ulp_diff_counts_adjacent_doubles() {
        assert_eq!(ulp_diff(1.0, 1.0), 0);
        assert_eq!(ulp_diff(1.0, next_up(1.0)), 1);
        assert_eq!(ulp_diff(next_up(next_up(1.0)), 1.0), 2);
    }

    #[test]
    fn test_ulp_diff_across_zero() {
        // +0 and -0 are distinct steps in the ordering, so crossing zero counts one more
        // than the number of doubles in between. Over-counting keeps the checks strict.
        assert_eq!(ulp_diff(0.0, -0.0), 1);
        assert_eq!(ulp_diff(5e-324, -5e-324), 3);
    }

    #[test]
    fn test_assert_ulp_le_passes_at_limit() {
        assert_ulp_le(1.0, next_up(next_up(1.0)), 2, "at limit");
        assert_ulp_le(0.0, -0.0, 0, "signed zeros");
    }

    #[test]
    #[should_panic(expected = "over: ULP=2 exceeds 1")]
    fn test_assert_ulp_le_panics_over_limit() {
        assert_ulp_le(1.0, next_up(next_up(1.0)), 1, "over");
    }

    #[test]
    #[should_panic(expected = "non-finite value in nan")]
    fn test_assert_ulp_le_rejects_nan() {
        assert_ulp_le(f64::NAN, 1.0, u64::MAX, "nan");
    }

    #[test]
    fn test_macro_forms() {
        crate::assert_ulp_le!(1.0, next_up(1.0), 1);
        crate::assert_ulp_le!(1.0, next_up(1.0), 1, "case {}", 7);
    }

    #[test]
    #[should_panic(expected = "case 7")]
    fn test_macro_with_context_panics_with_context() {
        crate::assert_ulp_le!(1.0, next_up(next_up(1.0)), 1, "case {}", 7);
    }
}
