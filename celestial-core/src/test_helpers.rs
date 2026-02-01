#[inline]
pub fn f64_to_ordered_u64(x: f64) -> u64 {
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
    assert!(
        a.is_finite() && b.is_finite(),
        "non-finite value in {}",
        ctx
    );
    let d = ulp_diff(a, b);
    assert!(
        d <= max_ulp,
        "{}: ULP={} exceeds {}, a={} (0x{:016x}) b={} (0x{:016x})",
        ctx,
        d,
        max_ulp,
        a,
        a.to_bits(),
        b,
        b.to_bits()
    );
}

#[track_caller]
pub fn assert_float_eq(a: f64, b: f64, max_ulp: u64) {
    if a == 0.0 && b == 0.0 {
        return;
    }
    assert!(a.is_finite() && b.is_finite());
    let d = ulp_diff(a, b);
    assert!(
        d <= max_ulp,
        "ULP={} exceeds {}, a={} (0x{:016x}) b={} (0x{:016x})",
        d,
        max_ulp,
        a,
        a.to_bits(),
        b,
        b.to_bits()
    );
}

#[macro_export]
macro_rules! assert_ulp_lt {
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
