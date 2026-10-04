use crate::matrix::Vector3;

#[inline]
pub fn fmod(x: f64, y: f64) -> f64 {
    libm::fmod(x, y)
}

// Horner's rule, lowest-order coefficient first. Each step is `c + acc * t`, the same
// operations as the nested `c0 + (c1 + (c2 + ...) * t) * t` form written out by hand.
pub(crate) fn polynomial(coeffs: &[f64], t: f64) -> f64 {
    coeffs
        .iter()
        .rev()
        .copied()
        .reduce(|acc, c| c + acc * t)
        .unwrap_or(0.0)
}

pub fn angular_separation(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let a = Vector3::from_spherical(lon1, lat1);
    let b = Vector3::from_spherical(lon2, lat2);
    libm::atan2(a.cross(&b).magnitude(), a.dot(&b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{HALF_PI, PI};

    // (lon1, lat1, lon2, lat2, separation): ERFA seps outputs, run with Rust libm for
    // sin/cos/atan2/sqrt.
    const ERFA_SEPS: [(f64, f64, f64, f64, f64); 7] = [
        (
            1.6752317257415452,
            0.4332597300282792,
            0.8654392245194842,
            0.3845526318877417,
            0.7411237456008262,
        ),
        (
            0.4982450033087944,
            -0.75293570115539,
            -2.9619095109238525,
            0.20420923832323146,
            2.5274010189113785,
        ),
        (1.0, 0.1, 0.2, -1.1, 1.34336844772573),
        (0.0, 0.0, PI, 0.0, PI),
        (0.0, 0.0, 0.0, 0.0, 0.0),
        (0.5, HALF_PI, 2.0, HALF_PI, 8.347667256373471e-17),
        (0.0, 0.0, 1e-9, 0.0, 1e-9),
    ];

    #[test]
    fn test_angular_separation_matches_erfa_seps() {
        for (lon1, lat1, lon2, lat2, expected) in ERFA_SEPS {
            let got = angular_separation(lon1, lat1, lon2, lat2);
            assert_eq!(got, expected, "({lon1}, {lat1}) to ({lon2}, {lat2})");
        }
    }
}
