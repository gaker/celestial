#[inline]
pub fn horner(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

pub fn chebyshev(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let (mut t_prev, mut t_curr) = (1.0, x);
            for _ in 2..=n {
                let t_next = 2.0 * x * t_curr - t_prev;
                t_prev = t_curr;
                t_curr = t_next;
            }
            t_curr
        }
    }
}

pub fn legendre(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let (mut p_prev, mut p_curr) = (1.0, x);
            for k in 1..n {
                let p_next = ((2 * k + 1) as f64 * x * p_curr - k as f64 * p_prev) / (k + 1) as f64;
                p_prev = p_curr;
                p_curr = p_next;
            }
            p_curr
        }
    }
}

#[inline]
pub fn power_term(x: f64, y: f64, p: u32, q: u32) -> f64 {
    x.powi(p as i32) * y.powi(q as i32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_horner_polynomial_degrees() {
        // Constant, linear, quadratic, cubic - all four exact-match.
        assert_eq!(horner(&[5.0], 10.0), 5.0);
        assert_eq!(horner(&[2.0, 3.0], 5.0), 2.0 + 3.0 * 5.0);
        assert_eq!(horner(&[1.0, 2.0, 3.0], 2.0), 1.0 + 2.0 * 2.0 + 3.0 * 4.0);
        assert_eq!(
            horner(&[1.0, -1.0, 2.0, -2.0], 3.0),
            1.0 - 3.0 + 2.0 * 9.0 - 2.0 * 27.0,
        );
    }

    #[test]
    fn test_chebyshev_polynomials_up_to_t5() {
        // T_0 and T_1 are byte-exact closed forms; T_2..T_5 are derived from
        // the recurrence and compared against the explicit polynomial form.
        assert_eq!(chebyshev(0, 0.5), 1.0);
        assert_eq!(chebyshev(0, -0.3), 1.0);
        assert_eq!(chebyshev(1, 0.5), 0.5);
        assert_eq!(chebyshev(1, -0.7), -0.7);

        let cases: &[(usize, f64, f64)] = &[
            (2, 0.5, 2.0 * 0.5 * 0.5 - 1.0),
            (3, 0.6, 4.0 * 0.6_f64.powi(3) - 3.0 * 0.6),
            (4, 0.4, 8.0 * 0.4_f64.powi(4) - 8.0 * 0.4_f64.powi(2) + 1.0),
            (5, 0.3,
                16.0 * 0.3_f64.powi(5) - 20.0 * 0.3_f64.powi(3) + 5.0 * 0.3),
        ];
        for &(n, x, expected) in cases {
            assert!(
                (chebyshev(n, x) - expected).abs() < 1e-14,
                "T_{}({}) mismatch", n, x,
            );
        }
    }

    #[test]
    fn test_legendre_polynomials_up_to_p5() {
        // P_0 and P_1 are exact; P_2..P_5 come from the three-term recurrence.
        assert_eq!(legendre(0, 0.5), 1.0);
        assert_eq!(legendre(0, -0.3), 1.0);
        assert_eq!(legendre(1, 0.5), 0.5);
        assert_eq!(legendre(1, -0.7), -0.7);

        let cases: &[(usize, f64, f64)] = &[
            (2, 0.5, (3.0 * 0.5 * 0.5 - 1.0) / 2.0),
            (3, 0.6, (5.0 * 0.6_f64.powi(3) - 3.0 * 0.6) / 2.0),
            (4, 0.4,
                (35.0 * 0.4_f64.powi(4) - 30.0 * 0.4_f64.powi(2) + 3.0) / 8.0),
            (5, 0.3,
                (63.0 * 0.3_f64.powi(5) - 70.0 * 0.3_f64.powi(3) + 15.0 * 0.3) / 8.0),
        ];
        for &(n, x, expected) in cases {
            assert!(
                (legendre(n, x) - expected).abs() < 1e-14,
                "P_{}({}) mismatch", n, x,
            );
        }
    }

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
