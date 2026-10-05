pub(super) fn value(coeffs: &[f64], tau: f64) -> f64 {
    let Some((&first, rest)) = coeffs.split_first() else {
        return 0.0;
    };
    let two_tau = 2.0 * tau;
    let (mut b_k, mut b_k1) = (0.0, 0.0);
    for &coeff in rest.iter().rev() {
        (b_k, b_k1) = (two_tau * b_k - b_k1 + coeff, b_k);
    }
    tau * b_k - b_k1 + first
}

pub(super) fn derivative(coeffs: &[f64], tau: f64, radius: f64) -> f64 {
    let [_, first, rest @ ..] = coeffs else {
        return 0.0;
    };
    let two_tau = 2.0 * tau;
    let (mut u_prev, mut u) = (1.0, two_tau);
    let mut sum = *first;
    for (i, &coeff) in rest.iter().enumerate() {
        sum += ((i + 2) as f64) * coeff * u;
        (u_prev, u) = (u, two_tau * u - u_prev);
    }
    sum / radius
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_series() {
        assert_eq!(value(&[5.0, 0.0, 0.0], 0.0), 5.0);
        assert_eq!(value(&[5.0, 0.0, 0.0], 0.5), 5.0);
        assert_eq!(derivative(&[5.0, 0.0, 0.0], 0.5, 2.0), 0.0);
    }

    #[test]
    fn linear_series() {
        assert_eq!(value(&[0.0, 1.0, 0.0], 0.5), 0.5);
        assert_eq!(value(&[0.0, 1.0, 0.0], -0.5), -0.5);
        assert_eq!(derivative(&[0.0, 3.0, 0.0], 0.0, 2.0), 1.5);
        assert_eq!(derivative(&[0.0, 3.0, 0.0], 0.0, 1.0), 3.0);
    }

    #[test]
    fn quadratic_series() {
        // 1·T0 + 1·T2 = 2τ², whose derivative is 4τ.
        assert_eq!(value(&[1.0, 0.0, 1.0], 0.5), 0.5);
        assert_eq!(value(&[1.0, 0.0, 1.0], 1.0), 2.0);
        assert_eq!(derivative(&[1.0, 0.0, 1.0], 0.5, 1.0), 2.0);
    }

    #[test]
    fn cubic_series() {
        // T3 = 4τ³ − 3τ, whose derivative is 12τ² − 3.
        assert_eq!(value(&[0.0, 0.0, 0.0, 1.0], 0.5), -1.0);
        assert_eq!(derivative(&[0.0, 0.0, 0.0, 1.0], 0.5, 1.0), 0.0);
        assert_eq!(derivative(&[0.0, 0.0, 0.0, 1.0], 1.0, 3.0), 3.0);
    }

    #[test]
    fn short_series() {
        assert_eq!(value(&[], 0.5), 0.0);
        assert_eq!(value(&[7.0], 0.5), 7.0);
        assert_eq!(derivative(&[7.0], 0.5, 1.0), 0.0);
        assert_eq!(derivative(&[], 0.5, 1.0), 0.0);
    }
}
