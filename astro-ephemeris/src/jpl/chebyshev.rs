use super::SpkError;

pub fn evaluate_chebyshev(
    coeffs: &[f64],
    n_coeffs: usize,
    t_normalized: f64,
) -> Result<f64, SpkError> {
    if n_coeffs == 0 {
        return Err(SpkError::InvalidData("Empty coefficient array".into()));
    }
    if coeffs.len() < n_coeffs {
        return Err(SpkError::InvalidData("Insufficient coefficients".into()));
    }
    let mut b_k1 = 0.0;
    let mut b_k = 0.0;
    let two_t = 2.0 * t_normalized;
    for i in (1..n_coeffs).rev() {
        let b_k_prev = b_k;
        b_k = two_t * b_k - b_k1 + coeffs[i];
        b_k1 = b_k_prev;
    }
    Ok(t_normalized * b_k - b_k1 + coeffs[0])
}

pub fn evaluate_chebyshev_derivative(
    coeffs: &[f64],
    n_coeffs: usize,
    t_normalized: f64,
    half_interval: f64,
) -> Result<f64, SpkError> {
    if n_coeffs < 2 {
        return Ok(0.0);
    }
    if coeffs.len() < n_coeffs {
        return Err(SpkError::InvalidData("Insufficient coefficients".into()));
    }
    let two_t = 2.0 * t_normalized;
    let mut u_prev = 1.0;
    let mut u_curr = two_t;
    let mut derivative = coeffs[1];
    for (i, &coeff) in coeffs.iter().enumerate().take(n_coeffs).skip(2) {
        derivative += (i as f64) * coeff * u_curr;
        let u_next = two_t * u_curr - u_prev;
        u_prev = u_curr;
        u_curr = u_next;
    }
    Ok(derivative / half_interval)
}

pub fn evaluate_position_velocity(
    coeffs_x: &[f64],
    coeffs_y: &[f64],
    coeffs_z: &[f64],
    n_coeffs: usize,
    t_normalized: f64,
    half_interval: f64,
) -> Result<([f64; 3], [f64; 3]), SpkError> {
    let px = evaluate_chebyshev(coeffs_x, n_coeffs, t_normalized)?;
    let py = evaluate_chebyshev(coeffs_y, n_coeffs, t_normalized)?;
    let pz = evaluate_chebyshev(coeffs_z, n_coeffs, t_normalized)?;
    let vx = evaluate_chebyshev_derivative(coeffs_x, n_coeffs, t_normalized, half_interval)?;
    let vy = evaluate_chebyshev_derivative(coeffs_y, n_coeffs, t_normalized, half_interval)?;
    let vz = evaluate_chebyshev_derivative(coeffs_z, n_coeffs, t_normalized, half_interval)?;
    Ok(([px, py, pz], [vx, vy, vz]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chebyshev_constant() {
        let coeffs = [5.0, 0.0, 0.0];
        let result = evaluate_chebyshev(&coeffs, 3, 0.0).unwrap();
        assert!((result - 5.0).abs() < 1e-14);
        let result = evaluate_chebyshev(&coeffs, 3, 0.5).unwrap();
        assert!((result - 5.0).abs() < 1e-14);
    }

    #[test]
    fn test_chebyshev_linear() {
        let coeffs = [0.0, 1.0, 0.0];
        let result = evaluate_chebyshev(&coeffs, 3, 0.5).unwrap();
        assert!((result - 0.5).abs() < 1e-14);
        let result = evaluate_chebyshev(&coeffs, 3, -0.5).unwrap();
        assert!((result - (-0.5)).abs() < 1e-14);
    }

    #[test]
    fn test_chebyshev_derivative_linear() {
        let coeffs = [0.0, 3.0, 0.0];
        let half_interval = 2.0;
        let result = evaluate_chebyshev_derivative(&coeffs, 3, 0.0, half_interval).unwrap();
        assert!((result - 1.5).abs() < 1e-14);
        let half_interval = 1.0;
        let result = evaluate_chebyshev_derivative(&coeffs, 3, 0.0, half_interval).unwrap();
        assert!((result - 3.0).abs() < 1e-14);
    }

    #[test]
    fn test_chebyshev_empty_coeffs() {
        let coeffs: [f64; 0] = [];
        let result = evaluate_chebyshev(&coeffs, 0, 0.0);
        assert!(result.is_err());
        match result.unwrap_err() {
            SpkError::InvalidData(msg) => assert!(msg.contains("Empty")),
            _ => panic!("Expected InvalidData error"),
        }
    }

    #[test]
    fn test_chebyshev_insufficient_coeffs() {
        let coeffs = [1.0, 2.0];
        let result = evaluate_chebyshev(&coeffs, 5, 0.0);
        assert!(result.is_err());
        match result.unwrap_err() {
            SpkError::InvalidData(msg) => assert!(msg.contains("Insufficient")),
            _ => panic!("Expected InvalidData error"),
        }
    }

    #[test]
    fn test_chebyshev_derivative_less_than_two_coeffs() {
        let coeffs = [1.0];
        let result = evaluate_chebyshev_derivative(&coeffs, 1, 0.0, 1.0).unwrap();
        assert!((result - 0.0).abs() < 1e-14);
    }

    #[test]
    fn test_chebyshev_derivative_insufficient_coeffs() {
        let coeffs = [1.0, 2.0];
        let result = evaluate_chebyshev_derivative(&coeffs, 5, 0.0, 1.0);
        assert!(result.is_err());
        match result.unwrap_err() {
            SpkError::InvalidData(msg) => assert!(msg.contains("Insufficient")),
            _ => panic!("Expected InvalidData error"),
        }
    }

    #[test]
    fn test_evaluate_position_velocity_basic() {
        let coeffs_x = [1.0, 0.5, 0.0];
        let coeffs_y = [2.0, 0.3, 0.0];
        let coeffs_z = [3.0, 0.1, 0.0];
        let n_coeffs = 3;
        let t_normalized = 0.0;
        let half_interval = 1.0;

        let (pos, vel) = evaluate_position_velocity(
            &coeffs_x,
            &coeffs_y,
            &coeffs_z,
            n_coeffs,
            t_normalized,
            half_interval,
        )
        .unwrap();

        // At t=0, position should be just the first coefficient
        assert!((pos[0] - 1.0).abs() < 1e-14);
        assert!((pos[1] - 2.0).abs() < 1e-14);
        assert!((pos[2] - 3.0).abs() < 1e-14);

        // Velocity is derivative of position
        assert!((vel[0] - 0.5).abs() < 1e-14);
        assert!((vel[1] - 0.3).abs() < 1e-14);
        assert!((vel[2] - 0.1).abs() < 1e-14);
    }

    #[test]
    fn test_evaluate_position_velocity_at_nonzero_t() {
        let coeffs_x = [0.0, 1.0, 0.0];
        let coeffs_y = [0.0, 2.0, 0.0];
        let coeffs_z = [0.0, 3.0, 0.0];
        let n_coeffs = 3;
        let t_normalized = 0.5;
        let half_interval = 1.0;

        let (pos, _vel) = evaluate_position_velocity(
            &coeffs_x,
            &coeffs_y,
            &coeffs_z,
            n_coeffs,
            t_normalized,
            half_interval,
        )
        .unwrap();

        // Linear: pos = t
        assert!((pos[0] - 0.5).abs() < 1e-14);
        assert!((pos[1] - 1.0).abs() < 1e-14);
        assert!((pos[2] - 1.5).abs() < 1e-14);
    }

    #[test]
    fn test_evaluate_position_velocity_error_propagation() {
        let coeffs_x = [1.0];
        let coeffs_y = [2.0, 0.3, 0.0];
        let coeffs_z = [3.0, 0.1, 0.0];

        // coeffs_x is too short for 3 coefficients
        let result = evaluate_position_velocity(&coeffs_x, &coeffs_y, &coeffs_z, 3, 0.0, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_chebyshev_quadratic() {
        // T_2(x) = 2x^2 - 1
        // So if coeffs = [a, b, c], f(x) = a*T_0 + b*T_1 + c*T_2 = a + b*x + c*(2x^2-1)
        // With [1, 0, 1]: f(x) = 1 + 0 + (2x^2 - 1) = 2x^2
        let coeffs = [1.0, 0.0, 1.0];
        let result = evaluate_chebyshev(&coeffs, 3, 0.5).unwrap();
        // Expected: 2 * 0.25 = 0.5
        assert!((result - 0.5).abs() < 1e-14);

        let result = evaluate_chebyshev(&coeffs, 3, 1.0).unwrap();
        // Expected: 2 * 1 = 2
        assert!((result - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_chebyshev_derivative_quadratic() {
        // f(x) = a + b*x + c*(2x^2-1)
        // f'(x) = b + 4cx
        // With [1, 0, 1] and half_interval=1: f'(x) = 0 + 4*1*x = 4x
        let coeffs = [1.0, 0.0, 1.0];
        let result = evaluate_chebyshev_derivative(&coeffs, 3, 0.5, 1.0).unwrap();
        // Expected: 4 * 0.5 = 2
        assert!((result - 2.0).abs() < 1e-14);
    }
}
