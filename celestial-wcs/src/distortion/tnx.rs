use crate::common::newton_raphson_2d;
use crate::error::{WcsError, WcsResult};

use super::polynomial::{chebyshev, legendre};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurfaceType {
    Chebyshev = 1,
    Legendre = 2,
    Polynomial = 3,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CrossTerms {
    None = 1,
    Half = 2,
    Full = 3,
}

#[derive(Debug, Clone)]
pub struct TnxSurface {
    surface_type: SurfaceType,
    x_order: u32,
    y_order: u32,
    cross_terms: CrossTerms,
    x_range: (f64, f64),
    y_range: (f64, f64),
    coefficients: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct TnxDistortion {
    lng_surface: TnxSurface,
    lat_surface: TnxSurface,
}

impl SurfaceType {
    fn from_value(val: u32) -> WcsResult<Self> {
        match val {
            1 => Ok(Self::Chebyshev),
            2 => Ok(Self::Legendre),
            3 => Ok(Self::Polynomial),
            _ => Err(WcsError::invalid_parameter(format!(
                "invalid TNX surface type: {}",
                val
            ))),
        }
    }
}

impl CrossTerms {
    fn from_value(val: u32) -> WcsResult<Self> {
        match val {
            1 => Ok(Self::None),
            2 => Ok(Self::Half),
            3 => Ok(Self::Full),
            _ => Err(WcsError::invalid_parameter(format!(
                "invalid TNX cross-terms type: {}",
                val
            ))),
        }
    }
}

impl TnxSurface {
    pub fn new(
        surface_type: SurfaceType,
        x_order: u32,
        y_order: u32,
        cross_terms: CrossTerms,
        x_range: (f64, f64),
        y_range: (f64, f64),
        coefficients: Vec<f64>,
    ) -> WcsResult<Self> {
        let expected = Self::expected_coeffs(x_order, y_order, cross_terms);
        if coefficients.len() != expected {
            return Err(WcsError::invalid_parameter(format!(
                "TNX surface expects {} coefficients, got {}",
                expected,
                coefficients.len()
            )));
        }
        Ok(Self {
            surface_type,
            x_order,
            y_order,
            cross_terms,
            x_range,
            y_range,
            coefficients,
        })
    }

    fn expected_coeffs(x_order: u32, y_order: u32, cross_terms: CrossTerms) -> usize {
        match cross_terms {
            CrossTerms::None => (x_order + y_order) as usize,
            CrossTerms::Half => Self::half_cross_count(x_order, y_order),
            CrossTerms::Full => (x_order * y_order) as usize,
        }
    }

    fn half_cross_count(x_order: u32, y_order: u32) -> usize {
        let max_order = x_order.max(y_order);
        let mut count = 0;
        for j in 0..y_order {
            for i in 0..x_order {
                if i + j <= max_order {
                    count += 1;
                }
            }
        }
        count
    }

    #[inline]
    fn normalize_x(&self, x: f64) -> f64 {
        let (xmin, xmax) = self.x_range;
        (2.0 * x - (xmax + xmin)) / (xmax - xmin)
    }

    #[inline]
    fn normalize_y(&self, y: f64) -> f64 {
        let (ymin, ymax) = self.y_range;
        (2.0 * y - (ymax + ymin)) / (ymax - ymin)
    }

    fn basis(&self, order: u32, x_norm: f64) -> f64 {
        match self.surface_type {
            SurfaceType::Chebyshev => chebyshev(order as usize, x_norm),
            SurfaceType::Legendre => legendre(order as usize, x_norm),
            SurfaceType::Polynomial => x_norm.powi(order as i32),
        }
    }

    pub fn evaluate(&self, x: f64, y: f64) -> f64 {
        let x_norm = self.normalize_x(x);
        let y_norm = self.normalize_y(y);

        match self.cross_terms {
            CrossTerms::None => self.evaluate_no_cross(x_norm, y_norm),
            CrossTerms::Half => self.evaluate_half_cross(x_norm, y_norm),
            CrossTerms::Full => self.evaluate_full_cross(x_norm, y_norm),
        }
    }

    fn evaluate_no_cross(&self, x_norm: f64, y_norm: f64) -> f64 {
        let mut result = 0.0;
        let mut idx = 0;

        for i in 0..self.x_order {
            result += self.coefficients[idx] * self.basis(i, x_norm);
            idx += 1;
        }
        for j in 0..self.y_order {
            result += self.coefficients[idx] * self.basis(j, y_norm);
            idx += 1;
        }
        result
    }

    fn evaluate_half_cross(&self, x_norm: f64, y_norm: f64) -> f64 {
        let max_order = self.x_order.max(self.y_order);
        let mut result = 0.0;
        let mut idx = 0;

        for j in 0..self.y_order {
            let y_basis = self.basis(j, y_norm);
            for i in 0..self.x_order {
                if i + j <= max_order {
                    let x_basis = self.basis(i, x_norm);
                    result += self.coefficients[idx] * x_basis * y_basis;
                    idx += 1;
                }
            }
        }
        result
    }

    fn evaluate_full_cross(&self, x_norm: f64, y_norm: f64) -> f64 {
        let mut result = 0.0;
        let mut idx = 0;

        for j in 0..self.y_order {
            let y_basis = self.basis(j, y_norm);
            for i in 0..self.x_order {
                let x_basis = self.basis(i, x_norm);
                result += self.coefficients[idx] * x_basis * y_basis;
                idx += 1;
            }
        }
        result
    }

    pub fn parse(content: &str) -> WcsResult<Self> {
        let tokens: Vec<&str> = content.split_whitespace().collect();

        if tokens.len() < 8 {
            return Err(WcsError::invalid_parameter(
                "TNX surface requires at least 8 values",
            ));
        }

        let surface_type = Self::parse_u32(tokens[0], "surface_type")?;
        let x_order = Self::parse_u32(tokens[1], "xorder")?;
        let y_order = Self::parse_u32(tokens[2], "yorder")?;
        let xterms = Self::parse_u32(tokens[3], "xterms")?;
        let xmin = Self::parse_f64(tokens[4], "xmin")?;
        let xmax = Self::parse_f64(tokens[5], "xmax")?;
        let ymin = Self::parse_f64(tokens[6], "ymin")?;
        let ymax = Self::parse_f64(tokens[7], "ymax")?;

        let coefficients: Result<Vec<f64>, _> = tokens[8..]
            .iter()
            .enumerate()
            .map(|(i, s)| Self::parse_f64(s, &format!("coefficient[{}]", i)))
            .collect();

        Self::new(
            SurfaceType::from_value(surface_type)?,
            x_order,
            y_order,
            CrossTerms::from_value(xterms)?,
            (xmin, xmax),
            (ymin, ymax),
            coefficients?,
        )
    }

    fn parse_u32(s: &str, name: &str) -> WcsResult<u32> {
        s.parse::<f64>()
            .map(|v| v as u32)
            .map_err(|_| WcsError::invalid_parameter(format!("invalid {}: '{}'", name, s)))
    }

    fn parse_f64(s: &str, name: &str) -> WcsResult<f64> {
        s.parse()
            .map_err(|_| WcsError::invalid_parameter(format!("invalid {}: '{}'", name, s)))
    }
}

impl TnxDistortion {
    pub fn new(lng_surface: TnxSurface, lat_surface: TnxSurface) -> Self {
        Self {
            lng_surface,
            lat_surface,
        }
    }

    pub fn parse(wat1: &str, wat2: &str) -> WcsResult<Self> {
        let lng_content = Self::extract_correction(wat1, "lngcor")?;
        let lat_content = Self::extract_correction(wat2, "latcor")?;

        let lng_surface = TnxSurface::parse(&lng_content)?;
        let lat_surface = TnxSurface::parse(&lat_content)?;

        Ok(Self::new(lng_surface, lat_surface))
    }

    fn extract_correction(wat: &str, key: &str) -> WcsResult<String> {
        let pattern = format!("{} = \"", key);
        let start = wat.find(&pattern).ok_or_else(|| {
            WcsError::missing_keyword(format!("TNX {} not found in WAT string", key))
        })?;

        let after_key = &wat[start + pattern.len()..];
        let end = after_key
            .find('"')
            .ok_or_else(|| WcsError::invalid_parameter(format!("unterminated {} string", key)))?;

        Ok(after_key[..end].to_string())
    }

    pub fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        let dx = self.lng_surface.evaluate(x, y);
        let dy = self.lat_surface.evaluate(x, y);
        (x + dx, y + dy)
    }

    pub fn apply_inverse(&self, x: f64, y: f64) -> WcsResult<(f64, f64)> {
        let distort_fn = |px: f64, py: f64| self.apply(px, py);

        newton_raphson_2d((x, y), (x, y), distort_fn, 20, 1e-12, "TNX inverse distortion")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zero_surface(cross_terms: CrossTerms) -> TnxSurface {
        let (x_order, y_order) = (3, 3);
        let n = TnxSurface::expected_coeffs(x_order, y_order, cross_terms);
        TnxSurface::new(
            SurfaceType::Chebyshev,
            x_order,
            y_order,
            cross_terms,
            (0.0, 100.0),
            (0.0, 100.0),
            vec![0.0; n],
        )
        .unwrap()
    }

    #[test]
    fn test_zero_surface_evaluates_to_zero_for_every_cross_term_mode() {
        // A zero-coefficient surface must produce 0 regardless of cross-term
        // selection.  The coefficient count differs per mode (Full=9, Half=6,
        // None=6 for x_order=y_order=3); zero_surface() handles that.
        for mode in [CrossTerms::Full, CrossTerms::Half, CrossTerms::None] {
            let surf = zero_surface(mode);
            assert_eq!(surf.evaluate(0.0, 0.0), 0.0, "{:?} at origin", mode);
            assert_eq!(surf.evaluate(50.0, 50.0), 0.0, "{:?} at center", mode);
            assert_eq!(surf.evaluate(100.0, 100.0), 0.0, "{:?} at corner", mode);
        }
    }

    #[test]
    fn test_basis_functions_match_explicit_polynomials() {
        // Surface.basis(n, x) must agree with the explicit polynomial forms
        // for n = 0..3 across all three basis types.
        struct Case {
            ty: SurfaceType,
            x: f64,
            expected: [f64; 4],
        }
        let x = 0.5;
        let cases = [
            // Chebyshev: T_0=1, T_1=x, T_2=2x^2-1, T_3=4x^3-3x.
            Case {
                ty: SurfaceType::Chebyshev,
                x,
                expected: [1.0, 0.5, 2.0 * 0.25 - 1.0, 4.0 * 0.125 - 3.0 * 0.5],
            },
            // Legendre: P_0=1, P_1=x, P_2=(3x^2-1)/2, P_3=(5x^3-3x)/2.
            Case {
                ty: SurfaceType::Legendre,
                x,
                expected: [
                    1.0, 0.5, (3.0 * 0.25 - 1.0) / 2.0, (5.0 * 0.125 - 3.0 * 0.5) / 2.0,
                ],
            },
            // Polynomial: x^n.
            Case {
                ty: SurfaceType::Polynomial,
                x,
                expected: [1.0, 0.5, 0.25, 0.125],
            },
        ];
        for case in &cases {
            let surf = TnxSurface::new(
                case.ty,
                3, 3,
                CrossTerms::Full,
                (0.0, 1.0), (0.0, 1.0),
                vec![0.0; 9],
            ).unwrap();
            for (n, &expected) in case.expected.iter().enumerate() {
                assert!(
                    (surf.basis(n as u32, case.x) - expected).abs() < 1e-14,
                    "{:?} basis {} at {} mismatch", case.ty, n, case.x,
                );
            }
        }
    }

    #[test]
    fn test_normalization_maps_range_to_unit_interval() {
        // normalize_x/normalize_y must map x_range/y_range onto [-1, +1].
        let surf = TnxSurface::new(
            SurfaceType::Chebyshev,
            2, 2,
            CrossTerms::Full,
            (10.0, 20.0), (30.0, 50.0),
            vec![0.0; 4],
        ).unwrap();
        // Centre of the asymmetric (10, 20) x (30, 50) range maps to 0.
        assert_eq!(surf.normalize_x(15.0), 0.0);
        assert_eq!(surf.normalize_y(40.0), 0.0);

        // Edges of the (0, 100) range map to -1 and +1.
        let surf = TnxSurface::new(
            SurfaceType::Chebyshev,
            2, 2,
            CrossTerms::Full,
            (0.0, 100.0), (0.0, 100.0),
            vec![0.0; 4],
        ).unwrap();
        assert_eq!(surf.normalize_x(0.0), -1.0);
        assert_eq!(surf.normalize_x(100.0), 1.0);
        assert_eq!(surf.normalize_y(0.0), -1.0);
        assert_eq!(surf.normalize_y(100.0), 1.0);
    }

    #[test]
    fn test_polynomial_surface_evaluation_known_terms() {
        // Constant surface (5) is invariant across the domain.
        let constant = TnxSurface::new(
            SurfaceType::Polynomial,
            2, 2,
            CrossTerms::Full,
            (0.0, 100.0), (0.0, 100.0),
            vec![5.0, 0.0, 0.0, 0.0],
        ).unwrap();
        assert_eq!(constant.evaluate(0.0, 0.0), 5.0);
        assert_eq!(constant.evaluate(50.0, 50.0), 5.0);
        assert_eq!(constant.evaluate(100.0, 100.0), 5.0);

        // Linear-in-x surface produces normalized x in [-1, 1] (independent of y).
        let linear_x = TnxSurface::new(
            SurfaceType::Polynomial,
            2, 2,
            CrossTerms::Full,
            (0.0, 100.0), (0.0, 100.0),
            vec![0.0, 1.0, 0.0, 0.0],
        ).unwrap();
        assert_eq!(linear_x.evaluate(0.0, 50.0), -1.0);
        assert_eq!(linear_x.evaluate(50.0, 50.0), 0.0);
        assert_eq!(linear_x.evaluate(100.0, 50.0), 1.0);
    }

    #[test]
    fn test_parse_surface() {
        let content = "3 2 2 3 0.0 100.0 0.0 100.0 1.0 0.0 0.0 0.0";
        let surf = TnxSurface::parse(content).unwrap();

        assert_eq!(surf.surface_type, SurfaceType::Polynomial);
        assert_eq!(surf.x_order, 2);
        assert_eq!(surf.y_order, 2);
        assert_eq!(surf.cross_terms, CrossTerms::Full);
        assert_eq!(surf.x_range, (0.0, 100.0));
        assert_eq!(surf.y_range, (0.0, 100.0));
        assert_eq!(surf.coefficients.len(), 4);
    }

    #[test]
    fn test_parse_wat_string() {
        let wat1 = "wtype=tnx axtype=ra lngcor = \"3 2 2 3 0 100 0 100 0 0 0 0\"";
        let wat2 = "wtype=tnx axtype=dec latcor = \"3 2 2 3 0 100 0 100 0 0 0 0\"";

        let tnx = TnxDistortion::parse(wat1, wat2).unwrap();
        let (x, y) = tnx.apply(50.0, 50.0);

        assert_eq!(x, 50.0);
        assert_eq!(y, 50.0);
    }

    #[test]
    fn test_roundtrip_across_surface_types() {
        // All three basis types share the same Newton-Raphson inverse, so a
        // single sweep covers the polynomial, Chebyshev and Legendre paths.
        let cases: &[(SurfaceType, [f64; 9], &[(f64, f64)])] = &[
            (
                SurfaceType::Polynomial,
                [0.0, 0.001, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0],
                &[(45.0, 55.0)],
            ),
            (
                SurfaceType::Chebyshev,
                [0.0, 0.0005, 0.0, 0.0, 0.0, 0.0003, 0.0, 0.0, 0.0],
                &[(25.0, 25.0), (50.0, 75.0), (80.0, 20.0)],
            ),
            (
                SurfaceType::Legendre,
                [0.0, 0.0004, 0.0, 0.0, 0.0, 0.0002, 0.0, 0.0, 0.0],
                &[(60.0, 40.0)],
            ),
        ];
        for (ty, coeffs, pts) in cases {
            let lng = TnxSurface::new(
                *ty, 3, 3, CrossTerms::Full,
                (0.0, 100.0), (0.0, 100.0),
                coeffs.to_vec(),
            ).unwrap();
            let lat = TnxSurface::new(
                *ty, 3, 3, CrossTerms::Full,
                (0.0, 100.0), (0.0, 100.0),
                coeffs.to_vec(),
            ).unwrap();
            let tnx = TnxDistortion::new(lng, lat);
            for &(x_orig, y_orig) in *pts {
                let (xd, yd) = tnx.apply(x_orig, y_orig);
                let (xb, yb) = tnx.apply_inverse(xd, yd).unwrap();
                assert!(
                    (xb - x_orig).abs() < 1e-10,
                    "{:?} x at ({}, {})", ty, x_orig, y_orig,
                );
                assert!(
                    (yb - y_orig).abs() < 1e-10,
                    "{:?} y at ({}, {})", ty, x_orig, y_orig,
                );
            }
        }
    }

    #[test]
    fn test_cross_terms_coefficient_count_and_evaluation() {
        // expected_coeffs covers every cross-term mode:
        //   Full = x_order * y_order
        //   None = max(x_order, y_order) (one-axis terms only)
        //   Half = upper-triangular total-degree polynomial
        assert_eq!(TnxSurface::expected_coeffs(2, 2, CrossTerms::Full), 4);
        assert_eq!(TnxSurface::expected_coeffs(3, 3, CrossTerms::Full), 9);
        assert_eq!(TnxSurface::expected_coeffs(4, 4, CrossTerms::Full), 16);
        assert_eq!(TnxSurface::expected_coeffs(2, 3, CrossTerms::Full), 6);
        assert_eq!(TnxSurface::expected_coeffs(2, 2, CrossTerms::None), 4);
        assert_eq!(TnxSurface::expected_coeffs(3, 3, CrossTerms::None), 6);
        assert_eq!(TnxSurface::expected_coeffs(4, 5, CrossTerms::None), 9);

        // A constant surface with Half cross terms still evaluates to its
        // constant - the cross-term mode affects only the term layout.
        let half_n = TnxSurface::expected_coeffs(3, 3, CrossTerms::Half);
        let mut half_coeffs = vec![0.0; half_n];
        half_coeffs[0] = 1.0;
        let half = TnxSurface::new(
            SurfaceType::Polynomial,
            3, 3,
            CrossTerms::Half,
            (0.0, 100.0), (0.0, 100.0),
            half_coeffs,
        ).unwrap();
        assert_eq!(half.evaluate(50.0, 50.0), 1.0);

        // None-cross polynomial: const + 0.5*x normalised, evaluated at the
        // right edge where normalised x = 1.
        let none = TnxSurface::new(
            SurfaceType::Polynomial,
            3, 3,
            CrossTerms::None,
            (0.0, 100.0), (0.0, 100.0),
            vec![1.0, 0.5, 0.0, 0.0, 0.0, 0.0],
        ).unwrap();
        assert_eq!(none.evaluate(100.0, 50.0), 1.0 + 0.5);
    }

    #[test]
    fn test_surface_construction_validates_inputs() {
        // Unknown surface-type codes and cross-term codes are rejected.
        assert!(SurfaceType::from_value(0).is_err());
        assert!(SurfaceType::from_value(4).is_err());
        assert!(CrossTerms::from_value(0).is_err());
        assert!(CrossTerms::from_value(4).is_err());

        // Wrong coefficient count (expected 4, given 3) is rejected.
        let bad_n = TnxSurface::new(
            SurfaceType::Polynomial,
            2, 2,
            CrossTerms::Full,
            (0.0, 100.0), (0.0, 100.0),
            vec![0.0; 3],
        );
        assert!(bad_n.is_err());
    }

    #[test]
    fn test_parse_rejects_malformed_wat_strings() {
        // Missing lngcor in WAT1.
        let wat1 = "wtype=tnx axtype=ra";
        let wat2 = "wtype=tnx axtype=dec latcor = \"3 2 2 3 0 100 0 100 0 0 0 0\"";
        assert!(TnxDistortion::parse(wat1, wat2).is_err());

        // Missing latcor in WAT2.
        let wat1 = "wtype=tnx axtype=ra lngcor = \"3 2 2 3 0 100 0 100 0 0 0 0\"";
        let wat2 = "wtype=tnx axtype=dec";
        assert!(TnxDistortion::parse(wat1, wat2).is_err());

        // Fewer than the 8 required header tokens.
        let err = TnxSurface::parse("1 2 3 4 5 6 7").unwrap_err().to_string();
        assert!(err.contains("at least 8"), "got: {}", err);

        // Unterminated quoted correction string.
        let wat1 = "wtype=tnx lngcor = \"1 2 3 4 5 6 7 8";
        let wat2 = "wtype=tnx latcor = \"1 2 3 4 5 6 7 8 0.0\"";
        let err = TnxDistortion::parse(wat1, wat2).unwrap_err().to_string();
        assert!(err.contains("unterminated"), "got: {}", err);
    }

    #[test]
    fn test_legendre_surface_evaluation() {
        let coeffs = vec![0.0, 0.1, 0.0, 0.0];
        let surface = TnxSurface::new(
            SurfaceType::Legendre,
            2,
            2,
            CrossTerms::Full,
            (-1.0, 1.0),
            (-1.0, 1.0),
            coeffs,
        )
        .unwrap();

        let val = surface.evaluate(0.5, 0.0);
        assert!((val - 0.05).abs() < 1e-10);
    }

}
