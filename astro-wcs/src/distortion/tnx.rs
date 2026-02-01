use crate::error::{WcsError, WcsResult};

use super::polynomial::{chebyshev, legendre, newton_raphson_2d};

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

        newton_raphson_2d((x, y), (x, y), distort_fn, 20, 1e-12).map_err(|msg| {
            WcsError::convergence_failure(format!("TNX inverse distortion: {}", msg))
        })
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
    fn test_identity_full_cross() {
        let surf = zero_surface(CrossTerms::Full);
        assert_eq!(surf.evaluate(50.0, 50.0), 0.0);
        assert_eq!(surf.evaluate(0.0, 0.0), 0.0);
        assert_eq!(surf.evaluate(100.0, 100.0), 0.0);
    }

    #[test]
    fn test_identity_half_cross() {
        let surf = zero_surface(CrossTerms::Half);
        assert_eq!(surf.evaluate(50.0, 50.0), 0.0);
    }

    #[test]
    fn test_identity_no_cross() {
        let surf = zero_surface(CrossTerms::None);
        assert_eq!(surf.evaluate(50.0, 50.0), 0.0);
    }

    #[test]
    fn test_chebyshev_basis() {
        let surf = TnxSurface::new(
            SurfaceType::Chebyshev,
            3,
            3,
            CrossTerms::Full,
            (0.0, 1.0),
            (0.0, 1.0),
            vec![0.0; 9],
        )
        .unwrap();

        assert_eq!(surf.basis(0, 0.5), 1.0);
        assert_eq!(surf.basis(1, 0.5), 0.5);
        assert_eq!(surf.basis(2, 0.5), 2.0 * 0.25 - 1.0);
    }

    #[test]
    fn test_legendre_basis() {
        let surf = TnxSurface::new(
            SurfaceType::Legendre,
            3,
            3,
            CrossTerms::Full,
            (0.0, 1.0),
            (0.0, 1.0),
            vec![0.0; 9],
        )
        .unwrap();

        assert_eq!(surf.basis(0, 0.5), 1.0);
        assert_eq!(surf.basis(1, 0.5), 0.5);
        let expected_p2 = (3.0 * 0.25 - 1.0) / 2.0;
        assert_eq!(surf.basis(2, 0.5), expected_p2);
    }

    #[test]
    fn test_polynomial_basis() {
        let surf = TnxSurface::new(
            SurfaceType::Polynomial,
            3,
            3,
            CrossTerms::Full,
            (0.0, 1.0),
            (0.0, 1.0),
            vec![0.0; 9],
        )
        .unwrap();

        assert_eq!(surf.basis(0, 0.5), 1.0);
        assert_eq!(surf.basis(1, 0.5), 0.5);
        assert_eq!(surf.basis(2, 0.5), 0.25);
        assert_eq!(surf.basis(3, 0.5), 0.125);
    }

    #[test]
    fn test_normalization_center() {
        let surf = TnxSurface::new(
            SurfaceType::Chebyshev,
            2,
            2,
            CrossTerms::Full,
            (10.0, 20.0),
            (30.0, 50.0),
            vec![0.0; 4],
        )
        .unwrap();

        assert_eq!(surf.normalize_x(15.0), 0.0);
        assert_eq!(surf.normalize_y(40.0), 0.0);
    }

    #[test]
    fn test_normalization_edges() {
        let surf = TnxSurface::new(
            SurfaceType::Chebyshev,
            2,
            2,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            vec![0.0; 4],
        )
        .unwrap();

        assert_eq!(surf.normalize_x(0.0), -1.0);
        assert_eq!(surf.normalize_x(100.0), 1.0);
        assert_eq!(surf.normalize_y(0.0), -1.0);
        assert_eq!(surf.normalize_y(100.0), 1.0);
    }

    #[test]
    fn test_constant_surface() {
        let coeffs = vec![5.0, 0.0, 0.0, 0.0];
        let surf = TnxSurface::new(
            SurfaceType::Polynomial,
            2,
            2,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs,
        )
        .unwrap();

        assert_eq!(surf.evaluate(0.0, 0.0), 5.0);
        assert_eq!(surf.evaluate(50.0, 50.0), 5.0);
        assert_eq!(surf.evaluate(100.0, 100.0), 5.0);
    }

    #[test]
    fn test_linear_x_surface() {
        let coeffs = vec![0.0, 1.0, 0.0, 0.0];
        let surf = TnxSurface::new(
            SurfaceType::Polynomial,
            2,
            2,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs,
        )
        .unwrap();

        assert_eq!(surf.evaluate(0.0, 50.0), -1.0);
        assert_eq!(surf.evaluate(50.0, 50.0), 0.0);
        assert_eq!(surf.evaluate(100.0, 50.0), 1.0);
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
    fn test_roundtrip_polynomial() {
        let coeffs = vec![0.0, 0.001, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0];
        let lng = TnxSurface::new(
            SurfaceType::Polynomial,
            3,
            3,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs.clone(),
        )
        .unwrap();
        let lat = TnxSurface::new(
            SurfaceType::Polynomial,
            3,
            3,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs,
        )
        .unwrap();

        let tnx = TnxDistortion::new(lng, lat);

        let (x_orig, y_orig) = (45.0, 55.0);
        let (x_dist, y_dist) = tnx.apply(x_orig, y_orig);
        let (x_back, y_back) = tnx.apply_inverse(x_dist, y_dist).unwrap();

        assert!((x_back - x_orig).abs() < 1e-10);
        assert!((y_back - y_orig).abs() < 1e-10);
    }

    #[test]
    fn test_roundtrip_chebyshev() {
        let coeffs = vec![0.0, 0.0005, 0.0, 0.0, 0.0, 0.0003, 0.0, 0.0, 0.0];
        let lng = TnxSurface::new(
            SurfaceType::Chebyshev,
            3,
            3,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs.clone(),
        )
        .unwrap();
        let lat = TnxSurface::new(
            SurfaceType::Chebyshev,
            3,
            3,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs,
        )
        .unwrap();

        let tnx = TnxDistortion::new(lng, lat);

        for (x_orig, y_orig) in [(25.0, 25.0), (50.0, 75.0), (80.0, 20.0)] {
            let (x_dist, y_dist) = tnx.apply(x_orig, y_orig);
            let (x_back, y_back) = tnx.apply_inverse(x_dist, y_dist).unwrap();

            assert!((x_back - x_orig).abs() < 1e-10);
            assert!((y_back - y_orig).abs() < 1e-10);
        }
    }

    #[test]
    fn test_roundtrip_legendre() {
        let coeffs = vec![0.0, 0.0004, 0.0, 0.0, 0.0, 0.0002, 0.0, 0.0, 0.0];
        let lng = TnxSurface::new(
            SurfaceType::Legendre,
            3,
            3,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs.clone(),
        )
        .unwrap();
        let lat = TnxSurface::new(
            SurfaceType::Legendre,
            3,
            3,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs,
        )
        .unwrap();

        let tnx = TnxDistortion::new(lng, lat);

        let (x_orig, y_orig) = (60.0, 40.0);
        let (x_dist, y_dist) = tnx.apply(x_orig, y_orig);
        let (x_back, y_back) = tnx.apply_inverse(x_dist, y_dist).unwrap();

        assert!((x_back - x_orig).abs() < 1e-10);
        assert!((y_back - y_orig).abs() < 1e-10);
    }

    #[test]
    fn test_half_cross_terms() {
        let n = TnxSurface::expected_coeffs(3, 3, CrossTerms::Half);
        let mut coeffs = vec![0.0; n];
        coeffs[0] = 1.0;

        let surf = TnxSurface::new(
            SurfaceType::Polynomial,
            3,
            3,
            CrossTerms::Half,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs,
        )
        .unwrap();

        assert_eq!(surf.evaluate(50.0, 50.0), 1.0);
    }

    #[test]
    fn test_no_cross_terms() {
        let n = TnxSurface::expected_coeffs(3, 3, CrossTerms::None);
        assert_eq!(n, 6);

        let coeffs = vec![1.0, 0.5, 0.0, 0.0, 0.0, 0.0];
        let surf = TnxSurface::new(
            SurfaceType::Polynomial,
            3,
            3,
            CrossTerms::None,
            (0.0, 100.0),
            (0.0, 100.0),
            coeffs,
        )
        .unwrap();

        let result = surf.evaluate(100.0, 50.0);
        let expected = 1.0 + 0.5 * 1.0;
        assert_eq!(result, expected);
    }

    #[test]
    fn test_full_cross_coefficient_count() {
        assert_eq!(TnxSurface::expected_coeffs(2, 2, CrossTerms::Full), 4);
        assert_eq!(TnxSurface::expected_coeffs(3, 3, CrossTerms::Full), 9);
        assert_eq!(TnxSurface::expected_coeffs(4, 4, CrossTerms::Full), 16);
        assert_eq!(TnxSurface::expected_coeffs(2, 3, CrossTerms::Full), 6);
    }

    #[test]
    fn test_no_cross_coefficient_count() {
        assert_eq!(TnxSurface::expected_coeffs(2, 2, CrossTerms::None), 4);
        assert_eq!(TnxSurface::expected_coeffs(3, 3, CrossTerms::None), 6);
        assert_eq!(TnxSurface::expected_coeffs(4, 5, CrossTerms::None), 9);
    }

    #[test]
    fn test_invalid_surface_type() {
        let result = SurfaceType::from_value(0);
        assert!(result.is_err());

        let result = SurfaceType::from_value(4);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_cross_terms() {
        let result = CrossTerms::from_value(0);
        assert!(result.is_err());

        let result = CrossTerms::from_value(4);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_coefficient_count() {
        let result = TnxSurface::new(
            SurfaceType::Polynomial,
            2,
            2,
            CrossTerms::Full,
            (0.0, 100.0),
            (0.0, 100.0),
            vec![0.0; 3],
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_missing_lngcor() {
        let wat1 = "wtype=tnx axtype=ra";
        let wat2 = "wtype=tnx axtype=dec latcor = \"3 2 2 3 0 100 0 100 0 0 0 0\"";

        let result = TnxDistortion::parse(wat1, wat2);
        assert!(result.is_err());
    }

    #[test]
    fn test_missing_latcor() {
        let wat1 = "wtype=tnx axtype=ra lngcor = \"3 2 2 3 0 100 0 100 0 0 0 0\"";
        let wat2 = "wtype=tnx axtype=dec";

        let result = TnxDistortion::parse(wat1, wat2);
        assert!(result.is_err());
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

    #[test]
    fn test_roundtrip_legendre_distortion() {
        let coeffs = vec![0.0, 1e-4, 0.0, 0.0];
        let surface = TnxSurface::new(
            SurfaceType::Legendre,
            2,
            2,
            CrossTerms::Full,
            (-10.0, 10.0),
            (-10.0, 10.0),
            coeffs.clone(),
        )
        .unwrap();

        let tnx = TnxDistortion::new(surface.clone(), surface);
        let (x_out, y_out) = tnx.apply(5.0, 5.0);
        let (x_back, y_back) = tnx.apply_inverse(x_out, y_out).unwrap();

        assert!((x_back - 5.0).abs() < 1e-9);
        assert!((y_back - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_half_cross_terms_surface() {
        // Create and evaluate a half cross-terms surface
        let n = TnxSurface::expected_coeffs(3, 3, CrossTerms::Half);
        let mut coeffs = vec![0.0; n];
        if n > 1 {
            coeffs[1] = 0.01;
        }

        let surface = TnxSurface::new(
            SurfaceType::Chebyshev,
            3,
            3,
            CrossTerms::Half,
            (-1.0, 1.0),
            (-1.0, 1.0),
            coeffs,
        )
        .unwrap();

        let _val = surface.evaluate(0.5, 0.5);
    }

    #[test]
    fn test_parse_insufficient_tokens() {
        let result = TnxSurface::parse("1 2 3 4 5 6 7");
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("at least 8"));
    }

    #[test]
    fn test_parse_unterminated_correction() {
        let wat1 = "wtype=tnx lngcor = \"1 2 3 4 5 6 7 8";
        let wat2 = "wtype=tnx latcor = \"1 2 3 4 5 6 7 8 0.0\"";
        let result = TnxDistortion::parse(wat1, wat2);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("unterminated"));
    }
}
