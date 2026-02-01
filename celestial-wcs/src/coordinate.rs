use celestial_core::Angle;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PixelCoord {
    x: f64,
    y: f64,
}

impl PixelCoord {
    #[inline]
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    #[inline]
    pub fn x(&self) -> f64 {
        self.x
    }

    #[inline]
    pub fn y(&self) -> f64 {
        self.y
    }

    #[inline]
    pub fn to_array_index(&self) -> (usize, usize) {
        let row = (self.y - 1.0).round() as usize;
        let col = (self.x - 1.0).round() as usize;
        (row, col)
    }

    #[inline]
    pub fn from_array_index(row: usize, col: usize) -> Self {
        Self {
            x: col as f64 + 1.0,
            y: row as f64 + 1.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntermediateCoord {
    x: f64,
    y: f64,
}

impl IntermediateCoord {
    #[inline]
    pub fn new(x_deg: f64, y_deg: f64) -> Self {
        Self { x: x_deg, y: y_deg }
    }

    #[inline]
    pub fn x_deg(&self) -> f64 {
        self.x
    }

    #[inline]
    pub fn y_deg(&self) -> f64 {
        self.y
    }

    #[inline]
    pub fn x_rad(&self) -> f64 {
        self.x.to_radians()
    }

    #[inline]
    pub fn y_rad(&self) -> f64 {
        self.y.to_radians()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NativeCoord {
    phi: Angle,
    theta: Angle,
}

impl NativeCoord {
    #[inline]
    pub fn new(phi: Angle, theta: Angle) -> Self {
        Self { phi, theta }
    }

    #[inline]
    pub fn phi(&self) -> Angle {
        self.phi
    }

    #[inline]
    pub fn theta(&self) -> Angle {
        self.theta
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CelestialCoord {
    alpha: Angle,
    delta: Angle,
}

impl CelestialCoord {
    #[inline]
    pub fn new(alpha: Angle, delta: Angle) -> Self {
        Self { alpha, delta }
    }

    #[inline]
    pub fn alpha(&self) -> Angle {
        self.alpha
    }

    #[inline]
    pub fn delta(&self) -> Angle {
        self.delta
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pixel_coord_new_and_accessors() {
        let p = PixelCoord::new(100.5, 200.5);
        assert_eq!(p.x(), 100.5);
        assert_eq!(p.y(), 200.5);
    }

    #[test]
    fn test_pixel_to_array_index() {
        let p = PixelCoord::new(1.0, 1.0);
        assert_eq!(p.to_array_index(), (0, 0));

        let p2 = PixelCoord::new(10.0, 20.0);
        assert_eq!(p2.to_array_index(), (19, 9));
    }

    #[test]
    fn test_array_index_to_pixel() {
        let p = PixelCoord::from_array_index(0, 0);
        assert_eq!(p.x(), 1.0);
        assert_eq!(p.y(), 1.0);

        let p2 = PixelCoord::from_array_index(19, 9);
        assert_eq!(p2.x(), 10.0);
        assert_eq!(p2.y(), 20.0);
    }

    #[test]
    fn test_pixel_roundtrip() {
        let p = PixelCoord::new(50.0, 100.0);
        let (row, col) = p.to_array_index();
        let p2 = PixelCoord::from_array_index(row, col);
        assert_eq!(p, p2);
    }

    #[test]
    fn test_intermediate_coord() {
        let c = IntermediateCoord::new(0.001, -0.002);
        assert_eq!(c.x_deg(), 0.001);
        assert_eq!(c.y_deg(), -0.002);
        assert!((c.x_rad() - 0.001_f64.to_radians()).abs() < 1e-15);
        assert!((c.y_rad() - (-0.002_f64).to_radians()).abs() < 1e-15);
    }

    #[test]
    fn test_native_coord() {
        let phi = Angle::from_degrees(45.0);
        let theta = Angle::from_degrees(30.0);
        let n = NativeCoord::new(phi, theta);
        assert_eq!(n.phi(), phi);
        assert_eq!(n.theta(), theta);
    }

    #[test]
    fn test_celestial_coord() {
        let alpha = Angle::from_degrees(180.0);
        let delta = Angle::from_degrees(45.0);
        let c = CelestialCoord::new(alpha, delta);
        assert_eq!(c.alpha(), alpha);
        assert_eq!(c.delta(), delta);
    }
}
