use celestial_core::angle::Angle;

use crate::error::{WcsError, WcsResult};

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
    pub fn to_array_index(&self) -> WcsResult<(usize, usize)> {
        Ok((array_index(self.y)?, array_index(self.x)?))
    }

    #[inline]
    pub fn from_array_index(row: usize, col: usize) -> Self {
        Self {
            x: col as f64 + 1.0,
            y: row as f64 + 1.0,
        }
    }
}

// FITS pixel n is centred on n and covers array index n − 1.
fn array_index(pixel: f64) -> WcsResult<usize> {
    let index = libm::round(pixel - 1.0);
    if !(0.0..usize::MAX as f64).contains(&index) {
        let message = format!("pixel {pixel} has no array index");
        return Err(WcsError::out_of_bounds(message));
    }
    Ok(index as usize)
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
}

// Paper II gives native and celestial coordinates in degrees, so they are
// stored that way; an `Angle` is made only when a caller asks for one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NativeCoord {
    phi: f64,
    theta: f64,
}

impl NativeCoord {
    #[inline]
    pub fn new(phi_deg: f64, theta_deg: f64) -> Self {
        Self {
            phi: phi_deg,
            theta: theta_deg,
        }
    }

    #[inline]
    pub fn phi_deg(&self) -> f64 {
        self.phi
    }

    #[inline]
    pub fn theta_deg(&self) -> f64 {
        self.theta
    }

    #[inline]
    pub fn phi(&self) -> Angle {
        Angle::from_degrees(self.phi)
    }

    #[inline]
    pub fn theta(&self) -> Angle {
        Angle::from_degrees(self.theta)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CelestialCoord {
    alpha: f64,
    delta: f64,
}

impl CelestialCoord {
    #[inline]
    pub fn new(alpha_deg: f64, delta_deg: f64) -> Self {
        Self {
            alpha: alpha_deg,
            delta: delta_deg,
        }
    }

    #[inline]
    pub fn alpha_deg(&self) -> f64 {
        self.alpha
    }

    #[inline]
    pub fn delta_deg(&self) -> f64 {
        self.delta
    }

    #[inline]
    pub fn alpha(&self) -> Angle {
        Angle::from_degrees(self.alpha)
    }

    #[inline]
    pub fn delta(&self) -> Angle {
        Angle::from_degrees(self.delta)
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
    fn test_pixel_array_index_conversion() {
        // Forward direction: 1-indexed FITS pixel -> 0-indexed (row, col).
        assert_eq!(PixelCoord::new(1.0, 1.0).to_array_index().unwrap(), (0, 0));
        assert_eq!(
            PixelCoord::new(10.0, 20.0).to_array_index().unwrap(),
            (19, 9)
        );
        assert_eq!(PixelCoord::new(0.6, 1.4).to_array_index().unwrap(), (0, 0));

        // Inverse direction: (row, col) -> FITS pixel.
        assert_eq!(
            PixelCoord::from_array_index(0, 0),
            PixelCoord::new(1.0, 1.0)
        );
        assert_eq!(
            PixelCoord::from_array_index(19, 9),
            PixelCoord::new(10.0, 20.0)
        );

        // Round trip preserves the original pixel.
        let p = PixelCoord::new(50.0, 100.0);
        let (row, col) = p.to_array_index().unwrap();
        assert_eq!(PixelCoord::from_array_index(row, col), p);
    }

    #[test]
    fn test_pixel_off_the_array_has_no_index() {
        let off = [
            (f64::NAN, 1.0),
            (1.0, f64::INFINITY),
            (0.4, 1.0),
            (1.0, -3.0),
            (1e20, 1.0),
        ];
        for (x, y) in off {
            let result = PixelCoord::new(x, y).to_array_index();
            assert!(
                matches!(result, Err(WcsError::OutOfBounds { .. })),
                "({x}, {y}): {result:?}"
            );
        }
    }

    #[test]
    fn test_intermediate_coord() {
        let c = IntermediateCoord::new(0.001, -0.002);
        assert_eq!(c.x_deg(), 0.001);
        assert_eq!(c.y_deg(), -0.002);
    }

    #[test]
    fn test_native_coord() {
        let n = NativeCoord::new(45.0, 30.0);
        assert_eq!(n.phi_deg(), 45.0);
        assert_eq!(n.theta_deg(), 30.0);
        assert_eq!(n.phi(), Angle::from_degrees(45.0));
        assert_eq!(n.theta(), Angle::from_degrees(30.0));
    }

    #[test]
    fn test_celestial_coord_keeps_degrees_exact() {
        // 123.456° does not survive a trip through radians.
        let c = CelestialCoord::new(123.456, -45.0);
        assert_eq!(c.alpha_deg(), 123.456);
        assert_eq!(c.delta_deg(), -45.0);
        assert_eq!(c.alpha(), Angle::from_degrees(123.456));
        assert_eq!(c.delta(), Angle::from_degrees(-45.0));
    }
}
