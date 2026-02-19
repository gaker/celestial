use crate::coordinate::{IntermediateCoord, PixelCoord};
use crate::error::{WcsError, WcsResult};

const DETERMINANT_THRESHOLD: f64 = 1e-15;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearTransform {
    crpix: [f64; 2],
    cd: [[f64; 2]; 2],
    cd_inverse: [[f64; 2]; 2],
    determinant: f64,
}

impl LinearTransform {
    pub fn from_cd(crpix: [f64; 2], cd: [[f64; 2]; 2]) -> WcsResult<Self> {
        let determinant = cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0];
        if determinant.abs() < DETERMINANT_THRESHOLD {
            return Err(WcsError::non_invertible_matrix(determinant));
        }
        let cd_inverse = compute_inverse(cd, determinant);
        Ok(Self {
            crpix,
            cd,
            cd_inverse,
            determinant,
        })
    }

    pub fn from_pc_cdelt(crpix: [f64; 2], pc: [[f64; 2]; 2], cdelt: [f64; 2]) -> WcsResult<Self> {
        let cd = [
            [cdelt[0] * pc[0][0], cdelt[0] * pc[0][1]],
            [cdelt[1] * pc[1][0], cdelt[1] * pc[1][1]],
        ];
        Self::from_cd(crpix, cd)
    }

    pub fn pixel_to_intermediate(&self, pixel: PixelCoord) -> IntermediateCoord {
        // Paper I Eq. 1: q[i] = sum over j of m[i][j] * (p[j] - r[j])
        // Row-wise: complete each output before moving to next
        let d0 = pixel.x() - self.crpix[0];
        let d1 = pixel.y() - self.crpix[1];
        let x = self.cd[0][0] * d0 + self.cd[0][1] * d1;
        let y = self.cd[1][0] * d0 + self.cd[1][1] * d1;
        IntermediateCoord::new(x, y)
    }

    pub fn intermediate_to_pixel(&self, inter: IntermediateCoord) -> PixelCoord {
        // Row-wise: for each output, sum over input then add crpix
        let x = inter.x_deg();
        let y = inter.y_deg();
        let px = self.cd_inverse[0][0] * x + self.cd_inverse[0][1] * y + self.crpix[0];
        let py = self.cd_inverse[1][0] * x + self.cd_inverse[1][1] * y + self.crpix[1];
        PixelCoord::new(px, py)
    }

    #[inline]
    pub fn crpix(&self) -> [f64; 2] {
        self.crpix
    }

    #[inline]
    pub fn cd_matrix(&self) -> [[f64; 2]; 2] {
        self.cd
    }

    #[inline]
    pub fn pixel_scale(&self) -> f64 {
        libm::sqrt(self.determinant.abs())
    }
}

fn compute_inverse(m: [[f64; 2]; 2], det: f64) -> [[f64; 2]; 2] {
    let inv_det = 1.0 / det;
    [
        [m[1][1] * inv_det, -m[0][1] * inv_det],
        [-m[1][0] * inv_det, m[0][0] * inv_det],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_pixel_intermediate_pixel() {
        let crpix = [512.0, 512.0];
        let cd = [[0.001, 0.0], [0.0, 0.001]];
        let transform = LinearTransform::from_cd(crpix, cd).unwrap();

        let original = PixelCoord::new(256.0, 768.0);
        let intermediate = transform.pixel_to_intermediate(original);
        let recovered = transform.intermediate_to_pixel(intermediate);

        assert_eq!(original.x(), recovered.x());
        assert_eq!(original.y(), recovered.y());
    }

    #[test]
    fn test_known_values() {
        let crpix = [512.0, 512.0];
        let cd = [[0.001, 0.0], [0.0, 0.001]];
        let transform = LinearTransform::from_cd(crpix, cd).unwrap();

        let pixel = PixelCoord::new(256.0, 256.0);
        let inter = transform.pixel_to_intermediate(pixel);

        assert_eq!(inter.x_deg(), -0.256);
        assert_eq!(inter.y_deg(), -0.256);
    }

    #[test]
    fn test_pc_cdelt_equivalence() {
        let crpix = [100.0, 100.0];
        let cd = [[0.002, 0.001], [-0.001, 0.002]];
        let transform_cd = LinearTransform::from_cd(crpix, cd).unwrap();

        let cdelt = [0.002, 0.002];
        let pc = [[1.0, 0.5], [-0.5, 1.0]];
        let transform_pc = LinearTransform::from_pc_cdelt(crpix, pc, cdelt).unwrap();

        assert_eq!(transform_cd.cd_matrix(), transform_pc.cd_matrix());

        let pixel = PixelCoord::new(150.0, 175.0);
        let inter_cd = transform_cd.pixel_to_intermediate(pixel);
        let inter_pc = transform_pc.pixel_to_intermediate(pixel);

        assert_eq!(inter_cd.x_deg(), inter_pc.x_deg());
        assert_eq!(inter_cd.y_deg(), inter_pc.y_deg());
    }

    #[test]
    fn test_non_invertible_matrix() {
        let crpix = [512.0, 512.0];
        let cd = [[1.0, 2.0], [2.0, 4.0]];
        let result = LinearTransform::from_cd(crpix, cd);

        assert!(result.is_err());
        match result {
            Err(WcsError::NonInvertibleMatrix { determinant }) => {
                assert_eq!(determinant, 0.0);
            }
            _ => panic!("Expected NonInvertibleMatrix error"),
        }
    }

    #[test]
    fn test_pixel_scale() {
        let crpix = [512.0, 512.0];
        let cd = [[0.001, 0.0], [0.0, 0.001]];
        let transform = LinearTransform::from_cd(crpix, cd).unwrap();

        assert_eq!(transform.pixel_scale(), 0.001);
    }

    #[test]
    fn test_rotated_matrix_roundtrip() {
        let crpix = [256.0, 256.0];
        let angle = celestial_core::constants::PI / 6.0;
        let scale = 0.0005;
        let (angle_s, angle_c) = angle.sin_cos();
        let cd = [
            [scale * angle_c, -scale * angle_s],
            [scale * angle_s, scale * angle_c],
        ];
        let transform = LinearTransform::from_cd(crpix, cd).unwrap();

        let original = PixelCoord::new(100.0, 400.0);
        let intermediate = transform.pixel_to_intermediate(original);
        let recovered = transform.intermediate_to_pixel(intermediate);

        assert_eq!(original.x(), recovered.x());
        assert_eq!(original.y(), recovered.y());
    }

    #[test]
    fn test_crpix_accessor() {
        let crpix = [123.456, 789.012];
        let cd = [[0.001, 0.0], [0.0, 0.001]];
        let transform = LinearTransform::from_cd(crpix, cd).unwrap();

        assert_eq!(transform.crpix(), crpix);
    }
}
