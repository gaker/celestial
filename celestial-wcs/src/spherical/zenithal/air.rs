use crate::common::{
    intermediate_to_polar, native_coord, pole_native_coord, radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::newton::bracketed_newton;
use crate::trig::{asind, sincosd, to_degrees, to_radians};

pub(crate) fn project_air(native: NativeCoord, k: &WcsResult<f64>) -> WcsResult<IntermediateCoord> {
    let theta = native.theta_deg();
    if theta == 90.0 {
        return Ok(IntermediateCoord::new(0.0, 0.0));
    }
    if theta == -90.0 {
        return Err(WcsError::singularity(
            "AIR projection diverges at theta = -90",
        ));
    }
    let r_theta = compute_air_r_theta(theta, k)?;
    Ok(radial_to_intermediate(r_theta, native.phi_deg()))
}

// Paper II eq. 71.
fn compute_air_r_theta(theta: f64, k: &WcsResult<f64>) -> WcsResult<f64> {
    let h = half_colatitude(theta)
        .ok_or_else(|| WcsError::singularity("AIR projection: cos(xi) <= 0"))?;
    Ok(-2.0 * to_degrees(h.ln_cos / h.tan + k.clone()? * h.tan))
}

// Eq. 71's ln(cos ξb) / tan² ξb, which the paper notes tends to −1/2 as θb
// approaches 90°. That limit stands for it at θb = 90°, where ξb = 0.
pub(crate) fn air_k(theta_b: f64) -> WcsResult<f64> {
    if theta_b == 90.0 {
        return Ok(-0.5);
    }
    match half_colatitude(theta_b) {
        Some(h) => Ok(h.ln_cos / (h.tan * h.tan)),
        None => Err(WcsError::singularity(
            "AIR projection: invalid theta_b parameter",
        )),
    }
}

// The functions of ξ = (90° − θ)/2 that eq. 71 and its derivative use. Where
// cos ξ is near 1, ln cos ξ is ½ log1p(−sin² ξ), which keeps the digits that
// the log of the rounded cosine loses.
struct HalfColatitude {
    ln_cos: f64,
    tan: f64,
    sin2: f64,
    cos2: f64,
}

fn half_colatitude(theta: f64) -> Option<HalfColatitude> {
    let (sin_xi, cos_xi) = sincosd((90.0 - theta) / 2.0);
    (cos_xi > 0.0).then(|| {
        let (sin2, cos2) = (sin_xi * sin_xi, cos_xi * cos_xi);
        let ln_cos = if sin2 < cos2 {
            libm::log1p(-sin2) / 2.0
        } else {
            libm::log(cos_xi)
        };
        HalfColatitude {
            ln_cos,
            tan: sin_xi / cos_xi,
            sin2,
            cos2,
        }
    })
}

pub(crate) fn deproject_air(
    inter: IntermediateCoord,
    k: &WcsResult<f64>,
) -> WcsResult<NativeCoord> {
    let (phi, r_theta, is_pole) = intermediate_to_polar(inter.x_deg(), inter.y_deg());
    if is_pole {
        return Ok(pole_native_coord());
    }
    native_coord(phi, solve_air_inverse(r_theta, k)?)
}

// Eq. 71 solved for θ. R runs from 0 at the north pole to infinity at the
// south, so (−90°, 90°] always holds the root, and r − R rises through it.
// The ZEA radius (eq. 70) seeds the iteration where it is defined.
fn solve_air_inverse(r: f64, k: &WcsResult<f64>) -> WcsResult<f64> {
    let rho = to_radians(r) / 2.0;
    let seed = if rho < 1.0 {
        90.0 - 2.0 * asind(rho)
    } else {
        0.0
    };
    let f_and_slope = |theta| {
        let r_theta = compute_air_r_theta(theta, k)?;
        Ok((r - r_theta, -compute_air_dr_dtheta(theta, k)?))
    };
    bracketed_newton((-90.0, 90.0), seed, f_and_slope, "AIR inverse")
}

// dRθ/dθ in degrees per degree, the same number as in radians per radian.
fn compute_air_dr_dtheta(theta: f64, k: &WcsResult<f64>) -> WcsResult<f64> {
    let h = half_colatitude(theta)
        .ok_or_else(|| WcsError::singularity("AIR derivative: cos(xi) <= 0"))?;
    Ok(-1.0 - h.ln_cos / h.sin2 + k.clone()? / h.cos2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::assert_pinned;
    use crate::spherical::Projection;

    // Against a 60-digit evaluation the projected values are within 3.53 ulp
    // and the deprojected φ within 0.44 ulp. The inverse finds θ to the ~2 ulp
    // to which eq. 71 evaluates R, so θ is 16.17 ulp low at θ = 10° for
    // θb = 30°: R ≈ 78 there is three binades above θ.
    #[test]
    fn test_air_pinned_both_ways() {
        let cases = [
            (
                90.0,
                (
                    (-180.0, 10.0),
                    (0.0, 84.47347480700826),
                    (180.0, 10.000000000000012),
                ),
            ),
            (
                90.0,
                (
                    (135.0, 30.0),
                    (43.57833528667381, 43.57833528667382),
                    (135.0, 29.999999999999993),
                ),
            ),
            (
                45.0,
                (
                    (30.0, 89.0),
                    (0.4807314393582714, -0.8326512777642426),
                    (30.000000000000004, 89.0),
                ),
            ),
            (
                45.0,
                ((-90.0, 60.0), (-28.99519052100856, -0.0), (-90.0, 60.0)),
            ),
            (
                30.0,
                (
                    (30.0, 10.0),
                    (38.944582975502854, -67.45399639315286),
                    (30.000000000000004, 10.000000000000009),
                ),
            ),
        ];
        for (theta_b, row) in cases {
            assert_pinned(&Projection::air(theta_b), &[row]);
        }
    }

    #[test]
    fn test_air_near_pole_keeps_the_radius() {
        // 1e-6° from the pole, cos ξ is within an ulp of 1 and its log rounds
        // to zero. Against a 60-digit evaluation, θb = 90 is correctly rounded;
        // the θb = 45 x is 0.79 ulp low and y 1.62 ulp high, from the rounding
        // of ln cos ξb / tan² ξb.
        let native = NativeCoord::new(37.5, 89.999999);
        let cases = [
            (90.0, (6.087614274717458e-7, -7.933533382882105e-7)),
            (45.0, (5.85298246638504e-7, -7.627755256987684e-7)),
        ];
        for (theta_b, expected) in cases {
            let inter = Projection::air(theta_b).project(native).unwrap();
            assert_eq!((inter.x_deg(), inter.y_deg()), expected, "θb = {theta_b}");
        }
    }

    #[test]
    fn test_air_near_origin_keeps_the_latitude() {
        // Against a 60-digit evaluation, φ is 0.81 ulp low and both θ are
        // correctly rounded.
        let inter = IntermediateCoord::new(1e-6, -2e-6);
        let cases = [(90.0, 89.99999776393203), (45.0, 89.99999767429351)];
        for (theta_b, expected) in cases {
            let native = Projection::air(theta_b).deproject(inter).unwrap();
            let got = (native.phi_deg(), native.theta_deg());
            assert_eq!(got, (26.565051177077986, expected), "θb = {theta_b}");
        }
    }

    #[test]
    fn test_air_inverse_reaches_the_far_south() {
        // These latitudes lie south of −89.43°, where an iteration that stalls
        // against a bound must not report success. This far south 90° − θ
        // keeps one bit fewer than θ, so against a 60-digit evaluation the
        // θb = 90 θ is 1.28 ulp low and the θb = 45 θ 0.76 ulp low.
        let inter = IntermediateCoord::new(20000.0, 0.0);
        let cases = [(90.0, -89.67168870721409), (45.0, -89.69699833037569)];
        for (theta_b, expected) in cases {
            let native = Projection::air(theta_b).deproject(inter).unwrap();
            let got = (native.phi_deg(), native.theta_deg());
            assert_eq!(got, (90.0, expected), "θb = {theta_b}");
        }
    }

    #[test]
    fn test_air_edge_cases() {
        // South pole: AIR diverges at theta = -90.
        let proj = Projection::air(90.0);
        let native = NativeCoord::new(0.0, -90.0);
        assert!(proj.project(native).is_err());

        // Invalid theta_b (outside valid range) errors.
        let native = NativeCoord::new(0.0, 45.0);
        let err = Projection::air(-100.0).project(native).unwrap_err();
        assert!(
            matches!(err, WcsError::Singularity { .. }),
            "invalid theta_b: {:?}",
            err
        );
    }
}
