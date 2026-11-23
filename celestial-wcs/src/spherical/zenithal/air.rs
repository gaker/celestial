use crate::common::{
    intermediate_to_polar, native_coord, pole_native_coord, radial_to_intermediate,
};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::newton::bracketed_newton;
use crate::trig::{acosd, asind, sincosd, to_degrees, to_radians};

// What eq. 71 derives from θb: its constant k and the θ where R turns.
#[derive(Debug, Clone)]
pub(crate) struct AirRadius {
    k: f64,
    turns: Option<AirTurns>,
}

// The southern and northern θ where dR/dθ is 0, and R at the northern one.
#[derive(Debug, Clone, Copy)]
struct AirTurns {
    south: f64,
    north: f64,
    r_north: f64,
}

pub(crate) fn air_radius(theta_b: f64) -> WcsResult<AirRadius> {
    let k = air_k(theta_b)?;
    let turns = air_turns(k)?;
    Ok(AirRadius { k, turns })
}

pub(crate) fn project_air(
    native: NativeCoord,
    radius: &WcsResult<AirRadius>,
) -> WcsResult<IntermediateCoord> {
    let k = radius.as_ref().map_err(Clone::clone)?.k;
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

// Paper II eq. 71. At the pole ln(cos ξ)/tan ξ is 0/0, and R's limit is 0.
fn compute_air_r_theta(theta: f64, k: f64) -> WcsResult<f64> {
    if theta == 90.0 {
        return Ok(0.0);
    }
    let h = half_colatitude(theta)
        .ok_or_else(|| WcsError::singularity("AIR projection: cos(xi) <= 0"))?;
    Ok(-2.0 * to_degrees(h.ln_cos / h.tan + k * h.tan))
}

// Eq. 71's ln(cos ξb) / tan² ξb, which the paper notes tends to −1/2 as θb
// approaches 90°. That limit stands for it at θb = 90°, where ξb = 0.
fn air_k(theta_b: f64) -> WcsResult<f64> {
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

// With u = cos² ξ, dR/dθ is (k + q(u))/u for q(u) = −u − u ln u / (2(1 − u)).
// q is 0 at u = 0 and −1/2 at u = 1, and q′ changes sign once, near
// u = 0.0699 (θ ≈ −59.33°), where q is about 0.0301. k + q is negative at both
// ends, so R turns at two θ when k + q is positive at that maximum, which is
// when θb is below about −76.47°, and nowhere otherwise. Small u is the south.
fn air_turns(k: f64) -> WcsResult<Option<AirTurns>> {
    let q = |u: f64| -u - u * libm::log(u) / (2.0 * (1.0 - u));
    let q_slope = |u: f64| {
        let one_minus_u = 1.0 - u;
        -1.0 - (libm::log(u) + one_minus_u) / (2.0 * one_minus_u * one_minus_u)
    };
    let peak = sign_change(|u| -q_slope(u), (0.0, 1.0));
    if k + q(peak) <= 0.0 {
        return Ok(None);
    }
    let theta = |u: f64| 90.0 - 2.0 * acosd(libm::sqrt(u));
    let south = theta(sign_change(|u| k + q(u), (0.0, peak)));
    let north = theta(sign_change(|u| -(k + q(u)), (peak, 1.0)));
    let r_north = compute_air_r_theta(north, k)?;
    Ok(Some(AirTurns {
        south,
        north,
        r_north,
    }))
}

// Bisection for the point in (lo, hi) where f goes from negative to positive.
// It stops when the midpoint rounds to an end, so it ends within the ~1100
// halvings that separate two doubles in (0, 1).
fn sign_change(f: impl Fn(f64) -> f64, (mut lo, mut hi): (f64, f64)) -> f64 {
    loop {
        let mid = lo + (hi - lo) / 2.0;
        if mid <= lo || mid >= hi {
            return mid;
        }
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
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
    radius: &WcsResult<AirRadius>,
) -> WcsResult<NativeCoord> {
    let (phi, r_theta, is_pole) = intermediate_to_polar(inter.x_deg(), inter.y_deg());
    if is_pole {
        return Ok(pole_native_coord());
    }
    let radius = radius.as_ref().map_err(Clone::clone)?;
    native_coord(phi, solve_air_inverse(r_theta, radius)?)
}

// Eq. 71 solved for θ. R runs from 0 at the north pole to infinity at the
// south. Where it turns, a radius up to R at the northern turn belongs to as
// many as three θ, and this takes the one nearest the pole, as the ZPN inverse
// does. A larger radius belongs to one θ, south of the southern turn. On each
// bracket r − R rises through the root. The ZEA radius (eq. 70) seeds the
// iteration where it is defined and inside the bracket.
fn solve_air_inverse(r: f64, radius: &AirRadius) -> WcsResult<f64> {
    let (lo, hi) = match radius.turns {
        Some(turns) if r <= turns.r_north => (turns.north, 90.0),
        Some(turns) => (-90.0, turns.south),
        None => (-90.0, 90.0),
    };
    let rho = to_radians(r) / 2.0;
    let seed = if rho < 1.0 {
        90.0 - 2.0 * asind(rho)
    } else {
        0.0
    };
    let seed = if (lo..=hi).contains(&seed) {
        seed
    } else {
        lo + (hi - lo) / 2.0
    };
    let k = radius.k;
    let f_and_slope = |theta| {
        let r_theta = compute_air_r_theta(theta, k)?;
        Ok((r - r_theta, -compute_air_dr_dtheta(theta, k)?))
    };
    bracketed_newton((lo, hi), seed, f_and_slope, "AIR inverse")
}

// dRθ/dθ in degrees per degree, the same number as in radians per radian.
// At the pole ln(cos ξ)/sin² ξ is 0/0 with limit −1/2, so dR/dθ is k − 1/2.
fn compute_air_dr_dtheta(theta: f64, k: f64) -> WcsResult<f64> {
    if theta == 90.0 {
        return Ok(k - 0.5);
    }
    let h = half_colatitude(theta)
        .ok_or_else(|| WcsError::singularity("AIR derivative: cos(xi) <= 0"))?;
    Ok(-1.0 - h.ln_cos / h.sin2 + k / h.cos2)
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
    // A point next to the reference point deprojects to just below the pole.
    // The θ are the 60-digit values, correctly rounded.
    #[test]
    fn test_air_inverse_next_to_the_pole() {
        let air = Projection::air(45.0);
        for (x, theta) in [
            (1e-15, 90.0),
            (1e-13, 89.9999999999999),
            (1e-9, 89.99999999895991),
        ] {
            let native = air.deproject(IntermediateCoord::new(x, 0.0)).unwrap();
            assert_eq!((native.phi_deg(), native.theta_deg()), (90.0, theta), "{x}");
        }
    }

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

    // Below θb ≈ −76.47° R turns twice. For θb = −80 it rises between
    // θ ≈ −74.69° and −45.20°, where R is 46.97° and 50.76°. Against 60 digits
    // the southern θ is 0.29 ulp high, the northern θ 3.36 ulp high and R there
    // 1.07 ulp high.
    #[test]
    fn test_air_turns_below_the_threshold() {
        assert!(air_turns(air_k(-76.47).unwrap()).unwrap().is_none());
        assert!(air_turns(air_k(-76.48).unwrap()).unwrap().is_some());
        let turns = air_turns(air_k(-80.0).unwrap()).unwrap().unwrap();
        assert_eq!(
            (turns.south, turns.north, turns.r_north),
            (-74.68610976650166, -45.19507584653405, 50.75987255064677)
        );
    }

    // (120.29°, −60.98°) lies between the turns of θb = −80, so its radius also
    // belongs to two θ further north, and the inverse returns the nearest the
    // pole, −28.42°. A radius above R at the northern turn belongs to one θ,
    // south of the turns. Against 60 digits x and y are within 2.06 ulp, φ is
    // 0.80 ulp low and θ = −83° is 0.32 ulp low. θ = −28.42° is 11.53 ulp
    // high: dR/dθ is −0.164 there and R is a binade above θ, so each ulp of R
    // moves θ 12.2 ulp.
    #[test]
    fn test_air_inverse_takes_the_root_nearest_the_pole() {
        let rows = [
            (
                (120.29, -60.98),
                (42.584879596534776, 24.874625037334816),
                (120.28999999999999, -28.42452060022634),
            ),
            ((180.0, -83.0), (0.0, 54.58911310938161), (180.0, -83.0)),
        ];
        assert_pinned(&Projection::air(-80.0), &rows);
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
            matches!(err, WcsError::InvalidParameter { .. }),
            "invalid theta_b: {:?}",
            err
        );
    }
}
