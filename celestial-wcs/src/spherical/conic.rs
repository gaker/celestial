use crate::common::{check_nonzero_param, native_coord, ZERO_LIMIT};
use crate::coordinate::{IntermediateCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{asind, atan2d, atand, cotd, sincosd, sind, tand, to_degrees, to_radians};

// Paper II eqs. 119 and 120.
fn standard_parallels(theta_a: f64, eta: f64) -> (f64, f64) {
    (theta_a - eta, theta_a + eta)
}

// Eqs. 114 and 115, with Rθ, Y0 and Cφ in degrees.
fn conic_to_intermediate(r_theta: f64, y0: f64, c: f64, phi: f64) -> IntermediateCoord {
    let (sin_c_phi, cos_c_phi) = sincosd(c * phi);
    IntermediateCoord::new(r_theta * sin_c_phi, -r_theta * cos_c_phi + y0)
}

// Eqs. 116 and 117, returning (φ, Rθ). Rθ takes the sign of θa, and dividing
// both arguments of arg() by it only flips their signs with θa's. φ has no
// value at the apex; 0 stands for it.
fn intermediate_to_conic(inter: IntermediateCoord, y0: f64, theta_a: f64, c: f64) -> (f64, f64) {
    let (x, y_offset) = (inter.x_deg(), y0 - inter.y_deg());
    let r_theta = libm::copysign(libm::sqrt(x * x + y_offset * y_offset), theta_a);
    if r_theta == 0.0 {
        return (0.0, r_theta);
    }
    let sign = libm::copysign(1.0, theta_a);
    (atan2d(sign * x, sign * y_offset) / c, r_theta)
}

// C and Y0 from eqs. 121 and 123, cos η for eqs. 122 and 124, and cot θa for
// eq. 124. Y0 is kept in radians too, as eq. 122's factor.
#[derive(Debug, Clone, Copy)]
pub(super) struct CopCone {
    c: f64,
    y0: f64,
    y0_rad: f64,
    cos_eta: f64,
    cos_theta_a: f64,
    cot_theta_a: f64,
}

pub(super) fn cop_cone(theta_a: f64, eta: f64) -> WcsResult<CopCone> {
    check_nonzero_param(theta_a, "COP projection: theta_a")?;
    let cos_eta = sincosd(eta).1;
    if cos_eta == 0.0 {
        return Err(WcsError::invalid_parameter("COP projection: eta = +/-90"));
    }
    let (sin_theta_a, cos_theta_a) = sincosd(theta_a);
    let cot_theta_a = cotd(theta_a);
    let y0_rad = cos_eta * cot_theta_a;
    Ok(CopCone {
        c: sin_theta_a,
        y0: to_degrees(y0_rad),
        y0_rad,
        cos_eta,
        cos_theta_a,
        cot_theta_a,
    })
}

// Eqs. 121 to 123.
pub(crate) fn project_cop(
    native: NativeCoord,
    theta_a: f64,
    cone: &WcsResult<CopCone>,
) -> WcsResult<IntermediateCoord> {
    let cone = cone.clone()?;
    let theta = native.theta_deg();
    let delta = cop_offset(theta, theta_a)?;
    let r_theta = cop_radius(&cone, theta, delta);
    Ok(conic_to_intermediate(
        r_theta,
        cone.y0,
        cone.c,
        native.phi_deg(),
    ))
}

// Eq. 122 with cot θa − tan(θ − θa) = cot θa cos θ / (cos θa cos(θ − θa)). The
// printed difference cancels near the pole. The ratio is exactly 1 at θ = θa,
// so Rθ there is Y0. At θa = ±90° cot θa is 0 and the ratio is 0/0, so the
// identity is written as cos θ / (sin θa cos(θ − θa)) instead.
fn cop_radius(cone: &CopCone, theta: f64, delta: f64) -> f64 {
    let cos_theta = sincosd(theta).1;
    let cos_delta = sincosd(delta).1;
    if cone.cos_theta_a == 0.0 {
        return to_degrees(cone.cos_eta * cos_theta / (cone.c * cos_delta));
    }
    to_degrees(cone.y0_rad * (cos_theta / (cone.cos_theta_a * cos_delta)))
}

// Rθ diverges at θ = θa ± 90° (Fig. 24). Past that, tan(θ − θa) comes back
// finite and would put the point on the wrong side of the apex.
fn cop_offset(theta: f64, theta_a: f64) -> WcsResult<f64> {
    let delta = theta - theta_a;
    if delta.abs() == 90.0 {
        return Err(WcsError::singularity(
            "COP projection: theta - theta_a = +/-90",
        ));
    }
    if delta.abs() > 90.0 {
        return Err(WcsError::out_of_bounds(
            "COP projection: point beyond theta_a +/- 90",
        ));
    }
    Ok(delta)
}

// Eq. 124.
pub(crate) fn deproject_cop(
    inter: IntermediateCoord,
    theta_a: f64,
    cone: &WcsResult<CopCone>,
) -> WcsResult<NativeCoord> {
    let cone = cone.clone()?;
    let (phi, r_theta) = intermediate_to_conic(inter, cone.y0, theta_a, cone.c);
    let theta = theta_a + atand(cone.cot_theta_a - to_radians(r_theta) / cone.cos_eta);
    native_coord(phi, theta)
}

// Eqs. 125 and 126.
pub(crate) fn project_coe(
    native: NativeCoord,
    cone: &WcsResult<CoeCone>,
) -> WcsResult<IntermediateCoord> {
    let CoeCone {
        gamma,
        sin_product,
        y0,
    } = cone.clone()?;
    let radicand = 1.0 + sin_product - gamma * sind(native.theta_deg());
    if radicand < 0.0 {
        return Err(WcsError::out_of_bounds(
            "COE projection: point outside valid region",
        ));
    }
    let r_theta = to_degrees(2.0 / gamma * libm::sqrt(radicand));
    Ok(conic_to_intermediate(
        r_theta,
        y0,
        gamma / 2.0,
        native.phi_deg(),
    ))
}

// Eq. 129.
pub(crate) fn deproject_coe(
    inter: IntermediateCoord,
    theta_a: f64,
    cone: &WcsResult<CoeCone>,
) -> WcsResult<NativeCoord> {
    let CoeCone {
        gamma,
        sin_product,
        y0,
    } = cone.clone()?;
    let (phi, r_theta) = intermediate_to_conic(inter, y0, theta_a, gamma / 2.0);
    let half_r = to_radians(r_theta) / 2.0;
    let sin_theta = 1.0 / gamma + sin_product / gamma - gamma * (half_r * half_r);
    if sin_theta.abs() > 1.0 + 1e-12 {
        return Err(WcsError::out_of_bounds(
            "COE deprojection: point outside valid region",
        ));
    }
    native_coord(phi, asind(sin_theta.clamp(-1.0, 1.0)))
}

// γ from eq. 128, the product sin θ1 sin θ2 that eqs. 126, 127 and 129
// share, and Y0.
#[derive(Debug, Clone, Copy)]
pub(super) struct CoeCone {
    gamma: f64,
    sin_product: f64,
    y0: f64,
}

pub(super) fn coe_cone(theta_a: f64, eta: f64) -> WcsResult<CoeCone> {
    check_nonzero_param(theta_a, "COE projection: theta_a")?;
    let (theta_1, theta_2) = standard_parallels(theta_a, eta);
    let (sin_t1, sin_t2) = (sind(theta_1), sind(theta_2));
    let gamma = sin_t1 + sin_t2;
    if gamma.abs() < ZERO_LIMIT {
        return Err(WcsError::invalid_parameter(
            "COE projection: sin(theta_1) + sin(theta_2) = 0",
        ));
    }
    let sin_product = sin_t1 * sin_t2;
    let y0 = coe_y0(gamma, sin_product, theta_a);
    Ok(CoeCone {
        gamma,
        sin_product,
        y0,
    })
}

// Eq. 127, with (θ1 + θ2)/2 = θa (eq. 111). The radicand is
// (cos η − sin²θa)² + (sin θa cos θa)², so only rounding can take it below 0.
fn coe_y0(gamma: f64, sin_product: f64, theta_a: f64) -> f64 {
    let radicand = 1.0 + sin_product - gamma * sind(theta_a);
    to_degrees(2.0 / gamma * libm::sqrt(radicand.max(0.0)))
}

// Eq. 131, or eq. 135 for one standard parallel.
pub(crate) fn project_cod(
    native: NativeCoord,
    theta_a: f64,
    cone: &WcsResult<CodCone>,
) -> WcsResult<IntermediateCoord> {
    let CodCone { c, y0 } = cone.clone()?;
    let r_theta = theta_a - native.theta_deg() + y0;
    Ok(conic_to_intermediate(r_theta, y0, c, native.phi_deg()))
}

// Eq. 133, or eq. 137 for one standard parallel.
pub(crate) fn deproject_cod(
    inter: IntermediateCoord,
    theta_a: f64,
    cone: &WcsResult<CodCone>,
) -> WcsResult<NativeCoord> {
    let CodCone { c, y0 } = cone.clone()?;
    let (phi, r_theta) = intermediate_to_conic(inter, y0, theta_a, c);
    native_coord(phi, theta_a + (y0 - r_theta))
}

// C and Y0 from eqs. 130 and 132. With one standard parallel (η = 0),
// η cot η and sin η / η take their limits, which gives eqs. 134 and 136.
#[derive(Debug, Clone, Copy)]
pub(super) struct CodCone {
    c: f64,
    y0: f64,
}

pub(super) fn cod_cone(theta_a: f64, eta: f64) -> WcsResult<CodCone> {
    check_nonzero_param(theta_a, "COD projection: theta_a")?;
    let (sin_theta_a, cot_theta_a) = (sind(theta_a), cotd(theta_a));
    if eta == 0.0 {
        let y0 = to_degrees(cot_theta_a);
        return Ok(CodCone { c: sin_theta_a, y0 });
    }
    let c = to_degrees(sin_theta_a * sind(eta) / eta);
    let y0 = eta * cotd(eta) * cot_theta_a;
    Ok(CodCone { c, y0 })
}

// Eq. 140. Rθ is 0 at the pole on θa's side and diverges at the other.
pub(crate) fn project_coo(
    native: NativeCoord,
    theta_a: f64,
    cone: &WcsResult<CooCone>,
) -> WcsResult<IntermediateCoord> {
    let CooCone { c, psi, y0 } = cone.clone()?;
    let theta = native.theta_deg();
    if theta.abs() == 90.0 && (theta > 0.0) != (theta_a > 0.0) {
        return Err(WcsError::singularity(
            "COO projection: singularity at opposite pole",
        ));
    }
    let r_theta = psi * libm::pow(tand((90.0 - theta) / 2.0), c);
    Ok(conic_to_intermediate(r_theta, y0, c, native.phi_deg()))
}

// Eq. 144 with 90° − 2 tan⁻¹ w = 2 tan⁻¹((1 − w)/(1 + w)). The printed
// difference cancels for small θ. At the apex Rθ = 0, and w is 0 or infinite
// with the sign of C, which atan2 turns into the pole on θa's side.
pub(crate) fn deproject_coo(
    inter: IntermediateCoord,
    theta_a: f64,
    cone: &WcsResult<CooCone>,
) -> WcsResult<NativeCoord> {
    let CooCone { c, psi, y0 } = cone.clone()?;
    let (phi, r_theta) = intermediate_to_conic(inter, y0, theta_a, c);
    let w = libm::pow(r_theta / psi, 1.0 / c);
    native_coord(phi, 2.0 * atan2d(1.0 - w, 1.0 + w))
}

// C, ψ and Y0 from eqs. 139, 142 and 141. Both standard parallels lie
// strictly between the poles, so the tangents of (90° − θi)/2 are positive.
#[derive(Debug, Clone, Copy)]
pub(super) struct CooCone {
    c: f64,
    psi: f64,
    y0: f64,
}

pub(super) fn coo_cone(theta_a: f64, eta: f64) -> WcsResult<CooCone> {
    check_nonzero_param(theta_a, "COO projection: theta_a")?;
    let (theta_1, theta_2) = standard_parallels(theta_a, eta);
    let (cos_t1, cos_t2) = (sincosd(theta_1).1, sincosd(theta_2).1);
    if cos_t1 <= 0.0 || cos_t2 <= 0.0 {
        return Err(WcsError::invalid_parameter(
            "COO projection: standard parallels reach the poles",
        ));
    }
    let tan_t1 = tand((90.0 - theta_1) / 2.0);
    let c = coo_cone_constant(theta_a, eta)?;
    let psi = to_degrees(cos_t1 / (c * libm::pow(tan_t1, c)));
    let y0 = psi * libm::pow(tand((90.0 - theta_a) / 2.0), c);
    Ok(CooCone { c, psi, y0 })
}

// Eq. 139. With θ1,2 = θa ∓ η, cos θ2 / cos θ1 = (1 − t)/(1 + t) for
// t = tan θa tan η, and ln tan((90° − θ)/2) = −atanh(sin θ), so
// C = atanh(tan θa tan η) / atanh(sin η / cos θa). The printed ratio of logs
// is 0/0 as η → 0 and loses C's relative precision as ε/η. At η = 0 C is
// sin θa (§5.4.4). The caller has checked cos θ1,2 > 0, so |θa| < 90°.
fn coo_cone_constant(theta_a: f64, eta: f64) -> WcsResult<f64> {
    let c = if eta == 0.0 {
        sind(theta_a)
    } else {
        let (sin_eta, cos_theta_a) = (sind(eta), sincosd(theta_a).1);
        libm::atanh(tand(theta_a) * tand(eta)) / libm::atanh(sin_eta / cos_theta_a)
    };
    if !c.is_finite() || c.abs() < ZERO_LIMIT {
        return Err(WcsError::invalid_parameter(
            "COO projection: degenerate cone (C is 0 or not finite)",
        ));
    }
    Ok(c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{assert_pinned, PinnedRow};
    use crate::spherical::Projection;

    fn assert_cases(make: fn(f64, f64) -> Projection, cases: &[(f64, f64, &[PinnedRow])]) {
        for &(theta_a, eta, rows) in cases {
            assert_pinned(&make(theta_a, eta), rows);
        }
    }

    // At θa, Rθ goes through the same steps as Y0, so y = Y0 − Rθ is exactly 0.
    #[test]
    fn test_conic_native_reference_and_origin() {
        let theta_a = 45.0;
        let projs = [
            Projection::cop(theta_a, 0.0),
            Projection::coe(theta_a, 0.0),
            Projection::cod(theta_a, 0.0),
            Projection::coo(theta_a, 0.0),
            Projection::cop(theta_a, 15.0),
            Projection::coe(theta_a, 15.0),
            Projection::cod(theta_a, 15.0),
            Projection::coo(theta_a, 15.0),
        ];
        for proj in &projs {
            let (phi0, theta0) = proj.native_reference();
            assert_eq!(phi0, 0.0);
            assert_eq!(theta0, theta_a);
            let inter = proj.project(NativeCoord::new(0.0, theta_a)).unwrap();
            assert_eq!((inter.x_deg(), inter.y_deg()), (0.0, 0.0), "{proj:?}");
        }
    }

    // y = Y0 − Rθ cos Cφ cancels for points near the arc through the reference
    // point, and Y0 − y cancels coming back, so offsets there are measured
    // against Y0 and Rθ rather than the result.

    // Against a 60-digit evaluation x is within 3.52 ulp, at (−90°, 25°) with
    // η = 25°. y at (−60°, 40°) with θa = 60° is −0.15 from 33.08 − 33.23, and
    // its 288 ulp are 1.1 ulp of Y0 and of Rθ; at (−60°, 60°) with θa = 75° it
    // is −0.94, and its 63 ulp are 3.9 ulp of Y0. The deprojected values are
    // within 4.72 ulp, φ at (30°, 85°), where Y0 − y cancels near the apex.
    // θa = 45° and η = 25° is Paper II Fig. 24.
    #[test]
    fn test_cop_pinned_both_ways() {
        let cases: [(f64, f64, &[PinnedRow]); 4] = [
            (
                45.0,
                0.0,
                &[
                    (
                        (-60.0, 20.0),
                        (-56.678915923679604, -4.718130475405566),
                        (-59.999999999999986, 20.000000000000004),
                    ),
                    ((0.0, 40.0), (0.0, -5.012731175275583), (0.0, 40.0)),
                    (
                        (180.0, 60.0),
                        (33.374095467600185, 82.70090443909527),
                        (180.0, 60.000000000000014),
                    ),
                ],
            ),
            (
                60.0,
                0.0,
                &[(
                    (-60.0, 40.0),
                    (-42.478021661903284, -0.15369538291923845),
                    (-60.0, 39.999999999999986),
                )],
            ),
            (
                75.0,
                0.0,
                &[(
                    (-60.0, 60.0),
                    (-26.026444699002543, -0.9388588807857126),
                    (-59.99999999999999, 60.00000000000001),
                )],
            ),
            (
                45.0,
                25.0,
                &[
                    (
                        (-180.0, 45.0),
                        (-41.318447155383694, 83.38015830046777),
                        (180.0, 45.000000000000014),
                    ),
                    (
                        (-90.0, 25.0),
                        (-63.46297465946199, 20.478983331117462),
                        (-89.99999999999999, 25.00000000000001),
                    ),
                    (
                        (30.0, 85.0),
                        (3.023230415139762, 44.138582096289916),
                        (29.999999999999993, 85.0),
                    ),
                    ((37.5, 90.0), (0.0, 51.92761113704146), (0.0, 90.0)),
                ],
            ),
        ];
        assert_cases(Projection::cop, &cases);
    }

    // Near the pole eq. 122's cot θa − tan(θ − θa) cancels. Against a 60-digit
    // evaluation x at (30°, 89.999999°) is 1.90 ulp low and at
    // (30°, 89.9999999999°) 1.69 ulp low.
    #[test]
    fn test_cop_near_pole_keeps_x() {
        let x = |theta| project_at(Projection::cop(45.0, 0.0), 30.0, theta).0;
        assert_eq!(x(89.999999), 7.236788022764762e-7);
        assert_eq!(x(89.9999999999), 7.236917317285143e-11);
    }

    // At θa = 90° cot θa is 0 and C is 1. Against a 60-digit evaluation every
    // value is within 1.90 ulp.
    #[test]
    fn test_cop_at_the_pole_pinned_both_ways() {
        let rows: [PinnedRow; 2] = [
            (
                (30.0, 60.0),
                (16.539866862653763, -28.64788975654116),
                (30.000000000000004, 60.0),
            ),
            (
                (-40.0, 20.0),
                (-99.64963675251275, -118.7578126005151),
                (-40.0, 20.0),
            ),
        ];
        assert_pinned(&Projection::cop(90.0, 0.0), &rows[..1]);
        assert_pinned(&Projection::cop(90.0, 10.0), &rows[1..]);
    }

    // Against a 60-digit evaluation x is within 4.93 ulp, at θ = 80° where eq.
    // 126's radicand cancels to 0.11, and y within 13.55 ulp, at (0°, 50°), 1.7
    // ulp of Y0. The deprojected values are within 7.35 ulp, θ at (−120°, 20°)
    // where eq. 129's sin θ = 0.34 is 1.06 − 0.72. The pole is an arc (eq. 126
    // with sin θ = 1), and the projection of (37.5°, 90°) rounds 1.6e-17 past it
    // in sin θ, which comes back as the pole.
    #[test]
    fn test_coe_pinned_both_ways() {
        let cases: [(f64, f64, &[PinnedRow]); 3] = [
            (
                45.0,
                0.0,
                &[
                    (
                        (-120.0, 20.0),
                        (-81.35720235078874, 49.96730180371959),
                        (-120.00000000000001, 19.999999999999996),
                    ),
                    (
                        (0.0, 50.0),
                        (0.0, 4.993201784217014),
                        (0.0, 50.00000000000001),
                    ),
                    (
                        (90.0, 80.0),
                        (23.77918874440753, 45.512170875723996),
                        (89.99999999999999, 80.00000000000001),
                    ),
                ],
            ),
            (
                45.0,
                10.0,
                &[(
                    (90.0, 20.0),
                    (72.85836293567, 19.653748923228733),
                    (89.99999999999999, 19.999999999999996),
                )],
            ),
            (
                45.0,
                25.0,
                &[
                    (
                        (-120.0, 20.0),
                        (-81.8277960736184, 38.56336501285425),
                        (-120.0, 20.000000000000004),
                    ),
                    (
                        (180.0, 50.0),
                        (47.078686415465384, 79.90957765076973),
                        (179.99999999999997, 50.0),
                    ),
                    (
                        (37.5, 90.0),
                        (7.252923041537599, 41.33533619266996),
                        (37.5, 90.0),
                    ),
                ],
            ),
        ];
        assert_cases(Projection::coe, &cases);
    }

    // Against a 60-digit evaluation x is within 3.19 ulp and y within 7.19 ulp,
    // at (90°, 20°) with η = 20°, 1.8 ulp of Rθ. The deprojected values are
    // within 2.63 ulp.
    // At the reference point Rθ = |Y0| exactly, so eq. 133 gives θa with no
    // rounding.
    #[test]
    fn test_cod_reference_point_is_exact() {
        for theta_a in [-60.0, -30.0, 10.0, 30.0, 45.0, 75.0] {
            for eta in [0.0, 5.0, 20.0] {
                let back =
                    Projection::cod(theta_a, eta).deproject(IntermediateCoord::new(0.0, 0.0));
                assert_eq!(back.unwrap().theta_deg(), theta_a, "({theta_a}, {eta})");
            }
        }
    }

    #[test]
    fn test_cod_pinned_both_ways() {
        let cases: [(f64, f64, &[PinnedRow]); 3] = [
            (
                45.0,
                0.0,
                &[
                    (
                        (-120.0, 20.0),
                        (-81.963923283894, 49.912649717541704),
                        (-120.0, 20.0),
                    ),
                    (
                        (90.0, 50.0),
                        (46.85800871277491, 34.07562502706671),
                        (90.0, 50.00000000000001),
                    ),
                ],
            ),
            (
                45.0,
                5.0,
                &[(
                    (-120.0, 50.0),
                    (-51.93108379398082, 52.37403790713023),
                    (-120.0, 50.0),
                )],
            ),
            (
                45.0,
                20.0,
                &[
                    (
                        (90.0, 20.0),
                        (70.82251377535216, 17.85370185393586),
                        (90.0, 19.999999999999993),
                    ),
                    (
                        (0.0, 20.0),
                        (0.0, -25.000000000000007),
                        (0.0, 19.999999999999993),
                    ),
                    (
                        (180.0, 80.0),
                        (16.39938865090083, 66.30932535897378),
                        (180.0, 80.0),
                    ),
                    (
                        (37.5, 90.0),
                        (4.358673377559067, 46.00552918435849),
                        (37.500000000000014, 90.0),
                    ),
                ],
            ),
        ];
        assert_cases(Projection::cod, &cases);
    }

    // Paper II eq. 140: with θa > 0, C > 0 and Rθ = ψ tan((90° − θ)/2)^C is 0 at
    // the north pole, which projects to (0, Y0).
    #[test]
    fn test_coo_projects_its_convergent_pole() {
        let proj = Projection::coo(45.0, 25.0);
        let pole = NativeCoord::new(37.5, 90.0);
        let inter = proj.project(pole).unwrap();
        assert_eq!(inter.x_deg(), 0.0);
        assert_eq!(proj.deproject(inter).unwrap().theta_deg(), 90.0);
    }

    // Against a 60-digit evaluation x is within 2.08 ulp and y within 188 ulp,
    // at (−72.25°, 10°) with η = 25°, where it is −0.54 and the offset is 1.5
    // ulp of Rθ. The deprojected values are within 7.21 ulp, θ at (−72.25°, 10°),
    // where eq. 144 turns a relative error in w into 5.6 times as much in θ. At
    // (−180°, 22.5°) the projected point is 5e-14° past the seam, and φ comes
    // back as 180°.
    #[test]
    fn test_coo_pinned_both_ways() {
        let cases: [(f64, f64, &[PinnedRow]); 3] = [
            (
                45.0,
                0.0,
                &[
                    (
                        (-90.0, 50.0),
                        (-46.85217987645359, 34.078513465640455),
                        (-90.0, 50.00000000000002),
                    ),
                    (
                        (0.0, 80.0),
                        (0.0, 38.213386083412345),
                        (0.0, 80.00000000000001),
                    ),
                ],
            ),
            (
                45.0,
                10.0,
                &[(
                    (90.0, 25.0),
                    (68.49310113518582, 22.6804811195149),
                    (90.0, 25.0),
                )],
            ),
            (
                45.0,
                25.0,
                &[
                    (
                        (0.0, 80.0),
                        (0.0, 34.00174564100206),
                        (0.0, 80.00000000000001),
                    ),
                    (
                        (-72.25, 10.0),
                        (-66.94763490803183, -0.5380267336768654),
                        (-72.24999999999999, 10.00000000000002),
                    ),
                    (
                        (-180.0, 22.5),
                        (-52.81811549810938, 97.42463424571493),
                        (180.0, 22.500000000000004),
                    ),
                ],
            ),
        ];
        assert_cases(Projection::coo, &cases);
    }

    // Negating θa negates C, Y0 and Rθ (Paper II §5.4), which mirrors the plane in
    // y with (φ, θ) ↦ (φ, −θ). Every step is odd or even in θa, so it is bit-exact.
    #[test]
    fn test_southern_conic_mirrors_northern() {
        let native = |phi: f64, theta: f64| NativeCoord::new(phi, theta);
        for (north, south) in [
            (Projection::cop(45.0, 25.0), Projection::cop(-45.0, 25.0)),
            (Projection::coe(45.0, 25.0), Projection::coe(-45.0, 25.0)),
            (Projection::cod(45.0, 20.0), Projection::cod(-45.0, 20.0)),
        ] {
            let inter = north.project(native(60.0, 30.0)).unwrap();
            let south_inter = south.project(native(60.0, -30.0)).unwrap();
            assert_eq!(
                (south_inter.x_deg(), south_inter.y_deg()),
                (inter.x_deg(), -inter.y_deg()),
                "{south:?} forward"
            );
            let n = north.deproject(inter).unwrap();
            let s = south.deproject(south_inter).unwrap();
            assert_eq!(
                (s.phi_deg(), s.theta_deg()),
                (n.phi_deg(), -n.theta_deg()),
                "{south:?} inverse"
            );
        }
    }

    // Eq. 144's printed 90° − 2 tan⁻¹ w cancels where θ is small: at
    // (0, 180°/π) with θa = −45°, η = 10°, 2 tan⁻¹ w is 83.8° and one ulp of it
    // is 16 ulp of θ. Against a 60-digit evaluation θ is 11.6 ulp high; a
    // relative error in w moves θ 9.2 times as much there.
    #[test]
    fn test_coo_inverse_keeps_small_latitudes() {
        let coo = Projection::coo(-45.0, 10.0);
        let inter = IntermediateCoord::new(0.0, 57.29577951308232);
        assert_eq!(coo.deproject(inter).unwrap().theta_deg(), 6.195262831128368);
    }

    // COO's C goes through tan((90° − θ)/2), which is not odd in θ, so the mirror is
    // not bit-exact for it. Paper II §5.4 evaluated to 60 digits rounds φ to
    // −120.00000000000004, and this is 1.20 ulp above that.
    #[test]
    fn test_southern_coo_inverse_longitude() {
        let coo = Projection::coo(-45.0, 25.0);
        let inter = IntermediateCoord::new(-45.449742896703405, -48.37048269969766);
        let native = coo.deproject(inter).unwrap();
        assert_eq!(native.phi_deg(), -120.00000000000003);
    }

    #[test]
    fn test_cop_rejects_points_beyond_its_divergence() {
        // COP diverges at θ = θa ± 90° (Paper II Fig. 24).
        let cases = [
            (45.0, -44.0, true),
            (45.0, -46.0, false),
            (45.0, -90.0, false),
            (-45.0, 44.0, true),
            (-45.0, 46.0, false),
            (-45.0, 90.0, false),
        ];
        for (theta_a, theta, inside) in cases {
            let native = NativeCoord::new(37.5, theta);
            let got = Projection::cop(theta_a, 25.0).project(native);
            let label = format!("θa = {theta_a} at θ = {theta}: {got:?}");
            assert_eq!(got.is_ok(), inside, "{label}");
            assert!(
                inside || matches!(got, Err(WcsError::OutOfBounds { .. })),
                "{label}"
            );
        }
    }

    #[test]
    fn test_conic_error_cases() {
        // theta_a = 0 is invalid for every conic (the cone degenerates to a plane).
        for proj in [
            Projection::cop(0.0, 0.0),
            Projection::coe(0.0, 0.0),
            Projection::cod(0.0, 0.0),
            Projection::coo(0.0, 0.0),
        ] {
            let native = NativeCoord::new(0.0, 45.0);
            assert!(
                proj.project(native).is_err(),
                "{:?} theta_a=0 should error",
                proj
            );
        }

        // COP: theta - theta_a = +/-90 is singular (cos(delta) = 0).  With
        // theta_a = 45, theta = -45 hits the singularity.
        let cop = Projection::cop(45.0, 0.0);
        let native = NativeCoord::new(0.0, -45.0);
        assert!(cop.project(native).is_err());

        // COE maps the sphere inside Rθ = (180°/π)(2/γ)√(1 + sin θ1 sin θ2 +
        // γ), 142.0° for (45°, 20°), and y = −500 lies past it.
        let far = Projection::coe(45.0, 20.0).deproject(IntermediateCoord::new(0.0, -500.0));
        assert!(matches!(far, Err(WcsError::OutOfBounds { .. })), "{far:?}");

        // COO is singular at the pole opposite theta_a (tan((90-theta)/2) <= 0).
        let coo = Projection::coo(45.0, 0.0);
        let native = NativeCoord::new(0.0, -90.0);
        assert!(coo.project(native).is_err());
    }

    fn project_at(proj: Projection, phi: f64, theta: f64) -> (f64, f64) {
        let inter = proj.project(NativeCoord::new(phi, theta)).unwrap();
        (inter.x_deg(), inter.y_deg())
    }

    // As η → 0 each conic meets its one-standard-parallel form. COP and COE
    // meet it bit for bit at η = 1e-8. COD's eqs. 130 and 132 there land an ulp
    // below eqs. 134 and 136 in x; against a 60-digit evaluation they are 1.95
    // and 0.95 ulp low, and y is 2.88 ulp low in both.
    #[test]
    fn test_small_eta_meets_one_standard_parallel() {
        let at = |proj| project_at(proj, 30.0, 55.0);
        assert_eq!(
            at(Projection::cop(45.0, 1e-8)),
            at(Projection::cop(45.0, 0.0))
        );
        assert_eq!(
            at(Projection::coe(45.0, 1e-8)),
            at(Projection::coe(45.0, 0.0))
        );
        let one_standard = (17.113476877273744, 13.204741093454295);
        assert_eq!(at(Projection::cod(45.0, 0.0)), one_standard);
        let small = (17.11347687727374, 13.204741093454295);
        assert_eq!(at(Projection::cod(45.0, 1e-8)), small);
    }

    // The printed eq. 139 is 0/0 as η → 0. Against a 60-digit evaluation at
    // η = 1e-8, x is 1.92 ulp high and y 6.70 ulp high; at η = 0, eq. 135's
    // C = sin θa gives 0.90 and 2.76 ulp.
    #[test]
    fn test_coo_small_eta_meets_one_standard_parallel() {
        let at = |eta| project_at(Projection::coo(45.0, eta), 30.0, 55.0);
        assert_eq!(at(1e-8), (17.09404806188135, 13.254797353328186));
        assert_eq!(at(0.0), (17.09404806188134, 13.254797353328179));
    }
}
