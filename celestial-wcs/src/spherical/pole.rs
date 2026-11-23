use celestial_core::constants::RAD_TO_DEG;

use crate::common::ZERO_LIMIT;
use crate::error::{WcsError, WcsResult};
use crate::trig::{acosd, atan2d, sincosd, sind};

// A δp solved from eq. 8 for a pole that is exactly ±90° lands within
// rounding of it; 1e-14 rad covers that.
const POLE_TOL: f64 = 1e-14 * RAD_TO_DEG;

// Eq. 8, with φp − φ0 as phi_offset.
pub(super) fn delta_p(delta_0: f64, theta_0: f64, phi_offset: f64, latpole: f64) -> WcsResult<f64> {
    let (sin_offset, cos_offset) = sincosd(phi_offset);
    if sin_offset != 0.0 {
        return oblique_delta_p(delta_0, theta_0, (sin_offset, cos_offset), latpole);
    }
    // The default LONPOLE puts the pole on the reference meridian, where
    // eq. 8 is exactly arg(±cos θ0, sin θ0) ± (90° − δ0). Taking those
    // terms in degrees keeps the fiducial point on CRVAL.
    let base = if cos_offset > 0.0 {
        theta_0
    } else {
        libm::copysign(180.0, theta_0) - theta_0
    };
    let half_width = 90.0 - delta_0;
    closest_valid_root(base + half_width, base - half_width, latpole)
}

fn oblique_delta_p(
    delta_0: f64,
    theta_0: f64,
    (sin_offset, cos_offset): (f64, f64),
    latpole: f64,
) -> WcsResult<f64> {
    let (sin_theta_0, cos_theta_0) = sincosd(theta_0);
    let cos_theta_0_sin_offset = cos_theta_0 * sin_offset;
    let denom_sq = 1.0 - cos_theta_0_sin_offset * cos_theta_0_sin_offset;
    let Some(cosine) = acos_argument(sind(delta_0), denom_sq)? else {
        return Ok(latpole);
    };
    let acos_term = acosd(cosine);
    let base = atan2d(sin_theta_0, cos_theta_0 * cos_offset);
    closest_valid_root(base + acos_term, base - acos_term, latpole)
}

// None when cos θ0 sin(φp − φ0) = ±1 and δ0 = 0, which every δp satisfies.
fn acos_argument(sin_delta_0: f64, denom_sq: f64) -> WcsResult<Option<f64>> {
    if denom_sq.abs() < ZERO_LIMIT {
        if sin_delta_0.abs() < ZERO_LIMIT {
            return Ok(None);
        }
        return Err(WcsError::invalid_parameter(
            "Invalid combination of θ₀, δ₀, and φₚ - no solution for δₚ",
        ));
    }
    let arg = sin_delta_0 / libm::sqrt(denom_sq);
    if arg.abs() > 1.0 + 1e-15 {
        return Err(WcsError::invalid_parameter(
            "Invalid combination of θ₀, δ₀, and φₚ - acos argument out of range",
        ));
    }
    Ok(Some(arg.clamp(-1.0, 1.0)))
}

// The roots of eq. 8 are angles, so −285° is the valid root 75°. The
// remainder by 360° is exact.
fn closest_valid_root(root_1: f64, root_2: f64, latpole: f64) -> WcsResult<f64> {
    let (root_1, root_2) = (
        libm::remainder(root_1, 360.0),
        libm::remainder(root_2, 360.0),
    );
    let valid = |v: f64| (-90.0 - POLE_TOL..=90.0 + POLE_TOL).contains(&v);
    let chosen = match (valid(root_1), valid(root_2)) {
        (true, false) => root_1,
        (false, true) => root_2,
        (true, true) if (root_1 - latpole).abs() <= (root_2 - latpole).abs() => root_1,
        (true, true) => root_2,
        (false, false) => {
            return Err(WcsError::invalid_parameter(
                "No valid solution for δₚ in range [-90°, 90°]",
            ))
        }
    };
    if near_pole(chosen) {
        return Ok(libm::copysign(90.0, chosen));
    }
    Ok(chosen.clamp(-90.0, 90.0))
}

// Paper II eqs. 9 and 10, with the notes that follow them for a pole at
// ±90°. All of them take φp only as φp − φ0, here phi_offset.
pub(super) fn alpha_p(
    alpha_0: f64,
    delta_0: f64,
    theta_0: f64,
    phi_offset: f64,
    delta_p: f64,
) -> f64 {
    if near_pole(delta_0) {
        return alpha_0;
    }
    if delta_p == 90.0 {
        return alpha_0 + phi_offset - 180.0;
    }
    if delta_p == -90.0 {
        return alpha_0 - phi_offset;
    }
    let (sin_term, cos_term) = offset_from_pole(delta_0, theta_0, phi_offset, delta_p);
    alpha_0 - atan2d(sin_term, cos_term)
}

// The sine and cosine of α0 − αp from eqs. 9 and 10, both scaled by
// cos δ0 cos δp, which is positive off the poles.
pub(super) fn offset_from_pole(
    delta_0: f64,
    theta_0: f64,
    phi_offset: f64,
    delta_p: f64,
) -> (f64, f64) {
    let (sin_theta_0, cos_theta_0) = sincosd(theta_0);
    let (sin_delta_p, cos_delta_p) = sincosd(delta_p);
    let sin_term = sind(phi_offset) * cos_theta_0 * cos_delta_p;
    let cos_term = sin_theta_0 - sin_delta_p * sind(delta_0);
    (sin_term, cos_term)
}

pub(super) fn near_pole(lat: f64) -> bool {
    (lat.abs() - 90.0).abs() < POLE_TOL
}

#[cfg(test)]
mod tests {
    use super::*;

    // With θ0 = 0, δ0 = 0 and φp − φ0 = 90°, eq. 8 holds for every δp, and
    // LATPOLE is the answer.
    #[test]
    fn test_every_delta_p_solves_eq_8_on_the_equator() {
        assert_eq!(delta_p(0.0, 0.0, 90.0, 37.0).unwrap(), 37.0);
    }

    // (δ0, θ0, φp − φ0): cos θ0 sin(φp − φ0) = 1 with δ0 ≠ 0; an acos
    // argument of sin 45° / cos 60° = √2; roots −130° and 110°, both past a
    // pole.
    #[test]
    fn test_delta_p_without_a_solution_is_refused() {
        for (delta_0, theta_0, offset) in
            [(10.0, 0.0, 90.0), (45.0, 0.0, 60.0), (30.0, 10.0, 180.0)]
        {
            let result = delta_p(delta_0, theta_0, offset, 90.0);
            assert!(
                matches!(result, Err(WcsError::InvalidParameter { .. })),
                "({delta_0}, {theta_0}, {offset}): {result:?}"
            );
        }
    }

    // On the equator with the default LONPOLE the roots are ±90°, and
    // LATPOLE chooses between them.
    #[test]
    fn test_latpole_chooses_either_root() {
        assert_eq!(delta_p(0.0, 0.0, 0.0, 90.0).unwrap(), 90.0);
        assert_eq!(delta_p(0.0, 0.0, 0.0, -90.0).unwrap(), -90.0);
    }

    // With θ0 = 30° and δ0 = 0 the roots are 120° and −60°, and only the
    // second is a latitude.
    #[test]
    fn test_delta_p_takes_the_one_valid_root() {
        assert_eq!(delta_p(0.0, 30.0, 0.0, 90.0).unwrap(), -60.0);
    }

    // A reference point on the celestial pole leaves αp at α0 (the note to
    // eqs. 9 and 10).
    #[test]
    fn test_alpha_p_of_a_reference_point_on_the_pole() {
        assert_eq!(alpha_p(123.0, 90.0, 0.0, 180.0, 90.0), 123.0);
        assert_eq!(alpha_p(-45.0, -90.0, 0.0, 180.0, -90.0), -45.0);
    }
}
