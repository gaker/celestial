mod coefficients;

use crate::{CoordError, CoordResult};
use astro_core::{
    constants::{DAYS_PER_JULIAN_YEAR, J2000_JD, SPEED_OF_LIGHT_AU_PER_DAY},
    Vector3,
};
use astro_time::TT;
use coefficients::Coefficients;

pub struct EarthState {
    pub barycentric_velocity: Vector3,
    pub heliocentric_position: Vector3,
}

pub fn compute_earth_state(tt: &TT) -> CoordResult<EarthState> {
    let jd = tt.to_julian_date();
    let t = ((jd.jd1() - J2000_JD) + jd.jd2()) / DAYS_PER_JULIAN_YEAR;

    if t.abs() > 100.0 {
        return Err(CoordError::invalid_coordinate("Epoch outside 1900-2100"));
    }

    let t2 = t * t;

    let e_coefs = [
        [Coefficients::E0X, Coefficients::E1X, Coefficients::E2X],
        [Coefficients::E0Y, Coefficients::E1Y, Coefficients::E2Y],
        [Coefficients::E0Z, Coefficients::E1Z, Coefficients::E2Z],
    ];
    let s_coefs = [
        [Coefficients::S0X, Coefficients::S1X, Coefficients::S2X],
        [Coefficients::S0Y, Coefficients::S1Y, Coefficients::S2Y],
        [Coefficients::S0Z, Coefficients::S1Z, Coefficients::S2Z],
    ];

    let mut ph = [0.0; 3];
    let mut vh = [0.0; 3];
    let mut vb = [0.0; 3];

    for i in 0..3 {
        let mut xyz = 0.0;
        let mut xyzd = 0.0;

        accumulate_terms(t, t2, e_coefs[i][0], &mut xyz, &mut xyzd, 0);
        accumulate_terms(t, t2, e_coefs[i][1], &mut xyz, &mut xyzd, 1);
        accumulate_terms(t, t2, e_coefs[i][2], &mut xyz, &mut xyzd, 2);

        ph[i] = xyz;
        vh[i] = xyzd / DAYS_PER_JULIAN_YEAR;

        accumulate_terms(t, t2, s_coefs[i][0], &mut xyz, &mut xyzd, 0);
        accumulate_terms(t, t2, s_coefs[i][1], &mut xyz, &mut xyzd, 1);
        accumulate_terms(t, t2, s_coefs[i][2], &mut xyz, &mut xyzd, 2);

        vb[i] = xyzd / DAYS_PER_JULIAN_YEAR;
    }

    let helio_pos_bcrs = rotate_to_bcrs(ph);
    let bary_vel_bcrs = rotate_to_bcrs(vb);

    Ok(EarthState {
        barycentric_velocity: Vector3::new(bary_vel_bcrs[0], bary_vel_bcrs[1], bary_vel_bcrs[2]),
        heliocentric_position: Vector3::new(
            helio_pos_bcrs[0],
            helio_pos_bcrs[1],
            helio_pos_bcrs[2],
        ),
    })
}

fn accumulate_terms(t: f64, t2: f64, coefs: &[[f64; 3]], xyz: &mut f64, xyzd: &mut f64, power: u8) {
    for &[a, b, c] in coefs {
        let ct = c * t;
        let p = b + ct;
        let (sin_p, cos_p) = p.sin_cos();

        match power {
            0 => {
                *xyz += a * cos_p;
                *xyzd -= a * c * sin_p;
            }
            1 => {
                *xyz += a * t * cos_p;
                *xyzd += a * (cos_p - ct * sin_p);
            }
            _ => {
                *xyz += a * t2 * cos_p;
                *xyzd += a * t * (2.0 * cos_p - ct * sin_p);
            }
        }
    }
}

fn rotate_to_bcrs(v: [f64; 3]) -> [f64; 3] {
    const AM12: f64 = 0.000000211284;
    const AM13: f64 = -0.000000091603;
    const AM21: f64 = -0.000000230286;
    const AM22: f64 = 0.917482137087;
    const AM23: f64 = -0.397776982902;
    const AM32: f64 = 0.397776982902;
    const AM33: f64 = 0.917482137087;

    let (x, y, z) = (v[0], v[1], v[2]);
    [
        x + AM12 * y + AM13 * z,
        AM21 * x + AM22 * y + AM23 * z,
        AM32 * y + AM33 * z,
    ]
}

const SCHWARZSCHILD_RADIUS_SUN_AU: f64 = 1.97412574336e-8;

/// Apply gravitational light deflection by the Sun.
///
/// This implements the relativistic bending of starlight as it passes near the Sun,
/// based on the algorithm in IERS Conventions.
///
/// # Arguments
/// * `star_direction` - Unit vector from observer to star (BCRS)
/// * `sun_to_observer` - Unit vector from Sun to observer (BCRS)
/// * `sun_observer_distance_au` - Distance from Sun to observer in AU
///
/// # Returns
/// Deflected star direction (unit vector)
pub fn apply_light_deflection(
    star_direction: Vector3,
    sun_to_observer: Vector3,
    sun_observer_distance_au: f64,
) -> Vector3 {
    // Deflection limiter: for nearby observers, use smaller limit
    let em2 = sun_observer_distance_au * sun_observer_distance_au;
    let em2_clamped = if em2 < 1.0 { 1.0 } else { em2 };
    let dlim = 1e-6 / em2_clamped;

    // For distant stars, the direction from Sun to star ≈ direction from observer to star
    // So we use star_direction for both p and q in ERFA's eraLd
    let q = star_direction;
    let e = sun_to_observer;

    // q + e
    let qpe = Vector3::new(q.x + e.x, q.y + e.y, q.z + e.z);

    // q . (q + e)
    let qdqpe = q.dot(&qpe);

    // Apply limiter to avoid division by zero when star is near Sun
    let qdqpe_limited = if qdqpe > dlim { qdqpe } else { dlim };

    // 2G*M / (c^2 * r * (q.(q+e))) = SRS / (em * (q.(q+e)))
    // where SRS = Schwarzschild radius of Sun in AU
    let w = SCHWARZSCHILD_RADIUS_SUN_AU / sun_observer_distance_au / qdqpe_limited;

    // e × q (cross product)
    let eq = Vector3::new(
        e.y * q.z - e.z * q.y,
        e.z * q.x - e.x * q.z,
        e.x * q.y - e.y * q.x,
    );

    // p × (e × q) (cross product)
    let p = star_direction;
    let peq = Vector3::new(
        p.y * eq.z - p.z * eq.y,
        p.z * eq.x - p.x * eq.z,
        p.x * eq.y - p.y * eq.x,
    );

    // Apply deflection: p1 = p + w * (p × (e × q))
    Vector3::new(p.x + w * peq.x, p.y + w * peq.y, p.z + w * peq.z)
}

/// Remove gravitational light deflection by the Sun (inverse operation).
pub fn remove_light_deflection(
    deflected_direction: Vector3,
    sun_to_observer: Vector3,
    sun_observer_distance_au: f64,
) -> Vector3 {
    // We compute: d = forward(estimate) - estimate, then refine: estimate = observed - d
    let mut d = Vector3::zeros();

    for _ in 0..5 {
        // Estimate the original direction
        let before = Vector3::new(
            deflected_direction.x - d.x,
            deflected_direction.y - d.y,
            deflected_direction.z - d.z,
        )
        .normalize();

        // Apply forward deflection to the estimate
        let after = apply_light_deflection(before, sun_to_observer, sun_observer_distance_au);

        // Update correction
        d = Vector3::new(after.x - before.x, after.y - before.y, after.z - before.z);
    }

    // Final result
    Vector3::new(
        deflected_direction.x - d.x,
        deflected_direction.y - d.y,
        deflected_direction.z - d.z,
    )
    .normalize()
}

/// Remove stellar aberration (inverse operation).
pub fn remove_aberration(
    apparent_direction: Vector3,
    velocity_au_day: Vector3,
    sun_earth_distance_au: f64,
) -> Vector3 {
    // Iterative approach matching ERFA's eraAticq
    let mut d = Vector3::zeros();

    for _ in 0..2 {
        // Estimate the original direction
        let before = Vector3::new(
            apparent_direction.x - d.x,
            apparent_direction.y - d.y,
            apparent_direction.z - d.z,
        )
        .normalize();

        // Apply forward aberration to the estimate
        let after = apply_aberration(before, velocity_au_day, sun_earth_distance_au);

        // Update correction
        d = Vector3::new(after.x - before.x, after.y - before.y, after.z - before.z);
    }

    // Final result
    Vector3::new(
        apparent_direction.x - d.x,
        apparent_direction.y - d.y,
        apparent_direction.z - d.z,
    )
    .normalize()
}

pub fn apply_aberration(
    direction: Vector3,
    velocity_au_day: Vector3,
    sun_earth_distance_au: f64,
) -> Vector3 {
    let v = Vector3::new(
        velocity_au_day.x / SPEED_OF_LIGHT_AU_PER_DAY,
        velocity_au_day.y / SPEED_OF_LIGHT_AU_PER_DAY,
        velocity_au_day.z / SPEED_OF_LIGHT_AU_PER_DAY,
    );

    let v2 = v.x * v.x + v.y * v.y + v.z * v.z;
    let bm1 = (1.0 - v2).sqrt();

    let pdv = direction.dot(&v);
    let w1 = 1.0 + pdv / (1.0 + bm1);
    let w2 = SCHWARZSCHILD_RADIUS_SUN_AU / sun_earth_distance_au;

    let p2 = Vector3::new(
        direction.x * bm1 + w1 * v.x + w2 * (v.x - pdv * direction.x),
        direction.y * bm1 + w1 * v.y + w2 * (v.y - pdv * direction.y),
        direction.z * bm1 + w1 * v.z + w2 * (v.z - pdv * direction.z),
    );

    let r = p2.magnitude();
    Vector3::new(p2.x / r, p2.y / r, p2.z / r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_earth_velocity_j2000() {
        let tt = TT::j2000();
        let state = compute_earth_state(&tt).unwrap();

        let v_mag = state.barycentric_velocity.magnitude();
        assert!(
            (v_mag - 0.017).abs() < 0.001,
            "Barycentric velocity should be ~0.017 AU/day, got {}",
            v_mag
        );
    }

    #[test]
    fn test_aberration_magnitude() {
        let tt = TT::j2000();
        let state = compute_earth_state(&tt).unwrap();

        let p = Vector3::new(1.0, 0.0, 0.0);
        let s = state.heliocentric_position.magnitude();
        let p_ab = apply_aberration(p, state.barycentric_velocity, s);

        let disp_sq = (p_ab.x - p.x).powi(2) + (p_ab.y - p.y).powi(2) + (p_ab.z - p.z).powi(2);
        let aberr_arcsec = disp_sq.sqrt() * 206264.806247;

        assert!(
            aberr_arcsec < 25.0,
            "Aberration should be < 25\", got {}",
            aberr_arcsec
        );
    }

    #[test]
    fn test_aberration_roundtrip() {
        let tt = TT::j2000();
        let state = compute_earth_state(&tt).unwrap();
        let sun_dist = state.heliocentric_position.magnitude();

        let directions = [
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0).normalize(),
        ];

        for dir in directions {
            let aberrated = apply_aberration(dir, state.barycentric_velocity, sun_dist);
            let recovered = remove_aberration(aberrated, state.barycentric_velocity, sun_dist);

            let diff = ((dir.x - recovered.x).powi(2)
                + (dir.y - recovered.y).powi(2)
                + (dir.z - recovered.z).powi(2))
            .sqrt();

            // Iterative inverse gives ~1e-13 precision, not machine epsilon
            assert!(diff < 1e-12, "Aberration roundtrip error: {:.2e}", diff);
        }
    }

    #[test]
    fn test_remove_aberration_is_inverse() {
        let velocity = Vector3::new(0.01, 0.005, 0.002);
        let sun_dist = 1.0;
        let original = Vector3::new(0.6, 0.7, 0.3).normalize();

        let aberrated = apply_aberration(original, velocity, sun_dist);
        let recovered = remove_aberration(aberrated, velocity, sun_dist);

        let diff = ((original.x - recovered.x).powi(2)
            + (original.y - recovered.y).powi(2)
            + (original.z - recovered.z).powi(2))
        .sqrt();

        // Iterative inverse gives ~1e-13 precision, not machine epsilon
        assert!(
            diff < 1e-12,
            "Aberration inverse should be within iterative precision: {:.2e}",
            diff
        );
    }
}
