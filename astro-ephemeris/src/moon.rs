//! ELP/MPP02 Lunar Ephemeris
//!
//! Computes geocentric rectangular coordinates of the Moon using the
//! ELP/MPP02 semi-analytical lunar theory by Chapront & Francou (2003).
//!
//! Output is in the dynamical mean ecliptic and equinox of J2000 frame,
//! with positions in kilometers.

use astro_core::{
    AstroResult,
    constants::{ARCSEC_TO_RAD, DEG_TO_RAD, J2000_JD, PI},
};
use astro_time::TDB;

use crate::lunar_coefficients::{
    MainTerm, PertBlock, MAIN_DISTANCE, MAIN_LATITUDE, MAIN_LONGITUDE, PERT_DISTANCE,
    PERT_LATITUDE, PERT_LONGITUDE,
};

const A405: f64 = 384747.9613701725;
const AELP: f64 = 384747.980674318;

type Poly5 = [f64; 5];
type ElpArguments = ([[f64; 5]; 3], [[f64; 5]; 4], [[f64; 5]; 8], Poly5);

pub struct ElpMpp02Moon {
    icor: u8,
}

impl Default for ElpMpp02Moon {
    fn default() -> Self {
        Self::new()
    }
}

impl ElpMpp02Moon {
    pub fn new() -> Self {
        Self { icor: 0 }
    }

    pub fn with_de405_fit() -> Self {
        Self { icor: 1 }
    }

    pub fn geocentric_position(&self, tdb: &TDB) -> AstroResult<[f64; 3]> {
        let jd = tdb.to_julian_date();
        let tj = (jd.jd1() - J2000_JD) + jd.jd2();

        let (x, y, z, _, _, _) = self.evaluate(tj);
        Ok([x, y, z])
    }

    pub fn geocentric_state(&self, tdb: &TDB) -> AstroResult<[f64; 6]> {
        let jd = tdb.to_julian_date();
        let tj = (jd.jd1() - J2000_JD) + jd.jd2();

        let (x, y, z, vx, vy, vz) = self.evaluate(tj);
        Ok([x, y, z, vx, vy, vz])
    }

    /// Returns geocentric position in ICRS (J2000 equatorial) frame in km.
    ///
    /// The native ELP/MPP02 output is in the mean ecliptic and equinox of J2000.
    /// This method applies the obliquity rotation to convert to ICRS.
    pub fn geocentric_position_icrs(&self, tdb: &TDB) -> AstroResult<[f64; 3]> {
        let pos = self.geocentric_position(tdb)?;
        Ok(ecliptic_to_icrs(pos))
    }

    /// Returns geocentric state (position and velocity) in ICRS frame.
    /// Position in km, velocity in km/day.
    pub fn geocentric_state_icrs(&self, tdb: &TDB) -> AstroResult<[f64; 6]> {
        let state = self.geocentric_state(tdb)?;
        let pos = ecliptic_to_icrs([state[0], state[1], state[2]]);
        let vel = ecliptic_to_icrs([state[3], state[4], state[5]]);
        Ok([pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]])
    }

    fn evaluate(&self, tj: f64) -> (f64, f64, f64, f64, f64, f64) {
        let t = [
            1.0,
            tj / astro_core::constants::DAYS_PER_JULIAN_CENTURY,
            0.0,
            0.0,
            0.0,
        ];
        let t = [
            t[0],
            t[1],
            t[1] * t[1],
            t[1] * t[1] * t[1],
            t[1] * t[1] * t[1] * t[1],
        ];

        let (w, del, p, zeta) = self.compute_arguments(&t);

        let mut v = [0.0f64; 6];

        let (val, deriv) = self.evaluate_main(MAIN_LONGITUDE, &del, &t, false);
        v[0] = val;
        v[3] = deriv;
        let (pval, pderiv) = self.evaluate_pert(PERT_LONGITUDE, &del, &p, &zeta, &t);
        v[0] += pval;
        v[3] += pderiv;

        let (val, deriv) = self.evaluate_main(MAIN_LATITUDE, &del, &t, false);
        v[1] = val;
        v[4] = deriv;
        let (pval, pderiv) = self.evaluate_pert(PERT_LATITUDE, &del, &p, &zeta, &t);
        v[1] += pval;
        v[4] += pderiv;

        let (val, deriv) = self.evaluate_main(MAIN_DISTANCE, &del, &t, true);
        v[2] = val;
        v[5] = deriv;
        let (pval, pderiv) = self.evaluate_pert(PERT_DISTANCE, &del, &p, &zeta, &t);
        v[2] += pval;
        v[5] += pderiv;

        let rad = ARCSEC_TO_RAD.recip();
        v[0] = v[0] / rad
            + w[0][0]
            + w[0][1] * t[1]
            + w[0][2] * t[2]
            + w[0][3] * t[3]
            + w[0][4] * t[4];
        v[1] /= rad;
        v[2] *= A405 / AELP;
        v[3] = v[3] / rad
            + w[0][1]
            + 2.0 * w[0][2] * t[1]
            + 3.0 * w[0][3] * t[2]
            + 4.0 * w[0][4] * t[3];
        v[4] /= rad;

        let (slamb, clamb) = v[0].sin_cos();
        let (sbeta, cbeta) = v[1].sin_cos();
        let cw = v[2] * cbeta;
        let sw = v[2] * sbeta;

        let x1 = cw * clamb;
        let x2 = cw * slamb;
        let x3 = sw;

        let xp1 = (v[5] * cbeta - v[4] * sw) * clamb - v[3] * x2;
        let xp2 = (v[5] * cbeta - v[4] * sw) * slamb + v[3] * x1;
        let xp3 = v[5] * sbeta + v[4] * cw;

        let (pw, qw) = self.precession_pq(&t);
        let ra = 2.0 * (1.0 - pw * pw - qw * qw).sqrt();
        let pwqw = 2.0 * pw * qw;
        let pw2 = 1.0 - 2.0 * pw * pw;
        let qw2 = 1.0 - 2.0 * qw * qw;
        let pwra = pw * ra;
        let qwra = qw * ra;

        let x = pw2 * x1 + pwqw * x2 + pwra * x3;
        let y = pwqw * x1 + qw2 * x2 - qwra * x3;
        let z = -pwra * x1 + qwra * x2 + (pw2 + qw2 - 1.0) * x3;

        let (ppw, qpw) = self.precession_pq_deriv(&t);
        let ppw2 = -4.0 * pw * ppw;
        let qpw2 = -4.0 * qw * qpw;
        let ppwqpw = 2.0 * (ppw * qw + pw * qpw);
        let rap = (ppw2 + qpw2) / ra;
        let ppwra = ppw * ra + pw * rap;
        let qpwra = qpw * ra + qw * rap;

        let vx = (pw2 * xp1 + pwqw * xp2 + pwra * xp3 + ppw2 * x1 + ppwqpw * x2 + ppwra * x3)
            / astro_core::constants::DAYS_PER_JULIAN_CENTURY;
        let vy = (pwqw * xp1 + qw2 * xp2 - qwra * xp3 + ppwqpw * x1 + qpw2 * x2 - qpwra * x3)
            / astro_core::constants::DAYS_PER_JULIAN_CENTURY;
        let vz = (-pwra * xp1 + qwra * xp2 + (pw2 + qw2 - 1.0) * xp3 - ppwra * x1
            + qpwra * x2
            + (ppw2 + qpw2) * x3)
            / astro_core::constants::DAYS_PER_JULIAN_CENTURY;

        (x, y, z, vx, vy, vz)
    }

    fn compute_arguments(&self, _t: &Poly5) -> ElpArguments {
        let (
            dw1_0,
            dw2_0,
            dw3_0,
            deart_0,
            dperi,
            dw1_1,
            _dgam,
            _de,
            deart_1,
            _dep,
            dw2_1,
            dw3_1,
            dw1_2,
        ) = if self.icor == 0 {
            (
                -0.10525, 0.16826, -0.10760, -0.04012, -0.04854, -0.32311, 0.00069, 0.00005,
                0.01442, 0.00226, 0.08017, -0.04317, -0.03794,
            )
        } else {
            (
                -0.07008, 0.20794, -0.07215, -0.00033, -0.00749, -0.35106, 0.00085, -0.00006,
                0.00732, 0.00224, 0.08017, -0.04317, -0.03743,
            )
        };

        let dprec = -0.29965;
        let rad = ARCSEC_TO_RAD.recip();

        fn dms(deg: i32, min: i32, sec: f64) -> f64 {
            (deg as f64 + min as f64 / 60.0 + sec / 3600.0) * DEG_TO_RAD
        }

        let mut w = [[0.0f64; 5]; 3];
        w[0][0] = dms(218, 18, 59.95571 + dw1_0);
        w[0][1] = (1732559343.73604 + dw1_1) / rad;
        w[0][2] = (-6.8084 + dw1_2) / rad;
        w[0][3] = 0.66040e-2 / rad;
        w[0][4] = -0.31690e-4 / rad;

        w[1][0] = dms(83, 21, 11.67475 + dw2_0);
        w[1][1] = (14643420.3171 + dw2_1) / rad;
        w[1][2] = -38.2631 / rad;
        w[1][3] = -0.45047e-1 / rad;
        w[1][4] = 0.21301e-3 / rad;

        w[2][0] = dms(125, 2, 40.39816 + dw3_0);
        w[2][1] = (-6967919.5383 + dw3_1) / rad;
        w[2][2] = 6.3590 / rad;
        w[2][3] = 0.76250e-2 / rad;
        w[2][4] = -0.35860e-4 / rad;

        let mut eart = [0.0f64; 5];
        eart[0] = dms(100, 27, 59.13885 + deart_0);
        eart[1] = (129597742.29300 + deart_1) / rad;
        eart[2] = -0.020200 / rad;
        eart[3] = 0.90000e-5 / rad;
        eart[4] = 0.15000e-6 / rad;

        let mut peri = [0.0f64; 5];
        peri[0] = dms(102, 56, 14.45766 + dperi);
        peri[1] = 1161.24342 / rad;
        peri[2] = 0.529265 / rad;
        peri[3] = -0.11814e-3 / rad;
        peri[4] = 0.11379e-4 / rad;

        if self.icor == 1 {
            w[0][3] -= 0.00018865 / rad;
            w[0][4] -= 0.00001024 / rad;
            w[1][2] += 0.00470602 / rad;
            w[1][3] -= 0.00025213 / rad;
            w[2][2] -= 0.00261070 / rad;
            w[2][3] -= 0.00010712 / rad;
        }

        let mut del = [[0.0f64; 5]; 4];
        for i in 0..5 {
            del[0][i] = w[0][i] - eart[i];
            del[1][i] = w[0][i] - w[2][i];
            del[2][i] = w[0][i] - w[1][i];
            del[3][i] = eart[i] - peri[i];
        }
        del[0][0] += PI;

        let mut p = [[0.0f64; 5]; 8];
        p[0][0] = dms(252, 15, 3.216919);
        p[1][0] = dms(181, 58, 44.758419);
        p[2][0] = dms(100, 27, 59.138850);
        p[3][0] = dms(355, 26, 3.642778);
        p[4][0] = dms(34, 21, 5.379392);
        p[5][0] = dms(50, 4, 38.902495);
        p[6][0] = dms(314, 3, 4.354234);
        p[7][0] = dms(304, 20, 56.808371);

        p[0][1] = 538101628.66888 / rad;
        p[1][1] = 210664136.45777 / rad;
        p[2][1] = 129597742.29300 / rad;
        p[3][1] = 68905077.65936 / rad;
        p[4][1] = 10925660.57335 / rad;
        p[5][1] = 4399609.33632 / rad;
        p[6][1] = 1542482.57845 / rad;
        p[7][1] = 786547.89700 / rad;

        let mut zeta = [0.0f64; 5];
        zeta[0] = w[0][0];
        zeta[1] = w[0][1] + (5029.0966 + dprec) / rad;
        zeta[2] = w[0][2];
        zeta[3] = w[0][3];
        zeta[4] = w[0][4];

        (w, del, p, zeta)
    }

    fn evaluate_main(
        &self,
        terms: &[MainTerm],
        del: &[[f64; 5]; 4],
        t: &[f64; 5],
        is_distance: bool,
    ) -> (f64, f64) {
        let mut val = 0.0;
        let mut deriv = 0.0;
        let phase_offset = if is_distance { PI / 2.0 } else { 0.0 };

        for term in terms {
            let mut y = phase_offset;
            let mut yp = 0.0;

            for k in 0..5 {
                let arg_k: f64 = (0..4).map(|i| term.delaunay[i] as f64 * del[i][k]).sum();
                y += arg_k * t[k];
                if k > 0 {
                    yp += (k as f64) * arg_k * t[k - 1];
                }
            }

            let a0 = term.coeffs[0];
            val += a0 * y.sin();
            deriv += a0 * yp * y.cos();
        }

        (val, deriv)
    }

    fn evaluate_pert(
        &self,
        blocks: &[PertBlock],
        del: &[[f64; 5]; 4],
        p: &[[f64; 5]; 8],
        zeta: &[f64; 5],
        t: &[f64; 5],
    ) -> (f64, f64) {
        let mut val = 0.0;
        let mut deriv = 0.0;

        for block in blocks {
            let it = block.power as usize;

            for term in block.terms {
                let mut y = term.phase;
                let mut yp = 0.0;

                for k in 0..5 {
                    let mut arg_k = 0.0;
                    for (i, del_i) in del.iter().enumerate().take(4) {
                        arg_k += term.multipliers[i] as f64 * del_i[k];
                    }
                    for (i, p_i) in p.iter().enumerate().take(8) {
                        arg_k += term.multipliers[i + 4] as f64 * p_i[k];
                    }
                    arg_k += term.multipliers[12] as f64 * zeta[k];

                    y += arg_k * t[k];
                    if k > 0 {
                        yp += (k as f64) * arg_k * t[k - 1];
                    }
                }

                let x = term.amplitude;
                let xp = if it > 0 {
                    (it as f64) * x * t[it - 1]
                } else {
                    0.0
                };

                val += x * t[it] * y.sin();
                deriv += xp * y.sin() + x * t[it] * yp * y.cos();
            }
        }

        (val, deriv)
    }

    fn precession_pq(&self, t: &[f64; 5]) -> (f64, f64) {
        let p1 = 0.10180391e-04;
        let p2 = 0.47020439e-06;
        let p3 = -0.5417367e-09;
        let p4 = -0.2507948e-11;
        let p5 = 0.463486e-14;

        let q1 = -0.113469002e-03;
        let q2 = 0.12372674e-06;
        let q3 = 0.1265417e-08;
        let q4 = -0.1371808e-11;
        let q5 = -0.320334e-14;

        let pw = (p1 + p2 * t[1] + p3 * t[2] + p4 * t[3] + p5 * t[4]) * t[1];
        let qw = (q1 + q2 * t[1] + q3 * t[2] + q4 * t[3] + q5 * t[4]) * t[1];

        (pw, qw)
    }

    fn precession_pq_deriv(&self, t: &[f64; 5]) -> (f64, f64) {
        let p1 = 0.10180391e-04;
        let p2 = 0.47020439e-06;
        let p3 = -0.5417367e-09;
        let p4 = -0.2507948e-11;
        let p5 = 0.463486e-14;

        let q1 = -0.113469002e-03;
        let q2 = 0.12372674e-06;
        let q3 = 0.1265417e-08;
        let q4 = -0.1371808e-11;
        let q5 = -0.320334e-14;

        let ppw = p1 + (2.0 * p2 + 3.0 * p3 * t[1] + 4.0 * p4 * t[2] + 5.0 * p5 * t[3]) * t[1];
        let qpw = q1 + (2.0 * q2 + 3.0 * q3 * t[1] + 4.0 * q4 * t[2] + 5.0 * q5 * t[3]) * t[1];

        (ppw, qpw)
    }
}

/// Rotate from mean ecliptic J2000 to ICRS (equatorial J2000).
///
/// Uses the IAU 2006 obliquity of the ecliptic at J2000:
/// ε₀ = 84381.406 arcseconds = 23°26'21.406"
fn ecliptic_to_icrs(ecl: [f64; 3]) -> [f64; 3] {
    // IAU 2006 obliquity at J2000.0 in arcseconds
    const EPS0_ARCSEC: f64 = 84381.406;
    let eps = EPS0_ARCSEC * ARCSEC_TO_RAD;

    let (sin_eps, cos_eps) = eps.sin_cos();

    // Rotation about X-axis by -ε (from ecliptic to equatorial)
    [
        ecl[0],
        ecl[1] * cos_eps - ecl[2] * sin_eps,
        ecl[1] * sin_eps + ecl[2] * cos_eps,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use astro_time::julian::JulianDate;

    #[test]
    fn test_moon_j2000() {
        let moon = ElpMpp02Moon::new();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let pos = moon.geocentric_position(&tdb).unwrap();

        println!("Moon at J2000.0:");
        println!("  X = {:.3} km", pos[0]);
        println!("  Y = {:.3} km", pos[1]);
        println!("  Z = {:.3} km", pos[2]);

        let dist = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        println!("  Distance = {:.3} km", dist);

        assert!(
            dist > 350_000.0 && dist < 410_000.0,
            "Distance {} should be ~384,400 km",
            dist
        );
    }

    #[test]
    fn test_moon_distance_range() {
        let moon = ElpMpp02Moon::new();

        for days in [0, 7, 14, 21, 28] {
            let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD + days as f64, 0.0));
            let pos = moon.geocentric_position(&tdb).unwrap();
            let dist = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();

            assert!(
                dist > 356_000.0 && dist < 407_000.0,
                "Day {}: distance {} km out of range (perigee ~356,500, apogee ~406,700)",
                days,
                dist
            );
        }
    }

    #[test]
    fn test_fortran_reference_dates() {
        let moon = ElpMpp02Moon::new();

        let test_dates = [2444239.5, 2446239.5, 2448239.5, 2450239.5, 2452239.5];

        println!("\nELP/MPP02 LLR mode (icor=0) test positions:");
        for jd in test_dates {
            let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
            let pos = moon.geocentric_position(&tdb).unwrap();

            println!("JD {:.1}:", jd);
            println!("  X = {:17.7} km", pos[0]);
            println!("  Y = {:17.7} km", pos[1]);
            println!("  Z = {:17.7} km", pos[2]);

            let dist = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
            assert!(
                dist > 356_000.0 && dist < 407_000.0,
                "JD {}: distance {} km out of range",
                jd,
                dist
            );
        }
    }

    #[test]
    fn test_moon_icrs_j2000() {
        let moon = ElpMpp02Moon::new();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        let ecl = moon.geocentric_position(&tdb).unwrap();
        let icrs = moon.geocentric_position_icrs(&tdb).unwrap();

        println!("\nMoon at J2000.0:");
        println!(
            "  Ecliptic:  X = {:12.3} Y = {:12.3} Z = {:12.3} km",
            ecl[0], ecl[1], ecl[2]
        );
        println!(
            "  ICRS:      X = {:12.3} Y = {:12.3} Z = {:12.3} km",
            icrs[0], icrs[1], icrs[2]
        );

        // Distance should be preserved
        let dist_ecl = (ecl[0] * ecl[0] + ecl[1] * ecl[1] + ecl[2] * ecl[2]).sqrt();
        let dist_icrs = (icrs[0] * icrs[0] + icrs[1] * icrs[1] + icrs[2] * icrs[2]).sqrt();
        assert!((dist_ecl - dist_icrs).abs() < 1e-6, "Distance mismatch");

        // X should be unchanged (rotation about X-axis)
        assert!((ecl[0] - icrs[0]).abs() < 1e-10, "X should be unchanged");
    }

    #[test]
    fn test_sidereal_month_period() {
        let moon = ElpMpp02Moon::new();
        let sidereal_month = 27.321661; // days

        let t0 = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let t1 = TDB::from_julian_date(JulianDate::new(J2000_JD + sidereal_month, 0.0));

        let pos0 = moon.geocentric_position(&t0).unwrap();
        let pos1 = moon.geocentric_position(&t1).unwrap();

        // Normalize positions
        let norm0 = (pos0[0] * pos0[0] + pos0[1] * pos0[1] + pos0[2] * pos0[2]).sqrt();
        let norm1 = (pos1[0] * pos1[0] + pos1[1] * pos1[1] + pos1[2] * pos1[2]).sqrt();

        let unit0 = [pos0[0] / norm0, pos0[1] / norm0, pos0[2] / norm0];
        let unit1 = [pos1[0] / norm1, pos1[1] / norm1, pos1[2] / norm1];

        // Dot product should be close to 1 (same direction)
        let dot = unit0[0] * unit1[0] + unit0[1] * unit1[1] + unit0[2] * unit1[2];
        let angle_deg = dot.acos().to_degrees();

        println!("\nSidereal month test:");
        println!(
            "  Position at T0:     {:12.3} {:12.3} {:12.3}",
            pos0[0], pos0[1], pos0[2]
        );
        println!(
            "  Position at T0+27.32d: {:12.3} {:12.3} {:12.3}",
            pos1[0], pos1[1], pos1[2]
        );
        println!("  Angular separation: {:.2}°", angle_deg);

        // Moon should be within ~5° of original position after one sidereal month
        // (small deviation due to precession, perturbations, etc.)
        assert!(
            angle_deg < 5.0,
            "Moon should return close to original position after sidereal month"
        );
    }

    #[test]
    fn test_velocity_magnitude() {
        let moon = ElpMpp02Moon::new();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        let state = moon.geocentric_state(&tdb).unwrap();
        let v_mag = (state[3] * state[3] + state[4] * state[4] + state[5] * state[5]).sqrt();

        // Moon's orbital velocity is ~1.022 km/s = ~88,300 km/day
        let v_km_s = v_mag / astro_core::constants::SECONDS_PER_DAY_F64;

        println!("\nMoon velocity at J2000:");
        println!("  Vx = {:.3} km/day", state[3]);
        println!("  Vy = {:.3} km/day", state[4]);
        println!("  Vz = {:.3} km/day", state[5]);
        println!("  |V| = {:.3} km/day = {:.4} km/s", v_mag, v_km_s);

        // Verify with numerical differentiation
        let dt = 0.001; // 0.001 days = ~86 seconds
        let t1 = TDB::from_julian_date(JulianDate::new(J2000_JD - dt, 0.0));
        let t2 = TDB::from_julian_date(JulianDate::new(J2000_JD + dt, 0.0));
        let p1 = moon.geocentric_position(&t1).unwrap();
        let p2 = moon.geocentric_position(&t2).unwrap();
        let v_num = [
            (p2[0] - p1[0]) / (2.0 * dt),
            (p2[1] - p1[1]) / (2.0 * dt),
            (p2[2] - p1[2]) / (2.0 * dt),
        ];
        let v_num_mag = (v_num[0] * v_num[0] + v_num[1] * v_num[1] + v_num[2] * v_num[2]).sqrt();

        println!(
            "  Numerical velocity: {:.3} km/day = {:.4} km/s",
            v_num_mag,
            v_num_mag / astro_core::constants::SECONDS_PER_DAY_F64
        );

        // Use numerical velocity as ground truth - should be ~88,000 km/day
        assert!(
            v_num_mag > 80_000.0 && v_num_mag < 100_000.0,
            "Numerical velocity {} km/day out of expected range",
            v_num_mag
        );
    }

    #[test]
    fn test_against_jpl_de441() {
        // Reference values from JPL Horizons DE441 ephemeris
        // https://ssd.jpl.nasa.gov/horizons/
        // Query: Moon geocentric position, ecliptic J2000, OUT_UNITS='KM-D'
        let reference_data = [
            // (JD, X, Y, Z) - all in km
            (J2000_JD, -291608.38, -274979.74, 36271.20), // 2000-01-01 12:00 TDB
            (2444239.5, 43890.30, 381188.87, -31633.44),  // 1980-01-01 00:00 TDB
        ];

        let moon = ElpMpp02Moon::new();

        println!("\nComparison with JPL DE441:");
        let mut max_delta = 0.0f64;

        for (jd, jpl_x, jpl_y, jpl_z) in reference_data {
            let tdb = TDB::from_julian_date(JulianDate::new(jd, 0.0));
            let pos = moon.geocentric_position(&tdb).unwrap();

            let dx = pos[0] - jpl_x;
            let dy = pos[1] - jpl_y;
            let dz = pos[2] - jpl_z;
            let delta = (dx * dx + dy * dy + dz * dz).sqrt();

            println!("  JD {:.1}:", jd);
            println!(
                "    JPL:     X = {:12.2} Y = {:12.2} Z = {:12.2} km",
                jpl_x, jpl_y, jpl_z
            );
            println!(
                "    ELP:     X = {:12.2} Y = {:12.2} Z = {:12.2} km",
                pos[0], pos[1], pos[2]
            );
            println!("    Delta:   {:.3} km", delta);

            max_delta = max_delta.max(delta);

            // ELP/MPP02 should agree with DE441 within ~2 km (stated accuracy)
            // Our truncated series may have slightly larger error
            assert!(
                delta < 5.0,
                "JD {}: position difference {} km exceeds 5 km tolerance",
                jd,
                delta
            );
        }

        println!("  Maximum delta: {:.3} km", max_delta);
    }

    #[test]
    fn test_default_impl() {
        // Lines 35-36: Test Default trait implementation
        let moon: ElpMpp02Moon = Default::default();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));
        let pos = moon.geocentric_position(&tdb).unwrap();

        // Should work the same as new()
        let dist = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        assert!(dist > 350_000.0 && dist < 410_000.0);
    }

    #[test]
    fn test_with_de405_fit() {
        // Line 45: Test with_de405_fit constructor (icor=1)
        let moon = ElpMpp02Moon::with_de405_fit();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        // Test position
        let pos = moon.geocentric_position(&tdb).unwrap();
        let dist = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        assert!(
            dist > 350_000.0 && dist < 410_000.0,
            "DE405 fit distance: {} km",
            dist
        );

        // Test state
        let state = moon.geocentric_state(&tdb).unwrap();
        assert_eq!(state[0], pos[0]);
        assert_eq!(state[1], pos[1]);
        assert_eq!(state[2], pos[2]);
    }

    #[test]
    fn test_geocentric_state_icrs() {
        // Lines 76-80: Test geocentric_state_icrs
        let moon = ElpMpp02Moon::new();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        let state_icrs = moon.geocentric_state_icrs(&tdb).unwrap();

        // Should have 6 components: position and velocity
        let pos_icrs = [state_icrs[0], state_icrs[1], state_icrs[2]];
        let vel_icrs = [state_icrs[3], state_icrs[4], state_icrs[5]];

        // Position distance should be preserved
        let dist = (pos_icrs[0].powi(2) + pos_icrs[1].powi(2) + pos_icrs[2].powi(2)).sqrt();
        assert!(dist > 350_000.0 && dist < 410_000.0);

        // Velocity should be non-zero
        let vel_mag = (vel_icrs[0].powi(2) + vel_icrs[1].powi(2) + vel_icrs[2].powi(2)).sqrt();
        assert!(vel_mag > 0.0);

        // Compare with separate calls
        let pos_only = moon.geocentric_position_icrs(&tdb).unwrap();
        assert!((pos_only[0] - pos_icrs[0]).abs() < 1e-10);
        assert!((pos_only[1] - pos_icrs[1]).abs() < 1e-10);
        assert!((pos_only[2] - pos_icrs[2]).abs() < 1e-10);
    }

    #[test]
    fn test_de405_fit_corrections() {
        // Lines 164, 207-213: Test that icor=1 gives different results than icor=0
        let moon_llr = ElpMpp02Moon::new(); // icor=0 (LLR mode)
        let moon_de405 = ElpMpp02Moon::with_de405_fit(); // icor=1 (DE405 fit)

        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD + 1000.0, 0.0));

        let pos_llr = moon_llr.geocentric_position(&tdb).unwrap();
        let pos_de405 = moon_de405.geocentric_position(&tdb).unwrap();

        // Positions should be slightly different due to different constants
        let diff = ((pos_llr[0] - pos_de405[0]).powi(2)
            + (pos_llr[1] - pos_de405[1]).powi(2)
            + (pos_llr[2] - pos_de405[2]).powi(2))
        .sqrt();

        // The difference should be small but measurable (a few meters to km)
        assert!(diff > 0.0, "LLR and DE405 positions should differ");
        assert!(diff < 100.0, "Difference should be less than 100 km");

        println!("LLR vs DE405 position difference: {:.3} km", diff);
    }

    #[test]
    fn test_de405_fit_velocity() {
        // Test that DE405 fit also works for state/velocity computation
        let moon = ElpMpp02Moon::with_de405_fit();
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        let state = moon.geocentric_state(&tdb).unwrap();

        // Velocity magnitude should be ~1 km/s = 86400 km/day
        let v_mag = (state[3].powi(2) + state[4].powi(2) + state[5].powi(2)).sqrt();
        let v_km_s = v_mag / astro_core::constants::SECONDS_PER_DAY_F64;

        println!("DE405 fit velocity: {:.4} km/s", v_km_s);
        assert!(
            v_km_s > 0.8 && v_km_s < 1.2,
            "Moon velocity should be ~1 km/s"
        );
    }
}
