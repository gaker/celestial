#[cfg(test)]
mod tests;

use astro_coords::{CartesianFrame, EclipticCartesian, Vector3};
use astro_core::constants::{DAYS_PER_JULIAN_MILLENNIUM, J2000_JD, TWOPI};
use astro_core::AstroResult;

use crate::earth::Vsop2013Earth;
use crate::planetary_coefficients::*;

use astro_time::julian::JulianDate;
use astro_time::TDB;

const DT_DAYS: f64 = 1.0 / astro_core::constants::SECONDS_PER_DAY_F64;

#[allow(clippy::excessive_precision)]
const LAMBDA0: [f64; 17] = [
    4.402608631669, // Mercury
    3.176134461576, // Venus
    1.753470369433, // Earth-Moon Barycenter
    6.203500014141, // Mars
    4.091360003050, // Vesta
    1.713740719173, // Iris
    5.598641292287, // Bamberga
    2.805136360408, // Ceres
    2.326989734620, // Pallas
    0.599546107035, // Jupiter
    0.874018510107, // Saturn
    5.481225395663, // Uranus
    5.311897933164, // Neptune
    0.0,            // Pluto mean longitude (unused)
    5.19846640063,  // Moon D
    1.62790513602,  // Moon F
    2.35555563875,  // Moon l
];

#[allow(clippy::excessive_precision)]
const LAMBDA_DOT: [f64; 17] = [
    26087.90314068555,  // Mercury (rad/millennium)
    10213.28554743445,  // Venus
    6283.075850353215,  // Earth-Moon Barycenter
    3340.612434145457,  // Mars
    1731.170452721855,  // Vesta
    1704.450855027201,  // Iris
    1428.948917844273,  // Bamberga
    1364.756513629990,  // Ceres
    1361.923207632842,  // Pallas
    529.6909615623250,  // Jupiter
    213.2990861084880,  // Saturn
    74.78165903077800,  // Uranus
    38.13297222612500,  // Neptune
    0.3595362285049309, // Pluto (mu from TOP2013)
    77713.7714481804,   // Moon D
    84334.6615717837,   // Moon F
    83286.9142477147,   // Moon l
];

fn precompute_lambdas(t: f64) -> [f64; 17] {
    let mut lambdas = [0.0; 17];
    for i in 0..17 {
        lambdas[i] = LAMBDA0[i] + LAMBDA_DOT[i] * t;
    }
    lambdas
}

fn evaluate_variable(blocks: &[TimeBlock], t: f64, lambdas: &[f64; 17]) -> f64 {
    let mut result = 0.0;
    for block in blocks {
        let mut block_sum = 0.0;
        for term in block.terms {
            let arg = compute_argument(&term.mult, lambdas);
            let (arg_sin, arg_cos) = arg.sin_cos();
            block_sum += term.s * arg_sin + term.c * arg_cos;
        }
        result += block_sum * t.powi(block.power as i32);
    }
    result
}

fn compute_argument(mult: &[i16; 17], lambdas: &[f64; 17]) -> f64 {
    let mut arg = 0.0;
    for i in 0..17 {
        if mult[i] != 0 {
            arg += (mult[i] as f64) * lambdas[i];
        }
    }
    arg
}

fn elements_to_cartesian(el: &OrbitalElements) -> AstroResult<Vector3> {
    let (xa, xl, xk, xh, xq, xp) = (el.a, el.lambda, el.k, el.h, el.q, el.p);

    let xfi = (1.0 - xk * xk - xh * xh).sqrt();
    let xki = (1.0 - xq * xq - xp * xp).sqrt();
    let u = 1.0 / (1.0 + xfi);

    let ex = (xk * xk + xh * xh).sqrt();
    let gl = xl % TWOPI;
    let gm = gl - xh.atan2(xk);

    let mut e_anom = gl
        + (ex - 0.125 * ex.powi(3)) * gm.sin()
        + 0.5 * ex.powi(2) * (2.0 * gm).sin()
        + 0.375 * ex.powi(3) * (3.0 * gm).sin();

    loop {
        let (sin_e, cos_e) = e_anom.sin_cos();
        let z3_real = xk * cos_e + xh * sin_e;
        let z3_imag = xk * sin_e - xh * cos_e;
        let dl = gl - e_anom + z3_imag;
        e_anom += dl / (1.0 - z3_real);
        if dl.abs() < 1e-15 {
            break;
        }
    }

    let (sin_e, cos_e) = e_anom.sin_cos();
    let z3_real = xk * cos_e + xh * sin_e;
    let z3_imag = xk * sin_e - xh * cos_e;
    let rsa = 1.0 - z3_real;

    let z1_real = u * xk * z3_imag;
    let z1_imag = u * xh * z3_imag;
    let zto_real = (-xk + cos_e + z1_imag) / rsa;
    let zto_imag = (-xh + sin_e - z1_real) / rsa;

    let xm = xp * zto_real - xq * zto_imag;
    let xr = xa * rsa;

    Ok(Vector3::new(
        xr * (zto_real - 2.0 * xp * xm),
        xr * (zto_imag + 2.0 * xq * xm),
        -2.0 * xr * xki * xm,
    ))
}

fn central_difference_velocity(p_minus: &Vector3, p_plus: &Vector3) -> Vector3 {
    let inv_2dt = 1.0 / (2.0 * DT_DAYS);
    Vector3::new(
        (p_plus.x - p_minus.x) * inv_2dt,
        (p_plus.y - p_minus.y) * inv_2dt,
        (p_plus.z - p_minus.z) * inv_2dt,
    )
}

fn offset_tdb(tdb: &TDB, dt_days: f64) -> TDB {
    let jd = tdb.to_julian_date();
    TDB::from_julian_date(JulianDate::new(jd.jd1(), jd.jd2() + dt_days))
}

macro_rules! impl_vsop2013_planet {
    ($name:ident, $coeffs:ident) => {
        pub struct $name;

        impl $name {
            pub fn heliocentric_position(&self, tdb: &TDB) -> AstroResult<Vector3> {
                let jd = tdb.to_julian_date();
                let t = (jd.jd1() + jd.jd2() - J2000_JD) / DAYS_PER_JULIAN_MILLENNIUM;
                let lambdas = precompute_lambdas(t);

                let a = evaluate_variable($coeffs::A, t, &lambdas);
                let mut lambda = evaluate_variable($coeffs::LAMBDA, t, &lambdas) % TWOPI;
                if lambda < 0.0 {
                    lambda += TWOPI;
                }
                let k = evaluate_variable($coeffs::K, t, &lambdas);
                let h = evaluate_variable($coeffs::H, t, &lambdas);
                let q = evaluate_variable($coeffs::Q, t, &lambdas);
                let p = evaluate_variable($coeffs::P, t, &lambdas);

                let el = OrbitalElements {
                    a,
                    lambda,
                    k,
                    h,
                    q,
                    p,
                };
                let pos_ecl = elements_to_cartesian(&el)?;
                Ok(EclipticCartesian::from_vector3(&pos_ecl).to_icrs())
            }

            pub fn heliocentric_state(&self, tdb: &TDB) -> AstroResult<(Vector3, Vector3)> {
                let pos = self.heliocentric_position(tdb)?;
                let t_minus = offset_tdb(tdb, -DT_DAYS);
                let t_plus = offset_tdb(tdb, DT_DAYS);
                let p_minus = self.heliocentric_position(&t_minus)?;
                let p_plus = self.heliocentric_position(&t_plus)?;
                let vel = central_difference_velocity(&p_minus, &p_plus);
                Ok((pos, vel))
            }

            pub fn geocentric_position(&self, tdb: &TDB) -> AstroResult<Vector3> {
                let planet_helio = self.heliocentric_position(tdb)?;
                let earth_helio = Vsop2013Earth::new().heliocentric_position(tdb)?;
                Ok(Vector3::new(
                    planet_helio.x - earth_helio.x,
                    planet_helio.y - earth_helio.y,
                    planet_helio.z - earth_helio.z,
                ))
            }

            pub fn geocentric_state(&self, tdb: &TDB) -> AstroResult<(Vector3, Vector3)> {
                let pos = self.geocentric_position(tdb)?;
                let t_minus = offset_tdb(tdb, -DT_DAYS);
                let t_plus = offset_tdb(tdb, DT_DAYS);
                let p_minus = self.geocentric_position(&t_minus)?;
                let p_plus = self.geocentric_position(&t_plus)?;
                let vel = central_difference_velocity(&p_minus, &p_plus);
                Ok((pos, vel))
            }
        }
    };
}

impl_vsop2013_planet!(Vsop2013Mercury, mercury);
impl_vsop2013_planet!(Vsop2013Venus, venus);
impl_vsop2013_planet!(Vsop2013Mars, mars);
impl_vsop2013_planet!(Vsop2013Jupiter, jupiter);
impl_vsop2013_planet!(Vsop2013Saturn, saturn);
impl_vsop2013_planet!(Vsop2013Uranus, uranus);
impl_vsop2013_planet!(Vsop2013Neptune, neptune);
impl_vsop2013_planet!(Vsop2013Pluto, pluto);

impl_vsop2013_planet!(Vsop2013Emb, emb);

struct OrbitalElements {
    a: f64,
    lambda: f64,
    k: f64,
    h: f64,
    q: f64,
    p: f64,
}

#[cfg(test)]
mod test {
    use super::*;
    use astro_core::constants::J2000_JD;

    #[test]
    fn test_elements_to_cartesian() {
        let el = OrbitalElements {
            a: 39.2648542648,
            lambda: 4.1726045776,
            k: -0.1758641167,
            h: -0.1701234143,
            q: -0.0517015914,
            p: 0.1398654514,
        };

        let pos_ecl = elements_to_cartesian(&el).unwrap();
        let expected = [-9.8753625435, -27.9588613710, 5.8504463318];

        for i in 0..3 {
            let diff = (pos_ecl[i] - expected[i]).abs();
            assert!(diff < 1e-8, "Component {} error {:.2e} AU", i, diff);
        }
    }

    #[test]
    fn test_ecliptic_to_icrs() {
        let ecl = Vector3::new(-9.8753625435, -27.9588613710, 5.8504463318);
        let icrs = EclipticCartesian::from_vector3(&ecl).to_icrs();

        assert!(icrs.x < 0.0, "X should be negative");
        assert!(icrs.y < 0.0, "Y should be negative");
        assert!(icrs.z < 0.0, "Z should be negative in ICRS");
    }

    #[test]
    fn test_mars_heliocentric_velocity() {
        let mars = Vsop2013Mars;
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        let (pos, vel) = mars.heliocentric_state(&tdb).unwrap();
        let vel_mag = (vel.x.powi(2) + vel.y.powi(2) + vel.z.powi(2)).sqrt();

        println!("Mars at J2000.0:");
        println!("  Position: ({:.6}, {:.6}, {:.6}) AU", pos.x, pos.y, pos.z);
        println!(
            "  Velocity: ({:.6}, {:.6}, {:.6}) AU/day",
            vel.x, vel.y, vel.z
        );
        println!("  |V| = {:.6} AU/day", vel_mag);

        assert!(
            vel_mag > 0.01 && vel_mag < 0.03,
            "Mars heliocentric velocity {} AU/day should be ~0.014-0.027 AU/day",
            vel_mag
        );

        let pos_mag = (pos.x.powi(2) + pos.y.powi(2) + pos.z.powi(2)).sqrt();
        let dot = pos.x * vel.x + pos.y * vel.y + pos.z * vel.z;
        let cos_angle = dot / (pos_mag * vel_mag);
        let angle_deg = cos_angle.acos().to_degrees();

        println!("  Angle between position and velocity: {:.1}°", angle_deg);
        assert!(
            angle_deg > 60.0 && angle_deg < 120.0,
            "Velocity should be roughly tangent to orbit (angle ~90°), got {}°",
            angle_deg
        );
    }

    #[test]
    fn test_mars_velocity_against_de432s() {
        use crate::jpl::{bodies, SpkFile};

        const AU_KM: f64 = 149597870.7;
        const SECONDS_PER_DAY: f64 = astro_core::constants::SECONDS_PER_DAY_F64;

        let spk_path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/de432s.bsp");
        let spk = SpkFile::open(&spk_path).expect("Failed to open de432s.bsp");

        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        let (mars_bary_pos_km, mars_bary_vel_kms) = spk
            .compute_state(
                bodies::MARS_BARYCENTER,
                bodies::SOLAR_SYSTEM_BARYCENTER,
                J2000_JD,
            )
            .expect("Failed to get Mars barycentric state");
        let (sun_bary_pos_km, sun_bary_vel_kms) = spk
            .compute_state(bodies::SUN, bodies::SOLAR_SYSTEM_BARYCENTER, J2000_JD)
            .expect("Failed to get Sun barycentric state");

        let de_helio_pos_au = [
            (mars_bary_pos_km[0] - sun_bary_pos_km[0]) / AU_KM,
            (mars_bary_pos_km[1] - sun_bary_pos_km[1]) / AU_KM,
            (mars_bary_pos_km[2] - sun_bary_pos_km[2]) / AU_KM,
        ];
        let de_helio_vel_au_day = [
            (mars_bary_vel_kms[0] - sun_bary_vel_kms[0]) * SECONDS_PER_DAY / AU_KM,
            (mars_bary_vel_kms[1] - sun_bary_vel_kms[1]) * SECONDS_PER_DAY / AU_KM,
            (mars_bary_vel_kms[2] - sun_bary_vel_kms[2]) * SECONDS_PER_DAY / AU_KM,
        ];

        let mars = Vsop2013Mars;
        let (vsop_pos, vsop_vel) = mars.heliocentric_state(&tdb).unwrap();

        let pos_error_au = [
            vsop_pos.x - de_helio_pos_au[0],
            vsop_pos.y - de_helio_pos_au[1],
            vsop_pos.z - de_helio_pos_au[2],
        ];
        let pos_error_km =
            (pos_error_au[0].powi(2) + pos_error_au[1].powi(2) + pos_error_au[2].powi(2)).sqrt()
                * AU_KM;

        let vel_error_au_day = [
            vsop_vel.x - de_helio_vel_au_day[0],
            vsop_vel.y - de_helio_vel_au_day[1],
            vsop_vel.z - de_helio_vel_au_day[2],
        ];
        let vel_error_mag = (vel_error_au_day[0].powi(2)
            + vel_error_au_day[1].powi(2)
            + vel_error_au_day[2].powi(2))
        .sqrt();

        assert!(
            pos_error_km < 1000.0,
            "Mars position error {:.1} km exceeds 1000 km tolerance",
            pos_error_km
        );
        assert!(
            vel_error_mag < 0.0001,
            "Mars velocity error {:.6} AU/day exceeds 0.0001 AU/day tolerance",
            vel_error_mag
        );
    }

    #[test]
    fn test_mars_geocentric_state() {
        let mars = Vsop2013Mars;
        let tdb = TDB::from_julian_date(JulianDate::new(J2000_JD, 0.0));

        let (pos, vel) = mars.geocentric_state(&tdb).unwrap();

        let dist_au = (pos.x.powi(2) + pos.y.powi(2) + pos.z.powi(2)).sqrt();
        assert!(
            dist_au > 0.5 && dist_au < 2.7,
            "Mars geocentric distance {} AU outside expected range",
            dist_au
        );

        let speed_au_day = (vel.x.powi(2) + vel.y.powi(2) + vel.z.powi(2)).sqrt();
        assert!(
            speed_au_day > 0.001 && speed_au_day < 0.05,
            "Mars geocentric velocity {} AU/day outside expected range",
            speed_au_day
        );
    }
}
