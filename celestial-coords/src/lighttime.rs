use crate::Distance;
use celestial_core::Vector3;

const C_AU_PER_DAY: f64 = 173.1446326846693;

pub struct LightTimeCorrection {
    light_time_days: f64,
}

impl LightTimeCorrection {
    pub fn from_distance(distance: Distance) -> Self {
        let distance_au = distance.au();
        let light_time_days = distance_au / C_AU_PER_DAY;

        Self { light_time_days }
    }

    pub fn from_position_vector(pos: Vector3) -> Self {
        let distance_au = libm::sqrt(pos.x.powi(2) + pos.y.powi(2) + pos.z.powi(2));
        let light_time_days = distance_au / C_AU_PER_DAY;

        Self { light_time_days }
    }

    pub fn light_time_days(&self) -> f64 {
        self.light_time_days
    }

    pub fn light_time_seconds(&self) -> f64 {
        self.light_time_days * celestial_core::constants::SECONDS_PER_DAY_F64
    }

    pub fn apply_proper_motion(
        &self,
        pm_ra_mas_per_year: f64,
        pm_dec_mas_per_year: f64,
    ) -> (f64, f64) {
        let years = self.light_time_days / 365.25;

        let delta_ra_mas = pm_ra_mas_per_year * years;
        let delta_dec_mas = pm_dec_mas_per_year * years;

        (delta_ra_mas, delta_dec_mas)
    }

    pub fn apply_radial_velocity(&self, radial_velocity_km_s: f64) -> f64 {
        radial_velocity_km_s * self.light_time_seconds()
    }

    /// Corrects a position vector for light-time delay.
    ///
    /// # Arguments
    /// * `pos` - Position vector in AU
    /// * `vel` - Velocity vector in AU/day
    ///
    /// # Returns
    /// Corrected position vector in AU
    pub fn correct_position_vector(pos: Vector3, vel: Vector3) -> Vector3 {
        let lt = Self::from_position_vector(pos);
        let t = lt.light_time_days();

        Vector3::new(pos.x - vel.x * t, pos.y - vel.y * t, pos.z - vel.z * t)
    }

    /// Iteratively corrects a position vector for light-time delay with convergence checking.
    ///
    /// # Arguments
    /// * `pos` - Position vector in AU
    /// * `vel` - Velocity vector in AU/day
    /// * `max_iterations` - Maximum number of iterations
    ///
    /// # Returns
    /// Converged position vector in AU
    pub fn iterate_correction(pos: Vector3, vel: Vector3, max_iterations: usize) -> Vector3 {
        let mut corrected = pos;

        for _ in 0..max_iterations {
            let lt = Self::from_position_vector(corrected);
            let t = lt.light_time_days();

            let new_corrected =
                Vector3::new(pos.x - vel.x * t, pos.y - vel.y * t, pos.z - vel.z * t);

            let delta = libm::sqrt(
                (new_corrected.x - corrected.x).powi(2)
                    + (new_corrected.y - corrected.y).powi(2)
                    + (new_corrected.z - corrected.z).powi(2),
            );

            corrected = new_corrected;

            if delta < 1e-12 {
                break;
            }
        }

        corrected
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_light_time_from_distance() {
        let distance = Distance::from_au(1.0).unwrap();
        let lt = LightTimeCorrection::from_distance(distance);

        assert!((lt.light_time_days() - 1.0 / C_AU_PER_DAY).abs() < 1e-12);
        assert!((lt.light_time_seconds() - 499.00478).abs() < 1.0);
    }

    #[test]
    fn test_proper_motion_correction() {
        let distance = Distance::from_parsecs(10.0).unwrap();
        let lt = LightTimeCorrection::from_distance(distance);

        let (delta_ra, delta_dec) = lt.apply_proper_motion(100.0, 50.0);

        assert!(delta_ra > 0.0);
        assert!(delta_dec > 0.0);
    }

    #[test]
    fn test_radial_velocity_correction() {
        let distance = Distance::from_au(1.0).unwrap();
        let lt = LightTimeCorrection::from_distance(distance);

        let distance_change_km = lt.apply_radial_velocity(30.0);

        assert!(distance_change_km > 0.0);
    }

    #[test]
    fn test_position_vector_correction() {
        let pos = Vector3::new(1.0, 0.0, 0.0);
        let vel = Vector3::new(0.1, 0.0, 0.0);

        let corrected = LightTimeCorrection::correct_position_vector(pos, vel);

        assert!(corrected.x < pos.x);
    }

    #[test]
    fn test_iterative_correction() {
        let pos = Vector3::new(5.0, 0.0, 0.0);
        let vel = Vector3::new(0.01, 0.0, 0.0);

        let corrected = LightTimeCorrection::iterate_correction(pos, vel, 10);

        assert!(corrected.x < pos.x);
    }

    #[test]
    fn test_jupiter_light_time() {
        let jupiter_distance = Distance::from_au(5.2).unwrap();
        let lt = LightTimeCorrection::from_distance(jupiter_distance);

        let expected_minutes = 5.2 * 499.0 / 60.0;
        let actual_minutes = lt.light_time_seconds() / 60.0;

        assert!((actual_minutes - expected_minutes).abs() < 1.0);
    }

    #[test]
    fn test_romer_jupiter_moons_1676() {
        let near_opposition = Distance::from_au(5.2 - 1.0).unwrap();
        let far_opposition = Distance::from_au(5.2 + 1.0).unwrap();

        let lt_near = LightTimeCorrection::from_distance(near_opposition);
        let lt_far = LightTimeCorrection::from_distance(far_opposition);

        let delay_seconds = lt_far.light_time_seconds() - lt_near.light_time_seconds();
        let delay_minutes = delay_seconds / 60.0;

        assert!((delay_minutes - 16.6).abs() < 1.0);
    }

    #[test]
    fn test_barnards_star_proper_motion() {
        let distance = Distance::from_parsecs(1.83).unwrap();
        let lt = LightTimeCorrection::from_distance(distance);

        let pm_ra = -798.58;
        let pm_dec = 10328.12;

        let (delta_ra, delta_dec) = lt.apply_proper_motion(pm_ra, pm_dec);

        let light_years = lt.light_time_days() / 365.25;
        assert!((delta_ra - (pm_ra * light_years)).abs() < 0.01);
        assert!((delta_dec - (pm_dec * light_years)).abs() < 0.01);
    }

    #[test]
    fn test_proxima_centauri_radial_velocity() {
        let distance = Distance::from_parsecs(1.30).unwrap();
        let lt = LightTimeCorrection::from_distance(distance);

        let rv_km_s = -22.2;

        let distance_change_km = lt.apply_radial_velocity(rv_km_s);

        let expected_km = rv_km_s * lt.light_time_seconds();
        assert!((distance_change_km - expected_km).abs() < 1e-6);
    }

    #[test]
    fn test_sun_earth_light_time() {
        let earth_sun = Distance::from_au(1.0).unwrap();
        let lt = LightTimeCorrection::from_distance(earth_sun);

        assert!((lt.light_time_seconds() - 499.0).abs() < 1.0);
        assert!((lt.light_time_days() - 0.00577).abs() < 0.0001);
    }

    #[test]
    fn test_mercury_transit_timing() {
        let mercury_near = Distance::from_au(0.31).unwrap();
        let mercury_far = Distance::from_au(0.47).unwrap();

        let lt_near = LightTimeCorrection::from_distance(mercury_near);
        let lt_far = LightTimeCorrection::from_distance(mercury_far);

        let timing_difference = (lt_far.light_time_seconds() - lt_near.light_time_seconds()).abs();

        assert!(timing_difference > 0.0);
        assert!(timing_difference < 100.0);
    }

    #[test]
    fn test_convergence_high_velocity_object() {
        let pos = Vector3::new(100.0, 0.0, 0.0);
        let vel = Vector3::new(1.0, 0.0, 0.0);

        let single = LightTimeCorrection::correct_position_vector(pos, vel);
        let iterated = LightTimeCorrection::iterate_correction(pos, vel, 10);

        let improvement = (iterated.x - single.x).abs();
        assert!(improvement > 0.0);
    }

    #[test]
    fn test_aberration_effect_simulation() {
        let distance = Distance::from_au(1.0).unwrap();
        let lt = LightTimeCorrection::from_distance(distance);

        let earth_orbital_velocity = 29.78;
        let distance_change_km = lt.apply_radial_velocity(earth_orbital_velocity);

        let distance_change_au = distance_change_km / 1.495978707e8;
        assert!(distance_change_au.abs() < 0.1);
    }
}
