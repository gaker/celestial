mod space_motion;

use crate::errors::{CoordError, CoordResult};
use celestial_core::angle::{wrap_0_2pi, Angle};
use celestial_core::constants::{
    ARCSEC_PER_RAD, AU_M, DAYS_PER_JULIAN_YEAR, MILLIARCSEC_TO_RAD, SECONDS_PER_DAY_F64,
};
use celestial_core::math::angular_separation;
use celestial_time::scales::tt::TT;
use space_motion::{SpaceMotion, SphericalMotion};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// As in SOFA's pmsafe, a small, zero or negative parallax (arcsec) is raised to the star's
// yearly motion (radians) times this factor, which holds the implied transverse speed near
// 1% of c, and to at least the floor.
const SAFE_PARALLAX_FACTOR: f64 = 326.0;
const PARALLAX_FLOOR_ARCSEC: f64 = 5e-7;

#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CatalogStar {
    pub ra: Angle,
    pub dec: Angle,
    pub pm_ra_cos_dec_mas_per_year: f64,
    pub pm_dec_mas_per_year: f64,
    pub parallax_mas: f64,
    pub radial_velocity_km_per_s: f64,
}

impl CatalogStar {
    pub fn propagate(&self, from: &TT, to: &TT) -> CoordResult<Self> {
        let observed = self.observed_motion()?;
        let days = interval_days(from, to)?;
        let moved = SpaceMotion::from_observed(&observed)?.propagate(days)?;
        Self::from_observed_motion(&moved.to_observed()?)
    }

    fn observed_motion(&self) -> CoordResult<SphericalMotion> {
        self.check_finite()?;
        let (ra, dec) = (self.ra.radians(), self.dec.validate_latitude()?.radians());
        let pm_ra = self.pm_ra_cos_dec_mas_per_year * MILLIARCSEC_TO_RAD / libm::cos(dec);
        let pm_dec = self.pm_dec_mas_per_year * MILLIARCSEC_TO_RAD;
        let parallax = safe_parallax(ra, dec, pm_ra, pm_dec, self.parallax_mas / 1e3);
        Ok(SphericalMotion {
            ra,
            dec,
            distance: ARCSEC_PER_RAD / parallax,
            ra_rate: pm_ra / DAYS_PER_JULIAN_YEAR,
            dec_rate: pm_dec / DAYS_PER_JULIAN_YEAR,
            radial_rate: SECONDS_PER_DAY_F64 * self.radial_velocity_km_per_s * 1e3 / AU_M,
        })
    }

    fn from_observed_motion(m: &SphericalMotion) -> CoordResult<Self> {
        let pm_ra = m.ra_rate * DAYS_PER_JULIAN_YEAR;
        let pm_dec = m.dec_rate * DAYS_PER_JULIAN_YEAR;
        let parallax = ARCSEC_PER_RAD / m.distance;
        Ok(Self {
            ra: Angle::from_radians(wrap_0_2pi(m.ra)?),
            dec: Angle::from_radians(m.dec),
            pm_ra_cos_dec_mas_per_year: pm_ra * libm::cos(m.dec) / MILLIARCSEC_TO_RAD,
            pm_dec_mas_per_year: pm_dec / MILLIARCSEC_TO_RAD,
            parallax_mas: parallax * 1e3,
            radial_velocity_km_per_s: 1e-3 * m.radial_rate * AU_M / SECONDS_PER_DAY_F64,
        })
    }

    fn check_finite(&self) -> CoordResult<()> {
        let values = [
            self.ra.radians(),
            self.pm_ra_cos_dec_mas_per_year,
            self.pm_dec_mas_per_year,
            self.parallax_mas,
            self.radial_velocity_km_per_s,
        ];
        if values.iter().all(|v| v.is_finite()) {
            Ok(())
        } else {
            Err(CoordError::invalid_coordinate(
                "Catalog star value is not finite",
            ))
        }
    }
}

fn safe_parallax(ra: f64, dec: f64, pm_ra: f64, pm_dec: f64, parallax: f64) -> f64 {
    let yearly = angular_separation(ra, dec, ra + pm_ra, dec + pm_dec);
    libm::fmax(
        libm::fmax(parallax, yearly * SAFE_PARALLAX_FACTOR),
        PARALLAX_FLOOR_ARCSEC,
    )
}

// SOFA takes TDB epochs. TT differs from TDB by under 2 ms, far below anything proper motion
// can show.
fn interval_days(from: &TT, to: &TT) -> CoordResult<f64> {
    let (from, to) = (from.to_julian_date(), to.to_julian_date());
    let days = (to.jd1() - from.jd1()) + (to.jd2() - from.jd2());
    if days.is_finite() {
        Ok(days)
    } else {
        Err(CoordError::invalid_coordinate("Epoch is not finite"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::{DAYS_PER_JULIAN_CENTURY, J2000_JD, SPEED_OF_LIGHT_M_PER_S};
    use celestial_time::julian::JulianDate;

    const J2016: f64 = 2457389.0;
    const BARNARD: [f64; 6] = [
        4.702763489559941,
        0.08271813456901925,
        -801.55,
        10362.39,
        546.98,
        -110.5,
    ];

    fn tt(jd1: f64, jd2: f64) -> TT {
        TT::from_julian_date(JulianDate::new(jd1, jd2))
    }

    fn star(v: [f64; 6]) -> CatalogStar {
        CatalogStar {
            ra: Angle::from_radians(v[0]),
            dec: Angle::from_radians(v[1]),
            pm_ra_cos_dec_mas_per_year: v[2],
            pm_dec_mas_per_year: v[3],
            parallax_mas: v[4],
            radial_velocity_km_per_s: v[5],
        }
    }

    fn propagated(v: [f64; 6], to: TT) -> [f64; 6] {
        let s = star(v).propagate(&tt(J2016, 0.0), &to).unwrap();
        [
            s.ra.radians(),
            s.dec.radians(),
            s.pm_ra_cos_dec_mas_per_year,
            s.pm_dec_mas_per_year,
            s.parallax_mas,
            s.radial_velocity_km_per_s,
        ]
    }

    // Expected values in these tests are eraPmsafe's, in the units propagate converts to.
    #[test]
    fn test_barnards_star_over_a_century_matches_erfa() {
        assert_eq!(
            propagated(BARNARD, tt(J2016, DAYS_PER_JULIAN_CENTURY)),
            [
                4.702370963609202,
                0.08777316646287313,
                -811.88367199483,
                10491.423017780007,
                550.3761078308636,
                -110.04171924720609
            ]
        );
    }

    // The raised parallax comes back out, and the light-time terms give the star a small
    // radial velocity.
    #[test]
    fn test_zero_parallax_is_raised_as_erfa_raises_it() {
        assert_eq!(
            propagated([1.3, -0.4, 35.7, -12.4, 0.0, 0.0], tt(2415020.0, 0.0)),
            [
                1.2999782022633735,
                -0.39999302635480083,
                35.699894726818606,
                -12.40030302870577,
                0.059730277472676487,
                -0.06374762723009973
            ]
        );
    }

    #[test]
    fn test_motion_across_the_pole_matches_erfa() {
        let start = [0.7, 1.5707945814656445, 0.0, 10000.0, 0.0, 0.0];
        assert_eq!(
            propagated(start, tt(J2016, DAYS_PER_JULIAN_CENTURY)),
            [
                3.8415926535897933,
                1.5659499732960585,
                -1.328212764943143e-12,
                -9999.764966918956,
                15.805123566051114,
                14.540620592250239
            ]
        );
    }

    #[test]
    fn test_star_whose_correction_never_settles_matches_erfa() {
        let start = [
            1.9997587983405145,
            0.7258225557123114,
            -0.8030870173270971,
            0.5068112021302534,
            0.0,
            0.0,
        ];
        assert_eq!(
            propagated(start, tt(J2000_JD, 0.0)),
            [
                1.99975888162857,
                0.7258225163988689,
                -0.8030869893092908,
                0.5068112465269244,
                0.001500890834188184,
                -0.00022094311847780304
            ]
        );
    }

    #[test]
    fn test_non_finite_star_is_err() {
        let epochs = (tt(J2016, 0.0), tt(J2016, DAYS_PER_JULIAN_YEAR));
        for field in 0..6 {
            for bad in [f64::NAN, f64::INFINITY] {
                let mut v = BARNARD;
                v[field] = bad;
                assert!(star(v).propagate(&epochs.0, &epochs.1).is_err());
            }
        }
    }

    #[test]
    fn test_non_finite_epoch_is_err() {
        let s = star(BARNARD);
        assert!(s.propagate(&tt(f64::NAN, 0.0), &tt(J2016, 0.0)).is_err());
        assert!(s
            .propagate(&tt(J2016, 0.0), &tt(f64::INFINITY, 0.0))
            .is_err());
    }

    #[test]
    fn test_declination_beyond_the_pole_is_err() {
        let mut v = BARNARD;
        v[1] = 1.6;
        assert!(star(v)
            .propagate(&tt(J2016, 0.0), &tt(J2016, DAYS_PER_JULIAN_YEAR))
            .is_err());
    }

    #[test]
    fn test_speed_beyond_half_of_light_is_err() {
        let mut v = BARNARD;
        v[5] = 0.6 * SPEED_OF_LIGHT_M_PER_S / 1e3;
        assert!(star(v)
            .propagate(&tt(J2016, 0.0), &tt(J2016, DAYS_PER_JULIAN_YEAR))
            .is_err());
    }
}
