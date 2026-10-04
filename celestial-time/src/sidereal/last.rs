use super::gast::GAST;
use super::gmst::GMST;
use super::lmst::LMST;
use crate::julian::JulianDate;
use crate::scales::tt::TT;
use crate::scales::ut1::UT1;
use crate::TimeResult;
use celestial_core::angle::{wrap_0_2pi, wrap_pm_pi};

local_sidereal_time!(LAST, GAST, to_last, to_gast);

impl LAST {
    pub fn to_lmst(&self, tt: &TT) -> TimeResult<LMST> {
        let lmst = wrap_0_2pi(self.radians() - equation_of_equinoxes(tt)?)?;
        LMST::from_radians(lmst, &self.location)
    }

    pub fn to_gmst(&self, tt: &TT) -> TimeResult<GMST> {
        self.to_lmst(tt)?.to_gmst()
    }
}

// GAST - GMST: the equation of the equinoxes that matches how LAST is built,
// complementary terms included. UT1 cancels in the difference; (0, 0) is the
// choice ERFA's ee06a makes.
fn equation_of_equinoxes(tt: &TT) -> TimeResult<f64> {
    let ut1 = UT1::from_julian_date(JulianDate::new(0.0, 0.0));
    let gast = GAST::from_ut1_and_tt(&ut1, tt)?;
    let gmst = GMST::from_ut1_and_tt(&ut1, tt)?;
    Ok(wrap_pm_pi(gast.radians() - gmst.radians())?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::constants::MJD_ZERO_POINT;
    use celestial_core::location::Location;

    #[test]
    fn test_last_to_lmst_conversion() {
        // anp(anp(eraGst06a + elong) - eraEe06a). The equation of the equinoxes
        // must include the complementary terms (up to ~2.6 mas) to get these.
        let cases = [
            (51544.5, 51544.50079861111, 0.0, 4.894961283639734),
            (60000.25, 60000.25079861111, -2.7144, 1.5590097514538646),
            (40000.75, 40000.75079861111, 1.2345, 3.894275327009728),
            (73000.125, 73000.12579861112, 3.0, 3.9274851400914272),
        ];
        for (ut1_mjd, tt_mjd, longitude, expected) in cases {
            let location = Location::new(0.5, longitude, 0.0).unwrap();
            let ut1 = UT1::from_julian_date(JulianDate::new(MJD_ZERO_POINT, ut1_mjd));
            let tt = TT::from_julian_date(JulianDate::new(MJD_ZERO_POINT, tt_mjd));
            let last = LAST::from_ut1_tt_and_location(&ut1, &tt, &location).unwrap();
            let lmst = last.to_lmst(&tt).unwrap();
            assert_eq!(lmst.radians(), expected, "UT1 MJD {ut1_mjd}");
        }
    }

    #[test]
    fn test_last_to_gast_matches_erfa_anp() {
        // anp(LAST - elong), all in radians.
        for (last, longitude, expected) in [(1.0, 2.9, 4.383185307179586), (2.5, -3.1, 5.6)] {
            let location = Location::new(0.5, longitude, 0.0).unwrap();
            let gast = LAST::from_radians(last, &location)
                .unwrap()
                .to_gast()
                .unwrap();
            assert_eq!(gast.radians(), expected, "LAST {last}, elong {longitude}");
        }
    }

    #[test]
    fn test_last_to_gmst_matches_erfa() {
        // anp(anp(LAST - eraEe06a) - elong), all in radians.
        let tt = TT::from_julian_date(JulianDate::new(MJD_ZERO_POINT, 51544.50079861111));
        let cases = [
            (4.894961283639734, -2.7136, 1.3254379368685765),
            (1.0, 2.9, 4.383247267588015),
        ];
        for (last, longitude, expected) in cases {
            let location = Location::new(0.5, longitude, 0.0).unwrap();
            let gmst = LAST::from_radians(last, &location)
                .unwrap()
                .to_gmst(&tt)
                .unwrap();
            assert_eq!(gmst.radians(), expected, "LAST {last}, elong {longitude}");
        }
    }
}
