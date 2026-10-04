use super::check_epoch;
use crate::eop::record::EopParameters;
use crate::errors::{CoordError, CoordResult};
use crate::frames::cirs::CIRSPosition;
use crate::frames::itrs::ITRSPosition;
use crate::frames::tirs::TIRSPosition;
use crate::frames::topocentric::HourAnglePosition;
use celestial_core::angle::{wrap_0_2pi, wrap_pm_pi, Angle};
use celestial_core::constants::{ARCSEC_TO_RAD, MJD_ZERO_POINT};
use celestial_core::location::Location;
use celestial_core::matrix::RotationMatrix3;
use celestial_core::utils::jd_to_centuries;
use celestial_time::scales::conversions::utc_ut1::ToUT1WithDUT1;
use celestial_time::scales::conversions::{ToTAI, ToUTC};
use celestial_time::scales::tt::TT;
use celestial_time::transforms::rotation::earth_rotation_angle;

// UT1-UTC changes by up to a few milliseconds a day, so a record more than a
// day from the epoch is treated as missing rather than silently used.
const MAX_EOP_AGE_DAYS: f64 = 1.0;

// The Earth's orientation at one epoch: the Earth rotation angle that turns
// CIRS into TIRS, and the polar motion that turns TIRS into ITRS.
#[derive(Debug, Clone, PartialEq)]
pub struct EarthRotation {
    epoch: TT,
    era: f64,
    cirs_to_tirs: RotationMatrix3,
    tirs_to_itrs: RotationMatrix3,
}

impl EarthRotation {
    pub fn new(epoch: &TT, eop: &EopParameters) -> CoordResult<Self> {
        check_eop(epoch, eop)?;
        let era = earth_rotation_angle_at(epoch, eop.ut1_utc)?;
        let mut cirs_to_tirs = RotationMatrix3::identity();
        cirs_to_tirs.rotate_z(era);
        Ok(Self {
            epoch: *epoch,
            era,
            cirs_to_tirs,
            tirs_to_itrs: polar_motion(epoch, eop.x_p * ARCSEC_TO_RAD, eop.y_p * ARCSEC_TO_RAD),
        })
    }

    pub fn epoch(&self) -> TT {
        self.epoch
    }

    pub fn era(&self) -> f64 {
        self.era
    }

    pub fn cirs_to_tirs(&self, cirs: &CIRSPosition) -> CoordResult<TIRSPosition> {
        check_epoch(self.epoch, cirs.epoch())?;
        let p = match cirs.distance() {
            Some(distance) => cirs.unit_vector() * distance.au(),
            None => cirs.unit_vector(),
        };
        TIRSPosition::from_position_vector(self.cirs_to_tirs * p, self.epoch)
    }

    pub fn tirs_to_cirs(&self, tirs: &TIRSPosition) -> CoordResult<CIRSPosition> {
        check_epoch(self.epoch, tirs.epoch())?;
        let p = self.cirs_to_tirs.transpose() * tirs.position_vector();
        CIRSPosition::from_unit_vector(p, self.epoch)
    }

    pub fn tirs_to_itrs(&self, tirs: &TIRSPosition) -> CoordResult<ITRSPosition> {
        check_epoch(self.epoch, tirs.epoch())?;
        let p = self.tirs_to_itrs * tirs.position_vector();
        ITRSPosition::from_position_vector(p, self.epoch)
    }

    pub fn itrs_to_tirs(&self, itrs: &ITRSPosition) -> CoordResult<TIRSPosition> {
        check_epoch(self.epoch, itrs.epoch())?;
        let p = self.tirs_to_itrs.transpose() * itrs.position_vector();
        TIRSPosition::from_position_vector(p, self.epoch)
    }

    // The hour angle is measured on the CIO-based equator from the observer's
    // meridian, which sits at ERA + east longitude.
    pub fn cirs_to_hour_angle(
        &self,
        cirs: &CIRSPosition,
        observer: &Location,
    ) -> CoordResult<HourAnglePosition> {
        check_epoch(self.epoch, cirs.epoch())?;
        let ha = wrap_pm_pi(self.era + observer.longitude() - cirs.ra().radians())?;
        let mut position =
            HourAnglePosition::new(Angle::from_radians(ha), cirs.dec(), *observer, self.epoch)?;
        if let Some(distance) = cirs.distance() {
            position.set_distance(distance);
        }
        Ok(position)
    }

    pub fn hour_angle_to_cirs(&self, position: &HourAnglePosition) -> CoordResult<CIRSPosition> {
        check_epoch(self.epoch, position.epoch())?;
        let longitude = position.observer().longitude();
        let ra = wrap_0_2pi(self.era + longitude - position.hour_angle().radians())?;
        let mut cirs =
            CIRSPosition::new(Angle::from_radians(ra), position.declination(), self.epoch)?;
        if let Some(distance) = position.distance() {
            cirs.set_distance(distance);
        }
        Ok(cirs)
    }
}

fn check_eop(epoch: &TT, eop: &EopParameters) -> CoordResult<()> {
    let age = epoch.to_julian_date().to_f64() - (eop.mjd + MJD_ZERO_POINT);
    if age.is_nan() || libm::fabs(age) > MAX_EOP_AGE_DAYS {
        return Err(CoordError::data_unavailable(format!(
            "EOP record at MJD {} is more than {MAX_EOP_AGE_DAYS} day from the epoch",
            eop.mjd
        )));
    }
    for (name, value) in [("UT1-UTC", eop.ut1_utc), ("x_p", eop.x_p), ("y_p", eop.y_p)] {
        if !value.is_finite() {
            return Err(CoordError::invalid_coordinate(format!(
                "EOP {name} {value} is not finite"
            )));
        }
    }
    Ok(())
}

fn earth_rotation_angle_at(epoch: &TT, ut1_utc: f64) -> CoordResult<f64> {
    let ut1 = epoch.to_tai()?.to_utc()?.to_ut1_with_dut1(ut1_utc)?;
    Ok(earth_rotation_angle(&ut1.to_julian_date())?)
}

// W = R1(-yp) R2(-xp) R3(s'), with the TIO locator s' from its IERS 2003
// linear model, -47 microarcseconds per Julian century of TT.
fn polar_motion(epoch: &TT, xp: f64, yp: f64) -> RotationMatrix3 {
    let jd = epoch.to_julian_date();
    let s_prime = -47e-6 * jd_to_centuries(jd.jd1(), jd.jd2()) * ARCSEC_TO_RAD;
    let mut w = RotationMatrix3::identity();
    w.rotate_z(s_prime);
    w.rotate_y(-xp);
    w.rotate_x(-yp);
    w
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eop::record::EopRecord;
    use celestial_core::constants::{J2000_JD, MJD_ZERO_POINT};
    use celestial_core::matrix::Vector3;

    fn eop_at_j2000(ut1_utc: f64) -> EopParameters {
        EopRecord::new(J2000_JD - MJD_ZERO_POINT, 0.0, 0.0, ut1_utc)
            .unwrap()
            .to_parameters()
    }

    #[test]
    fn era_matches_erfa_at_j2000() {
        // eraTttai, eraTaiutc, eraUtcut1 with DUT1 = 0.3 s, then eraEra00.
        let rotation = EarthRotation::new(&TT::j2000(), &eop_at_j2000(0.3)).unwrap();
        assert_eq!(rotation.era(), 4.890302717983435);
    }

    #[test]
    fn cirs_origin_turns_westward_by_era() {
        let epoch = TT::j2000();
        let rotation = EarthRotation::new(&epoch, &eop_at_j2000(0.3)).unwrap();
        let origin =
            CIRSPosition::new(Angle::from_radians(0.0), Angle::from_radians(0.0), epoch).unwrap();

        let tirs = rotation.cirs_to_tirs(&origin).unwrap();

        let (s, c) = libm::sincos(rotation.era());
        assert_eq!(tirs.position_vector(), Vector3::new(c, -s, 0.0));
    }

    #[test]
    fn stale_eop_is_rejected() {
        let eop = EopRecord::new(51542.0, 0.0, 0.0, 0.3)
            .unwrap()
            .to_parameters();
        assert!(EarthRotation::new(&TT::j2000(), &eop).is_err());
    }

    #[test]
    fn position_at_another_epoch_is_rejected() {
        let rotation = EarthRotation::new(&TT::j2000(), &eop_at_j2000(0.3)).unwrap();
        let later = TT::from_julian_date(celestial_time::julian::JulianDate::new(J2000_JD, 0.25));
        let tirs = TIRSPosition::new(1.0, 0.0, 0.0, later).unwrap();
        assert!(rotation.tirs_to_cirs(&tirs).is_err());
        assert!(rotation.tirs_to_itrs(&tirs).is_err());
    }
}
