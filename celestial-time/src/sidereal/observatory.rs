use super::{GAST, GMST, LAST, LMST};
use crate::scales::{TT, UT1};
use crate::TimeResult;
use celestial_core::Location;

#[derive(Debug, Clone, Copy)]
pub struct ObservatoryContext<'a> {
    ut1: &'a UT1,
    tt: &'a TT,
    location: &'a Location,
}

impl<'a> ObservatoryContext<'a> {
    pub fn new(ut1: &'a UT1, tt: &'a TT, location: &'a Location) -> Self {
        Self { ut1, tt, location }
    }

    /// Get location for a famous observatory
    ///
    /// Provides quick access to well-known observatory locations.
    /// Returns an owned Location that can be used with `new()`.
    ///
    /// # Supported Observatories
    /// * `"mauna_kea"` - Mauna Kea Observatory, Hawaii
    /// * `"greenwich"` - Royal Observatory Greenwich, UK
    /// * `"palomar"` - Palomar Observatory, California
    /// * `"vlt"` - Very Large Telescope, Chile
    /// * `"keck"` - W. M. Keck Observatory, Hawaii
    ///
    /// # Examples
    /// ```
    /// use celestial_time::{UT1, TT};
    /// use celestial_time::sidereal::ObservatoryContext;
    ///
    /// let ut1 = UT1::j2000();
    /// let tt = TT::j2000();
    /// let location = ObservatoryContext::observatory_location("mauna_kea").unwrap();
    /// let observatory = ObservatoryContext::new(&ut1, &tt, &location);
    /// ```
    pub fn observatory_location(observatory_name: &str) -> TimeResult<Location> {
        match observatory_name {
            "mauna_kea" | "keck" => Ok(Location::from_degrees(19.8283, -155.4783, 4145.0)
                .expect("Keck Observatory coordinates are valid")),
            "greenwich" => Ok(Location::greenwich()),
            "palomar" => Ok(Location::from_degrees(33.3563, -116.8650, 1712.0)
                .expect("Palomar coordinates are valid")),
            "vlt" => Ok(Location::from_degrees(-24.6275, -70.4044, 2635.0)
                .expect("VLT coordinates are valid")),
            _ => Err(crate::TimeError::CalculationError(format!(
                "Unknown observatory: {}",
                observatory_name
            ))),
        }
    }

    /// Get the UT1 time
    pub fn ut1(&self) -> &UT1 {
        self.ut1
    }

    /// Get the TT time
    pub fn tt(&self) -> &TT {
        self.tt
    }

    /// Get the observer location
    pub fn location(&self) -> &Location {
        self.location
    }

    /// Calculate Greenwich Mean Sidereal Time (GMST)
    ///
    /// # Examples
    /// ```
    /// use celestial_time::{UT1, TT};
    /// use celestial_time::sidereal::ObservatoryContext;
    /// use celestial_core::Location;
    ///
    /// let ut1 = UT1::j2000();
    /// let tt = TT::j2000();
    /// let location = Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap();
    /// let observatory = ObservatoryContext::new(&ut1, &tt, &location);
    /// let gmst = observatory.gmst().unwrap();
    /// ```
    pub fn gmst(&self) -> TimeResult<GMST> {
        GMST::from_ut1_and_tt(self.ut1, self.tt)
    }

    /// Calculate Greenwich Apparent Sidereal Time (GAST)
    ///
    /// # Examples
    /// ```
    /// use celestial_time::{UT1, TT};
    /// use celestial_time::sidereal::ObservatoryContext;
    /// use celestial_core::Location;
    ///
    /// let ut1 = UT1::j2000();
    /// let tt = TT::j2000();
    /// let location = Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap();
    /// let observatory = ObservatoryContext::new(&ut1, &tt, &location);
    /// let gast = observatory.gast().unwrap();
    /// ```
    pub fn gast(&self) -> TimeResult<GAST> {
        GAST::from_ut1_and_tt(self.ut1, self.tt)
    }

    /// Calculate Local Mean Sidereal Time (LMST)
    ///
    /// # Examples
    /// ```
    /// use celestial_time::{UT1, TT};
    /// use celestial_time::sidereal::ObservatoryContext;
    /// use celestial_core::Location;
    ///
    /// let ut1 = UT1::j2000();
    /// let tt = TT::j2000();
    /// let location = Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap();
    /// let observatory = ObservatoryContext::new(&ut1, &tt, &location);
    /// let lmst = observatory.lmst().unwrap();
    /// ```
    pub fn lmst(&self) -> TimeResult<LMST> {
        LMST::from_ut1_tt_and_location(self.ut1, self.tt, self.location)
    }

    /// Calculate Local Apparent Sidereal Time (LAST)
    ///
    /// # Examples
    /// ```
    /// use celestial_time::{UT1, TT};
    /// use celestial_time::sidereal::ObservatoryContext;
    /// use celestial_core::Location;
    ///
    /// let ut1 = UT1::j2000();
    /// let tt = TT::j2000();
    /// let location = Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap();
    /// let observatory = ObservatoryContext::new(&ut1, &tt, &location);
    /// let last = observatory.last().unwrap();
    /// ```
    pub fn last(&self) -> TimeResult<LAST> {
        LAST::from_ut1_tt_and_location(self.ut1, self.tt, self.location)
    }

    /// Get all sidereal times at once
    ///
    /// Returns a tuple of (GMST, GAST, LMST, LAST) for convenience.
    ///
    /// # Examples
    /// ```
    /// use celestial_time::{UT1, TT};
    /// use celestial_time::sidereal::ObservatoryContext;
    /// use celestial_core::Location;
    ///
    /// let ut1 = UT1::j2000();
    /// let tt = TT::j2000();
    /// let location = Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap();
    /// let observatory = ObservatoryContext::new(&ut1, &tt, &location);
    /// let (gmst, gast, lmst, last) = observatory.all_sidereal_times().unwrap();
    /// ```
    pub fn all_sidereal_times(&self) -> TimeResult<(GMST, GAST, LMST, LAST)> {
        let gmst = self.gmst()?;
        let gast = self.gast()?;
        let lmst = self.lmst()?;
        let last = self.last()?;
        Ok((gmst, gast, lmst, last))
    }

    /// Calculate hour angle to a target right ascension
    ///
    /// Uses Local Apparent Sidereal Time for the most accurate hour angle calculation.
    ///
    /// # Arguments
    /// * `target_ra_hours` - Target right ascension in hours
    ///
    /// # Returns
    /// Hour angle in hours, normalized to [-12, 12) range
    ///
    /// # Examples
    /// ```
    /// use celestial_time::{UT1, TT};
    /// use celestial_time::sidereal::ObservatoryContext;
    /// use celestial_core::Location;
    ///
    /// let ut1 = UT1::j2000();
    /// let tt = TT::j2000();
    /// let location = Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap();
    /// let observatory = ObservatoryContext::new(&ut1, &tt, &location);
    /// let hour_angle = observatory.hour_angle_to_target(6.0).unwrap(); // RA = 6h
    /// ```
    pub fn hour_angle_to_target(&self, target_ra_hours: f64) -> TimeResult<f64> {
        let last = self.last()?;
        Ok(last.hour_angle_to_target(target_ra_hours))
    }

    /// Get observatory information as a formatted string
    ///
    /// Useful for logging and debugging.
    ///
    /// # Examples
    /// ```
    /// use celestial_time::{UT1, TT};
    /// use celestial_time::sidereal::ObservatoryContext;
    /// use celestial_core::Location;
    ///
    /// let ut1 = UT1::j2000();
    /// let tt = TT::j2000();
    /// let location = Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap();
    /// let observatory = ObservatoryContext::new(&ut1, &tt, &location);
    /// println!("{}", observatory.info());
    /// ```
    pub fn info(&self) -> String {
        let lat_deg = self.location.latitude.to_degrees();
        let lon_deg = self.location.longitude.to_degrees();
        let height_m = self.location.height;

        format!(
            "Observatory at ({:.4}°, {:.4}°, {:.0}m) - UT1: {}, TT: {}",
            lat_deg,
            lon_deg,
            height_m,
            self.ut1.to_julian_date().jd1() + self.ut1.to_julian_date().jd2(),
            self.tt.to_julian_date().jd1() + self.tt.to_julian_date().jd2()
        )
    }
}

impl<'a> std::fmt::Display for ObservatoryContext<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.info())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mauna_kea() -> Location {
        Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap()
    }

    fn greenwich() -> Location {
        Location::greenwich()
    }

    #[test]
    fn test_observatory_context_creation() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();
        let observatory = ObservatoryContext::new(&ut1, &tt, &location);

        let obs_ut1_jd = observatory.ut1().to_julian_date();
        let ut1_jd = ut1.to_julian_date();
        let obs_tt_jd = observatory.tt().to_julian_date();
        let tt_jd = tt.to_julian_date();
        assert_eq!(
            obs_ut1_jd.jd1() + obs_ut1_jd.jd2(),
            ut1_jd.jd1() + ut1_jd.jd2()
        );
        assert_eq!(obs_tt_jd.jd1() + obs_tt_jd.jd2(), tt_jd.jd1() + tt_jd.jd2());
        assert_eq!(observatory.location().latitude, location.latitude);
        assert_eq!(observatory.location().longitude, location.longitude);
        assert_eq!(observatory.location().height, location.height);
    }

    #[test]
    fn test_all_sidereal_times() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();
        let observatory = ObservatoryContext::new(&ut1, &tt, &location);

        let (gmst, gast, lmst, last) = observatory.all_sidereal_times().unwrap();

        assert!(gmst.hours() >= 0.0 && gmst.hours() < 24.0);
        assert!(gast.hours() >= 0.0 && gast.hours() < 24.0);
        assert!(lmst.hours() >= 0.0 && lmst.hours() < 24.0);
        assert!(last.hours() >= 0.0 && last.hours() < 24.0);

        assert!((gast.hours() - gmst.hours()).abs() < 1.0);

        let expected_offset = -155.4783 / 15.0;
        let actual_offset = lmst.hours() - gmst.hours();
        let normalized_offset = if actual_offset > 12.0 {
            actual_offset - 24.0
        } else if actual_offset < -12.0 {
            actual_offset + 24.0
        } else {
            actual_offset
        };
        assert!((normalized_offset - expected_offset).abs() < 1e-10);

        let last_offset = last.hours() - gast.hours();
        let normalized_last_offset = if last_offset > 12.0 {
            last_offset - 24.0
        } else if last_offset < -12.0 {
            last_offset + 24.0
        } else {
            last_offset
        };
        assert!((normalized_last_offset - expected_offset).abs() < 1e-10);
    }

    #[test]
    fn test_hour_angle_calculation() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();
        let observatory = ObservatoryContext::new(&ut1, &tt, &location);

        let target_ra = 6.0;
        let hour_angle = observatory.hour_angle_to_target(target_ra).unwrap();

        let last = observatory.last().unwrap();
        let expected_ha = last.hour_angle_to_target(target_ra);
        assert!((hour_angle - expected_ha).abs() < 1e-12);
    }

    #[test]
    fn test_famous_observatories() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();

        let observatories = ["mauna_kea", "greenwich", "palomar", "vlt", "keck"];

        for name in observatories {
            let location = ObservatoryContext::observatory_location(name).unwrap();
            let observatory = ObservatoryContext::new(&ut1, &tt, &location);
            let (gmst, gast, lmst, last) = observatory.all_sidereal_times().unwrap();

            assert!(
                gmst.hours() >= 0.0 && gmst.hours() < 24.0,
                "Invalid GMST for {}",
                name
            );
            assert!(
                gast.hours() >= 0.0 && gast.hours() < 24.0,
                "Invalid GAST for {}",
                name
            );
            assert!(
                lmst.hours() >= 0.0 && lmst.hours() < 24.0,
                "Invalid LMST for {}",
                name
            );
            assert!(
                last.hours() >= 0.0 && last.hours() < 24.0,
                "Invalid LAST for {}",
                name
            );
        }
    }

    #[test]
    fn test_unknown_observatory() {
        let result = ObservatoryContext::observatory_location("unknown_observatory");
        assert!(result.is_err());
    }

    #[test]
    fn test_individual_sidereal_calculations() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();
        let observatory = ObservatoryContext::new(&ut1, &tt, &location);

        let (gmst_batch, gast_batch, lmst_batch, last_batch) =
            observatory.all_sidereal_times().unwrap();

        let gmst_individual = observatory.gmst().unwrap();
        let gast_individual = observatory.gast().unwrap();
        let lmst_individual = observatory.lmst().unwrap();
        let last_individual = observatory.last().unwrap();

        assert!((gmst_individual.hours() - gmst_batch.hours()).abs() < 1e-12);
        assert!((gast_individual.hours() - gast_batch.hours()).abs() < 1e-12);
        assert!((lmst_individual.hours() - lmst_batch.hours()).abs() < 1e-12);
        assert!((last_individual.hours() - last_batch.hours()).abs() < 1e-12);
    }

    #[test]
    fn test_display_and_info() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = mauna_kea();
        let observatory = ObservatoryContext::new(&ut1, &tt, &location);

        let info = observatory.info();
        let display = format!("{}", observatory);

        assert!(info.contains("19.8283"));
        assert!(info.contains("-155.4783"));
        assert!(info.contains("4145"));

        assert_eq!(info, display);
    }

    #[test]
    fn test_observatory_context_hour_angle() {
        let ut1 = UT1::j2000();
        let tt = TT::j2000();
        let location = greenwich();
        let observatory = ObservatoryContext::new(&ut1, &tt, &location);

        let target_ra = 12.0;
        let hour_angle = observatory.hour_angle_to_target(target_ra).unwrap();

        let last = observatory.last().unwrap();
        let expected_ha = last.hour_angle_to_target(target_ra);
        assert_eq!(hour_angle, expected_ha);

        assert!((-12.0..12.0).contains(&hour_angle));
    }
}
