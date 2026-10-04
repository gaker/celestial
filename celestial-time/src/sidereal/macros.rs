// GMST and GAST are a sidereal angle at Greenwich, and LMST and LAST add the
// observer's location. The two in each pair differ only in how the Greenwich
// angle is computed, so these macros write everything else once.

macro_rules! sidereal_angle_accessors {
    ($field:tt) => {
        pub fn angle(&self) -> $crate::sidereal::angle::SiderealAngle {
            self.$field
        }

        pub fn hours(&self) -> f64 {
            self.$field.hours()
        }

        pub fn degrees(&self) -> f64 {
            self.$field.degrees()
        }

        pub fn radians(&self) -> f64 {
            self.$field.radians()
        }

        pub fn hour_angle_to_target(&self, target_ra_hours: f64) -> $crate::TimeResult<f64> {
            self.$field.hour_angle_to_target(target_ra_hours)
        }
    };
}

macro_rules! greenwich_sidereal_time {
    ($name:ident, $calculate:ident, $local:ident, $to_local:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        pub struct $name($crate::sidereal::angle::SiderealAngle);

        impl $name {
            pub fn from_ut1_and_tt(
                ut1: &$crate::scales::ut1::UT1,
                tt: &$crate::scales::tt::TT,
            ) -> $crate::TimeResult<Self> {
                $crate::sidereal::angle::SiderealAngle::from_radians($calculate(ut1, tt)?).map(Self)
            }

            pub fn from_hours(hours: f64) -> $crate::TimeResult<Self> {
                $crate::sidereal::angle::SiderealAngle::from_hours(hours).map(Self)
            }

            pub fn from_degrees(degrees: f64) -> $crate::TimeResult<Self> {
                $crate::sidereal::angle::SiderealAngle::from_degrees(degrees).map(Self)
            }

            pub fn from_radians(radians: f64) -> $crate::TimeResult<Self> {
                $crate::sidereal::angle::SiderealAngle::from_radians(radians).map(Self)
            }

            sidereal_angle_accessors!(0);

            pub fn $to_local(
                &self,
                location: &celestial_core::location::Location,
            ) -> $crate::TimeResult<$local> {
                let local =
                    celestial_core::angle::wrap_0_2pi(self.radians() + location.longitude())?;
                $local::from_radians(local, location)
            }
        }

        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, concat!(stringify!($name), " {}"), self.0)
            }
        }
    };
}

macro_rules! local_sidereal_time {
    ($name:ident, $greenwich:ident, $to_local:ident, $to_greenwich:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq)]
        #[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
        pub struct $name {
            angle: $crate::sidereal::angle::SiderealAngle,
            location: celestial_core::location::Location,
        }

        impl $name {
            pub fn from_ut1_tt_and_location(
                ut1: &$crate::scales::ut1::UT1,
                tt: &$crate::scales::tt::TT,
                location: &celestial_core::location::Location,
            ) -> $crate::TimeResult<Self> {
                $greenwich::from_ut1_and_tt(ut1, tt)?.$to_local(location)
            }

            pub fn from_hours(
                hours: f64,
                location: &celestial_core::location::Location,
            ) -> $crate::TimeResult<Self> {
                let angle = $crate::sidereal::angle::SiderealAngle::from_hours(hours)?;
                Ok(Self {
                    angle,
                    location: *location,
                })
            }

            pub fn from_degrees(
                degrees: f64,
                location: &celestial_core::location::Location,
            ) -> $crate::TimeResult<Self> {
                let angle = $crate::sidereal::angle::SiderealAngle::from_degrees(degrees)?;
                Ok(Self {
                    angle,
                    location: *location,
                })
            }

            pub fn from_radians(
                radians: f64,
                location: &celestial_core::location::Location,
            ) -> $crate::TimeResult<Self> {
                let angle = $crate::sidereal::angle::SiderealAngle::from_radians(radians)?;
                Ok(Self {
                    angle,
                    location: *location,
                })
            }

            pub fn location(&self) -> celestial_core::location::Location {
                self.location
            }

            sidereal_angle_accessors!(angle);

            pub fn $to_greenwich(&self) -> $crate::TimeResult<$greenwich> {
                let greenwich =
                    celestial_core::angle::wrap_0_2pi(self.radians() - self.location.longitude())?;
                $greenwich::from_radians(greenwich)
            }
        }

        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let lat_deg = self.location.latitude() * celestial_core::constants::RAD_TO_DEG;
                let lon_deg = self.location.longitude() * celestial_core::constants::RAD_TO_DEG;
                write!(
                    f,
                    concat!(stringify!($name), " {} at ({:.4}°, {:.4}°)"),
                    self.angle, lat_deg, lon_deg
                )
            }
        }
    };
}

#[cfg(test)]
mod tests {
    use crate::scales::tt::TT;
    use crate::scales::ut1::UT1;
    use crate::sidereal::angle::SiderealAngle;
    use crate::sidereal::gast::GAST;
    use crate::sidereal::gmst::GMST;
    use crate::sidereal::last::LAST;
    use crate::sidereal::lmst::LMST;
    use crate::TimeError;
    use celestial_core::constants::PI;
    use celestial_core::location::Location;

    fn not_finite(name: &str, shown: &str) -> TimeError {
        TimeError::ConversionError(format!("{} must be finite, got {}", name, shown))
    }

    macro_rules! check_greenwich_sidereal_time {
        ($name:ident) => {{
            let label = stringify!($name);
            assert_eq!($name::from_hours(f64::NAN), Err(not_finite("hours", "NaN")));
            assert_eq!(
                $name::from_degrees(f64::INFINITY),
                Err(not_finite("degrees", "inf"))
            );
            assert_eq!(
                $name::from_radians(-f64::INFINITY),
                Err(not_finite("radians", "-inf"))
            );

            let from_degrees = $name::from_degrees(180.0).unwrap();
            assert_eq!(
                from_degrees.angle(),
                SiderealAngle::from_hours(12.0).unwrap()
            );
            assert_eq!(
                (from_degrees.hours(), from_degrees.degrees()),
                (12.0, 180.0)
            );
            let from_radians = $name::from_radians(PI).unwrap();
            assert_eq!((from_radians.hours(), from_radians.radians()), (12.0, PI));

            assert_eq!(from_degrees.hour_angle_to_target(6.0), Ok(6.0));
            assert_eq!(
                from_degrees.hour_angle_to_target(f64::NAN),
                Err(not_finite("target_ra_hours", "NaN"))
            );
            assert_eq!(from_degrees.to_string(), format!("{label} 12.000000h"));
        }};
    }

    macro_rules! check_local_sidereal_time {
        ($name:ident, $greenwich:ident, $to_greenwich:ident) => {{
            let label = stringify!($name);
            let mauna_kea = Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap();
            assert_eq!(
                $name::from_hours(f64::NAN, &mauna_kea),
                Err(not_finite("hours", "NaN"))
            );
            assert_eq!(
                $name::from_degrees(f64::INFINITY, &mauna_kea),
                Err(not_finite("degrees", "inf"))
            );
            assert_eq!(
                $name::from_radians(-f64::INFINITY, &mauna_kea),
                Err(not_finite("radians", "-inf"))
            );

            let from_degrees = $name::from_degrees(180.0, &mauna_kea).unwrap();
            assert_eq!(from_degrees.location(), mauna_kea);
            assert_eq!(
                from_degrees.angle(),
                SiderealAngle::from_hours(12.0).unwrap()
            );
            assert_eq!(
                (from_degrees.hours(), from_degrees.degrees()),
                (12.0, 180.0)
            );
            let from_radians = $name::from_radians(PI, &mauna_kea).unwrap();
            assert_eq!((from_radians.hours(), from_radians.radians()), (12.0, PI));

            assert_eq!(from_degrees.hour_angle_to_target(6.0), Ok(6.0));
            assert_eq!(
                from_degrees.hour_angle_to_target(f64::NAN),
                Err(not_finite("target_ra_hours", "NaN"))
            );
            assert_eq!(
                from_degrees.to_string(),
                format!("{label} 12.000000h at (19.8283°, -155.4783°)")
            );

            let (ut1, tt) = (UT1::j2000(), TT::j2000());
            let greenwich = $greenwich::from_ut1_and_tt(&ut1, &tt).unwrap();
            let at_greenwich =
                $name::from_ut1_tt_and_location(&ut1, &tt, &Location::greenwich()).unwrap();
            assert_eq!(at_greenwich.angle(), greenwich.angle());
            assert_eq!(at_greenwich.$to_greenwich(), Ok(greenwich));
        }};
    }

    #[test]
    fn test_gmst_shared_api() {
        check_greenwich_sidereal_time!(GMST);
    }

    #[test]
    fn test_gast_shared_api() {
        check_greenwich_sidereal_time!(GAST);
    }

    #[test]
    fn test_lmst_shared_api() {
        check_local_sidereal_time!(LMST, GMST, to_gmst);
    }

    #[test]
    fn test_last_shared_api() {
        check_local_sidereal_time!(LAST, GAST, to_gast);
    }
}
