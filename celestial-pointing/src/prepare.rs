//! Pre-fit helpers for preparing observation coordinates.
//!
//! The INDAT format expects telescope coordinates to be refraction-removed
//! before being written to a dat file. This module exposes the helpers
//! callers need to do that preparation — they are NOT auto-invoked inside
//! the fit. The fit assumes the inputs already meet the format contract.
//!
//! See `book/src/pointing/indat-format.md` for the full input-coordinate
//! contract.

use crate::observation::SiteParams;
use celestial_core::{Angle, Location};
use celestial_coords::frames::HourAnglePosition;
use celestial_time::TT;

/// Adds atmospheric refraction to a sky direction.
///
/// Takes a true topocentric direction (where a star actually is in the
/// sky) and returns the refracted direction an observer would see through
/// the atmosphere.
///
/// Per the TPOINT pipeline, catalog coordinates step through this
/// transform on the way to Observed Place, where they can be compared
/// against telescope encoder readings.
///
/// When the site has no atmospheric data (`pressure <= 0.0`) the inputs
/// are returned unchanged, matching the TPOINT convention. The inverse
/// transform (`refraction_shift` with `apply=false`) is used by
/// round-trip tests within this module.
pub(crate) fn apply_refraction(
    ra: Angle,
    dec: Angle,
    lst: Angle,
    location: &Location,
    site: &SiteParams,
) -> (Angle, Angle) {
    refraction_shift(ra, dec, lst, location, site, true)
}

fn refraction_shift(
    ra: Angle,
    dec: Angle,
    lst: Angle,
    location: &Location,
    site: &SiteParams,
    apply: bool,
) -> (Angle, Angle) {
    if site.pressure <= 0.0 {
        return (ra, dec);
    }
    let ha = (lst - ra).wrapped();
    let epoch = TT::j2000();
    let Ok(hourangle) = HourAnglePosition::new(ha, dec, *location, epoch) else {
        return (ra, dec);
    };
    let Ok(topo) = hourangle.to_topocentric() else {
        return (ra, dec);
    };
    let shifted = if apply {
        topo.with_refraction(
            site.pressure,
            site.temperature,
            site.humidity,
            site.wavelength,
        )
    } else {
        topo.without_refraction(
            site.pressure,
            site.temperature,
            site.humidity,
            site.wavelength,
        )
    };
    let Ok(shifted_ha) = shifted.to_hour_angle() else {
        return (ra, dec);
    };
    let shifted_ra = (lst - shifted_ha.hour_angle()).wrapped();
    (shifted_ra, shifted_ha.declination())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn standard_site() -> SiteParams {
        SiteParams {
            latitude: Angle::from_degrees(39.0),
            longitude: Angle::from_degrees(0.0),
            temperature: 10.0,
            pressure: 1013.25,
            elevation: 0.0,
            humidity: 0.5,
            wavelength: 0.55,
            lapse_rate: 0.0065,
        }
    }

    fn location_at(lat_deg: f64) -> Location {
        Location::from_degrees(lat_deg, 0.0, 0.0).unwrap()
    }

    fn arcsec(a: Angle, b: Angle) -> f64 {
        (a - b).wrapped().arcseconds().abs()
    }

    #[test]
    fn zenith_correction_is_near_zero() {
        let lat = 39.0;
        let lst = Angle::from_hours(0.0);
        let ra = lst;
        let dec = Angle::from_degrees(lat);
        let (ra2, dec2) = refraction_shift(ra, dec, lst, &location_at(lat), &standard_site(), false);
        assert!(arcsec(ra, ra2) < 1e-6);
        assert!(arcsec(dec, dec2) < 1e-6);
    }

    #[test]
    fn low_altitude_correction_is_on_order_of_hundred_arcsec() {
        let lat = 39.0;
        let lst = Angle::from_hours(0.0);
        let ra = lst;
        let dec = Angle::from_degrees(lat - 60.0);
        let (ra2, dec2) = refraction_shift(ra, dec, lst, &location_at(lat), &standard_site(), false);

        let epoch = TT::j2000();
        let apparent = HourAnglePosition::new(
            (lst - ra).wrapped(),
            dec,
            location_at(lat),
            epoch,
        )
        .unwrap()
        .to_topocentric()
        .unwrap();
        let true_hap = HourAnglePosition::new(
            (lst - ra2).wrapped(),
            dec2,
            location_at(lat),
            epoch,
        )
        .unwrap();
        let true_topo = true_hap.to_topocentric().unwrap();
        let elev_delta = (apparent.elevation() - true_topo.elevation())
            .wrapped()
            .arcseconds();
        assert!(
            (70.0..160.0).contains(&elev_delta),
            "expected ~100 arcsec refraction at 30 deg alt, got {}",
            elev_delta
        );
        assert!(true_topo.elevation() < apparent.elevation());
    }

    #[test]
    fn zero_pressure_disables_correction() {
        let mut site = standard_site();
        site.pressure = 0.0;
        let lat = 39.0;
        let lst = Angle::from_hours(0.0);
        let ra = lst;
        let dec = Angle::from_degrees(lat - 60.0);
        let (ra2, dec2) = refraction_shift(ra, dec, lst, &location_at(lat), &site, false);
        assert_eq!(ra, ra2);
        assert_eq!(dec, dec2);
    }

    #[test]
    fn negative_pressure_disables_correction() {
        let mut site = standard_site();
        site.pressure = -1.0;
        let lat = 39.0;
        let lst = Angle::from_hours(0.0);
        let ra = lst;
        let dec = Angle::from_degrees(lat - 60.0);
        let (ra2, dec2) = refraction_shift(ra, dec, lst, &location_at(lat), &site, false);
        assert_eq!(ra, ra2);
        assert_eq!(dec, dec2);
    }

    #[test]
    fn higher_altitude_smaller_correction_than_lower() {
        let lat = 39.0;
        let lst = Angle::from_hours(0.0);
        let site = standard_site();
        let loc = location_at(lat);

        let ra = lst;
        let dec_high = Angle::from_degrees(lat - 10.0);
        let (_, dec_high_true) = refraction_shift(ra, dec_high, lst, &loc, &site, false);
        let shift_high = arcsec(dec_high, dec_high_true);

        let dec_low = Angle::from_degrees(lat - 60.0);
        let (_, dec_low_true) = refraction_shift(ra, dec_low, lst, &loc, &site, false);
        let shift_low = arcsec(dec_low, dec_low_true);

        assert!(
            shift_low > shift_high,
            "low altitude shift {} should exceed high altitude shift {}",
            shift_low,
            shift_high
        );
    }

    #[test]
    fn apply_and_remove_are_approximately_inverses() {
        // apply() computes refraction at the true elevation; remove() computes
        // it at the apparent (refracted) elevation. They are not exact
        // algebraic inverses because the zenith distance differs by the
        // refraction amount itself. Round-trip error should be a small
        // fraction of an arcsecond at moderate altitude.
        let lat = 39.0;
        let lst = Angle::from_hours(0.0);
        let ra = lst;
        let dec = Angle::from_degrees(lat - 30.0);
        let loc = location_at(lat);
        let site = standard_site();
        let (ra2, dec2) = apply_refraction(ra, dec, lst, &loc, &site);
        let (ra3, dec3) = refraction_shift(ra2, dec2, lst, &loc, &site, false);
        assert!(arcsec(ra, ra3) < 0.5, "ra drift {}", arcsec(ra, ra3));
        assert!(arcsec(dec, dec3) < 0.5, "dec drift {}", arcsec(dec, dec3));
    }
}
