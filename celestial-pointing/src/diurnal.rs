//! Diurnal aberration correction for telescope coordinates.
//!
//! Default convention assumes observations have had diurnal
//! aberration removed before being written. The INDAT option `:NODA`
//! signals that this step was skipped. When `:NODA` is absent, the
//! fitter is expected to apply the correction itself.

use celestial_core::{Angle, Location, Vector3};
use celestial_coords::aberration::{apply_aberration, remove_aberration};

const EARTH_ROTATION_RATE_RAD_SEC: f64 = 7.292_115_0e-5;
const SECONDS_PER_DAY: f64 = 86_400.0;
const AU_KM: f64 = celestial_core::constants::AU_KM;

/// Adds diurnal aberration to a sky direction.
///
/// Takes a geocentric-apparent direction (as from a catalog, already
/// precessed/nutated to date) and returns the topocentric direction a
/// rotating observer sees — i.e., the direction shifted by the observer's
/// diurnal velocity.
///
/// The inverse transform (`apply_or_remove` with `apply=false`) is used by
/// round-trip tests within this module.
pub(crate) fn apply_diurnal(
    cat_ra: Angle,
    cat_dec: Angle,
    lst: Angle,
    latitude: Angle,
    longitude: Angle,
    height_m: f64,
) -> (Angle, Angle) {
    apply_or_remove(cat_ra, cat_dec, lst, latitude, longitude, height_m, true)
}

fn apply_or_remove(
    ra: Angle,
    dec: Angle,
    lst: Angle,
    latitude: Angle,
    longitude: Angle,
    height_m: f64,
    apply: bool,
) -> (Angle, Angle) {
    let loc = match Location::new(latitude.radians(), longitude.radians(), height_m) {
        Ok(l) => l,
        Err(_) => return (ra, dec),
    };
    let (u_km, _v_km) = match loc.to_geocentric_km() {
        Ok(uv) => uv,
        Err(_) => return (ra, dec),
    };

    let velocity = diurnal_velocity_au_day(u_km, lst);
    let direction = spherical_to_cartesian(ra, dec);
    let corrected = if apply {
        apply_aberration(direction, velocity, 1.0)
    } else {
        remove_aberration(direction, velocity, 1.0)
    };
    cartesian_to_spherical(corrected)
}

fn diurnal_velocity_au_day(u_km: f64, lst: Angle) -> Vector3 {
    let speed_km_sec = EARTH_ROTATION_RATE_RAD_SEC * u_km;
    let speed_au_day = speed_km_sec * SECONDS_PER_DAY / AU_KM;
    let (sin_lst, cos_lst) = libm::sincos(lst.radians());
    Vector3::new(-speed_au_day * sin_lst, speed_au_day * cos_lst, 0.0)
}

fn spherical_to_cartesian(ra: Angle, dec: Angle) -> Vector3 {
    let (sin_ra, cos_ra) = libm::sincos(ra.radians());
    let (sin_dec, cos_dec) = libm::sincos(dec.radians());
    Vector3::new(cos_dec * cos_ra, cos_dec * sin_ra, sin_dec)
}

fn cartesian_to_spherical(v: Vector3) -> (Angle, Angle) {
    let ra = libm::atan2(v.y, v.x);
    let dec = libm::asin(v.z / v.magnitude());
    let ra_wrapped = if ra < 0.0 { ra + 2.0 * libm::atan(1.0) * 4.0 } else { ra };
    (Angle::from_radians(ra_wrapped), Angle::from_radians(dec))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn arcsec(a: Angle, b: Angle) -> f64 {
        (a - b).wrapped().arcseconds().abs()
    }

    #[test]
    fn zero_correction_at_pole() {
        let ra = Angle::from_hours(5.0);
        let dec = Angle::from_degrees(89.999);
        let lst = Angle::from_hours(12.0);
        let lat = Angle::from_degrees(90.0 - 1e-6);
        let (ra2, dec2) = apply_or_remove(ra, dec, lst, lat, Angle::from_radians(0.0), 0.0, false);
        assert!(arcsec(ra, ra2) < 1e-3);
        assert!(arcsec(dec, dec2) < 1e-3);
    }

    #[test]
    fn maximum_correction_on_meridian_equator() {
        // At LST=0, observer velocity points toward RA=6h (east of meridian).
        // A star on the meridian (RA=0h, dec=0) is perpendicular to velocity
        // and therefore sees maximum aberration.
        let lat = Angle::from_degrees(0.0);
        let lst = Angle::from_hours(0.0);
        let tel_ra = Angle::from_hours(0.0);
        let tel_dec = Angle::from_degrees(0.0);
        let (ra2, _dec2) = apply_or_remove(tel_ra, tel_dec, lst, lat, Angle::from_radians(0.0), 0.0, false);
        let shift = arcsec(tel_ra, ra2);
        assert!(
            shift > 0.25 && shift < 0.40,
            "expected ~0.32 arcsec, got {}",
            shift
        );
    }

    #[test]
    fn correction_smaller_at_mid_latitude_than_equator() {
        // Same meridian star, evaluated at two latitudes.
        let lst = Angle::from_hours(0.0);
        let tel_ra = Angle::from_hours(0.0);
        let tel_dec = Angle::from_degrees(0.0);
        let (ra_eq, _) = apply_or_remove(
            tel_ra,
            tel_dec,
            lst,
            Angle::from_degrees(0.0),
            Angle::from_radians(0.0),
            0.0,
            false,
        );
        let (ra_mid, _) = apply_or_remove(
            tel_ra,
            tel_dec,
            lst,
            Angle::from_degrees(60.0),
            Angle::from_radians(0.0),
            0.0,
            false,
        );
        assert!(arcsec(tel_ra, ra_eq) > arcsec(tel_ra, ra_mid));
    }

    #[test]
    fn apply_and_remove_are_inverses() {
        let lat = Angle::from_degrees(39.0);
        let lst = Angle::from_hours(3.0);
        let ra = Angle::from_hours(7.5);
        let dec = Angle::from_degrees(30.0);
        let lon = Angle::from_radians(0.0);
        let (ra2, dec2) = apply_diurnal(ra, dec, lst, lat, lon, 0.0);
        let (ra3, dec3) = apply_or_remove(ra2, dec2, lst, lat, lon, 0.0, false);
        assert!(arcsec(ra, ra3) < 1e-3);
        assert!(arcsec(dec, dec3) < 1e-3);
    }

    #[test]
    fn correction_direction_moves_tel_westward() {
        // Observer moves east (+RA). Stars appear shifted east (larger RA).
        // Removing diurnal aberration moves the direction back west
        // (smaller RA). Meridian star at LST=0: tel_ra reported as slightly
        // larger than truth, so corrected RA should be smaller.
        let lat = Angle::from_degrees(0.0);
        let lst = Angle::from_hours(0.0);
        let tel_ra = Angle::from_hours(0.0);
        let tel_dec = Angle::from_degrees(0.0);
        let (ra2, _) = apply_or_remove(tel_ra, tel_dec, lst, lat, Angle::from_radians(0.0), 0.0, false);
        let delta_rad = (ra2 - tel_ra).wrapped().radians();
        assert!(delta_rad < 0.0, "expected westward (negative) shift, got {}", delta_rad);
    }
}
