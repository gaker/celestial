use crate::common::{require_on_sphere, wrap_longitude};
use crate::coordinate::{CelestialCoord, NativeCoord};
use crate::error::{WcsError, WcsResult};
use crate::trig::{atan2d, sincosd};

use super::pole;

// Angles are in degrees, as Paper II states its equations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphericalRotation {
    alpha_p: f64,
    delta_p: f64,
    phi_p: f64,
    sin_delta_p: f64,
    cos_delta_p: f64,
    // Only used when δp = ±90°: a meridian (α_r, φ_r) common to both frames, from
    // which α − α_r = φ − φ_r (north) or φ_r − φ (south). Kept as a pair, not a
    // difference, so that φ_r itself maps to α_r exactly.
    pole_meridian: (f64, f64),
    // Only used when |δp| < 90°: where each frame measures longitude from.
    celestial_origin: Origin,
    native_origin: Origin,
}

// Eqs. 2 and 5 measure longitude from αp and φp. Measuring it from a point in
// the field instead, as α − αp = (α − α_r) − (αp − α_r), keeps the offsets
// small near that point, where αp and the angle eq. 2 gives can each exceed
// 100° and round at twice the ulp of the result.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Origin {
    lon: f64,
    // The sine and cosine of the other frame's pole longitude (αp or φp) less
    // `lon`, and the same up to a positive factor, which atan2 takes without
    // the rounding of normalizing it.
    to_pole: (f64, f64),
    to_pole_unscaled: (f64, f64),
}

impl Origin {
    fn new(lon: f64, to_pole: (f64, f64), to_pole_unscaled: (f64, f64)) -> Self {
        Self {
            lon,
            to_pole,
            to_pole_unscaled,
        }
    }

    fn at_pole(lon: f64) -> Self {
        Self::new(lon, (0.0, 1.0), (0.0, 1.0))
    }
}

impl SphericalRotation {
    pub fn new(alpha_p: f64, delta_p: f64, phi_p: f64) -> Self {
        let (sin_delta_p, cos_delta_p) = sincosd(delta_p);
        // Eqs. 3 and 4 put φp at αp − 180° and αp respectively.
        let alpha_r = if delta_p > 0.0 {
            alpha_p - 180.0
        } else {
            alpha_p
        };
        Self {
            alpha_p,
            delta_p,
            phi_p,
            sin_delta_p,
            cos_delta_p,
            pole_meridian: (alpha_r, phi_p),
            celestial_origin: Origin::at_pole(alpha_p),
            native_origin: Origin::at_pole(phi_p),
        }
    }

    // Paper II §2.2 gives 0 or 180° (φ0 is 0 there) as the condition that
    // celestial latitude increase with native latitude at the fiducial point.
    // For any φ0 that puts the celestial pole on the meridian φ0 when δ0 ≥ θ0,
    // and on the opposite one otherwise.
    pub(crate) fn default_lonpole(delta_0: f64, phi_0: f64, theta_0: f64) -> f64 {
        if delta_0 >= theta_0 {
            phi_0
        } else {
            phi_0 + 180.0
        }
    }

    pub fn from_crval(
        alpha_0: f64,
        delta_0: f64,
        phi_0: f64,
        theta_0: f64,
        lonpole: Option<f64>,
        latpole: Option<f64>,
    ) -> WcsResult<Self> {
        let given = [alpha_0, delta_0, phi_0, theta_0]
            .into_iter()
            .chain(lonpole);
        Self::require_valid(given.chain(latpole), theta_0)?;
        let phi_p = lonpole.unwrap_or_else(|| Self::default_lonpole(delta_0, phi_0, theta_0));
        if theta_0 == 90.0 {
            return Ok(Self::new(wrap_longitude(alpha_0)?, delta_0, phi_p));
        }
        let phi_offset = phi_p - phi_0;
        let delta_p = pole::delta_p(delta_0, theta_0, phi_offset, latpole.unwrap_or(90.0))?;
        let alpha_p = pole::alpha_p(alpha_0, delta_0, theta_0, phi_offset, delta_p);
        let rotation = Self::new(wrap_longitude(alpha_p)?, delta_p, phi_p);
        let crval = (wrap_longitude(alpha_0)?, delta_0);
        Ok(rotation.measured_from(crval, (phi_0, theta_0)))
    }

    // With coincident poles the fiducial point pins the shared meridian
    // exactly, as (α0, φ0); otherwise longitude is measured from it. A
    // fiducial point on the celestial pole leaves both to the poles.
    fn measured_from(self, (alpha_0, delta_0): (f64, f64), fiducial: (f64, f64)) -> Self {
        if pole::near_pole(delta_0) {
            return self;
        }
        if self.delta_p.abs() == 90.0 {
            return Self {
                pole_meridian: (alpha_0, fiducial.0),
                ..self
            };
        }
        self.with_origins((alpha_0, delta_0), fiducial)
    }

    fn with_origins(self, (alpha_0, delta_0): (f64, f64), (phi_0, theta_0): (f64, f64)) -> Self {
        let phi_offset = self.phi_p - phi_0;
        let (sin_term, cos_term) =
            pole::offset_from_pole(delta_0, theta_0, phi_offset, self.delta_p);
        let norm = libm::hypot(sin_term, cos_term);
        // Both vanish only where cos δ0 cos δp rounds away. αp is then α0,
        // and the pole serves as the origin.
        if norm == 0.0 {
            return self;
        }
        let to_alpha_p = (-sin_term / norm, cos_term / norm);
        let to_phi_p = sincosd(phi_offset);
        Self {
            celestial_origin: Origin::new(alpha_0, to_alpha_p, (-sin_term, cos_term)),
            native_origin: Origin::new(phi_0, to_phi_p, to_phi_p),
            ..self
        }
    }

    fn require_valid(mut angles: impl Iterator<Item = f64>, theta_0: f64) -> WcsResult<()> {
        if !angles.all(f64::is_finite) {
            return Err(WcsError::invalid_parameter(
                "rotation reference angles must be finite",
            ));
        }
        if theta_0.abs() > 90.0 {
            return Err(WcsError::invalid_parameter(
                "fiducial native latitude θ₀ lies beyond a pole",
            ));
        }
        Ok(())
    }

    pub fn native_to_celestial(&self, native: NativeCoord) -> WcsResult<CelestialCoord> {
        let (alpha, delta) = self.rotate(native.phi_deg(), native.theta_deg(), false)?;
        Ok(CelestialCoord::new(alpha, delta))
    }

    pub fn celestial_to_native(&self, celestial: CelestialCoord) -> WcsResult<NativeCoord> {
        let (phi, theta) = self.rotate(celestial.alpha_deg(), celestial.delta_deg(), true)?;
        Ok(NativeCoord::new(phi, theta))
    }

    // Eqs. 2 and 5 are the same rotation with the longitude origins swapped.
    fn rotate(&self, lon: f64, lat: f64, to_native: bool) -> WcsResult<(f64, f64)> {
        require_on_sphere(lon, lat)?;
        let (lon, lat) = if self.delta_p.abs() == 90.0 {
            self.shift_about_shared_pole(lon, lat, to_native)
        } else {
            self.rotate_oblique(lon, lat, to_native)
        };
        Ok((wrap_longitude(lon)?, lat))
    }

    fn rotate_oblique(&self, lon: f64, lat: f64, to_native: bool) -> (f64, f64) {
        let (from, to, pole_lon) = self.origins(to_native);
        if libm::fabs(lat) == 90.0 {
            return self.input_pole_image(lat, pole_lon);
        }
        let (sin_turn, cos_turn) = from.to_pole;
        let offset = sincosd(longitude_offset(lon, from.lon));
        let (sin_dlon, cos_dlon) = add_angles(offset, (-sin_turn, cos_turn));
        if sin_dlon == 0.0 {
            let (arg, lat) = self.rotate_on_meridian(lat, cos_dlon);
            return (pole_lon + arg, lat);
        }
        let (sin_lat, cos_lat) = sincosd(lat);
        let x = sin_lat * self.cos_delta_p - cos_lat * self.sin_delta_p * cos_dlon;
        let y = -cos_lat * sin_dlon;
        let z = sin_lat * self.sin_delta_p + cos_lat * self.cos_delta_p * cos_dlon;
        let (y_to, x_to) = add_angles((y, x), to.to_pole_unscaled);
        // atan2 keeps the latitude error near an ulp at every latitude, where
        // the paper's asin(z) degrades toward the poles.
        (to.lon + atan2d(y_to, x_to), atan2d(z, libm::hypot(x, y)))
    }

    // The origins of input and output longitude, and the output longitude of
    // the input frame's pole.
    fn origins(&self, to_native: bool) -> (&Origin, &Origin, f64) {
        if to_native {
            (&self.celestial_origin, &self.native_origin, self.phi_p)
        } else {
            (&self.native_origin, &self.celestial_origin, self.alpha_p)
        }
    }

    // A pole of the input frame lies on every meridian, at latitude ±δp in
    // the output frame.
    fn input_pole_image(&self, lat: f64, pole_lon: f64) -> (f64, f64) {
        if lat > 0.0 {
            (pole_lon, self.delta_p)
        } else {
            (pole_lon + 180.0, -self.delta_p)
        }
    }

    // On the meridian through both poles eq. 2 turns the point along that
    // meridian alone: x and z are the sine and cosine of lat − δp, or of
    // lat + δp with the sign of z flipped, so the output latitude is exact
    // in degrees.
    fn rotate_on_meridian(&self, lat: f64, cos_dlon: f64) -> (f64, f64) {
        let (along, lat) = if cos_dlon > 0.0 {
            let along = lat - self.delta_p;
            (along, 90.0 - along.abs())
        } else {
            let along = lat + self.delta_p;
            (along, along.abs() - 90.0)
        };
        (self.arg(along, 0.0, lat), lat)
    }

    // arg(x, y) has no value at a pole of the output frame. There the
    // longitude is its limit along the input meridian through the pole,
    // approached from the input equator, or from the pole's own hemisphere
    // when the pole lies on that equator.
    fn arg(&self, x: f64, y: f64, z: f64) -> f64 {
        if x != 0.0 || y != 0.0 {
            return atan2d(y, x);
        }
        if (z > 0.0) == (self.delta_p > 0.0) {
            180.0
        } else {
            0.0
        }
    }

    // With the poles coincident the rotation is a longitude shift. The general
    // form only reaches it through rounded cos(±90°) products, and going
    // through αp and φp would add and remove 180°.
    fn shift_about_shared_pole(&self, lon: f64, lat: f64, to_native: bool) -> (f64, f64) {
        let (alpha_r, phi_r) = self.pole_meridian;
        match (self.delta_p > 0.0, to_native) {
            (true, false) => ((lon - phi_r) + alpha_r, lat),
            (true, true) => ((lon - alpha_r) + phi_r, lat),
            (false, false) => (alpha_r - (lon - phi_r), -lat),
            (false, true) => (phi_r - (lon - alpha_r), -lat),
        }
    }

    #[inline]
    pub fn phi_p_deg(&self) -> f64 {
        self.phi_p
    }

    #[inline]
    pub fn delta_p_deg(&self) -> f64 {
        self.delta_p
    }
}

// A pair straddling the ±180° seam differs by nearly ±360°, and rounding that
// difference loses most of the small offset between them. The reduction by
// 360° is exact, and carrying the rounding error past it leaves a single
// rounding.
fn longitude_offset(lon: f64, origin: f64) -> f64 {
    let raw = lon - origin;
    if libm::fabs(raw) <= 180.0 {
        return raw;
    }
    libm::remainder(raw, 360.0) + difference_error(lon, origin, raw)
}

// The sine and cosine of a + b from those of a and b. A pair may carry a
// positive factor, and the result then carries it too.
fn add_angles((sin_a, cos_a): (f64, f64), (sin_b, cos_b): (f64, f64)) -> (f64, f64) {
    (sin_a * cos_b + cos_a * sin_b, cos_a * cos_b - sin_a * sin_b)
}

// Knuth's TwoSum: the exact amount by which `raw = a − b` was rounded.
fn difference_error(a: f64, b: f64, raw: f64) -> f64 {
    let a_part = raw + b;
    let b_part = a_part - raw;
    (a - a_part) + (b_part - b)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn crval_rotation(crval: (f64, f64), theta_0: f64, lonpole: Option<f64>) -> SphericalRotation {
        SphericalRotation::from_crval(crval.0, crval.1, 0.0, theta_0, lonpole, None).unwrap()
    }

    fn to_celestial(rot: &SphericalRotation, phi: f64, theta: f64) -> (f64, f64) {
        let native = NativeCoord::new(phi, theta);
        let c = rot.native_to_celestial(native).unwrap();
        (c.alpha_deg(), c.delta_deg())
    }

    fn to_native(rot: &SphericalRotation, alpha: f64, delta: f64) -> (f64, f64) {
        let sky = CelestialCoord::new(alpha, delta);
        let n = rot.celestial_to_native(sky).unwrap();
        (n.phi_deg(), n.theta_deg())
    }

    #[test]
    fn test_from_crval_southern_cylindrical() {
        // δ0 < θ0 takes the default LONPOLE of 180°, and αp then depends on δp.
        let rot = crval_rotation((37.5, -22.5), 0.0, None);
        assert_eq!(rot.delta_p_deg(), 67.5);
        assert_eq!(to_celestial(&rot, 0.0, 0.0).0, 37.5);
        // 60 digits give −7.954231362487325.
        assert_eq!(to_native(&rot, 50.0, -30.0).1, -7.954231362487324);
    }

    #[test]
    fn test_from_crval_meridian_pole_is_exact() {
        // With φp − φ0 a multiple of 180°, eq. 8 is θ0 ± (90° − δ0) or
        // ±180° − θ0 ± (90° − δ0), exact in degrees, so the fiducial point
        // returns CRVAL.
        let cases = [
            ((150.0, 40.0), (30.0, 22.5), None, 72.5),
            ((150.0, 10.0), (0.0, 22.5), None, 77.5),
            ((37.5, -22.5), (0.0, 0.0), None, 67.5),
            ((10.0, 5.0), (0.0, 0.0), Some(0.0), 85.0),
        ];
        for (crval, fiducial, lonpole, delta_p) in cases {
            let rot = fiducial_rotation(crval, fiducial, lonpole, None);
            assert_eq!(rot.delta_p_deg(), delta_p, "{crval:?} {fiducial:?}");
            let back = to_celestial(&rot, fiducial.0, fiducial.1);
            assert_eq!(back, crval, "{crval:?} {fiducial:?}");
        }
    }

    #[test]
    fn test_from_crval_takes_eq8_roots_as_angles() {
        // Eq. 8 gives −150° ± 135°; −285° is the angle 75°, a valid root.
        let northern = fiducial_rotation((10.0, -45.0), (0.0, -30.0), None, None);
        assert_eq!(northern.delta_p_deg(), 75.0);
        assert_eq!(to_celestial(&northern, 0.0, -30.0), (10.0, -45.0));
        let southern = fiducial_rotation((10.0, -45.0), (0.0, -30.0), None, Some(-90.0));
        assert_eq!(southern.delta_p_deg(), -15.0);
        assert_eq!(to_celestial(&southern, 0.0, -30.0), (10.0, -45.0));
    }

    #[test]
    fn test_from_crval_lonpole_off_meridian() {
        // 60 digits give δp = 24.627121690838486; the inverse cosine in eq. 8
        // has slope 2.4 here and magnifies the rounding of its argument.
        let rot = crval_rotation((100.0, 40.0), 0.0, Some(45.0));
        assert_eq!(rot.delta_p_deg(), 24.627121690838507);
        assert_eq!(
            to_celestial(&rot, 10.0, 5.0),
            (101.91313532571672, 51.08934021953901)
        );
        let (alpha, delta) = to_celestial(&rot, 45.0, 30.0);
        assert_eq!(
            (alpha.rem_euclid(360.0), delta),
            (347.3779636918118, 84.62712169083851)
        );
    }

    // Here αp is −92.56° and the angle eq. 2 gives near CRVAL about 130°;
    // both round at twice the ulp of α, which is near 37.5°.
    #[test]
    fn test_from_crval_lonpole_longitude_near_crval() {
        let rot = crval_rotation((37.5, 22.5), 0.0, Some(45.0));
        assert_eq!(to_celestial(&rot, 0.0, 0.0), (37.5, 22.5));
        assert_eq!(to_native(&rot, 37.5, 22.5), (0.0, 0.0));
        // 60 digits give 39.822657095051306 and 20.72537425691167.
        assert_eq!(
            to_celestial(&rot, 1.25, -2.5),
            (39.822657095051305, 20.725374256911667)
        );
        // 60 digits give 1.6116739364859836 and −2.085211299351829. Rounding
        // the inputs alone moves φ by 70 ulp.
        assert_eq!(
            to_native(&rot, 40.0, 21.25),
            (1.611673936485985, -2.0852112993518337)
        );
    }

    // A pole of the input frame lies on every meridian, so it lands on
    // latitude ±δp exactly, as the meridian through both poles does.
    #[test]
    fn test_input_pole_lands_on_delta_p() {
        let given = SphericalRotation::new(100.0, -25.0, 180.0);
        assert_eq!(to_native(&given, 13.0, 90.0), (180.0, -25.0));
        assert_eq!(to_native(&given, 13.0, -90.0), (0.0, 25.0));
        let rot = crval_rotation((37.5, 22.5), 0.0, Some(45.0));
        let delta_p = rot.delta_p_deg();
        assert_eq!(to_native(&rot, 10.0, 90.0), (45.0, delta_p));
        assert_eq!(to_native(&rot, 10.0, -90.0), (-135.0, -delta_p));
    }

    #[test]
    fn test_rotation_latitude_near_native_pole() {
        let rot = crval_rotation((200.0, -30.0), 90.0, None);
        assert_eq!(to_native(&rot, 200.0, -30.0).1, 90.0);
        // 60 digits give 89.99986771246301 and −29.999799999949634.
        assert_eq!(to_native(&rot, 200.0001, -30.0001).1, 89.999867712463);
        assert_eq!(
            to_celestial(&rot, 153.43494882292202, 89.99977639320227).1,
            -29.999799999949634
        );
    }

    #[test]
    fn test_rotation_longitude_seam_is_one_meridian() {
        // α = ±180° differ from αp by 0 or exactly one wrap period, whose sine
        // is not 0 in floating point.
        let rot = crval_rotation((180.0, 45.0), 90.0, None);
        assert_eq!(to_native(&rot, -180.0, 45.0).1, 90.0);
        assert_eq!(to_native(&rot, 180.0, 45.0).1, 90.0);
    }

    #[test]
    fn test_rotation_same_meridian_either_side_of_seam() {
        // α and α − 360° name one meridian. Near the native pole, the raw
        // difference of the pair straddling the seam from αp rounds away the
        // small offset.
        let rot = crval_rotation((200.0, -30.0), 90.0, None);
        let alpha = 200.0001;
        let native = |a: f64| to_native(&rot, a, -30.0001);
        assert_eq!(native(alpha - 360.0), native(alpha));
    }

    #[test]
    fn test_from_crval_zenithal_pole_is_crval() {
        let rot = crval_rotation((123.456, 37.5), 90.0, None);
        assert_eq!(rot.delta_p_deg(), 37.5);
        assert_eq!(to_celestial(&rot, 0.0, 90.0).1, 37.5);
        assert_eq!(to_native(&rot, 123.456, 37.5).1, 90.0);
        // 60 digits give 37.74832107683319.
        assert_eq!(
            to_celestial(&rot, 116.56505117707799, 89.44097413616623).1,
            37.748321076833186
        );
    }

    #[test]
    fn test_rotation_about_shared_pole_is_longitude_shift() {
        // CRVAL on the equator with θ0 = 0 puts the native pole on the celestial pole.
        let rot = crval_rotation((37.5, 0.0), 0.0, None);
        assert_eq!(rot.delta_p_deg(), 90.0);
        assert_eq!(to_native(&rot, 0.0, 90.0), (-37.5, 90.0));
        assert_eq!(to_native(&rot, 0.0, -90.0), (-37.5, -90.0));
        assert_eq!(to_native(&rot, 10.0, 45.0), (-27.5, 45.0));
        assert_eq!(to_celestial(&rot, 10.0, 45.0), (47.5, 45.0));
        assert_eq!(to_celestial(&rot, 0.0, 0.0), (37.5, 0.0));
    }

    // At a pole of the output frame longitude has no value, and the one that
    // comes back stands for it.
    #[test]
    fn test_spherical_rotation_known_points() {
        let coincident = SphericalRotation::new(0.0, 90.0, 180.0);
        assert_eq!(to_celestial(&coincident, 0.0, 90.0), (0.0, 90.0));
        assert_eq!(to_celestial(&coincident, 90.0, 0.0), (90.0, 0.0));
        let oblique = SphericalRotation::new(180.0, 45.0, 180.0);
        assert_eq!(to_celestial(&oblique, 0.0, 90.0), (180.0, 45.0));
        assert_eq!(to_native(&oblique, 180.0, 45.0), (0.0, 90.0));
    }

    // Each direction pinned one way against 60 digits: within 2.60 ulp, at
    // most half the summed half-ulp rounding of the steps as evaluated.
    #[test]
    fn test_spherical_rotation_pinned_both_ways() {
        let north = SphericalRotation::new(120.0, 35.0, 180.0);
        let sky = (141.18557055266368, 11.953500173493248);
        assert_eq!(to_celestial(&north, 45.0, 60.0), sky);
        assert_eq!(to_native(&north, sky.0, sky.1), (45.0, 59.99999999999998));
        let south = SphericalRotation::new(100.0, -25.0, 180.0);
        let native = (58.34492397467059, 79.82426722795314);
        assert_eq!(to_native(&south, 110.0, -30.0), native);
        assert_eq!(
            to_celestial(&south, native.0, native.1),
            (110.0, -30.000000000000004)
        );
    }

    // With θ0 = 90° the native pole is CRVAL, the celestial pole included.
    #[test]
    fn test_from_crval_zenithal_at_pole() {
        for crval in [(180.0, 45.0), (0.0, 90.0)] {
            let rot = crval_rotation(crval, 90.0, Some(180.0));
            assert_eq!(to_celestial(&rot, 0.0, 90.0), crval);
        }
    }

    #[test]
    fn test_from_crval_roundtrip() {
        let rot = SphericalRotation::from_crval(120.0, 35.0, 0.0, 90.0, Some(180.0), None).unwrap();

        // 60 digits give the same celestial pair. Rounding it to doubles moves
        // φ by 7 ulp, so the trip back cannot return 30° exactly; 60 digits
        // give 30.000000000000025 and 59.99999999999999 from it.
        let celestial = to_celestial(&rot, 30.0, 60.0);
        assert_eq!(celestial, (134.62907023771479, 8.165245288585398));
        let recovered = to_native(&rot, celestial.0, celestial.1);
        assert_eq!(recovered, (30.00000000000003, 59.99999999999999));
    }

    // Within 1.75 ulp of 60 digits, at most 0.66 of the summed half-ulp
    // rounding of the steps.
    #[test]
    fn test_explicit_lonpole_overrides_default() {
        let default = crval_rotation((180.0, 45.0), 90.0, None);
        let rotated = crval_rotation((180.0, 45.0), 90.0, Some(90.0));
        assert_eq!(
            to_celestial(&default, 30.0, 60.0),
            (-164.77484298932217, 17.829543848069374)
        );
        assert_eq!(
            to_celestial(&rotated, 30.0, 60.0),
            (-135.17038375556047, 52.10606741594716)
        );
    }

    // θ0 = 0, δ0 = 45° and LONPOLE 0 give eq. 8 the roots ±45°, both valid.
    // LATPOLE picks the nearer, and 90° stands in for it when absent. Within
    // 3.11 ulp of 60 digits, at most 0.55 of the summed half-ulp rounding.
    #[test]
    fn test_latpole_picks_between_two_poles() {
        let crval = (180.0, 45.0);
        let rotation = |latpole| fiducial_rotation(crval, (0.0, 0.0), Some(0.0), latpole);
        let (north, south) = (rotation(Some(90.0)), rotation(Some(-90.0)));
        assert_eq!(rotation(None), north);
        assert_eq!((north.delta_p, south.delta_p), (45.0, -45.0));
        assert_eq!(to_celestial(&north, 0.0, 0.0), crval);
        assert_eq!(to_celestial(&south, 0.0, 0.0), crval);
        let north_sky = (-112.20765429859648, 62.1144331639063);
        assert_eq!(to_celestial(&north, 30.0, 30.0), north_sky);
        let south_sky = (153.8998612177136, 10.182067403158909);
        assert_eq!(to_celestial(&south, 30.0, 30.0), south_sky);
    }

    // With θ0 = 90° CRVAL fixes δp, leaving LATPOLE nothing to choose. The
    // latitude is 7.88 ulp low against 60 digits, 0.35 of the summed half-ulp
    // rounding: z cancels to 0.047 from terms near 0.35.
    #[test]
    fn test_latpole_is_unused_at_a_zenithal_fiducial() {
        let given = fiducial_rotation((180.0, 45.0), (0.0, 90.0), Some(180.0), Some(45.0));
        let absent = fiducial_rotation((180.0, 45.0), (0.0, 90.0), Some(180.0), None);
        assert_eq!(given, absent);
        assert_eq!(
            to_celestial(&given, 60.0, 30.0),
            (-131.33656448563264, 2.7149549662644747)
        );
    }

    // Here δp is 89.9999998°, and sin θ0, sin δ0 and sin δp all round to 1,
    // so both terms of eqs. 9 and 10 vanish. αp is then α0, and longitude is
    // measured from the poles. 60 digits give −139.9999999 and
    // 45.00000017320508; δp rounding 0.84 ulp high moves the latitude 1.2.
    #[test]
    fn test_fiducial_next_to_pole_measures_from_poles() {
        let rot = fiducial_rotation((10.0, 89.9999999), (0.0, 89.9999999), None, Some(-90.0));
        assert_eq!(rot.alpha_p, 10.0);
        assert_eq!(
            to_celestial(&rot, 30.0, 45.0),
            (-139.9999999, 45.00000017320506)
        );
    }

    // CRVAL at the celestial pole exercises Paper II rule 1 (δ0 = ±90° gives
    // αp = α0); the poles coincide and the rotation is a shift of longitude.
    #[test]
    fn test_from_crval_pole_crval() {
        let rot = SphericalRotation::from_crval(0.0, 90.0, 0.0, 90.0, None, None).unwrap();
        assert_eq!((rot.alpha_p, rot.delta_p), (0.0, 90.0));
        assert_eq!(to_celestial(&rot, 45.0, 80.0), (-135.0, 80.0));
    }

    fn fiducial_rotation(
        crval: (f64, f64),
        fiducial: (f64, f64),
        lonpole: Option<f64>,
        latpole: Option<f64>,
    ) -> SphericalRotation {
        SphericalRotation::from_crval(crval.0, crval.1, fiducial.0, fiducial.1, lonpole, latpole)
            .unwrap()
    }

    // The default is φ0 when δ0 ≥ θ0 and φ0 + 180° otherwise.
    #[test]
    fn test_default_lonpole_is_measured_from_phi_0() {
        let cases = [
            ((150.0, 40.0), (30.0, 22.5), 30.0),
            ((150.0, 10.0), (30.0, 22.5), 210.0),
            ((180.0, 45.0), (0.0, 90.0), 180.0),
            ((180.0, 45.0), (0.0, 0.0), 0.0),
        ];
        for (crval, fiducial, lonpole) in cases {
            let default = fiducial_rotation(crval, fiducial, None, None);
            let explicit = fiducial_rotation(crval, fiducial, Some(lonpole), None);
            assert_eq!(default, explicit, "{crval:?} {fiducial:?}");
        }
    }

    #[test]
    fn test_coincident_poles_put_crval_at_phi_0() {
        // With θ0 = δ0 = 0 and the default LONPOLE, eq. 8 gives δp = ±90° and
        // LATPOLE picks the sign; the rotation is then a shift of longitude.
        let cases = [
            (
                90.0,
                [(30.0, 0.0), (150.0, 0.0)],
                [(40.0, 10.0), (160.0, 10.0)],
            ),
            (
                -90.0,
                [(30.0, 0.0), (150.0, 0.0)],
                [(40.0, 10.0), (140.0, -10.0)],
            ),
        ];
        for (latpole, at_fiducial, off_fiducial) in cases {
            let rot = fiducial_rotation((150.0, 0.0), (30.0, 0.0), None, Some(latpole));
            for [native, celestial] in [at_fiducial, off_fiducial] {
                let got = to_celestial(&rot, native.0, native.1);
                assert_eq!(got, celestial, "LATPOLE {latpole}, {native:?}");
            }
        }
    }

    #[test]
    fn test_from_crval_rejects_a_fiducial_latitude_beyond_a_pole() {
        let rot = SphericalRotation::from_crval(150.0, 40.0, 0.0, 90.5, None, None);
        assert!(
            matches!(rot, Err(WcsError::InvalidParameter { .. })),
            "{rot:?}"
        );
    }

    // CRVAL comes back exactly, its longitude in the normal range.
    #[test]
    fn test_from_crval_theta0_zero_maps_fiducial_to_crval() {
        let alphas = [
            (0.0, 0.0),
            (45.0, 45.0),
            (90.0, 90.0),
            (180.0, 180.0),
            (270.0, -90.0),
            (315.0, -45.0),
        ];
        for (alpha_0, alpha) in alphas {
            let rot = crval_rotation((alpha_0, 0.0), 0.0, None);
            assert_eq!(to_celestial(&rot, 0.0, 0.0), (alpha, 0.0), "α0 = {alpha_0}");
        }
    }

    // A CAR reference point on the celestial pole: θ0 = 0, so δp = 0 and
    // αp = α0 (eqs. 8-10). Eq. 2 then takes (90°, 0°) to α0 + arg(0, −1).
    #[test]
    fn test_from_crval_reference_point_on_the_pole() {
        let rot = SphericalRotation::from_crval(30.0, 90.0, 0.0, 0.0, None, None).unwrap();
        assert_eq!((rot.phi_p_deg(), rot.delta_p_deg()), (0.0, 0.0));
        let at_crval = rot.native_to_celestial(NativeCoord::new(0.0, 0.0)).unwrap();
        assert_eq!(at_crval, CelestialCoord::new(30.0, 90.0));
        let quarter = rot
            .native_to_celestial(NativeCoord::new(90.0, 0.0))
            .unwrap();
        assert_eq!(quarter, CelestialCoord::new(-60.0, 0.0));
    }

    // With δp = −90° eq. 2 reduces to α = αp + φp − φ and δ = −θ.
    #[test]
    fn test_rotation_about_a_shared_south_pole() {
        let rot = SphericalRotation::new(0.0, -90.0, 180.0);
        let celestial = rot
            .native_to_celestial(NativeCoord::new(30.0, 20.0))
            .unwrap();
        assert_eq!(celestial, CelestialCoord::new(150.0, -20.0));
        let native = rot.celestial_to_native(celestial).unwrap();
        assert_eq!(native, NativeCoord::new(30.0, 20.0));
    }

    #[test]
    fn test_from_crval_refuses_non_finite_angles() {
        for (alpha_0, lonpole) in [(f64::NAN, None), (0.0, Some(f64::INFINITY))] {
            let result = SphericalRotation::from_crval(alpha_0, 45.0, 0.0, 90.0, lonpole, None);
            assert!(
                matches!(result, Err(WcsError::InvalidParameter { .. })),
                "{alpha_0} {lonpole:?}: {result:?}"
            );
        }
    }
}
