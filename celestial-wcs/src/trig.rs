use celestial_core::angle::Angle;

// Paper II states every angle in degrees. Reducing by whole quadrants in
// degrees is exact, so multiples of 90° give exact zeros and ones and the
// pole and seam tests downstream can compare against them.
pub(crate) fn sincosd(deg: f64) -> (f64, f64) {
    let (r, quadrant) = libm::remquo(deg, 90.0);
    let (s, c) = Angle::from_degrees(r).sin_cos();
    let (s, c) = match quadrant & 3 {
        0 => (s, c),
        1 => (c, -s),
        2 => (-s, -c),
        _ => (-c, s),
    };
    // An exact zero is +0 whichever quadrant produced it.
    (s + 0.0, c + 0.0)
}

pub(crate) fn sind(deg: f64) -> f64 {
    sincosd(deg).0
}

pub(crate) fn tand(deg: f64) -> f64 {
    let (s, c) = sincosd(deg);
    odd_multiple_of_45(deg, s, c).unwrap_or(s / c)
}

pub(crate) fn cotd(deg: f64) -> f64 {
    let (s, c) = sincosd(deg);
    odd_multiple_of_45(deg, s, c).unwrap_or(c / s)
}

// tan and cot are exactly ±1 at odd multiples of 45°, but sin 45° and cos 45°
// round an ulp apart and their ratio misses it. Mercator's reference point
// (eq. 87 at θ = 0) needs tan 45° = 1 to land on the origin.
fn odd_multiple_of_45(deg: f64, s: f64, c: f64) -> Option<f64> {
    let odd = libm::fabs(libm::remainder(deg, 90.0)) == 45.0;
    odd.then(|| libm::copysign(1.0, s * c))
}

pub(crate) fn atan2d(y: f64, x: f64) -> f64 {
    to_degrees(libm::atan2(y, x))
}

pub(crate) fn atand(v: f64) -> f64 {
    to_degrees(libm::atan(v))
}

pub(crate) fn asind(v: f64) -> f64 {
    to_degrees(libm::asin(v))
}

pub(crate) fn acosd(v: f64) -> f64 {
    to_degrees(libm::acos(v))
}

// Paper II marks each conversion between radians and degrees, as a factor of
// 180°/π or π/180°. Each of these rounds once.
pub(crate) fn to_degrees(rad: f64) -> f64 {
    Angle::from_radians(rad).degrees()
}

pub(crate) fn to_radians(deg: f64) -> f64 {
    Angle::from_degrees(deg).radians()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sincosd_exact_at_quadrants() {
        let cases = [
            (0.0, (0.0, 1.0)),
            (90.0, (1.0, 0.0)),
            (180.0, (0.0, -1.0)),
            (270.0, (-1.0, 0.0)),
            (360.0, (0.0, 1.0)),
            (-90.0, (-1.0, 0.0)),
            (-180.0, (0.0, -1.0)),
            (-270.0, (1.0, 0.0)),
            (450.0, (1.0, 0.0)),
            (-720.0, (0.0, 1.0)),
            (1.0e6 * 90.0 + 90.0, (1.0, 0.0)),
        ];
        for (deg, expected) in cases {
            let got = sincosd(deg);
            assert_eq!(got, expected, "{deg}");
            assert!(got.0.is_sign_positive() || got.0 != 0.0, "{deg}: -0 sine");
            assert!(got.1.is_sign_positive() || got.1 != 0.0, "{deg}: -0 cosine");
        }
    }

    #[test]
    fn test_sincosd_reduces_by_whole_quadrants() {
        for deg in [10.0, 30.0, 60.0, 123.5, -37.25] {
            let (s, c) = sincosd(deg);
            assert_eq!(sincosd(deg + 90.0), (c, -s), "{deg} + 90");
            assert_eq!(sincosd(deg + 180.0), (-s, -c), "{deg} + 180");
            assert_eq!(sincosd(deg - 90.0), (-c, s), "{deg} - 90");
            assert_eq!(sincosd(-deg), (-s, c), "-{deg}");
        }
    }

    #[test]
    fn test_sincosd_matches_angle_inside_first_octant() {
        for deg in [1.0, 30.0, 44.999] {
            assert_eq!(sincosd(deg), Angle::from_degrees(deg).sin_cos(), "{deg}");
        }
        assert_eq!(sind(30.0), 0.5);
        assert_eq!(sincosd(60.0).1, 0.5);
        assert_eq!(sind(150.0), 0.5);
    }

    #[test]
    fn test_inverse_functions_exact_at_quadrants() {
        assert_eq!(atan2d(0.0, 1.0), 0.0);
        assert_eq!(atan2d(1.0, 0.0), 90.0);
        assert_eq!(atan2d(0.0, -1.0), 180.0);
        assert_eq!(atan2d(-1.0, 0.0), -90.0);
        assert_eq!(atan2d(1.0, 1.0), 45.0);
        assert_eq!(acosd(1.0), 0.0);
        assert_eq!(acosd(0.0), 90.0);
        assert_eq!(acosd(-1.0), 180.0);
        assert_eq!(asind(1.0), 90.0);
        assert_eq!(asind(-1.0), -90.0);
        assert_eq!(atand(1.0), 45.0);
        assert_eq!(atand(f64::INFINITY), 90.0);
    }

    #[test]
    fn test_tand_exact_where_representable() {
        assert_eq!(tand(0.0), 0.0);
        assert_eq!(tand(180.0), 0.0);
        assert_eq!(tand(90.0), f64::INFINITY);
        assert_eq!(tand(-90.0), f64::NEG_INFINITY);
        assert_eq!(tand(45.0), 1.0);
        assert_eq!(tand(-45.0), -1.0);
        assert_eq!(tand(135.0), -1.0);
        assert_eq!(tand(225.0), 1.0);
        assert_eq!(tand(30.0), sind(30.0) / sincosd(30.0).1);
    }

    #[test]
    fn test_cotd_exact_where_representable() {
        assert_eq!(cotd(90.0), 0.0);
        assert_eq!(cotd(-90.0), 0.0);
        assert_eq!(cotd(0.0), f64::INFINITY);
        assert_eq!(cotd(45.0), 1.0);
        assert_eq!(cotd(-45.0), -1.0);
        assert_eq!(cotd(135.0), -1.0);
        assert_eq!(cotd(30.0), sincosd(30.0).1 / sind(30.0));
    }

    #[test]
    fn test_conversions_round_once() {
        assert_eq!(to_radians(180.0), celestial_core::constants::PI);
        assert_eq!(to_degrees(celestial_core::constants::HALF_PI), 90.0);
    }
}
