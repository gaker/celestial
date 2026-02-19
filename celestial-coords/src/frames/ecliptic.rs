use crate::{transforms::CoordinateFrame, CoordResult, Distance, ICRSPosition};
use celestial_core::{matrix::RotationMatrix3, Angle};
use celestial_time::{transforms::PrecessionCalculator, TT};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EclipticPosition {
    lambda: Angle,
    beta: Angle,
    epoch: TT,
    distance: Option<Distance>,
}

impl EclipticPosition {
    pub fn new(lambda: Angle, beta: Angle, epoch: TT) -> CoordResult<Self> {
        let lambda = lambda.validate_longitude(true)?;
        let beta = beta.validate_latitude()?;

        Ok(Self {
            lambda,
            beta,
            epoch,
            distance: None,
        })
    }

    pub fn with_distance(
        lambda: Angle,
        beta: Angle,
        epoch: TT,
        distance: Distance,
    ) -> CoordResult<Self> {
        let mut pos = Self::new(lambda, beta, epoch)?;
        pos.distance = Some(distance);
        Ok(pos)
    }

    pub fn from_degrees(lambda_deg: f64, beta_deg: f64, epoch: TT) -> CoordResult<Self> {
        Self::new(
            Angle::from_degrees(lambda_deg),
            Angle::from_degrees(beta_deg),
            epoch,
        )
    }

    pub fn lambda(&self) -> Angle {
        self.lambda
    }

    pub fn beta(&self) -> Angle {
        self.beta
    }

    pub fn epoch(&self) -> TT {
        self.epoch
    }

    pub fn distance(&self) -> Option<Distance> {
        self.distance
    }

    pub fn set_distance(&mut self, distance: Distance) {
        self.distance = Some(distance);
    }

    pub fn mean_obliquity(&self) -> Angle {
        let jd = self.epoch.to_julian_date();
        Angle::from_radians(celestial_core::obliquity::iau_2006_mean_obliquity(
            jd.jd1(),
            jd.jd2(),
        ))
    }

    pub fn true_obliquity(&self) -> CoordResult<Angle> {
        use celestial_time::transforms::nutation::NutationCalculator;

        let nutation =
            self.epoch
                .nutation_iau2006a()
                .map_err(|e| crate::CoordError::CoreError {
                    message: format!("Nutation calculation failed: {}", e),
                })?;

        let jd = self.epoch.to_julian_date();
        let mean_obliquity = celestial_core::obliquity::iau_2006_mean_obliquity(jd.jd1(), jd.jd2());

        let true_obliquity = mean_obliquity + nutation.nutation_obliquity();

        Ok(Angle::from_radians(true_obliquity))
    }

    pub fn vernal_equinox(epoch: TT) -> Self {
        Self {
            lambda: Angle::ZERO,
            beta: Angle::ZERO,
            epoch,
            distance: None,
        }
    }

    pub fn summer_solstice(epoch: TT) -> Self {
        Self {
            lambda: Angle::HALF_PI,
            beta: Angle::ZERO,
            epoch,
            distance: None,
        }
    }

    pub fn autumnal_equinox(epoch: TT) -> Self {
        Self {
            lambda: Angle::PI,
            beta: Angle::ZERO,
            epoch,
            distance: None,
        }
    }

    pub fn winter_solstice(epoch: TT) -> Self {
        Self {
            lambda: Angle::from_degrees(270.0),
            beta: Angle::ZERO,
            epoch,
            distance: None,
        }
    }

    pub fn north_ecliptic_pole(epoch: TT) -> Self {
        Self {
            lambda: Angle::ZERO,
            beta: Angle::HALF_PI,
            epoch,
            distance: None,
        }
    }

    pub fn south_ecliptic_pole(epoch: TT) -> Self {
        Self {
            lambda: Angle::ZERO,
            beta: -Angle::HALF_PI,
            epoch,
            distance: None,
        }
    }

    pub fn is_near_ecliptic_plane(&self) -> bool {
        self.beta.abs().degrees() < 5.0
    }

    pub fn is_near_ecliptic_pole(&self) -> bool {
        self.beta.abs().degrees() > 85.0
    }

    pub fn season_index(&self) -> u8 {
        let lambda_deg = self.lambda.degrees();
        if lambda_deg < 90.0 {
            0
        } else if lambda_deg < 180.0 {
            1
        } else if lambda_deg < 270.0 {
            2
        } else {
            3
        }
    }

    pub fn angular_separation(&self, other: &Self) -> Angle {
        let (sin_b1, cos_b1) = self.beta.sin_cos();
        let (sin_b2, cos_b2) = other.beta.sin_cos();
        let delta_lambda = (self.lambda - other.lambda).radians();

        let angle_rad = celestial_core::math::vincenty_angular_separation(
            sin_b1,
            cos_b1,
            sin_b2,
            cos_b2,
            delta_lambda,
        );

        Angle::from_radians(angle_rad)
    }
}

fn ecm06_matrix(epoch: &TT) -> CoordResult<RotationMatrix3> {
    let precession = epoch.precession()?;
    let bias_precession_matrix = precession.bias_precession_matrix;

    let jd = epoch.to_julian_date();
    let mean_obliquity = celestial_core::obliquity::iau_2006_mean_obliquity(jd.jd1(), jd.jd2());

    let mut ecliptic_rotation = RotationMatrix3::identity();
    ecliptic_rotation.rotate_x(mean_obliquity);

    Ok(ecliptic_rotation.multiply(&bias_precession_matrix))
}

impl CoordinateFrame for EclipticPosition {
    fn to_icrs(&self, epoch: &TT) -> CoordResult<ICRSPosition> {
        let lambda = self.lambda.radians();
        let beta = self.beta.radians();

        let (sin_beta, cos_beta) = libm::sincos(beta);
        let (sin_lambda, cos_lambda) = libm::sincos(lambda);

        let ecliptic_cartesian = [cos_beta * cos_lambda, cos_beta * sin_lambda, sin_beta];

        let icrs_to_ecliptic_matrix = ecm06_matrix(epoch)?;
        let icrs_cartesian = icrs_to_ecliptic_matrix
            .transpose()
            .apply_to_vector(ecliptic_cartesian);

        let x = icrs_cartesian[0];
        let y = icrs_cartesian[1];
        let z = icrs_cartesian[2];

        let rxy2 = x * x + y * y;
        let ra = if rxy2 == 0.0 { 0.0 } else { libm::atan2(y, x) };
        let dec = if z == 0.0 {
            0.0
        } else {
            libm::atan2(z, libm::sqrt(rxy2))
        };

        let d2pi = celestial_core::constants::TWOPI;
        let mut ra_normalized = ra % d2pi;
        if ra_normalized < 0.0 {
            ra_normalized += d2pi;
        }

        let dpi = celestial_core::constants::PI;
        let mut dec_normalized = dec % d2pi;
        if dec_normalized.abs() >= dpi {
            dec_normalized -= libm::copysign(d2pi, dec);
        }

        let mut icrs = ICRSPosition::new(
            Angle::from_radians(ra_normalized),
            Angle::from_radians(dec_normalized),
        )?;

        if let Some(distance) = self.distance {
            icrs.set_distance(distance);
        }

        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let ra = icrs.ra().radians();
        let dec = icrs.dec().radians();

        let (sin_dec, cos_dec) = libm::sincos(dec);
        let (sin_ra, cos_ra) = libm::sincos(ra);

        let icrs_cartesian = [cos_dec * cos_ra, cos_dec * sin_ra, sin_dec];

        let icrs_to_ecliptic_matrix = ecm06_matrix(epoch)?;
        let ecliptic_cartesian = icrs_to_ecliptic_matrix.apply_to_vector(icrs_cartesian);

        let x = ecliptic_cartesian[0];
        let y = ecliptic_cartesian[1];
        let z = ecliptic_cartesian[2];

        let rxy2 = x * x + y * y;
        let lambda = if rxy2 != 0.0 { libm::atan2(y, x) } else { 0.0 };
        let beta = if rxy2 != 0.0 || z != 0.0 {
            libm::atan2(z, libm::sqrt(rxy2))
        } else {
            0.0
        };

        let d2pi = celestial_core::constants::TWOPI;
        let mut lambda_normalized = lambda % d2pi;
        if lambda_normalized < 0.0 {
            lambda_normalized += d2pi;
        }

        let mut ecliptic = Self::new(
            Angle::from_radians(lambda_normalized),
            Angle::from_radians(beta),
            *epoch,
        )?;

        if let Some(distance) = icrs.distance() {
            ecliptic.set_distance(distance);
        }

        Ok(ecliptic)
    }
}

impl std::fmt::Display for EclipticPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Ecliptic(λ={:.6}°, β={:.6}°, epoch=J{:.1}",
            self.lambda.degrees(),
            self.beta.degrees(),
            self.epoch.julian_year()
        )?;

        if let Some(distance) = self.distance {
            write!(f, ", d={}", distance)?;
        }

        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Distance;

    mod erfa_reference {
        use super::*;

        #[test]
        fn test_obliquity_at_j2000() {
            let epoch = TT::j2000();
            let pos = EclipticPosition::from_degrees(0.0, 0.0, epoch).unwrap();
            let mean_obliquity = pos.mean_obliquity();
            assert_eq!(mean_obliquity.radians(), 4.09092600600582889658e-01);
        }

        #[test]
        fn test_ecm06_matrix_at_j2000() {
            let epoch = TT::j2000();
            let matrix = ecm06_matrix(&epoch).unwrap();
            let m = matrix.elements();

            assert_eq!(m[0][0], 9.99999999999994115818e-01);
            assert_eq!(m[0][1], -7.07836896097155612759e-08);
            assert_eq!(m[0][2], 8.05621397761318608390e-08);
            assert_eq!(m[1][0], 3.28970040774196464850e-08);
            assert_eq!(m[1][1], 9.17482129914958366435e-01);
            assert_eq!(m[1][2], 3.97776999444047929533e-01);
            assert_eq!(m[2][0], -1.02070447254843554005e-07);
            assert_eq!(m[2][1], -3.97776999444043044551e-01);
            assert_eq!(m[2][2], 9.17482129914955590877e-01);
        }

        #[test]
        fn test_north_ecliptic_pole_to_icrs() {
            let epoch = TT::j2000();
            let north_pole = EclipticPosition::north_ecliptic_pole(epoch);
            let icrs = north_pole.to_icrs(&epoch).unwrap();

            assert_eq!(icrs.ra().radians(), 4.71238872378250484019e+00);
            assert_eq!(icrs.dec().radians(), 1.16170369313486876450e+00);
        }

        #[test]
        fn test_south_ecliptic_pole_to_icrs() {
            let epoch = TT::j2000();
            let south_pole = EclipticPosition::south_ecliptic_pole(epoch);
            let icrs = south_pole.to_icrs(&epoch).unwrap();

            assert_eq!(icrs.ra().radians(), 1.57079607019271128010e+00);
            assert_eq!(icrs.dec().radians(), -1.16170369313486876450e+00);
        }

        #[test]
        fn test_roundtrip_decimal_coords() {
            let epoch = TT::j2000();

            // (123.456789°, 45.678901°)
            let original = EclipticPosition::from_degrees(123.456789, 45.678901, epoch).unwrap();
            let icrs = original.to_icrs(&epoch).unwrap();
            let roundtrip = EclipticPosition::from_icrs(&icrs, &epoch).unwrap();

            assert_eq!(
                original.lambda().radians(),
                roundtrip.lambda().radians(),
                "Lambda mismatch: orig={}, round={}",
                original.lambda().degrees(),
                roundtrip.lambda().degrees()
            );
        }

        #[test]
        fn test_roundtrip_negative_beta() {
            let epoch = TT::j2000();

            // (267.314159°, -23.271828°) roundtrip: ERFA shows lambda diff = 0 exactly
            let original = EclipticPosition::from_degrees(267.314159, -23.271828, epoch).unwrap();
            let icrs = original.to_icrs(&epoch).unwrap();
            let roundtrip = EclipticPosition::from_icrs(&icrs, &epoch).unwrap();

            assert_eq!(
                original.lambda().radians(),
                roundtrip.lambda().radians(),
                "Lambda mismatch: orig={}, round={}",
                original.lambda().degrees(),
                roundtrip.lambda().degrees()
            );
        }

        #[test]
        fn test_roundtrip_high_latitude() {
            let epoch = TT::j2000();

            // (45.123456°, 67.890123°) roundtrip
            let original = EclipticPosition::from_degrees(45.123456, 67.890123, epoch).unwrap();
            let icrs = original.to_icrs(&epoch).unwrap();
            let roundtrip = EclipticPosition::from_icrs(&icrs, &epoch).unwrap();

            // ERFA shows lambda diff = -1.11e-16 rad which wraps to 0
            // Our implementation should also produce 0 or very close
            let lambda_diff = (original.lambda().radians() - roundtrip.lambda().radians()).abs();
            assert!(
                lambda_diff < 1e-14,
                "Lambda diff too large: {} rad",
                lambda_diff
            );
        }
    }

    #[test]
    fn test_constructor_with_distance() {
        let epoch = TT::j2000();
        let distance = Distance::from_parsecs(10.0).unwrap();

        let pos = EclipticPosition::with_distance(
            Angle::from_degrees(180.0),
            Angle::from_degrees(45.0),
            epoch,
            distance,
        )
        .unwrap();

        assert_eq!(pos.lambda().degrees(), 180.0);
        assert_eq!(pos.beta().degrees(), 45.0);
        assert_eq!(pos.epoch(), epoch);
        assert_eq!(pos.distance().unwrap(), distance);
    }

    #[test]
    fn test_accessor_methods() {
        let epoch = TT::j2000();
        let distance = Distance::from_parsecs(5.0).unwrap();

        let mut pos = EclipticPosition::from_degrees(90.0, -30.0, epoch).unwrap();

        let expected_lambda = Angle::from_degrees(90.0).degrees();
        let expected_beta = Angle::from_degrees(-30.0).degrees();

        assert_eq!(pos.lambda().degrees(), expected_lambda);
        assert_eq!(pos.beta().degrees(), expected_beta);
        assert_eq!(pos.epoch(), epoch);
        assert_eq!(pos.distance(), None);

        pos.set_distance(distance);
        assert_eq!(pos.distance().unwrap(), distance);
    }

    #[test]
    fn test_obliquity_calculations() {
        let epoch = TT::j2000();
        let pos = EclipticPosition::from_degrees(0.0, 0.0, epoch).unwrap();

        let mean_obliquity = pos.mean_obliquity();
        let true_obliquity = pos.true_obliquity().unwrap();

        // IAU 2006 mean obliquity at J2000.0: 84381.406 arcseconds
        let expected_mean_obliquity_arcsec = 84381.406;
        let expected_mean_obliquity_deg = expected_mean_obliquity_arcsec / 3600.0;
        assert_eq!(mean_obliquity.degrees(), expected_mean_obliquity_deg);

        // True obliquity differs from mean by nutation in obliquity
        // At J2000.0, nutation is small but non-zero
        assert_ne!(true_obliquity.radians(), mean_obliquity.radians());
    }

    #[test]
    fn test_special_position_constructors() {
        let epoch = TT::j2000();

        let vernal = EclipticPosition::vernal_equinox(epoch);
        assert_eq!(vernal.lambda().degrees(), 0.0);
        assert_eq!(vernal.beta().degrees(), 0.0);
        assert_eq!(vernal.epoch(), epoch);
        assert_eq!(vernal.distance(), None);

        let summer = EclipticPosition::summer_solstice(epoch);
        assert_eq!(summer.lambda().degrees(), 90.0);
        assert_eq!(summer.beta().degrees(), 0.0);

        let autumn = EclipticPosition::autumnal_equinox(epoch);
        assert_eq!(autumn.lambda().degrees(), 180.0);
        assert_eq!(autumn.beta().degrees(), 0.0);

        let winter = EclipticPosition::winter_solstice(epoch);
        assert_eq!(winter.lambda().degrees(), 270.0);
        assert_eq!(winter.beta().degrees(), 0.0);

        let north_pole = EclipticPosition::north_ecliptic_pole(epoch);
        assert_eq!(north_pole.lambda().degrees(), 0.0);
        assert_eq!(north_pole.beta().degrees(), 90.0);

        let south_pole = EclipticPosition::south_ecliptic_pole(epoch);
        assert_eq!(south_pole.lambda().degrees(), 0.0);
        assert_eq!(south_pole.beta().degrees(), -90.0);
    }

    #[test]
    fn test_angular_separation() {
        let epoch = TT::j2000();

        let vernal = EclipticPosition::vernal_equinox(epoch);
        let summer = EclipticPosition::summer_solstice(epoch);
        let north_pole = EclipticPosition::north_ecliptic_pole(epoch);

        let sep_vernal_summer = vernal.angular_separation(&summer);
        assert_eq!(sep_vernal_summer.degrees(), 90.0);

        let sep_pole_vernal = north_pole.angular_separation(&vernal);
        assert_eq!(sep_pole_vernal.degrees(), 90.0);

        let sep_self = vernal.angular_separation(&vernal);
        assert_eq!(sep_self.degrees(), 0.0);
    }

    #[test]
    fn test_coordinate_transformations_roundtrip() {
        let epoch = TT::j2000();

        // Test roundtrip coordinate transformations
        // Lambda diff can be ~1e-15 rad due to numerical paths, but angular separation is ~0
        let test_cases = [(90.0, 0.0), (180.0, 0.0), (270.0, 0.0)];

        for (lambda_deg, beta_deg) in test_cases {
            let original = EclipticPosition::from_degrees(lambda_deg, beta_deg, epoch).unwrap();

            let icrs = original.to_icrs(&epoch).unwrap();
            let roundtrip = EclipticPosition::from_icrs(&icrs, &epoch).unwrap();

            // Verify angular separation is essentially zero
            // Accounts for ~1e-15 rad lambda drift and ~1e-16 rad beta drift
            let separation = original.angular_separation(&roundtrip);
            assert!(
                separation.radians() < 1e-14,
                "Separation too large for ({}, {}): {} rad",
                lambda_deg,
                beta_deg,
                separation.radians()
            );
        }
    }

    #[test]
    fn test_coordinate_transformations_roundtrip_zero_boundary() {
        // The 0/360 boundary: lambda 0° can become 6.28... due to atan2 returning near-2π
        // ERFA shows same behavior: lambda roundtrip from 0° goes through RA ~360° back to λ ~360°
        let epoch = TT::j2000();

        let original = EclipticPosition::from_degrees(0.0, 0.0, epoch).unwrap();
        let icrs = original.to_icrs(&epoch).unwrap();
        let roundtrip = EclipticPosition::from_icrs(&icrs, &epoch).unwrap();

        // Angular separation should be essentially zero even if lambda differs by ~2π
        // ERFA shows wrapped lambda diff = 0
        let separation = original.angular_separation(&roundtrip);
        assert!(
            separation.radians() < 1e-14,
            "Angular separation too large: {} rad",
            separation.radians()
        );
    }

    #[test]
    fn test_coordinate_transformations_with_distance() {
        let epoch = TT::j2000();
        let distance = Distance::from_parsecs(10.0).unwrap();

        let original = EclipticPosition::with_distance(
            Angle::from_degrees(45.0),
            Angle::from_degrees(30.0),
            epoch,
            distance,
        )
        .unwrap();

        let icrs = original.to_icrs(&epoch).unwrap();
        assert_eq!(icrs.distance().unwrap(), distance);

        let roundtrip = EclipticPosition::from_icrs(&icrs, &epoch).unwrap();
        assert_eq!(roundtrip.distance().unwrap(), distance);
    }

    #[test]
    fn test_display_formatting() {
        let epoch = TT::j2000();
        let distance = Distance::from_parsecs(5.0).unwrap();

        let pos_no_dist = EclipticPosition::from_degrees(45.123456, -30.987654, epoch).unwrap();
        let display_no_dist = format!("{}", pos_no_dist);
        assert!(display_no_dist.contains("λ=45.123456°"));
        assert!(display_no_dist.contains("β=-30.987654°"));
        assert!(display_no_dist.contains("epoch=J2000.0"));
        assert!(!display_no_dist.contains("d="));

        let mut pos_with_dist = pos_no_dist.clone();
        pos_with_dist.set_distance(distance);
        let display_with_dist = format!("{}", pos_with_dist);
        assert!(display_with_dist.contains("λ=45.123456°"));
        assert!(display_with_dist.contains("β=-30.987654°"));
        assert!(display_with_dist.contains("epoch=J2000.0"));
        assert!(display_with_dist.contains("d=5"));
    }

    #[test]
    fn test_seasonal_classification() {
        let epoch = TT::j2000();

        let spring = EclipticPosition::from_degrees(45.0, 0.0, epoch).unwrap();
        assert_eq!(spring.season_index(), 0);

        let summer = EclipticPosition::from_degrees(135.0, 0.0, epoch).unwrap();
        assert_eq!(summer.season_index(), 1);

        let autumn = EclipticPosition::from_degrees(225.0, 0.0, epoch).unwrap();
        assert_eq!(autumn.season_index(), 2);

        let winter = EclipticPosition::from_degrees(315.0, 0.0, epoch).unwrap();
        assert_eq!(winter.season_index(), 3);
    }

    #[test]
    fn test_ecliptic_plane_classification() {
        let epoch = TT::j2000();

        let on_plane = EclipticPosition::from_degrees(45.0, 2.0, epoch).unwrap();
        assert!(on_plane.is_near_ecliptic_plane());
        assert!(!on_plane.is_near_ecliptic_pole());

        let off_plane = EclipticPosition::from_degrees(45.0, 45.0, epoch).unwrap();
        assert!(!off_plane.is_near_ecliptic_plane());

        let near_pole = EclipticPosition::from_degrees(0.0, 87.0, epoch).unwrap();
        assert!(near_pole.is_near_ecliptic_pole());
    }

    #[test]
    fn test_coordinate_edge_cases() {
        let epoch = TT::j2000();

        let wrapped_lambda = EclipticPosition::from_degrees(370.0, 0.0, epoch).unwrap();
        let expected_wrapped = Angle::from_degrees(370.0)
            .validate_longitude(true)
            .unwrap()
            .degrees();
        assert_eq!(wrapped_lambda.lambda().degrees(), expected_wrapped);

        let negative_lambda = EclipticPosition::from_degrees(-90.0, 0.0, epoch).unwrap();
        let expected_negative = Angle::from_degrees(-90.0)
            .validate_longitude(true)
            .unwrap()
            .degrees();
        assert_eq!(negative_lambda.lambda().degrees(), expected_negative);

        assert!(EclipticPosition::from_degrees(0.0, 95.0, epoch).is_err());
        assert!(EclipticPosition::from_degrees(0.0, -95.0, epoch).is_err());

        let max_beta = EclipticPosition::from_degrees(0.0, 90.0, epoch).unwrap();
        assert_eq!(max_beta.beta().degrees(), 90.0);

        let min_beta = EclipticPosition::from_degrees(0.0, -90.0, epoch).unwrap();
        assert_eq!(min_beta.beta().degrees(), -90.0);
    }

    #[test]
    fn test_pole_angular_separation_edge_cases() {
        let epoch = TT::j2000();

        let north_pole = EclipticPosition::north_ecliptic_pole(epoch);
        let south_pole = EclipticPosition::south_ecliptic_pole(epoch);

        let pole_separation = north_pole.angular_separation(&south_pole);
        assert_eq!(pole_separation.degrees(), 180.0);

        // Two points at the same pole with different longitudes should have zero separation.
        // At poles, longitude is undefined (singularity), so we test with identical coordinates.
        let same_pole = EclipticPosition::north_ecliptic_pole(epoch);
        let pole_separation_same = north_pole.angular_separation(&same_pole);
        assert_eq!(pole_separation_same.degrees(), 0.0);
    }

    #[test]
    fn test_pole_singularity_different_longitudes() {
        // At the poles, longitude is mathematically undefined (singularity).
        // Two points at beta=90° with different lambda values represent the same point.
        // Due to cos(90°) ≈ 6e-17 (not exactly 0), Vincenty formula returns tiny non-zero.
        // This test documents this known floating-point limitation.
        let epoch = TT::j2000();

        let north_pole = EclipticPosition::north_ecliptic_pole(epoch);
        let same_pole_diff_lon = EclipticPosition::from_degrees(123.456, 90.0, epoch).unwrap();

        let separation = north_pole.angular_separation(&same_pole_diff_lon);

        // The separation should be essentially zero (within floating-point noise)
        // cos(π/2) ≈ 6.12e-17, which propagates through Vincenty formula
        assert!(
            separation.degrees() < 1e-12,
            "Pole singularity: expected ~0, got {} degrees",
            separation.degrees()
        );
    }

    #[test]
    fn test_coordinate_transformations_at_poles() {
        let epoch = TT::j2000();

        let north_pole = EclipticPosition::north_ecliptic_pole(epoch);
        let icrs_north = north_pole.to_icrs(&epoch).unwrap();
        let roundtrip_north = EclipticPosition::from_icrs(&icrs_north, &epoch).unwrap();

        assert_eq!(roundtrip_north.beta().degrees(), 90.0);

        let south_pole = EclipticPosition::south_ecliptic_pole(epoch);
        let icrs_south = south_pole.to_icrs(&epoch).unwrap();
        let roundtrip_south = EclipticPosition::from_icrs(&icrs_south, &epoch).unwrap();

        assert_eq!(roundtrip_south.beta().degrees(), -90.0);
    }

    #[test]
    fn test_seasonal_boundary_cases() {
        let epoch = TT::j2000();

        let exactly_90 = EclipticPosition::from_degrees(90.0, 0.0, epoch).unwrap();
        assert_eq!(exactly_90.season_index(), 1);

        let exactly_180 = EclipticPosition::from_degrees(180.0, 0.0, epoch).unwrap();
        assert_eq!(exactly_180.season_index(), 2);

        let exactly_270 = EclipticPosition::from_degrees(270.0, 0.0, epoch).unwrap();
        assert_eq!(exactly_270.season_index(), 3);

        let exactly_0 = EclipticPosition::from_degrees(0.0, 0.0, epoch).unwrap();
        assert_eq!(exactly_0.season_index(), 0);

        let almost_360 = EclipticPosition::from_degrees(359.9, 0.0, epoch).unwrap();
        assert_eq!(almost_360.season_index(), 3);
    }

    #[test]
    fn test_plane_classification_boundary_cases() {
        let epoch = TT::j2000();

        let exactly_5_deg = EclipticPosition::from_degrees(0.0, 5.0, epoch).unwrap();
        assert!(!exactly_5_deg.is_near_ecliptic_plane());

        let just_under_5_deg = EclipticPosition::from_degrees(0.0, 4.99, epoch).unwrap();
        assert!(just_under_5_deg.is_near_ecliptic_plane());

        let exactly_85_deg = EclipticPosition::from_degrees(0.0, 85.0, epoch).unwrap();
        assert!(!exactly_85_deg.is_near_ecliptic_pole());

        let just_over_85_deg = EclipticPosition::from_degrees(0.0, 85.01, epoch).unwrap();
        assert!(just_over_85_deg.is_near_ecliptic_pole());

        let neg_85_deg = EclipticPosition::from_degrees(0.0, -85.01, epoch).unwrap();
        assert!(neg_85_deg.is_near_ecliptic_pole());
    }
}
