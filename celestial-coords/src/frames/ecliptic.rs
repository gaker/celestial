use super::direction::spherical_angles;
use crate::distance::Distance;
use crate::errors::CoordResult;
use crate::frames::icrs::ICRSPosition;
use crate::transforms::CoordinateFrame;
use celestial_core::angle::Angle;
use celestial_core::matrix::{RotationMatrix3, Vector3};
use celestial_time::scales::tt::TT;
use celestial_time::transforms::precession::PrecessionCalculator;

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
        let lambda = lambda.normalized()?;
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

    pub fn mean_obliquity(&self) -> CoordResult<Angle> {
        let jd = self.epoch.to_julian_date();
        let obliquity = celestial_core::obliquity::iau_2006_mean_obliquity(jd.jd1(), jd.jd2())?;
        Ok(Angle::from_radians(obliquity))
    }

    pub fn true_obliquity(&self) -> CoordResult<Angle> {
        use celestial_time::transforms::nutation::NutationCalculator;

        let nutation = self.epoch.nutation_iau2006a()?;

        let true_obliquity = self.mean_obliquity()?.radians() + nutation.nutation_obliquity();

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
        Angle::from_radians(celestial_core::math::angular_separation(
            self.lambda.radians(),
            self.beta.radians(),
            other.lambda.radians(),
            other.beta.radians(),
        ))
    }
}

pub(crate) fn ecm06_matrix(epoch: &TT) -> CoordResult<RotationMatrix3> {
    let precession = epoch.precession()?;
    let bias_precession_matrix = precession.bias_precession_matrix;

    let jd = epoch.to_julian_date();
    let mean_obliquity = celestial_core::obliquity::iau_2006_mean_obliquity(jd.jd1(), jd.jd2())?;

    let mut ecliptic_rotation = RotationMatrix3::identity();
    ecliptic_rotation.rotate_x(mean_obliquity);

    Ok(ecliptic_rotation.multiply(&bias_precession_matrix))
}

impl CoordinateFrame for EclipticPosition {
    fn to_icrs(&self, _epoch: &TT) -> CoordResult<ICRSPosition> {
        let ecliptic = Vector3::from_spherical(self.lambda.radians(), self.beta.radians());
        let (ra, dec) = spherical_angles(ecm06_matrix(&self.epoch)?.transpose() * ecliptic)?;
        let mut icrs = ICRSPosition::new(ra, dec)?;
        if let Some(distance) = self.distance {
            icrs.set_distance(distance);
        }
        Ok(icrs)
    }

    fn from_icrs(icrs: &ICRSPosition, epoch: &TT) -> CoordResult<Self> {
        let direction = Vector3::from_spherical(icrs.ra().radians(), icrs.dec().radians());
        let (lambda, beta) = spherical_angles(ecm06_matrix(epoch)? * direction)?;
        let mut ecliptic = Self::new(lambda, beta, *epoch)?;
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
    use crate::distance::Distance;

    // Expected values are ERFA's (eraObl06, eraEcm06, eraEceq06).
    mod erfa_reference {
        use super::*;
        use celestial_core::constants::{HALF_PI, PI};

        #[test]
        fn test_obliquity_at_j2000() {
            let epoch = TT::j2000();
            let pos = EclipticPosition::from_degrees(0.0, 0.0, epoch).unwrap();
            let mean_obliquity = pos.mean_obliquity().unwrap();
            assert_eq!(mean_obliquity.radians(), 0.4090926006005829);
        }

        #[test]
        fn test_ecm06_matrix_at_j2000() {
            let epoch = TT::j2000();
            let matrix = ecm06_matrix(&epoch).unwrap();
            let m = matrix.elements();

            assert_eq!(m[0][0], 0.9999999999999941);
            assert_eq!(m[0][1], -7.078368960971556e-8);
            assert_eq!(m[0][2], 8.056213977613186e-8);
            assert_eq!(m[1][0], 3.2897004077419646e-8);
            assert_eq!(m[1][1], 0.9174821299149584);
            assert_eq!(m[1][2], 0.39777699944404793);
            assert_eq!(m[2][0], -1.0207044725484355e-7);
            assert_eq!(m[2][1], -0.39777699944404304);
            assert_eq!(m[2][2], 0.9174821299149556);
        }

        #[test]
        fn test_north_ecliptic_pole_to_icrs() {
            let epoch = TT::j2000();
            let north_pole = EclipticPosition::north_ecliptic_pole(epoch);
            let icrs = north_pole.to_icrs(&epoch).unwrap();

            assert_eq!(icrs.ra().radians(), 4.712388723782505);
            assert_eq!(icrs.dec().radians(), 1.1617036931348688);
        }

        #[test]
        fn test_south_ecliptic_pole_to_icrs() {
            let epoch = TT::j2000();
            let south_pole = EclipticPosition::south_ecliptic_pole(epoch);
            let icrs = south_pole.to_icrs(&epoch).unwrap();

            assert_eq!(icrs.ra().radians(), 1.5707960701927113);
            assert_eq!(icrs.dec().radians(), -1.1617036931348688);
        }

        // eraEceq06 then eraEqec06 at J2000: (lambda, beta) in degrees, the ICRS place
        // (ra, dec) and the ecliptic place recovered from it, in radians.
        const ROUND_TRIPS: [(f64, f64, [f64; 2], [f64; 2]); 8] = [
            (
                123.456789,
                45.678901,
                [2.565485497265132, 1.0935584854954725],
                [2.1547274519899178, 0.7972472211425305],
            ),
            (
                267.314159,
                -23.271828,
                [4.649606254950768, -0.8146774167786783],
                [4.6655122117496335, -0.4061700215578071],
            ),
            (
                45.123456,
                67.890123,
                [5.846910916176382, 1.2734169633587697],
                [0.7875528770787902, 1.1849061759339308],
            ),
            (
                90.0,
                0.0,
                [1.5707962909391529, 0.40909263366001886],
                [HALF_PI, 0.0],
            ),
            (
                180.0,
                0.0,
                [3.1415925828061035, -8.056213972741833e-8],
                [PI, 9.248752751069707e-18],
            ),
            (
                270.0,
                0.0,
                [4.7123889445289455, -0.40909263366001886],
                [4.712388980384689, 0.0],
            ),
            // The equinox comes back just below 2 pi, not at 0.
            (
                0.0,
                0.0,
                [6.283185236395896, 8.056213977613197e-8],
                [6.283185307179585, 2.608072413413319e-16],
            ),
            (
                0.0,
                90.0,
                [4.712388723782505, 1.1617036931348688],
                [5.423206103920281, HALF_PI],
            ),
        ];

        #[test]
        fn test_round_trips() {
            let epoch = TT::j2000();
            for (lambda, beta, icrs, back) in ROUND_TRIPS {
                let original = EclipticPosition::from_degrees(lambda, beta, epoch).unwrap();
                let place = original.to_icrs(&epoch).unwrap();
                let radians = [place.ra().radians(), place.dec().radians()];
                assert_eq!(radians, icrs, "({lambda}, {beta})");
                let recovered = EclipticPosition::from_icrs(&place, &epoch).unwrap();
                let radians = [recovered.lambda().radians(), recovered.beta().radians()];
                assert_eq!(radians, back, "({lambda}, {beta})");
            }
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

        let mean_obliquity = pos.mean_obliquity().unwrap();
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
        let expected_wrapped = Angle::from_degrees(370.0).normalized().unwrap().degrees();
        assert_eq!(wrapped_lambda.lambda().degrees(), expected_wrapped);

        let negative_lambda = EclipticPosition::from_degrees(-90.0, 0.0, epoch).unwrap();
        let expected_negative = Angle::from_degrees(-90.0).normalized().unwrap().degrees();
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
        // The cosine of the double nearest pi/2 is 6.1e-17, not 0, so two poles with different
        // longitudes sit a hair apart; eraSeps gives the same.
        let epoch = TT::j2000();

        let north_pole = EclipticPosition::north_ecliptic_pole(epoch);
        let same_pole_diff_lon = EclipticPosition::from_degrees(123.456, 90.0, epoch).unwrap();

        let separation = north_pole.angular_separation(&same_pole_diff_lon);
        assert_eq!(separation.radians(), 1.0785573740246409e-16);
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

    #[test]
    fn test_nutation_failure_is_an_epoch_error() {
        let epoch = TT::from_julian_date(celestial_time::julian::JulianDate::new(f64::NAN, 0.0));
        let pos = EclipticPosition::new(Angle::ZERO, Angle::ZERO, epoch).unwrap();
        let result = pos.true_obliquity();
        assert!(
            matches!(result, Err(crate::errors::CoordError::EpochError(_))),
            "{result:?}"
        );
    }
}
