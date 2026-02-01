use crate::{frames::ITRSPosition, CoordResult};
use astro_core::constants::ARCSEC_TO_RAD;
use astro_core::Vector3;
use astro_time::{scales::conversions::ToUT1WithDeltaT, transforms::earth_rotation_angle, TT};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TIRSPosition {
    x: f64,
    y: f64,
    z: f64,
    epoch: TT,
}

impl TIRSPosition {
    pub fn new(x: f64, y: f64, z: f64, epoch: TT) -> Self {
        Self { x, y, z, epoch }
    }

    pub fn x(&self) -> f64 {
        self.x
    }

    pub fn y(&self) -> f64 {
        self.y
    }

    pub fn z(&self) -> f64 {
        self.z
    }

    pub fn epoch(&self) -> TT {
        self.epoch
    }

    pub fn position_vector(&self) -> Vector3 {
        Vector3::new(self.x, self.y, self.z)
    }

    pub fn from_position_vector(pos: Vector3, epoch: TT) -> Self {
        Self::new(pos.x, pos.y, pos.z, epoch)
    }

    pub fn geocentric_distance(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn distance_to(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    fn polar_motion_matrix(xp: f64, yp: f64, sp: f64) -> [[f64; 3]; 3] {
        let mut matrix = [[0.0; 3]; 3];
        matrix[0][0] = 1.0;
        matrix[1][1] = 1.0;
        matrix[2][2] = 1.0;

        Self::apply_rotation_z(&mut matrix, sp);
        Self::apply_rotation_y(&mut matrix, -xp);
        Self::apply_rotation_x(&mut matrix, -yp);

        matrix
    }

    fn apply_rotation_z(matrix: &mut [[f64; 3]; 3], angle: f64) {
        let (s, c) = angle.sin_cos();
        let temp = *matrix;

        for j in 0..3 {
            matrix[0][j] = c * temp[0][j] + s * temp[1][j];
            matrix[1][j] = -s * temp[0][j] + c * temp[1][j];
        }
    }

    fn apply_rotation_y(matrix: &mut [[f64; 3]; 3], angle: f64) {
        let (s, c) = angle.sin_cos();
        let temp = *matrix;

        for j in 0..3 {
            matrix[0][j] = c * temp[0][j] - s * temp[2][j];
            matrix[2][j] = s * temp[0][j] + c * temp[2][j];
        }
    }

    fn apply_rotation_x(matrix: &mut [[f64; 3]; 3], angle: f64) {
        let (s, c) = angle.sin_cos();
        let temp = *matrix;

        for j in 0..3 {
            matrix[1][j] = c * temp[1][j] + s * temp[2][j];
            matrix[2][j] = -s * temp[1][j] + c * temp[2][j];
        }
    }

    /// Computes ΔT (TT - UT1) using EOP parameters.
    ///
    /// # EOP Freshness Check
    /// Rejects EOP data >1 day from the target epoch. This ensures UT1-UTC is current, as
    /// Earth's rotation is irregular. For sparse datasets, consider:
    /// - Using an interpolator to fill gaps (see `EopManager`)
    /// - Relaxing this check if lower precision is acceptable (modify this threshold)
    /// - Pre-fetching/caching EOP data for your observation window
    ///
    /// Typical use cases handle this via interpolation, so sparse raw data should be rare.
    pub fn compute_delta_t(epoch: &TT, eop: &crate::eop::EopParameters) -> CoordResult<f64> {
        use astro_time::scales::conversions::{ToTAI, ToUTC};

        let epoch_jd = epoch.to_julian_date().to_f64();
        let eop_jd = eop.mjd + astro_core::constants::MJD_ZERO_POINT;

        if (epoch_jd - eop_jd).abs() > 1.0 {
            return Err(crate::CoordError::data_unavailable(
                "EOP data is more than 1 day from epoch - Delta-T computation may be inaccurate",
            ));
        }

        let tai = epoch
            .to_tai()
            .map_err(|e| crate::CoordError::external_library("TT to TAI", &e.to_string()))?;
        let utc = tai
            .to_utc()
            .map_err(|e| crate::CoordError::external_library("TAI to UTC", &e.to_string()))?;

        let utc_jd = utc.to_julian_date().to_f64();
        let ut1_jd = utc_jd + eop.ut1_utc / astro_core::constants::SECONDS_PER_DAY_F64;

        let delta_t_days = epoch_jd - ut1_jd;
        let delta_t_seconds = delta_t_days * astro_core::constants::SECONDS_PER_DAY_F64;

        Ok(delta_t_seconds)
    }

    pub fn to_itrs(
        &self,
        epoch: &TT,
        eop: &crate::eop::EopParameters,
    ) -> CoordResult<ITRSPosition> {
        let delta_t_seconds = Self::compute_delta_t(epoch, eop)?;
        let ut1 = epoch.to_ut1_with_delta_t(delta_t_seconds)?;
        let era = earth_rotation_angle(&ut1.to_julian_date())?;

        let (sin_era, cos_era) = era.sin_cos();

        let x_after_era = cos_era * self.x + sin_era * self.y;
        let y_after_era = -sin_era * self.x + cos_era * self.y;
        let z_after_era = self.z;

        let xp_rad = eop.x_p * ARCSEC_TO_RAD;
        let yp_rad = eop.y_p * ARCSEC_TO_RAD;

        let polar_motion_matrix = Self::polar_motion_matrix(xp_rad, yp_rad, eop.s_prime);

        let x_itrs = polar_motion_matrix[0][0] * x_after_era
            + polar_motion_matrix[0][1] * y_after_era
            + polar_motion_matrix[0][2] * z_after_era;
        let y_itrs = polar_motion_matrix[1][0] * x_after_era
            + polar_motion_matrix[1][1] * y_after_era
            + polar_motion_matrix[1][2] * z_after_era;
        let z_itrs = polar_motion_matrix[2][0] * x_after_era
            + polar_motion_matrix[2][1] * y_after_era
            + polar_motion_matrix[2][2] * z_after_era;

        Ok(ITRSPosition::new(x_itrs, y_itrs, z_itrs, *epoch))
    }

    pub fn from_itrs(
        itrs: &ITRSPosition,
        epoch: &TT,
        eop: &crate::eop::EopParameters,
    ) -> CoordResult<Self> {
        let xp_rad = eop.x_p * ARCSEC_TO_RAD;
        let yp_rad = eop.y_p * ARCSEC_TO_RAD;

        let polar_motion_matrix = Self::polar_motion_matrix(xp_rad, yp_rad, eop.s_prime);

        let x_before_pm = polar_motion_matrix[0][0] * itrs.x()
            + polar_motion_matrix[1][0] * itrs.y()
            + polar_motion_matrix[2][0] * itrs.z();
        let y_before_pm = polar_motion_matrix[0][1] * itrs.x()
            + polar_motion_matrix[1][1] * itrs.y()
            + polar_motion_matrix[2][1] * itrs.z();
        let z_before_pm = polar_motion_matrix[0][2] * itrs.x()
            + polar_motion_matrix[1][2] * itrs.y()
            + polar_motion_matrix[2][2] * itrs.z();

        let delta_t_seconds = Self::compute_delta_t(epoch, eop)?;
        let ut1 = epoch.to_ut1_with_delta_t(delta_t_seconds)?;
        let era = earth_rotation_angle(&ut1.to_julian_date())?;

        let (sin_era, cos_era) = era.sin_cos();

        let x_tirs = cos_era * x_before_pm - sin_era * y_before_pm;
        let y_tirs = sin_era * x_before_pm + cos_era * y_before_pm;
        let z_tirs = z_before_pm;

        Ok(Self::new(x_tirs, y_tirs, z_tirs, *epoch))
    }

    pub fn from_cirs(
        cirs_vec: Vector3,
        epoch: &TT,
        eop: &crate::eop::EopParameters,
    ) -> CoordResult<Self> {
        let delta_t_seconds = Self::compute_delta_t(epoch, eop)?;
        let ut1 = epoch.to_ut1_with_delta_t(delta_t_seconds)?;
        let era = earth_rotation_angle(&ut1.to_julian_date())?;

        let (sin_era, cos_era) = era.sin_cos();

        let x_tirs = cos_era * cirs_vec.x - sin_era * cirs_vec.y;
        let y_tirs = sin_era * cirs_vec.x + cos_era * cirs_vec.y;
        let z_tirs = cirs_vec.z;

        Ok(Self::new(x_tirs, y_tirs, z_tirs, *epoch))
    }

    pub fn to_cirs(
        &self,
        eop: &crate::eop::EopParameters,
    ) -> CoordResult<crate::frames::CIRSPosition> {
        use crate::frames::CIRSPosition;

        let delta_t_seconds = Self::compute_delta_t(&self.epoch, eop)?;
        let ut1 = self.epoch.to_ut1_with_delta_t(delta_t_seconds)?;
        let era = earth_rotation_angle(&ut1.to_julian_date())?;

        let (sin_era, cos_era) = era.sin_cos();

        let x_cirs = cos_era * self.x + sin_era * self.y;
        let y_cirs = -sin_era * self.x + cos_era * self.y;
        let z_cirs = self.z;

        let cirs_vec = Vector3::new(x_cirs, y_cirs, z_cirs);
        CIRSPosition::from_unit_vector(cirs_vec, self.epoch)
    }
}

impl std::fmt::Display for TIRSPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TIRS(X={:.3}m, Y={:.3}m, Z={:.3}m, epoch=J{:.1})",
            self.x,
            self.y,
            self.z,
            self.epoch.julian_year()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use astro_core::constants::TWOPI;

    #[test]
    fn test_tirs_creation() {
        let epoch = TT::j2000();
        let pos = TIRSPosition::new(1000000.0, 2000000.0, 3000000.0, epoch);

        assert_eq!(pos.x(), 1000000.0);
        assert_eq!(pos.y(), 2000000.0);
        assert_eq!(pos.z(), 3000000.0);
        assert_eq!(pos.epoch(), epoch);
    }

    #[test]
    fn test_vector_operations() {
        let epoch = TT::j2000();
        let original = TIRSPosition::new(1000.0, 2000.0, 3000.0, epoch);

        let vec = original.position_vector();
        assert_eq!(vec.x, 1000.0);
        assert_eq!(vec.y, 2000.0);
        assert_eq!(vec.z, 3000.0);

        let recovered = TIRSPosition::from_position_vector(vec, epoch);
        assert_eq!(recovered.x(), original.x());
        assert_eq!(recovered.y(), original.y());
        assert_eq!(recovered.z(), original.z());
    }

    #[test]
    fn test_distance_calculations() {
        let epoch = TT::j2000();

        let pos1 = TIRSPosition::new(1000.0, 0.0, 0.0, epoch);
        let pos2 = TIRSPosition::new(2000.0, 0.0, 0.0, epoch);

        // Distance between positions
        assert_eq!(pos1.distance_to(&pos2), 1000.0);
        assert_eq!(pos2.distance_to(&pos1), 1000.0);

        // Distance from origin
        assert_eq!(pos1.geocentric_distance(), 1000.0);
        assert_eq!(pos2.geocentric_distance(), 2000.0);
    }

    #[test]
    fn test_itrs_transformation_roundtrip() {
        use crate::eop::EopRecord;

        let epoch = TT::j2000();
        let original_tirs = TIRSPosition::new(4000000.0, 3000000.0, 5000000.0, epoch);

        let eop = EopRecord::new(51544.5, 0.0, 0.0, 0.3, 0.0)
            .unwrap()
            .to_parameters();

        // Transform to ITRS and back
        let itrs = original_tirs.to_itrs(&epoch, &eop).unwrap();
        let recovered_tirs = TIRSPosition::from_itrs(&itrs, &epoch, &eop).unwrap();

        // Should roundtrip with floating-point precision (rotation with EOP-based Delta-T)
        // Allow 2 ULP due to additional operations in Delta-T calculation
        astro_core::test_helpers::assert_ulp_le(
            recovered_tirs.x(),
            original_tirs.x(),
            2,
            "X roundtrip",
        );
        astro_core::test_helpers::assert_ulp_le(
            recovered_tirs.y(),
            original_tirs.y(),
            2,
            "Y roundtrip",
        );
        astro_core::test_helpers::assert_ulp_le(
            recovered_tirs.z(),
            original_tirs.z(),
            2,
            "Z roundtrip",
        );
    }

    #[test]
    fn test_itrs_transformation_properties() {
        use crate::eop::EopRecord;

        let epoch = TT::j2000();
        let tirs = TIRSPosition::new(1000000.0, 0.0, 0.0, epoch);

        let eop = EopRecord::new(51544.5, 0.0, 0.0, 0.3, 0.0)
            .unwrap()
            .to_parameters();

        let itrs = tirs.to_itrs(&epoch, &eop).unwrap();

        // Z coordinate should be unchanged (rotation around Z-axis)
        assert_eq!(itrs.z(), tirs.z());

        // Distance from origin should be preserved
        astro_core::test_helpers::assert_ulp_le(
            itrs.geocentric_distance(),
            tirs.geocentric_distance(),
            1,
            "Distance preservation",
        );

        // X and Y should change due to Earth rotation (unless ERA = 0)
        let delta_t_seconds = 69.184;
        let ut1 = epoch.to_ut1_with_delta_t(delta_t_seconds).unwrap();
        let era = earth_rotation_angle(&ut1.to_julian_date()).unwrap();
        if era != 0.0 {
            // Only test if ERA is non-zero (it should be for J2000)
            assert!(itrs.x() != tirs.x() || itrs.y() != tirs.y());
        }
    }

    #[test]
    fn test_display_formatting() {
        let epoch = TT::j2000();
        let pos = TIRSPosition::new(1234567.89, -987654.32, 555666.77, epoch);

        let display = format!("{}", pos);
        assert!(display.contains("TIRS"));
        assert!(display.contains("1234567.890m"));
        assert!(display.contains("-987654.320m"));
        assert!(display.contains("555666.770m"));
        assert!(display.contains("J2000.0"));
    }

    #[test]
    fn test_itrs_consistency() {
        use crate::eop::EopRecord;

        let epoch = TT::j2000();
        let itrs_original = ITRSPosition::new(2000000.0, 1000000.0, 6000000.0, epoch);

        let eop = EopRecord::new(51544.5, 0.0, 0.0, 0.3, 0.0)
            .unwrap()
            .to_parameters();

        // Transform ITRS -> TIRS -> ITRS
        let tirs = TIRSPosition::from_itrs(&itrs_original, &epoch, &eop).unwrap();
        let itrs_recovered = tirs.to_itrs(&epoch, &eop).unwrap();

        // Should roundtrip with floating-point precision
        astro_core::test_helpers::assert_ulp_le(
            itrs_recovered.x(),
            itrs_original.x(),
            1,
            "ITRS X roundtrip",
        );
        astro_core::test_helpers::assert_ulp_le(
            itrs_recovered.y(),
            itrs_original.y(),
            1,
            "ITRS Y roundtrip",
        );
        astro_core::test_helpers::assert_ulp_le(
            itrs_recovered.z(),
            itrs_original.z(),
            1,
            "ITRS Z roundtrip",
        );
    }

    #[test]
    fn test_earth_rotation_angle_usage() {
        use crate::eop::EopRecord;

        let epoch = TT::j2000();

        let eop = EopRecord::new(51544.5, 0.0, 0.0, 0.3, 0.0)
            .unwrap()
            .to_parameters();

        let delta_t_seconds = TIRSPosition::compute_delta_t(&epoch, &eop).unwrap();
        let ut1 = epoch.to_ut1_with_delta_t(delta_t_seconds).unwrap();
        let era = earth_rotation_angle(&ut1.to_julian_date()).unwrap();

        assert!(era >= 0.0);
        assert!(era < TWOPI);

        let tirs_x_axis = TIRSPosition::new(6378137.0, 0.0, 0.0, epoch);
        let itrs = tirs_x_axis.to_itrs(&epoch, &eop).unwrap();

        let expected_x = 6378137.0 * era.cos();
        let expected_y = -6378137.0 * era.sin();

        astro_core::test_helpers::assert_ulp_le(itrs.x(), expected_x, 1, "ERA X transformation");
        astro_core::test_helpers::assert_ulp_le(itrs.y(), expected_y, 1, "ERA Y transformation");
        assert_eq!(itrs.z(), 0.0);
    }

    #[test]
    fn test_eop_timing_check() {
        use crate::eop::EopRecord;

        let epoch = TT::j2000();
        let eop_old = EopRecord::new(51542.0, 0.0, 0.0, 0.3, 0.0)
            .unwrap()
            .to_parameters();

        let result = TIRSPosition::compute_delta_t(&epoch, &eop_old);
        assert!(result.is_err());
    }

    #[test]
    fn test_cirs_transformation() {
        use crate::eop::EopRecord;

        let epoch = TT::j2000();
        let eop = EopRecord::new(51544.5, 0.0, 0.0, 0.3, 0.0)
            .unwrap()
            .to_parameters();

        let tirs = TIRSPosition::new(6378137.0, 0.0, 0.0, epoch);

        let cirs = tirs.to_cirs(&eop).unwrap();

        assert!(cirs.ra().degrees() >= 0.0);
        assert!(cirs.ra().degrees() < 360.0);
        assert!(cirs.dec().degrees() >= -90.0);
        assert!(cirs.dec().degrees() <= 90.0);
    }
}
