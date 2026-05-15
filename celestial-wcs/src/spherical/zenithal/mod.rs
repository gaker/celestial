mod tan;
mod sin;
mod arc;
mod stg;
mod zea;
mod azp;
mod szp;
mod zpn;
mod air;

pub(super) use tan::{project_tan, deproject_tan};
pub(super) use sin::{project_sin, deproject_sin};
pub(super) use arc::{project_arc, deproject_arc};
pub(super) use stg::{project_stg, deproject_stg};
pub(super) use zea::{project_zea, deproject_zea};
pub(super) use azp::{project_azp, deproject_azp};
pub(super) use szp::{project_szp, deproject_szp};
pub(super) use zpn::{project_zpn, deproject_zpn};
pub(super) use air::{project_air, deproject_air};


#[cfg(test)]
mod tests {
    use crate::coordinate::IntermediateCoord;
    use crate::Projection;

    #[test]
    fn test_deproject_origin_returns_pole() {
        let origin = IntermediateCoord::new(0.0, 0.0);

        let tan_result = Projection::tan().deproject(origin).unwrap();
        assert_eq!(tan_result.phi().degrees(), 0.0);
        assert_eq!(tan_result.theta().degrees(), 90.0);

        let arc_result = Projection::arc().deproject(origin).unwrap();
        assert_eq!(arc_result.phi().degrees(), 0.0);
        assert_eq!(arc_result.theta().degrees(), 90.0);

        let stg_result = Projection::stg().deproject(origin).unwrap();
        assert_eq!(stg_result.phi().degrees(), 0.0);
        assert_eq!(stg_result.theta().degrees(), 90.0);

        let zea_result = Projection::zea().deproject(origin).unwrap();
        assert_eq!(zea_result.phi().degrees(), 0.0);
        assert_eq!(zea_result.theta().degrees(), 90.0);

        let azp_result = Projection::azp(2.0, 0.0).deproject(origin).unwrap();
        assert_eq!(azp_result.phi().degrees(), 0.0);
        assert_eq!(azp_result.theta().degrees(), 90.0);
    }

    #[test]
    fn test_all_projections_native_reference() {
        let projections = [
            Projection::tan(),
            Projection::sin(),
            Projection::arc(),
            Projection::stg(),
            Projection::zea(),
            Projection::azp(2.0, 0.0),
        ];

        for proj in projections {
            let (phi0, theta0) = proj.native_reference();
            assert_eq!(phi0, 0.0);
            assert_eq!(theta0, 90.0);
        }
    }

}
