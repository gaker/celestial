mod air;
mod arc;
mod azp;
mod sin;
mod stg;
mod szp;
mod tan;
mod zea;
mod zpn;

pub(super) use air::{air_k, deproject_air, project_air};
pub(super) use arc::{deproject_arc, project_arc};
pub(super) use azp::{deproject_azp, project_azp};
pub(super) use sin::{deproject_sin, project_sin};
pub(super) use stg::{deproject_stg, project_stg};
pub(super) use szp::{deproject_szp, project_szp, projection_point};
pub(super) use tan::{deproject_tan, project_tan};
pub(super) use zea::{deproject_zea, project_zea};
pub(super) use zpn::{deproject_zpn, project_zpn};

#[cfg(test)]
mod tests {
    use crate::coordinate::IntermediateCoord;
    use crate::spherical::Projection;

    #[test]
    fn test_deproject_origin_returns_pole() {
        let origin = IntermediateCoord::new(0.0, 0.0);

        let tan_result = Projection::tan().deproject(origin).unwrap();
        assert_eq!(tan_result.phi_deg(), 0.0);
        assert_eq!(tan_result.theta_deg(), 90.0);

        let arc_result = Projection::arc().deproject(origin).unwrap();
        assert_eq!(arc_result.phi_deg(), 0.0);
        assert_eq!(arc_result.theta_deg(), 90.0);

        let stg_result = Projection::stg().deproject(origin).unwrap();
        assert_eq!(stg_result.phi_deg(), 0.0);
        assert_eq!(stg_result.theta_deg(), 90.0);

        let zea_result = Projection::zea().deproject(origin).unwrap();
        assert_eq!(zea_result.phi_deg(), 0.0);
        assert_eq!(zea_result.theta_deg(), 90.0);

        let azp_result = Projection::azp(2.0, 0.0).deproject(origin).unwrap();
        assert_eq!(azp_result.phi_deg(), 0.0);
        assert_eq!(azp_result.theta_deg(), 90.0);

        for sin in [Projection::sin(), Projection::sin_with_params(0.1, -0.2)] {
            let sin_result = sin.deproject(origin).unwrap();
            assert_eq!(sin_result.phi_deg(), 0.0, "{sin:?}");
            assert_eq!(sin_result.theta_deg(), 90.0, "{sin:?}");
        }
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
