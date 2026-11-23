use crate::error::{WcsError, WcsResult};

use super::{Kind, Projection};

// What a projection derives follows from its parameters, so they alone
// decide equality; a derived NaN would otherwise make a projection unequal to
// itself.
impl PartialEq for Projection {
    fn eq(&self, other: &Self) -> bool {
        self.parameters() == other.parameters()
    }
}

impl Projection {
    // Paper II §5.1.7: ZPN, the projection with the most parameters, reads
    // PV2_0 through PV2_29.
    pub(crate) const PV_LAST_INDEX: u8 = 29;

    // Paper II §2.4: the projection the code names, with each parameter m read
    // from `pv` or, when that has none, given its default.
    pub(crate) fn from_code(code: &str, pv: impl Fn(u8) -> Option<f64>) -> WcsResult<Self> {
        let (first, defaults) = parameter_defaults(code);
        let values = (first..)
            .zip(defaults)
            .map(|(m, default)| {
                pv(m).or(*default).ok_or_else(|| {
                    WcsError::missing_keyword(format!("{code} projection requires PV2_{m}"))
                })
            })
            .collect::<WcsResult<Vec<_>>>()?;
        Self::with_parameters(code, &values).ok_or_else(|| WcsError::unsupported_projection(code))
    }

    pub(crate) fn code(&self) -> &'static str {
        self.parameters().0
    }

    // The parameters other than their defaults, as (m, PV2_m).
    pub(crate) fn pv_values(&self) -> Vec<(usize, f64)> {
        let (code, values) = self.parameters();
        let (first, defaults) = parameter_defaults(code);
        let default = |i: usize| defaults.get(i).copied().flatten();
        let differs = |&(i, value): &(usize, f64)| default(i) != Some(value);
        let at_m = |(i, value)| (usize::from(first) + i, value);
        values
            .into_iter()
            .enumerate()
            .filter(differs)
            .map(at_m)
            .collect()
    }

    fn parameters(&self) -> (&'static str, Vec<f64>) {
        match &self.0 {
            Kind::Tan => ("TAN", vec![]),
            Kind::Sin { xi, eta } => ("SIN", vec![*xi, *eta]),
            Kind::Arc => ("ARC", vec![]),
            Kind::Stg => ("STG", vec![]),
            Kind::Zea => ("ZEA", vec![]),
            Kind::Azp { mu, gamma, .. } => ("AZP", vec![*mu, *gamma]),
            Kind::Szp {
                mu, phi_c, theta_c, ..
            } => ("SZP", vec![*mu, *phi_c, *theta_c]),
            Kind::Zpn { coeffs } => ("ZPN", coeffs.clone()),
            Kind::Air { theta_b, .. } => ("AIR", vec![*theta_b]),
            Kind::Car => ("CAR", vec![]),
            Kind::Mer => ("MER", vec![]),
            Kind::Cea { lambda } => ("CEA", vec![*lambda]),
            Kind::Cyp { mu, lambda } => ("CYP", vec![*mu, *lambda]),
            Kind::Sfl => ("SFL", vec![]),
            Kind::Par => ("PAR", vec![]),
            Kind::Mol => ("MOL", vec![]),
            Kind::Ait => ("AIT", vec![]),
            Kind::Cop { theta_a, eta, .. } => ("COP", vec![*theta_a, *eta]),
            Kind::Coe { theta_a, eta, .. } => ("COE", vec![*theta_a, *eta]),
            Kind::Cod { theta_a, eta, .. } => ("COD", vec![*theta_a, *eta]),
            Kind::Coo { theta_a, eta, .. } => ("COO", vec![*theta_a, *eta]),
            Kind::Bon { theta_1, .. } => ("BON", vec![*theta_1]),
            Kind::Pco => ("PCO", vec![]),
            Kind::Tsc => ("TSC", vec![]),
            Kind::Csc => ("CSC", vec![]),
            Kind::Qsc => ("QSC", vec![]),
        }
    }

    fn with_parameters(code: &str, values: &[f64]) -> Option<Self> {
        let projection = match (code, values) {
            ("TAN", []) => Self::tan(),
            ("SIN", &[xi, eta]) => Self::sin_with_params(xi, eta),
            ("ARC", []) => Self::arc(),
            ("STG", []) => Self::stg(),
            ("ZEA", []) => Self::zea(),
            ("AZP", &[mu, gamma]) => Self::azp(mu, gamma),
            ("SZP", &[mu, phi_c, theta_c]) => Self::szp(mu, phi_c, theta_c),
            ("ZPN", coeffs) => Self::zpn(coeffs.to_vec()),
            ("AIR", &[theta_b]) => Self::air(theta_b),
            ("CAR", []) => Self::car(),
            ("MER", []) => Self::mer(),
            ("CEA", &[lambda]) => Self::cea_with_lambda(lambda),
            ("CYP", &[mu, lambda]) => Self::cyp(mu, lambda),
            ("SFL", []) => Self::sfl(),
            ("PAR", []) => Self::par(),
            ("MOL", []) => Self::mol(),
            ("AIT", []) => Self::ait(),
            ("COP", &[theta_a, eta]) => Self::cop(theta_a, eta),
            ("COE", &[theta_a, eta]) => Self::coe(theta_a, eta),
            ("COD", &[theta_a, eta]) => Self::cod(theta_a, eta),
            ("COO", &[theta_a, eta]) => Self::coo(theta_a, eta),
            ("BON", &[theta_1]) => Self::bon(theta_1),
            ("PCO", []) => Self::pco(),
            ("TSC", []) => Self::tsc(),
            ("CSC", []) => Self::csc(),
            ("QSC", []) => Self::qsc(),
            _ => return None,
        };
        Some(projection)
    }
}

// Paper II §5: the first m of the code's PV2_m and each parameter's default
// from there, where None marks one the projection requires.
fn parameter_defaults(code: &str) -> (u8, &'static [Option<f64>]) {
    match code {
        "SIN" | "AZP" => (1, &[Some(0.0), Some(0.0)]),
        "SZP" => (1, &[Some(0.0), Some(0.0), Some(90.0)]),
        "ZPN" => (0, &[Some(0.0); Projection::PV_LAST_INDEX as usize + 1]),
        "AIR" => (1, &[Some(90.0)]),
        "CEA" => (1, &[Some(1.0)]),
        "CYP" => (1, &[Some(1.0), Some(1.0)]),
        "COP" | "COE" | "COD" | "COO" => (1, &[None, Some(0.0)]),
        "BON" => (1, &[None]),
        _ => (1, &[]),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;

    fn every_code() -> [(&'static str, Projection); 26] {
        [
            ("TAN", Projection::tan()),
            ("SIN", Projection::sin_with_params(0.1, -0.2)),
            ("ARC", Projection::arc()),
            ("STG", Projection::stg()),
            ("ZEA", Projection::zea()),
            ("AZP", Projection::azp(2.0, 30.0)),
            ("SZP", Projection::szp(2.0, 45.0, 60.0)),
            ("ZPN", Projection::zpn(vec![0.0, 1.0, 0.0, 0.01])),
            ("AIR", Projection::air(45.0)),
            ("CAR", Projection::car()),
            ("MER", Projection::mer()),
            ("CEA", Projection::cea_with_lambda(0.5)),
            ("CYP", Projection::cyp(1.0, 0.5)),
            ("SFL", Projection::sfl()),
            ("PAR", Projection::par()),
            ("MOL", Projection::mol()),
            ("AIT", Projection::ait()),
            ("COP", Projection::cop(45.0, 15.0)),
            ("COE", Projection::coe(40.0, 10.0)),
            ("COD", Projection::cod(35.0, 5.0)),
            ("COO", Projection::coo(50.0, 20.0)),
            ("BON", Projection::bon(30.0)),
            ("PCO", Projection::pco()),
            ("TSC", Projection::tsc()),
            ("CSC", Projection::csc()),
            ("QSC", Projection::qsc()),
        ]
    }

    fn read_back(projection: &Projection) -> WcsResult<Projection> {
        let pv: HashMap<_, _> = projection.pv_values().into_iter().collect();
        Projection::from_code(projection.code(), |m| pv.get(&usize::from(m)).copied())
    }

    // What pv_values writes, from_code reads back as the same projection.
    #[test]
    fn test_every_code_round_trips_through_its_pv_values() {
        for (code, projection) in every_code() {
            assert_eq!(projection.code(), code);
            assert_eq!(read_back(&projection).unwrap(), projection, "{code}");
        }
    }

    // sin γ and cos γ of an infinite γ are NaN.
    #[test]
    fn test_equal_parameters_make_equal_projections_whatever_they_derive() {
        let azp = || Projection::azp(0.0, f64::INFINITY);
        assert_eq!(azp(), azp());
        assert_ne!(
            Projection::azp(0.0, f64::NAN),
            Projection::azp(0.0, f64::NAN)
        );
    }

    // Each parameter goes out under its own m, and one at its Paper II
    // default does not go out at all.
    #[test]
    fn test_pv_values_omit_defaults() {
        let zpn = Projection::zpn(vec![0.0, 1.0, 0.0, 0.01]);
        let cases: [(Projection, &[(usize, f64)]); 17] = [
            (Projection::tan(), &[]),
            (
                Projection::sin_with_params(0.1, -0.2),
                &[(1, 0.1), (2, -0.2)],
            ),
            (Projection::sin(), &[]),
            (Projection::azp(2.0, 30.0), &[(1, 2.0), (2, 30.0)]),
            (Projection::azp(0.0, 0.0), &[]),
            (
                Projection::szp(2.0, 45.0, 60.0),
                &[(1, 2.0), (2, 45.0), (3, 60.0)],
            ),
            (Projection::szp(0.0, 0.0, 90.0), &[]),
            (zpn, &[(1, 1.0), (3, 0.01)]),
            (Projection::air(45.0), &[(1, 45.0)]),
            (Projection::air(90.0), &[]),
            (Projection::cea(), &[]),
            (Projection::cea_with_lambda(0.5), &[(1, 0.5)]),
            (Projection::cyp(1.0, 0.5), &[(2, 0.5)]),
            (Projection::cyp(1.0, 1.0), &[]),
            (Projection::cop(45.0, 15.0), &[(1, 45.0), (2, 15.0)]),
            (Projection::coo(50.0, 0.0), &[(1, 50.0)]),
            (Projection::bon(30.0), &[(1, 30.0)]),
        ];
        for (projection, expected) in cases {
            assert_eq!(projection.pv_values(), expected, "{projection:?}");
        }
    }

    #[test]
    fn test_from_code_gives_absent_parameters_their_defaults() {
        let none = |_| None;
        let cases = [
            ("SIN", Projection::sin()),
            ("AZP", Projection::azp(0.0, 0.0)),
            ("SZP", Projection::szp(0.0, 0.0, 90.0)),
            ("AIR", Projection::air(90.0)),
            ("CEA", Projection::cea()),
            // Paper II §5.2.1 gives µ and λ the default value 1.
            ("CYP", Projection::cyp(1.0, 1.0)),
            // Paper II Table 13: PV2_m default to 0, so ZPN without them is
            // not ARC.
            ("ZPN", Projection::zpn(vec![])),
        ];
        for (code, expected) in cases {
            assert_eq!(Projection::from_code(code, none).unwrap(), expected);
        }
        let eta_absent = Projection::from_code("COD", |m| (m == 1).then_some(30.0));
        assert_eq!(eta_absent.unwrap(), Projection::cod(30.0, 0.0));
    }

    // Paper II §5.1.7: ZPN reads PV2_0 through PV2_29.
    #[test]
    fn test_from_code_reads_zpn_coefficients_up_to_29() {
        let pv = HashMap::from([(0, 0.0), (1, 1.0), (29, 1e-9), (30, 1.0)]);
        let projection = Projection::from_code("ZPN", |m| pv.get(&m).copied()).unwrap();
        let mut coeffs = vec![0.0; 30];
        coeffs[1] = 1.0;
        coeffs[29] = 1e-9;
        assert_eq!(projection, Projection::zpn(coeffs));
    }

    #[test]
    fn test_from_code_requires_the_conic_and_bon_parameter() {
        for code in ["COP", "COE", "COD", "COO", "BON"] {
            let err = Projection::from_code(code, |m| (m == 2).then_some(5.0)).unwrap_err();
            let expected = format!("{code} projection requires PV2_1");
            assert!(err.to_string().contains(&expected), "{code}: {err}");
        }
    }

    #[test]
    fn test_from_code_rejects_an_unknown_code() {
        for code in ["XYZ", "tan", ""] {
            let err = Projection::from_code(code, |_| None).unwrap_err();
            assert!(
                matches!(err, WcsError::UnsupportedProjection { .. }),
                "{code}: {err}"
            );
        }
    }
}
