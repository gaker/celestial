use super::{GAST, GMST, LAST, LMST};

impl From<LMST> for GMST {
    fn from(lmst: LMST) -> GMST {
        lmst.to_gmst()
    }
}

impl From<LAST> for GAST {
    fn from(last: LAST) -> GAST {
        last.to_gast()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use astro_core::Location;

    fn mauna_kea() -> Location {
        Location::from_degrees(19.8283, -155.4783, 4145.0).unwrap()
    }

    #[test]
    fn test_lmst_to_gmst_conversion() {
        let location = mauna_kea();
        let lmst = LMST::from_hours(12.0, &location);
        let gmst: GMST = lmst.into();

        let expected_hours = 12.0 - (-155.4783 / 15.0);
        assert!((gmst.hours() - expected_hours).abs() < 1e-12);
    }

    #[test]
    fn test_last_to_gast_conversion() {
        let location = mauna_kea();
        let last = LAST::from_hours(12.0, &location);
        let gast: GAST = last.into();

        let expected_hours = 12.0 - (-155.4783 / 15.0);
        assert!((gast.hours() - expected_hours).abs() < 1e-12);
    }

    #[test]
    fn test_gmst_to_lmst_conversion() {
        let location = mauna_kea();
        let gmst = GMST::from_hours(12.0);
        let lmst = gmst.to_lmst(&location);

        let expected_hours = 12.0 + (-155.4783 / 15.0);
        let expected_normalized = if expected_hours < 0.0 {
            expected_hours + 24.0
        } else {
            expected_hours
        };
        assert!((lmst.hours() - expected_normalized).abs() < 1e-12);
    }

    #[test]
    fn test_gast_to_last_conversion() {
        let location = mauna_kea();
        let gast = GAST::from_hours(12.0);
        let last = gast.to_last(&location);

        let expected_hours = 12.0 + (-155.4783 / 15.0);
        let expected_normalized = if expected_hours < 0.0 {
            expected_hours + 24.0
        } else {
            expected_hours
        };
        assert!((last.hours() - expected_normalized).abs() < 1e-12);
    }
}
