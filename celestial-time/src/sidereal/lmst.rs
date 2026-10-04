use super::gmst::GMST;

local_sidereal_time!(LMST, GMST, to_lmst, to_gmst);

#[cfg(test)]
mod tests {
    use super::*;
    use celestial_core::location::Location;

    #[test]
    fn test_lmst_to_gmst_matches_erfa_anp() {
        // anp(LMST - elong), all in radians.
        let cases = [
            (4.894961283639734, -2.7136, 1.3253759764601476),
            (0.3, 1.2, 5.383185307179586),
            (6.0, -0.5, 0.21681469282041377),
        ];
        for (lmst, longitude, expected) in cases {
            let location = Location::new(0.5, longitude, 0.0).unwrap();
            let gmst = LMST::from_radians(lmst, &location)
                .unwrap()
                .to_gmst()
                .unwrap();
            assert_eq!(gmst.radians(), expected, "LMST {lmst}, elong {longitude}");
        }
    }
}
