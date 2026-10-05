use crate::planetary_coefficients::pluto::A;

// VSOP2013p9.dat line 76: the mu multiplier 72490 is outside i16.
#[test]
fn table_keeps_multipliers_beyond_i16() {
    let term = A[0]
        .terms
        .iter()
        .find(|t| t.c == 0.1200203132340981e-3)
        .unwrap();
    assert_eq!(term.s, -0.3107290842332163e-4);
    assert_eq!(term.mult, [72490, 0, 0, 0, 0, 0]);
    assert_eq!(term.index, [13, 0, 0, 0, 0, 0]);
}
