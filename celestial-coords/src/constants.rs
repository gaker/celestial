// Rotation from ICRS to Galactic coordinates: the rows are the Galactic x, y and z axes in ICRS.
// The Galactic frame is the one the Hipparcos Catalogue (ESA 1997, vol. 1 §1.5.3) ties to ICRS,
// with the north Galactic pole at RA 192.85948°, Dec 27.12825°.
pub(crate) const ICRS_TO_GALACTIC: [[f64; 3]; 3] = [
    [
        -0.05487556041621537,
        -0.873437090234885,
        -0.4838350155487132,
    ],
    [0.49410942787558365, -0.4448296299600112, 0.7469822444972188],
    [
        -0.8676661490190047,
        -0.19807637343120152,
        0.4559837761750669,
    ],
];

// Light time for one astronomical unit, in days: multiplying a velocity in
// AU/day by it gives the velocity in units of c.
pub(crate) const AU_LIGHT_TIME_DAYS: f64 = (celestial_core::constants::AU_M
    / celestial_core::constants::SPEED_OF_LIGHT_M_PER_S)
    / celestial_core::constants::SECONDS_PER_DAY_F64;
