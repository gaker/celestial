use super::Angle;
use crate::constants::{
    ARCMIN_TO_RAD, ARCMIN_TO_RAD_LO, ARCSEC_PER_RAD, ARCSEC_PER_RAD_LO, ARCSEC_TO_RAD,
    ARCSEC_TO_RAD_LO, DEG_TO_RAD, DEG_TO_RAD_LO, HOURS_TO_RAD, HOURS_TO_RAD_LO, RAD_TO_ARCMIN,
    RAD_TO_ARCMIN_LO, RAD_TO_DEG, RAD_TO_DEG_LO, RAD_TO_HOURS, RAD_TO_HOURS_LO,
};

impl Angle {
    /// Creates an angle from degrees.
    ///
    /// # Example
    ///
    /// ```
    /// use celestial_core::angle::Angle;
    ///
    /// let angle = Angle::from_degrees(180.0);
    /// assert_eq!(angle.radians(), celestial_core::constants::PI);
    /// ```
    #[inline]
    pub fn from_degrees(deg: f64) -> Self {
        Self::from_radians(mul_split(deg, DEG_TO_RAD, DEG_TO_RAD_LO))
    }

    /// Creates an angle from hours.
    ///
    /// In astronomy, right ascension is measured in hours where 24h = 360 degrees.
    /// Each hour equals 15 degrees.
    ///
    /// # Example
    ///
    /// ```
    /// use celestial_core::angle::Angle;
    ///
    /// let ra = Angle::from_hours(6.0);  // 6h = 90 degrees
    /// assert_eq!(ra.degrees(), 90.0);
    ///
    /// let ra_24h = Angle::from_hours(24.0);  // Full circle
    /// assert_eq!(ra_24h.degrees(), 360.0);
    /// ```
    #[inline]
    pub fn from_hours(h: f64) -> Self {
        Self::from_radians(mul_split(h, HOURS_TO_RAD, HOURS_TO_RAD_LO))
    }

    /// Creates an angle from arcseconds.
    ///
    /// One arcsecond = 1/3600 of a degree. Commonly used for:
    /// - Parallax measurements
    /// - Proper motion
    /// - Small angular separations
    ///
    /// # Example
    ///
    /// ```
    /// use celestial_core::angle::Angle;
    ///
    /// let angle = Angle::from_arcseconds(3600.0);  // 1 degree
    /// assert_eq!(angle.degrees(), 1.0);
    ///
    /// // Proxima Centauri's parallax is about 0.77 arcseconds
    /// let parallax = Angle::from_arcseconds(0.77);
    /// ```
    #[inline]
    pub fn from_arcseconds(arcsec: f64) -> Self {
        Self::from_radians(mul_split(arcsec, ARCSEC_TO_RAD, ARCSEC_TO_RAD_LO))
    }

    /// Creates an angle from arcminutes.
    ///
    /// One arcminute = 1/60 of a degree. Commonly used for:
    /// - Field of view specifications
    /// - Object sizes (e.g., the Moon is about 31 arcminutes)
    ///
    /// # Example
    ///
    /// ```
    /// use celestial_core::angle::Angle;
    ///
    /// let angle = Angle::from_arcminutes(60.0);  // 1 degree
    /// assert_eq!(angle.degrees(), 1.0);
    ///
    /// // Full Moon's apparent diameter
    /// let moon_diameter = Angle::from_arcminutes(31.0);
    /// ```
    #[inline]
    pub fn from_arcminutes(arcmin: f64) -> Self {
        Self::from_radians(mul_split(arcmin, ARCMIN_TO_RAD, ARCMIN_TO_RAD_LO))
    }

    /// Returns the angle in degrees.
    #[inline]
    pub fn degrees(self) -> f64 {
        mul_split(self.radians(), RAD_TO_DEG, RAD_TO_DEG_LO)
    }

    /// Returns the angle in hours.
    ///
    /// Useful for right ascension where 24h = 360 degrees.
    #[inline]
    pub fn hours(self) -> f64 {
        mul_split(self.radians(), RAD_TO_HOURS, RAD_TO_HOURS_LO)
    }

    /// Returns the angle in arcseconds.
    #[inline]
    pub fn arcseconds(self) -> f64 {
        mul_split(self.radians(), ARCSEC_PER_RAD, ARCSEC_PER_RAD_LO)
    }

    /// Returns the angle in arcminutes.
    #[inline]
    pub fn arcminutes(self) -> f64 {
        mul_split(self.radians(), RAD_TO_ARCMIN, RAD_TO_ARCMIN_LO)
    }
}

// head + tail holds the constant to ~106 bits and the fma rounds once, which gives the
// correctly rounded product for every x whose x * tail is not subnormal.
fn mul_split(x: f64, head: f64, tail: f64) -> f64 {
    libm::fma(x, head, x * tail)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_arcseconds() {
        let angle = Angle::from_arcseconds(3600.0);
        assert_eq!(angle.degrees(), 1.0);
    }

    #[test]
    fn test_from_arcminutes() {
        let angle = Angle::from_arcminutes(60.0);
        assert_eq!(angle.degrees(), 1.0);
    }

    #[test]
    fn test_arcseconds_getter() {
        let angle = Angle::from_degrees(1.0);
        assert_eq!(angle.arcseconds(), 3600.0);
    }

    #[test]
    fn test_arcminutes_getter() {
        let angle = Angle::from_degrees(1.0);
        assert_eq!(angle.arcminutes(), 60.0);
    }

    #[test]
    fn test_arc_unit_getters_invert_constructors() {
        assert_eq!(Angle::from_arcseconds(0.371).arcseconds(), 0.371);
        assert_eq!(Angle::from_arcminutes(0.371).arcminutes(), 0.371);
    }

    // Expected values are the exact results rounded once to f64.
    #[test]
    fn test_hour_conversions_are_correctly_rounded() {
        assert_eq!(Angle::from_hours(2.5).radians(), 0.6544984694978736);
        assert_eq!(Angle::from_radians(1.0).hours(), 3.819718634205488);
    }

    // (input, exact product rounded once) for inputs whose exact product lies closest to a
    // rounding midpoint, found by lattice reduction over every f64 mantissa.
    const RAD_TO_HOURS_CASES: [(f64, f64); 6] = [
        (1.8004476283930817, 6.877203356084132),
        (1.409763821700704, 5.384901139578922),
        (1.89553258904251, 7.240401152109449),
        (1.885188624947846, 7.200890119705508),
        (-2.13002350284922e-271, -8.136090465128812e-271),
        (1.1916328275158352e271, 4.55170211639321e271),
    ];
    const RAD_TO_ARCMIN_CASES: [(f64, f64); 6] = [
        (1.8084124859401745, 6216.864183787999),
        (1.2953833795781495, 4453.200030073264),
        (1.6416729883771952, 5643.6560144785635),
        (1.987962597176241, 6834.111998883863),
        (-2.139446344982815e-271, -7.354874763732313e-268),
        (1.0949503282482752e271, 3.764161955105417e274),
    ];
    // In the HOURS_TO_RAD, RAD_TO_DEG and DEG_TO_RAD cases the first input is the one a
    // round-to-nearest tail gets wrong.
    const HOURS_TO_RAD_CASES: [(f64, f64); 7] = [
        (1.06888805882525, 0.27983423942627167),
        (1.5327630619578818, 0.40127643126172324),
        (1.9561772196378415, 0.5121259985278291),
        (1.6235625109642111, 0.42504767142408034),
        (1.936421703435871, 0.5069540164804972),
        (-1.8133386913249516e-271, -4.747309592613832e-272),
        (1.6535003633069484e271, 4.328853828394302e270),
    ];
    const ARCMIN_TO_RAD_CASES: [(f64, f64); 6] = [
        (1.372469773336267, 0.0003992352738136357),
        (1.663108441403967, 0.0004837786353368402),
        (1.8837754091165113, 0.0005479680542864389),
        (1.2325263726259557, 0.00035852738866642395),
        (-1.6237033657932884e-271, -4.723161634801126e-275),
        (1.4057777508468412e271, 4.0892417172596474e267),
    ];

    const RAD_TO_DEG_CASES: [(f64, f64); 6] = [
        (1.8557083013070463, 106.32425367228505),
        (1.0897770731007919, 62.43962689879516),
        (1.2208272251910588, 69.948247518115),
        (1.3827501476383275, 79.22574758076763),
        (-2.195399762749133e-271, -1.2578714074954756e-269),
        (9.211572286018887e270, 5.277842146685578e272),
    ];
    const DEG_TO_RAD_CASES: [(f64, f64); 7] = [
        (1.5888359943899804, 0.02773041937630331),
        (1.7337969328689278, 0.030260465039541887),
        (1.8844922002263165, 0.032890593622101456),
        (1.0204608982310732, 0.01781040256199101),
        (1.1654218367100206, 0.020340448225229582),
        (-2.0511722517998337e-271, -3.5799709319453287e-273),
        (1.592907077355808e271, 2.7801473178178865e269),
    ];
    const RAD_TO_ARCSEC_CASES: [(f64, f64); 6] = [
        (1.11411846919191, 229803.43018418088),
        (1.7348776422955747, 357844.200750516),
        (1.7511178542690082, 361193.98492662807),
        (1.1303586811653434, 233153.2143602929),
        (-1.3180602906262695e-271, -2.7186945046801897e-266),
        (1.4664421929828624e271, 3.024754148081772e276),
    ];
    const ARCSEC_TO_RAD_CASES: [(f64, f64); 6] = [
        (1.559164163816219, 7.559041177138128e-6),
        (1.1709704244292953, 5.6770248193796245e-6),
        (1.9665255220087154, 9.533984773208973e-6),
        (1.5719425763532675, 7.620992669346355e-6),
        (-1.8445725725956692e-271, -8.942740189937932e-277),
        (9.897876341561688e270, 4.7986258643195094e265),
    ];

    fn assert_converts(cases: &[(f64, f64)], convert: fn(f64) -> f64) {
        for &(input, expected) in cases {
            assert_eq!(convert(input), expected, "{input:e}");
        }
    }

    #[test]
    fn test_unit_getters_are_correctly_rounded_near_midpoints() {
        assert_converts(&RAD_TO_DEG_CASES, |r| Angle::from_radians(r).degrees());
        assert_converts(&RAD_TO_HOURS_CASES, |r| Angle::from_radians(r).hours());
        assert_converts(&RAD_TO_ARCMIN_CASES, |r| {
            Angle::from_radians(r).arcminutes()
        });
        assert_converts(&RAD_TO_ARCSEC_CASES, |r| {
            Angle::from_radians(r).arcseconds()
        });
    }

    #[test]
    fn test_unit_constructors_are_correctly_rounded_near_midpoints() {
        assert_converts(&DEG_TO_RAD_CASES, |d| Angle::from_degrees(d).radians());
        assert_converts(&HOURS_TO_RAD_CASES, |h| Angle::from_hours(h).radians());
        assert_converts(&ARCMIN_TO_RAD_CASES, |m| {
            Angle::from_arcminutes(m).radians()
        });
        assert_converts(&ARCSEC_TO_RAD_CASES, |s| {
            Angle::from_arcseconds(s).radians()
        });
    }
}
