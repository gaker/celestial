use super::keywords::Keyword;

// The cards the metadata setters of `Image` and `AstroImage` write, so the
// two write the same name and value type.

pub(crate) fn exposure(seconds: f64) -> Keyword {
    Keyword::real("EXPTIME", seconds).with_comment("Exposure time in seconds")
}

pub(crate) fn temperature(celsius: f64) -> Keyword {
    Keyword::real("CCD-TEMP", celsius).with_comment("CCD temperature in Celsius")
}

pub(crate) fn gain(gain: f64) -> Keyword {
    Keyword::real("GAIN", gain)
}

pub(crate) fn binning(x: u32, y: u32) -> [Keyword; 2] {
    [
        Keyword::integer("XBINNING", i64::from(x)),
        Keyword::integer("YBINNING", i64::from(y)),
    ]
}

pub(crate) fn filter(name: &str) -> Keyword {
    Keyword::string("FILTER", name)
}

pub(crate) fn object(name: &str) -> Keyword {
    Keyword::string("OBJECT", name)
}

pub(crate) fn observer(name: &str) -> Keyword {
    Keyword::string("OBSERVER", name)
}

pub(crate) fn telescope(name: &str) -> Keyword {
    Keyword::string("TELESCOP", name)
}

pub(crate) fn instrument(name: &str) -> Keyword {
    Keyword::string("INSTRUME", name)
}

pub(crate) fn date_obs(iso8601: &str) -> Keyword {
    Keyword::string("DATE-OBS", iso8601)
}
