use super::Image;
use crate::fits::header::Keyword;

impl Image {
    pub fn exposure(&mut self, seconds: f64) -> &mut Self {
        self.set_keyword(Keyword::real("EXPTIME", seconds));
        self
    }

    pub fn temperature(&mut self, celsius: f64) -> &mut Self {
        self.set_keyword(Keyword::real("CCD-TEMP", celsius));
        self
    }

    pub fn gain(&mut self, gain: f64) -> &mut Self {
        self.set_keyword(Keyword::real("GAIN", gain));
        self
    }

    pub fn binning(&mut self, x: u32, y: u32) -> &mut Self {
        self.set_keyword(Keyword::integer("XBINNING", x as i64));
        self.set_keyword(Keyword::integer("YBINNING", y as i64));
        self
    }

    pub fn filter(&mut self, name: &str) -> &mut Self {
        self.set_keyword(Keyword::string("FILTER", name));
        self
    }

    pub fn object(&mut self, name: &str) -> &mut Self {
        self.set_keyword(Keyword::string("OBJECT", name));
        self
    }

    pub fn observer(&mut self, name: &str) -> &mut Self {
        self.set_keyword(Keyword::string("OBSERVER", name));
        self
    }

    pub fn telescope(&mut self, name: &str) -> &mut Self {
        self.set_keyword(Keyword::string("TELESCOP", name));
        self
    }

    pub fn instrument(&mut self, name: &str) -> &mut Self {
        self.set_keyword(Keyword::string("INSTRUME", name));
        self
    }

    pub fn date_obs(&mut self, iso8601: &str) -> &mut Self {
        self.set_keyword(Keyword::string("DATE-OBS", iso8601));
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::KeywordValue;
    use crate::formats::pixel_data::PixelData;

    fn img() -> Image {
        Image::new(PixelData::U8(vec![0; 4]), vec![2usize, 2])
    }

    #[test]
    fn exposure_sets_exptime() {
        let mut i = img();
        i.exposure(30.0);
        assert_eq!(
            i.get_keyword("EXPTIME").unwrap().value,
            Some(KeywordValue::Real(30.0))
        );
    }

    #[test]
    fn temperature_sets_ccd_temp() {
        let mut i = img();
        i.temperature(-10.5);
        assert_eq!(
            i.get_keyword("CCD-TEMP").unwrap().value,
            Some(KeywordValue::Real(-10.5))
        );
    }

    #[test]
    fn gain_sets_gain() {
        let mut i = img();
        i.gain(1.5);
        assert_eq!(
            i.get_keyword("GAIN").unwrap().value,
            Some(KeywordValue::Real(1.5))
        );
    }

    #[test]
    fn binning_sets_both_axes() {
        let mut i = img();
        i.binning(2, 3);
        assert_eq!(
            i.get_keyword("XBINNING").unwrap().value,
            Some(KeywordValue::Integer(2))
        );
        assert_eq!(
            i.get_keyword("YBINNING").unwrap().value,
            Some(KeywordValue::Integer(3))
        );
    }

    #[test]
    fn string_setters_write_correct_keywords() {
        let mut i = img();
        i.filter("V")
            .object("M31")
            .observer("Hubble")
            .telescope("Keck")
            .instrument("WFC3")
            .date_obs("2026-04-21T20:00:00");

        assert_eq!(
            i.get_keyword("FILTER").unwrap().value,
            Some(KeywordValue::String("V".to_string()))
        );
        assert_eq!(
            i.get_keyword("OBJECT").unwrap().value,
            Some(KeywordValue::String("M31".to_string()))
        );
        assert_eq!(
            i.get_keyword("OBSERVER").unwrap().value,
            Some(KeywordValue::String("Hubble".to_string()))
        );
        assert_eq!(
            i.get_keyword("TELESCOP").unwrap().value,
            Some(KeywordValue::String("Keck".to_string()))
        );
        assert_eq!(
            i.get_keyword("INSTRUME").unwrap().value,
            Some(KeywordValue::String("WFC3".to_string()))
        );
        assert_eq!(
            i.get_keyword("DATE-OBS").unwrap().value,
            Some(KeywordValue::String("2026-04-21T20:00:00".to_string()))
        );
    }

    #[test]
    fn setters_return_self_for_chaining() {
        let mut i = img();
        let ret = i.exposure(10.0).filter("R").gain(2.0).binning(1, 1);
        // Chaining compiles and returns an &mut Image — no need to assert the address.
        let _: &mut Image = ret;
    }

    #[test]
    fn setters_overwrite_previous_values() {
        let mut i = img();
        i.exposure(10.0);
        i.exposure(30.0);
        assert_eq!(
            i.get_keyword("EXPTIME").unwrap().value,
            Some(KeywordValue::Real(30.0))
        );
        assert_eq!(i.keywords.iter().filter(|k| k.name == "EXPTIME").count(), 1);
    }
}
