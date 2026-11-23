use super::Image;
use crate::fits::header::keywords::Keyword;
use crate::fits::header::metadata;

impl Image {
    pub fn exposure(&mut self, seconds: f64) -> &mut Self {
        self.with(metadata::exposure(seconds))
    }

    pub fn temperature(&mut self, celsius: f64) -> &mut Self {
        self.with(metadata::temperature(celsius))
    }

    pub fn gain(&mut self, gain: f64) -> &mut Self {
        self.with(metadata::gain(gain))
    }

    pub fn binning(&mut self, x: u32, y: u32) -> &mut Self {
        let [x, y] = metadata::binning(x, y);
        self.with(x).with(y)
    }

    pub fn filter(&mut self, name: &str) -> &mut Self {
        self.with(metadata::filter(name))
    }

    pub fn object(&mut self, name: &str) -> &mut Self {
        self.with(metadata::object(name))
    }

    pub fn observer(&mut self, name: &str) -> &mut Self {
        self.with(metadata::observer(name))
    }

    pub fn telescope(&mut self, name: &str) -> &mut Self {
        self.with(metadata::telescope(name))
    }

    pub fn instrument(&mut self, name: &str) -> &mut Self {
        self.with(metadata::instrument(name))
    }

    pub fn date_obs(&mut self, iso8601: &str) -> &mut Self {
        self.with(metadata::date_obs(iso8601))
    }

    fn with(&mut self, kw: Keyword) -> &mut Self {
        self.set_keyword(kw);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::keywords::KeywordValue;
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
    fn get_last_returns_the_last_of_a_repeated_name() {
        let mut i = img();
        i.keywords.push(Keyword::real("CRVAL1", 83.6));
        i.keywords.push(Keyword::real("HIERARCH ESO DET", 1.0));
        i.keywords.push(Keyword::real("CRVAL1", 84.0));
        i.keywords.push(Keyword::real("ESO DET", 2.0));
        assert_eq!(i.get("CRVAL1").as_f64(), Some(83.6));
        assert_eq!(i.get_last("CRVAL1").as_f64(), Some(84.0));
        assert_eq!(i.get_last("ESO DET").as_f64(), Some(2.0));
        assert_eq!(i.get_last("HIERARCH ESO DET").as_f64(), Some(1.0));
        assert!(!i.get_last("CRVAL2").exists());
    }

    #[test]
    fn set_keyword_appends_commentary_cards() {
        let mut i = img();
        i.set_keyword(Keyword::history("a"));
        i.set_keyword(Keyword::history("b"));
        let history: Vec<_> = i
            .keywords
            .iter()
            .filter(|k| k.name == "HISTORY")
            .map(|k| k.comment.as_deref())
            .collect();
        assert_eq!(history, [Some("a"), Some("b")]);
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
