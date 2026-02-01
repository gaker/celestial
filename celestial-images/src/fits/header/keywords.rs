use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameType {
    Light,
    Dark,
    Bias,
    Flat,
    Tricolor,
}

impl FrameType {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Light => "Light Frame",
            Self::Dark => "Dark Frame",
            Self::Bias => "Bias Frame",
            Self::Flat => "Flat Frame",
            Self::Tricolor => "Tricolor Image",
        }
    }
}

impl fmt::Display for FrameType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Keyword {
    pub name: String,
    pub value: Option<KeywordValue>,
    pub comment: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum KeywordValue {
    Logical(bool),
    Integer(i64),
    Real(f64),
    String(String),
    Complex(f64, f64),
}

impl Keyword {
    pub fn new(name: String) -> Self {
        Self {
            name,
            value: None,
            comment: None,
        }
    }

    pub fn with_value(mut self, value: impl Into<KeywordValue>) -> Self {
        self.value = Some(value.into());
        self
    }

    pub fn with_comment<S: Into<String>>(mut self, comment: S) -> Self {
        self.comment = Some(comment.into());
        self
    }

    pub fn logical<S: Into<String>>(name: S, value: bool) -> Self {
        Self {
            name: name.into(),
            value: Some(KeywordValue::Logical(value)),
            comment: None,
        }
    }

    pub fn integer<S: Into<String>>(name: S, value: i64) -> Self {
        Self {
            name: name.into(),
            value: Some(KeywordValue::Integer(value)),
            comment: None,
        }
    }

    pub fn real<S: Into<String>>(name: S, value: f64) -> Self {
        Self {
            name: name.into(),
            value: Some(KeywordValue::Real(value)),
            comment: None,
        }
    }

    pub fn string<S: Into<String>>(name: S, value: S) -> Self {
        Self {
            name: name.into(),
            value: Some(KeywordValue::String(value.into())),
            comment: None,
        }
    }

    /// Create a HISTORY keyword (no value, just text in comment position).
    pub fn history<S: Into<String>>(text: S) -> Self {
        Self {
            name: "HISTORY".to_string(),
            value: None,
            comment: Some(text.into()),
        }
    }

    /// Create a COMMENT keyword (no value, just text in comment position).
    pub fn comment<S: Into<String>>(text: S) -> Self {
        Self {
            name: "COMMENT".to_string(),
            value: None,
            comment: Some(text.into()),
        }
    }

    pub fn is_mandatory(&self) -> bool {
        matches!(
            self.name.as_str(),
            "SIMPLE" | "BITPIX" | "NAXIS" | "EXTEND" | "END"
        ) || (self.name.starts_with("NAXIS")
            && self.name.len() > 5
            && self.name[5..].chars().all(|c| c.is_ascii_digit()))
    }
}

impl KeywordValue {
    pub fn as_logical(&self) -> Option<bool> {
        match self {
            Self::Logical(b) => Some(*b),
            _ => None,
        }
    }

    pub fn as_integer(&self) -> Option<i64> {
        match self {
            Self::Integer(i) => Some(*i),
            _ => None,
        }
    }

    pub fn as_real(&self) -> Option<f64> {
        match self {
            Self::Real(f) => Some(*f),
            Self::Integer(i) => Some(*i as f64),
            _ => None,
        }
    }

    pub fn as_string(&self) -> Option<&str> {
        match self {
            Self::String(s) => Some(s),
            _ => None,
        }
    }
}

impl fmt::Display for KeywordValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Logical(b) => write!(f, "{}", if *b { "T" } else { "F" }),
            Self::Integer(i) => write!(f, "{}", i),
            Self::Real(r) => write!(f, "{}", r),
            Self::String(s) => write!(f, "'{}'", s),
            Self::Complex(real, imag) => write!(f, "({}, {})", real, imag),
        }
    }
}

impl From<bool> for KeywordValue {
    fn from(value: bool) -> Self {
        Self::Logical(value)
    }
}

impl From<i64> for KeywordValue {
    fn from(value: i64) -> Self {
        Self::Integer(value)
    }
}

impl From<i32> for KeywordValue {
    fn from(value: i32) -> Self {
        Self::Integer(value as i64)
    }
}

impl From<f64> for KeywordValue {
    fn from(value: f64) -> Self {
        Self::Real(value)
    }
}

impl From<f32> for KeywordValue {
    fn from(value: f32) -> Self {
        Self::Real(value as f64)
    }
}

impl From<String> for KeywordValue {
    fn from(value: String) -> Self {
        Self::String(value)
    }
}

impl From<&str> for KeywordValue {
    fn from(value: &str) -> Self {
        Self::String(value.to_string())
    }
}

impl From<(f64, f64)> for KeywordValue {
    fn from(value: (f64, f64)) -> Self {
        Self::Complex(value.0, value.1)
    }
}

#[derive(Debug, Clone)]
pub struct KeywordBuilder {
    keywords: Vec<Keyword>,
}

impl KeywordBuilder {
    pub fn new() -> Self {
        Self {
            keywords: Vec::new(),
        }
    }

    pub fn from_image(dimensions: impl AsRef<[usize]>, bitpix: crate::core::BitPix) -> Self {
        let mut builder = Self::new();
        builder.add_mandatory(dimensions.as_ref(), bitpix);
        builder
    }

    fn add_mandatory(&mut self, dimensions: &[usize], bitpix: crate::core::BitPix) {
        self.keywords.push(
            Keyword::logical("SIMPLE", true).with_comment("file does conform to FITS standard"),
        );
        self.keywords.push(
            Keyword::integer("BITPIX", bitpix.value() as i64)
                .with_comment("number of bits per data pixel"),
        );

        self.keywords.push(
            Keyword::integer("NAXIS", dimensions.len() as i64).with_comment("number of data axes"),
        );

        for (i, &dim) in dimensions.iter().enumerate() {
            self.keywords.push(
                Keyword::integer(format!("NAXIS{}", i + 1), dim as i64)
                    .with_comment(format!("length of data axis {}", i + 1)),
            );
        }

        self.keywords.push(Keyword::comment(
            "FITS (Flexible Image Transport System) format is defined in 'Astronomy",
        ));

        self.keywords.push(Keyword::comment(
            "and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H",
        ));
    }

    pub fn date<S: Into<String>>(&mut self, iso_date: S) -> &mut Self {
        let date_str: String = iso_date.into();
        self.keywords.push(Keyword::string("DATE", &date_str));
        self
    }

    pub fn date_obs<S: Into<String>>(&mut self, iso_date: S) -> &mut Self {
        let date_str: String = iso_date.into();
        self.keywords.push(Keyword::string("DATE-OBS", &date_str));
        self
    }

    pub fn date_from_utc(&mut self, utc: &celestial_time::UTC) -> &mut Self {
        self.date(utc.to_iso8601())
    }

    pub fn date_obs_from_utc(&mut self, utc: &celestial_time::UTC) -> &mut Self {
        self.date_obs(utc.to_iso8601())
    }

    pub fn object<S: Into<String>>(&mut self, name: S) -> &mut Self {
        let s: String = name.into();
        self.keywords.push(Keyword::string("OBJECT", &s));
        self
    }

    pub fn observer<S: Into<String>>(&mut self, name: S) -> &mut Self {
        let s: String = name.into();
        self.keywords.push(Keyword::string("OBSERVER", &s));
        self
    }

    pub fn telescope<S: Into<String>>(&mut self, name: S) -> &mut Self {
        let s: String = name.into();
        self.keywords.push(Keyword::string("TELESCOP", &s));
        self
    }

    pub fn instrument<S: Into<String>>(&mut self, name: S) -> &mut Self {
        let s: String = name.into();
        self.keywords.push(Keyword::string("INSTRUME", &s));
        self
    }

    pub fn exposure(&mut self, seconds: f64) -> &mut Self {
        self.keywords
            .push(Keyword::real("EXPTIME", seconds).with_comment("Exposure time in seconds"));
        self
    }

    pub fn filter<S: Into<String>>(&mut self, name: S) -> &mut Self {
        let s: String = name.into();
        self.keywords.push(Keyword::string("FILTER", &s));
        self
    }

    pub fn temperature(&mut self, celsius: f64) -> &mut Self {
        self.keywords
            .push(Keyword::real("CCD-TEMP", celsius).with_comment("CCD temperature in Celsius"));
        self
    }

    pub fn software<S: Into<String>>(&mut self, name: S) -> &mut Self {
        let s: String = name.into();
        self.keywords
            .push(Keyword::string("SWCREATE", &s).with_comment("SBIGFITSEXT Name & vers"));

        self
    }

    pub fn gain(&mut self, value: i64) -> &mut Self {
        self.keywords
            .push(Keyword::integer("GAINRAW", value).with_comment("Your gain value (integer)"));
        self
    }

    pub fn offset(&mut self, value: i64) -> &mut Self {
        self.keywords
            .push(Keyword::integer("OFFSET", value).with_comment("camera offset"));
        self
    }

    pub fn frame_type(&mut self, frame_type: FrameType) -> &mut Self {
        self.keywords.push(
            Keyword::string("IMAGETYP", frame_type.as_str())
                .with_comment("SBIGFITSEXT Light, Dark, Bias or Flat"),
        );

        let ft = match frame_type {
            FrameType::Light => 1,
            FrameType::Bias => 2,
            _ => 0,
        };

        self.keywords.push(
            Keyword::integer("PICTTYPE".to_string(), ft)
                .with_comment("Image type as index 0= Unknown 1=Light, 2=Bias,"),
        );
        self
    }

    pub fn binning(&mut self, x: u32, y: u32) -> &mut Self {
        self.keywords.push(
            Keyword::integer("XBINNING", x as i64)
                .with_comment("SBIGFITSEXT Binning factor in width"),
        );
        self.keywords.push(
            Keyword::integer("YBINNING", y as i64)
                .with_comment("SBIGFITSEXT Binning factor in height"),
        );
        self
    }

    pub fn bscale(&mut self, value: f64) -> &mut Self {
        self.keywords.push(Keyword::real("BSCALE", value));
        self
    }

    pub fn bzero(&mut self, value: f64) -> &mut Self {
        self.keywords.push(Keyword::real("BZERO", value));
        self
    }

    pub fn keyword(&mut self, kw: Keyword) -> &mut Self {
        self.keywords.push(kw);
        self
    }

    pub fn history<S: Into<String>>(&mut self, text: S) -> &mut Self {
        self.keywords.push(Keyword::history(text));
        self
    }

    pub fn comment<S: Into<String>>(&mut self, text: S) -> &mut Self {
        self.keywords.push(Keyword::comment(text));
        self
    }

    pub fn build(self) -> Vec<Keyword> {
        self.keywords
    }

    pub fn keywords(&self) -> &[Keyword] {
        &self.keywords
    }
}

impl Default for KeywordBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl From<KeywordBuilder> for Vec<Keyword> {
    fn from(builder: KeywordBuilder) -> Self {
        builder.keywords
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keyword_new_creates_empty_keyword() {
        let keyword = Keyword::new("TEST".to_string());
        assert_eq!(keyword.name, "TEST");
        assert_eq!(keyword.value, None);
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_with_value_sets_value() {
        let keyword = Keyword::new("TEST".to_string()).with_value(KeywordValue::Integer(42));

        assert_eq!(keyword.name, "TEST");
        assert_eq!(keyword.value, Some(KeywordValue::Integer(42)));
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_with_comment_sets_comment() {
        let keyword = Keyword::new("TEST".to_string()).with_comment("This is a comment");

        assert_eq!(keyword.name, "TEST");
        assert_eq!(keyword.value, None);
        assert_eq!(keyword.comment, Some("This is a comment".to_string()));
    }

    #[test]
    fn keyword_chaining() {
        let keyword = Keyword::new("TEST".to_string())
            .with_value(KeywordValue::String("hello".to_string()))
            .with_comment("A test keyword");

        assert_eq!(keyword.name, "TEST");
        assert_eq!(
            keyword.value,
            Some(KeywordValue::String("hello".to_string()))
        );
        assert_eq!(keyword.comment, Some("A test keyword".to_string()));
    }

    #[test]
    fn keyword_logical_creates_logical_keyword() {
        let keyword = Keyword::logical("SIMPLE", true);

        assert_eq!(keyword.name, "SIMPLE");
        assert_eq!(keyword.value, Some(KeywordValue::Logical(true)));
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_logical_false() {
        let keyword = Keyword::logical("EXTEND", false);

        assert_eq!(keyword.name, "EXTEND");
        assert_eq!(keyword.value, Some(KeywordValue::Logical(false)));
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_integer_creates_integer_keyword() {
        let keyword = Keyword::integer("NAXIS", 2);

        assert_eq!(keyword.name, "NAXIS");
        assert_eq!(keyword.value, Some(KeywordValue::Integer(2)));
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_integer_negative() {
        let keyword = Keyword::integer("TEST", -42);

        assert_eq!(keyword.name, "TEST");
        assert_eq!(keyword.value, Some(KeywordValue::Integer(-42)));
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_real_creates_real_keyword() {
        let keyword = Keyword::real("CRVAL1", celestial_core::constants::PI);

        assert_eq!(keyword.name, "CRVAL1");
        assert_eq!(
            keyword.value,
            Some(KeywordValue::Real(celestial_core::constants::PI))
        );
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_real_zero() {
        let keyword = Keyword::real("TEST", 0.0);

        assert_eq!(keyword.name, "TEST");
        assert_eq!(keyword.value, Some(KeywordValue::Real(0.0)));
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_string_creates_string_keyword() {
        let keyword = Keyword::string("OBJECT", "M31");

        assert_eq!(keyword.name, "OBJECT");
        assert_eq!(keyword.value, Some(KeywordValue::String("M31".to_string())));
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_string_empty() {
        let keyword = Keyword::string("TEST", "");

        assert_eq!(keyword.name, "TEST");
        assert_eq!(keyword.value, Some(KeywordValue::String("".to_string())));
        assert_eq!(keyword.comment, None);
    }

    #[test]
    fn keyword_is_mandatory() {
        let keyword = Keyword::logical("SIMPLE", true);
        assert!(keyword.is_mandatory());
    }

    #[test]
    fn keyword_is_mandatory_bitpix() {
        let keyword = Keyword::integer("BITPIX", 16);
        assert!(keyword.is_mandatory());
    }

    #[test]
    fn keyword_is_mandatory_naxis() {
        let keyword = Keyword::integer("NAXIS", 2);
        assert!(keyword.is_mandatory());
    }

    #[test]
    fn keyword_is_mandatory_naxis_numbered() {
        for i in 1..=10 {
            let keyword = Keyword::integer(format!("NAXIS{}", i), 100);
            assert!(keyword.is_mandatory());
        }
    }

    #[test]
    fn keyword_is_mandatory_extend() {
        let keyword = Keyword::logical("EXTEND", false);
        assert!(keyword.is_mandatory());
    }

    #[test]
    fn keyword_is_mandatory_end() {
        let keyword = Keyword::new("END".to_string());
        assert!(keyword.is_mandatory());
    }

    #[test]
    fn keyword_is_not_mandatory_custom() {
        let keyword = Keyword::string("OBJECT", "M31");
        assert!(!keyword.is_mandatory());
    }

    #[test]
    fn keyword_is_not_mandatory_similar_names() {
        let keyword = Keyword::string("NAXIS_TEST", "test");
        assert!(!keyword.is_mandatory());

        let keyword2 = Keyword::string("SIMPLE_TEST", "test");
        assert!(!keyword2.is_mandatory());
    }

    #[test]
    fn keyword_value_as_logical_returns_bool() {
        let value = KeywordValue::Logical(true);
        assert_eq!(value.as_logical(), Some(true));

        let value = KeywordValue::Logical(false);
        assert_eq!(value.as_logical(), Some(false));
    }

    #[test]
    fn keyword_value_as_logical_returns_none_for_other_types() {
        let value = KeywordValue::Integer(42);
        assert_eq!(value.as_logical(), None);

        let value = KeywordValue::String("true".to_string());
        assert_eq!(value.as_logical(), None);
    }

    #[test]
    fn keyword_value_as_integer_returns_integer() {
        let value = KeywordValue::Integer(42);
        assert_eq!(value.as_integer(), Some(42));

        let value = KeywordValue::Integer(-123);
        assert_eq!(value.as_integer(), Some(-123));
    }

    #[test]
    fn keyword_value_as_integer_returns_none_for_other_types() {
        let value = KeywordValue::Logical(true);
        assert_eq!(value.as_integer(), None);

        let value = KeywordValue::Real(celestial_core::constants::PI);
        assert_eq!(value.as_integer(), None);
    }

    #[test]
    fn keyword_value_as_real_returns_real() {
        let value = KeywordValue::Real(celestial_core::constants::PI);
        assert_eq!(value.as_real(), Some(celestial_core::constants::PI));

        let value = KeywordValue::Real(-2.71);
        assert_eq!(value.as_real(), Some(-2.71));
    }

    #[test]
    fn keyword_value_as_real_converts_integer() {
        let value = KeywordValue::Integer(42);
        assert_eq!(value.as_real(), Some(42.0));

        let value = KeywordValue::Integer(-5);
        assert_eq!(value.as_real(), Some(-5.0));
    }

    #[test]
    fn keyword_value_as_real_returns_none_for_other_types() {
        let value = KeywordValue::Logical(true);
        assert_eq!(value.as_real(), None);

        let value = KeywordValue::String("3.14".to_string());
        assert_eq!(value.as_real(), None);
    }

    #[test]
    fn keyword_value_as_string_returns_string() {
        let value = KeywordValue::String("hello".to_string());
        assert_eq!(value.as_string(), Some("hello"));

        let value = KeywordValue::String("".to_string());
        assert_eq!(value.as_string(), Some(""));
    }

    #[test]
    fn keyword_value_as_string_returns_none_for_other_types() {
        let value = KeywordValue::Integer(42);
        assert_eq!(value.as_string(), None);

        let value = KeywordValue::Logical(true);
        assert_eq!(value.as_string(), None);
    }

    #[test]
    fn keyword_value_display_logical() {
        let value = KeywordValue::Logical(true);
        assert_eq!(format!("{}", value), "T");

        let value = KeywordValue::Logical(false);
        assert_eq!(format!("{}", value), "F");
    }

    #[test]
    fn keyword_value_display_integer() {
        let value = KeywordValue::Integer(42);
        assert_eq!(format!("{}", value), "42");

        let value = KeywordValue::Integer(-123);
        assert_eq!(format!("{}", value), "-123");

        let value = KeywordValue::Integer(0);
        assert_eq!(format!("{}", value), "0");
    }

    #[test]
    fn keyword_value_display_real() {
        let value = KeywordValue::Real(celestial_core::constants::PI);
        assert_eq!(
            format!("{}", value),
            format!("{}", celestial_core::constants::PI)
        );

        let value = KeywordValue::Real(-2.71);
        assert_eq!(format!("{}", value), "-2.71");

        let value = KeywordValue::Real(0.0);
        assert_eq!(format!("{}", value), "0");
    }

    #[test]
    fn keyword_value_display_string() {
        let value = KeywordValue::String("hello".to_string());
        assert_eq!(format!("{}", value), "'hello'");

        let value = KeywordValue::String("".to_string());
        assert_eq!(format!("{}", value), "''");

        let value = KeywordValue::String("test with spaces".to_string());
        assert_eq!(format!("{}", value), "'test with spaces'");
    }

    #[test]
    fn keyword_value_display_complex() {
        let value = KeywordValue::Complex(1.0, 2.0);
        assert_eq!(format!("{}", value), "(1, 2)");

        let value = KeywordValue::Complex(-1.5, celestial_core::constants::PI);
        assert_eq!(
            format!("{}", value),
            format!("(-1.5, {})", celestial_core::constants::PI)
        );

        let value = KeywordValue::Complex(0.0, 0.0);
        assert_eq!(format!("{}", value), "(0, 0)");
    }

    #[test]
    fn keyword_value_equality() {
        let val1 = KeywordValue::Integer(42);
        let val2 = KeywordValue::Integer(42);
        let val3 = KeywordValue::Integer(43);

        assert_eq!(val1, val2);
        assert_ne!(val1, val3);

        let val4 = KeywordValue::String("test".to_string());
        let val5 = KeywordValue::String("test".to_string());
        let val6 = KeywordValue::String("other".to_string());

        assert_eq!(val4, val5);
        assert_ne!(val4, val6);
        assert_ne!(val1, val4);
    }

    #[test]
    fn keyword_equality() {
        let kw1 = Keyword::integer("NAXIS", 2);
        let kw2 = Keyword::integer("NAXIS", 2);
        let kw3 = Keyword::integer("NAXIS", 3);
        let kw4 = Keyword::integer("BITPIX", 2);

        assert_eq!(kw1, kw2);
        assert_ne!(kw1, kw3);
        assert_ne!(kw1, kw4);
    }

    #[test]
    fn keyword_clone() {
        let original = Keyword::string("OBJECT", "M31").with_comment("Andromeda Galaxy");
        let cloned = original.clone();

        assert_eq!(original, cloned);
        assert_eq!(cloned.name, "OBJECT");
        assert_eq!(cloned.value, Some(KeywordValue::String("M31".to_string())));
        assert_eq!(cloned.comment, Some("Andromeda Galaxy".to_string()));
    }

    #[test]
    fn keyword_value_clone() {
        let original = KeywordValue::Complex(1.0, 2.0);
        let cloned = original.clone();

        assert_eq!(original, cloned);
        if let KeywordValue::Complex(r, i) = cloned {
            assert_eq!(r, 1.0);
            assert_eq!(i, 2.0);
        } else {
            panic!("Expected Complex value");
        }
    }

    #[test]
    fn keyword_comprehensive_example() {
        let keyword = Keyword::new("CRVAL1".to_string())
            .with_value(KeywordValue::Real(180.0))
            .with_comment("Reference coordinate value");

        assert_eq!(keyword.name, "CRVAL1");
        assert_eq!(keyword.value.as_ref().unwrap().as_real().unwrap(), 180.0);
        assert_eq!(
            keyword.comment.as_ref().unwrap(),
            "Reference coordinate value"
        );
        assert!(!keyword.is_mandatory());
        assert_eq!(format!("{}", keyword.value.as_ref().unwrap()), "180");
    }

    #[test]
    fn keyword_builder_from_image() {
        use crate::core::BitPix;
        let mut builder = KeywordBuilder::from_image(&[1024, 768], BitPix::I16);
        builder.object("M31");
        builder.exposure(300.0);
        builder.filter("Ha");
        let keywords = builder.build();

        assert!(keywords.iter().any(|k| k.name == "SIMPLE"));
        assert!(keywords
            .iter()
            .any(|k| k.name == "BITPIX" && k.value == Some(KeywordValue::Integer(16))));
        assert!(keywords
            .iter()
            .any(|k| k.name == "NAXIS" && k.value == Some(KeywordValue::Integer(2))));
        assert!(keywords
            .iter()
            .any(|k| k.name == "NAXIS1" && k.value == Some(KeywordValue::Integer(1024))));
        assert!(keywords
            .iter()
            .any(|k| k.name == "NAXIS2" && k.value == Some(KeywordValue::Integer(768))));
        assert!(keywords.iter().any(|k| k.name == "OBJECT"));
        assert!(keywords.iter().any(|k| k.name == "EXPTIME"));
        assert!(keywords.iter().any(|k| k.name == "FILTER"));
    }

    #[test]
    fn keyword_builder_into_vec() {
        use crate::core::BitPix;
        let builder = KeywordBuilder::from_image(&[100, 100], BitPix::F32);
        let keywords: Vec<Keyword> = builder.into();
        assert!(!keywords.is_empty());
    }

    #[test]
    fn keyword_builder_date_from_utc() {
        use crate::core::BitPix;
        let utc = celestial_time::UTC::j2000();
        let mut builder = KeywordBuilder::from_image(&[100, 100], BitPix::F32);
        builder.date_from_utc(&utc);
        let keywords = builder.build();

        let date_kw = keywords.iter().find(|k| k.name == "DATE").unwrap();
        let date_val = date_kw.value.as_ref().unwrap().as_string().unwrap();
        assert!(date_val.starts_with("2000-01-01T12:"));
    }
}
