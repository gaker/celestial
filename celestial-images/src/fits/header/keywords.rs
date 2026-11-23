use crate::fits::io::card;
use std::convert::Infallible;
use std::fmt;

// FITS 4.0 §4.1.2.1: a keyword name has at most 8 characters, each one of
// A-Z, 0-9, '-' and '_'.
pub(crate) fn is_standard_name(name: &str) -> bool {
    name.len() <= 8
        && name
            .bytes()
            .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'-' || b == b'_')
}

// The other spelling of a HIERARCH keyword, which a lookup also tries.
// "ESO DET" and "HIERARCH ESO DET" name the same card, since the writer
// stores a non-standard name under the HIERARCH convention. A standard name
// has no alias, so "NAXIS1" never matches "HIERARCH NAXIS1".
pub(crate) fn hierarch_alias(name: &str) -> Option<String> {
    match name.strip_prefix("HIERARCH ") {
        Some(rest) => (!is_standard_name(rest)).then(|| rest.to_string()),
        None => (!is_standard_name(name)).then(|| format!("HIERARCH {}", name)),
    }
}

// The first keyword named `name`, or else the first named with its HIERARCH
// alias. astropy and cfitsio return the first card of a repeated name.
pub(crate) fn find_keyword(keywords: &[Keyword], name: &str) -> Option<usize> {
    let find = |name: &str| keywords.iter().position(|k| k.name == name);
    find(name).or_else(|| hierarch_alias(name).and_then(|full| find(&full)))
}

// The last keyword named `name`, or else the last named with its HIERARCH
// alias. Readers apply the last of a repeated scaling or WCS card.
pub(crate) fn find_last_keyword(keywords: &[Keyword], name: &str) -> Option<usize> {
    let find = |name: &str| keywords.iter().rposition(|k| k.name == name);
    find(name).or_else(|| hierarch_alias(name).and_then(|full| find(&full)))
}

// Replaces the first card with the keyword's name or its HIERARCH alias and
// removes the later ones, or appends the keyword. Data and WCS readers use
// the last of a repeated card, so a later copy would hide the new value.
// The replacement keeps the stored spelling of the name. A header may hold
// many commentary cards, so those are always appended.
pub(crate) fn set_keyword(keywords: &mut Vec<Keyword>, keyword: Keyword) {
    if keyword.is_commentary() {
        keywords.push(keyword);
        return;
    }
    let alias = hierarch_alias(&keyword.name);
    let named = |k: &Keyword| k.name == keyword.name || alias.as_deref() == Some(&k.name);
    let Some(first) = keywords.iter().position(named) else {
        keywords.push(keyword);
        return;
    };
    let mut index = 0;
    keywords.retain(|k| {
        index += 1;
        index - 1 == first || !named(k)
    });
    if let Some(existing) = keywords.get_mut(first) {
        let name = std::mem::take(&mut existing.name);
        *existing = Keyword { name, ..keyword };
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Keyword {
    pub name: String,
    pub value: Option<KeywordValue>,
    pub comment: Option<String>,
    pub(super) commentary: bool,
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
            commentary: false,
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
            commentary: false,
        }
    }

    pub fn integer<S: Into<String>>(name: S, value: i64) -> Self {
        Self {
            name: name.into(),
            value: Some(KeywordValue::Integer(value)),
            comment: None,
            commentary: false,
        }
    }

    pub fn real<S: Into<String>>(name: S, value: f64) -> Self {
        Self {
            name: name.into(),
            value: Some(KeywordValue::Real(value)),
            comment: None,
            commentary: false,
        }
    }

    pub fn string<S: Into<String>>(name: S, value: S) -> Self {
        Self {
            name: name.into(),
            value: Some(KeywordValue::String(value.into())),
            comment: None,
            commentary: false,
        }
    }

    /// Create a HISTORY keyword (no value, just text in comment position).
    pub fn history<S: Into<String>>(text: S) -> Self {
        Self {
            name: "HISTORY".to_string(),
            value: None,
            comment: Some(text.into()),
            commentary: true,
        }
    }

    /// Create a COMMENT keyword (no value, just text in comment position).
    pub fn comment<S: Into<String>>(text: S) -> Self {
        Self {
            name: "COMMENT".to_string(),
            value: None,
            comment: Some(text.into()),
            commentary: true,
        }
    }

    // A card without the value indicator "= " in columns 9-10 holds text in
    // columns 9-80 (FITS 4.0 §4.1.2.2). The writer writes the name and the
    // text back unchanged.
    pub fn commentary<S: Into<String>>(name: S, text: S) -> Self {
        Self {
            name: name.into(),
            value: None,
            comment: Some(text.into()),
            commentary: true,
        }
    }

    pub fn is_commentary(&self) -> bool {
        self.commentary || matches!(self.name.as_str(), "COMMENT" | "HISTORY" | "")
    }

    // XTENSION, PCOUNT and GCOUNT are the extension keywords of FITS 4.0
    // §4.4.1.2. Writers generate them, so an image read from an extension
    // does not carry them into a primary header.
    pub fn is_mandatory(&self) -> bool {
        matches!(
            self.name.as_str(),
            "SIMPLE" | "BITPIX" | "NAXIS" | "EXTEND" | "END" | "XTENSION" | "PCOUNT" | "GCOUNT"
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

// The value field as a FITS card holds it. A NaN or infinite real has no
// FITS form, so it prints as Rust prints it.
impl fmt::Display for KeywordValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Ok(text) = card::value_field(self, |r| Ok::<_, Infallible>(r.to_string()));
        f.write_str(&text)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_keyword_replaces_the_first_card_and_removes_later_copies() {
        let mut keywords = vec![
            Keyword::real("CRVAL1", 83.6),
            Keyword::string("ESO DET", "a"),
            Keyword::real("CRVAL1", 84.0),
            Keyword::string("HIERARCH ESO DET", "b"),
            Keyword::history("h"),
        ];
        set_keyword(&mut keywords, Keyword::real("CRVAL1", 50.0));
        set_keyword(&mut keywords, Keyword::string("HIERARCH ESO DET", "c"));
        set_keyword(&mut keywords, Keyword::history("h"));
        set_keyword(&mut keywords, Keyword::new(String::new()));
        set_keyword(&mut keywords, Keyword::new(String::new()));
        assert_eq!(
            keywords,
            vec![
                Keyword::real("CRVAL1", 50.0),
                Keyword::string("ESO DET", "c"),
                Keyword::history("h"),
                Keyword::history("h"),
                Keyword::new(String::new()),
                Keyword::new(String::new()),
            ]
        );
    }

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
    fn keyword_is_mandatory_extension_keywords() {
        assert!(Keyword::string("XTENSION", "IMAGE").is_mandatory());
        assert!(Keyword::integer("PCOUNT", 0).is_mandatory());
        assert!(Keyword::integer("GCOUNT", 1).is_mandatory());
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
        assert_eq!(format!("{}", value), "0.0");

        let value = KeywordValue::Real(512.0);
        assert_eq!(format!("{}", value), "512.0");

        let value = KeywordValue::Real(1e25);
        assert_eq!(format!("{}", value), "1.0E+25");

        let value = KeywordValue::Real(4.57e-19);
        assert_eq!(format!("{}", value), "4.57E-19");
    }

    // FITS cannot hold these, so they print as Rust prints them.
    #[test]
    fn keyword_value_display_non_finite_real() {
        assert_eq!(KeywordValue::Real(f64::NAN).to_string(), "NaN");
        assert_eq!(
            KeywordValue::Complex(1.0, f64::INFINITY).to_string(),
            "(1.0, inf)"
        );
    }

    #[test]
    fn keyword_value_display_string() {
        let value = KeywordValue::String("hello".to_string());
        assert_eq!(format!("{}", value), "'hello'");

        let value = KeywordValue::String("".to_string());
        assert_eq!(format!("{}", value), "''");

        let value = KeywordValue::String("test with spaces".to_string());
        assert_eq!(format!("{}", value), "'test with spaces'");

        let value = KeywordValue::String("O'Brien  ".to_string());
        assert_eq!(format!("{}", value), "'O''Brien'");
    }

    #[test]
    fn keyword_value_display_complex() {
        let value = KeywordValue::Complex(1.0, 2.0);
        assert_eq!(format!("{}", value), "(1.0, 2.0)");

        let value = KeywordValue::Complex(-1.5, celestial_core::constants::PI);
        assert_eq!(
            format!("{}", value),
            format!("(-1.5, {})", celestial_core::constants::PI)
        );

        let value = KeywordValue::Complex(0.0, 0.0);
        assert_eq!(format!("{}", value), "(0.0, 0.0)");
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
    fn keyword_example() {
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
        assert_eq!(format!("{}", keyword.value.as_ref().unwrap()), "180.0");
    }
}
