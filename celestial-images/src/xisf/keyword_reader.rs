use crate::fits::header::keywords::{is_standard_name, Keyword};
use crate::fits::header::parser::HeaderCard;
use crate::xisf::errors::{Result, XisfError};
use crate::xisf::reader::attribute_text;
use quick_xml::events::BytesStart;
use std::collections::{HashMap, HashSet};

// The FITSKeyword and Reference elements of a header, collected while it is
// parsed. XISF 1.0 section 11.6 gives an image its FITSKeyword children and
// the root FITSKeyword elements its Reference children name by uid.
#[derive(Default)]
pub(super) struct KeywordScan {
    // The FITSKeyword children of the open Image element.
    pub(super) image: Vec<Keyword>,
    references: Vec<Reference>,
    // The uids of `references`, so a repeated Reference adds nothing.
    image_uids: HashSet<String>,
    images: Vec<ImageKeywords>,
    root: Vec<(Option<String>, Keyword)>,
    // The ref attribute of every Reference element in the header.
    referenced: HashSet<String>,
}

// The ref attribute of a Reference child of an Image element, and the
// number of FITSKeyword children before it.
struct Reference {
    position: usize,
    uid: String,
}

struct ImageKeywords {
    own: Vec<Keyword>,
    references: Vec<Reference>,
}

// The root keywords, split by how they are associated with images.
struct RootKeywords<'a> {
    // The keywords with a uid that a Reference element names.
    by_uid: HashMap<&'a str, Vec<&'a Keyword>>,
    // The keywords that apply to every image.
    shared: Vec<&'a Keyword>,
}

impl KeywordScan {
    pub(super) fn open_image(&mut self) {
        self.image.clear();
        self.references.clear();
        self.image_uids.clear();
    }

    pub(super) fn finish_image(&mut self) {
        self.images.push(ImageKeywords {
            own: std::mem::take(&mut self.image),
            references: std::mem::take(&mut self.references),
        });
        self.image_uids.clear();
    }

    // A Reference child of the open Image element.
    pub(super) fn add_reference(&mut self, element: &BytesStart) -> Result<()> {
        if let Some(uid) = attribute(element, b"ref")? {
            if self.image_uids.insert(uid.clone()) {
                let position = self.image.len();
                self.references.push(Reference {
                    position,
                    uid: uid.clone(),
                });
            }
            self.referenced.insert(uid);
        }
        Ok(())
    }

    // A Reference element anywhere else, as in a Thumbnail element.
    pub(super) fn mark_referenced(&mut self, element: &BytesStart) -> Result<()> {
        if let Some(uid) = attribute(element, b"ref")? {
            self.referenced.insert(uid);
        }
        Ok(())
    }

    pub(super) fn add_root_keyword(&mut self, element: &BytesStart) -> Result<()> {
        if let Some(keyword) = fits_keyword_element(element)? {
            self.root.push((attribute(element, b"uid")?, keyword));
        }
        Ok(())
    }

    // Returns the keywords of each image and the root keywords. A root
    // keyword that no Reference element names applies to every image that
    // has no keyword of its name, so a lookup of the first card and of the
    // last card of a name agree. COMMENT and HISTORY cards are always added.
    // `header_len` is the length of the XML header in bytes.
    pub(super) fn finish(self, header_len: usize) -> Result<(Vec<Vec<Keyword>>, Vec<Keyword>)> {
        let root = self.split_root();
        let referenced_copies = self
            .images
            .iter()
            .flat_map(|i| &i.references)
            .fold(0usize, |n, r| {
                n.saturating_add(root.by_uid.get(r.uid.as_str()).map_or(0, Vec::len))
            });
        check_copies(
            self.images.len(),
            root.shared.len(),
            referenced_copies,
            header_len,
        )?;
        let images = self
            .images
            .iter()
            .map(|image| image.resolve(&root))
            .collect();
        let keywords = self.root.into_iter().map(|(_, keyword)| keyword).collect();
        Ok((images, keywords))
    }

    fn split_root(&self) -> RootKeywords<'_> {
        let mut by_uid: HashMap<&str, Vec<&Keyword>> = HashMap::new();
        let mut shared = Vec::new();
        for (uid, keyword) in &self.root {
            match uid {
                Some(uid) if self.referenced.contains(uid) => {
                    by_uid.entry(uid.as_str()).or_default().push(keyword)
                }
                _ => shared.push(keyword),
            }
        }
        RootKeywords { by_uid, shared }
    }
}

// Each image gets its own copy of a root keyword, while the header holds the
// keyword once, so a small header could make the reader copy keywords
// without bound. A file asking for more copies than its header has bytes is
// refused; a real file asks for far fewer.
const MAX_KEYWORD_COPIES: usize = 100_000;

fn check_copies(
    images: usize,
    shared: usize,
    referenced_copies: usize,
    header_len: usize,
) -> Result<()> {
    let copies = images
        .saturating_mul(shared)
        .saturating_add(referenced_copies);
    // A copied Keyword takes about 100 bytes, so a bound by header length
    // alone lets a 4 MB header allocate about 400 MB. The fixed limit keeps
    // the copies near 10 MB; real files hold one image and tens of keywords.
    let bound = if header_len < MAX_KEYWORD_COPIES {
        format!("the {header_len} bytes of the header")
    } else {
        format!("the limit of {MAX_KEYWORD_COPIES}")
    };
    if copies <= header_len.min(MAX_KEYWORD_COPIES) {
        return Ok(());
    }
    Err(XisfError::InvalidFormat(format!(
        "{images} images and {shared} root FITSKeyword elements without a Reference, plus \
         {referenced_copies} referenced copies, give {copies} keyword copies, more than {bound}"
    )))
}

impl ImageKeywords {
    fn resolve(&self, root: &RootKeywords) -> Vec<Keyword> {
        let mut keywords: Vec<&Keyword> = Vec::new();
        let mut references = self.references.iter().peekable();
        for position in 0..=self.own.len() {
            while let Some(r) = references.next_if(|r| r.position == position) {
                let named = root.by_uid.get(r.uid.as_str());
                keywords.extend(named.into_iter().flatten().copied());
            }
            keywords.extend(self.own.get(position));
        }
        let mut names: HashSet<&str> = keywords.iter().map(|k| lookup_name(&k.name)).collect();
        for &keyword in &root.shared {
            if keyword.is_commentary() || names.insert(lookup_name(&keyword.name)) {
                keywords.push(keyword);
            }
        }
        keywords.into_iter().cloned().collect()
    }
}

// `find_keyword` matches a name and its HIERARCH alias, so both spellings
// give the same lookup name here: a non-standard name without the HIERARCH
// prefix.
fn lookup_name(name: &str) -> &str {
    match name.strip_prefix("HIERARCH ") {
        Some(rest) if !is_standard_name(rest) => rest,
        _ => name,
    }
}

fn attribute(element: &BytesStart, name: &[u8]) -> Result<Option<String>> {
    for attr in element.attributes() {
        let attr = attr.map_err(|e| XisfError::XmlParse(e.to_string()))?;
        if attr.key.as_ref() == name {
            return attribute_text(&attr).map(Some);
        }
    }
    Ok(None)
}

// Returns the keyword a FITSKeyword element holds, or None when the element
// has no name attribute.
pub(super) fn fits_keyword_element(element: &BytesStart) -> Result<Option<Keyword>> {
    let mut name = None;
    let mut value = String::new();
    let mut comment = String::new();
    for attr in element.attributes() {
        let attr = attr.map_err(|e| XisfError::XmlParse(e.to_string()))?;
        let text = printable(&attribute_text(&attr)?);
        match attr.key.as_ref() {
            b"name" => name = Some(text),
            b"value" => value = text,
            b"comment" => comment = text,
            _ => {}
        }
    }
    name.map(|name| fits_keyword(name, &value, comment))
        .transpose()
}

// The FITS header reader maps bytes outside 0x20-0x7E to '?', so a keyword
// read here holds text that a FITS header can hold too.
fn printable(text: &str) -> String {
    text.bytes()
        .map(|b| {
            if (b' '..=b'~').contains(&b) {
                char::from(b)
            } else {
                '?'
            }
        })
        .collect()
}

// XISF 1.0 section 11.6.1 gives COMMENT and HISTORY an empty value, and a
// blank FITS name has no value either. Any other empty value is the FITS
// undefined value (FITS 4.0 section 4.1.2.3). The value is FITS value field
// text, so it is typed by the FITS header parser.
fn fits_keyword(name: String, value: &str, comment: String) -> Result<Keyword> {
    let comment = (!comment.is_empty()).then_some(comment);
    if matches!(name.as_str(), "COMMENT" | "HISTORY" | "") {
        let mut keyword = Keyword::commentary(name, String::new());
        keyword.comment = comment;
        return Ok(keyword);
    }
    let mut keyword = Keyword::new(name);
    keyword.comment = comment;
    if !value.trim().is_empty() {
        let parsed = HeaderCard::parse_value(value).map_err(XisfError::Keyword)?;
        keyword.value = Some(parsed);
    }
    Ok(keyword)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keyword_copies_up_to_the_header_length_are_accepted() {
        assert!(check_copies(10, 10, 0, 100).is_ok());
        assert!(check_copies(10, 9, 10, 100).is_ok());
        assert!(check_copies(0, usize::MAX, 0, 0).is_ok());
    }

    #[test]
    fn keyword_copies_past_the_header_length_are_refused() {
        let err = check_copies(10, 10, 1, 100).unwrap_err();
        assert_eq!(
            err.to_string(),
            "Invalid XISF format: 10 images and 10 root FITSKeyword elements without a \
             Reference, plus 1 referenced copies, give 101 keyword copies, more than the 100 \
             bytes of the header"
        );
        assert!(check_copies(usize::MAX, 2, 0, usize::MAX - 1).is_err());
        assert!(check_copies(1, 0, usize::MAX, usize::MAX - 1).is_err());
    }

    #[test]
    fn keyword_copies_past_the_fixed_limit_are_refused_in_a_long_header() {
        assert!(check_copies(1000, 100, 0, 4_000_000).is_ok());
        let err = check_copies(1000, 100, 1, 4_000_000).unwrap_err();
        assert_eq!(
            err.to_string(),
            "Invalid XISF format: 1000 images and 100 root FITSKeyword elements without a \
             Reference, plus 1 referenced copies, give 100001 keyword copies, more than the \
             limit of 100000"
        );
    }

    #[test]
    fn both_spellings_of_a_hierarch_name_have_one_lookup_name() {
        assert_eq!(lookup_name("HIERARCH LONG-NAME"), "LONG-NAME");
        assert_eq!(lookup_name("LONG-NAME"), "LONG-NAME");
        assert_eq!(lookup_name("HIERARCH NAXIS1"), "HIERARCH NAXIS1");
        assert_eq!(lookup_name("NAXIS1"), "NAXIS1");
    }
}
