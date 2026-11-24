//! Lists the Gaia DR3 source files in the S3-compatible bucket that mirrors the ESA CDN.
//!
//! The bucket answers with S3 ListObjects (version 1) pages. A page with
//! `<IsTruncated>true</IsTruncated>` continues at its `NextMarker`, or at its last key when it
//! has no `NextMarker`.

use crate::http::Downloader;
use anyhow::Context;
use quick_xml::events::{BytesRef, Event};
use quick_xml::Reader;

const LISTING_URL: &str =
    "https://gaia.eu-1.cdn77-storage.com/?prefix=Gaia/gdr3/gaia_source/&delimiter=/";

pub(crate) struct RemoteFile {
    pub(crate) key: String,
    pub(crate) filename: String,
    pub(crate) size: u64,
    pub(crate) etag: String,
}

struct ListingPage {
    files: Vec<RemoteFile>,
    // None when this page ends the listing.
    next_marker: Option<String>,
}

pub(crate) fn list_remote_files(
    http: &Downloader,
    verbose: bool,
) -> anyhow::Result<Vec<RemoteFile>> {
    let mut files = Vec::new();
    let mut marker: Option<String> = None;
    loop {
        let url = listing_url(marker.as_deref())?;
        if verbose {
            eprintln!("Listing: {}", url);
        }
        let page = parse_listing(&http.get_text(&url)?)?;
        files.extend(page.files);
        if verbose {
            eprintln!("  Total keys so far: {}", files.len());
        }
        match page.next_marker {
            Some(next) if marker.as_ref() == Some(&next) => {
                anyhow::bail!("The bucket listing repeated marker {:?}", next)
            }
            Some(next) => marker = Some(next),
            None => return Ok(files),
        }
    }
}

fn listing_url(marker: Option<&str>) -> anyhow::Result<String> {
    let mut url = reqwest::Url::parse(LISTING_URL).context("Invalid listing URL")?;
    if let Some(marker) = marker {
        url.query_pairs_mut().append_pair("marker", marker);
    }
    Ok(url.into())
}

#[derive(Default)]
struct ListingParser {
    // Names of the open elements, outermost first.
    open: Vec<String>,
    // Text of the innermost open element.
    text: String,
    key: Option<String>,
    size: Option<u64>,
    etag: String,
    last_key: Option<String>,
    files: Vec<RemoteFile>,
    next_marker: Option<String>,
    is_truncated: Option<bool>,
    complete: bool,
}

fn parse_listing(xml: &str) -> anyhow::Result<ListingPage> {
    let mut reader = Reader::from_str(xml);
    let mut parser = ListingParser::default();
    loop {
        match reader
            .read_event()
            .context("The bucket listing is not valid XML")?
        {
            Event::Start(e) => parser.open_element(element_name(e.name().as_ref())?),
            Event::End(_) => parser.close_element()?,
            Event::Text(e) => parser.text.push_str(&e.decode()?),
            Event::GeneralRef(e) => parser.text.push(resolve_reference(&e)?),
            Event::CData(_) => parser.reject_cdata()?,
            Event::Eof => return parser.finish(),
            _ => {}
        }
    }
}

fn element_name(raw: &[u8]) -> anyhow::Result<String> {
    Ok(std::str::from_utf8(raw)?.to_string())
}

// S3 escapes the quotes around an ETag as `&quot;`.
fn resolve_reference(reference: &BytesRef) -> anyhow::Result<char> {
    if let Some(c) = reference.resolve_char_ref()? {
        return Ok(c);
    }
    let name = reference.decode()?;
    quick_xml::escape::resolve_predefined_entity(&name)
        .and_then(|s| s.chars().next())
        .with_context(|| format!("Unknown entity &{};", name))
}

impl ListingParser {
    // Text before a child element belongs to no field, so each element starts with none.
    fn open_element(&mut self, name: String) {
        self.text.clear();
        self.open.push(name);
    }

    // The parser reads no CDATA, so CDATA in an element it reads would leave that field empty.
    fn reject_cdata(&self) -> anyhow::Result<()> {
        match self.open.last().map(String::as_str) {
            Some(name @ ("Key" | "Size" | "ETag" | "NextMarker" | "IsTruncated")) => {
                anyhow::bail!("The bucket listing has CDATA in {}", name)
            }
            _ => Ok(()),
        }
    }

    fn close_element(&mut self) -> anyhow::Result<()> {
        let name = self.open.pop().context("Unbalanced end tag")?;
        let parent = self.open.last().map(String::as_str);
        let text = std::mem::take(&mut self.text);
        match (parent, name.as_str()) {
            (Some("Contents"), "Key") => self.key = Some(text),
            (Some("Contents"), "Size") => self.size = Some(parse_size(&text)?),
            (Some("Contents"), "ETag") => self.etag = text.trim().trim_matches('"').to_string(),
            (Some(_), "Contents") => self.close_contents()?,
            (Some(_), "NextMarker") => self.next_marker = Some(text.trim().to_string()),
            (Some(_), "IsTruncated") => self.is_truncated = Some(parse_bool(&text)?),
            (None, "ListBucketResult") => self.complete = true,
            _ => {}
        }
        Ok(())
    }

    fn close_contents(&mut self) -> anyhow::Result<()> {
        let key = self.key.take().context("A Contents element has no Key")?;
        let size = self
            .size
            .take()
            .with_context(|| format!("{} has no Size", key))?;
        let etag = std::mem::take(&mut self.etag);
        if key.ends_with(".csv.gz") {
            let filename = key.rsplit('/').next().unwrap_or(&key).to_string();
            let file = RemoteFile {
                key: key.clone(),
                filename,
                size,
                etag,
            };
            self.files.push(file);
        }
        self.last_key = Some(key);
        Ok(())
    }

    fn finish(self) -> anyhow::Result<ListingPage> {
        if !self.complete {
            anyhow::bail!("The bucket listing ends before </ListBucketResult>");
        }
        let next_marker =
            match self.is_truncated {
                None => anyhow::bail!("The bucket listing has no IsTruncated element"),
                Some(false) => None,
                Some(true) => Some(self.next_marker.or(self.last_key).context(
                    "The bucket listing is truncated but gives no NextMarker and no keys",
                )?),
            };
        Ok(ListingPage {
            files: self.files,
            next_marker,
        })
    }
}

fn parse_size(text: &str) -> anyhow::Result<u64> {
    text.trim()
        .parse()
        .with_context(|| format!("Invalid object size {:?} in the bucket listing", text))
}

fn parse_bool(text: &str) -> anyhow::Result<bool> {
    match text.trim() {
        "true" => Ok(true),
        "false" => Ok(false),
        other => anyhow::bail!("Invalid IsTruncated value {:?}", other),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PREFIX: &str = "Gaia/gdr3/gaia_source/";

    fn contents(name: &str, size: &str) -> String {
        format!(
            "<Contents>\n  <Key>{PREFIX}{name}</Key>\n  <Size>{size}</Size>\n  \
             <ETag>&quot;abc{name}&quot;</ETag>\n</Contents>\n"
        )
    }

    fn listing(body: &str) -> String {
        format!(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
             <ListBucketResult xmlns=\"http://s3.amazonaws.com/doc/2006-03-01/\">\n\
             <Name>gaia</Name>\n{body}</ListBucketResult>\n"
        )
    }

    fn parse_err(xml: &str) -> String {
        match parse_listing(xml) {
            Ok(page) => panic!("parsed {} files", page.files.len()),
            Err(e) => format!("{:#}", e),
        }
    }

    #[test]
    fn test_listing_page_is_parsed() {
        let body = contents("a.csv.gz", "12") + &contents("_MD5SUM.txt", "5");
        let xml = listing(&(body + "<IsTruncated>false</IsTruncated>\n"));
        let page = parse_listing(&xml).unwrap();
        assert_eq!(page.files.len(), 1);
        let file = &page.files[0];
        assert_eq!(file.key, format!("{PREFIX}a.csv.gz"));
        assert_eq!(file.filename, "a.csv.gz");
        assert_eq!(file.size, 12);
        assert_eq!(file.etag, "abca.csv.gz");
        assert_eq!(page.next_marker, None);
    }

    #[test]
    fn test_unparseable_size_is_an_error() {
        let xml = listing(&(contents("a.csv.gz", "12x") + "<IsTruncated>false</IsTruncated>"));
        assert!(parse_err(&xml).contains("Invalid object size \"12x\""));
    }

    #[test]
    fn test_truncated_page_continues_at_next_marker() {
        let tail = "<IsTruncated>true</IsTruncated><NextMarker>m/x</NextMarker>";
        let page = parse_listing(&listing(&(contents("a.csv.gz", "1") + tail))).unwrap();
        assert_eq!(page.next_marker.as_deref(), Some("m/x"));
    }

    #[test]
    fn test_truncated_page_without_next_marker_continues_at_last_key() {
        let body = contents("a.csv.gz", "1") + &contents("b.txt", "1");
        let page = parse_listing(&listing(&(body + "<IsTruncated>true</IsTruncated>"))).unwrap();
        assert_eq!(page.next_marker, Some(format!("{PREFIX}b.txt")));
    }

    #[test]
    fn test_truncated_page_with_no_keys_is_an_error() {
        let xml = listing("<IsTruncated>true</IsTruncated>");
        assert!(parse_err(&xml).contains("no NextMarker and no keys"));
    }

    #[test]
    fn test_page_without_is_truncated_is_an_error() {
        let xml = listing(&contents("a.csv.gz", "1"));
        assert!(parse_err(&xml).contains("no IsTruncated"));
    }

    #[test]
    fn test_cut_off_page_is_an_error() {
        let xml = listing(&(contents("a.csv.gz", "1") + "<IsTruncated>false</IsTruncated>"));
        let cut = &xml[..xml.find("</ListBucketResult>").unwrap()];
        assert!(parse_err(cut).contains("ends before </ListBucketResult>"));
    }

    #[test]
    fn test_unknown_entity_is_an_error() {
        let body = contents("a&bogus;.csv.gz", "1") + "<IsTruncated>false</IsTruncated>";
        assert!(parse_err(&listing(&body)).contains("&bogus;"));
    }

    #[test]
    fn test_marker_is_url_encoded() {
        let url = listing_url(Some("Gaia/gdr3/a b&c=d.csv.gz")).unwrap();
        assert_eq!(
            url,
            format!("{LISTING_URL}&marker=Gaia%2Fgdr3%2Fa+b%26c%3Dd.csv.gz")
        );
        assert_eq!(listing_url(None).unwrap(), LISTING_URL);
    }

    #[test]
    fn test_cdata_key_is_an_error() {
        let body = "<Contents><Key><![CDATA[Gaia/a.csv.gz]]></Key><Size>1</Size><ETag>e</ETag>\
                    </Contents><IsTruncated>false</IsTruncated>";
        assert!(parse_err(&listing(body)).contains("CDATA in Key"));
    }

    #[test]
    fn test_cdata_size_is_an_error() {
        let body = "<Contents><Key>Gaia/a.csv.gz</Key><Size><![CDATA[1]]></Size><ETag>e</ETag>\
                    </Contents><IsTruncated>false</IsTruncated>";
        assert!(parse_err(&listing(body)).contains("CDATA in Size"));
    }
}
