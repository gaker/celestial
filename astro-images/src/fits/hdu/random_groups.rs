use super::{HduTrait, HduType};
use crate::fits::header::Header;
use crate::fits::io::reader::HduInfo;

#[derive(Debug)]
pub struct RandomGroupsHdu {
    header: Header,
    info: HduInfo,
}

impl RandomGroupsHdu {
    pub fn new(header: Header, info: HduInfo) -> Self {
        Self { header, info }
    }

    pub fn group_count(&self) -> Option<i64> {
        self.header
            .get_keyword_value("GCOUNT")
            .and_then(|v| v.as_integer())
    }

    pub fn parameter_count(&self) -> Option<i64> {
        self.header
            .get_keyword_value("PCOUNT")
            .and_then(|v| v.as_integer())
    }

    pub fn extension_name(&self) -> Option<&str> {
        self.header
            .get_keyword_value("EXTNAME")
            .and_then(|v| v.as_string())
    }

    pub fn extension_version(&self) -> Option<i64> {
        self.header
            .get_keyword_value("EXTVER")
            .and_then(|v| v.as_integer())
    }
}

impl HduTrait for RandomGroupsHdu {
    fn header(&self) -> &Header {
        &self.header
    }

    fn info(&self) -> &HduInfo {
        &self.info
    }

    fn hdu_type(&self) -> HduType {
        HduType::RandomGroups
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fits::header::{Header, Keyword};
    use crate::fits::io::reader::HduInfo;

    fn create_test_header(extname: Option<&str>, groups: bool) -> Header {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 3));
        header.add_keyword(Keyword::integer("NAXIS1", 0));
        header.add_keyword(Keyword::integer("NAXIS2", 10));
        header.add_keyword(Keyword::integer("NAXIS3", 5));
        header.add_keyword(Keyword::integer("BITPIX", -32));

        if groups {
            header.add_keyword(Keyword::logical("GROUPS", true));
            header.add_keyword(Keyword::integer("GCOUNT", 10));
            header.add_keyword(Keyword::integer("PCOUNT", 3));
        }

        if let Some(name) = extname {
            header.add_keyword(Keyword::string("EXTNAME", name));
            header.add_keyword(Keyword::integer("EXTVER", 1));
        }

        header
    }

    fn create_test_hdu_info() -> HduInfo {
        HduInfo {
            index: 0,
            header_start: 0,
            header_size: 2880,
            data_start: 2880,
            data_size: 1200,
        }
    }

    #[test]
    fn new_creates_random_groups_hdu() {
        let header = create_test_header(Some("UVDATA"), true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.info.index, 0);
        assert_eq!(hdu.group_count(), Some(10));
    }

    #[test]
    fn header_returns_header_reference() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        let header_ref = hdu.header();
        assert!(header_ref
            .get_keyword_value("SIMPLE")
            .unwrap()
            .as_logical()
            .unwrap());
    }

    #[test]
    fn info_returns_info_reference() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        let info_ref = hdu.info();
        assert_eq!(info_ref.index, 0);
        assert_eq!(info_ref.data_start, 2880);
    }

    #[test]
    fn hdu_type_returns_random_groups() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::RandomGroups);
    }

    #[test]
    fn group_count_returns_gcount_value() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.group_count(), Some(10));
    }

    #[test]
    fn group_count_returns_none_when_missing() {
        let header = create_test_header(None, false);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.group_count(), None);
    }

    #[test]
    fn parameter_count_returns_pcount_value() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.parameter_count(), Some(3));
    }

    #[test]
    fn parameter_count_returns_none_when_missing() {
        let header = create_test_header(None, false);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.parameter_count(), None);
    }

    #[test]
    fn extension_name_returns_extname_value() {
        let header = create_test_header(Some("VISIBILITY"), true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.extension_name(), Some("VISIBILITY"));
    }

    #[test]
    fn extension_name_returns_none_when_missing() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.extension_name(), None);
    }

    #[test]
    fn extension_version_returns_extver_value() {
        let header = create_test_header(Some("TEST"), true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.extension_version(), Some(1));
    }

    #[test]
    fn extension_version_returns_none_when_missing() {
        let header = create_test_header(None, true);
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.extension_version(), None);
    }

    #[test]
    fn all_methods_work_together() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 4));
        header.add_keyword(Keyword::integer("NAXIS1", 0));
        header.add_keyword(Keyword::integer("NAXIS2", 2));
        header.add_keyword(Keyword::integer("NAXIS3", 1024));
        header.add_keyword(Keyword::integer("NAXIS4", 1024));
        header.add_keyword(Keyword::integer("BITPIX", -64));
        header.add_keyword(Keyword::logical("GROUPS", true));
        header.add_keyword(Keyword::integer("GCOUNT", 100));
        header.add_keyword(Keyword::integer("PCOUNT", 5));
        header.add_keyword(Keyword::string("EXTNAME", "RADIODATA"));
        header.add_keyword(Keyword::integer("EXTVER", 2));
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::RandomGroups);
        assert_eq!(hdu.group_count(), Some(100));
        assert_eq!(hdu.parameter_count(), Some(5));
        assert_eq!(hdu.extension_name(), Some("RADIODATA"));
        assert_eq!(hdu.extension_version(), Some(2));
    }

    #[test]
    fn minimal_valid_random_groups() {
        let mut header = Header::new();
        header.add_keyword(Keyword::logical("SIMPLE", true));
        header.add_keyword(Keyword::integer("NAXIS", 1));
        header.add_keyword(Keyword::integer("NAXIS1", 0));
        header.add_keyword(Keyword::integer("BITPIX", 16));
        header.add_keyword(Keyword::logical("GROUPS", true));
        header.add_keyword(Keyword::integer("GCOUNT", 1));
        header.add_keyword(Keyword::integer("PCOUNT", 0));
        let info = create_test_hdu_info();
        let hdu = RandomGroupsHdu::new(header, info);

        assert_eq!(hdu.hdu_type(), HduType::RandomGroups);
        assert_eq!(hdu.group_count(), Some(1));
        assert_eq!(hdu.parameter_count(), Some(0));
        assert_eq!(hdu.extension_name(), None);
        assert_eq!(hdu.extension_version(), None);
    }
}
