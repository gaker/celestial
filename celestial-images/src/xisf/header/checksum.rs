use crate::xisf::errors::{Result, XisfError};
use crate::xisf::header::text_encoding::decode_hex_digits;
use sha1::Sha1;
use sha2::{Digest, Sha256, Sha512};

// The checksum attribute of a data block (XISF 1.0 section 10.5). The
// digest covers the block as stored: the compressed bytes of a compressed
// block (section 10.6.2) and the decoded bytes of an inline or embedded one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Checksum {
    algorithm: Algorithm,
    digest: Vec<u8>,
}

// Section 10.5 requires decoders to support SHA-1, SHA-256 and SHA-512. The
// optional SHA3-256 and SHA3-512 are not supported.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Algorithm {
    Sha1,
    Sha256,
    Sha512,
}

impl Algorithm {
    fn parse(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "sha-1" | "sha1" => Some(Self::Sha1),
            "sha-256" | "sha256" => Some(Self::Sha256),
            "sha-512" | "sha512" => Some(Self::Sha512),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Sha1 => "sha-1",
            Self::Sha256 => "sha-256",
            Self::Sha512 => "sha-512",
        }
    }

    fn digest_len(self) -> usize {
        match self {
            Self::Sha1 => 20,
            Self::Sha256 => 32,
            Self::Sha512 => 64,
        }
    }

    fn digest(self, data: &[u8]) -> Vec<u8> {
        match self {
            Self::Sha1 => Sha1::digest(data).to_vec(),
            Self::Sha256 => Sha256::digest(data).to_vec(),
            Self::Sha512 => Sha512::digest(data).to_vec(),
        }
    }
}

impl Checksum {
    // A value is algorithm:digest, with an algorithm name of letters, digits
    // and hyphens (Table 9) and a digest of hexadecimal digits. A value of
    // that form with a name this reader does not support is unsupported; any
    // other value, including one with whitespace, is invalid.
    pub(crate) fn parse(value: &str) -> Result<Self> {
        let invalid = |problem: String| {
            XisfError::InvalidFormat(format!("Invalid checksum attribute \"{value}\": {problem}"))
        };
        let not_algorithm_digest = || invalid("it is not algorithm:digest".to_string());
        let (name, digits) = value.split_once(':').ok_or_else(not_algorithm_digest)?;
        let is_name = |c: char| c.is_ascii_alphanumeric() || c == '-';
        if name.is_empty() || !name.chars().all(is_name) {
            return Err(not_algorithm_digest());
        }
        let Some(algorithm) = Algorithm::parse(name) else {
            if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_hexdigit()) {
                return Err(not_algorithm_digest());
            }
            return Err(XisfError::UnsupportedChecksum(name.to_string()));
        };
        let wrong_digest = || {
            let digits = algorithm.digest_len() * 2;
            invalid(format!(
                "a {} digest is {digits} hexadecimal digits",
                algorithm.name()
            ))
        };
        let digest = decode_hex_digits(digits).map_err(|_| wrong_digest())?;
        if digest.len() != algorithm.digest_len() {
            return Err(wrong_digest());
        }
        Ok(Self { algorithm, digest })
    }

    // Section 10.5 requires SHA-1 of every encoder that writes checksums.
    pub(crate) fn sha1(stored: &[u8]) -> Self {
        let algorithm = Algorithm::Sha1;
        Self {
            algorithm,
            digest: algorithm.digest(stored),
        }
    }

    pub(crate) fn verify(&self, stored: &[u8]) -> Result<()> {
        let actual = self.algorithm.digest(stored);
        if actual == self.digest {
            return Ok(());
        }
        Err(XisfError::ChecksumMismatch {
            algorithm: self.algorithm.name(),
            expected: lowercase_hex(&self.digest),
            actual: lowercase_hex(&actual),
        })
    }
}

// The checksum attribute value, with the digest in lowercase Base16 as
// section 10.5 requires.
impl std::fmt::Display for Checksum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}",
            self.algorithm.name(),
            lowercase_hex(&self.digest)
        )
    }
}

fn lowercase_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // FIPS 180 test vectors for the message "abc".
    const SHA1_ABC: &str = "a9993e364706816aba3e25717850c26c9cd0d89d";
    const SHA256_ABC: &str = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
    const SHA512_ABC: &str = "ddaf35a193617abacc417349ae20413112e6fa4e89a97ea20a9eeee64b55d39a2192992a274fc1a836ba3c23a3feebbd454d4423643ce80e2a9ac94fa54ca49f";

    #[test]
    fn fips_180_vectors_verify() {
        for (name, digest) in [
            ("sha1", SHA1_ABC),
            ("sha256", SHA256_ABC),
            ("sha512", SHA512_ABC),
        ] {
            let checksum = Checksum::parse(&format!("{name}:{digest}")).unwrap();
            assert!(checksum.verify(b"abc").is_ok());
            assert!(matches!(
                checksum.verify(b"abd"),
                Err(XisfError::ChecksumMismatch { expected, .. }) if expected == digest
            ));
        }
    }

    #[test]
    fn sha1_checksum_formats_as_an_attribute_value_that_parses_back() {
        let checksum = Checksum::sha1(b"abc");
        assert_eq!(checksum.to_string(), format!("sha-1:{SHA1_ABC}"));
        assert_eq!(Checksum::parse(&checksum.to_string()).unwrap(), checksum);
    }

    #[test]
    fn algorithm_names_ignore_case() {
        let checksum = Checksum::parse(&format!("SHA-1:{SHA1_ABC}")).unwrap();
        assert_eq!(checksum.algorithm, Algorithm::Sha1);
    }

    // Whitespace or any other malformed shape is invalid. Only a well-formed
    // name:digest with a name this reader does not support is unsupported.
    #[test]
    fn malformed_values_are_invalid_and_unknown_names_unsupported() {
        for value in [
            format!(" sha1:{SHA1_ABC}"),
            format!("sha1 :{SHA1_ABC}"),
            format!("sha1: {SHA1_ABC}"),
            format!("sha1:{SHA1_ABC} "),
            "".to_string(),
            ":00".to_string(),
            "sha1".to_string(),
            "sha 1:00".to_string(),
            "md5:".to_string(),
            "md5: 00".to_string(),
            "md5:zz".to_string(),
            " md5:00".to_string(),
        ] {
            let result = Checksum::parse(&value);
            assert!(
                matches!(result, Err(XisfError::InvalidFormat(_))),
                "{value:?}: {result:?}"
            );
        }
        for (value, name) in [
            ("md5:00", "md5"),
            ("sha3-512:abc", "sha3-512"),
            ("SHA3-256:00", "SHA3-256"),
        ] {
            assert!(matches!(
                Checksum::parse(value),
                Err(XisfError::UnsupportedChecksum(n)) if n == name
            ));
        }
    }

    #[test]
    fn sha3_is_unsupported() {
        assert!(matches!(
            Checksum::parse("sha3-256:00"),
            Err(XisfError::UnsupportedChecksum(name)) if name == "sha3-256"
        ));
    }
}
