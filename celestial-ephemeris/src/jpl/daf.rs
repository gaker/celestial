use super::SpkError;
use std::collections::HashSet;
use std::io::{ErrorKind, Read};

const RECORD_BYTES: usize = 1024;
const RECORD_WORDS: usize = 128;
const SUMMARY_WORDS: usize = 5;
const MAX_SUMMARIES: usize = (RECORD_WORDS - 3) / SUMMARY_WORDS;
const CHUNK_BYTES: usize = 64 * 1024;
const FTP_MARKER: &[u8] = b"FTPSTR:";
const FTP_STRING: &[u8] = b"FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP";

pub(super) struct Daf {
    pub(super) words: Vec<f64>,
    pub(super) summaries: Vec<Summary>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Summary {
    pub(super) start: f64,
    pub(super) end: f64,
    pub(super) body: i32,
    pub(super) center: i32,
    pub(super) frame: i32,
    pub(super) data_type: i32,
    pub(super) begin: i32,
    pub(super) end_address: i32,
}

#[derive(Clone, Copy)]
enum Endian {
    Little,
    Big,
}

impl Endian {
    fn word(self, bytes: [u8; 8]) -> f64 {
        match self {
            Self::Little => f64::from_le_bytes(bytes),
            Self::Big => f64::from_be_bytes(bytes),
        }
    }

    fn int(self, bytes: [u8; 4]) -> i32 {
        match self {
            Self::Little => i32::from_le_bytes(bytes),
            Self::Big => i32::from_be_bytes(bytes),
        }
    }

    // A summary packs two i32s into each double word, in file byte order.
    fn ints(self, word: f64) -> [i32; 2] {
        let bits = word.to_bits();
        let (high, low) = ((bits >> 32) as u32 as i32, bits as u32 as i32);
        match self {
            Self::Little => [low, high],
            Self::Big => [high, low],
        }
    }
}

impl Daf {
    pub(super) fn read(mut reader: impl Read, size_hint: u64) -> Result<Self, SpkError> {
        let mut record = [0u8; RECORD_BYTES];
        reader.read_exact(&mut record).map_err(short_file)?;
        let (endian, fward) = parse_file_record(&record)?;
        let mut words = Vec::new();
        let capacity = usize::try_from(size_hint / 8).unwrap_or(usize::MAX);
        words
            .try_reserve_exact(capacity)
            .map_err(|e| std::io::Error::new(ErrorKind::OutOfMemory, e))?;
        push_words(&mut words, endian, &record);
        read_words(&mut reader, endian, &mut words)?;
        let summaries = Walker {
            words: &words,
            endian,
        }
        .summaries(fward)?;
        Ok(Self { words, summaries })
    }
}

pub(super) fn is_whole(x: f64) -> bool {
    libm::trunc(x) == x
}

fn short_file(error: std::io::Error) -> SpkError {
    if error.kind() == ErrorKind::UnexpectedEof {
        return SpkError::InvalidFormat("the file is shorter than one 1024-byte DAF record".into());
    }
    error.into()
}

fn bytes<const N: usize>(record: &[u8; RECORD_BYTES], at: usize) -> [u8; N] {
    let mut out = [0; N];
    out.copy_from_slice(&record[at..at + N]);
    out
}

fn parse_file_record(record: &[u8; RECORD_BYTES]) -> Result<(Endian, i32), SpkError> {
    check_id_word(&record[..8])?;
    let endian = binary_format(record)?;
    let (nd, ni) = (endian.int(bytes(record, 8)), endian.int(bytes(record, 12)));
    if (nd, ni) != (2, 6) {
        return Err(SpkError::InvalidFormat(format!(
            "summary format ND = {}, NI = {}; an SPK has ND = 2, NI = 6",
            nd, ni
        )));
    }
    check_ftp(record)?;
    Ok((endian, endian.int(bytes(record, 76))))
}

fn check_id_word(id: &[u8]) -> Result<(), SpkError> {
    if id == b"DAF/SPK " || id == b"NAIF/DAF" {
        return Ok(());
    }
    Err(SpkError::InvalidFormat(format!(
        "ID word {:?} is not DAF/SPK or NAIF/DAF",
        String::from_utf8_lossy(id)
    )))
}

fn binary_format(record: &[u8; RECORD_BYTES]) -> Result<Endian, SpkError> {
    match &record[88..96] {
        b"LTL-IEEE" => Ok(Endian::Little),
        b"BIG-IEEE" => Ok(Endian::Big),
        // Files written before the format field existed leave it blank; an
        // SPK's ND is always 2, so its byte order shows which one was used.
        id if id.iter().all(|&b| b == b' ' || b == 0) => {
            if i32::from_le_bytes(bytes(record, 8)) == 2 {
                Ok(Endian::Little)
            } else {
                Ok(Endian::Big)
            }
        }
        id => Err(SpkError::InvalidFormat(format!(
            "binary format {:?} is not LTL-IEEE or BIG-IEEE",
            String::from_utf8_lossy(id)
        ))),
    }
}

// Files older than the FTP string have none, so only a damaged one is an error.
fn check_ftp(record: &[u8; RECORD_BYTES]) -> Result<(), SpkError> {
    let Some(at) = record
        .windows(FTP_MARKER.len())
        .position(|w| w == FTP_MARKER)
    else {
        return Ok(());
    };
    if record[at..].starts_with(FTP_STRING) {
        return Ok(());
    }
    Err(SpkError::InvalidFormat(
        "FTP validation string is damaged; the file was probably transferred in text mode".into(),
    ))
}

fn push_words(words: &mut Vec<f64>, endian: Endian, bytes: &[u8]) {
    let (whole, _) = bytes.as_chunks::<8>();
    words.extend(whole.iter().map(|&w| endian.word(w)));
}

fn read_words(
    reader: &mut impl Read,
    endian: Endian,
    words: &mut Vec<f64>,
) -> Result<(), SpkError> {
    let mut chunk = vec![0u8; CHUNK_BYTES];
    let mut filled = 0;
    loop {
        match reader.read(&mut chunk[filled..]) {
            Ok(0) => return Ok(()),
            Ok(n) => filled += n,
            Err(e) if e.kind() == ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
        let whole = filled - filled % 8;
        push_words(words, endian, &chunk[..whole]);
        chunk.copy_within(whole..filled, 0);
        filled -= whole;
    }
}

struct Walker<'a> {
    words: &'a [f64],
    endian: Endian,
}

impl Walker<'_> {
    fn summaries(&self, fward: i32) -> Result<Vec<Summary>, SpkError> {
        let mut seen = HashSet::new();
        let mut summaries = Vec::new();
        let mut number = self.record_number(f64::from(fward), "first summary record")?;
        loop {
            if !seen.insert(number) {
                return Err(SpkError::InvalidFormat(format!(
                    "summary record {} is reached twice; the record chain loops",
                    number
                )));
            }
            let next = self.summary_record(number, &mut summaries)?;
            if next == 0.0 {
                return Ok(summaries);
            }
            number = self.record_number(next, "next summary record")?;
        }
    }

    // Returns the record's pointer to the next one, 0 for the last.
    fn summary_record(&self, number: usize, summaries: &mut Vec<Summary>) -> Result<f64, SpkError> {
        let record = &self.words[(number - 1) * RECORD_WORDS..number * RECORD_WORDS];
        let count = summary_count(record[2], number)?;
        let packed = record[3..]
            .as_chunks::<SUMMARY_WORDS>()
            .0
            .iter()
            .take(count);
        summaries.extend(packed.map(|words| self.summary(words)));
        Ok(record[0])
    }

    fn record_number(&self, value: f64, what: &str) -> Result<usize, SpkError> {
        let last = self.words.len() / RECORD_WORDS;
        if is_whole(value) && value >= 2.0 && value <= last as f64 {
            return Ok(value as usize);
        }
        Err(SpkError::InvalidFormat(format!(
            "{} {} is not a record of this file (2 to {})",
            what, value, last
        )))
    }

    fn summary(&self, &[start, end, ids, kinds, range]: &[f64; SUMMARY_WORDS]) -> Summary {
        let [body, center] = self.endian.ints(ids);
        let [frame, data_type] = self.endian.ints(kinds);
        let [begin, end_address] = self.endian.ints(range);
        Summary {
            start,
            end,
            body,
            center,
            frame,
            data_type,
            begin,
            end_address,
        }
    }
}

fn summary_count(nsum: f64, record: usize) -> Result<usize, SpkError> {
    if is_whole(nsum) && nsum >= 0.0 && nsum <= MAX_SUMMARIES as f64 {
        return Ok(nsum as usize);
    }
    Err(SpkError::InvalidFormat(format!(
        "summary record {} holds {} summaries; expected a whole number from 0 to {}",
        record, nsum, MAX_SUMMARIES
    )))
}
