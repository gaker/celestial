use super::chain::{self, Sum};
use super::daf::{Daf, Summary};
use super::type2::Type2;
use super::SpkError;
use celestial_core::constants::{J2000_JD, SECONDS_PER_DAY_F64};
use celestial_core::matrix::Vector3;
use celestial_time::julian::JulianDate;
use celestial_time::scales::tdb::TDB;
use std::fmt;
use std::fs::File;
use std::io::Read;
use std::path::Path;

const J2000_FRAME: i32 = 1;
const CHEBYSHEV_POSITION: i32 = 2;

pub struct SpkFile {
    words: Vec<f64>,
    segments: Vec<SpkSegment>,
}

#[derive(Debug, Clone)]
pub struct SpkSegment {
    index: usize,
    body: i32,
    center: i32,
    frame: i32,
    data_type: i32,
    start: f64,
    end: f64,
    type2: Option<Type2>,
}

struct Label {
    index: usize,
    body: i32,
    center: i32,
}

#[derive(Clone, Copy)]
struct Epoch {
    seconds: f64,
    jd: f64,
}

impl SpkFile {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, SpkError> {
        let file = File::open(path)?;
        let size = file.metadata()?.len();
        Self::from_reader(file, size)
    }

    pub(super) fn from_reader(reader: impl Read, size_hint: u64) -> Result<Self, SpkError> {
        let Daf { words, summaries } = Daf::read(reader, size_hint)?;
        let segments = summaries
            .iter()
            .enumerate()
            .map(|(index, summary)| SpkSegment::new(index, summary, &words))
            .collect::<Result<_, _>>()?;
        Ok(Self { words, segments })
    }

    pub fn segments(&self) -> &[SpkSegment] {
        &self.segments
    }

    pub fn compute_state(
        &self,
        body: i32,
        center: i32,
        tdb: &TDB,
    ) -> Result<(Vector3, Vector3), SpkError> {
        let epoch = Epoch::new(tdb)?;
        self.chain(body, center, epoch, |s| s.state(&self.words, epoch.seconds))
    }

    pub fn compute_position(&self, body: i32, center: i32, tdb: &TDB) -> Result<Vector3, SpkError> {
        let epoch = Epoch::new(tdb)?;
        self.chain(body, center, epoch, |s| {
            s.position(&self.words, epoch.seconds)
        })
    }

    fn chain<T: Sum>(
        &self,
        body: i32,
        center: i32,
        epoch: Epoch,
        link: impl Fn(&SpkSegment) -> Result<T, SpkError>,
    ) -> Result<T, SpkError> {
        chain::relative(&self.segments, (body, center), epoch.seconds, link)?.ok_or(
            SpkError::SegmentNotFound {
                body,
                center,
                jd: epoch.jd,
            },
        )
    }
}

impl SpkSegment {
    fn new(index: usize, summary: &Summary, words: &[f64]) -> Result<Self, SpkError> {
        let mut segment = Self {
            index,
            body: summary.body,
            center: summary.center,
            frame: summary.frame,
            data_type: summary.data_type,
            start: summary.start,
            end: summary.end,
            type2: None,
        };
        segment.type2 = segment.validate(summary, words).map_err(|reason| {
            SpkError::InvalidFormat(format!("{}: {}", segment.label(), reason))
        })?;
        Ok(segment)
    }

    fn validate(&self, summary: &Summary, words: &[f64]) -> Result<Option<Type2>, String> {
        let (begin, end) = data_range(summary, words.len())?;
        let ordered = self.start.is_finite() && self.end.is_finite() && self.start <= self.end;
        if !ordered {
            return Err(format!(
                "coverage {} to {} s is not a finite, ordered interval",
                self.start, self.end
            ));
        }
        if self.data_type != CHEBYSHEV_POSITION {
            return Ok(None);
        }
        Type2::new(words, begin, end, [self.start, self.end]).map(Some)
    }

    pub fn body(&self) -> i32 {
        self.body
    }

    pub fn center(&self) -> i32 {
        self.center
    }

    pub fn frame(&self) -> i32 {
        self.frame
    }

    pub fn data_type(&self) -> i32 {
        self.data_type
    }

    pub fn start(&self) -> TDB {
        tdb_at(self.start)
    }

    pub fn end(&self) -> TDB {
        tdb_at(self.end)
    }

    pub(super) fn covers(&self, t: f64) -> bool {
        self.start <= t && t <= self.end
    }

    fn label(&self) -> Label {
        Label {
            index: self.index,
            body: self.body,
            center: self.center,
        }
    }

    fn chebyshev(&self) -> Result<&Type2, SpkError> {
        let type2 = self
            .type2
            .as_ref()
            .ok_or(SpkError::UnsupportedType(self.data_type))?;
        if self.frame != J2000_FRAME {
            return Err(SpkError::UnsupportedFrame(self.frame));
        }
        Ok(type2)
    }

    fn state(&self, words: &[f64], t: f64) -> Result<(Vector3, Vector3), SpkError> {
        self.chebyshev()?
            .state(words, t, &self.label())
            .map_err(SpkError::InvalidData)
    }

    fn position(&self, words: &[f64], t: f64) -> Result<Vector3, SpkError> {
        self.chebyshev()?
            .position(words, t, &self.label())
            .map_err(SpkError::InvalidData)
    }
}

impl fmt::Debug for SpkFile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SpkFile")
            .field("words", &self.words.len())
            .field("segments", &self.segments)
            .finish()
    }
}

impl fmt::Display for Label {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "segment {} (body {}, center {})",
            self.index, self.body, self.center
        )
    }
}

impl Epoch {
    fn new(tdb: &TDB) -> Result<Self, SpkError> {
        let jd = tdb.to_julian_date();
        let days = (jd.jd1() - J2000_JD) + jd.jd2();
        let seconds = days * SECONDS_PER_DAY_F64;
        let jd = J2000_JD + days;
        if seconds.is_finite() {
            return Ok(Self { seconds, jd });
        }
        Err(SpkError::InvalidEpoch { jd })
    }
}

fn tdb_at(seconds: f64) -> TDB {
    TDB::from_julian_date(JulianDate::new(J2000_JD, seconds / SECONDS_PER_DAY_F64))
}

fn data_range(summary: &Summary, file_words: usize) -> Result<(usize, usize), String> {
    match (
        usize::try_from(summary.begin),
        usize::try_from(summary.end_address),
    ) {
        (Ok(begin), Ok(end)) if 1 <= begin && begin <= end && end <= file_words => Ok((begin, end)),
        _ => Err(format!(
            "data addresses {} to {} are not an ordered range inside 1 to {}",
            summary.begin, summary.end_address, file_words
        )),
    }
}
