use super::chebyshev;
use super::daf::is_whole;
use celestial_core::matrix::Vector3;
use std::fmt::Display;

const DIRECTORY_WORDS: usize = 4;
// τ can land a few ULPs past ±1 at a record boundary; anything further means
// the record does not cover the epoch.
const TAU_LIMIT: f64 = 1.0 + 1.0 / 1_048_576.0;

#[derive(Debug, Clone)]
pub(super) struct Type2 {
    init: f64,
    intlen: f64,
    rsize: usize,
    n_records: usize,
    n_coeffs: usize,
    first: usize,
}

struct Located<'w> {
    index: usize,
    record: &'w [f64],
    tau: f64,
}

impl Type2 {
    pub(super) fn new(
        words: &[f64],
        begin: usize,
        end: usize,
        coverage: [f64; 2],
    ) -> Result<Self, String> {
        let [init, intlen, rsize, n] = directory(&words[begin - 1..end])?;
        check_timing(init, intlen)?;
        let rsize = record_size(rsize)?;
        let n_records = record_count(n, words.len())?;
        check_length(n_records, rsize, end - begin + 1)?;
        let segment = Self {
            init,
            intlen,
            rsize,
            n_records,
            n_coeffs: (rsize - 2) / 3,
            first: begin - 1,
        };
        segment.check_coverage(coverage)?;
        Ok(segment)
    }

    pub(super) fn state(
        &self,
        words: &[f64],
        t: f64,
        label: &dyn Display,
    ) -> Result<(Vector3, Vector3), String> {
        let at = self.locate(words, t, label)?;
        let position = self.position_at(&at);
        let velocity = self.velocity_at(&at);
        if is_finite(&position) && is_finite(&velocity) {
            return Ok((position, velocity));
        }
        Err(format!(
            "record {} of {} gives a non-finite state",
            at.index, label
        ))
    }

    pub(super) fn position(
        &self,
        words: &[f64],
        t: f64,
        label: &dyn Display,
    ) -> Result<Vector3, String> {
        let at = self.locate(words, t, label)?;
        let position = self.position_at(&at);
        if is_finite(&position) {
            return Ok(position);
        }
        Err(format!(
            "record {} of {} gives a non-finite position",
            at.index, label
        ))
    }

    fn locate<'w>(
        &self,
        words: &'w [f64],
        t: f64,
        label: &dyn Display,
    ) -> Result<Located<'w>, String> {
        let index = (libm::floor((t - self.init) / self.intlen) as usize).min(self.n_records - 1);
        let start = self.first + index * self.rsize;
        let record = &words[start..start + self.rsize];
        let tau = (t - record[0]) / record[1];
        let at = Located { index, record, tau };
        at.check(t, label)?;
        Ok(at)
    }

    fn series<'r>(&self, record: &'r [f64]) -> [&'r [f64]; 3] {
        let n = self.n_coeffs;
        [
            &record[2..2 + n],
            &record[2 + n..2 + 2 * n],
            &record[2 + 2 * n..2 + 3 * n],
        ]
    }

    fn position_at(&self, at: &Located) -> Vector3 {
        let [x, y, z] = self.series(at.record);
        Vector3::new(
            chebyshev::value(x, at.tau),
            chebyshev::value(y, at.tau),
            chebyshev::value(z, at.tau),
        )
    }

    fn velocity_at(&self, at: &Located) -> Vector3 {
        let radius = at.record[1];
        let [x, y, z] = self.series(at.record);
        Vector3::new(
            chebyshev::derivative(x, at.tau, radius),
            chebyshev::derivative(y, at.tau, radius),
            chebyshev::derivative(z, at.tau, radius),
        )
    }

    fn check_coverage(&self, [start, end]: [f64; 2]) -> Result<(), String> {
        let last = self.init + self.n_records as f64 * self.intlen;
        if self.init <= start && end <= last {
            return Ok(());
        }
        Err(format!(
            "records from {} s for {} x {} s do not cover the segment's {} to {} s",
            self.init, self.n_records, self.intlen, start, end
        ))
    }
}

impl Located<'_> {
    fn check(&self, t: f64, label: &dyn Display) -> Result<(), String> {
        let radius = self.record[1];
        let usable = radius > 0.0 && radius.is_finite();
        if !usable {
            return Err(format!(
                "record {} of {} has radius {} s",
                self.index, label, radius
            ));
        }
        let distance = libm::fabs(self.tau);
        if distance <= TAU_LIMIT {
            return Ok(());
        }
        Err(format!(
            "TDB {} s is {} radii from the midpoint of record {} of {}",
            t, distance, self.index, label
        ))
    }
}

fn directory(segment: &[f64]) -> Result<[f64; DIRECTORY_WORDS], String> {
    match segment.last_chunk() {
        Some(&directory) => Ok(directory),
        None => Err(format!(
            "{} words cannot hold a type 2 directory",
            segment.len()
        )),
    }
}

fn check_timing(init: f64, intlen: f64) -> Result<(), String> {
    if !init.is_finite() {
        return Err(format!("INIT {} is not finite", init));
    }
    let usable = intlen > 0.0 && intlen.is_finite();
    if !usable {
        return Err(format!("INTLEN {} is not a positive finite number", intlen));
    }
    Ok(())
}

fn record_size(rsize: f64) -> Result<usize, String> {
    let size = rsize as usize;
    if is_whole(rsize) && rsize >= 5.0 && (size - 2).is_multiple_of(3) {
        return Ok(size);
    }
    Err(format!("RSIZE {} is not 2 + 3n for a whole n >= 1", rsize))
}

fn record_count(n: f64, file_words: usize) -> Result<usize, String> {
    if is_whole(n) && n >= 1.0 && n <= file_words as f64 {
        return Ok(n as usize);
    }
    Err(format!(
        "N {} is not a whole number from 1 to {}",
        n, file_words
    ))
}

fn check_length(n_records: usize, rsize: usize, length: usize) -> Result<(), String> {
    let words = n_records
        .checked_mul(rsize)
        .and_then(|w| w.checked_add(DIRECTORY_WORDS));
    if words == Some(length) {
        return Ok(());
    }
    Err(format!(
        "{} records of {} words plus 4 is not the segment's {} words",
        n_records, rsize, length
    ))
}

fn is_finite(v: &Vector3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}
