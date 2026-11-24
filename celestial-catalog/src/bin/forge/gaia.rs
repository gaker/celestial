//! Gaia DR3 ECSV parser
//!
//! Handles gzipped ECSV format with `#` metadata lines.

use crate::bp_rp_excess::is_excess_suspect;
use crate::intermediate;
use celestial_catalog::query::catalog::{
    StarRecord, FLAG_5_OR_6_PARAM, FLAG_BP_RP_EXCESS_SUSPECT, FLAG_DUPLICATED_SOURCE,
    FLAG_HAS_PARALLAX, FLAG_HAS_PROPER_MOTION, FLAG_NO_5PARAM, FLAG_RUWE_SUSPECT,
};
use std::collections::HashMap;
use std::io::BufRead;
use std::str::FromStr;

// Gaia DR3 positions are given at epoch J2016.0.
pub(crate) const GAIA_DR3_EPOCH: f64 = 2016.0;

// No star is brighter than Sirius, V = -1.46, so a smaller magnitude is a corrupt value.
const MIN_MAG: f32 = -2.0;

const COLUMNS: [&str; 12] = [
    "source_id",
    "ra",
    "dec",
    "pmra",
    "pmdec",
    "parallax",
    "phot_g_mean_mag",
    "bp_rp",
    "ruwe",
    "astrometric_params_solved",
    "duplicated_source",
    "phot_bp_rp_excess_factor",
];

// Positions of the COLUMNS fields in a data row, in COLUMNS order.
struct ColumnIndices([usize; COLUMNS.len()]);

#[derive(Clone, Copy)]
enum Col {
    SourceId,
    Ra,
    Dec,
    Pmra,
    Pmdec,
    Parallax,
    GMag,
    BpRp,
    Ruwe,
    ParamsSolved,
    DuplicatedSource,
    ExcessFactor,
}

// Data rows read so far, by outcome. Rows fainter than the limit are the remainder.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(crate) struct RowCounts {
    pub(crate) read: u64,
    pub(crate) kept: u64,
    // phot_g_mean_mag missing, unparseable, not finite or below MIN_MAG.
    pub(crate) bad_magnitude: u64,
    // source_id, ra or dec missing or unparseable, or a place outside RA [0, 360) or
    // Dec [-90, 90].
    pub(crate) bad_position: u64,
}

impl RowCounts {
    pub(crate) fn add(&mut self, other: RowCounts) {
        self.read += other.read;
        self.kept += other.kept;
        self.bad_magnitude += other.bad_magnitude;
        self.bad_position += other.bad_position;
    }
}

enum Row {
    Kept(StarRecord),
    Fainter,
    BadMagnitude,
    BadPosition,
}

pub(crate) struct GaiaParser<R: BufRead> {
    reader: R,
    indices: ColumnIndices,
    mag_limit: f32,
    line_buf: String,
    counts: RowCounts,
}

impl<R: BufRead> GaiaParser<R> {
    pub(crate) fn new(mut reader: R, mag_limit: f32) -> anyhow::Result<Self> {
        let indices = Self::parse_header(&mut reader)?;
        Ok(Self {
            reader,
            indices,
            mag_limit,
            line_buf: String::with_capacity(4096),
            counts: RowCounts::default(),
        })
    }

    fn parse_header(reader: &mut R) -> anyhow::Result<ColumnIndices> {
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                anyhow::bail!("EOF before finding header");
            }
            if !line.starts_with('#') {
                break;
            }
        }
        build_column_indices(&line)
    }

    pub(crate) fn counts(&self) -> RowCounts {
        self.counts
    }

    fn count(&mut self, row: Row) -> Option<StarRecord> {
        match row {
            Row::Kept(star) => {
                self.counts.kept += 1;
                return Some(star);
            }
            Row::Fainter => {}
            Row::BadMagnitude => self.counts.bad_magnitude += 1,
            Row::BadPosition => self.counts.bad_position += 1,
        }
        None
    }
}

fn build_column_indices(header_line: &str) -> anyhow::Result<ColumnIndices> {
    let col_map: HashMap<&str, usize> = header_line
        .trim()
        .split(',')
        .enumerate()
        .map(|(idx, col)| (col, idx))
        .collect();
    let mut indices = [0; COLUMNS.len()];
    for (slot, name) in indices.iter_mut().zip(COLUMNS) {
        *slot = *col_map
            .get(name)
            .ok_or_else(|| anyhow::anyhow!("Missing column: {}", name))?;
    }
    Ok(ColumnIndices(indices))
}

impl<R: BufRead> Iterator for GaiaParser<R> {
    type Item = anyhow::Result<StarRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_buf.clear();
            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(e) => return Some(Err(e.into())),
            }
            if self.line_buf.starts_with('#') {
                continue;
            }
            self.counts.read += 1;
            let row = Fields::new(&self.line_buf, &self.indices).parse(self.mag_limit);
            if let Some(star) = self.count(row) {
                return Some(Ok(star));
            }
        }
    }
}

struct Fields<'a> {
    values: Vec<&'a str>,
    indices: &'a ColumnIndices,
}

impl<'a> Fields<'a> {
    fn new(line: &'a str, indices: &'a ColumnIndices) -> Self {
        Self {
            values: line.trim().split(',').collect(),
            indices,
        }
    }

    fn get<T: FromStr>(&self, col: Col) -> Option<T> {
        parse_field(self.values.get(self.indices.0[col as usize]).copied())
    }

    // Real DR3 files write a missing value as `null`, which does not parse. "NaN" and "inf"
    // parse as floats, so a non-finite value also counts as missing.
    fn finite(&self, col: Col) -> Option<f64> {
        self.get::<f64>(col).filter(|v| v.is_finite())
    }

    // Booleans are quoted in DR3 files: "True" or "False".
    fn is_true(&self, col: Col) -> bool {
        let value = self.values.get(self.indices.0[col as usize]).copied();
        value.is_some_and(|v| v.trim_matches('"').eq_ignore_ascii_case("true"))
    }

    fn parse(&self, mag_limit: f32) -> Row {
        let mag = self.get::<f32>(Col::GMag);
        let Some(mag) = mag.filter(|m| m.is_finite() && *m >= MIN_MAG) else {
            return Row::BadMagnitude;
        };
        if mag > mag_limit {
            return Row::Fainter;
        }
        match self.record(mag) {
            Some(star) => Row::Kept(star),
            None => Row::BadPosition,
        }
    }

    fn record(&self, mag: f32) -> Option<StarRecord> {
        let ra = self
            .finite(Col::Ra)
            .filter(|ra| (0.0..360.0).contains(ra))?;
        let dec = self
            .finite(Col::Dec)
            .filter(|d| (-90.0..=90.0).contains(d))?;
        let pm = self.finite(Col::Pmra).zip(self.finite(Col::Pmdec));
        let parallax = self.finite(Col::Parallax);
        let mut star = intermediate::zeroed_record();
        star.source_id = self.get(Col::SourceId)?;
        (star.ra, star.dec, star.mag) = (ra, dec, mag);
        (star.pmra, star.pmdec) = pm.unwrap_or((0.0, 0.0));
        star.parallax = parallax.unwrap_or(0.0);
        star.flags = self.quality_flags(mag);
        if pm.is_some() {
            star.flags |= FLAG_HAS_PROPER_MOTION;
        }
        if parallax.is_some() {
            star.flags |= FLAG_HAS_PARALLAX;
        }
        Some(star)
    }

    fn quality_flags(&self, mag: f32) -> u16 {
        let mut flags = astrometry_flags(self.get::<u8>(Col::ParamsSolved));
        if self.finite(Col::Ruwe).is_some_and(|ruwe| ruwe > 1.4) {
            flags |= FLAG_RUWE_SUSPECT;
        }
        if self.is_true(Col::DuplicatedSource) {
            flags |= FLAG_DUPLICATED_SOURCE;
        }
        let excess = self.finite(Col::ExcessFactor).zip(self.finite(Col::BpRp));
        if excess.is_some_and(|(c, bp_rp)| is_excess_suspect(c, bp_rp, f64::from(mag))) {
            flags |= FLAG_BP_RP_EXCESS_SUSPECT;
        }
        flags
    }
}

// DR3 astrometric_params_solved is 3 (position only), 31 (5-parameter solution) or 95
// (6-parameter solution, with the pseudocolour solved too).
fn astrometry_flags(params_solved: Option<u8>) -> u16 {
    match params_solved {
        Some(31) => FLAG_5_OR_6_PARAM,
        Some(95) => FLAG_5_OR_6_PARAM | FLAG_NO_5PARAM,
        Some(_) => FLAG_NO_5PARAM,
        None => 0,
    }
}

fn parse_field<T: FromStr>(value: Option<&str>) -> Option<T> {
    value.filter(|v| !v.is_empty())?.parse().ok()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) const HEADER: &str = "source_id,ra,dec,pmra,pmdec,parallax,phot_g_mean_mag,bp_rp,\
ruwe,astrometric_params_solved,duplicated_source,phot_bp_rp_excess_factor";

    // A DR3 row with the given ID, place, G and astrometric_params_solved; booleans quoted as
    // in DR3 files, and a 2-parameter solution with null motion and parallax.
    pub(crate) fn row(id: &str, ra: &str, dec: &str, g: &str, params: u8) -> String {
        let motion = if params == 3 {
            "null,null,null"
        } else {
            "1.5,-2.5,3.0"
        };
        format!("{id},{ra},{dec},{motion},{g},1.0,1.0,{params},\"False\",1.216844")
    }

    pub(crate) fn csv(rows: &[String]) -> String {
        format!("# %ECSV 1.0\n{HEADER}\n{}\n", rows.join("\n"))
    }

    fn parse(rows: &[String]) -> (Vec<StarRecord>, RowCounts) {
        let text = csv(rows);
        let mut parser = GaiaParser::new(std::io::Cursor::new(text), 18.5).unwrap();
        let stars = parser.by_ref().map(|s| s.unwrap()).collect();
        (stars, parser.counts())
    }

    fn with_field(row: String, index: usize, value: &str) -> String {
        let mut fields: Vec<&str> = row.split(',').collect();
        fields[index] = value;
        fields.join(",")
    }

    #[test]
    fn test_duplicated_source_rows_are_kept_and_flagged() {
        let base = || row("1", "10.0", "20.0", "12.0", 31);
        let rows = [
            with_field(base(), 10, "\"True\""),
            with_field(base(), 10, "True"),
            base(),
        ];
        let (stars, counts) = parse(&rows);
        let dup: Vec<bool> = stars
            .iter()
            .map(|s| s.flags & FLAG_DUPLICATED_SOURCE != 0)
            .collect();
        assert_eq!(dup, [true, true, false]);
        assert_eq!((counts.read, counts.kept), (3, 3));
    }

    #[test]
    fn test_rows_without_valid_id_place_or_magnitude_are_skipped_and_counted() {
        let rows = [
            row("null", "10.0", "20.0", "12.0", 31),
            row("x7", "10.0", "20.0", "12.0", 31),
            row("2", "null", "20.0", "12.0", 31),
            row("3", "10.0", "NaN", "12.0", 31),
            row("4", "360.0", "20.0", "12.0", 31),
            row("5", "10.0", "-90.5", "12.0", 31),
            row("6", "10.0", "20.0", "null", 31),
            row("7", "10.0", "20.0", "-inf", 31),
            row("8", "10.0", "20.0", "-nan", 31),
            row("9", "10.0", "20.0", "-3.0", 31),
            row("10", "10.0", "20.0", "19.0", 31),
            row("11", "359.5", "-90.0", "18.5", 31),
        ];
        let (stars, counts) = parse(&rows);
        let kept: Vec<(i64, f64, f64)> = stars.iter().map(|s| (s.source_id, s.ra, s.dec)).collect();
        assert_eq!(kept, [(11, 359.5, -90.0)]);
        let expected = RowCounts {
            read: 12,
            kept: 1,
            bad_magnitude: 4,
            bad_position: 6,
        };
        assert_eq!(counts, expected);
    }

    #[test]
    fn test_non_finite_motion_counts_as_missing() {
        let rows = [with_field(row("1", "10.0", "20.0", "12.0", 31), 3, "NaN")];
        let (stars, _) = parse(&rows);
        assert_eq!((stars[0].pmra, stars[0].pmdec), (0.0, 0.0));
        assert_eq!(stars[0].flags, FLAG_HAS_PARALLAX | FLAG_5_OR_6_PARAM);
    }

    #[test]
    fn test_astrometry_flags_by_solution() {
        let rows = [31, 95, 3].map(|params| row("1", "10.0", "20.0", "12.0", params));
        let (stars, _) = parse(&rows);
        let flags: Vec<u16> = stars.iter().map(|s| s.flags).collect();
        let full = FLAG_HAS_PROPER_MOTION | FLAG_HAS_PARALLAX | FLAG_5_OR_6_PARAM;
        assert_eq!(flags, [full, full | FLAG_NO_5PARAM, FLAG_NO_5PARAM]);
    }

    // At BP-RP = 1 and G = 10, f = 1.216844 and the 3-sigma limit on |C*| is 0.019068.
    #[test]
    fn test_excess_flag_uses_corrected_excess() {
        let excess = |c: &str, bp_rp: &str| {
            let star = with_field(row("1", "10.0", "20.0", "10.0", 31), 11, c);
            with_field(star, 7, bp_rp)
        };
        let rows = [
            excess("1.216844", "1.0"),
            excess("1.25", "1.0"),
            excess("1.19", "1.0"),
            excess("3.5", "null"),
            excess("2.0", "6.0"),
        ];
        let (stars, _) = parse(&rows);
        let flagged: Vec<bool> = stars
            .iter()
            .map(|s| s.flags & FLAG_BP_RP_EXCESS_SUSPECT != 0)
            .collect();
        assert_eq!(flagged, [false, true, true, false, true]);
    }
}
