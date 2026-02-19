//! Gaia DR3 ECSV parser
//!
//! Handles gzipped ECSV format with `#` metadata lines.

use std::collections::HashMap;
use std::io::BufRead;

pub struct GaiaStarRaw {
    pub source_id: i64,
    pub ra: f64,
    pub dec: f64,
    pub pmra: f64,
    pub pmdec: f64,
    pub parallax: f64,
    pub mag: f32,
    pub flags: u16,
}

pub const FLAG_HAS_PROPER_MOTION: u16 = 1 << 0;
pub const FLAG_HAS_PARALLAX: u16 = 1 << 1;
pub const FLAG_RUWE_SUSPECT: u16 = 1 << 2;
pub const FLAG_NO_5PARAM: u16 = 1 << 3;
pub const FLAG_BP_RP_EXCESS_SUSPECT: u16 = 1 << 4;
pub const FLAG_SOURCE_HIPPARCOS: u16 = 1 << 5;

struct ColumnIndices {
    source_id: usize,
    ra: usize,
    dec: usize,
    pmra: usize,
    pmdec: usize,
    parallax: usize,
    phot_g_mean_mag: usize,
    ruwe: usize,
    astrometric_params_solved: usize,
    duplicated_source: usize,
    phot_bp_rp_excess_factor: usize,
}

pub struct GaiaParser<R: BufRead> {
    reader: R,
    indices: ColumnIndices,
    mag_limit: f32,
    line_buf: String,
}

impl<R: BufRead> GaiaParser<R> {
    pub fn new(mut reader: R, mag_limit: f32) -> anyhow::Result<Self> {
        let indices = Self::parse_header(&mut reader)?;
        Ok(Self {
            reader,
            indices,
            mag_limit,
            line_buf: String::with_capacity(4096),
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
        Self::build_column_indices(&line)
    }

    fn build_column_indices(header_line: &str) -> anyhow::Result<ColumnIndices> {
        let mut col_map: HashMap<&str, usize> = HashMap::new();
        for (idx, col) in header_line.trim().split(',').enumerate() {
            col_map.insert(col, idx);
        }
        Ok(ColumnIndices {
            source_id: Self::require_column(&col_map, "source_id")?,
            ra: Self::require_column(&col_map, "ra")?,
            dec: Self::require_column(&col_map, "dec")?,
            pmra: Self::require_column(&col_map, "pmra")?,
            pmdec: Self::require_column(&col_map, "pmdec")?,
            parallax: Self::require_column(&col_map, "parallax")?,
            phot_g_mean_mag: Self::require_column(&col_map, "phot_g_mean_mag")?,
            ruwe: Self::require_column(&col_map, "ruwe")?,
            astrometric_params_solved: Self::require_column(&col_map, "astrometric_params_solved")?,
            duplicated_source: Self::require_column(&col_map, "duplicated_source")?,
            phot_bp_rp_excess_factor: Self::require_column(&col_map, "phot_bp_rp_excess_factor")?,
        })
    }

    fn require_column(col_map: &HashMap<&str, usize>, name: &str) -> anyhow::Result<usize> {
        col_map
            .get(name)
            .copied()
            .ok_or_else(|| anyhow::anyhow!("Missing column: {}", name))
    }
}

impl<R: BufRead> Iterator for GaiaParser<R> {
    type Item = anyhow::Result<GaiaStarRaw>;

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
            match self.parse_row() {
                Ok(Some(star)) => return Some(Ok(star)),
                Ok(None) => continue,
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

impl<R: BufRead> GaiaParser<R> {
    fn parse_row(&self) -> anyhow::Result<Option<GaiaStarRaw>> {
        let fields: Vec<&str> = self.line_buf.trim().split(',').collect();
        if self.should_skip_row(&fields) {
            return Ok(None);
        }
        let mag = parse_f32(fields.get(self.indices.phot_g_mean_mag).copied());
        if mag.is_none() || mag.unwrap() > self.mag_limit {
            return Ok(None);
        }
        Ok(Some(self.build_star(&fields, mag.unwrap())))
    }

    fn should_skip_row(&self, fields: &[&str]) -> bool {
        let dup = fields
            .get(self.indices.duplicated_source)
            .copied()
            .unwrap_or("");
        dup.eq_ignore_ascii_case("true")
    }

    fn build_star(&self, fields: &[&str], mag: f32) -> GaiaStarRaw {
        let pmra = parse_f64(fields.get(self.indices.pmra).copied());
        let pmdec = parse_f64(fields.get(self.indices.pmdec).copied());
        let parallax = parse_f64(fields.get(self.indices.parallax).copied());
        let flags = self.compute_flags(fields, pmra, pmdec, parallax);
        GaiaStarRaw {
            source_id: parse_i64(fields.get(self.indices.source_id).copied()).unwrap_or(0),
            ra: parse_f64(fields.get(self.indices.ra).copied()).unwrap_or(0.0),
            dec: parse_f64(fields.get(self.indices.dec).copied()).unwrap_or(0.0),
            pmra: pmra.unwrap_or(0.0),
            pmdec: pmdec.unwrap_or(0.0),
            parallax: parallax.unwrap_or(0.0),
            mag,
            flags,
        }
    }

    fn compute_flags(
        &self,
        fields: &[&str],
        pmra: Option<f64>,
        pmdec: Option<f64>,
        parallax: Option<f64>,
    ) -> u16 {
        let mut flags = 0u16;
        if pmra.is_some() && pmdec.is_some() {
            flags |= FLAG_HAS_PROPER_MOTION;
        }
        if parallax.is_some() {
            flags |= FLAG_HAS_PARALLAX;
        }
        flags |= self.compute_quality_flags(fields);
        flags
    }

    fn compute_quality_flags(&self, fields: &[&str]) -> u16 {
        let mut flags = 0u16;
        if let Some(ruwe) = parse_f32(fields.get(self.indices.ruwe).copied()) {
            if ruwe > 1.4 {
                flags |= FLAG_RUWE_SUSPECT;
            }
        }
        if let Some(params) = parse_i8(fields.get(self.indices.astrometric_params_solved).copied())
        {
            if params != 31 {
                flags |= FLAG_NO_5PARAM;
            }
        }
        if let Some(excess) = parse_f32(fields.get(self.indices.phot_bp_rp_excess_factor).copied())
        {
            if excess > 3.0 {
                flags |= FLAG_BP_RP_EXCESS_SUSPECT;
            }
        }
        flags
    }
}

fn parse_i64(s: Option<&str>) -> Option<i64> {
    s.and_then(|v| if v.is_empty() { None } else { v.parse().ok() })
}

fn parse_i8(s: Option<&str>) -> Option<i8> {
    s.and_then(|v| if v.is_empty() { None } else { v.parse().ok() })
}

fn parse_f64(s: Option<&str>) -> Option<f64> {
    s.and_then(|v| if v.is_empty() { None } else { v.parse().ok() })
}

fn parse_f32(s: Option<&str>) -> Option<f32> {
    s.and_then(|v| if v.is_empty() { None } else { v.parse().ok() })
}
