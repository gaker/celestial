use super::quantize::Method;
use super::spec::{Quantization, Spec};
use super::table::Table;
use crate::core::types::BitPix;
use crate::fits::errors::{FitsError, Result};
use crate::fits::header::parser::Header;

// The columns of FITS 4.0 §10.1.3 and §10.2 that this reader uses.
pub(super) struct Columns {
    pub(super) compressed: usize,
    pub(super) gzip: Option<usize>,
    pub(super) raw: Option<usize>,
    pub(super) quantized: Option<Quantized>,
    pub(super) blank: Option<usize>,
    pub(super) mask: Option<usize>,
}

// The ZSCALE and ZZERO columns of quantized floats, and the ZQUANTIZ method.
#[derive(Clone, Copy)]
pub(super) struct Quantized {
    pub(super) scale: usize,
    pub(super) zero: usize,
    pub(super) method: Method,
}

impl Columns {
    pub(super) fn find(header: &Header, table: &Table, spec: &Spec) -> Result<Self> {
        let array = |name: &str| match table.column(name) {
            Some(c) if !table.is_array(c) => Err(invalid(format!("{name} is not an array column"))),
            found => Ok(found),
        };
        let compressed = array("COMPRESSED_DATA")?
            .ok_or_else(|| invalid("the table has no COMPRESSED_DATA column".to_string()))?;
        Ok(Self {
            compressed,
            gzip: array("GZIP_COMPRESSED_DATA")?,
            raw: array("UNCOMPRESSED_DATA")?,
            quantized: quantized(header, table, spec)?,
            blank: table.column("ZBLANK"),
            mask: array("NULL_PIXEL_MASK")?,
        })
    }

    // FITS 4.0 §10.1.3: BLANK gives the undefined value of an integer image.
    // Without a BLANK card, ZBLANK gives it, from the column when there is
    // one. The restored image holds one BLANK, so every tile must agree.
    pub(super) fn integer_blank(
        &self,
        header: &Header,
        table: &Table,
        spec: &Spec,
    ) -> Result<Option<i64>> {
        let float = matches!(spec.bitpix, BitPix::F32 | BitPix::F64);
        if float || header.get_keyword("BLANK").is_some() {
            return Ok(None);
        }
        let Some(column) = self.blank else {
            return Ok(spec.blank);
        };
        let mut values = (0..table.rows()).map(|row| table.integer(row, column));
        let first = values.next().transpose()?;
        for value in values {
            if Some(value?) != first {
                return Err(FitsError::UnsupportedCompression(
                    "integer tiles whose ZBLANK differs between tiles".to_string(),
                ));
            }
        }
        Ok(first)
    }
}

// FITS 4.0 §10.2: ZSCALE and ZZERO columns mark quantized floats, and their
// absence marks lossless floats. The keyword form of ZSCALE is not standard.
fn quantized(header: &Header, table: &Table, spec: &Spec) -> Result<Option<Quantized>> {
    let Some((scale, zero)) = scale_columns(header, table)? else {
        return Ok(None);
    };
    if !matches!(spec.bitpix, BitPix::F32 | BitPix::F64) {
        return Err(FitsError::UnsupportedCompression(format!(
            "quantized tiles of BITPIX {}",
            spec.bitpix.value()
        )));
    }
    match spec.quantization {
        Quantization::Quantized(method) => Ok(Some(Quantized {
            scale,
            zero,
            method,
        })),
        Quantization::Lossless => Err(invalid("ZQUANTIZ = 'NONE' with ZSCALE columns".to_string())),
    }
}

// Returns the ZSCALE and ZZERO column indexes, or `None` when the table has
// neither.
fn scale_columns(header: &Header, table: &Table) -> Result<Option<(usize, usize)>> {
    match (table.column("ZSCALE"), table.column("ZZERO")) {
        (Some(scale), Some(zero)) => Ok(Some((scale, zero))),
        (None, None)
            if ["ZSCALE", "ZZERO"]
                .iter()
                .any(|k| header.get_keyword(k).is_some()) =>
        {
            Err(FitsError::UnsupportedCompression(
                "ZSCALE or ZZERO given as a keyword".to_string(),
            ))
        }
        (None, None) => Ok(None),
        _ => Err(invalid(
            "ZSCALE and ZZERO columns must come as a pair".to_string(),
        )),
    }
}

fn invalid(message: String) -> FitsError {
    FitsError::InvalidFormat(message)
}
