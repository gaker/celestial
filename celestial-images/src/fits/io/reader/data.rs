use super::{data_size, FitsFile};
use crate::core::types::BitPix;
use crate::fits::data::array::DataArray;
use crate::fits::errors::{FitsError, Result};
use crate::fits::hdu::primary::PrimaryHdu;
use crate::fits::hdu::{read_array, Hdu, HduTrait};
use std::io::{Read, Seek, SeekFrom};

impl<R: Read + Seek> FitsFile<R> {
    pub fn get_image_info(&mut self, index: usize) -> Result<(Vec<usize>, BitPix)> {
        let hdu = self.read_hdu(index)?;
        let (dimensions, bitpix) = match hdu {
            Hdu::Primary(primary) => (primary.data_dimensions()?, primary.bitpix()),
            Hdu::Image(image) => (image.data_dimensions()?, image.bitpix()),
            _ => {
                return Err(FitsError::InvalidFormat(
                    "HDU type does not support image data".to_string(),
                ))
            }
        };

        let bitpix = bitpix.ok_or_else(|| FitsError::KeywordNotFound {
            keyword: "BITPIX".to_string(),
        })?;

        Ok((dimensions, bitpix))
    }

    pub fn primary_hdu_with_data<T>(&mut self) -> Result<(PrimaryHdu, Vec<T>)>
    where
        T: DataArray,
    {
        let primary = self.primary_hdu()?;
        match primary {
            Hdu::Primary(hdu) => {
                let data = hdu.read_data(&mut self.reader)?;
                Ok((hdu, data))
            }
            _ => Err(FitsError::InvalidFormat(
                "the primary HDU holds random groups, not an image".to_string(),
            )),
        }
    }

    // The data unit without its block padding. The scan has checked that the
    // file holds this many bytes after the header.
    pub(crate) fn read_data_unit(&mut self, index: usize) -> Result<Vec<u8>> {
        let header = self.get_header(index)?;
        let start = self
            .hdus
            .get(index)
            .ok_or(FitsError::HduNotFound(index))?
            .data_start;
        let len = data_size::data_unit_len(&header)?;
        let mut data = Vec::new();
        data.try_reserve_exact(len)
            .map_err(|_| FitsError::InvalidFormat(format!("cannot allocate {len} data bytes")))?;
        data.resize(len, 0);
        self.reader.seek(SeekFrom::Start(start))?;
        self.reader.read_exact(&mut data)?;
        Ok(data)
    }

    // The image array of a primary or IMAGE HDU, as stored.
    pub(crate) fn read_image_array(&mut self, index: usize) -> Result<Vec<u8>> {
        match self.read_hdu(index)? {
            Hdu::Primary(hdu) => read_array(&hdu, &mut self.reader),
            Hdu::Image(hdu) => read_array(hdu.as_ref(), &mut self.reader),
            _ => Err(FitsError::InvalidFormat(
                "HDU type does not support image data reading".to_string(),
            )),
        }
    }

    pub fn read_hdu_with_data<T>(&mut self, index: usize) -> Result<(Hdu, Vec<T>)>
    where
        T: DataArray,
    {
        let hdu = self.read_hdu(index)?;
        let data = match &hdu {
            Hdu::Primary(primary) => primary.read_data(&mut self.reader)?,
            Hdu::Image(image) => image.read_data(&mut self.reader)?,
            _ => {
                return Err(FitsError::InvalidFormat(
                    "HDU type does not support image data reading".to_string(),
                ))
            }
        };
        Ok((hdu, data))
    }
}
