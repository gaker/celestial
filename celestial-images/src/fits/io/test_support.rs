use crate::fits::compression::CompressionAlgorithm;
use crate::fits::data::array::DataArray;
use crate::fits::errors::Result;
use crate::fits::header::keywords::Keyword;
use crate::fits::io::writer::{compressed_image_bytes, primary_image_bytes};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

// The file is created at the first write, after every header card has been
// formatted, so a value that cannot be written leaves an existing file intact.
pub(crate) struct FitsWriter {
    path: PathBuf,
    writer: Option<BufWriter<File>>,
}

impl FitsWriter {
    pub(crate) fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        Ok(Self {
            path: path.as_ref().to_path_buf(),
            writer: None,
        })
    }

    fn output(&mut self) -> Result<&mut BufWriter<File>> {
        let writer = match self.writer.take() {
            Some(writer) => writer,
            None => BufWriter::new(File::create(&self.path)?),
        };
        Ok(self.writer.insert(writer))
    }

    fn emit(&mut self, bytes: &[u8]) -> Result<()> {
        let output = self.output()?;
        output.write_all(bytes)?;
        output.flush()?;
        Ok(())
    }

    pub(crate) fn write_primary_image<T>(
        &mut self,
        data: &[T],
        dimensions: &[usize],
        keywords: &[Keyword],
    ) -> Result<()>
    where
        T: DataArray,
    {
        self.emit(&primary_image_bytes(data, dimensions, keywords, false)?)
    }

    pub(crate) fn write_primary_image_with_checksum<T>(
        &mut self,
        data: &[T],
        dimensions: &[usize],
        keywords: &[Keyword],
    ) -> Result<()>
    where
        T: DataArray,
    {
        self.emit(&primary_image_bytes(data, dimensions, keywords, true)?)
    }

    // Every tile is compressed and both headers are formatted before the file
    // is created, so an error leaves no file behind.
    pub(crate) fn write_compressed_image<T>(
        &mut self,
        data: &[T],
        dimensions: &[usize],
        tile_size: Option<(usize, usize)>,
        algorithm: CompressionAlgorithm,
        keywords: &[Keyword],
    ) -> Result<()>
    where
        T: DataArray,
    {
        let bytes = compressed_image_bytes(data, dimensions, tile_size, algorithm, keywords)?;
        self.emit(&bytes)
    }
}
