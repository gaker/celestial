// These tests compile as a separate crate, so they reach the library only
// through its public paths.

use celestial_images::core::types::ByteOrder;
use celestial_images::fits::data::array::DataArray;
use celestial_images::fits::io::reader::FitsFile;
use celestial_images::formats::astro_image::AstroImage;
use celestial_images::xisf::header::color_space::ColorSpace;
use celestial_images::xisf::header::image_info::ImageInfo;
use celestial_images::xisf::header::pixel_storage::PixelStorage;
use celestial_images::xisf::header::sample_format::SampleFormat;
use celestial_images::xisf::reader::XisfFile;
use celestial_images::xisf::writer::data_type::XisfDataType;
use celestial_images::xisf::writer::XisfWriter;
use std::io::Cursor;

fn two_image_xisf() -> Vec<u8> {
    let mut writer = XisfWriter::new(Cursor::new(Vec::new()));
    writer.add_image(&[1u8, 2, 3, 4], 2, 2, 1).unwrap();
    writer
        .add_image(&[10u16, 20, 30, 40, 50, 60], 3, 2, 1)
        .unwrap();
    writer.write_to_vec().unwrap()
}

#[test]
fn second_xisf_image_is_described_and_decoded() {
    let mut file = XisfFile::new(Cursor::new(two_image_xisf())).unwrap();
    let info: &ImageInfo = file.image_info(1).unwrap();
    assert_eq!(info.geometry, vec![3, 2]);
    assert_eq!(info.sample_format, SampleFormat::UInt16);
    assert_eq!(info.pixel_storage, PixelStorage::Planar);
    assert!(matches!(info.color_space, ColorSpace::Gray));
    assert_eq!(info.byte_order, ByteOrder::LittleEndian);
    let width = info.sample_format.bytes_per_sample();

    let raw = file.read_image_data_raw(1).unwrap();
    let samples: Vec<u16> = raw
        .chunks_exact(width)
        .map(|b| u16::from_le_bytes([b[0], b[1]]))
        .collect();
    assert_eq!(samples, vec![10, 20, 30, 40, 50, 60]);
    assert_eq!(file.read_image_data_typed::<u16>(1).unwrap(), samples);
}

fn xisf_of<T: XisfDataType>(samples: &[T], width: usize, height: usize) -> Vec<u8> {
    let mut writer = XisfWriter::new(Cursor::new(Vec::new()));
    writer.add_image(samples, width, height, 1).unwrap();
    writer.write_to_vec().unwrap()
}

#[test]
fn xisf_data_type_is_a_nameable_bound() {
    let bytes = xisf_of(&[0.25f32, 0.5], 2, 1);
    let mut file = XisfFile::new(Cursor::new(bytes)).unwrap();
    let info = file.image_info(0).unwrap();
    assert_eq!(info.sample_format, SampleFormat::Float32);
    assert_eq!(info.bounds, (0.25, 0.5));
    assert_eq!(
        file.read_image_data_typed::<f32>(0).unwrap(),
        vec![0.25, 0.5]
    );
}

fn primary_samples<T: DataArray>(path: &std::path::Path) -> Vec<T> {
    let mut fits = FitsFile::open(path).unwrap();
    fits.primary_hdu_with_data::<T>().unwrap().1
}

fn write_fits<T: DataArray + Clone>(samples: &[T], path: &std::path::Path) {
    AstroImage::new(samples, [2, 1]).write_fits(path).unwrap();
}

#[test]
fn data_array_is_a_nameable_bound() {
    let file = tempfile::Builder::new().suffix(".fits").tempfile().unwrap();
    write_fits(&[-3i16, 7], file.path());
    assert_eq!(primary_samples::<i16>(file.path()), vec![-3, 7]);
}
