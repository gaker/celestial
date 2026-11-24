// The setters and `run` of `Solver`.

use celestial_coords::frames::icrs::ICRSPosition;
use celestial_time::julian::JulianDate;

use crate::error::SolveError;
use crate::metadata::image_metadata::{metadata_from_image, ImageMetadata};

use super::pipeline::run_pipeline;
use super::{SolveParams, SolveResult, Solver};

impl<'a> Solver<'a> {
    /// Overrides the position hint. Bypasses `CRVAL1/2`, `RA/DEC`, and `OBJCTRA/DEC`
    /// header lookups.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_coords::frames::icrs::ICRSPosition;
    ///
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog)
    ///     .hint(ICRSPosition::from_degrees(83.633, 22.014)?)
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn hint(mut self, pos: ICRSPosition) -> Self {
        self.overrides.hint = Some(pos);
        self
    }

    /// Overrides the observation epoch. Bypasses `DATE-OBS` and `MJD-OBS` header lookups.
    ///
    /// The epoch is used for catalog proper-motion propagation; nearby stars have moved
    /// significantly since J2000.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_time::scales::utc::utc_from_calendar;
    ///
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let jd = utc_from_calendar(2026, 4, 21, 20, 0, 0.0)?.to_julian_date();
    /// let result = celestial_solver::solve::solve(&img, &catalog).epoch(jd).run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn epoch(mut self, jd: JulianDate) -> Self {
        self.overrides.epoch = Some(jd);
        self
    }

    /// Overrides telescope focal length in mm. Bypasses the `FOCALLEN` header.
    ///
    /// Combined with [`Solver::pixel_size_um`] to compute the plate scale. Use this
    /// when the header reflects the bare OTA but the rig has a focal reducer, flattener,
    /// or Barlow inline.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog)
    ///     .focal_length_mm(500)    // integer literal auto-converts to f64
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn focal_length_mm<T: Into<f64>>(mut self, mm: T) -> Self {
        self.overrides.focal_mm = Some(mm.into());
        self
    }

    /// Overrides pixel size in microns. Bypasses the `XPIXSZ` header.
    ///
    /// Combined with [`Solver::focal_length_mm`] to compute the plate scale. Use this
    /// when the header reports the unbinned pixel size but the frame was binned, or
    /// when the sensor's spec differs from the header.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog)
    ///     .pixel_size_um(3.76)
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn pixel_size_um<T: Into<f64>>(mut self, um: T) -> Self {
        self.overrides.pixel_um = Some(um.into());
        self
    }

    /// Overrides the plate scale directly in arcseconds per pixel.
    ///
    /// Takes priority over [`Solver::focal_length_mm`] / [`Solver::pixel_size_um`] and
    /// skips the `206.2648 * pixel_um / focal_mm` calculation. Use this when you know
    /// the plate scale from a prior solve and don't have the optics parameters.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog)
    ///     .plate_scale_arcsec(1.55)
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn plate_scale_arcsec<T: Into<f64>>(mut self, arcsec: T) -> Self {
        self.overrides.scale_arcsec = Some(arcsec.into());
        self
    }

    /// Replaces the full metadata. Skips header parsing entirely.
    ///
    /// The escape hatch for images without FITS headers: PNG, raw TIFF, anything where
    /// the position / epoch / scale live outside the file. When set, the individual
    /// override setters are ignored.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_solver::metadata::image_metadata::ImageMetadata;
    /// use celestial_coords::frames::icrs::ICRSPosition;
    /// use celestial_time::scales::utc::utc_from_calendar;
    ///
    /// # let img = celestial_images::formats::image::Image::open("frame.png")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let meta = ImageMetadata {
    ///     hint: ICRSPosition::from_degrees(83.633, 22.014)?,
    ///     scale_arcsec: 1.5,
    ///     epoch: utc_from_calendar(2026, 4, 21, 20, 0, 0.0)?.to_julian_date(),
    ///     focal_mm: None,
    ///     pixel_um: None,
    /// };
    /// let result = celestial_solver::solve::solve(&img, &catalog)
    ///     .with_metadata(meta)
    ///     .run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn with_metadata(mut self, meta: ImageMetadata) -> Self {
        self.metadata = Some(meta);
        self
    }

    /// Replaces the solver tuning parameters.
    ///
    /// See [`SolveParams`] for the individual knobs.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use celestial_solver::solve::SolveParams;
    ///
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let params = SolveParams { sip_order: None, ..SolveParams::default() };
    /// let result = celestial_solver::solve::solve(&img, &catalog).params(params).run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn params(mut self, params: SolveParams) -> Self {
        self.params = params;
        self
    }

    /// Runs the solve pipeline and returns the [`SolveResult`].
    ///
    /// # Errors
    ///
    /// Returns [`SolveError::Metadata`] when header parsing fails (missing hint, missing
    /// plate scale) or the metadata is not finite, [`SolveError::InvalidParameter`] for a
    /// parameter out of its range, [`SolveError::NoStars`] when no stars are detected,
    /// [`SolveError::NoMatch`] when no quad match passes verification, the initial WCS
    /// disagrees with the hint, or neither the refined nor the initial solution has
    /// log-odds above ln(1e9) and at least [`SolveParams::min_refined_pairs`] fitted
    /// pairs, [`SolveError::FitFailed`] when the initial WCS fit fails,
    /// [`SolveError::Catalog`] when a catalog query fails, and
    /// [`SolveError::Verification`] when the returned solution cannot be verified, such
    /// as for a catalog error during its verification.
    /// A failed refinement keeps the initial WCS, and so does a refined WCS that fails
    /// verification or does not build; a SIP distortion whose inverse cannot be fitted is
    /// left out. [`SolveResult::refinement`] records each of these.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # let img = celestial_images::formats::image::Image::open("f.fits")?;
    /// # let catalog = celestial_catalog::query::catalog::Catalog::open("cat.bin")?;
    /// let result = celestial_solver::solve::solve(&img, &catalog).run()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn run(self) -> Result<SolveResult, SolveError> {
        let meta = match self.metadata {
            Some(m) => m,
            None => metadata_from_image(self.image, &self.overrides)?,
        };
        meta.validate()?;
        run_pipeline(self.image, self.catalog, &meta, &self.params)
    }
}

#[cfg(test)]
mod tests {
    use celestial_coords::frames::icrs::ICRSPosition;
    use celestial_time::julian::JulianDate;

    use crate::solve::test_support::{build_catalog_keep_alive, meta_at, zero_image};
    use crate::solve::{solve, SolveParams};

    #[test]
    fn solver_builder_chain_compiles_and_runs_to_no_stars() {
        // End-to-end exercise of every builder setter: hint, epoch, focal_length_mm,
        // pixel_size_um, plate_scale_arcsec, params. The image is zeroed so the
        // detector returns no stars and we exit at the `ensure!(!stars.is_empty())`.
        let (cat, _file) = build_catalog_keep_alive();
        let img = zero_image(1024, 1024);

        let result = solve(&img, &cat)
            .hint(ICRSPosition::from_degrees(10.0, 20.0).unwrap())
            .epoch(JulianDate::new(2451545.0, 0.0))
            .focal_length_mm(500)
            .pixel_size_um(3.76_f64)
            .plate_scale_arcsec(1.55_f64)
            .params(SolveParams::default())
            .run();
        let err = result.err().expect("zero image must error");
        assert!(
            err.to_string().contains("no stars detected"),
            "expected 'no stars detected', got: {err}"
        );
    }

    #[test]
    fn solver_with_metadata_skips_header_parse_and_reaches_detector() {
        // with_metadata is the escape hatch: if it works, the missing-headers
        // image still reaches the detection stage instead of erroring early.
        let (cat, _file) = build_catalog_keep_alive();
        let img = zero_image(1024, 1024);
        let meta = meta_at(0.0, 0.0, 1.5);

        let err = solve(&img, &cat)
            .with_metadata(meta)
            .run()
            .err()
            .expect("zero image must error");
        assert!(err.to_string().contains("no stars detected"));
    }

    #[test]
    fn run_without_metadata_or_headers_returns_metadata_error() {
        // No overrides, no FITS headers → metadata_from_image fails before
        // we ever reach the detector.
        let (cat, _file) = build_catalog_keep_alive();
        let img = zero_image(1024, 1024);

        let err = solve(&img, &cat)
            .run()
            .err()
            .expect("no headers, no overrides must error");
        // The exact wording comes from MetadataError; just confirm it isn't
        // the detection-stage error, so the solve failed earlier.
        assert!(!err.to_string().contains("no stars detected"));
    }
}
