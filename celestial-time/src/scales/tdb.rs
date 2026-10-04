//! Barycentric Dynamical Time (TDB) scale.
//!
//! TDB is the independent time argument for barycentric ephemerides of the solar system.
//! It differs from TT by small periodic terms (max ~1.7 ms) due to relativistic effects
//! from Earth's orbital motion around the solar system barycenter.
//!
//! # Background
//!
//! TDB has no secular rate relative to TT; it differs from TT by periodic terms from
//! gravitational time dilation and velocity effects. The difference TDB-TT is dominated by a ~1.7 ms amplitude term with a period
//! of one year, plus smaller terms. For most applications, TDB ≈ TT to within 2 ms.
//!
//! The IAU recommends using TCB (Barycentric Coordinate Time) for rigorous relativistic
//! work. TDB is defined as a linear transformation of TCB that keeps TDB-TT bounded.
//!
//! # Usage
//!
//! ```
//! use celestial_time::julian::JulianDate;
//! use celestial_time::scales::tdb::TDB;
//!
//! let tdb = TDB::j2000();
//! let tdb_plus_day = tdb.add_days(1.0);
//!
//! let tdb_from_cal = celestial_time::scales::tdb::tdb_from_calendar(2000, 1, 1, 12, 0, 0.0).unwrap();
//! ```
//!
//! # Precision
//!
//! Internally stores time as a split Julian Date for sub-microsecond precision.
//! Arithmetic operations preserve precision by operating on the underlying JulianDate.

time_scale! {
    /// Barycentric Dynamical Time.
    ///
    /// A time scale for solar system barycentric ephemerides. Wraps a split Julian Date
    /// for high-precision arithmetic. TDB tracks TT to within ~2 ms over centuries.
    TDB
}

uniform_day_scale!(TDB, tdb_from_calendar);
