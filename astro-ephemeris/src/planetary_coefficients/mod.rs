//! VSOP2013 planetary coefficients
//!
//! Truncated coefficient tables for analytical planetary ephemeris.
//! These coefficients are used to compute heliocentric positions of planets.

pub mod emb;
pub mod jupiter;
pub mod mars;
pub mod mercury;
pub mod neptune;
pub mod pluto;
pub mod saturn;
pub mod uranus;
pub mod venus;

/// A single Fourier term in the VSOP2013 series
#[derive(Debug, Clone, Copy)]
pub struct Term {
    /// Sine coefficient
    pub s: f64,
    /// Cosine coefficient
    pub c: f64,
    /// Argument multipliers for 17 fundamental arguments
    pub mult: [i16; 17],
}

/// Terms grouped by power of T (time)
#[derive(Debug, Clone, Copy)]
pub struct TimeBlock {
    /// Power of T (0, 1, 2, ...)
    pub power: u8,
    /// Terms for this power, sorted by amplitude descending
    pub terms: &'static [Term],
}
