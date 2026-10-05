//! VSOP2013 planetary coefficients
//!
//! Truncated coefficient tables for analytical planetary ephemeris.
//! These coefficients are used to compute heliocentric positions of planets.

pub(crate) mod emb;
pub(crate) mod jupiter;
pub(crate) mod mars;
pub(crate) mod mercury;
pub(crate) mod neptune;
pub(crate) mod pluto;
pub(crate) mod saturn;
pub(crate) mod uranus;
pub(crate) mod venus;

/// A single Fourier term in the VSOP2013 series
#[derive(Debug, Clone, Copy)]
pub(crate) struct Term {
    /// Sine coefficient
    pub(crate) s: f64,
    /// Cosine coefficient
    pub(crate) c: f64,
    // The nonzero multipliers, zero-padded, and which of the 17 fundamental
    // arguments each one multiplies.
    pub(crate) mult: [i32; 6],
    pub(crate) index: [u8; 6],
}

/// Terms grouped by power of T (time)
#[derive(Debug, Clone, Copy)]
pub(crate) struct TimeBlock {
    /// Power of T (0, 1, 2, ...)
    pub(crate) power: u8,
    /// Terms for this power, sorted by amplitude descending
    pub(crate) terms: &'static [Term],
}
