# celestial-pointing

Telescope pointing model fitting and correction.

[![Crates.io](https://img.shields.io/crates/v/celestial-pointing)](https://crates.io/crates/celestial-pointing)
[![Documentation](https://docs.rs/celestial-pointing/badge.svg)](https://docs.rs/celestial-pointing)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/celestial-pointing)](https://github.com/gaker/celestial)

Build, fit, and apply telescope pointing models using standard equatorial and harmonic terms. Interactive REPL with TPOINT-compatible workflow: load observations, fit models via least-squares, analyze residuals with plots, and export corrections. No runtime FFI.

## Installation

```toml
[dependencies]
celestial-pointing = "0.1"
```

Or run the interactive CLI:

```bash
cargo install celestial-pointing
pointing
```

## Modules

| Module        | Purpose                                                    |
|---------------|------------------------------------------------------------|
| `observation` | Observation struct (catalog, observed, commanded coords)   |
| `model`       | PointingModel with coefficient management and application  |
| `terms`       | 6 base + 8 physical + 96 harmonic term implementations     |
| `solver`      | Weighted least-squares fitting via nalgebra SVD            |
| `session`     | Interactive session state (observations, model, fit)       |
| `commands`    | 30+ REPL commands (FIT, SHOW, OPTIMAL, plots, etc.)        |
| `plot`        | SVG and terminal ASCII residual visualization              |

## Pointing Terms

| Category  | Terms                                                        |
|-----------|--------------------------------------------------------------|
| Base      | IH, ID, CH, NP, MA, ME (index, collimation, polar alignment) |
| Physical  | TF, TX, DAF, FO, HCES, HCEC, DCES, DCEC (flexure, fork)      |
| Harmonic  | HxSy, HxCy patterns for periodic errors (x=H/D/X, y=H/D, 1-8)|

## Example Session

```sh
pointing> INDAT observations.dat
Loaded 47 observations
pointing> USE IH ID CH NP MA ME
6 terms active
pointing> FIT
RMS = 4.23"  (was 127.8")
pointing> OPTIMAL
Base: IH ID CH NP MA ME (BIC=312.4, RMS=4.23")
+ TF (dBIC=-18.2, RMS=3.41")
+ HHSH (dBIC=-7.1, RMS=3.12")
Final model: 8 terms, RMS=3.12"
pointing> GSCAT residuals.svg
Written to residuals.svg
pointing> OUTMOD model.dat
```

## Commands

| Command   | Purpose                                      |
|-----------|----------------------------------------------|
| `INDAT`   | Load observation data file                   |
| `USE`     | Add terms to active model                    |
| `FIT`     | Fit model to observations                    |
| `SHOW`    | Display current model coefficients           |
| `OPTIMAL` | Auto-build model using BIC selection         |
| `OUTMOD`  | Export model coefficients                    |
| `CORRECT` | Compute correction for target coordinates    |
| `PREDICT` | Show per-term correction breakdown           |
| `GSCAT`   | Scatter plot (dX vs dDec)                    |
| `GDIST`   | Histogram of residual distributions          |
| `GMAP`    | Sky map with residual vectors                |
| `GHA`     | Residuals vs hour angle                      |
| `GDEC`    | Residuals vs declination                     |
| `GHYST`   | Hysteresis plot by pier side                 |

## Design Notes

- **OPTIMAL uses BIC**: Forward stepwise selection with Bayesian Information Criterion prevents overfitting. Terms must improve BIC by at least -6.0 to be added.
- **Parallel harmonic search**: Candidate evaluation uses rayon for fast model selection across 96 harmonic terms.
- **Dual output modes**: All plot commands support terminal ASCII (no args) or SVG file output (with path argument).
- **Pier-side aware**: Observations track East/West pier side for hysteresis analysis on German equatorial mounts.

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT License

## Contributing

See the [repository](https://github.com/gaker/celestial) for contribution guidelines.
