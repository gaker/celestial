# Celestial

Pure-Rust astronomy library: time scales, IAU coordinate chain, ephemeris, WCS, FITS/XISF/SER I/O, telescope pointing models, star catalog, plate solving. Results are expected to match ERFA bit-for-bit.

## Crate layering

`core` → `time` → `coords` → `ephemeris` · `pointing` · `catalog` · `wcs` → `images` → `solver`

`celestial-core` depends on nothing internal. Keep it that way.

## Rules

- **Math:** `libm::sin(x)`, never `x.sin()`. Constants come from `celestial_core::constants`, never `std::f64::consts`. Check `celestial_core::math` before writing a helper.
- **Precision:** tests use exact `assert_eq!` against reference values. A mismatch is a bug, so find the cause. Adding or widening any tolerance (including `assert_ulp_le!`) needs approval first.
- **Imports:** no prelude, and crate roots re-export nothing. Callers import from the defining module, e.g. `celestial_core::matrix::{RotationMatrix3, Vector3}`. This overrides the global "use a prelude" rule.
- **Visibility:** private, then `pub(super)` / `pub(crate)`. Use `pub` only for real API.
- **Errors:** `thiserror`, and no panics in library code. NaN or non-finite input returns an error rather than garbage.
- **Docs:**
  - Don't add `///` comments. User docs live in `book/` (mdBook, `rust,ignore` examples).
  - Edit existing doc comments only when they're wrong.

## Commands

```sh
cargo test --workspace
cargo clippy --workspace -- -D warnings            # CI runs both
cargo llvm-cov --features serde -p celestial-core  # coverage
rustfmt --edition 2021 path/to/file.rs             # format one file
```

The toolchain is pinned in `rust-toolchain.toml`. Optional features: `serde` (most crates), `test-support` (core), `cli` (catalog, ephemeris), `parallel`/`simd`/`standard-formats` (images).

## Reference material

`references/` holds the IERS conventions and format specs. `vsop2013/` holds the planetary series data.
