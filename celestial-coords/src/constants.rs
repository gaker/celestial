/// Rotation matrix for transforming Galactic coordinates to ICRS.
///
/// This matrix represents the transformation from the IAU 1958 Galactic coordinate system
/// to the International Celestial Reference System (ICRS). The matrix is derived from
/// the IAU-defined Galactic pole and zero-point:
/// - North Galactic Pole (NGP): RA = 192.859508°, Dec = 27.128336° (J2000/ICRS)
/// - Galactic Center direction: l=0°, b=0° points toward RA = 266.405°, Dec = -28.936° (J2000)
///
/// Reference: Liu, J.-C., Zhu, Z., & Zhang, H. (2011). "Reconsidering the Galactic
/// coordinate system". Astronomy & Astrophysics, 526, A16.
/// See also: ERFA function eraG2icrs documentation
#[allow(clippy::excessive_precision)]
pub const GALACTIC_TO_ICRS: [[f64; 3]; 3] = [
    [
        -0.054875560416215368492398900454,
        -0.873437090234885048760383168409,
        -0.483835015548713226831774175116,
    ],
    [
        0.494109427875583673525222371358,
        -0.444829629960011178146614061616,
        0.746982244497218890527388004556,
    ],
    [
        -0.867666149019004701181616534570,
        -0.198076373431201528180486091412,
        0.455983776175066922272100478348,
    ],
];
