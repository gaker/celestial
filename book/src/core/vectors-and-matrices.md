# Vectors & Rotation Matrices

Two types handle 3D geometry: `Vector3` for positions and directions, `RotationMatrix3`
for frame transformations. Both are re-exported at the crate root.

## Vector3

A 3D Cartesian vector with public `x`, `y`, `z` fields.

### Construction

```rust,ignore
use celestial_core::Vector3;

// From components
let v = Vector3::new(1.0, 2.0, 3.0);

// Unit vectors along axes
let x = Vector3::x_axis();  // [1, 0, 0]
let y = Vector3::y_axis();  // [0, 1, 0]
let z = Vector3::z_axis();  // [0, 0, 1]

// Zero vector
let zero = Vector3::zeros(); // [0, 0, 0]

// From/to arrays
let v = Vector3::from_array([1.0, 2.0, 3.0]);
let arr: [f64; 3] = v.to_array();
```

### Spherical Coordinates

Convert between Cartesian and spherical (RA/Dec-style) representations.
`from_spherical` always produces a unit vector. Both angles are in radians.

```rust,ignore
use celestial_core::Vector3;
use celestial_core::constants::HALF_PI;

// RA=0, Dec=0 → points along +X
let v = Vector3::from_spherical(0.0, 0.0);
// v = [1, 0, 0]

// RA=90°, Dec=0 → points along +Y
let v = Vector3::from_spherical(HALF_PI, 0.0);
// v ≈ [0, 1, 0]

// Dec=90° → points along +Z (north celestial pole)
let v = Vector3::from_spherical(0.0, HALF_PI);
// v ≈ [0, 0, 1]
```

`to_spherical` returns `(theta, phi)` where theta is the azimuthal angle
(like RA, range `(-π, π]`) and phi is the elevation (like Dec, range `[-π/2, π/2]`).
The vector does not need to be normalized:

```rust,ignore
use celestial_core::Vector3;

let v = Vector3::new(0.0, 0.0, 1.0); // north pole
let (theta, phi) = v.to_spherical();
// theta = 0.0, phi = π/2
```

The convention matches standard astronomy: theta is azimuthal from +X toward +Y
(like RA), phi is elevation from the XY plane (like Dec). This differs from the
physics convention where the names are swapped.

### Magnitude and Normalization

```rust,ignore
use celestial_core::Vector3;

let v = Vector3::new(3.0, 4.0, 0.0);

v.magnitude();         // 5.0
v.magnitude_squared(); // 25.0 (cheaper, good for comparisons)

let unit = v.normalize(); // [0.6, 0.8, 0.0], magnitude = 1.0
```

`normalize` returns the zero vector unchanged (avoids NaN).

### Dot and Cross Products

For unit vectors, the dot product is the cosine of the angle between them.
The cross product gives the perpendicular axis (right-hand rule):

```rust,ignore
use celestial_core::Vector3;

let a = Vector3::x_axis();
let b = Vector3::y_axis();

a.dot(&b);    // 0.0 (perpendicular)
a.cross(&b);  // Vector3 { x: 0, y: 0, z: 1 } (= z_axis)

let c = Vector3::new(1.0, 2.0, 3.0);
let d = Vector3::new(4.0, 5.0, 6.0);
c.dot(&d); // 32.0 (1*4 + 2*5 + 3*6)
```

### Element Access

Fields are public, but index-based access is also available:

```rust,ignore
use celestial_core::Vector3;

let mut v = Vector3::new(1.0, 2.0, 3.0);

// Direct field access
v.x; // 1.0

// Index access (panics if index > 2)
v[0]; // 1.0 (x)
v[1]; // 2.0 (y)
v[2]; // 3.0 (z)

// Mutable indexing
v[0] = 10.0;

// Checked access (returns Result)
v.get(0).unwrap();      // 10.0
v.set(1, 20.0).unwrap();
v.get(3); // Err — index out of bounds
```

### Arithmetic

```rust,ignore
use celestial_core::Vector3;

let a = Vector3::new(1.0, 2.0, 3.0);
let b = Vector3::new(4.0, 5.0, 6.0);

let sum  = a + b;     // [5, 7, 9]
let diff = b - a;     // [3, 3, 3]
let scaled = a * 2.0; // [2, 4, 6]
let also   = 3.0 * a; // [3, 6, 9]  (scalar * vector works too)
let half = a / 2.0;   // [0.5, 1.0, 1.5]
let neg  = -a;         // [-1, -2, -3]

// In-place division
let mut v = Vector3::new(10.0, 20.0, 30.0);
v /= 2.0; // [5, 10, 15]
```

### Display

Formats as `Vector3(x, y, z)` with 9 decimal places:

```rust,ignore
let v = Vector3::new(1.0, -2.5, 3.0);
println!("{}", v);
// Vector3(1.000000000, -2.500000000, 3.000000000)
```

---

## RotationMatrix3

A 3x3 rotation matrix stored row-major. Elements are private; access via methods
or indexing.

### Construction

Start from identity and build up rotations, or provide elements directly:

```rust,ignore
use celestial_core::RotationMatrix3;

// Identity matrix (no rotation)
let m = RotationMatrix3::identity();

// From a 3x3 array (row-major)
let m = RotationMatrix3::from_array([
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
]);
```

### Building Rotations

Three in-place methods apply rotations about the principal axes. Each modifies
the matrix to become `R_axis(angle) * self`. Angles are in radians:

```rust,ignore
use celestial_core::RotationMatrix3;

let mut m = RotationMatrix3::identity();
m.rotate_z(0.1);   // Apply Rz(0.1 rad)
m.rotate_y(-0.05); // Then Ry(-0.05 rad) — now m = Ry * Rz
m.rotate_x(0.02);  // Then Rx(0.02 rad) — now m = Rx * Ry * Rz
```

Rotations follow the ERFA convention: positive angles rotate counterclockwise
when looking from the positive axis toward the origin. This is the passive/alias
convention where we rotate the coordinate frame, not the vector.

What this means concretely:

```rust,ignore
use celestial_core::RotationMatrix3;
use celestial_core::constants::HALF_PI;

let mut m = RotationMatrix3::identity();
m.rotate_z(HALF_PI);

// [1, 0, 0] rotates to [0, -1, 0]
let v = m.apply_to_vector([1.0, 0.0, 0.0]);
// v ≈ [0.0, -1.0, 0.0]
```

### Applying to Vectors

Two ways to transform a vector:

```rust,ignore
use celestial_core::{RotationMatrix3, Vector3};

let mut m = RotationMatrix3::identity();
m.rotate_z(0.5);

// Method: takes [f64; 3], returns [f64; 3]
let result = m.apply_to_vector([1.0, 0.0, 0.0]);

// Operator: takes Vector3, returns Vector3
let v = Vector3::new(1.0, 0.0, 0.0);
let result = m * v;
let result = &m * v; // also works with references
```

### Composing Rotations

Matrix multiplication composes rotations. The rightmost matrix acts first:

```rust,ignore
use celestial_core::RotationMatrix3;

let mut rx = RotationMatrix3::identity();
rx.rotate_x(0.1);

let mut rz = RotationMatrix3::identity();
rz.rotate_z(0.2);

// "First apply rx, then rz" = rz * rx
let combined = rz * rx;

// .multiply() does the same thing
let combined = rz.multiply(&rx);

// All reference combinations work: m*m, m*&m, &m*m, &m*&m
```

### Transforming Spherical Coordinates

The common case — transform RA/Dec (or lon/lat) through a rotation without
manually converting to/from Cartesian:

```rust,ignore
use celestial_core::RotationMatrix3;
use celestial_core::constants::QUARTER_PI;

let mut m = RotationMatrix3::identity();
m.rotate_z(QUARTER_PI); // 45° rotation

let (ra, dec) = (0.0, 0.0); // equator at RA=0
let (new_ra, new_dec) = m.transform_spherical(ra, dec);
// new_ra ≈ -π/4 (-45°), new_dec ≈ 0.0
```

Input and output angles are in radians. Output RA is in `(-π, π]`, output
Dec is in `[-π/2, π/2]`.

### Transpose (Inverse)

For a proper rotation matrix, the transpose equals the inverse. This is exact
and numerically stable:

```rust,ignore
use celestial_core::RotationMatrix3;

let mut m = RotationMatrix3::identity();
m.rotate_z(0.5);
m.rotate_x(0.3);

let m_inv = m.transpose();

// Round-trip: m * m_inv returns to original
let v = [1.0, 2.0, 3.0];
let rotated  = m.apply_to_vector(v);
let restored = m_inv.apply_to_vector(rotated);
// restored ≈ [1.0, 2.0, 3.0]
```

### Validation and Comparison

```rust,ignore
use celestial_core::RotationMatrix3;

let mut m = RotationMatrix3::identity();
m.rotate_z(0.5);

// Check if it's a proper rotation (det=+1, orthogonal)
m.is_rotation_matrix(1e-14); // true

m.determinant(); // 1.0 for proper rotations

// Compare two matrices element-wise
let a = RotationMatrix3::identity();
let b = RotationMatrix3::identity();
a.max_difference(&b); // 0.0
```

### Element Access

```rust,ignore
use celestial_core::RotationMatrix3;

let mut m = RotationMatrix3::identity();

// Method access
m.get(0, 0); // 1.0
m.set(0, 1, 0.5);

// Tuple indexing
m[(0, 0)]; // 1.0
m[(0, 1)] = 0.5;

// Direct array reference
let elems: &[[f64; 3]; 3] = m.elements();
```

### Display

Formats as a labeled multi-line matrix with 9 decimal places:

```rust,ignore
use celestial_core::RotationMatrix3;

let m = RotationMatrix3::identity();
println!("{}", m);
// RotationMatrix3:
//   [ 1.000000000  0.000000000  0.000000000]
//   [ 0.000000000  1.000000000  0.000000000]
//   [ 0.000000000  0.000000000  1.000000000]
```
