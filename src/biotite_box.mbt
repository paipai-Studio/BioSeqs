/// biotite_box.mbt — Periodic boundary conditions / unit cell operations

///|
fn biotite_sqrt(x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.pow(x, 0.5)
  }
}

///|
fn biotite_dabs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}
/// Inspired by biotite.structure.Box and biotite.structure.pbc

///|
/// A box (unit cell) represented by three edge vectors.
/// Box vectors: [a, b, c] where each is a 3D vector.
pub struct BiotiteBox {
  // Three edge vectors, stored as flat array: [ax, ay, az, bx, by, bz, cx, cy, cz]
  vectors : Array[Double]
}

///|
/// Create box from three edge vectors.
pub fn BiotiteBox::from_vectors(
  a : Array[Double],
  b : Array[Double],
  c : Array[Double],
) -> BiotiteBox {
  { vectors: [a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]], }
}

///|
/// Create box from unit cell parameters (a, b, c, alpha, beta, gamma).
/// a, b, c: edge lengths; alpha, beta, gamma: angles in degrees.
pub fn BiotiteBox::from_unitcell(
  a : Double,
  b : Double,
  c : Double,
  alpha : Double,
  beta : Double,
  gamma : Double,
) -> BiotiteBox {
  let pi = 3.14159265358979323846
  let deg2rad = pi / 180.0
  let cos_a = @math.cos(alpha * deg2rad)
  let cos_b = @math.cos(beta * deg2rad)
  let cos_g = @math.cos(gamma * deg2rad)
  let sin_g = @math.sin(gamma * deg2rad)
  // Standard crystallographic setup
  let ax = a
  let ay = 0.0
  let az = 0.0
  let bx = b * cos_g
  let by = b * sin_g
  let bz = 0.0
  let cx = c * cos_b
  let cy = c * (cos_a - cos_b * cos_g) / sin_g
  let cz_val = c * c - cx * cx - cy * cy
  let cz = if cz_val > 0.0 { biotite_sqrt(cz_val) } else { 0.0 }
  { vectors: [ax, ay, az, bx, by, bz, cx, cy, cz], }
}

///|
/// Extract unit cell parameters from box vectors.
/// Returns (a, b, c, alpha, beta, gamma) with angles in degrees.
pub fn BiotiteBox::to_unitcell(
  self : BiotiteBox,
) -> (Double, Double, Double, Double, Double, Double) {
  let v = self.vectors
  let ax = v[0]
  let ay = v[1]
  let az = v[2]
  let bx = v[3]
  let by = v[4]
  let bz = v[5]
  let cx = v[6]
  let cy = v[7]
  let cz = v[8]
  let a = biotite_sqrt(ax * ax + ay * ay + az * az)
  let b = biotite_sqrt(bx * bx + by * by + bz * bz)
  let c = biotite_sqrt(cx * cx + cy * cy + cz * cz)
  let pi = 3.14159265358979323846
  let rad2deg = 180.0 / pi
  // cos(alpha) = (b·c) / (|b||c|)
  let cos_alpha = if b > 0.0 && c > 0.0 {
    (bx * cx + by * cy + bz * cz) / (b * c)
  } else {
    0.0
  }
  let cos_beta = if a > 0.0 && c > 0.0 {
    (ax * cx + ay * cy + az * cz) / (a * c)
  } else {
    0.0
  }
  let cos_gamma = if a > 0.0 && b > 0.0 {
    (ax * bx + ay * by + az * bz) / (a * b)
  } else {
    0.0
  }
  // Clamp to [-1, 1] for numerical safety
  let clamp = fn(x : Double) -> Double {
    if x > 1.0 {
      1.0
    } else if x < -1.0 {
      -1.0
    } else {
      x
    }
  }
  let alpha = @math.acos(clamp(cos_alpha)) * rad2deg
  let beta = @math.acos(clamp(cos_beta)) * rad2deg
  let gamma = @math.acos(clamp(cos_gamma)) * rad2deg
  (a, b, c, alpha, beta, gamma)
}

///|
/// Calculate box volume (scalar triple product).
pub fn BiotiteBox::volume(self : BiotiteBox) -> Double {
  let v = self.vectors
  // a · (b × c)
  let cross_x = v[4] * v[8] - v[5] * v[7]
  let cross_y = v[5] * v[6] - v[3] * v[8]
  let cross_z = v[3] * v[7] - v[4] * v[6]
  v[0] * cross_x + v[1] * cross_y + v[2] * cross_z
}

///|
/// Remove periodic boundary conditions from coordinates.
/// Wraps coordinates into the box [0, box) using fractional coordinates.
pub fn BiotiteBox::remove_pbc(
  self : BiotiteBox,
  coords : Array[Double],
) -> Array[Double] {
  let n = coords.length()
  let result : Array[Double] = Array::make(n, 0.0)
  let inv = self.inverse_matrix()
  for i = 0; i < n / 3; i = i + 1 {
    let x = coords[i * 3]
    let y = coords[i * 3 + 1]
    let z = coords[i * 3 + 2]
    // Convert to fractional
    let fx = inv[0] * x + inv[1] * y + inv[2] * z
    let fy = inv[3] * x + inv[4] * y + inv[5] * z
    let fz = inv[6] * x + inv[7] * y + inv[8] * z
    // Wrap to [0, 1)
    let wfx = fx - @math.floor(fx)
    let wfy = fy - @math.floor(fy)
    let wfz = fz - @math.floor(fz)
    // Convert back to Cartesian
    result[i * 3] = self.vectors[0] * wfx +
      self.vectors[3] * wfy +
      self.vectors[6] * wfz
    result[i * 3 + 1] = self.vectors[1] * wfx +
      self.vectors[4] * wfy +
      self.vectors[7] * wfz
    result[i * 3 + 2] = self.vectors[2] * wfx +
      self.vectors[5] * wfy +
      self.vectors[8] * wfz
  }
  result
}

///|
/// Move coordinates inside the box (alias for remove_pbc).
pub fn BiotiteBox::move_inside_box(
  self : BiotiteBox,
  coords : Array[Double],
) -> Array[Double] {
  self.remove_pbc(coords)
}

///|
/// Convert Cartesian coordinates to fractional coordinates.
pub fn BiotiteBox::coord_to_fraction(
  self : BiotiteBox,
  coords : Array[Double],
) -> Array[Double] {
  let inv = self.inverse_matrix()
  let n = coords.length()
  let result : Array[Double] = Array::make(n, 0.0)
  for i = 0; i < n / 3; i = i + 1 {
    let x = coords[i * 3]
    let y = coords[i * 3 + 1]
    let z = coords[i * 3 + 2]
    result[i * 3] = inv[0] * x + inv[1] * y + inv[2] * z
    result[i * 3 + 1] = inv[3] * x + inv[4] * y + inv[5] * z
    result[i * 3 + 2] = inv[6] * x + inv[7] * y + inv[8] * z
  }
  result
}

///|
/// Convert fractional coordinates to Cartesian coordinates.
pub fn BiotiteBox::fraction_to_coord(
  self : BiotiteBox,
  frac : Array[Double],
) -> Array[Double] {
  let n = frac.length()
  let result : Array[Double] = Array::make(n, 0.0)
  for i = 0; i < n / 3; i = i + 1 {
    let fx = frac[i * 3]
    let fy = frac[i * 3 + 1]
    let fz = frac[i * 3 + 2]
    result[i * 3] = self.vectors[0] * fx +
      self.vectors[3] * fy +
      self.vectors[6] * fz
    result[i * 3 + 1] = self.vectors[1] * fx +
      self.vectors[4] * fy +
      self.vectors[7] * fz
    result[i * 3 + 2] = self.vectors[2] * fx +
      self.vectors[5] * fy +
      self.vectors[8] * fz
  }
  result
}

///|
/// Generate repeat box: scale box vectors by integer factors.
pub fn BiotiteBox::repeat_box(
  self : BiotiteBox,
  na : Int32,
  nb : Int32,
  nc : Int32,
) -> BiotiteBox {
  {
    vectors: [
      self.vectors[0] * na.to_double(),
      self.vectors[1] * na.to_double(),
      self.vectors[2] * na.to_double(),
      self.vectors[3] * nb.to_double(),
      self.vectors[4] * nb.to_double(),
      self.vectors[5] * nb.to_double(),
      self.vectors[6] * nc.to_double(),
      self.vectors[7] * nc.to_double(),
      self.vectors[8] * nc.to_double(),
    ],
  }
}

///|
/// Check if the box is orthogonal (all angles 90 degrees).
pub fn BiotiteBox::is_orthogonal(self : BiotiteBox) -> Bool {
  let (_, _, _, alpha, beta, gamma) = self.to_unitcell()
  let tol = 0.001
  biotite_dabs(alpha - 90.0) < tol &&
  biotite_dabs(beta - 90.0) < tol &&
  biotite_dabs(gamma - 90.0) < tol
}

///|
/// Compute inverse of the 3x3 box matrix.
fn BiotiteBox::inverse_matrix(self : BiotiteBox) -> Array[Double] {
  let v = self.vectors
  let det = self.volume()
  if biotite_dabs(det) < 0.0000000001 {
    return Array::make(9, 0.0)
  }
  // Cofactor matrix / determinant
  let inv_det = 1.0 / det
  [
    (v[4] * v[8] - v[5] * v[7]) * inv_det,
    (v[2] * v[7] - v[1] * v[8]) * inv_det,
    (v[1] * v[5] - v[2] * v[4]) * inv_det,
    (v[5] * v[6] - v[3] * v[8]) * inv_det,
    (v[0] * v[8] - v[2] * v[6]) * inv_det,
    (v[2] * v[3] - v[0] * v[5]) * inv_det,
    (v[3] * v[7] - v[4] * v[6]) * inv_det,
    (v[1] * v[6] - v[0] * v[7]) * inv_det,
    (v[0] * v[4] - v[1] * v[3]) * inv_det,
  ]
}
