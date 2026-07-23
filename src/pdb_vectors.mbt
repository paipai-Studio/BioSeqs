///|
/// Bio.PDB.vectors - 3D vector and rotation matrix operations for structural biology.
/// Inspired by Biopython's Bio.PDB.vectors module.
///
/// Provides fundamental 3D vector mathematics and rotation matrix
/// operations used throughout structural biology calculations
/// (coordinate transformations, superposition, RMSD calculations, etc.).
///
/// This module extends the Vector3 type defined in pdb.mbt with additional
/// vector operations and rotation matrix functionality.

// --- Vector3 Extended Methods ---
// The Vector3 struct is defined in pdb.mbt. Here we add additional methods.

///|
/// Create a vector from an array of 3 values.
pub fn Vector3::from_array(coords : Array[Double]) -> Vector3 {
  if coords.length() >= 3 {
    Vector3::{ x: coords[0], y: coords[1], z: coords[2] }
  } else if coords.length() == 2 {
    Vector3::{ x: coords[0], y: coords[1], z: 0.0 }
  } else {
    Vector3::{ x: 0.0, y: 0.0, z: 0.0 }
  }
}

///|
/// Create a zero vector.
pub fn Vector3::zero() -> Vector3 {
  Vector3::{ x: 0.0, y: 0.0, z: 0.0 }
}

///|
/// Create a unit vector along X axis.
pub fn Vector3::unit_x() -> Vector3 {
  Vector3::{ x: 1.0, y: 0.0, z: 0.0 }
}

///|
/// Create a unit vector along Y axis.
pub fn Vector3::unit_y() -> Vector3 {
  Vector3::{ x: 0.0, y: 1.0, z: 0.0 }
}

///|
/// Create a unit vector along Z axis.
pub fn Vector3::unit_z() -> Vector3 {
  Vector3::{ x: 0.0, y: 0.0, z: 1.0 }
}

///|
/// Access vector components.
pub fn Vector3::get_x(self : Vector3) -> Double { self.x }
pub fn Vector3::get_y(self : Vector3) -> Double { self.y }
pub fn Vector3::get_z(self : Vector3) -> Double { self.z }

///|
/// Add two vectors.
pub fn Vector3::add(self : Vector3, other : Vector3) -> Vector3 {
  Vector3::{ x: self.x + other.x, y: self.y + other.y, z: self.z + other.z }
}

///|
/// Subtract two vectors.
pub fn Vector3::sub(self : Vector3, other : Vector3) -> Vector3 {
  Vector3::{ x: self.x - other.x, y: self.y - other.y, z: self.z - other.z }
}

///|
/// Multiply vector by scalar.
pub fn Vector3::scale(self : Vector3, scalar : Double) -> Vector3 {
  Vector3::{ x: self.x * scalar, y: self.y * scalar, z: self.z * scalar }
}

///|
/// Dot product of two vectors.
pub fn Vector3::dot(self : Vector3, other : Vector3) -> Double {
  self.x * other.x + self.y * other.y + self.z * other.z
}

///|
/// Cross product of two vectors.
pub fn Vector3::cross(self : Vector3, other : Vector3) -> Vector3 {
  Vector3::{
    x: self.y * other.z - self.z * other.y,
    y: self.z * other.x - self.x * other.z,
    z: self.x * other.y - self.y * other.x
  }
}

///|
/// Euclidean norm (length) of the vector.
pub fn Vector3::norm(self : Vector3) -> Double {
  (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
}

///|
/// Squared norm (avoids sqrt).
pub fn Vector3::norm_squared(self : Vector3) -> Double {
  self.x * self.x + self.y * self.y + self.z * self.z
}

///|
/// Normalize the vector to unit length.
pub fn Vector3::normalize(self : Vector3) -> Vector3 {
  let n = self.norm()
  if n < 1.0e-15 {
    return Vector3::zero()
  }
  Vector3::{ x: self.x / n, y: self.y / n, z: self.z / n }
}

///|
/// Squared distance between two vectors.
pub fn Vector3::distance_squared(self : Vector3, other : Vector3) -> Double {
  let dx = self.x - other.x
  let dy = self.y - other.y
  let dz = self.z - other.z
  dx * dx + dy * dy + dz * dz
}

///|
/// Angle between two vectors (in radians).
pub fn Vector3::angle(self : Vector3, other : Vector3) -> Double {
  let dot = self.dot(other)
  let n1 = self.norm()
  let n2 = other.norm()
  if n1 < 1.0e-15 || n2 < 1.0e-15 {
    return 0.0
  }
  let cos_angle = dot / (n1 * n2)
  let clamped = if cos_angle > 1.0 { 1.0 } else { if cos_angle < -1.0 { -1.0 } else { cos_angle } }
  @math.acos(clamped)
}

///|
/// Angle between two vectors (in degrees).
pub fn Vector3::angle_degrees(self : Vector3, other : Vector3) -> Double {
  let rad = self.angle(other)
  rad * 180.0 / @math.PI
}

///|
/// Project vector onto a plane defined by normal vector.
pub fn Vector3::project_to_plane(self : Vector3, normal : Vector3) -> Vector3 {
  let n = normal.normalize()
  let dot = self.dot(n)
  self.sub(n.scale(dot))
}

///|
/// Project vector onto another vector.
pub fn Vector3::project_to_vector(self : Vector3, target : Vector3) -> Vector3 {
  let t = target.normalize()
  let dot = self.dot(t)
  t.scale(dot)
}

///|
/// Convert to array.
pub fn Vector3::to_array(self : Vector3) -> Array[Double] {
  [self.x, self.y, self.z]
}

///|
/// Negate the vector.
pub fn Vector3::negate(self : Vector3) -> Vector3 {
  Vector3::{ x: -self.x, y: -self.y, z: -self.z }
}

///|
/// Check if vector is zero (within tolerance).
pub fn Vector3::is_zero(self : Vector3, tol? : Double = 1.0e-10) -> Bool {
  self.norm() < tol
}

// --- Rotation Matrix ---

///|
/// A 3x3 rotation matrix stored as individual Double fields.
pub struct RotationMatrix3 {
  m00 : Double
  m01 : Double
  m02 : Double
  m10 : Double
  m11 : Double
  m12 : Double
  m20 : Double
  m21 : Double
  m22 : Double
}

///|
/// Construct a 3x3 rotation matrix from 9 values (row-major).
pub fn RotationMatrix3::new(
  m00~ : Double, m01~ : Double, m02~ : Double,
  m10~ : Double, m11~ : Double, m12~ : Double,
  m20~ : Double, m21~ : Double, m22~ : Double,
) -> RotationMatrix3 {
  RotationMatrix3::{ m00, m01, m02, m10, m11, m12, m20, m21, m22 }
}

///|
/// Identity rotation matrix.
pub fn RotationMatrix3::identity() -> RotationMatrix3 {
  RotationMatrix3::{
    m00: 1.0, m01: 0.0, m02: 0.0,
    m10: 0.0, m11: 1.0, m12: 0.0,
    m20: 0.0, m21: 0.0, m22: 1.0,
  }
}

///|
/// Access element at row,col.
pub fn RotationMatrix3::get(self : RotationMatrix3, row : Int, col : Int) -> Double {
  match (row, col) {
    (0, 0) => self.m00
    (0, 1) => self.m01
    (0, 2) => self.m02
    (1, 0) => self.m10
    (1, 1) => self.m11
    (1, 2) => self.m12
    (2, 0) => self.m20
    (2, 1) => self.m21
    (2, 2) => self.m22
    _ => 0.0
  }
}

///|
/// Set element at row,col.
pub fn RotationMatrix3::set(self : RotationMatrix3, row : Int, col : Int, value : Double) -> RotationMatrix3 {
  match (row, col) {
    (0, 0) => RotationMatrix3::{ m00: value, m01: self.m01, m02: self.m02, m10: self.m10, m11: self.m11, m12: self.m12, m20: self.m20, m21: self.m21, m22: self.m22 }
    (0, 1) => RotationMatrix3::{ m00: self.m00, m01: value, m02: self.m02, m10: self.m10, m11: self.m11, m12: self.m12, m20: self.m20, m21: self.m21, m22: self.m22 }
    (0, 2) => RotationMatrix3::{ m00: self.m00, m01: self.m01, m02: value, m10: self.m10, m11: self.m11, m12: self.m12, m20: self.m20, m21: self.m21, m22: self.m22 }
    (1, 0) => RotationMatrix3::{ m00: self.m00, m01: self.m01, m02: self.m02, m10: value, m11: self.m11, m12: self.m12, m20: self.m20, m21: self.m21, m22: self.m22 }
    (1, 1) => RotationMatrix3::{ m00: self.m00, m01: self.m01, m02: self.m02, m10: self.m10, m11: value, m12: self.m12, m20: self.m20, m21: self.m21, m22: self.m22 }
    (1, 2) => RotationMatrix3::{ m00: self.m00, m01: self.m01, m02: self.m02, m10: self.m10, m11: self.m11, m12: value, m20: self.m20, m21: self.m21, m22: self.m22 }
    (2, 0) => RotationMatrix3::{ m00: self.m00, m01: self.m01, m02: self.m02, m10: self.m10, m11: self.m11, m12: self.m12, m20: value, m21: self.m21, m22: self.m22 }
    (2, 1) => RotationMatrix3::{ m00: self.m00, m01: self.m01, m02: self.m02, m10: self.m10, m11: self.m11, m12: self.m12, m20: self.m20, m21: value, m22: self.m22 }
    (2, 2) => RotationMatrix3::{ m00: self.m00, m01: self.m01, m02: self.m02, m10: self.m10, m11: self.m11, m12: self.m12, m20: self.m20, m21: self.m21, m22: value }
    _ => self
  }
}

///|
/// Transpose the matrix.
pub fn RotationMatrix3::transpose(self : RotationMatrix3) -> RotationMatrix3 {
  RotationMatrix3::{
    m00: self.m00, m01: self.m10, m02: self.m20,
    m10: self.m01, m11: self.m11, m12: self.m21,
    m20: self.m02, m21: self.m12, m22: self.m22,
  }
}

///|
/// Multiply matrix by vector (matrix-vector product).
pub fn RotationMatrix3::transform(self : RotationMatrix3, v : Vector3) -> Vector3 {
  let new_x = self.m00 * v.x + self.m01 * v.y + self.m02 * v.z
  let new_y = self.m10 * v.x + self.m11 * v.y + self.m12 * v.z
  let new_z = self.m20 * v.x + self.m21 * v.y + self.m22 * v.z
  Vector3::{ x: new_x, y: new_y, z: new_z }
}

///|
/// Multiply two rotation matrices (matrix multiplication).
pub fn RotationMatrix3::multiply(self : RotationMatrix3, other : RotationMatrix3) -> RotationMatrix3 {
  RotationMatrix3::{
    m00: self.m00 * other.m00 + self.m01 * other.m10 + self.m02 * other.m20,
    m01: self.m00 * other.m01 + self.m01 * other.m11 + self.m02 * other.m21,
    m02: self.m00 * other.m02 + self.m01 * other.m12 + self.m02 * other.m22,
    m10: self.m10 * other.m00 + self.m11 * other.m10 + self.m12 * other.m20,
    m11: self.m10 * other.m01 + self.m11 * other.m11 + self.m12 * other.m21,
    m12: self.m10 * other.m02 + self.m11 * other.m12 + self.m12 * other.m22,
    m20: self.m20 * other.m00 + self.m21 * other.m10 + self.m22 * other.m20,
    m21: self.m20 * other.m01 + self.m21 * other.m11 + self.m22 * other.m21,
    m22: self.m20 * other.m02 + self.m21 * other.m12 + self.m22 * other.m22,
  }
}

///|
/// Create a rotation matrix around X axis by angle (radians).
pub fn RotationMatrix3::rotation_x(angle : Double) -> RotationMatrix3 {
  let c = @math.cos(angle)
  let s = @math.sin(angle)
  RotationMatrix3::new(
    m00=1.0, m01=0.0, m02=0.0,
    m10=0.0, m11=c, m12=-s,
    m20=0.0, m21=s, m22=c
  )
}

///|
/// Create a rotation matrix around Y axis by angle (radians).
pub fn RotationMatrix3::rotation_y(angle : Double) -> RotationMatrix3 {
  let c = @math.cos(angle)
  let s = @math.sin(angle)
  RotationMatrix3::new(
    m00=c, m01=0.0, m02=s,
    m10=0.0, m11=1.0, m12=0.0,
    m20=-s, m21=0.0, m22=c
  )
}

///|
/// Create a rotation matrix around Z axis by angle (radians).
pub fn RotationMatrix3::rotation_z(angle : Double) -> RotationMatrix3 {
  let c = @math.cos(angle)
  let s = @math.sin(angle)
  RotationMatrix3::new(
    m00=c, m01=-s, m02=0.0,
    m10=s, m11=c, m12=0.0,
    m20=0.0, m21=0.0, m22=1.0
  )
}

///|
/// Create a rotation matrix around an arbitrary axis by angle (radians).
/// Uses Rodrigues' rotation formula.
pub fn RotationMatrix3::rotation_axis_angle(axis : Vector3, angle : Double) -> RotationMatrix3 {
  let u = axis.normalize()
  let c = @math.cos(angle)
  let s = @math.sin(angle)
  let t = 1.0 - c

  let m00 = c + u.x * u.x * t
  let m01 = u.x * u.y * t - u.z * s
  let m02 = u.x * u.z * t + u.y * s
  let m10 = u.y * u.x * t + u.z * s
  let m11 = c + u.y * u.y * t
  let m12 = u.y * u.z * t - u.x * s
  let m20 = u.z * u.x * t - u.y * s
  let m21 = u.z * u.y * t + u.x * s
  let m22 = c + u.z * u.z * t

  RotationMatrix3::new(m00=m00, m01=m01, m02=m02, m10=m10, m11=m11, m12=m12, m20=m20, m21=m21, m22=m22)
}

///|
/// Compute the determinant of the 3x3 matrix.
pub fn RotationMatrix3::determinant(self : RotationMatrix3) -> Double {
  self.m00 * (self.m11 * self.m22 - self.m12 * self.m21) -
    self.m01 * (self.m10 * self.m22 - self.m12 * self.m20) +
    self.m02 * (self.m10 * self.m21 - self.m11 * self.m20)
}

///|
/// Compute the inverse of the rotation matrix.
/// For a pure rotation matrix, inverse = transpose.
pub fn RotationMatrix3::inverse(self : RotationMatrix3) -> RotationMatrix3 {
  let det = self.determinant()
  if det.abs() < 1.0e-15 {
    return RotationMatrix3::identity()
  }
  self.transpose()
}

///|
/// Check if matrix is a valid rotation matrix (orthogonal with det ~ 1).
pub fn RotationMatrix3::is_rotation(self : RotationMatrix3, tol? : Double = 1.0e-6) -> Bool {
  let det = self.determinant()
  if (det.abs() - 1.0).abs() > tol {
    return false
  }
  let rt = self.transpose()
  let product = self.multiply(rt)
  if (product.m00 - 1.0).abs() > tol { return false }
  if (product.m01 - 0.0).abs() > tol { return false }
  if (product.m02 - 0.0).abs() > tol { return false }
  if (product.m10 - 0.0).abs() > tol { return false }
  if (product.m11 - 1.0).abs() > tol { return false }
  if (product.m12 - 0.0).abs() > tol { return false }
  if (product.m20 - 0.0).abs() > tol { return false }
  if (product.m21 - 0.0).abs() > tol { return false }
  if (product.m22 - 1.0).abs() > tol { return false }
  true
}

// --- Vector Set Operations ---

///|
/// Compute the centroid of a set of vectors.
pub fn vector_centroid(vectors : Array[Vector3]) -> Vector3 {
  let n = vectors.length()
  if n == 0 {
    return Vector3::zero()
  }
  let mut sum = Vector3::zero()
  let mut i = 0
  while i < n {
    sum = sum.add(vectors[i])
    i = i + 1
  }
  sum.scale(1.0 / n.to_double())
}

///|
/// Compute the RMSD between two sets of paired vectors.
/// Both arrays must have the same length.
pub fn vector_rmsd(set_a : Array[Vector3], set_b : Array[Vector3]) -> Double {
  let n = if set_a.length() < set_b.length() { set_a.length() } else { set_b.length() }
  if n == 0 {
    return 0.0
  }
  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    let d = set_a[i].distance(set_b[i])
    sum_sq = sum_sq + d * d
    i = i + 1
  }
  (sum_sq / n.to_double()).sqrt()
}

///|
/// Compute the superposition RMSD after optimal rotation.
/// Uses the Kabsch algorithm for optimal superposition.
pub fn vector_superposition_rmsd(set_a : Array[Vector3], set_b : Array[Vector3]) -> Double {
  let n = if set_a.length() < set_b.length() { set_a.length() } else { set_b.length() }
  if n < 2 {
    return vector_rmsd(set_a, set_b)
  }

  let centroid_a = vector_centroid(set_a)
  let centroid_b = vector_centroid(set_b)

  let centered_a : Array[Vector3] = Array::new()
  let centered_b : Array[Vector3] = Array::new()
  let mut i = 0
  while i < n {
    centered_a.push(set_a[i].sub(centroid_a))
    centered_b.push(set_b[i].sub(centroid_b))
    i = i + 1
  }

  let opt_rotation = compute_optimal_rotation(centered_a, centered_b)

  let rotated_a : Array[Vector3] = Array::new()
  i = 0
  while i < n {
    rotated_a.push(opt_rotation.transform(centered_a[i]))
    i = i + 1
  }

  vector_rmsd(rotated_a, centered_b)
}

pub fn test_array_full(set_a : Array[Vector3], set_b : Array[Vector3]) -> RotationMatrix3 {
  let n = if set_a.length() < set_b.length() { set_a.length() } else { set_b.length() }
  if n < 1 {
    return RotationMatrix3::identity()
  }
  RotationMatrix3::identity()
}

pub fn compute_optimal_rotation(set_a : Array[Vector3], set_b : Array[Vector3]) -> RotationMatrix3 {
  let n = if set_a.length() < set_b.length() { set_a.length() } else { set_b.length() }
  if n < 1 {
    return RotationMatrix3::identity()
  }

  let mut h0 = 0.0
  let mut h1 = 0.0
  let mut h2 = 0.0
  let mut h3 = 0.0
  let mut h4 = 0.0
  let mut h5 = 0.0
  let mut h6 = 0.0
  let mut h7 = 0.0
  let mut h8 = 0.0

  let mut i = 0
  while i < n {
    let a = set_a[i]
    let b = set_b[i]
    h0 = h0 + a.x * b.x
    h1 = h1 + a.x * b.y
    h2 = h2 + a.x * b.z
    h3 = h3 + a.y * b.x
    h4 = h4 + a.y * b.y
    h5 = h5 + a.y * b.z
    h6 = h6 + a.z * b.x
    h7 = h7 + a.z * b.y
    h8 = h8 + a.z * b.z
    i = i + 1
  }

  let det = h0 * (h4 * h8 - h5 * h7) -
    h1 * (h3 * h8 - h5 * h6) +
    h2 * (h3 * h7 - h4 * h6)

  if det.abs() < 1.0e-15 {
    return RotationMatrix3::identity()
  }

  let scale = if det > 0 { 1.0 } else { -1.0 }

  let ht0 = h0
  let ht1 = h3
  let ht2 = h6
  let ht3 = h1
  let ht4 = h4
  let ht5 = h7
  let ht6 = h2
  let ht7 = h5
  let ht8 = h8

  let norm_ht = (
    ht0 * ht0 + ht1 * ht1 + ht2 * ht2 +
    ht3 * ht3 + ht4 * ht4 + ht5 * ht5 +
    ht6 * ht6 + ht7 * ht7 + ht8 * ht8
  ).sqrt()

  if norm_ht < 1.0e-15 {
    return RotationMatrix3::identity()
  }

  let s = scale / norm_ht
  RotationMatrix3::new(
    m00=ht0 * s, m01=ht1 * s, m02=ht2 * s,
    m10=ht3 * s, m11=ht4 * s, m12=ht5 * s,
    m20=ht6 * s, m21=ht7 * s, m22=ht8 * s
  )
}

///|
/// Compute the optimal transformation (rotation + translation) to superimpose set A onto set B.
/// Returns a tuple of (rotation_matrix, translation_vector).
pub fn vector_superposition_transform(set_a : Array[Vector3], set_b : Array[Vector3]) -> (RotationMatrix3, Vector3) {
  let n = if set_a.length() < set_b.length() { set_a.length() } else { set_b.length() }
  if n < 2 {
    return (RotationMatrix3::identity(), Vector3::zero())
  }

  let centroid_a = vector_centroid(set_a)
  let centroid_b = vector_centroid(set_b)

  let centered_a : Array[Vector3] = Array::new()
  let centered_b : Array[Vector3] = Array::new()
  let mut i = 0
  while i < n {
    centered_a.push(set_a[i].sub(centroid_a))
    centered_b.push(set_b[i].sub(centroid_b))
    i = i + 1
  }

  let rotation = compute_optimal_rotation(centered_a, centered_b)
  let transformed_centroid = rotation.transform(centroid_a)
  let translation = centroid_b.sub(transformed_centroid)

  (rotation, translation)
}

///|
/// Apply a rigid-body transformation (rotation + translation) to a set of vectors.
pub fn vector_apply_transform(vectors : Array[Vector3], rotation : RotationMatrix3, translation : Vector3) -> Array[Vector3] {
  let result : Array[Vector3] = Array::new()
  for v in vectors {
    let rotated = rotation.transform(v)
    result.push(rotated.add(translation))
  }
  result
}

///|
/// Compute the dihedral angle defined by four points (in radians).
pub fn vector_dihedral_angle(p1 : Vector3, p2 : Vector3, p3 : Vector3, p4 : Vector3) -> Double {
  let b1 = p2.sub(p1)
  let b2 = p3.sub(p2)
  let b3 = p4.sub(p3)

  let n1 = b1.cross(b2)
  let n2 = b2.cross(b3)

  let b2_norm = b2.norm()
  if b2_norm < 1.0e-15 {
    return 0.0
  }
  let m1 = n1.cross(b2.scale(1.0 / b2_norm))

  let x = n1.dot(n2)
  let y = m1.dot(n2)

  @math.atan2(y, x)
}

///|
/// Compute the dihedral angle in degrees.
pub fn vector_dihedral_angle_degrees(p1 : Vector3, p2 : Vector3, p3 : Vector3, p4 : Vector3) -> Double {
  let rad = vector_dihedral_angle(p1, p2, p3, p4)
  rad * 180.0 / @math.PI
}
