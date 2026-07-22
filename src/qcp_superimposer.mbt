///|
/// QCP Superimposer module - Protein structure superposition using Quaternion
/// Characteristic Polynomial algorithm.
/// Inspired by Biopython Bio.PDB.QCPSuperimposer.
///
/// The QCP algorithm provides more accurate and numerically stable structure
/// alignment compared to SVD-based methods.

///|
/// Represents 3D coordinates of an atom.
pub struct QCPAtomCoordinate {
  x : Double
  y : Double
  z : Double
}

///|
/// Create a new QCPAtomCoordinate.
pub fn QCPAtomCoordinate::new(x : Double, y : Double, z : Double) -> QCPAtomCoordinate {
  QCPAtomCoordinate::{ x, y, z }
}

///|
/// Convert to array.
pub fn QCPAtomCoordinate::to_array(self : QCPAtomCoordinate) -> Array[Double] {
  [self.x, self.y, self.z]
}

///|
/// Calculate distance between two atoms.
pub fn QCPAtomCoordinate::distance(self : QCPAtomCoordinate, other : QCPAtomCoordinate) -> Double {
  let dx = self.x - other.x
  let dy = self.y - other.y
  let dz = self.z - other.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// QCP Superimposer struct.
pub struct QCPSuperimposer {
  rot : Array[Array[Double]]
  tran : Array[Double]
  rmsd : Double
}

///|
/// Create a new QCPSuperimposer.
pub fn QCPSuperimposer::new() -> QCPSuperimposer {
  QCPSuperimposer::{
    rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    tran: [0.0, 0.0, 0.0],
    rmsd: 0.0,
  }
}

///|
/// Set the rotation matrix.
pub fn QCPSuperimposer::set_rot(self : QCPSuperimposer, rot : Array[Array[Double]]) -> QCPSuperimposer {
  QCPSuperimposer::{ ..self, rot }
}

///|
/// Set the translation vector.
pub fn QCPSuperimposer::set_tran(self : QCPSuperimposer, tran : Array[Double]) -> QCPSuperimposer {
  QCPSuperimposer::{ ..self, tran }
}

///|
/// Set the RMSD value.
pub fn QCPSuperimposer::set_rmsd(self : QCPSuperimposer, rmsd : Double) -> QCPSuperimposer {
  QCPSuperimposer::{ ..self, rmsd }
}

///|
/// Get the rotation matrix.
pub fn QCPSuperimposer::get_rot(self : QCPSuperimposer) -> Array[Array[Double]] {
  self.rot
}

///|
/// Get the translation vector.
pub fn QCPSuperimposer::get_tran(self : QCPSuperimposer) -> Array[Double] {
  self.tran
}

///|
/// Get the RMSD value.
pub fn QCPSuperimposer::get_rmsd(self : QCPSuperimposer) -> Double {
  self.rmsd
}

///|
/// Calculate the centroid of a set of coordinates.
pub fn calculate_centroid(coords : Array[QCPAtomCoordinate]) -> QCPAtomCoordinate {
  if coords.length() == 0 {
    return QCPAtomCoordinate::new(0.0, 0.0, 0.0)
  }
  
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_z = 0.0
  
  for coord in coords {
    sum_x = sum_x + coord.x
    sum_y = sum_y + coord.y
    sum_z = sum_z + coord.z
  }
  
  let n = coords.length().to_double()
  QCPAtomCoordinate::new(sum_x / n, sum_y / n, sum_z / n)
}

///|
/// Center coordinates around centroid.
pub fn center_coordinates(coords : Array[QCPAtomCoordinate]) -> Array[QCPAtomCoordinate] {
  let centroid = calculate_centroid(coords)
  
  let result : Array[QCPAtomCoordinate] = Array::new()
  for coord in coords {
    result.push(QCPAtomCoordinate::new(
      coord.x - centroid.x,
      coord.y - centroid.y,
      coord.z - centroid.z
    ))
  }
  
  result
}

///|
/// Build the correlation matrix.
pub fn build_correlation_matrix(
  fixed : Array[QCPAtomCoordinate],
  moving : Array[QCPAtomCoordinate]
) -> Array[Array[Double]] {
  let mut c00 = 0.0
  let mut c01 = 0.0
  let mut c02 = 0.0
  let mut c10 = 0.0
  let mut c11 = 0.0
  let mut c12 = 0.0
  let mut c20 = 0.0
  let mut c21 = 0.0
  let mut c22 = 0.0
  
  let len = fixed.length()
  let mut i = 0
  while i < len {
    c00 = c00 + fixed[i].x * moving[i].x
    c01 = c01 + fixed[i].x * moving[i].y
    c02 = c02 + fixed[i].x * moving[i].z
    c10 = c10 + fixed[i].y * moving[i].x
    c11 = c11 + fixed[i].y * moving[i].y
    c12 = c12 + fixed[i].y * moving[i].z
    c20 = c20 + fixed[i].z * moving[i].x
    c21 = c21 + fixed[i].z * moving[i].y
    c22 = c22 + fixed[i].z * moving[i].z
    i = i + 1
  }
  
  [[c00, c01, c02], [c10, c11, c12], [c20, c21, c22]]
}

///|
/// Calculate the trace of a 3x3 matrix.
pub fn matrix_trace(m : Array[Array[Double]]) -> Double {
  m[0][0] + m[1][1] + m[2][2]
}

///|
/// Calculate the determinant of a 3x3 matrix.
pub fn matrix_determinant_3x3(m : Array[Array[Double]]) -> Double {
  let term1 = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
  let term2 = m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
  let term3 = m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
  term1 - term2 + term3
}

///|
/// Build the 4x4 QCP matrix.
pub fn build_qcp_matrix(c : Array[Array[Double]]) -> Array[Array[Double]] {
  let q00 = c[0][0] + c[1][1] + c[2][2]
  let q01 = c[1][2] - c[2][1]
  let q02 = c[2][0] - c[0][2]
  let q03 = c[0][1] - c[1][0]
  
  let q10 = q01
  let q11 = c[0][0] - c[1][1] - c[2][2]
  let q12 = c[0][1] + c[1][0]
  let q13 = c[0][2] + c[2][0]
  
  let q20 = q02
  let q21 = q12
  let q22 = -c[0][0] + c[1][1] - c[2][2]
  let q23 = c[1][2] + c[2][1]
  
  let q30 = q03
  let q31 = q13
  let q32 = q23
  let q33 = -c[0][0] - c[1][1] + c[2][2]
  
  [[q00, q01, q02, q03], [q10, q11, q12, q13], [q20, q21, q22, q23], [q30, q31, q32, q33]]
}

///|
/// Solve the quartic equation and find the maximum real root.
pub fn solve_quartic(q : Array[Array[Double]]) -> Double {
  let a = -q[0][0] - q[1][1] - q[2][2] - q[3][3]
  
  let b1 = q[0][0] * q[1][1] + q[0][0] * q[2][2] + q[0][0] * q[3][3]
  let b2 = q[1][1] * q[2][2] + q[1][1] * q[3][3] + q[2][2] * q[3][3]
  let b3 = q[0][1] * q[1][0] + q[0][2] * q[2][0] + q[0][3] * q[3][0]
  let b4 = q[1][2] * q[2][1] + q[1][3] * q[3][1] + q[2][3] * q[3][2]
  let b = b1 + b2 - b3 - b4
  
  let c1 = q[0][0] * q[1][1] * q[2][2] + q[0][0] * q[1][1] * q[3][3] + q[0][0] * q[2][2] * q[3][3] + q[1][1] * q[2][2] * q[3][3]
  let c2 = q[0][0] * q[1][2] * q[2][1] + q[0][0] * q[1][3] * q[3][1] + q[0][0] * q[2][3] * q[3][2]
  let c3 = q[1][1] * q[0][2] * q[2][0] + q[1][1] * q[0][3] * q[3][0] + q[1][1] * q[2][3] * q[3][2]
  let c4 = q[2][2] * q[0][1] * q[1][0] + q[2][2] * q[0][3] * q[3][0] + q[2][2] * q[1][3] * q[3][1]
  let c5 = q[3][3] * q[0][1] * q[1][0] + q[3][3] * q[0][2] * q[2][0] + q[3][3] * q[1][2] * q[2][1]
  let c6 = q[0][1] * q[1][2] * q[2][0] + q[0][1] * q[1][3] * q[3][0] + q[0][2] * q[2][3] * q[3][0]
  let c7 = q[1][2] * q[2][3] * q[3][1] + q[0][1] * q[2][1] * q[1][0] + q[0][2] * q[3][2] * q[2][0] + q[1][3] * q[3][1] * q[1][0]
  let c = -(c1) + c2 + c3 + c4 + c5 - c6 - c7
  
  let d1 = q[0][0] * q[1][1] * q[2][2] * q[3][3]
  let d2 = q[0][1] * q[1][2] * q[2][3] * q[3][0] + q[0][2] * q[1][3] * q[2][0] * q[3][1] + q[0][3] * q[1][0] * q[2][1] * q[3][2]
  let d3 = q[0][3] * q[1][2] * q[2][1] * q[3][0] + q[0][2] * q[1][1] * q[2][3] * q[3][0] + q[0][1] * q[1][3] * q[2][0] * q[3][1]
  let d4 = q[0][0] * q[1][2] * q[2][3] * q[3][1] + q[0][3] * q[1][1] * q[2][0] * q[3][2] + q[0][2] * q[1][0] * q[2][1] * q[3][2]
  let d5 = q[0][1] * q[1][0] * q[2][3] * q[3][2] + q[0][0] * q[1][3] * q[2][1] * q[3][2]
  let d = d1 + d2 - d3 - d4 - d5
  
  let mut max_root = 0.0
  
  for i in 0..<100 {
    let t = i.to_double() * 0.02 - 1.0
    let val = t * t * t * t + a * t * t * t + b * t * t + c * t + d
    if val > 0.0 && t > max_root {
      max_root = t
    }
  }
  
  max_root
}

///|
/// Compute the eigenvector corresponding to the maximum eigenvalue.
pub fn compute_eigenvector(q : Array[Array[Double]], lambda : Double) -> Array[Double] {
  let m00 = q[0][0] - lambda
  let m01 = q[0][1]
  let m02 = q[0][2]
  let m03 = q[0][3]
  
  let m10 = q[1][0]
  let m11 = q[1][1] - lambda
  let m12 = q[1][2]
  let m13 = q[1][3]
  
  let m20 = q[2][0]
  let m21 = q[2][1]
  let m22 = q[2][2] - lambda
  let m23 = q[2][3]
  
  let m30 = q[3][0]
  let m31 = q[3][1]
  let m32 = q[3][2]
  let m33 = q[3][3] - lambda
  
  let mut v0 = 1.0
  let mut v1 = 0.0
  let mut v2 = 0.0
  let mut v3 = 0.0
  
  let mut iter = 0
  while iter < 50 {
    let nv0 = m00 * v0 + m01 * v1 + m02 * v2 + m03 * v3
    let nv1 = m10 * v0 + m11 * v1 + m12 * v2 + m13 * v3
    let nv2 = m20 * v0 + m21 * v1 + m22 * v2 + m23 * v3
    let nv3 = m30 * v0 + m31 * v1 + m32 * v2 + m33 * v3
    
    let norm = (nv0 * nv0 + nv1 * nv1 + nv2 * nv2 + nv3 * nv3).sqrt()
    
    if norm > 0.000001 {
      v0 = nv0 / norm
      v1 = nv1 / norm
      v2 = nv2 / norm
      v3 = nv3 / norm
    }
    
    iter = iter + 1
  }
  
  [v0, v1, v2, v3]
}

///|
/// Convert quaternion to rotation matrix.
pub fn quaternion_to_rotation(q : Array[Double]) -> Array[Array[Double]] {
  let w = q[0]
  let x = q[1]
  let y = q[2]
  let z = q[3]
  
  let ww = w * w
  let xx = x * x
  let yy = y * y
  let zz = z * z
  let wx = w * x
  let wy = w * y
  let wz = w * z
  let xy = x * y
  let xz = x * z
  let yz = y * z
  
  let rot = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  
  rot[0][0] = ww + xx - yy - zz
  rot[0][1] = 2.0 * (xy - wz)
  rot[0][2] = 2.0 * (xz + wy)
  
  rot[1][0] = 2.0 * (xy + wz)
  rot[1][1] = ww - xx + yy - zz
  rot[1][2] = 2.0 * (yz - wx)
  
  rot[2][0] = 2.0 * (xz - wy)
  rot[2][1] = 2.0 * (yz + wx)
  rot[2][2] = ww - xx - yy + zz
  
  rot
}

///|
/// Calculate RMSD between two sets of coordinates.
pub fn qcp_calculate_rmsd(
  fixed : Array[QCPAtomCoordinate],
  moving : Array[QCPAtomCoordinate],
  rot : Array[Array[Double]],
  tran : Array[Double]
) -> Double {
  let mut sum_sq = 0.0
  let len = fixed.length()
  
  let mut i = 0
  while i < len {
    let mx = moving[i].x
    let my = moving[i].y
    let mz = moving[i].z
    
    let rx = rot[0][0] * mx + rot[0][1] * my + rot[0][2] * mz + tran[0]
    let ry = rot[1][0] * mx + rot[1][1] * my + rot[1][2] * mz + tran[1]
    let rz = rot[2][0] * mx + rot[2][1] * my + rot[2][2] * mz + tran[2]
    
    let dx = fixed[i].x - rx
    let dy = fixed[i].y - ry
    let dz = fixed[i].z - rz
    
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
    i = i + 1
  }
  
  (sum_sq / len.to_double()).sqrt()
}

///|
/// Set the coordinates for superposition.
/// fixed: reference structure coordinates
/// moving: structure to be aligned coordinates
pub fn QCPSuperimposer::set(self : QCPSuperimposer, fixed : Array[QCPAtomCoordinate], moving : Array[QCPAtomCoordinate]) -> QCPSuperimposer {
  if fixed.length() != moving.length() {
    return self
  }
  
  if fixed.length() == 0 {
    return self.set_rmsd(0.0)
  }
  
  let fixed_centered = center_coordinates(fixed)
  let moving_centered = center_coordinates(moving)
  
  let centroid_fixed = calculate_centroid(fixed)
  let centroid_moving = calculate_centroid(moving)
  
  let c = build_correlation_matrix(fixed_centered, moving_centered)
  let q = build_qcp_matrix(c)
  
  let lambda_max = solve_quartic(q)
  let eigenvec = compute_eigenvector(q, lambda_max)
  
  let rot = quaternion_to_rotation(eigenvec)
  
  let t0 = rot[0][0] * centroid_moving.x + rot[0][1] * centroid_moving.y + rot[0][2] * centroid_moving.z
  let t1 = rot[1][0] * centroid_moving.x + rot[1][1] * centroid_moving.y + rot[1][2] * centroid_moving.z
  let t2 = rot[2][0] * centroid_moving.x + rot[2][1] * centroid_moving.y + rot[2][2] * centroid_moving.z
  
  let tran = [centroid_fixed.x - t0, centroid_fixed.y - t1, centroid_fixed.z - t2]
  
  let rmsd = qcp_calculate_rmsd(fixed, moving, rot, tran)
  
  QCPSuperimposer::{ rot, tran, rmsd }
}

///|
/// Apply rotation and translation to coordinates.
pub fn QCPSuperimposer::apply(self : QCPSuperimposer, coords : Array[QCPAtomCoordinate]) -> Array[QCPAtomCoordinate] {
  let result : Array[QCPAtomCoordinate] = Array::new()
  
  let mut i = 0
  while i < coords.length() {
    let coord = coords[i]
    let x = self.rot[0][0] * coord.x + self.rot[0][1] * coord.y + self.rot[0][2] * coord.z + self.tran[0]
    let y = self.rot[1][0] * coord.x + self.rot[1][1] * coord.y + self.rot[1][2] * coord.z + self.tran[1]
    let z = self.rot[2][0] * coord.x + self.rot[2][1] * coord.y + self.rot[2][2] * coord.z + self.tran[2]
    result.push(QCPAtomCoordinate::new(x, y, z))
    i = i + 1
  }
  
  result
}

///|
/// Bio API functions.

pub fn bio_qcp_superimpose(fixed : Array[QCPAtomCoordinate], moving : Array[QCPAtomCoordinate]) -> QCPSuperimposer {
  QCPSuperimposer::new().set(fixed, moving)
}

pub fn bio_qcp_calculate_rmsd(fixed : Array[QCPAtomCoordinate], moving : Array[QCPAtomCoordinate]) -> Double {
  let superimposer = bio_qcp_superimpose(fixed, moving)
  superimposer.get_rmsd()
}

pub fn bio_qcp_apply_transform(coords : Array[QCPAtomCoordinate], rot : Array[Array[Double]], tran : Array[Double]) -> Array[QCPAtomCoordinate] {
  let superimposer = QCPSuperimposer::new().set_rot(rot).set_tran(tran)
  superimposer.apply(coords)
}