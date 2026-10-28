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
pub fn QCPAtomCoordinate::new(
  x : Double,
  y : Double,
  z : Double,
) -> QCPAtomCoordinate {
  QCPAtomCoordinate::{ x, y, z, }
}

///|
/// Convert to array.
pub fn QCPAtomCoordinate::to_array(self : QCPAtomCoordinate) -> Array[Double] {
  [self.x, self.y, self.z]
}

///|
/// Calculate distance between two atoms.
pub fn QCPAtomCoordinate::distance(
  self : QCPAtomCoordinate,
  other : QCPAtomCoordinate,
) -> Double {
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
pub fn QCPSuperimposer::set_rot(
  self : QCPSuperimposer,
  rot : Array[Array[Double]],
) -> QCPSuperimposer {
  QCPSuperimposer::{ ..self, rot, }
}

///|
/// Set the translation vector.
pub fn QCPSuperimposer::set_tran(
  self : QCPSuperimposer,
  tran : Array[Double],
) -> QCPSuperimposer {
  QCPSuperimposer::{ ..self, tran, }
}

///|
/// Set the RMSD value.
pub fn QCPSuperimposer::set_rmsd(
  self : QCPSuperimposer,
  rmsd : Double,
) -> QCPSuperimposer {
  QCPSuperimposer::{ ..self, rmsd, }
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
pub fn calculate_centroid(
  coords : Array[QCPAtomCoordinate],
) -> QCPAtomCoordinate {
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
pub fn center_coordinates(
  coords : Array[QCPAtomCoordinate],
) -> Array[QCPAtomCoordinate] {
  let centroid = calculate_centroid(coords)

  let result : Array[QCPAtomCoordinate] = Array(capacity=0)
  for coord in coords {
    result.push(
      QCPAtomCoordinate::new(
        coord.x - centroid.x,
        coord.y - centroid.y,
        coord.z - centroid.z,
      ),
    )
  }

  result
}

///|
/// Build the correlation matrix.
pub fn build_correlation_matrix(
  fixed : Array[QCPAtomCoordinate],
  moving : Array[QCPAtomCoordinate],
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

  [
    [q00, q01, q02, q03],
    [q10, q11, q12, q13],
    [q20, q21, q22, q23],
    [q30, q31, q32, q33],
  ]
}

///|
/// Compute the eigenvector for the largest eigenvalue of a symmetric 4x4
/// matrix using Jacobi rotations.
fn qcp_largest_eigenvector(q : Array[Array[Double]]) -> Array[Double] {
  let matrix : Array[Array[Double]] = []
  let vectors : Array[Array[Double]] = []
  let mut row = 0
  while row < 4 {
    matrix.push(q[row].copy())
    let vector_row = Array::make(4, 0.0)
    vector_row[row] = 1.0
    vectors.push(vector_row)
    row = row + 1
  }
  let mut sweep = 0
  while sweep < 64 {
    let mut pivot_row = 0
    let mut pivot_column = 1
    let mut largest = matrix[0][1].abs()
    let mut current_row = 0
    while current_row < 3 {
      let mut current_column = current_row + 1
      while current_column < 4 {
        let magnitude = matrix[current_row][current_column].abs()
        if magnitude > largest {
          largest = magnitude
          pivot_row = current_row
          pivot_column = current_column
        }
        current_column = current_column + 1
      }
      current_row = current_row + 1
    }
    if largest < 1.0e-12 {
      break
    }
    let diagonal_difference = matrix[pivot_column][pivot_column] -
      matrix[pivot_row][pivot_row]
    let ratio = diagonal_difference / (2.0 * matrix[pivot_row][pivot_column])
    let tangent = if ratio >= 0.0 {
      1.0 / (ratio + (1.0 + ratio * ratio).sqrt())
    } else {
      -1.0 / (-ratio + (1.0 + ratio * ratio).sqrt())
    }
    let cosine = 1.0 / (1.0 + tangent * tangent).sqrt()
    let sine = tangent * cosine
    let pivot_value = matrix[pivot_row][pivot_column]
    let row_diagonal = matrix[pivot_row][pivot_row]
    let column_diagonal = matrix[pivot_column][pivot_column]
    matrix[pivot_row][pivot_row] = row_diagonal - tangent * pivot_value
    matrix[pivot_column][pivot_column] = column_diagonal + tangent * pivot_value
    matrix[pivot_row][pivot_column] = 0.0
    matrix[pivot_column][pivot_row] = 0.0
    let mut index = 0
    while index < 4 {
      if index != pivot_row && index != pivot_column {
        let row_value = matrix[index][pivot_row]
        let column_value = matrix[index][pivot_column]
        let rotated_row = cosine * row_value - sine * column_value
        let rotated_column = sine * row_value + cosine * column_value
        matrix[index][pivot_row] = rotated_row
        matrix[pivot_row][index] = rotated_row
        matrix[index][pivot_column] = rotated_column
        matrix[pivot_column][index] = rotated_column
      }
      let vector_row = vectors[index][pivot_row]
      let vector_column = vectors[index][pivot_column]
      vectors[index][pivot_row] = cosine * vector_row - sine * vector_column
      vectors[index][pivot_column] = sine * vector_row + cosine * vector_column
      index = index + 1
    }
    sweep = sweep + 1
  }
  let mut largest_index = 0
  let mut index = 1
  while index < 4 {
    if matrix[index][index] > matrix[largest_index][largest_index] {
      largest_index = index
    }
    index = index + 1
  }
  [
    vectors[0][largest_index],
    vectors[1][largest_index],
    vectors[2][largest_index],
    vectors[3][largest_index],
  ]
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
  tran : Array[Double],
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
pub fn QCPSuperimposer::set(
  self : QCPSuperimposer,
  fixed : Array[QCPAtomCoordinate],
  moving : Array[QCPAtomCoordinate],
) -> QCPSuperimposer {
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

  let c = build_correlation_matrix(moving_centered, fixed_centered)
  let q = build_qcp_matrix(c)
  let eigenvec = qcp_largest_eigenvector(q)

  let rot = quaternion_to_rotation(eigenvec)

  let t0 = rot[0][0] * centroid_moving.x +
    rot[0][1] * centroid_moving.y +
    rot[0][2] * centroid_moving.z
  let t1 = rot[1][0] * centroid_moving.x +
    rot[1][1] * centroid_moving.y +
    rot[1][2] * centroid_moving.z
  let t2 = rot[2][0] * centroid_moving.x +
    rot[2][1] * centroid_moving.y +
    rot[2][2] * centroid_moving.z

  let tran = [
    centroid_fixed.x - t0,
    centroid_fixed.y - t1,
    centroid_fixed.z - t2,
  ]

  let rmsd = qcp_calculate_rmsd(fixed, moving, rot, tran)

  QCPSuperimposer::{ rot, tran, rmsd, }
}

///|
/// Apply rotation and translation to coordinates.
pub fn QCPSuperimposer::apply(
  self : QCPSuperimposer,
  coords : Array[QCPAtomCoordinate],
) -> Array[QCPAtomCoordinate] {
  let result : Array[QCPAtomCoordinate] = Array(capacity=0)

  let mut i = 0
  while i < coords.length() {
    let coord = coords[i]
    let x = self.rot[0][0] * coord.x +
      self.rot[0][1] * coord.y +
      self.rot[0][2] * coord.z +
      self.tran[0]
    let y = self.rot[1][0] * coord.x +
      self.rot[1][1] * coord.y +
      self.rot[1][2] * coord.z +
      self.tran[1]
    let z = self.rot[2][0] * coord.x +
      self.rot[2][1] * coord.y +
      self.rot[2][2] * coord.z +
      self.tran[2]
    result.push(QCPAtomCoordinate::new(x, y, z))
    i = i + 1
  }

  result
}

///|
/// Bio API functions.
pub fn bio_qcp_superimpose(
  fixed : Array[QCPAtomCoordinate],
  moving : Array[QCPAtomCoordinate],
) -> QCPSuperimposer {
  QCPSuperimposer::new().set(fixed, moving)
}

///|
pub fn bio_qcp_calculate_rmsd(
  fixed : Array[QCPAtomCoordinate],
  moving : Array[QCPAtomCoordinate],
) -> Double {
  let superimposer = bio_qcp_superimpose(fixed, moving)
  superimposer.get_rmsd()
}

///|
pub fn bio_qcp_apply_transform(
  coords : Array[QCPAtomCoordinate],
  rot : Array[Array[Double]],
  tran : Array[Double],
) -> Array[QCPAtomCoordinate] {
  let superimposer = QCPSuperimposer::new().set_rot(rot).set_tran(tran)
  superimposer.apply(coords)
}
