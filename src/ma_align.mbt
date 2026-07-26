///|
/// Bio.PDB.MAalign - Multiple Protein Structure Alignment.
///
/// Provides functions for multiple protein structure alignment using
/// the Kabsch algorithm for optimal superposition and iterative
/// refinement to align multiple structures simultaneously.
///
/// References:
/// - Kabsch W (1978). "A discussion of the solution for the best
///   rotation to fit two sets of vectors." Acta Cryst. A34:827-828.
/// - Kabsch W (1978). "The solution for the best rotation to fit
///   two sets of vectors." Acta Cryst. A34:828.

///|
/// Represents a protein structure for multiple alignment.
pub struct StructureData {
  structure_id : String
  ca_coordinates : Array[(Double, Double, Double)]
  sequence : String
  residue_names : Array[String]
} derive(Eq)

///|
/// Create a new StructureData.
pub fn StructureData::new(
  structure_id : String,
  ca_coordinates : Array[(Double, Double, Double)],
  sequence : String,
  residue_names : Array[String],
) -> StructureData {
  StructureData::{ structure_id, ca_coordinates, sequence, residue_names }
}

///|
/// Get the structure ID.
pub fn StructureData::get_id(self : StructureData) -> String {
  self.structure_id
}

///|
/// Get the CA atom coordinates.
pub fn StructureData::get_coordinates(self : StructureData) -> Array[(Double, Double, Double)] {
  self.ca_coordinates
}

///|
/// Get the amino acid sequence.
pub fn StructureData::get_sequence(self : StructureData) -> String {
  self.sequence
}

///|
/// Get the residue names.
pub fn StructureData::get_residue_names(self : StructureData) -> Array[String] {
  self.residue_names
}

///|
/// Get the number of residues (CA atoms).
pub fn StructureData::length(self : StructureData) -> Int {
  self.ca_coordinates.length()
}

///|
/// Result of a multiple structure alignment.
pub struct AlignmentResult {
  structure_ids : Array[String]
  aligned_sequences : Array[String]
  rmsd_matrix : Array[Array[Double]]
  aligned_coordinates : Array[Array[(Double, Double, Double)]]
  conservation : Array[Double]
  total_rmsd : Double
}

///|
/// Create a new AlignmentResult.
pub fn AlignmentResult::new(
  structure_ids : Array[String],
  aligned_sequences : Array[String],
  rmsd_matrix : Array[Array[Double]],
  aligned_coordinates : Array[Array[(Double, Double, Double)]],
  conservation : Array[Double],
  total_rmsd : Double,
) -> AlignmentResult {
  AlignmentResult::{
    structure_ids,
    aligned_sequences,
    rmsd_matrix,
    aligned_coordinates,
    conservation,
    total_rmsd,
  }
}

///|
/// Get the structure IDs.
pub fn AlignmentResult::get_structure_ids(self : AlignmentResult) -> Array[String] {
  self.structure_ids
}

///|
/// Get the aligned sequences.
pub fn AlignmentResult::get_aligned_sequences(self : AlignmentResult) -> Array[String] {
  self.aligned_sequences
}

///|
/// Get the pairwise RMSD matrix.
pub fn AlignmentResult::get_rmsd_matrix(self : AlignmentResult) -> Array[Array[Double]] {
  self.rmsd_matrix
}

///|
/// Get the aligned coordinates.
pub fn AlignmentResult::get_aligned_coordinates(self : AlignmentResult) -> Array[Array[(Double, Double, Double)]] {
  self.aligned_coordinates
}

///|
/// Get the per-column conservation scores.
pub fn AlignmentResult::get_conservation(self : AlignmentResult) -> Array[Double] {
  self.conservation
}

///|
/// Get the total RMSD.
pub fn AlignmentResult::get_total_rmsd(self : AlignmentResult) -> Double {
  self.total_rmsd
}

///|
/// Get the alignment length (number of columns).
pub fn AlignmentResult::alignment_length(self : AlignmentResult) -> Int {
  if self.aligned_sequences.length() > 0 {
    self.aligned_sequences[0].length()
  } else {
    0
  }
}

///|
/// MAaligner - Multiple structure aligner with configurable parameters.
pub struct MAaligner {
  gap_opening : Double
  gap_extension : Double
  max_iterations : Int
  rmsd_threshold : Double
}

///|
/// Create a new MAaligner with default parameters.
pub fn ma_aligner_new() -> MAaligner {
  MAaligner::{
    gap_opening: -10.0,
    gap_extension: -1.0,
    max_iterations: 10,
    rmsd_threshold: 2.0,
  }
}

///|
/// Create a new MAaligner with custom parameters.
///
/// `gap_open` — gap opening penalty (typically negative).
/// `gap_ext` — gap extension penalty (typically negative).
/// `max_iter` — maximum number of iterations for refinement.
/// `rmsd_thresh` — RMSD convergence threshold in Angstroms.
pub fn ma_aligner_with_params(
  gap_open : Double,
  gap_ext : Double,
  max_iter : Int,
  rmsd_thresh : Double,
) -> MAaligner {
  MAaligner::{
    gap_opening: gap_open,
    gap_extension: gap_ext,
    max_iterations: max_iter,
    rmsd_threshold: rmsd_thresh,
  }
}

///|
/// Get the gap opening penalty.
pub fn MAaligner::get_gap_opening(self : MAaligner) -> Double {
  self.gap_opening
}

///|
/// Get the gap extension penalty.
pub fn MAaligner::get_gap_extension(self : MAaligner) -> Double {
  self.gap_extension
}

///|
/// Get the maximum number of iterations.
pub fn MAaligner::get_max_iterations(self : MAaligner) -> Int {
  self.max_iterations
}

///|
/// Get the RMSD convergence threshold.
pub fn MAaligner::get_rmsd_threshold(self : MAaligner) -> Double {
  self.rmsd_threshold
}

// ============================================================================
// Internal helper functions (not part of public API)
// ============================================================================

///|
/// Compute the centroid of a set of coordinates.
fn ma_compute_centroid(coords : Array[(Double, Double, Double)]) -> (Double, Double, Double) {
  let n = coords.length()
  if n == 0 {
    return (0.0, 0.0, 0.0)
  }
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sz = 0.0
  let mut i = 0
  while i < n {
    let (x, y, z) = coords[i]
    sx = sx + x
    sy = sy + y
    sz = sz + z
    i = i + 1
  }
  (sx / n.to_double(), sy / n.to_double(), sz / n.to_double())
}

///|
/// Compute the Euclidean distance between two 3D points.
fn ma_point_distance(p1 : (Double, Double, Double), p2 : (Double, Double, Double)) -> Double {
  let dx = p1.0 - p2.0
  let dy = p1.1 - p2.1
  let dz = p1.2 - p2.2
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// Compute the cross-covariance matrix H = P * Q^T for centered coordinates.
fn ma_compute_covariance(
  p : Array[(Double, Double, Double)],
  q : Array[(Double, Double, Double)],
) -> Array[Array[Double]] {
  let h : Array[Array[Double]] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let n = p.length()
  let mut i = 0
  while i < n {
    let (px, py, pz) = p[i]
    let (qx, qy, qz) = q[i]
    h[0][0] = h[0][0] + px * qx
    h[0][1] = h[0][1] + px * qy
    h[0][2] = h[0][2] + px * qz
    h[1][0] = h[1][0] + py * qx
    h[1][1] = h[1][1] + py * qy
    h[1][2] = h[1][2] + py * qz
    h[2][0] = h[2][0] + pz * qx
    h[2][1] = h[2][1] + pz * qy
    h[2][2] = h[2][2] + pz * qz
    i = i + 1
  }
  h
}

///|
/// Multiply two 3x3 matrices.
fn ma_matmul_3x3(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let mut i = 0
  while i < 3 {
    let mut j = 0
    while j < 3 {
      let mut sum = 0.0
      let mut k = 0
      while k < 3 {
        sum = sum + a[i][k] * b[k][j]
        k = k + 1
      }
      result[i][j] = sum
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Transpose a 3x3 matrix.
fn ma_transpose_3x3(m : Array[Array[Double]]) -> Array[Array[Double]] {
  [[m[0][0], m[1][0], m[2][0]],
   [m[0][1], m[1][1], m[2][1]],
   [m[0][2], m[1][2], m[2][2]]]
}

///|
/// Compute the determinant of a 3x3 matrix.
fn ma_det_3x3(m : Array[Array[Double]]) -> Double {
  let a = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
  let b = m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
  let c = m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
  a - b + c
}

///|
/// Compute the inverse of a 3x3 matrix using the adjugate method.
fn ma_invert_3x3(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let det = ma_det_3x3(m)
  if det.abs() < 1.0e-15 {
    return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  }
  let inv_det = 1.0 / det
  let result : Array[Array[Double]] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  result[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det
  result[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det
  result[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det
  result[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det
  result[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det
  result[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det
  result[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det
  result[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det
  result[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det
  result
}

///|
/// Compute eigenvalues and eigenvectors of a 3x3 symmetric matrix.
/// Uses the analytical solution for 3x3 symmetric matrices.
/// Returns (eigenvalues, eigenvectors_as_columns).
fn ma_eigen_3x3(m : Array[Array[Double]]) -> (Array[Double], Array[Array[Double]]) {
  let eps = 1.0e-12
  let a = m

  let trace = a[0][0] + a[1][1] + a[2][2]
  let q = trace / 3.0

  let b00 = a[0][0] - q
  let b11 = a[1][1] - q
  let b22 = a[2][2] - q
  let b01 = a[0][1]
  let b02 = a[0][2]
  let b12 = a[1][2]

  let p1 = (b01 * b01 + b02 * b02 + b12 * b12).sqrt()

  let eigenvalues : Array[Double] = [0.0, 0.0, 0.0]
  let eigenvectors : Array[Array[Double]] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

  if p1.abs() < eps {
    eigenvalues[0] = a[0][0]
    eigenvalues[1] = a[1][1]
    eigenvalues[2] = a[2][2]
    eigenvectors[0][0] = 1.0; eigenvectors[0][1] = 0.0; eigenvectors[0][2] = 0.0
    eigenvectors[1][0] = 0.0; eigenvectors[1][1] = 1.0; eigenvectors[1][2] = 0.0
    eigenvectors[2][0] = 0.0; eigenvectors[2][1] = 0.0; eigenvectors[2][2] = 1.0
    return (eigenvalues, eigenvectors)
  }

  let p2 = b00 * b00 + b01 * b01 + b02 * b02
  let p = p2.sqrt()

  let det_b_a = b00 * (b11 * b22 - b12 * b12)
  let det_b_b = b01 * (b01 * b22 - b12 * b02)
  let det_b_c = b02 * (b01 * b12 - b11 * b02)
  let det_b = det_b_a - det_b_b + det_b_c

  let r = det_b / (p1 * p2 + eps)

  let phi = if r <= -1.0 {
    @math.PI / 3.0
  } else if r >= 1.0 {
    0.0
  } else {
    @math.acos(r) / 3.0
  }

  eigenvalues[0] = q + 2.0 * p * @math.cos(phi)
  eigenvalues[2] = q + 2.0 * p * @math.cos(phi + 2.0 * @math.PI / 3.0)
  eigenvalues[1] = 3.0 * q - eigenvalues[0] - eigenvalues[2]

  let sorted_eigs : Array[Double] = [eigenvalues[0], eigenvalues[1], eigenvalues[2]]

  let v0_norm_sq = b01 * b01 + b02 * b02
  if v0_norm_sq < eps {
    eigenvectors[0][0] = 1.0; eigenvectors[0][1] = 0.0; eigenvectors[0][2] = 0.0
  } else {
    let v0n = v0_norm_sq.sqrt()
    eigenvectors[0][0] = b01 / v0n
    eigenvectors[0][1] = b02 / v0n
    eigenvectors[0][2] = 0.0
  }

  let l1 = eigenvalues[1]
  let e1 = (l1 - a[0][0]) / (l1 - eigenvalues[0] + eps)
  let e2 = (l1 - a[1][1]) / (l1 - eigenvalues[0] + eps)
  let e3 = (l1 - a[2][2]) / (l1 - eigenvalues[0] + eps)
  let e_norm = (e1 * e1 + e2 * e2 + e3 * e3).sqrt()
  if e_norm < eps {
    eigenvectors[1][0] = 0.0; eigenvectors[1][1] = 1.0; eigenvectors[1][2] = 0.0
  } else {
    eigenvectors[1][0] = e1 / e_norm
    eigenvectors[1][1] = e2 / e_norm
    eigenvectors[1][2] = e3 / e_norm
  }

  let cross_x = eigenvectors[0][1] * eigenvectors[1][2] - eigenvectors[0][2] * eigenvectors[1][1]
  let cross_y = eigenvectors[0][2] * eigenvectors[1][0] - eigenvectors[0][0] * eigenvectors[1][2]
  let cross_z = eigenvectors[0][0] * eigenvectors[1][1] - eigenvectors[0][1] * eigenvectors[1][0]
  let cross_norm = (cross_x * cross_x + cross_y * cross_y + cross_z * cross_z).sqrt()
  if cross_norm < eps {
    eigenvectors[2][0] = 0.0; eigenvectors[2][1] = 0.0; eigenvectors[2][2] = 1.0
  } else {
    eigenvectors[2][0] = cross_x / cross_norm
    eigenvectors[2][1] = cross_y / cross_norm
    eigenvectors[2][2] = cross_z / cross_norm
  }

  (sorted_eigs, eigenvectors)
}

///|
/// Compute the square root of a 3x3 symmetric positive-definite matrix
/// using eigendecomposition: sqrt(A) = V * diag(sqrt(lambda)) * V^T.
fn ma_matrix_sqrt_3x3(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let (eigenvalues, eigenvectors) = ma_eigen_3x3(m)

  let sqrt_eig0 = if eigenvalues[0] >= 0.0 { eigenvalues[0].sqrt() } else { 0.0 }
  let sqrt_eig1 = if eigenvalues[1] >= 0.0 { eigenvalues[1].sqrt() } else { 0.0 }
  let sqrt_eig2 = if eigenvalues[2] >= 0.0 { eigenvalues[2].sqrt() } else { 0.0 }

  let result : Array[Array[Double]] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let mut i = 0
  while i < 3 {
    let mut j = 0
    while j < 3 {
      let v0 = eigenvectors[i][0] * sqrt_eig0 * eigenvectors[j][0]
      let v1 = eigenvectors[i][1] * sqrt_eig1 * eigenvectors[j][1]
      let v2 = eigenvectors[i][2] * sqrt_eig2 * eigenvectors[j][2]
      let val = v0 + v1 + v2
      result[i][j] = val
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Build a 3x3 rotation matrix from three rotation angles (around x, y, z axes).
fn ma_build_rotation_matrix(angles : (Double, Double, Double)) -> Array[Array[Double]] {
  let rx = angles.0
  let ry = angles.1
  let rz = angles.2

  let cx = @math.cos(rx)
  let sx = @math.sin(rx)
  let cy = @math.cos(ry)
  let sy = @math.sin(ry)
  let cz = @math.cos(rz)
  let sz = @math.sin(rz)

  [[cy * cz, cx * sz + sx * sy * cz, sx * sz - cx * sy * cz],
   [-cy * sz, cx * cz - sx * sy * sz, sx * cz + cx * sy * sz],
   [sy, -sx * cy, cx * cy]]
}

///|
/// Extract Euler angles (rotation around x, y, z) from a rotation matrix.
fn ma_rotation_to_angles(rot : Array[Array[Double]]) -> (Double, Double, Double) {
  let angle_x = @math.atan2(rot[2][1], rot[2][2])
  let angle_y = @math.atan2(-rot[2][0], (rot[2][1] * rot[2][1] + rot[2][2] * rot[2][2]).sqrt())
  let angle_z = @math.atan2(rot[1][0], rot[0][0])
  (angle_x, angle_y, angle_z)
}

// ============================================================================
// Public API
// ============================================================================

///|
/// Center coordinates at the origin (subtract centroid).
pub fn ma_center_coordinates(coords : Array[(Double, Double, Double)]) -> Array[(Double, Double, Double)] {
  let centroid = ma_compute_centroid(coords)
  let result : Array[(Double, Double, Double)] = Array::new()
  let mut i = 0
  while i < coords.length() {
    let (x, y, z) = coords[i]
    result.push((x - centroid.0, y - centroid.1, z - centroid.2))
    i = i + 1
  }
  result
}

///|
/// Compute internal distance matrix between two coordinate sets.
/// Returns an n×m matrix where element [i][j] is the distance
/// between coords1[i] and coords2[j].
pub fn ma_compute_distance_matrix(
  coords1 : Array[(Double, Double, Double)],
  coords2 : Array[(Double, Double, Double)],
) -> Array[Array[Double]] {
  let n = coords1.length()
  let m = coords2.length()
  let matrix : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < m {
      row.push(ma_point_distance(coords1[i], coords2[j]))
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

///|
/// Find corresponding atom pairs by distance comparison.
/// For each atom in dist_mat1, find the best matching atom in dist_mat2
/// based on distance matrix comparison.
pub fn find_corresponding_atoms(
  dist_mat1 : Array[Array[Double]],
  dist_mat2 : Array[Array[Double]],
  threshold : Double,
) -> Array[(Int, Int)] {
  let n = dist_mat1.length()
  let m = dist_mat2.length()
  let correspondences : Array[(Int, Int)] = Array::new()

  let mut i = 0
  while i < n {
    let mut best_j = -1
    let mut best_dist = threshold
    let mut j = 0
    while j < m {
      let mut dist_sum = 0.0
      let mut k = 0
      while k < n && k != i {
        let d1 = dist_mat1[i][k]
        let d2 = dist_mat2[j][k]
        let diff = d1 - d2
        dist_sum = dist_sum + diff * diff
        k = k + 1
      }
      let dist_score = dist_sum.sqrt()
      if dist_score < best_dist {
        best_dist = dist_score
        best_j = j
      }
      j = j + 1
    }
    if best_j >= 0 {
      correspondences.push((i, best_j))
    }
    i = i + 1
  }
  correspondences
}

///|
/// Compute the RMSD between two coordinate sets.
/// Both arrays should have the same length.
pub fn compute_rmsd(
  coords1 : Array[(Double, Double, Double)],
  coords2 : Array[(Double, Double, Double)],
) -> Double {
  let n = if coords1.length() < coords2.length() { coords1.length() } else { coords2.length() }
  if n == 0 {
    return 0.0
  }
  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    let dx = coords1[i].0 - coords2[i].0
    let dy = coords1[i].1 - coords2[i].1
    let dz = coords1[i].2 - coords2[i].2
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
    i = i + 1
  }
  (sum_sq / n.to_double()).sqrt()
}

///|
/// Compute optimal rotation using the Kabsch algorithm (simplified).
/// Returns rotation matrix as 3 angles (rotation around x, y, z axes).
///
/// The Kabsch algorithm steps:
/// 1. Center both coordinate sets at origin
/// 2. Compute cross-covariance matrix H = P * Q^T
/// 3. Compute optimal rotation R = H * (H^T H)^(-1/2)
/// 4. Handle reflection case (det(R) < 0)
pub fn compute_svd_rotation(
  coords1 : Array[(Double, Double, Double)],
  coords2 : Array[(Double, Double, Double)],
) -> (Double, Double, Double) {
  let n = if coords1.length() < coords2.length() { coords1.length() } else { coords2.length() }
  if n < 2 {
    return (0.0, 0.0, 0.0)
  }

  let centered1 = ma_center_coordinates(coords1)
  let centered2 = ma_center_coordinates(coords2)

  let h = ma_compute_covariance(centered1, centered2)
  let h_t = ma_transpose_3x3(h)
  let ht_h = ma_matmul_3x3(h_t, h)

  let det_ht_h = ma_det_3x3(ht_h)
  let final_rot = if det_ht_h.abs() < 1.0e-15 {
    [[1.0, 0.0, 0.0],
     [0.0, 1.0, 0.0],
     [0.0, 0.0, 1.0]]
  } else {
    let sqrt_ht_h = ma_matrix_sqrt_3x3(ht_h)
    let sqrt_ht_h_inv = ma_invert_3x3(sqrt_ht_h)
    let rot = ma_matmul_3x3(h, sqrt_ht_h_inv)

    let det_rot = ma_det_3x3(rot)
    if det_rot < 0.0 {
      [[rot[0][0], rot[0][1], -rot[0][2]],
       [rot[1][0], rot[1][1], -rot[1][2]],
       [rot[2][0], rot[2][1], -rot[2][2]]]
    } else {
      rot
    }
  }

  ma_rotation_to_angles(final_rot)
}

///|
/// Apply rotation to a set of coordinates.
/// The rotation is specified as three angles (around x, y, z axes).
pub fn apply_rotation(
  coords : Array[(Double, Double, Double)],
  rotation : (Double, Double, Double),
) -> Array[(Double, Double, Double)] {
  let rot = ma_build_rotation_matrix(rotation)
  let result : Array[(Double, Double, Double)] = Array::new()
  let mut i = 0
  while i < coords.length() {
    let (x, y, z) = coords[i]
    let nx = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z
    let ny = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z
    let nz = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z
    result.push((nx, ny, nz))
    i = i + 1
  }
  result
}

///|
/// SVD-based pairwise alignment.
/// Returns the rotation matrix (as 3 angles) and RMSD.
///
/// Aligns coords2 to coords1 using the Kabsch algorithm.
pub fn ma_align_pairwise(
  ca1 : Array[(Double, Double, Double)],
  ca2 : Array[(Double, Double, Double)],
) -> ((Double, Double, Double), Double) {
  let n = if ca1.length() < ca2.length() { ca1.length() } else { ca2.length() }
  if n < 2 {
    return ((0.0, 0.0, 0.0), 0.0)
  }

  let centered1 = ma_center_coordinates(ca1)
  let centered2 = ma_center_coordinates(ca2)

  let h = ma_compute_covariance(centered1, centered2)
  let h_t = ma_transpose_3x3(h)
  let ht_h = ma_matmul_3x3(h_t, h)

  let det_ht_h = ma_det_3x3(ht_h)
  let final_rot = if det_ht_h.abs() < 1.0e-15 {
    [[1.0, 0.0, 0.0],
     [0.0, 1.0, 0.0],
     [0.0, 0.0, 1.0]]
  } else {
    let sqrt_ht_h = ma_matrix_sqrt_3x3(ht_h)
    let sqrt_ht_h_inv = ma_invert_3x3(sqrt_ht_h)
    let rot = ma_matmul_3x3(h, sqrt_ht_h_inv)

    let det_rot = ma_det_3x3(rot)
    if det_rot < 0.0 {
      [[rot[0][0], rot[0][1], -rot[0][2]],
       [rot[1][0], rot[1][1], -rot[1][2]],
       [rot[2][0], rot[2][1], -rot[2][2]]]
    } else {
      rot
    }
  }

  let rotation = ma_rotation_to_angles(final_rot)
  let rotated2 = apply_rotation(centered2, rotation)
  let rmsd = compute_rmsd(centered1, rotated2)

  (rotation, rmsd)
}

///|
/// Compute the average structure (centroid) from multiple aligned coordinate sets.
fn ma_compute_average_structure(
  coords_list : Array[Array[(Double, Double, Double)]],
) -> Array[(Double, Double, Double)] {
  if coords_list.length() == 0 {
    return []
  }
  let n_structures = coords_list.length()
  let n_atoms = coords_list[0].length()
  let avg : Array[(Double, Double, Double)] = Array::new()
  let mut i = 0
  while i < n_atoms {
    let mut sx = 0.0
    let mut sy = 0.0
    let mut sz = 0.0
    let mut j = 0
    while j < n_structures {
      let (x, y, z) = coords_list[j][i]
      sx = sx + x
      sy = sy + y
      sz = sz + z
      j = j + 1
    }
    avg.push((sx / n_structures.to_double(), sy / n_structures.to_double(), sz / n_structures.to_double()))
    i = i + 1
  }
  avg
}

///|
/// Build per-column conservation profile from aligned sequences.
/// Each column's conservation score is the fraction of non-gap characters.
pub fn build_conservation_profile(
  aligned_sequences : Array[String],
) -> Array[Double] {
  if aligned_sequences.length() == 0 {
    return []
  }
  let n_seqs = aligned_sequences.length()
  let aln_len = aligned_sequences[0].length()
  let conservation : Array[Double] = Array::new()

  let mut col = 0
  while col < aln_len {
    let counts : Map[Char, Int] = Map([], capacity=20)
    let mut row = 0
    while row < n_seqs {
      let ch = aligned_sequences[row].unsafe_get(col).unsafe_to_char()
      if ch != '-' {
        let count = match counts.get(ch) {
          Some(c) => c
          None => 0
        }
        counts[ch] = count + 1
      }
      row = row + 1
    }
    let mut max_count = 0
    for key in counts.keys() {
      let count = counts[key]
      if count > max_count {
        max_count = count
      }
    }
    let score = if n_seqs > 0 {
      max_count.to_double() / n_seqs.to_double()
    } else {
      0.0
    }
    conservation.push(score)
    col = col + 1
  }
  conservation
}

///|
/// Build aligned sequences from aligned coordinate sets.
fn ma_build_aligned_sequences(
  structures : Array[StructureData],
  aligned_coords_list : Array[Array[(Double, Double, Double)]],
) -> Array[String] {
  let result : Array[String] = Array::new()
  let n_structures = structures.length()

  let mut i = 0
  while i < n_structures {
    let seq = structures[i].sequence
    let coords_len = aligned_coords_list[i].length()
    let orig_len = structures[i].length()

    let sb = StringBuilder::new()
    let mut j = 0
    while j < coords_len {
      if j < orig_len && j < seq.length() {
        sb.write_char(seq.unsafe_get(j).unsafe_to_char())
      } else {
        sb.write_char('-')
      }
      j = j + 1
    }
    result.push(sb.to_string())
    i = i + 1
  }
  result
}

///|
/// Compute pairwise RMSD matrix from coordinate sets.
fn ma_compute_pairwise_rmsd(
  coords_list : Array[Array[(Double, Double, Double)]],
) -> Array[Array[Double]] {
  let n = coords_list.length()
  let matrix : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      if i == j {
        row.push(0.0)
      } else {
        let (_, rmsd) = ma_align_pairwise(coords_list[i], coords_list[j])
        row.push(rmsd)
      }
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

///|
/// Get pairwise RMSD matrix from an alignment result.
pub fn get_rmsd_matrix(result : AlignmentResult) -> Array[Array[Double]] {
  result.rmsd_matrix
}

///|
/// Get per-column conservation scores from an alignment result.
pub fn get_conservation(result : AlignmentResult) -> Array[Double] {
  result.conservation
}

///|
/// Single iteration of iterative alignment.
/// Aligns all structures to the current reference and updates the
/// reference as the average of aligned structures.
///
/// Returns the new reference coordinates and the maximum RMSD change.
pub fn iterate_alignment(
  _aligner : MAaligner,
  structures : Array[StructureData],
  _iteration : Int,
) -> (Array[(Double, Double, Double)], Double) {
  let n = structures.length()
  if n == 0 {
    return ([], 0.0)
  }

  let ref_coords = structures[0].ca_coordinates.copy()
  let aligned_all : Array[Array[(Double, Double, Double)]] = Array::new()

  aligned_all.push(ref_coords.copy())

  let mut max_rmsd = 0.0
  let mut i = 1
  while i < n {
    let (rotation, rmsd) = ma_align_pairwise(ref_coords, structures[i].ca_coordinates)
    let aligned = apply_rotation(structures[i].ca_coordinates, rotation)
    aligned_all.push(aligned)
    if rmsd > max_rmsd {
      max_rmsd = rmsd
    }
    i = i + 1
  }

  let new_ref = ma_compute_average_structure(aligned_all)
  let new_centered = ma_center_coordinates(new_ref)

  (new_centered, max_rmsd)
}

///|
/// Iterative multiple structure alignment.
///
/// Algorithm:
/// 1. Start with first structure as reference
/// 2. Align each subsequent structure to the reference
/// 3. Update reference as the average of all aligned structures
/// 4. Repeat until convergence or max iterations
pub fn align_multiple(
  aligner : MAaligner,
  structures : Array[StructureData],
) -> AlignmentResult {
  let n = structures.length()
  if n == 0 {
    return AlignmentResult::new([], [], [], [], [], 0.0)
  }

  if n == 1 {
    let ids : Array[String] = [structures[0].structure_id]
    let seqs : Array[String] = [structures[0].sequence]
    let coords : Array[Array[(Double, Double, Double)]] = [structures[0].ca_coordinates.copy()]
    let rmsd_mat : Array[Array[Double]] = [[0.0]]
    let cons = build_conservation_profile(seqs)
    return AlignmentResult::new(ids, seqs, rmsd_mat, coords, cons, 0.0)
  }

  let ref_coords = structures[0].ca_coordinates.copy()
  let mut aligned_all : Array[Array[(Double, Double, Double)]] = Array::new()
  aligned_all.push(ref_coords.copy())

  let mut max_rmsd = 0.0
  let mut i = 1
  while i < n {
    let (rotation, rmsd) = ma_align_pairwise(ref_coords, structures[i].ca_coordinates)
    let aligned = apply_rotation(structures[i].ca_coordinates, rotation)
    aligned_all.push(aligned)
    if rmsd > max_rmsd {
      max_rmsd = rmsd
    }
    i = i + 1
  }

  let mut iter = 0
  while iter < aligner.max_iterations {
    let new_ref = ma_compute_average_structure(aligned_all)
    let new_aligned : Array[Array[(Double, Double, Double)]] = Array::new()
    new_aligned.push(new_ref.copy())

    let mut iter_max_rmsd = 0.0
    let mut j = 1
    while j < n {
      let (rotation, rmsd) = ma_align_pairwise(new_ref, structures[j].ca_coordinates)
      let aligned = apply_rotation(structures[j].ca_coordinates, rotation)
      new_aligned.push(aligned)
      if rmsd > iter_max_rmsd {
        iter_max_rmsd = rmsd
      }
      j = j + 1
    }

    aligned_all = new_aligned

    if iter_max_rmsd < aligner.rmsd_threshold {
      max_rmsd = iter_max_rmsd
      break
    }
    max_rmsd = iter_max_rmsd
    iter = iter + 1
  }

  let ids : Array[String] = Array::new()
  let mut k = 0
  while k < n {
    ids.push(structures[k].structure_id)
    k = k + 1
  }

  let seqs = ma_build_aligned_sequences(structures, aligned_all)
  let rmsd_mat = ma_compute_pairwise_rmsd(aligned_all)
  let cons = build_conservation_profile(seqs)

  let total_rmsd = if n > 1 {
    let mut sum = 0.0
    let mut count = 0
    let mut a = 0
    while a < n {
      let mut b = a + 1
      while b < n {
        sum = sum + rmsd_mat[a][b]
        count = count + 1
        b = b + 1
      }
      a = a + 1
    }
    if count > 0 { sum / count.to_double() } else { 0.0 }
  } else {
    0.0
  }

  AlignmentResult::new(ids, seqs, rmsd_mat, aligned_all, cons, total_rmsd)
}

///|
/// Add a structure to the aligner.
/// Creates a StructureData from the provided components.
pub fn add_structure(
  _aligner : MAaligner,
  structure_id : String,
  ca_coords : Array[(Double, Double, Double)],
  sequence : String,
  residue_names : Array[String],
) -> StructureData {
  StructureData::new(structure_id, ca_coords, sequence, residue_names)
}

///|
/// Main entry point for multiple structure alignment.
///
/// `aligner` — MAaligner with alignment parameters.
/// `structures` — array of StructureData to align.
///
/// Returns an AlignmentResult containing aligned sequences,
/// RMSD matrix, aligned coordinates, and conservation scores.
pub fn align_structures(
  aligner : MAaligner,
  structures : Array[StructureData],
) -> AlignmentResult {
  align_multiple(aligner, structures)
}