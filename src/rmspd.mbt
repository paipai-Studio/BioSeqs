///|
/// RMSD of pair distances (rmspd) — internal distance-based comparison.
///
/// Unlike standard RMSD which requires optimal superimposition, rmspd
/// compares the set of all pairwise distances between corresponding
/// atoms. This is rotation/translation invariant by construction.
///
/// Biotite provides this as an alternative to RMSD for structures that
/// cannot be meaningfully superimposed (e.g., intrinsically disordered).

///|
/// Calculate RMSD of all pairwise atom distances.
///
/// For two structures with n atoms each, computes:
///   rmspd = sqrt( (2 / (n*(n-1))) * sum_{i<j} (d1_ij - d2_ij)^2 )
///
/// where d1_ij = |pos1[i] - pos1[j]| and d2_ij = |pos2[i] - pos2[j]|.
///
/// This is invariant to rotation and translation, making it useful for
/// comparing conformations without alignment.
pub fn rmspd(coords1 : Array[Vector3], coords2 : Array[Vector3]) -> Double {
  let n = coords1.length()
  if n < 2 || coords2.length() != n {
    return 0.0
  }
  let mut sum_sq = 0.0
  let mut count = 0
  for i in 0..<n {
    for j in (i + 1)..<n {
      let d1 = coords1[i].distance(coords1[j])
      let d2 = coords2[i].distance(coords2[j])
      let diff = d1 - d2
      sum_sq = sum_sq + diff * diff
      count = count + 1
    }
  }
  if count == 0 {
    return 0.0
  }
  (sum_sq / count.to_double()).sqrt()
}

///|
/// Calculate weighted RMSD of pair distances.
///
/// Weights can be used to emphasize certain atom pairs (e.g., C-alpha
/// pairs within a certain sequence separation).
pub fn rmspd_weighted(
  coords1 : Array[Vector3],
  coords2 : Array[Vector3],
  weights : Array[Double],
) -> Double {
  let n = coords1.length()
  if n < 2 || coords2.length() != n {
    return 0.0
  }
  let mut sum_sq = 0.0
  let mut sum_w = 0.0
  for i in 0..<n {
    for j in (i + 1)..<n {
      let d1 = coords1[i].distance(coords1[j])
      let d2 = coords2[i].distance(coords2[j])
      let diff = d1 - d2
      let w = if i < weights.length() && j < weights.length() {
        weights[i] * weights[j]
      } else {
        1.0
      }
      sum_sq = sum_sq + w * diff * diff
      sum_w = sum_w + w
    }
  }
  if sum_w == 0.0 {
    return 0.0
  }
  (sum_sq / sum_w).sqrt()
}

///|
/// Calculate RMSD of pair distances for a subset of atoms.
///
/// Only considers pairs where both atoms are in the selected indices.
pub fn rmspd_subset(
  coords1 : Array[Vector3],
  coords2 : Array[Vector3],
  indices : Array[Int32],
) -> Double {
  let n = indices.length()
  if n < 2 {
    return 0.0
  }
  let mut sum_sq = 0.0
  let mut count = 0
  for i in 0..<n {
    for j in (i + 1)..<n {
      let idx_i = indices[i]
      let idx_j = indices[j]
      if idx_i >= coords1.length() || idx_j >= coords1.length() {
        continue
      }
      if idx_i >= coords2.length() || idx_j >= coords2.length() {
        continue
      }
      let d1 = coords1[idx_i].distance(coords1[idx_j])
      let d2 = coords2[idx_i].distance(coords2[idx_j])
      let diff = d1 - d2
      sum_sq = sum_sq + diff * diff
      count = count + 1
    }
  }
  if count == 0 {
    return 0.0
  }
  (sum_sq / count.to_double()).sqrt()
}

///|
/// Calculate the distance difference matrix.
///
/// Returns a matrix where element [i][j] = d1_ij - d2_ij.
/// Useful for identifying which atom pairs contribute most to rmspd.
pub fn distance_difference_matrix(
  coords1 : Array[Vector3],
  coords2 : Array[Vector3],
) -> Array[Array[Double]] {
  let n = coords1.length()
  let result : Array[Array[Double]] = []
  for i in 0..<n {
    let row : Array[Double] = []
    for j in 0..<n {
      if i == j {
        row.push(0.0)
      } else {
        let d1 = coords1[i].distance(coords1[j])
        let d2 = coords2[i].distance(coords2[j])
        row.push(d1 - d2)
      }
    }
    result.push(row)
  }
  result
}
