///|
/// MechStiff: mechanical stiffness matrix from ANM normal modes.
///
/// The stiffness matrix K is reconstructed from the ANM eigendecomposition:
///   K = V Λ V^T
/// where V are eigenvectors and Λ are eigenvalues.
/// The effective spring constant between residues i and j is the
/// Frobenius norm of the 3×3 block K[3i:3i+3, 3j:3j+3].
///
/// Ported from ProDy's MechStiff suite.

///|
/// Mechanical stiffness analysis result.
pub struct MechStiff {
  n : Int
  stiffness_matrix : Array[Array[Double]]
  residue_stiffness : Array[Double]
  effective_springs : Array[Array[Double]]
}

///|
/// Build MechStiff from a protein chain via ANM.
///
/// # Parameters
/// - chain: protein chain (uses CA atoms)
/// - cutoff: ANM cutoff distance (default 15.0 Å)
/// - gamma: spring constant (default 1.0)
pub fn MechStiff::new(
  chain : Chain,
  cutoff? : Double = 15.0,
  gamma? : Double = 1.0,
) -> MechStiff? {
  let anm = ANM::new(chain, cutoff~, gamma~)
  match anm {
    Some(a) => {
      let n = a.num_nodes()
      let dim = 3 * n
      let eigenvalues = a.get_eigenvalues()
      let eigenvectors = a.get_eigenvectors()
      // Reconstruct stiffness matrix: K = V Λ V^T
      // K[i][j] = Σ_k λ_k * V[i][k] * V[j][k]
      // Skip zero eigenvalues (rigid body modes)
      let n_modes = eigenvalues.length()
      let k_matrix : Array[Array[Double]] = Array::makei(dim, { _ =>
        Array::make(dim, 0.0)
      })
      let mut mode = 0
      while mode < n_modes {
        let lambda = eigenvalues[mode]
        if lambda < 0.0000000001 { mode = mode + 1; continue }
        let mut i = 0
        while i < dim {
          let vi = eigenvectors[i][mode]
          let mut j = 0
          while j < dim {
            k_matrix[i][j] = k_matrix[i][j] + lambda * vi * eigenvectors[j][mode]
            j = j + 1
          }
          i = i + 1
        }
        mode = mode + 1
      }
      // Effective spring constants between residues
      // S[i][j] = ||K[3i:3i+3, 3j:3j+3]||_F (Frobenius norm)
      let springs : Array[Array[Double]] = Array::makei(n, { _ =>
        Array::make(n, 0.0)
      })
      let mut ri = 0
      while ri < n {
        let mut rj = 0
        while rj < n {
          let mut norm_sq = 0.0
          let mut bi = 0
          while bi < 3 {
            let mut bj = 0
            while bj < 3 {
              let val = k_matrix[3 * ri + bi][3 * rj + bj]
              norm_sq = norm_sq + val * val
              bj = bj + 1
            }
            bi = bi + 1
          }
          springs[ri][rj] = cav_sqrt(norm_sq)
          rj = rj + 1
        }
        ri = ri + 1
      }
      // Per-residue stiffness: average spring to all other residues
      let residue_stiff = Array::make(n, 0.0)
      let mut r = 0
      while r < n {
        let mut sum = 0.0
        let mut count = 0
        let mut s = 0
        while s < n {
          if s != r {
            sum = sum + springs[r][s]
            count = count + 1
          }
          s = s + 1
        }
        residue_stiff[r] = if count > 0 { sum / count.to_double() } else { 0.0 }
        r = r + 1
      }
      Some(MechStiff::{
        n, stiffness_matrix: k_matrix, residue_stiffness: residue_stiff,
        effective_springs: springs,
      })
    }
    None => None
  }
}

///|
/// Get the 3×3 stiffness block between residues i and j.
pub fn MechStiff::stiffness_between(
  self : MechStiff,
  i : Int,
  j : Int,
) -> Array[Array[Double]] {
  let block = Array::makei(3, { _ => Array::make(3, 0.0) })
  let mut bi = 0
  while bi < 3 {
    let mut bj = 0
    while bj < 3 {
      block[bi][bj] = self.stiffness_matrix[3 * i + bi][3 * j + bj]
      bj = bj + 1
    }
    bi = bi + 1
  }
  block
}

///|
/// Get the effective spring constant between residues i and j.
pub fn MechStiff::effective_spring(
  self : MechStiff,
  i : Int,
  j : Int,
) -> Double {
  self.effective_springs[i][j]
}

///|
/// Get the per-residue stiffness profile.
pub fn MechStiff::stiffness_profile(self : MechStiff) -> Array[Double] {
  self.residue_stiffness
}

///|
/// Get the number of residues.
pub fn MechStiff::num_residues(self : MechStiff) -> Int {
  self.n
}

///|
/// Get the full stiffness matrix (3N × 3N).
pub fn MechStiff::get_matrix(self : MechStiff) -> Array[Array[Double]] {
  self.stiffness_matrix
}

///|
/// Get the effective springs matrix (N × N).
pub fn MechStiff::get_springs(self : MechStiff) -> Array[Array[Double]] {
  self.effective_springs
}

///|
/// Compute the stiffness-based distance profile.
/// For each residue pair, the mechanical distance is 1/S[i][j].
/// Returns the mean mechanical distance for each residue.
pub fn MechStiff::mechanical_distance_profile(self : MechStiff) -> Array[Double] {
  let n = self.n
  let result = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let mut sum = 0.0
    let mut count = 0
    let mut j = 0
    while j < n {
      if i != j {
        let s = self.effective_springs[i][j]
        if s > 0.0000000001 {
          sum = sum + 1.0 / s
          count = count + 1
        }
      }
      j = j + 1
    }
    result[i] = if count > 0 { sum / count.to_double() } else { 0.0 }
    i = i + 1
  }
  result
}

///|
/// Identify stiff residues (above average stiffness).
/// Returns indices of residues with above-average stiffness.
pub fn MechStiff::stiff_residues(self : MechStiff) -> Array[Int] {
  let n = self.n
  let mut mean = 0.0
  for s in self.residue_stiffness { mean = mean + s }
  if n > 0 { mean = mean / n.to_double() }
  let result = Array(capacity=0)
  let mut i = 0
  while i < n {
    if self.residue_stiffness[i] > mean {
      result.push(i)
    }
    i = i + 1
  }
  result
}
