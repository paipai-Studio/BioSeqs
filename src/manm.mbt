///|
/// mANM — Membrane ANM.
///
/// Membrane proteins experience an asymmetric environment: the lipid
/// bilayer restricts motion in the membrane-normal direction (z) while
/// allowing lateral diffusion (x-y).  mANM models this by adding a
/// position restraint to membrane-embedded residues in the z-direction.
///
/// For each membrane residue i:
///   H[3i+2][3i+2] += k_membrane        (normal restraint)
///   H[3i][3i]   += k_lateral           (optional lateral restraint)
///   H[3i+1][3i+1] += k_lateral
///
/// The normal restraint represents the membrane's confining potential;
/// the lateral restraint (if enabled) represents lipid-protein friction.
///
/// Reference: Lezon, Bahar (2010) "Relaxation of functional residues
///            by network rigidity."  Also: Zheng, Tekpinar, Bahar
///            (2010) Biophys J 99:1778-1786.

///|
pub struct MANM {
  n : Int
  cutoff : Double
  gamma : Double
  k_membrane : Double
  k_lateral : Double
  n_membrane : Int
  hessian : Array[Array[Double]]
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]]
}

///|
/// Build an mANM from a protein chain and a list of membrane residue indices.
///
/// * `membrane_residues` — 0-based indices of residues in the membrane
/// * `cutoff`     — interaction range in Å (default 15.0)
/// * `gamma`      — intra-protein spring constant (default 1.0)
/// * `k_membrane` — normal-direction restraint strength (default 1.0)
/// * `k_lateral`  — lateral lipid contact strength (default 0.0 = disabled)
pub fn MANM::new(
  chain : Chain,
  membrane_residues : Array[Int],
  cutoff? : Double = 15.0,
  gamma? : Double = 1.0,
  k_membrane? : Double = 1.0,
  k_lateral? : Double = 0.0,
) -> MANM? {
  let residues = chain.get_residues()
  let ca_coords : Array[Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => ca_coords.push(ca.coord)
      None => ()
    }
  }
  let n = ca_coords.length()
  if n < 3 {
    return None
  }
  let dim = 3 * n
  let h : Array[Array[Double]] = Array::makei(dim, { _ => Array::make(dim, 0.0) })

  // Step 1: Standard ANM Hessian (intra-protein contacts)
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let dx = ca_coords[j].x - ca_coords[i].x
      let dy = ca_coords[j].y - ca_coords[i].y
      let dz = ca_coords[j].z - ca_coords[i].z
      let d2 = dx * dx + dy * dy + dz * dz
      let d = d2.sqrt()
      if d <= cutoff && d > 0.001 {
        let g = gamma / d2
        let hxx = g * dx * dx
        let hxy = g * dx * dy
        let hxz = g * dx * dz
        let hyy = g * dy * dy
        let hyz = g * dy * dz
        let hzz = g * dz * dz
        let i3 = 3 * i
        let j3 = 3 * j
        h[i3][j3] = h[i3][j3] - hxx
        h[i3][j3 + 1] = h[i3][j3 + 1] - hxy
        h[i3][j3 + 2] = h[i3][j3 + 2] - hxz
        h[i3 + 1][j3] = h[i3 + 1][j3] - hxy
        h[i3 + 1][j3 + 1] = h[i3 + 1][j3 + 1] - hyy
        h[i3 + 1][j3 + 2] = h[i3 + 1][j3 + 2] - hyz
        h[i3 + 2][j3] = h[i3 + 2][j3] - hxz
        h[i3 + 2][j3 + 1] = h[i3 + 2][j3 + 1] - hyz
        h[i3 + 2][j3 + 2] = h[i3 + 2][j3 + 2] - hzz
        h[j3][i3] = h[j3][i3] - hxx
        h[j3][i3 + 1] = h[j3][i3 + 1] - hxy
        h[j3][i3 + 2] = h[j3][i3 + 2] - hxz
        h[j3 + 1][i3] = h[j3 + 1][i3] - hxy
        h[j3 + 1][i3 + 1] = h[j3 + 1][i3 + 1] - hyy
        h[j3 + 1][i3 + 2] = h[j3 + 1][i3 + 2] - hyz
        h[j3 + 2][i3] = h[j3 + 2][i3] - hxz
        h[j3 + 2][i3 + 1] = h[j3 + 2][i3 + 1] - hyz
        h[j3 + 2][i3 + 2] = h[j3 + 2][i3 + 2] - hzz
        h[i3][i3] = h[i3][i3] + hxx
        h[i3][i3 + 1] = h[i3][i3 + 1] + hxy
        h[i3][i3 + 2] = h[i3][i3 + 2] + hxz
        h[i3 + 1][i3] = h[i3 + 1][i3] + hxy
        h[i3 + 1][i3 + 1] = h[i3 + 1][i3 + 1] + hyy
        h[i3 + 1][i3 + 2] = h[i3 + 1][i3 + 2] + hyz
        h[i3 + 2][i3] = h[i3 + 2][i3] + hxz
        h[i3 + 2][i3 + 1] = h[i3 + 2][i3 + 1] + hyz
        h[i3 + 2][i3 + 2] = h[i3 + 2][i3 + 2] + hzz
        h[j3][j3] = h[j3][j3] + hxx
        h[j3][j3 + 1] = h[j3][j3 + 1] + hxy
        h[j3][j3 + 2] = h[j3][j3 + 2] + hxz
        h[j3 + 1][j3] = h[j3 + 1][j3] + hxy
        h[j3 + 1][j3 + 1] = h[j3 + 1][j3 + 1] + hyy
        h[j3 + 1][j3 + 2] = h[j3 + 1][j3 + 2] + hyz
        h[j3 + 2][j3] = h[j3 + 2][j3] + hxz
        h[j3 + 2][j3 + 1] = h[j3 + 2][j3 + 1] + hyz
        h[j3 + 2][j3 + 2] = h[j3 + 2][j3 + 2] + hzz
      }
      j = j + 1
    }
    i = i + 1
  }

  // Step 2: Membrane restraints on membrane-embedded residues
  for idx in membrane_residues {
    if idx >= 0 && idx < n {
      let i3 = 3 * idx
      // Normal restraint (z-direction)
      h[i3 + 2][i3 + 2] = h[i3 + 2][i3 + 2] + k_membrane
      // Lateral restraint (x-y, optional)
      if k_lateral > 0.0 {
        h[i3][i3] = h[i3][i3] + k_lateral
        h[i3 + 1][i3 + 1] = h[i3 + 1][i3 + 1] + k_lateral
      }
    }
  }

  let max_iter = if dim > 100 { dim * 20 } else { 200 }
  let (eigenvalues, eigenvectors) = linalg_jacobi(h, dim, max_iter)
  Some(MANM::{
    n,
    cutoff,
    gamma,
    k_membrane,
    k_lateral,
    n_membrane: membrane_residues.length(),
    hessian: h,
    eigenvalues,
    eigenvectors,
  })
}

///|
/// Identify membrane residues by z-coordinate range.
/// Returns 0-based indices of residues whose CA z ∈ [z_min, z_max].
pub fn membrane_residues_by_z(
  chain : Chain,
  z_min : Double,
  z_max : Double,
) -> Array[Int] {
  let residues = chain.get_residues()
  let result : Array[Int] = Array(capacity=0)
  for i = 0; i < residues.length(); i = i + 1 {
    match residues[i].get_ca() {
      Some(ca) => {
        if ca.coord.z >= z_min && ca.coord.z <= z_max {
          result.push(i)
        }
      }
      None => ()
    }
  }
  result
}

///|
/// Get the Hessian matrix (3N × 3N).
pub fn MANM::hessian_matrix(self : MANM) -> Array[Array[Double]] {
  self.hessian
}

///|
/// Get sorted eigenvalues (ascending, first 6 ~0).
pub fn MANM::get_eigenvalues(self : MANM) -> Array[Double] {
  let dim = 3 * self.n
  let vals = Array::make(dim, 0.0)
  for i = 0; i < dim; i = i + 1 {
    vals[i] = self.eigenvalues[i]
  }
  for i = 1; i < dim; i = i + 1 {
    let key = vals[i]
    let mut j = i - 1
    while j >= 0 && vals[j] > key {
      vals[j + 1] = vals[j]
      j = j - 1
    }
    vals[j + 1] = key
  }
  vals
}

///|
/// Get eigenvector matrix.
pub fn MANM::get_eigenvectors(self : MANM) -> Array[Array[Double]] {
  self.eigenvectors
}

///|
/// Predict B-factors: B_i = sum_k (1/λ_k) * |v_k[3i:3i+3]|²
pub fn MANM::predict_bfactors(self : MANM) -> Array[Double] {
  let n = self.n
  let dim = 3 * n
  let bfactors = Array::make(n, 0.0)
  let indices = Array::makei(dim, { i => i })
  for i = 1; i < dim; i = i + 1 {
    let key = indices[i]
    let key_val = self.eigenvalues[key]
    let mut j = i - 1
    while j >= 0 && self.eigenvalues[indices[j]] > key_val {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
  }
  for k = 6; k < dim; k = k + 1 {
    let lambda = self.eigenvalues[indices[k]]
    if lambda.abs() < 1.0e-10 {
      continue
    }
    for i = 0; i < n; i = i + 1 {
      let vx = self.eigenvectors[3 * i][indices[k]]
      let vy = self.eigenvectors[3 * i + 1][indices[k]]
      let vz = self.eigenvectors[3 * i + 2][indices[k]]
      bfactors[i] = bfactors[i] + (vx * vx + vy * vy + vz * vz) / lambda
    }
  }
  for i = 0; i < n; i = i + 1 {
    if bfactors[i] < 0.0 {
      bfactors[i] = -bfactors[i]
    }
  }
  bfactors
}

///|
/// Mean-square fluctuations.
pub fn MANM::fluctuations(self : MANM) -> Array[Double] {
  self.predict_bfactors()
}

///|
/// Get the k slowest non-trivial modes (skip 6 zero modes).
pub fn MANM::slow_modes(
  self : MANM,
  k : Int,
) -> Array[(Int, Double, Array[Double])] {
  let n = self.n
  let dim = 3 * n
  let indices = Array::makei(dim, { i => i })
  for i = 1; i < dim; i = i + 1 {
    let key = indices[i]
    let key_val = self.eigenvalues[key]
    let mut j = i - 1
    while j >= 0 && self.eigenvalues[indices[j]] > key_val {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
  }
  let result : Array[(Int, Double, Array[Double])] = Array(capacity=0)
  let count = if k < dim - 6 { k } else { dim - 6 }
  for m = 6; m < 6 + count; m = m + 1 {
    let idx = indices[m]
    let eigval = self.eigenvalues[idx]
    let eigvec = Array::make(dim, 0.0)
    for i = 0; i < dim; i = i + 1 {
      eigvec[i] = self.eigenvectors[i][idx]
    }
    result.push((idx, eigval, eigvec))
  }
  result
}

///|
/// Cross-correlation matrix between residues.
pub fn MANM::cross_correlation(self : MANM) -> Array[Array[Double]] {
  let n = self.n
  let dim = 3 * n
  let corr : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  let indices = Array::makei(dim, { i => i })
  for i = 1; i < dim; i = i + 1 {
    let key = indices[i]
    let key_val = self.eigenvalues[key]
    let mut j = i - 1
    while j >= 0 && self.eigenvalues[indices[j]] > key_val {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
  }
  for k = 6; k < dim; k = k + 1 {
    let lambda = self.eigenvalues[indices[k]]
    if lambda.abs() < 1.0e-10 {
      continue
    }
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        let vix = self.eigenvectors[3 * i][indices[k]]
        let viy = self.eigenvectors[3 * i + 1][indices[k]]
        let viz = self.eigenvectors[3 * i + 2][indices[k]]
        let vjx = self.eigenvectors[3 * j][indices[k]]
        let vjy = self.eigenvectors[3 * j + 1][indices[k]]
        let vjz = self.eigenvectors[3 * j + 2][indices[k]]
        corr[i][j] = corr[i][j] + (vix * vjx + viy * vjy + viz * vjz) / lambda
      }
    }
  }
  let mut max_val = 1.0e-10
  for i = 0; i < n; i = i + 1 {
    if corr[i][i].abs() > max_val {
      max_val = corr[i][i].abs()
    }
  }
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      corr[i][j] = corr[i][j] / max_val
    }
  }
  corr
}

///|
/// Number of nodes.
pub fn MANM::num_nodes(self : MANM) -> Int {
  self.n
}

///|
/// Number of membrane-embedded residues.
pub fn MANM::num_membrane(self : MANM) -> Int {
  self.n_membrane
}

///|
/// Membrane restraint strength (normal direction).
pub fn MANM::membrane_restraint(self : MANM) -> Double {
  self.k_membrane
}

///|
/// Lateral restraint strength.
pub fn MANM::lateral_restraint(self : MANM) -> Double {
  self.k_lateral
}
