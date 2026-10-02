///|
/// ProDy-style Anisotropic Network Model (ANM).
///
/// ANM extends GNM by using a 3N×3N Hessian matrix that captures
/// directional information about protein motions.  Each node interacts
/// with neighbors via a spring of constant γ/d².
///
/// Reference: Atilgan, Durell, Atilgan, Bahar (2001) Biophys J 80:505-515.

///|
/// Result of an ANM calculation.
pub struct ANM {
  n : Int
  cutoff : Double
  gamma : Double
  hessian : Array[Array[Double]]
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]]
}

///|
/// Build an ANM from a protein chain using CA atoms.
pub fn ANM::new(
  chain : Chain,
  cutoff? : Double = 15.0,
  gamma? : Double = 1.0,
) -> ANM? {
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
  // Build Hessian: for each pair (i,j) in contact,
  // H_ij = -gamma/d² * [dxdx dydy dzdz] (3×3 super-element)
  // H_ii = -sum_j H_ij
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
        // Off-diagonal block H_ij (negative)
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
        // Symmetric block H_ji
        h[j3][i3] = h[j3][i3] - hxx
        h[j3][i3 + 1] = h[j3][i3 + 1] - hxy
        h[j3][i3 + 2] = h[j3][i3 + 2] - hxz
        h[j3 + 1][i3] = h[j3 + 1][i3] - hxy
        h[j3 + 1][i3 + 1] = h[j3 + 1][i3 + 1] - hyy
        h[j3 + 1][i3 + 2] = h[j3 + 1][i3 + 2] - hyz
        h[j3 + 2][i3] = h[j3 + 2][i3] - hxz
        h[j3 + 2][i3 + 1] = h[j3 + 2][i3 + 1] - hyz
        h[j3 + 2][i3 + 2] = h[j3 + 2][i3 + 2] - hzz
        // Diagonal block H_ii (positive)
        h[i3][i3] = h[i3][i3] + hxx
        h[i3][i3 + 1] = h[i3][i3 + 1] + hxy
        h[i3][i3 + 2] = h[i3][i3 + 2] + hxz
        h[i3 + 1][i3] = h[i3 + 1][i3] + hxy
        h[i3 + 1][i3 + 1] = h[i3 + 1][i3 + 1] + hyy
        h[i3 + 1][i3 + 2] = h[i3 + 1][i3 + 2] + hyz
        h[i3 + 2][i3] = h[i3 + 2][i3] + hxz
        h[i3 + 2][i3 + 1] = h[i3 + 2][i3 + 1] + hyz
        h[i3 + 2][i3 + 2] = h[i3 + 2][i3 + 2] + hzz
        // Diagonal block H_jj (positive)
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
  // Eigen-decompose (6 zero modes expected: 3 translation + 3 rotation)
  let max_iter = if dim > 100 { dim * 20 } else { 200 }
  let (eigenvalues, eigenvectors) = linalg_jacobi(h, dim, max_iter)
  Some(ANM::{ n, cutoff, gamma, hessian: h, eigenvalues, eigenvectors })
}

///|
/// Get the Hessian matrix (3N × 3N).
pub fn ANM::hessian_matrix(self : ANM) -> Array[Array[Double]] {
  self.hessian
}

///|
/// Get sorted eigenvalues (ascending, first 6 ~0).
pub fn ANM::get_eigenvalues(self : ANM) -> Array[Double] {
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
pub fn ANM::get_eigenvectors(self : ANM) -> Array[Array[Double]] {
  self.eigenvectors
}

///|
/// Predict B-factors from ANM: B_i = sum_k (1/lambda_k) * |v_k[3i:3i+3]|²
pub fn ANM::predict_bfactors(self : ANM) -> Array[Double] {
  let n = self.n
  let dim = 3 * n
  let bfactors = Array::make(n, 0.0)
  // Sort eigenvalue indices ascending
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
  // Skip first 6 zero modes (3 trans + 3 rot)
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
/// Mean-square fluctuations (same as B-factors up to scaling).
pub fn ANM::fluctuations(self : ANM) -> Array[Double] {
  self.predict_bfactors()
}

///|
/// Get the k slowest non-trivial modes.
/// Returns array of (mode_index, eigenvalue, eigenvector).
pub fn ANM::slow_modes(self : ANM, k : Int) -> Array[(Int, Double, Array[Double])] {
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
/// C_ij = sum_k (1/lambda_k) * v_k[3i:3i+3] · v_k[3j:3j+3]
pub fn ANM::cross_correlation(self : ANM) -> Array[Array[Double]] {
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
  // Normalize to [-1, 1]
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
/// Number of nodes (CA atoms).
pub fn ANM::num_nodes(self : ANM) -> Int {
  self.n
}
