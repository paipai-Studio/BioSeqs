///|
/// betaANM — distance-dependent force constant ANM.
///
/// In standard ANM, all pairs within the cutoff share the same force
/// constant γ.  In betaANM, the force constant decays with distance:
///   γ_ij = γ / d_ij^β
/// so closer residues interact more strongly.  Typical β = 2.
///
/// The Hessian super-element becomes:
///   H_ij = -(γ / d^β / d²) * [dx*dx dx*dy dx*dz; ...]
///
/// Reference: Kondrashov, Van Wynsberghe, Bhatt, Phillips, Phillips (2006)
///            Proteins 64:495-505.

///|
pub struct BetaANM {
  n : Int
  cutoff : Double
  gamma : Double
  beta : Double
  hessian : Array[Array[Double]]
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]]
}

///|
/// Build a betaANM from a protein chain using CA atoms.
///
/// * `cutoff` — interaction range in Å (default 15.0)
/// * `gamma`  — overall spring constant (default 1.0)
/// * `beta`   — distance exponent, higher = more local (default 2.0)
pub fn BetaANM::new(
  chain : Chain,
  cutoff? : Double = 15.0,
  gamma? : Double = 1.0,
  beta? : Double = 2.0,
) -> BetaANM? {
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
        let d_beta = @math.exp(beta * @math.ln(d))
        let g = gamma / (d_beta * d2)
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
  let max_iter = if dim > 100 { dim * 20 } else { 200 }
  let (eigenvalues, eigenvectors) = linalg_jacobi(h, dim, max_iter)
  Some(BetaANM::{ n, cutoff, gamma, beta, hessian: h, eigenvalues, eigenvectors })
}

///|
/// Get the Hessian matrix (3N × 3N).
pub fn BetaANM::hessian_matrix(self : BetaANM) -> Array[Array[Double]] {
  self.hessian
}

///|
/// Get sorted eigenvalues (ascending, first 6 ~0).
pub fn BetaANM::get_eigenvalues(self : BetaANM) -> Array[Double] {
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
pub fn BetaANM::get_eigenvectors(self : BetaANM) -> Array[Array[Double]] {
  self.eigenvectors
}

///|
/// Predict B-factors: B_i = sum_k (1/λ_k) * |v_k[3i:3i+3]|²
pub fn BetaANM::predict_bfactors(self : BetaANM) -> Array[Double] {
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
pub fn BetaANM::fluctuations(self : BetaANM) -> Array[Double] {
  self.predict_bfactors()
}

///|
/// Get the k slowest non-trivial modes (skip 6 zero modes).
pub fn BetaANM::slow_modes(
  self : BetaANM,
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
pub fn BetaANM::cross_correlation(self : BetaANM) -> Array[Array[Double]] {
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
pub fn BetaANM::num_nodes(self : BetaANM) -> Int {
  self.n
}

///|
/// The beta exponent used.
pub fn BetaANM::beta_exponent(self : BetaANM) -> Double {
  self.beta
}

// =========================================================================
// pfANM — Parameter-Free ANM
// =========================================================================

///|
/// pfANM — Parameter-Free Anisotropic Network Model.
///
/// Eliminates the need for a cutoff distance by using a steeply
/// decaying force constant:
///   γ_ij = γ / d_ij^6
/// The inverse-sixth-power law mimics van der Waals attraction,
/// so distant pairs contribute negligibly without an explicit cutoff.
///
/// Reference: Yang, Song, Jernigan (2009) PMC Biophys 2:11.

///|
pub struct PfANM {
  n : Int
  gamma : Double
  hessian : Array[Array[Double]]
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]]
}

///|
/// Build a pfANM from a protein chain using CA atoms.
///
/// * `gamma`  — overall scale (default 1.0)
/// * `cutoff` — practical cutoff to skip negligible pairs (default 30.0 Å)
pub fn PfANM::new(
  chain : Chain,
  gamma? : Double = 1.0,
  cutoff? : Double = 30.0,
) -> PfANM? {
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
        // Force constant: γ_ij = γ / d^6
        // Hessian element: g = γ / (d^6 * d²) = γ / d^8
        let d6 = d2 * d2 * d2
        let g = gamma / (d6 * d2)
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
  let max_iter = if dim > 100 { dim * 20 } else { 200 }
  let (eigenvalues, eigenvectors) = linalg_jacobi(h, dim, max_iter)
  Some(PfANM::{ n, gamma, hessian: h, eigenvalues, eigenvectors })
}

///|
/// Get the Hessian matrix (3N × 3N).
pub fn PfANM::hessian_matrix(self : PfANM) -> Array[Array[Double]] {
  self.hessian
}

///|
/// Get sorted eigenvalues (ascending, first 6 ~0).
pub fn PfANM::get_eigenvalues(self : PfANM) -> Array[Double] {
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
pub fn PfANM::get_eigenvectors(self : PfANM) -> Array[Array[Double]] {
  self.eigenvectors
}

///|
/// Predict B-factors: B_i = sum_k (1/λ_k) * |v_k[3i:3i+3]|²
pub fn PfANM::predict_bfactors(self : PfANM) -> Array[Double] {
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
pub fn PfANM::fluctuations(self : PfANM) -> Array[Double] {
  self.predict_bfactors()
}

///|
/// Get the k slowest non-trivial modes.
pub fn PfANM::slow_modes(
  self : PfANM,
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
pub fn PfANM::cross_correlation(self : PfANM) -> Array[Array[Double]] {
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
pub fn PfANM::num_nodes(self : PfANM) -> Int {
  self.n
}
