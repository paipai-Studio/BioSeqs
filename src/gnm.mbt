///|
/// ProDy-style Gaussian Network Model (GNM).
///
/// GNM models protein fluctuations using a Kirchhoff (Laplacian) matrix
/// built from CA-CA contacts.  Eigen-decomposition yields slow modes
/// (large-scale motions) and B-factor predictions.
///
/// Reference: Bahar, Atilgan, Erman (1997) Folding & Design 2:173-181.

///|
/// Result of a GNM calculation.
pub struct GNM {
  n : Int
  cutoff : Double
  gamma : Double
  kirchhoff : Array[Array[Double]]
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]]
}

///|
/// Build a GNM from a protein chain using CA atoms.
/// cutoff: contact distance threshold (default 10.0 Å)
/// gamma: spring constant (default 1.0)
pub fn GNM::new(
  chain : Chain,
  cutoff? : Double = 10.0,
  gamma? : Double = 1.0,
) -> GNM? {
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
  // Build Kirchhoff matrix: off-diagonal = -gamma if in contact, diagonal = sum of contacts
  let kirchhoff : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let d = ca_coords[i].distance(ca_coords[j])
      if d <= cutoff {
        kirchhoff[i][j] = -gamma
        kirchhoff[j][i] = -gamma
      }
      j = j + 1
    }
    i = i + 1
  }
  // Diagonal = negative sum of off-diagonal row elements
  for i = 0; i < n; i = i + 1 {
    let mut degree = 0.0
    for j = 0; j < n; j = j + 1 {
      if i != j {
        degree = degree - kirchhoff[i][j]
      }
    }
    kirchhoff[i][i] = degree
  }
  // Eigen-decompose (skip first zero eigenvalue mode)
  let (eigenvalues, eigenvectors) = linalg_jacobi(kirchhoff, n, 200)
  Some(GNM::{ n, cutoff, gamma, kirchhoff, eigenvalues, eigenvectors })
}

///|
/// Get the Kirchhoff matrix.
pub fn GNM::kirchhoff_matrix(self : GNM) -> Array[Array[Double]] {
  self.kirchhoff
}

///|
/// Get eigenvalues (sorted ascending, first is ~0).
pub fn GNM::get_eigenvalues(self : GNM) -> Array[Double] {
  let vals = Array::make(self.n, 0.0)
  for i = 0; i < self.n; i = i + 1 {
    vals[i] = self.eigenvalues[i]
  }
  // Simple insertion sort ascending
  for i = 1; i < self.n; i = i + 1 {
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
/// Get eigenvector matrix (columns are eigenvectors).
pub fn GNM::get_eigenvectors(self : GNM) -> Array[Array[Double]] {
  self.eigenvectors
}

///|
/// Predict B-factors from GNM: B_i = sum_k (1/lambda_k) * v_k[i]^2
/// Skip the first (zero) eigenvalue.
pub fn GNM::predict_bfactors(self : GNM) -> Array[Double] {
  let n = self.n
  let bfactors = Array::make(n, 0.0)
  // Sort eigenvalue indices ascending, skip first (near-zero)
  let indices = Array::makei(n, { i => i })
  // Sort indices by eigenvalue
  for i = 1; i < n; i = i + 1 {
    let key = indices[i]
    let key_val = self.eigenvalues[key]
    let mut j = i - 1
    while j >= 0 && self.eigenvalues[indices[j]] > key_val {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
  }
  // Sum over modes 1..n-1 (skip mode 0 = zero eigenvalue)
  for k = 1; k < n; k = k + 1 {
    let lambda = self.eigenvalues[indices[k]]
    if lambda.abs() < 1.0e-10 {
      continue
    }
    for i = 0; i < n; i = i + 1 {
      let v = self.eigenvectors[i][indices[k]]
      bfactors[i] = bfactors[i] + v * v / lambda
    }
  }
  // Scale to positive
  for i = 0; i < n; i = i + 1 {
    if bfactors[i] < 0.0 {
      bfactors[i] = -bfactors[i]
    }
  }
  bfactors
}

///|
/// Mean-square fluctuations: <delta r_i^2> = sum_k (1/lambda_k) * v_k[i]^2
/// Same as B-factors up to a constant.
pub fn GNM::fluctuations(self : GNM) -> Array[Double] {
  self.predict_bfactors()
}

///|
/// Get the k slowest modes (smallest non-zero eigenvalues).
/// Returns array of (mode_index, eigenvalue, eigenvector).
pub fn GNM::slow_modes(self : GNM, k : Int) -> Array[(Int, Double, Array[Double])] {
  let n = self.n
  let indices = Array::makei(n, { i => i })
  for i = 1; i < n; i = i + 1 {
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
  let count = if k < n - 1 { k } else { n - 1 }
  for m = 1; m <= count; m = m + 1 {
    let idx = indices[m]
    let eigval = self.eigenvalues[idx]
    let eigvec = Array::make(n, 0.0)
    for i = 0; i < n; i = i + 1 {
      eigvec[i] = self.eigenvectors[i][idx]
    }
    result.push((idx, eigval, eigvec))
  }
  result
}

///|
/// Number of residues (CA atoms).
pub fn GNM::num_nodes(self : GNM) -> Int {
  self.n
}

///|
/// Contact count (number of non-zero off-diagonal elements / 2).
pub fn GNM::num_contacts(self : GNM) -> Int {
  let mut count = 0
  for i = 0; i < self.n; i = i + 1 {
    for j = i + 1; j < self.n; j = j + 1 {
      if self.kirchhoff[i][j].abs() > 1.0e-10 {
        count = count + 1
      }
    }
  }
  count
}

///|
/// Build Kirchhoff matrix directly from CA coordinates.
pub fn gnm_kirchhoff(
  coords : Array[Vector3],
  cutoff? : Double = 10.0,
  gamma? : Double = 1.0,
) -> Array[Array[Double]] {
  let n = coords.length()
  let mat : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let d = coords[i].distance(coords[j])
      if d <= cutoff {
        mat[i][j] = -gamma
        mat[j][i] = -gamma
      }
      j = j + 1
    }
    i = i + 1
  }
  for i = 0; i < n; i = i + 1 {
    let mut degree = 0.0
    for j = 0; j < n; j = j + 1 {
      if i != j {
        degree = degree - mat[i][j]
      }
    }
    mat[i][i] = degree
  }
  mat
}
