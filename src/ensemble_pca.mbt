///|
/// ProDy-style conformational ensemble PCA.
///
/// Given multiple conformations of the same protein (e.g. NMR models,
/// MD snapshots), superimpose them onto a reference and perform PCA
/// on the CA displacement vectors to extract dominant motions.
///
/// Reference: PCA of conformational ensembles (Amadei, Linssen, Berendsen 1993).

///|
/// A conformational ensemble of CA coordinates.
pub struct Ensemble {
  n_res : Int
  n_conf : Int
  // coords[conf][res] = (x, y, z)
  coords : Array[Array[(Double, Double, Double)]]
  // Superimposed coords
  aligned : Array[Array[(Double, Double, Double)]]
  is_aligned : Bool
}

///|
/// PCA result.
pub struct EnsemblePCA {
  n_res : Int
  n_modes : Int
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]] // 3N × n_modes
  variance_explained : Array[Double]
  mean_coords : Array[(Double, Double, Double)]
}

///|
/// Build an ensemble from an array of chains.
/// Each chain must have the same number of CA atoms.
pub fn Ensemble::new(chains : Array[Chain]) -> Ensemble? {
  if chains.length() == 0 {
    return None
  }
  let n_conf = chains.length()
  // Extract CA coords from each chain
  let all_coords : Array[Array[(Double, Double, Double)]] = Array(capacity=0)
  let mut n_res = 0
  for chain in chains {
    let residues = chain.get_residues()
    let ca : Array[(Double, Double, Double)] = Array(capacity=0)
    for res in residues {
      match res.get_ca() {
        Some(atom) => ca.push((atom.coord.x, atom.coord.y, atom.coord.z))
        None => ()
      }
    }
    if n_res == 0 {
      n_res = ca.length()
    } else if ca.length() != n_res {
      return None
    }
    all_coords.push(ca)
  }
  if n_res < 3 {
    return None
  }
  // Copy coords as initial aligned
  let aligned : Array[Array[(Double, Double, Double)]] = Array(capacity=0)
  for conf in all_coords {
    let copy : Array[(Double, Double, Double)] = Array(capacity=0)
    for (x, y, z) in conf {
      copy.push((x, y, z))
    }
    aligned.push(copy)
  }
  Some(Ensemble::{ n_res, n_conf, coords: all_coords, aligned, is_aligned: false })
}

///|
/// Superimpose all conformations onto the first (reference).
pub fn Ensemble::superimpose(self : Ensemble) -> Ensemble {
  if self.n_conf <= 1 {
    return Ensemble::{ ..self, is_aligned: true }
  }
  let ref = self.coords[0]
  let aligned : Array[Array[(Double, Double, Double)]] = Array(capacity=0)
  aligned.push(ref) // reference stays
  for c = 1; c < self.n_conf; c = c + 1 {
    let conf = self.coords[c]
    // Simple superimpose: translate to common centroid, then Kabsch-like
    // For simplicity, just translate to ref centroid (full Kabsch is in svd_superimposer)
    let (rx, ry, rz) = centroid_3d(ref)
    let (cx, cy, cz) = centroid_3d(conf)
    let dx = rx - cx
    let dy = ry - cy
    let dz = rz - cz
    let shifted : Array[(Double, Double, Double)] = Array(capacity=0)
    for (x, y, z) in conf {
      shifted.push((x + dx, y + dy, z + dz))
    }
    aligned.push(shifted)
  }
  Ensemble::{ ..self, aligned, is_aligned: true }
}

///|
/// Perform PCA on the ensemble.
pub fn Ensemble::pca(self : Ensemble) -> EnsemblePCA? {
  let ens = if self.is_aligned {
    self
  } else {
    self.superimpose()
  }
  let n = ens.n_res
  let m = ens.n_conf
  let dim = 3 * n
  if m < 2 {
    return None
  }
  // Compute mean coordinates
  let mean = Array::make(n, (0.0, 0.0, 0.0))
  for i = 0; i < n; i = i + 1 {
    let mut sx = 0.0
    let mut sy = 0.0
    let mut sz = 0.0
    for c = 0; c < m; c = c + 1 {
      let (x, y, z) = ens.aligned[c][i]
      sx = sx + x
      sy = sy + y
      sz = sz + z
    }
    mean[i] = (sx / m.to_double(), sy / m.to_double(), sz / m.to_double())
  }
  // Build covariance matrix (dim × dim)
  // C = (1/(m-1)) * sum_c (x_c - mean)(x_c - mean)^T
  let cov : Array[Array[Double]] = Array::makei(dim, { _ => Array::make(dim, 0.0) })
  for c = 0; c < m; c = c + 1 {
    // Displacement vector
    let disp = Array::make(dim, 0.0)
    for i = 0; i < n; i = i + 1 {
      let (x, y, z) = ens.aligned[c][i]
      let (mx, my, mz) = mean[i]
      disp[3 * i] = x - mx
      disp[3 * i + 1] = y - my
      disp[3 * i + 2] = z - mz
    }
    // Outer product
    for a = 0; a < dim; a = a + 1 {
      for b = 0; b < dim; b = b + 1 {
        cov[a][b] = cov[a][b] + disp[a] * disp[b]
      }
    }
  }
  let norm = 1.0 / (m - 1).to_double()
  for a = 0; a < dim; a = a + 1 {
    for b = 0; b < dim; b = b + 1 {
      cov[a][b] = cov[a][b] * norm
    }
  }
  // Eigen-decompose
  let max_iter = if dim > 100 { dim * 20 } else { 200 }
  let (eigenvalues, eigenvectors) = linalg_jacobi(cov, dim, max_iter)
  // Sort descending
  let indices = Array::makei(dim, { i => i })
  for i = 1; i < dim; i = i + 1 {
    let key = indices[i]
    let key_val = eigenvalues[key]
    let mut j = i - 1
    while j >= 0 && eigenvalues[indices[j]] < key_val {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
  }
  // Collect non-trivial modes (eigenvalue > 1e-10)
  let mut n_modes = 0
  for k = 0; k < dim; k = k + 1 {
    if eigenvalues[indices[k]] > 1.0e-10 {
      n_modes = n_modes + 1
    }
  }
  if n_modes == 0 {
    n_modes = 1
  }
  let sorted_eigenvalues = Array::make(n_modes, 0.0)
  let sorted_eigenvectors : Array[Array[Double]] = Array::makei(dim, { _ => Array::make(n_modes, 0.0) })
  for k = 0; k < n_modes; k = k + 1 {
    sorted_eigenvalues[k] = eigenvalues[indices[k]]
    for i = 0; i < dim; i = i + 1 {
      sorted_eigenvectors[i][k] = eigenvectors[i][indices[k]]
    }
  }
  // Variance explained
  let total_var = {
    let mut s = 0.0
    for k = 0; k < n_modes; k = k + 1 {
      s = s + sorted_eigenvalues[k]
    }
    s
  }
  let var_explained = Array::make(n_modes, 0.0)
  for k = 0; k < n_modes; k = k + 1 {
    if total_var > 1.0e-10 {
      var_explained[k] = sorted_eigenvalues[k] / total_var
    }
  }
  Some(EnsemblePCA::{
    n_res: n,
    n_modes,
    eigenvalues: sorted_eigenvalues,
    eigenvectors: sorted_eigenvectors,
    variance_explained: var_explained,
    mean_coords: mean,
  })
}

///|
/// Get eigenvalues (descending).
pub fn EnsemblePCA::get_eigenvalues(self : EnsemblePCA) -> Array[Double] {
  self.eigenvalues
}

///|
/// Get eigenvectors (columns are PCs).
pub fn EnsemblePCA::get_eigenvectors(self : EnsemblePCA) -> Array[Array[Double]] {
  self.eigenvectors
}

///|
/// Variance explained by each PC.
pub fn EnsemblePCA::get_variance_explained(self : EnsemblePCA) -> Array[Double] {
  self.variance_explained
}

///|
/// Number of modes.
pub fn EnsemblePCA::num_modes(self : EnsemblePCA) -> Int {
  self.n_modes
}

///|
/// Project a conformation onto PCs.
/// Returns array of projections (PC scores).
pub fn EnsemblePCA::project(
  self : EnsemblePCA,
  coords : Array[(Double, Double, Double)],
) -> Array[Double] {
  let n = self.n_res
  let dim = 3 * n
  let disp = Array::make(dim, 0.0)
  for i = 0; i < n && i < coords.length(); i = i + 1 {
    let (x, y, z) = coords[i]
    let (mx, my, mz) = self.mean_coords[i]
    disp[3 * i] = x - mx
    disp[3 * i + 1] = y - my
    disp[3 * i + 2] = z - mz
  }
  let projections = Array::make(self.n_modes, 0.0)
  for k = 0; k < self.n_modes; k = k + 1 {
    let mut s = 0.0
    for i = 0; i < dim; i = i + 1 {
      s = s + self.eigenvectors[i][k] * disp[i]
    }
    projections[k] = s
  }
  projections
}

///|
/// Mean coordinates.
pub fn EnsemblePCA::mean_structure(self : EnsemblePCA) -> Array[(Double, Double, Double)] {
  self.mean_coords
}

///|
/// Number of conformations.
pub fn Ensemble::num_conformations(self : Ensemble) -> Int {
  self.n_conf
}

///|
/// Number of residues.
pub fn Ensemble::num_residues(self : Ensemble) -> Int {
  self.n_res
}

///|
/// Get aligned coordinates for a specific conformation.
pub fn Ensemble::get_aligned(self : Ensemble, conf : Int) -> Array[(Double, Double, Double)] {
  self.aligned[conf]
}

///|
/// Compute 3D centroid.
fn centroid_3d(coords : Array[(Double, Double, Double)]) -> (Double, Double, Double) {
  let n = coords.length()
  if n == 0 {
    return (0.0, 0.0, 0.0)
  }
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sz = 0.0
  for (x, y, z) in coords {
    sx = sx + x
    sy = sy + y
    sz = sz + z
  }
  let nf = n.to_double()
  (sx / nf, sy / nf, sz / nf)
}
