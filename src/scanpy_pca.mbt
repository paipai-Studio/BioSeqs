///|
/// scanpy `pp.pca` — principal component analysis of an expression matrix.
/// Implements the default sklearn convention (zero-center columns, eigendecompose
/// the covariance) to produce cell coordinates, gene loadings and the explained
/// variance ratio.

///|
pub struct ScanpyPca {
  coords : Array[Array[Double]]
  loadings : Array[Array[Double]]
  variance_ratio : Array[Double]
}

///|
/// Principal component analysis of a dense `n_obs × n_vars` matrix (zero-centred
/// columns), matching `scanpy.pp.pca`'s default (sklearn PCA). Returns cell
/// coordinates (`n_obs × n_comps`), gene loadings (`n_comps × n_vars`) and the
/// explained-variance ratio per component.
pub fn scanpy_pca(x : Array[Array[Double]], n_comps : Int) -> ScanpyPca {
  let n_obs = x.length()
  let n_vars = if n_obs > 0 { x[0].length() } else { 0 }
  if n_obs == 0 || n_vars == 0 {
    return ScanpyPca::{ coords: [], loadings: [], variance_ratio: [] }
  }

  // column means
  let means = Array::make(n_vars, 0.0)
  for j = 0; j < n_vars; j = j + 1 {
    for i = 0; i < n_obs; i = i + 1 {
      means[j] = means[j] + x[i][j]
    }
    means[j] = means[j] / n_obs.to_double()
  }

  // centred matrix + covariance (Xcᵀ Xc), eigendecomposed (n_vars × n_vars)
  let cov : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < n_vars; j = j + 1 {
    let row = Array::make(n_vars, 0.0)
    for k = 0; k < n_vars; k = k + 1 {
      let mut acc = 0.0
      for i = 0; i < n_obs; i = i + 1 {
        acc = acc + (x[i][j] - means[j]) * (x[i][k] - means[k])
      }
      row[k] = acc
    }
    cov.push(row)
  }

  let (eigenvalues, eigenvectors) = linalg_jacobi(cov, n_vars, 100 * n_vars * n_vars)
  let (evals, evecs) = linalg_sort_eigen(eigenvalues, eigenvectors, n_vars)

  let nc = if n_comps < n_vars { n_comps } else { n_vars }

  let mut var_sum = 0.0
  for c = 0; c < nc; c = c + 1 {
    if evals[c] > 0.0 {
      var_sum = var_sum + evals[c]
    }
  }

  let ratio : Array[Double] = Array(capacity=0)
  for c = 0; c < nc; c = c + 1 {
    ratio.push(if var_sum > 0.0 { evals[c] / var_sum } else { 0.0 })
  }

  let loadings : Array[Array[Double]] = Array(capacity=0)
  for c = 0; c < nc; c = c + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < n_vars; j = j + 1 {
      row.push(evecs[j][c])
    }
    loadings.push(row)
  }

  let coords : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n_obs; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for c = 0; c < nc; c = c + 1 {
      let mut acc = 0.0
      for j = 0; j < n_vars; j = j + 1 {
        acc = acc + (x[i][j] - means[j]) * evecs[j][c]
      }
      row.push(acc)
    }
    coords.push(row)
  }

  ScanpyPca::{ coords, loadings, variance_ratio: ratio }
}