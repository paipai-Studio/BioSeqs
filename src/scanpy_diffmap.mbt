///|
/// scanpy `tl.diffmap` — diffusion map embedding of a symmetric weighted
/// adjacency (the `pp.neighbors` connectivities). Computes the symmetrized
/// transition `D⁻½ A D⁻½`, eigendecomposes it, and returns the top diffusion
/// coordinates (eigenvectors) and eigenvalues. This basis underlies `tl.dpt`
/// (diffusion pseudotime).

///|
pub struct ScanpyDiffmap {
  eigen_basis : Array[Array[Double]]
  eigen_values : Array[Double]
}

///|
pub fn scanpy_diffmap(
  graph : SparseGraph,
  n_comps : Int,
) -> ScanpyDiffmap {
  let n = graph.num_vertices()
  if n == 0 {
    return ScanpyDiffmap::{ eigen_basis: [], eigen_values: [] }
  }
  // degrees (diagonal self-weight + incident edge weights)
  let deg = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    deg[i] = graph.degree(i)
  }
  // symmetrized transition P_sym = D^-1/2 A D^-1/2 (densified for the eigensolver)
  let psym : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    let di = deg[i]
    if di > 0.0 {
      row[i] = graph.diag[i] / di
    }
    psym.push(row)
  }
  for i = 0; i < n; i = i + 1 {
    let di = deg[i]
    for nb in graph.adj[i] {
      let j = nb.0
      let dj = deg[j]
      if di > 0.0 && dj > 0.0 {
        psym[i][j] = nb.1 / (di * dj).sqrt()
      }
    }
  }
  let (eigenvalues, eigenvectors) = linalg_jacobi(psym, n, 100 * n * n)
  let (evals, evecs) = linalg_sort_eigen(eigenvalues, eigenvectors, n)

  let nc = if n_comps < n { n_comps } else { n }
  let basis : Array[Array[Double]] = Array(capacity=0)
  let vals : Array[Double] = Array(capacity=0)
  for t = 0; t < nc; t = t + 1 {
    vals.push(evals[t])
  }
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for t = 0; t < nc; t = t + 1 {
      row.push(evecs[i][t])
    }
    basis.push(row)
  }
  ScanpyDiffmap::{ eigen_basis: basis, eigen_values: vals }
}
///|
/// Diffusion pseudotime from a rooted diffusion map: the cumulative diffusion
/// distance from the root cell, `Σ_t λ_t/(1−λ_t)·(ψ_t(i)−ψ_t(root))²` over the
/// non-trivial diffusion components (`t ≥ 1`).
pub fn scanpy_dpt(diffmap : ScanpyDiffmap, root : Int) -> Array[Double] {
  let n = diffmap.eigen_basis.length()
  let nc = diffmap.eigen_values.length()
  let out = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut acc = 0.0
    for t = 1; t < nc; t = t + 1 {
      let lambda = diffmap.eigen_values[t]
      let denom = 1.0 - lambda
      if denom > 1.0e-12 {
        let diff = diffmap.eigen_basis[i][t] - diffmap.eigen_basis[root][t]
        acc = acc + (lambda / denom) * diff * diff
      }
    }
    out[i] = acc
  }
  out
}
