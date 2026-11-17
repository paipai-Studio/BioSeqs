///|
/// scanpy-referenced imputation: MAGIC (Markov affinity-based graph
/// imputation of cells). Builds an adaptive Gaussian-kernel kNN affinity over
/// cells, forms the symmetric diffusion operator T = D^-1/2 W D^-1/2 and
/// diffuses the expression matrix `T^t X` (van Dijk et al., 2018; the
/// `scanpy.external.pp.magic` counterpart).

///|
/// MAGIC imputation of the expression matrix. `t` is the number of diffusion
/// steps, `k` the neighbor count for the kernel bandwidth (sigma_i = distance
/// to the k-th nearest neighbor). The result is stored back onto `adata` —
/// as the layer named `layer` when non-empty, otherwise replacing `adata.x` —
/// and also returned.
pub fn scanpy_magic(
  adata : AnnData,
  t~ : Int32 = 3,
  k~ : Int32 = 10,
  layer~ : String = "",
) -> Array[Array[Double]] {
  let x = if layer.length() > 0 && adata.layers.contains(layer) {
    adata.layers[layer]
  } else {
    adata.x_matrix()
  }
  let result = scanpy_magic_matrix(x, t~, k~)
  if layer.length() > 0 {
    let _ = adata.set_layer(layer, result)
  } else {
    adata.x = result
  }
  result
}

///|
/// MAGIC on a plain cell × gene matrix (no AnnData involvement).
pub fn scanpy_magic_matrix(
  x : Array[Array[Double]],
  t~ : Int32 = 3,
  k~ : Int32 = 10,
) -> Array[Array[Double]] {
  let n = x.length()
  if n == 0 {
    return []
  }
  let m = x[0].length()
  if n == 1 || m == 0 {
    return x.map(fn(row) { row.copy() })
  }
  let kk = if k < 1 { 1 } else if k > n - 1 { n - 1 } else { k }
  let steps = if t < 0 { 0 } else { t }
  // pairwise squared distances
  let d2 : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      let mut s = 0.0
      for g = 0; g < m; g = g + 1 {
        let diff = x[i][g] - x[j][g]
        s = s + diff * diff
      }
      row.push(s)
    }
    d2.push(row)
  }
  // adaptive bandwidth: k-th nearest neighbor distance (excluding self)
  let sigma = Array::make(n, 1.0)
  for i = 0; i < n; i = i + 1 {
    let mut best = @double.infinity
    // find k-th smallest via partial selection over d2[i]
    let dists : Array[Double] = Array(capacity=n - 1)
    for j = 0; j < n; j = j + 1 {
      if j != i {
        dists.push(d2[i][j].sqrt())
      }
    }
    dists.sort()
    let idx = if kk - 1 < dists.length() { kk - 1 } else { dists.length() - 1 }
    if idx >= 0 && idx < dists.length() {
      best = dists[idx]
    }
    sigma[i] = if best > 0.0 && best != @double.infinity { best } else { 1.0 }
  }
  // adaptive Gaussian kernel affinity (asymmetric, then symmetrized)
  let w : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    let s2 = sigma[i] * sigma[i]
    for j = 0; j < n; j = j + 1 {
      if i != j {
        row[j] = @math.pow(2.718281828459045, -d2[i][j] / s2)
      }
    }
    w.push(row)
  }
  // symmetrize and compute degrees
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < i; j = j + 1 {
      let v = (w[i][j] + w[j][i]) / 2.0
      w[i][j] = v
      w[j][i] = v
    }
  }
  let deg = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < n; j = j + 1 {
      s = s + w[i][j]
    }
    deg[i] = s
  }
  // diffuse: Y = T^t X with T = D^-1/2 W D^-1/2
  let mut y : Array[Array[Double]] = x.map(fn(row) { row.copy() })
  for _step in 0..<steps {
    let out : Array[Array[Double]] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      let row = Array::make(m, 0.0)
      let di = deg[i]
      for j = 0; j < n; j = j + 1 {
        let dj = deg[j]
        let tij = if di > 0.0 && dj > 0.0 {
          w[i][j] / (di * dj).sqrt()
        } else {
          0.0
        }
        if tij != 0.0 {
          for g = 0; g < m; g = g + 1 {
            row[g] = row[g] + tij * y[j][g]
          }
        }
      }
      out.push(row)
    }
    y = out
  }
  y
}
