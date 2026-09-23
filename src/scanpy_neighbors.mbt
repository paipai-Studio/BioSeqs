///|
/// scanpy `pp.neighbors` — k-nearest-neighbour graph and UMAP connectivities
/// (`method='umap'`), computing the fuzzy simplicial set from Euclidean
/// distances on a dense matrix.

///|
/// Pairwise Euclidean distance matrix (n × n).
fn scanpy_euclidean_distances(x : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = x.length()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      let mut acc = 0.0
      for d = 0; d < x[i].length(); d = d + 1 {
        let diff = x[i][d] - x[j][d]
        acc = acc + diff * diff
      }
      row.push(acc.sqrt())
    }
    out.push(row)
  }
  out
}

///|
/// Arg-sort a distance row ascending (indices, not values).
fn scanpy_argsort(dist_row : Array[Double]) -> Array[Int] {
  let n = dist_row.length()
  let idx = Array::make(n, 0)
  for t = 0; t < n; t = t + 1 {
    idx[t] = t
  }
  let mut i = 1
  while i < n {
    let key = idx[i]
    let kv = dist_row[key]
    let mut j = i - 1
    while j >= 0 && dist_row[idx[j]] > kv {
      idx[j + 1] = idx[j]
      j = j - 1
    }
    idx[j + 1] = key
    i = i + 1
  }
  idx
}

///|
/// kNN indices and distances (each row sorted ascending; index 0 is self at
/// distance 0).
fn scanpy_knn(
  dist : Array[Array[Double]],
  k : Int,
) -> (Array[Array[Int]], Array[Array[Double]]) {
  let n = dist.length()
  let idx : Array[Array[Int]] = Array(capacity=0)
  let d : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let sorted = scanpy_argsort(dist[i])
    let idx_row : Array[Int] = Array(capacity=0)
    let d_row : Array[Double] = Array(capacity=0)
    for t = 0; t < k; t = t + 1 {
      idx_row.push(sorted[t])
      d_row.push(dist[i][sorted[t]])
    }
    idx.push(idx_row)
    d.push(d_row)
  }
  (idx, d)
}

///|
/// UMAP `smooth_knn_dist`: estimate per-point `rho` (distance to the nearest
/// neighbour) and `sigma` (binary-searched bandwidth so the log2(k) perplexity
/// is satisfied).
fn scanpy_smooth_knn_dist(
  knn_dists : Array[Array[Double]],
  k : Int,
) -> (Array[Double], Array[Double]) {
  let n = knn_dists.length()
  let target = @math.log2(k.to_double())
  let rhos : Array[Double] = Array(capacity=0)
  let sigmas : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let rho = knn_dists[i][1]
    let mut lo = 0.0
    let mut hi = 1.0e30
    let mut mid = 1.0
    for _iter = 0; _iter < 64; _iter = _iter + 1 {
      let mut psum = 0.0
      for j = 1; j < k; j = j + 1 {
        let d = knn_dists[i][j] - rho
        if d > 0.0 {
          psum = psum + @math.exp(-(d / mid))
        } else {
          psum = psum + 1.0
        }
      }
      if psum > target {
        hi = mid
        mid = (lo + hi) / 2.0
      } else {
        lo = mid
        mid = if hi > 1.0e29 { mid * 2.0 } else { (lo + hi) / 2.0 }
      }
    }
    rhos.push(rho)
    sigmas.push(mid)
  }
  (rhos, sigmas)
}

///|
/// Compute the UMAP neighborhood connectivities (fuzzy simplicial set with
/// `local_connectivity=1`, `set_op_mix_ratio=1`). Returns a symmetric dense
/// `n_obs × n_obs` matrix whose `[i, j]` entry is the fuzzy union of the
/// membership strengths between cells `i` and `j`.
pub fn scanpy_neighbors(x : Array[Array[Double]], n_neighbors : Int) -> Array[Array[Double]] {
  let n = x.length()
  if n == 0 {
    return []
  }
  let k = if n_neighbors < n { n_neighbors } else { n }
  let dist = scanpy_euclidean_distances(x)
  let (knn_idx, knn_dists) = scanpy_knn(dist, k)
  let (rhos, sigmas) = scanpy_smooth_knn_dist(knn_dists, k)

  // directed membership strengths d (self = 1)
  let d : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    row[i] = 1.0
    d.push(row)
  }
  for i = 0; i < n; i = i + 1 {
    for j = 1; j < k; j = j + 1 {
      let nb = knn_idx[i][j]
      let delta = knn_dists[i][j] - rhos[i]
      let w = if delta > 0.0 {
        @math.exp(-(delta / sigmas[i]))
      } else {
        1.0
      }
      d[i][nb] = w
    }
  }

  // fuzzy union: conn = d + dᵀ - d ∘ dᵀ
  let conn : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      row.push(d[i][j] + d[j][i] - d[i][j] * d[j][i])
    }
    conn.push(row)
  }
  conn
}