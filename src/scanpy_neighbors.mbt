///|
/// scanpy `pp.neighbors` — k-nearest-neighbour graph and UMAP connectivities
/// (`method='umap'`), computing the fuzzy simplicial set from Euclidean
/// distances. The connectivities are returned as a `SparseGraph` (O(n·k))
/// instead of a dense `n × n` matrix.

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
/// distance 0). Computed row-by-row so no `n × n` distance matrix is ever
/// materialized: peak transient memory is `O(n)` per row instead of `O(n²)`.
fn scanpy_knn_distances(
  x : Array[Array[Double]],
  k : Int,
) -> (Array[Array[Int]], Array[Array[Double]]) {
  let n = x.length()
  let idx : Array[Array[Int]] = Array(capacity=n)
  let d : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      let mut acc = 0.0
      for dim = 0; dim < x[i].length(); dim = dim + 1 {
        let diff = x[i][dim] - x[j][dim]
        acc = acc + diff * diff
      }
      row.push(acc.sqrt())
    }
    let sorted = scanpy_argsort(row)
    let idx_row : Array[Int] = Array(capacity=k)
    let d_row : Array[Double] = Array(capacity=k)
    for t = 0; t < k; t = t + 1 {
      idx_row.push(sorted[t])
      d_row.push(row[sorted[t]])
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
/// `local_connectivity=1`, `set_op_mix_ratio=1`) and return them as a sparse
/// symmetric graph (O(n·k) edges). Diagonal self-weights are 1.
pub fn scanpy_neighbors(x : Array[Array[Double]], n_neighbors : Int) -> SparseGraph {
  let n = x.length()
  if n == 0 {
    return SparseGraph::{ n: 0, diag: [], adj: [] }
  }
  let k = if n_neighbors < n { n_neighbors } else { n }
  let (knn_idx, knn_dists) = scanpy_knn_distances(x, k)
  let (rhos, sigmas) = scanpy_smooth_knn_dist(knn_dists, k)

  // directed membership strengths d (self = 1), stored per-vertex as Map[j -> w]
  let directed : Array[Map[Int, Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Map[Int, Double] = Map([], capacity=k)
    for j = 1; j < k; j = j + 1 {
      let nb = knn_idx[i][j]
      let delta = knn_dists[i][j] - rhos[i]
      let w = if delta > 0.0 {
        @math.exp(-(delta / sigmas[i]))
      } else {
        1.0
      }
      row.set(nb, w)
    }
    directed.push(row)
  }

  // fuzzy union conn = d + dᵀ - d ∘ dᵀ, deduplicated undirected edges by a
  // canonical Int key (a * n + b with a < b).
  let undirected : Map[Int, Double] = Map([], capacity=n * k)
  for i = 0; i < n; i = i + 1 {
    for nb in directed[i] {
      let j = nb.0
      if j == i {
        continue
      }
      let a = if i < j { i } else { j }
      let b = if i < j { j } else { i }
      let key = a * n + b
      if undirected.contains(key) {
        continue
      }
      let w_ij = nb.1
      let w_ji = match directed[j].get(i) {
        None => 0.0
        Some(x) => x
      }
      let conn = w_ij + w_ji - w_ij * w_ji
      if conn > 0.0 {
        undirected.set(key, conn)
      }
    }
  }

  // diagonal self-weight = 1 (fuzzy union of a self-loop with itself)
  let diag : Array[Double] = Array::make(n, 1.0)
  let adj : Array[Array[(Int, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    adj.push([])
  }
  for e in undirected {
    let a = e.0 / n
    let b = e.0 % n
    let w = e.1
    adj[a].push((b, w))
    adj[b].push((a, w))
  }
  SparseGraph::{ n, diag, adj }
}