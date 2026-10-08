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
/// distance 0).
///
/// **Dispatch**: uses a kd-tree when `n > 64` and `d ≤ 20` (typical
/// PCA-reduced single-cell data), giving O(n·log n) total query time.
/// Otherwise falls back to brute-force O(n²).  Both paths produce
/// bit-identical output (same distances, same tie-breaking by index
/// ascending).
fn scanpy_knn_distances(
  x : Array[Array[Double]],
  k : Int,
) -> (Array[Array[Int]], Array[Array[Double]]) {
  let n = x.length()
  if n == 0 {
    return (Array(capacity=0), Array(capacity=0))
  }
  let d = x[0].length()
  // kd-tree: beneficial for low-dimensional data with enough points
  if n > 64 && d > 0 && d <= 20 {
    return scanpy_knn_distances_kdtree(x, k, n, d)
  }
  // Brute-force fallback (high-dim or small-n)
  scanpy_knn_distances_brute(x, k)
}

///|
/// Brute-force kNN: O(n²) distance computation, row-by-row (no n×n matrix
/// materialized).  Used as fallback when kd-tree is not beneficial
/// (d > 20 or n ≤ 64).
fn scanpy_knn_distances_brute(
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

// ============================================================================
// Extended neighbors API
//
// scanpy's `sc.pp.neighbors()` stores two matrices in `adata.obsp`:
//   - 'distances'      — the raw kNN distance matrix (sparse, symmetric)
//   - 'connectivities' — the UMAP fuzzy simplicial set (sparse, symmetric)
//
// The original `scanpy_neighbors` returns only connectivities.  The functions
// below add:
//   - scanpy_neighbors_with_distances  — returns both as a pair
//   - scanpy_neighbors_gauss           — Gaussian kernel (scanpy method='gauss')
// ============================================================================

///|
/// Result of `scanpy_neighbors_with_distances`: the kNN distance graph and
/// the UMAP connectivities graph, both as sparse symmetric matrices.
pub struct ScanpyNeighborsResult {
  distances : SparseGraph       // raw kNN distances (symmetric, 0 on diagonal)
  connectivities : SparseGraph  // UMAP fuzzy simplicial set (1 on diagonal)
}

///|
/// Compute both the kNN distance graph and the UMAP connectivities graph.
///
/// `distances` is the symmetrized kNN distance matrix: `dist[i,j]` is the
/// distance from i to its nearest neighbor j (0 if j is not among i's k
/// nearest).  The matrix is symmetrized by taking the minimum of the two
/// directions.
///
/// `connectivities` is the same UMAP fuzzy simplicial set returned by
/// `scanpy_neighbors`.
pub fn scanpy_neighbors_with_distances(
  x : Array[Array[Double]],
  n_neighbors : Int,
) -> ScanpyNeighborsResult {
  let n = x.length()
  if n == 0 {
    return ScanpyNeighborsResult::{
      distances: SparseGraph::{ n: 0, diag: [], adj: [] },
      connectivities: SparseGraph::{ n: 0, diag: [], adj: [] },
    }
  }
  let k = if n_neighbors < n { n_neighbors } else { n }
  let (knn_idx, knn_dists) = scanpy_knn_distances(x, k)
  let (rhos, sigmas) = scanpy_smooth_knn_dist(knn_dists, k)

  // --- Build distances graph (symmetrized by min) ---
  // For each vertex i, store kNN neighbors with distances. Symmetrize by
  // taking the minimum distance for each undirected edge.
  let dist_edges : Map[Int, Double] = Map([], capacity=n * k)
  for i = 0; i < n; i = i + 1 {
    for j = 1; j < k; j = j + 1 {
      let nb = knn_idx[i][j]
      let d = knn_dists[i][j]
      let a = if i < nb { i } else { nb }
      let b = if i < nb { nb } else { i }
      let key = a * n + b
      match dist_edges.get(key) {
        None => dist_edges.set(key, d)
        Some(existing) => {
          if d < existing {
            dist_edges.set(key, d)
          }
        }
      }
    }
  }
  let dist_diag : Array[Double] = Array::make(n, 0.0)
  let dist_adj : Array[Array[(Int, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    dist_adj.push([])
  }
  for e in dist_edges {
    let a = e.0 / n
    let b = e.0 % n
    let w = e.1
    dist_adj[a].push((b, w))
    dist_adj[b].push((a, w))
  }
  let distances = SparseGraph::{ n, diag: dist_diag, adj: dist_adj }

  // --- Build connectivities graph (UMAP fuzzy union) ---
  // directed membership strengths
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
  // fuzzy union: conn = w_ij + w_ji - w_ij * w_ji
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
  let conn_diag : Array[Double] = Array::make(n, 1.0)
  let conn_adj : Array[Array[(Int, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    conn_adj.push([])
  }
  for e in undirected {
    let a = e.0 / n
    let b = e.0 % n
    let w = e.1
    conn_adj[a].push((b, w))
    conn_adj[b].push((a, w))
  }
  let connectivities = SparseGraph::{ n, diag: conn_diag, adj: conn_adj }

  ScanpyNeighborsResult::{ distances, connectivities }
}

///|
/// Compute neighbors using a Gaussian kernel (scanpy's `method='gauss'`).
///
/// The Gaussian connectivities are: `conn[i,j] = exp(-dist² / (2 * sigma²))`
/// where `sigma` is the mean distance to the k-th nearest neighbor across
/// all points.  This is simpler than the UMAP kernel and is scanpy's
/// alternative to `method='umap'`.
pub fn scanpy_neighbors_gauss(
  x : Array[Array[Double]],
  n_neighbors : Int,
) -> SparseGraph {
  let n = x.length()
  if n == 0 {
    return SparseGraph::{ n: 0, diag: [], adj: [] }
  }
  let k = if n_neighbors < n { n_neighbors } else { n }
  let (knn_idx, knn_dists) = scanpy_knn_distances(x, k)

  // sigma = mean of k-th neighbor distances (per-point bandwidth)
  let mut sigma_sum = 0.0
  for i = 0; i < n; i = i + 1 {
    sigma_sum = sigma_sum + knn_dists[i][k - 1]
  }
  let sigma = if n > 0 { sigma_sum / n.to_double() } else { 1.0 }
  // Guard against sigma == 0 (all points identical)
  let sigma_sq = if sigma > 0.0 { 2.0 * sigma * sigma } else { 1.0 }

  // Build connectivities: exp(-dist² / (2*sigma²))
  let edges : Map[Int, Double] = Map([], capacity=n * k)
  for i = 0; i < n; i = i + 1 {
    for j = 1; j < k; j = j + 1 {
      let nb = knn_idx[i][j]
      if nb == i {
        continue
      }
      let d = knn_dists[i][j]
      let w = @math.exp(-(d * d / sigma_sq))
      let a = if i < nb { i } else { nb }
      let b = if i < nb { nb } else { i }
      let key = a * n + b
      // symmetrize by taking the max (both directions should give same w
      // since distance is symmetric, but kNN membership may not be)
      match edges.get(key) {
        None => {
          if w > 0.0 {
            edges.set(key, w)
          }
        }
        Some(existing) => {
          if w > existing {
            edges.set(key, w)
          }
        }
      }
    }
  }

  let diag : Array[Double] = Array::make(n, 1.0)
  let adj : Array[Array[(Int, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    adj.push([])
  }
  for e in edges {
    let a = e.0 / n
    let b = e.0 % n
    let w = e.1
    adj[a].push((b, w))
    adj[b].push((a, w))
  }
  SparseGraph::{ n, diag, adj }
}
