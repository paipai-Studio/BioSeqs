///|
/// scanpy `pp.neighbors` with pluggable distance metrics.  The core pipeline
/// (kNN -> smooth-knn bandwidth -> fuzzy simplicial set, or Gaussian kernel)
/// is metric-agnostic: `scanpy_knn_distances_metric` dispatches the pairwise
/// distance and the same UMAP/Gauss kernels are applied to the resulting kNN
/// graph.  Supported metrics: `euclidean` (kd-tree accelerated, identical to
/// `scanpy_neighbors_with_distances`), `cosine`, `correlation` (Pearson,
/// computed as cosine on mean-centred vectors) and `manhattan`.

///|
/// Pairwise distance between two vectors under a scanpy `pp.neighbors`
/// `metric` choice.  Zero-norm vectors under `cosine`/`correlation` yield the
/// maximal distance 1.0 (scikit-learn convention).
fn scanpy_pair_distance(
  a : Array[Double],
  b : Array[Double],
  metric : String,
) -> Double {
  let d = if a.length() < b.length() { a.length() } else { b.length() }
  match metric {
    "manhattan" => {
      let mut s = 0.0
      for t = 0; t < d; t = t + 1 {
        let diff = a[t] - b[t]
        s = s + (if diff < 0.0 { -diff } else { diff })
      }
      s
    }
    "cosine" => {
      let mut dot = 0.0
      let mut na = 0.0
      let mut nb = 0.0
      for t = 0; t < d; t = t + 1 {
        dot = dot + a[t] * b[t]
        na = na + a[t] * a[t]
        nb = nb + b[t] * b[t]
      }
      let denom = na.sqrt() * nb.sqrt()
      if denom <= 0.0 {
        1.0
      } else {
        let c = dot / denom
        if c > 1.0 {
          0.0
        } else if c < 0.0 {
          1.0
        } else {
          1.0 - c
        }
      }
    }
    "correlation" => {
      // cosine of the mean-centred vectors (Pearson distance)
      let mut ma = 0.0
      let mut mb = 0.0
      for t = 0; t < d; t = t + 1 {
        ma = ma + a[t]
        mb = mb + b[t]
      }
      ma = ma / d.to_double()
      mb = mb / d.to_double()
      let mut dot = 0.0
      let mut na = 0.0
      let mut nb = 0.0
      for t = 0; t < d; t = t + 1 {
        let da = a[t] - ma
        let db = b[t] - mb
        dot = dot + da * db
        na = na + da * da
        nb = nb + db * db
      }
      let denom = na.sqrt() * nb.sqrt()
      if denom <= 0.0 {
        1.0
      } else {
        let c = dot / denom
        if c > 1.0 {
          0.0
        } else if c < 0.0 {
          1.0
        } else {
          1.0 - c
        }
      }
    }
    _ => {
      let mut s = 0.0
      for t = 0; t < d; t = t + 1 {
        let diff = a[t] - b[t]
        s = s + diff * diff
      }
      s.sqrt()
    }
  }
}

///|
/// kNN indices and distances under an arbitrary metric.  `euclidean` reuses
/// the shared kd-tree/brute-force dispatch (`scanpy_knn_distances`); every
/// other metric uses the brute-force O(n^2) kernel with the generic distance
/// (row-by-row, no dense matrix materialized).  Rows are sorted ascending by
/// distance, index 0 is the self point at distance 0.
fn scanpy_knn_distances_metric(
  x : Array[Array[Double]],
  k : Int32,
  metric : String,
) -> (Array[Array[Int32]], Array[Array[Double]]) {
  if metric == "euclidean" {
    return scanpy_knn_distances(x, k)
  }
  let n = x.length()
  if n == 0 {
    return (Array(capacity=0), Array(capacity=0))
  }
  let idx : Array[Array[Int32]] = Array(capacity=n)
  let dists : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      row.push(scanpy_pair_distance(x[i], x[j], metric))
    }
    let sorted = scanpy_argsort(row)
    let idx_row : Array[Int32] = Array(capacity=k)
    let d_row : Array[Double] = Array(capacity=k)
    for t = 0; t < k; t = t + 1 {
      idx_row.push(sorted[t])
      d_row.push(row[sorted[t]])
    }
    idx.push(idx_row)
    dists.push(d_row)
  }
  (idx, dists)
}

///|
/// Symmetrized kNN distance graph (scanpy `obsp["distances"]`): edge weight is
/// the minimum of the two directed kNN distances, 0 diagonal.
fn scanpy_dist_graph_from_knn(
  knn_idx : Array[Array[Int32]],
  knn_dists : Array[Array[Double]],
  k : Int32,
  n : Int32,
) -> SparseGraph {
  let edges : Map[Int32, Double] = Map([], capacity=n * k)
  for i = 0; i < n; i = i + 1 {
    for j = 1; j < k; j = j + 1 {
      let nb = knn_idx[i][j]
      let d = knn_dists[i][j]
      let a = if i < nb { i } else { nb }
      let b = if i < nb { nb } else { i }
      let key = a * n + b
      match edges.get(key) {
        None => edges.set(key, d)
        Some(existing) => if d < existing { edges.set(key, d) }
      }
    }
  }
  let diag : Array[Double] = Array::make(n, 0.0)
  let adj : Array[Array[(Int32, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    adj.push([])
  }
  for e in edges {
    let a = e.0 / n
    let b = e.0 % n
    adj[a].push((b, e.1))
    adj[b].push((a, e.1))
  }
  SparseGraph::{ n, diag, adj, }
}

///|
/// UMAP fuzzy simplicial set connectivities from a kNN graph
/// (scanpy `obsp["connectivities"]`, `method='umap'`).  Diagonal self-weight
/// is 1; off-diagonal weights come from the fuzzy union
/// `w_ij + w_ji - w_ij * w_ji` of the directed membership strengths.
fn scanpy_fuzzy_connectivities(
  knn_idx : Array[Array[Int32]],
  knn_dists : Array[Array[Double]],
  k : Int32,
  n : Int32,
) -> SparseGraph {
  let (rhos, sigmas) = scanpy_smooth_knn_dist(knn_dists, k)
  let directed : Array[Map[Int32, Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Map[Int32, Double] = Map([], capacity=k)
    for j = 1; j < k; j = j + 1 {
      let nb = knn_idx[i][j]
      let delta = knn_dists[i][j] - rhos[i]
      let w = if delta > 0.0 { @math.exp(-(delta / sigmas[i])) } else { 1.0 }
      row.set(nb, w)
    }
    directed.push(row)
  }
  let undirected : Map[Int32, Double] = Map([], capacity=n * k)
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
  let diag : Array[Double] = Array::make(n, 1.0)
  let adj : Array[Array[(Int32, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    adj.push([])
  }
  for e in undirected {
    let a = e.0 / n
    let b = e.0 % n
    adj[a].push((b, e.1))
    adj[b].push((a, e.1))
  }
  SparseGraph::{ n, diag, adj, }
}

///|
/// Gaussian-kernel connectivities from a kNN graph (scanpy
/// `pp.neighbors(method='gauss')`): `exp(-d^2 / (2 * sigma^2))` with the
/// global bandwidth `sigma` = mean k-th neighbour distance.  Diagonal = 1.
fn scanpy_gauss_connectivities(
  knn_idx : Array[Array[Int32]],
  knn_dists : Array[Array[Double]],
  k : Int32,
  n : Int32,
) -> SparseGraph {
  let mut sigma_sum = 0.0
  for i = 0; i < n; i = i + 1 {
    sigma_sum = sigma_sum + knn_dists[i][k - 1]
  }
  let sigma = if n > 0 { sigma_sum / n.to_double() } else { 1.0 }
  let sigma_sq = if sigma > 0.0 { 2.0 * sigma * sigma } else { 1.0 }
  let edges : Map[Int32, Double] = Map([], capacity=n * k)
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
      match edges.get(key) {
        None => if w > 0.0 { edges.set(key, w) }
        Some(existing) => if w > existing { edges.set(key, w) }
      }
    }
  }
  let diag : Array[Double] = Array::make(n, 1.0)
  let adj : Array[Array[(Int32, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    adj.push([])
  }
  for e in edges {
    let a = e.0 / n
    let b = e.0 % n
    adj[a].push((b, e.1))
    adj[b].push((a, e.1))
  }
  SparseGraph::{ n, diag, adj, }
}

///|
/// Neighbors under an arbitrary distance metric (`pp.neighbors(metric=...)`).
/// Returns both the symmetrized kNN distance graph and the connectivities
/// graph.  For `metric="euclidean"` and `method="umap"` the output is
/// identical to `scanpy_neighbors_with_distances` (same kd-tree path, same
/// tie-breaking).  `method="gauss"` selects the Gaussian kernel instead of
/// the UMAP fuzzy simplicial set.
pub fn scanpy_neighbors_metric(
  x : Array[Array[Double]],
  n_neighbors : Int32,
  metric~ : String = "euclidean",
  method_~ : String = "umap",
) -> ScanpyNeighborsResult {
  let n = x.length()
  if n == 0 {
    let empty = SparseGraph::{ n: 0, diag: [], adj: [], }
    return ScanpyNeighborsResult::{ distances: empty, connectivities: empty, }
  }
  let k = if n_neighbors < n { n_neighbors } else { n }
  let (knn_idx, knn_dists) = scanpy_knn_distances_metric(x, k, metric)
  let distances = scanpy_dist_graph_from_knn(knn_idx, knn_dists, k, n)
  let connectivities = match method_ {
    "gauss" => scanpy_gauss_connectivities(knn_idx, knn_dists, k, n)
    _ => scanpy_fuzzy_connectivities(knn_idx, knn_dists, k, n)
  }
  ScanpyNeighborsResult::{ distances, connectivities, }
}
