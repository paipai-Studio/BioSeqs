///|
/// scanpy `tl.kmeans` — K-Means clustering with k-means++ initialization and
/// Lloyd iterations, brought into the AnnData workflow.  The matrix-level
/// entry point returns labels, centers, inertia and the iteration count; the
/// AnnData wrapper clusters the expression matrix and stores the labels as a
/// string column in `obs` (scanpy stores cluster labels as categoricals).

///|
/// K-Means clustering result: per-point labels (dense `0..k-1`), cluster
/// centroids, within-cluster sum of squared distances (inertia) and the
/// number of Lloyd iterations of the best run.
pub struct ScanpyKMeans {
  labels : Array[Int32]
  centers : Array[Array[Double]]
  inertia : Double
  n_iter : Int32
  n_clusters : Int32
}

///|
/// Small deterministic LCG used for seeding and tie-breaking (same
/// parameters as `scanpy_draw_graph`).
fn scanpy_kmeans_next_rand(state : Int32) -> Int32 {
  (state * 1103515245 + 12345) & 0x7FFFFFFF
}

///|
/// Squared Euclidean distance between two vectors of equal length.
fn scanpy_kmeans_dist2(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  let d = if a.length() < b.length() { a.length() } else { b.length() }
  for t = 0; t < d; t = t + 1 {
    let diff = a[t] - b[t]
    s = s + diff * diff
  }
  s
}

///|
/// k-means++ seeding: pick the first center uniformly at random, then each
/// subsequent center with probability proportional to the squared distance to
/// the closest already-chosen center.
fn scanpy_kmeans_pp_seed(
  x : Array[Array[Double]],
  k : Int32,
  seed : Int32,
) -> Array[Array[Double]] {
  let n = x.length()
  let d = if n > 0 { x[0].length() } else { 0 }
  let centers : Array[Array[Double]] = Array(capacity=k)
  if n == 0 || k == 0 {
    return centers
  }
  let mut rng = if seed == 0 { 1 } else { seed }
  rng = scanpy_kmeans_next_rand(rng)
  let first = (rng % n).abs()
  let first_row : Array[Double] = Array(capacity=d)
  for t = 0; t < d; t = t + 1 {
    first_row.push(x[first][t])
  }
  centers.push(first_row)
  let min_d2 = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    min_d2[i] = scanpy_kmeans_dist2(x[i], first_row)
  }
  while centers.length() < k {
    let mut total = 0.0
    for i = 0; i < n; i = i + 1 {
      total = total + min_d2[i]
    }
    let mut chosen = n - 1
    if total > 0.0 {
      rng = scanpy_kmeans_next_rand(rng)
      let r = rng.to_double() / 2147483647.0 * total
      let mut acc = 0.0
      chosen = n - 1
      for i = 0; i < n; i = i + 1 {
        acc = acc + min_d2[i]
        if acc >= r {
          chosen = i
          break
        }
      }
    }
    let row : Array[Double] = Array(capacity=d)
    for t = 0; t < d; t = t + 1 {
      row.push(x[chosen][t])
    }
    centers.push(row)
    for i = 0; i < n; i = i + 1 {
      let d2 = scanpy_kmeans_dist2(x[i], row)
      if d2 < min_d2[i] {
        min_d2[i] = d2
      }
    }
  }
  centers
}

///|
/// One K-Means run (k-means++ seed + Lloyd iterations).  Empty clusters are
/// re-seeded with the currently farthest point.  Returns labels, centers,
/// inertia and the number of iterations performed.
fn scanpy_kmeans_single(
  x : Array[Array[Double]],
  k : Int32,
  max_iter : Int32,
  seed : Int32,
) -> (Array[Int32], Array[Array[Double]], Double, Int32) {
  let n = x.length()
  let d = if n > 0 { x[0].length() } else { 0 }
  if n == 0 || k == 0 {
    return (Array::make(n, 0), [], 0.0, 0)
  }
  let centers = scanpy_kmeans_pp_seed(x, k, seed)
  let labels = Array::make(n, 0)
  let mut inertia = 0.0
  let mut iters = 0
  for iter = 0; iter < max_iter; iter = iter + 1 {
    iters = iter + 1
    // assignment
    inertia = 0.0
    for i = 0; i < n; i = i + 1 {
      let mut best = 0
      let mut best_d2 = scanpy_kmeans_dist2(x[i], centers[0])
      for c = 1; c < k; c = c + 1 {
        let d2 = scanpy_kmeans_dist2(x[i], centers[c])
        if d2 < best_d2 {
          best_d2 = d2
          best = c
        }
      }
      labels[i] = best
      inertia = inertia + best_d2
    }
    // update
    let counts = Array::make(k, 0)
    let sums : Array[Array[Double]] = Array(capacity=k)
    for c = 0; c < k; c = c + 1 {
      let row = Array::make(d, 0.0)
      sums.push(row)
    }
    for i = 0; i < n; i = i + 1 {
      let c = labels[i]
      counts[c] = counts[c] + 1
      for t = 0; t < d; t = t + 1 {
        sums[c][t] = sums[c][t] + x[i][t]
      }
    }
    let mut max_shift = 0.0
    for c = 0; c < k; c = c + 1 {
      if counts[c] == 0 {
        // re-seed empty cluster with the farthest point from its assigned
        // center (falls back to the global farthest point)
        let mut far = 0
        let mut far_d2 = -1.0
        for i = 0; i < n; i = i + 1 {
          let d2 = scanpy_kmeans_dist2(x[i], centers[labels[i]])
          if d2 > far_d2 {
            far_d2 = d2
            far = i
          }
        }
        for t = 0; t < d; t = t + 1 {
          let old = centers[c][t]
          centers[c][t] = x[far][t]
          let shift = (centers[c][t] - old).abs()
          if shift > max_shift {
            max_shift = shift
          }
        }
      } else {
        for t = 0; t < d; t = t + 1 {
          let next = sums[c][t] / counts[c].to_double()
          let shift = (next - centers[c][t]).abs()
          if shift > max_shift {
            max_shift = shift
          }
          centers[c][t] = next
        }
      }
    }
    if max_shift < 1.0e-10 {
      break
    }
  }
  // recompute inertia after the final center update
  inertia = 0.0
  for i = 0; i < n; i = i + 1 {
    inertia = inertia + scanpy_kmeans_dist2(x[i], centers[labels[i]])
  }
  (labels, centers, inertia, iters)
}

///|
/// K-Means clustering of a `n_obs x n_vars` matrix.  Runs `n_init`
/// independent k-means++ seeded Lloyd optimizations and keeps the run with
/// the lowest inertia (sklearn semantics, scanpy `tl.kmeans` defaults).
pub fn scanpy_kmeans(
  x : Array[Array[Double]],
  n_clusters : Int32,
  max_iter~ : Int32 = 300,
  n_init~ : Int32 = 10,
  seed~ : Int32 = 0,
) -> ScanpyKMeans {
  let n = x.length()
  if n == 0 {
    return ScanpyKMeans::{
      labels: [],
      centers: [],
      inertia: 0.0,
      n_iter: 0,
      n_clusters: 0,
    }
  }
  let k = if n_clusters < 1 {
    1
  } else if n_clusters > n {
    n
  } else {
    n_clusters
  }
  let mut best_labels = Array::make(n, 0)
  let mut best_centers : Array[Array[Double]] = []
  let mut best_inertia = 1.0e300
  let mut best_iter = 0
  for init = 0; init < (if n_init < 1 { 1 } else { n_init }); init = init + 1 {
    let (labels, centers, inertia, iters) = scanpy_kmeans_single(
      x,
      k,
      max_iter,
      seed + init * 7919,
    )
    if inertia < best_inertia {
      best_inertia = inertia
      best_labels = labels
      best_centers = centers
      best_iter = iters
    }
  }
  ScanpyKMeans::{
    labels: best_labels,
    centers: best_centers,
    inertia: best_inertia,
    n_iter: best_iter,
    n_clusters: k,
  }
}

///|
/// K-Means on the AnnData expression matrix (dense materialization of sparse
/// storage on demand), storing the labels as a string column under
/// `adata.obs[obs_key]` and the clustering summary under
/// `adata.uns[uns_key]` (method, k, inertia).  Returns the matrix-level
/// result.
pub fn scanpy_kmeans_anndata(
  adata : AnnData,
  n_clusters : Int32,
  obs_key~ : String = "kmeans",
  uns_key~ : String = "kmeans",
  max_iter~ : Int32 = 300,
  n_init~ : Int32 = 10,
  seed~ : Int32 = 0,
) -> ScanpyKMeans {
  let x = adata.x_matrix()
  let res = scanpy_kmeans(x, n_clusters, max_iter~, n_init~, seed~)
  let strs : Array[String] = Array(capacity=res.labels.length())
  for l in res.labels {
    strs.push(l.to_string())
  }
  adata.obs[obs_key] = strs
  let info : Map[String, JsonValue] = Map([], capacity=4)
  info["method"] = JsonStr("kmeans")
  info["n_clusters"] = JsonNum(res.n_clusters.to_double())
  info["inertia"] = JsonNum(res.inertia)
  adata.uns[uns_key] = JsonObj(info)
  res
}
