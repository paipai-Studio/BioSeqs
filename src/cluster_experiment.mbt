///|
/// clusterExperiment - A unified framework for consensus and sequential
/// clustering of high-dimensional biological data.
///
/// Inspired by Bioconductor `clusterExperiment` (Risso et al. 2018,
/// Bioinformatics 34:1628). This pure-MoonBit implementation provides the
/// core RSEC (Resampling-based Sequential Ensembling Clustering) workflow:
///   - clusterMany: parameter-grid ensemble of clusterings
///   - makeConsensus: co-clustering matrix + consensus clustering
///   - makeDendrogram: hierarchical linkage of cluster medoids
///   - mergeClusters: merge by significance (proportion test + q-value)
///   - sequential clustering: iteratively isolate small clusters
///
/// References:
///   https://www.bioconductor.org/packages/release/bioc/html/clusterExperiment.html

// ============================================================================
// Core data structures
// ============================================================================

///|
/// A single clustering result: labels per sample plus the parameter key.
pub struct CeClustering {
  labels : Array[Int] // -1 = unassigned, otherwise 0-based cluster id
  k : Int
  reduced_dim : Int // 0 = no reduction, else PCA target dim
  subsample : Double // subsample fraction (1.0 = no subsampling)
}

///|
/// Construct a CeClustering.
pub fn CeClustering::new(
  labels~ : Array[Int],
  k~ : Int,
  reduced_dim~ : Int,
  subsample~ : Double,
) -> CeClustering {
  CeClustering::{ labels, k, reduced_dim, subsample }
}

///|
/// A node in a cluster dendrogram.
pub struct CeDendroNode {
  left : CeDendroNode?
  right : CeDendroNode?
  height : Double
  cluster_id : Int // leaf cluster id, -1 for internal
  members : Array[Int] // cluster ids under this node
}

///|
/// Construct a CeDendroNode.
pub fn CeDendroNode::new(
  left~ : CeDendroNode?,
  right~ : CeDendroNode?,
  height~ : Double,
  cluster_id~ : Int,
  members~ : Array[Int],
) -> CeDendroNode {
  CeDendroNode::{ left, right, height, cluster_id, members }
}

///|
/// The RSEC workflow result container.
pub struct ClusterExperiment {
  labels : Array[Int] // final per-sample cluster labels
  cocluster : Array[Array[Double]] // n×n co-clustering matrix [0,1]
  dendro : CeDendroNode? // cluster dendrogram (over current clusters)
  n_clusters : Int
  merge_cutoff : Double // q-value cutoff used for merging
  merge_info : Array[CeMergeInfo] // per-merge diagnostics
}

///|
/// Construct a ClusterExperiment.
pub fn ClusterExperiment::new(
  labels~ : Array[Int],
  cocluster~ : Array[Array[Double]],
  dendro~ : CeDendroNode?,
  n_clusters~ : Int,
  merge_cutoff~ : Double,
  merge_info~ : Array[CeMergeInfo],
) -> ClusterExperiment {
  ClusterExperiment::{
    labels,
    cocluster,
    dendro,
    n_clusters,
    merge_cutoff,
    merge_info,
  }
}

///|
/// Diagnostics for one dendrogram-merge step.
pub struct CeMergeInfo {
  node_id : Int
  left_cluster : Int
  right_cluster : Int
  proportion : Double // proportion of significant features
  q_value : Double // Storey q-value of the merge
  merged : Bool // whether the merge passed the cutoff
}

///|
/// Construct a CeMergeInfo.
pub fn CeMergeInfo::new(
  node_id~ : Int,
  left_cluster~ : Int,
  right_cluster~ : Int,
  proportion~ : Double,
  q_value~ : Double,
  merged~ : Bool,
) -> CeMergeInfo {
  CeMergeInfo::{
    node_id,
    left_cluster,
    right_cluster,
    proportion,
    q_value,
    merged,
  }
}

///|
/// Parameters for the RSEC workflow.
pub struct CeRsecParams {
  ks : Array[Int] // candidate k values for clusterMany
  reduced_dims : Array[Int] // candidate PCA dims (0 = no reduction)
  subsamples : Array[Double] // candidate subsample fractions
  consensus_proportion : Double // co-clustering threshold (default 0.7)
  merge_cutoff : Double // q-value cutoff for merging (default 0.05)
  sequential_alpha : Double // significance for sequential isolation (0.05)
  sequential_beta : Double // min proportion to keep a cluster (0.7)
  n_init : Int // k-means restarts
  max_iter : Int // k-means iterations
}

///|
/// Default RSEC parameters.
pub fn ce_default_params() -> CeRsecParams {
  CeRsecParams::{
    ks: [2, 3, 4],
    reduced_dims: [0, 5],
    subsamples: [1.0, 0.8],
    consensus_proportion: 0.7,
    merge_cutoff: 0.05,
    sequential_alpha: 0.05,
    sequential_beta: 0.7,
    n_init: 5,
    max_iter: 100,
  }
}

// ============================================================================
// Numeric helpers
// ============================================================================

///|
/// Euclidean distance between two equal-length points.
pub fn ce_euclidean(a : Array[Double], b : Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<a.length() {
    let d = a[i] - b[i]
    sum = sum + d * d
  }
  sum.sqrt()
}

///|
/// Mean of an array of doubles. Returns 0.0 if empty.
pub fn ce_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in values {
    s = s + v
  }
  s / values.length().to_double()
}

///|
/// Variance (population, 1/n) of an array of doubles.
pub fn ce_variance(values : Array[Double]) -> Double {
  let n = values.length()
  if n < 2 {
    return 0.0
  }
  let m = ce_mean(values)
  let mut s = 0.0
  for v in values {
    let d = v - m
    s = s + d * d
  }
  s / n.to_double()
}

///|
/// Standard deviation (population) of an array of doubles.
pub fn ce_sd(values : Array[Double]) -> Double {
  ce_variance(values).sqrt()
}

///|
/// Median of an array of doubles (does not mutate input).
pub fn ce_median(values : Array[Double]) -> Double {
  let n = values.length()
  if n == 0 {
    return 0.0
  }
  let sorted = Array::new()
  for v in values {
    sorted.push(v)
  }
  // simple insertion sort (n usually small for per-feature stats)
  for i in 1..<n {
    let key = sorted[i]
    let mut j = i - 1
    while j >= 0 && sorted[j] > key {
      sorted[j + 1] = sorted[j]
      j = j - 1
    }
    sorted[j + 1] = key
  }
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

///|
/// Create an n×m matrix filled with a constant.
pub fn ce_make_matrix(n : Int, m : Int, fill : Double) -> Array[Array[Double]] {
  let out = Array::new()
  for _ in 0..<n {
    let row = Array::make(m, fill)
    out.push(row)
  }
  out
}

///|
/// Transpose a (rows × cols) matrix.
pub fn ce_transpose(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = m.length()
  if n == 0 {
    return []
  }
  let p = m[0].length()
  let out = Array::new()
  for j in 0..<p {
    let row = Array::new()
    for i in 0..<n {
      row.push(m[i][j])
    }
    out.push(row)
  }
  out
}

// ============================================================================
// PCA via power iteration on the covariance matrix (simple, deterministic)
// ============================================================================

///|
/// Center each column of a matrix in place (returns a new matrix).
pub fn ce_center_columns(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = m.length()
  if n == 0 {
    return []
  }
  let p = m[0].length()
  let out = Array::new()
  for row in m {
    out.push(Array::make(p, 0.0))
  }
  for j in 0..<p {
    let col = Array::new()
    for i in 0..<n {
      col.push(m[i][j])
    }
    let mu = ce_mean(col)
    for i in 0..<n {
      out[i][j] = m[i][j] - mu
    }
  }
  out
}

///|
/// Compute the leading `n_components` principal components of a data matrix
/// (samples × features) and return the projected samples.
///
/// Uses deflated power iteration on the sample covariance matrix. Suitable
/// for moderate dimensions; deterministic given the input.
pub fn ce_pca(
  data : Array[Array[Double]],
  n_components : Int,
) -> Array[Array[Double]] {
  let n = data.length()
  if n == 0 || n_components <= 0 {
    return []
  }
  let centered = ce_center_columns(data)
  let p = centered[0].length()
  let k = if n_components < p { n_components } else { p }
  // Project by accumulating leading eigenvectors of (1/n) X^T X
  // For each component: power-iterate on residual covariance.
  let components : Array[Array[Double]] = Array::new()
  let residual = centered
  for _ in 0..<k {
    let comp = Array::make(p, 0.0)
    for j in 0..<p {
      comp[j] = 1.0 / (j + 1).to_double() // deterministic-ish init
    }
    // normalize init
    let mut norm = 0.0
    for v in comp {
      norm = norm + v * v
    }
    norm = norm.sqrt()
    if norm > 0.0 {
      for j in 0..<p {
        comp[j] = comp[j] / norm
      }
    }
    // power iteration (30 steps)
    for _ in 0..<30 {
      // multiply X^T (X comp)
      let xc = Array::make(n, 0.0)
      for i in 0..<n {
        let mut s = 0.0
        for j in 0..<p {
          s = s + residual[i][j] * comp[j]
        }
        xc[i] = s
      }
      let next = Array::make(p, 0.0)
      for j in 0..<p {
        let mut s = 0.0
        for i in 0..<n {
          s = s + residual[i][j] * xc[i]
        }
        next[j] = s
      }
      let mut nn = 0.0
      for v in next {
        nn = nn + v * v
      }
      nn = nn.sqrt()
      if nn > 0.0 {
        for j in 0..<p {
          next[j] = next[j] / nn
        }
      }
      for j in 0..<p {
        comp[j] = next[j]
      }
    }
    // deflate: residual -= (X comp)(X comp)^T contribution via column space
    // i.e. residual = residual - (X comp) comp^T
    let proj = Array::make(n, 0.0)
    for i in 0..<n {
      let mut s = 0.0
      for j in 0..<p {
        s = s + residual[i][j] * comp[j]
      }
      proj[i] = s
    }
    for i in 0..<n {
      for j in 0..<p {
        residual[i][j] = residual[i][j] - proj[i] * comp[j]
      }
    }
    components.push(comp)
  }
  // Project original centered data onto components
  let out = Array::new()
  for i in 0..<n {
    let row = Array::make(k, 0.0)
    for c in 0..<k {
      let mut s = 0.0
      for j in 0..<p {
        s = s + centered[i][j] * components[c][j]
      }
      row[c] = s
    }
    out.push(row)
  }
  out
}

// ============================================================================
// K-means (Lloyd's algorithm) with k++ init
// ============================================================================

///|
/// Run k-means on `data` returning integer labels (0..k-1) and centroids.
pub fn ce_kmeans(
  data : Array[Array[Double]],
  k : Int,
  n_init : Int,
  max_iter : Int,
) -> (Array[Int], Array[Array[Double]]) {
  let n = data.length()
  if n == 0 || k <= 0 {
    return ([], [])
  }
  if k == 1 {
    return (Array::make(n, 0), [ce_centroid(data)])
  }
  let mut best_labels = Array::make(n, 0)
  let mut best_centroids : Array[Array[Double]] = Array::new()
  let mut best_inertia = 1.0e300
  for _ in 0..<n_init {
    let (labels, centroids, inertia) = ce_kmeans_single(data, k, max_iter)
    if inertia < best_inertia {
      best_inertia = inertia
      best_labels = labels
      best_centroids = centroids
    }
  }
  (best_labels, best_centroids)
}

///|
/// Compute the centroid (mean point) of a set of points.
pub fn ce_centroid(points : Array[Array[Double]]) -> Array[Double] {
  let n = points.length()
  if n == 0 {
    return []
  }
  let p = points[0].length()
  let out = Array::make(p, 0.0)
  for pt in points {
    for j in 0..<p {
      out[j] = out[j] + pt[j]
    }
  }
  for j in 0..<p {
    out[j] = out[j] / n.to_double()
  }
  out
}

///|
/// Single k-means run with k-means++ initialization.
fn ce_kmeans_single(
  data : Array[Array[Double]],
  k : Int,
  max_iter : Int,
) -> (Array[Int], Array[Array[Double]], Double) {
  let n = data.length()
  let centroids = ce_kmeans_pp_init(data, k)
  let labels = Array::make(n, 0)
  let mut inertia = 0.0
  for _ in 0..<max_iter {
    // assign
    let mut new_inertia = 0.0
    for i in 0..<n {
      let mut best = 0
      let mut best_d = 1.0e300
      for c in 0..<k {
        let d = ce_euclidean(data[i], centroids[c])
        if d < best_d {
          best_d = d
          best = c
        }
      }
      labels[i] = best
      new_inertia = new_inertia + best_d * best_d
    }
    // update
    let sums : Array[Array[Double]] = Array::new()
    let counts = Array::make(k, 0)
    for c in 0..<k {
      sums.push(Array::make(data[0].length(), 0.0))
    }
    for i in 0..<n {
      let c = labels[i]
      counts[c] = counts[c] + 1
      for j in 0..<data[0].length() {
        sums[c][j] = sums[c][j] + data[i][j]
      }
    }
    for c in 0..<k {
      if counts[c] > 0 {
        for j in 0..<data[0].length() {
          sums[c][j] = sums[c][j] / counts[c].to_double()
        }
        centroids[c] = sums[c]
      }
    }
    if (new_inertia - inertia).abs() < 1.0e-12 {
      inertia = new_inertia
      break
    }
    inertia = new_inertia
  }
  (labels, centroids, inertia)
}

///|
/// k-means++ seeding using a deterministic pseudo-random pick based on indices
/// (the project avoids real RNG; this gives reproducible spreads).
fn ce_kmeans_pp_init(data : Array[Array[Double]], k : Int) -> Array[Array[Double]] {
  let n = data.length()
  let centroids : Array[Array[Double]] = Array::new()
  // pick the point closest to the global centroid as first seed
  let global_center = ce_centroid(data)
  let mut best = 0
  let mut best_d = 1.0e300
  for i in 0..<n {
    let d = ce_euclidean(data[i], global_center)
    if d < best_d {
      best_d = d
      best = i
    }
  }
  centroids.push(data[best])
  while centroids.length() < k {
    // distances to nearest existing centroid
    let dists = Array::make(n, 0.0)
    let mut total = 0.0
    for i in 0..<n {
      let mut dmin = 1.0e300
      for c in centroids {
        let d = ce_euclidean(data[i], c)
        if d < dmin {
          dmin = d
        }
      }
      dists[i] = dmin * dmin
      total = total + dists[i]
    }
    if total <= 0.0 {
      // all points coincide with existing centroids; just pick by index
      centroids.push(data[centroids.length() % n])
      continue
    }
    // deterministic pick: choose the point with largest cumulative share
    // (a proxy for the "furthest" point, reproducible)
    let target = 0.5
    let mut cum = 0.0
    let mut picked = n - 1
    for i in 0..<n {
      cum = cum + dists[i] / total
      if cum >= target {
        picked = i
        break
      }
    }
    centroids.push(data[picked])
  }
  centroids
}

// ============================================================================
// clusterMany + makeConsensus
// ============================================================================

///|
/// Run clusterMany: an ensemble of clusterings over a grid of parameters.
/// Returns one CeClustering per parameter combination.
pub fn ce_cluster_many(
  data : Array[Array[Double]],
  params : CeRsecParams,
) -> Array[CeClustering] {
  let results : Array[CeClustering] = Array::new()
  for k in params.ks {
    for rd in params.reduced_dims {
      for ss in params.subsamples {
        let work = if rd > 0 {
          ce_pca(data, rd)
        } else {
          let copy : Array[Array[Double]] = Array::new()
          for row in data {
            copy.push(row)
          }
          copy
        }
        let subset = ce_subsample(work, ss)
        let (labels_full, _) = ce_kmeans(subset.data, k, params.n_init, params.max_iter)
        // map subsample labels back to full data
        let labels = Array::make(data.length(), -1)
        for idx in 0..<subset.indices.length() {
          labels[subset.indices[idx]] = labels_full[idx]
        }
        results.push(
          CeClustering::new(labels=labels, k=k, reduced_dim=rd, subsample=ss),
        )
      }
    }
  }
  results
}

///|
/// Subsample result.
struct CeSubsample {
  data : Array[Array[Double]]
  indices : Array[Int]
}

///|
/// Deterministic subsampling: take an evenly-spaced fraction of the data.
fn ce_subsample(data : Array[Array[Double]], fraction : Double) -> CeSubsample {
  let n = data.length()
  if fraction >= 1.0 {
    let copy : Array[Array[Double]] = Array::new()
    let idx : Array[Int] = Array::new()
    for i in 0..<n {
      copy.push(data[i])
      idx.push(i)
    }
    return CeSubsample::{ data: copy, indices: idx }
  }
  let keep_count = if n.to_double() * fraction < 1.0 {
    1
  } else {
    (n.to_double() * fraction).to_int()
  }
  let step = n.to_double() / keep_count.to_double()
  let out : Array[Array[Double]] = Array::new()
  let idx : Array[Int] = Array::new()
  let mut pos = 0.0
  for _ in 0..<keep_count {
    let i = pos.to_int()
    let ii = if i >= n { n - 1 } else { i }
    out.push(data[ii])
    idx.push(ii)
    pos = pos + step
  }
  CeSubsample::{ data: out, indices: idx }
}

///|
/// Build the n×n co-clustering matrix from an ensemble of clusterings.
/// Entry (i,j) = fraction of clusterings where i and j share a label
/// (counting only runs where both were assigned).
pub fn ce_cocluster_matrix(
  ensemble : Array[CeClustering],
  n_samples : Int,
) -> Array[Array[Double]] {
  let m = ce_make_matrix(n_samples, n_samples, 0.0)
  let counts = Array::make(n_samples, 0)
  for cl in ensemble {
    // mark which samples were assigned in this run
    for i in 0..<n_samples {
      if cl.labels[i] >= 0 {
        counts[i] = counts[i] + 1
      }
    }
    for i in 0..<n_samples {
      if cl.labels[i] < 0 {
        continue
      }
      for j in (i + 1)..<n_samples {
        if cl.labels[j] < 0 {
          continue
        }
        if cl.labels[i] == cl.labels[j] {
          m[i][j] = m[i][j] + 1.0
          m[j][i] = m[j][i] + 1.0
        }
      }
    }
  }
  for i in 0..<n_samples {
    m[i][i] = counts[i].to_double()
    for j in 0..<n_samples {
      if i != j && counts[i] > 0 && counts[j] > 0 {
        let denom = (if counts[i] < counts[j] {
          counts[i]
        } else {
          counts[j]
        }).to_double()
        if denom > 0.0 {
          m[i][j] = m[i][j] / denom
        }
      }
    }
  }
  m
}

///|
/// Derive a consensus clustering from the co-clustering matrix using a
/// distance of (1 - co-clustering) and average-linkage hierarchical clustering,
/// then cut at the requested number of clusters.
pub fn ce_make_consensus(
  cocluster : Array[Array[Double]],
  k : Int,
) -> Array[Int] {
  let n = cocluster.length()
  if n == 0 {
    return []
  }
  // distance matrix
  let dist = ce_make_matrix(n, n, 0.0)
  for i in 0..<n {
    for j in 0..<n {
      dist[i][j] = 1.0 - cocluster[i][j]
    }
  }
  ce_hclust_cut(dist, k)
}

///|
/// Agglomerative hierarchical clustering (average linkage) on a distance
/// matrix, cut to produce exactly `k` clusters. Returns 0-based labels.
pub fn ce_hclust_cut(dist : Array[Array[Double]], k : Int) -> Array[Int] {
  let n = dist.length()
  if n == 0 {
    return []
  }
  if k >= n {
    let out = Array::new()
    for i in 0..<n {
      out.push(i)
    }
    return out
  }
  // each sample starts in its own cluster
  let clusters : Array[Array[Int]] = Array::new()
  for i in 0..<n {
    clusters.push([i])
  }
  // active cluster indices
  let active : Array[Int] = Array::new()
  for i in 0..<n {
    active.push(i)
  }
  while active.length() > k {
    // find the closest pair of active clusters (average linkage)
    let mut best_i = 0
    let mut best_j = 1
    let mut best_d = 1.0e300
    for a in 0..<active.length() {
      for b in (a + 1)..<active.length() {
        let ci = active[a]
        let cj = active[b]
        // average linkage
        let mut s = 0.0
        for x in clusters[ci] {
          for y in clusters[cj] {
            s = s + dist[x][y]
          }
        }
        let d = s / (clusters[ci].length() * clusters[cj].length()).to_double()
        if d < best_d {
          best_d = d
          best_i = a
          best_j = b
        }
      }
    }
    // merge active[best_j] into active[best_i]
    let ci = active[best_i]
    let cj = active[best_j]
    for x in clusters[cj] {
      clusters[ci].push(x)
    }
    clusters[cj] = [] // mark empty
    // remove active[best_j] by shifting subsequent elements left
    let _ = best_i
    let alen = active.length()
    let mut p = best_j
    while p < alen - 1 {
      active[p] = active[p + 1]
      p = p + 1
    }
    ignore(active.pop())
  }
  // assign labels
  let labels = Array::make(n, -1)
  let mut next_id = 0
  for ci in active {
    if clusters[ci].length() == 0 {
      continue
    }
    for x in clusters[ci] {
      labels[x] = next_id
    }
    next_id = next_id + 1
  }
  labels
}

// ============================================================================
// makeDendrogram: hierarchical clustering of cluster medoids
// ============================================================================

///|
/// Build a dendrogram over clusters using medoids and average linkage.
/// `labels` are 0-based cluster ids per sample; returns the root node.
pub fn ce_make_dendrogram(
  data : Array[Array[Double]],
  labels : Array[Int],
) -> CeDendroNode? {
  // collect cluster ids
  let ids : Array[Int] = Array::new()
  for l in labels {
    if l >= 0 {
      let mut found = false
      for id in ids {
        if id == l {
          found = true
          break
        }
      }
      if !found {
        ids.push(l)
      }
    }
  }
  if ids.length() < 2 {
    return None
  }
  // compute medoids
  let medoids : Array[(Int, Array[Double])] = Array::new()
  for id in ids {
    let pts : Array[Array[Double]] = Array::new()
    for i in 0..<data.length() {
      if labels[i] == id {
        pts.push(data[i])
      }
    }
    medoids.push((id, ce_centroid(pts)))
  }
  // distance matrix between medoids
  let m = medoids.length()
  let dist = ce_make_matrix(m, m, 0.0)
  for i in 0..<m {
    for j in 0..<m {
      if i != j {
        dist[i][j] = ce_euclidean(medoids[i].1, medoids[j].1)
      }
    }
  }
  // agglomerative merge, recording nodes
  let nodes : Array[CeDendroNode] = Array::new()
  for i in 0..<m {
    nodes.push(
      CeDendroNode::new(left=None, right=None, height=0.0, cluster_id=i, members=[i]),
    )
  }
  let active : Array[Int] = Array::new()
  for i in 0..<m {
    active.push(i)
  }
  while active.length() > 1 {
    let mut best_a = 0
    let mut best_b = 1
    let mut best_d = 1.0e300
    for a in 0..<active.length() {
      for b in (a + 1)..<active.length() {
        let na = active[a]
        let nb = active[b]
        let mut s = 0.0
        for x in nodes[na].members {
          for y in nodes[nb].members {
            s = s + dist[x][y]
          }
        }
        let d = s / (nodes[na].members.length() * nodes[nb].members.length()).to_double()
        if d < best_d {
          best_d = d
          best_a = a
          best_b = b
        }
      }
    }
    let na = active[best_a]
    let nb = active[best_b]
    let merged_members = Array::new()
    for x in nodes[na].members {
      merged_members.push(x)
    }
    for x in nodes[nb].members {
      merged_members.push(x)
    }
    let parent = CeDendroNode::new(
      left=Some(nodes[na]),
      right=Some(nodes[nb]),
      height=best_d,
      cluster_id=-1,
      members=merged_members,
    )
    nodes.push(parent)
    let new_idx = nodes.length() - 1
    // remove active[best_b] by shifting left
    let alen = active.length()
    let mut p = best_b
    while p < alen - 1 {
      active[p] = active[p + 1]
      p = p + 1
    }
    ignore(active.pop())
    active[best_a] = new_idx
  }
  if active.length() == 1 {
    Some(nodes[active[0]])
  } else {
    None
  }
}

// ============================================================================
// mergeClusters: merge dendrogram nodes by significance
// ============================================================================

///|
/// Merge clusters by testing per-feature proportion differences.
///
/// For each internal dendrogram node, compare the two child clusters
/// feature-by-feature using a two-sample proportion-style z-test on the
/// fraction of samples above the cluster's per-feature median. Cluster pairs
/// are merged when the Storey q-value of the node exceeds `merge_cutoff`.
///
/// Returns updated labels and per-node merge diagnostics.
pub fn ce_merge_clusters(
  data : Array[Array[Double]],
  labels : Array[Int],
  dendro : CeDendroNode?,
  merge_cutoff : Double,
) -> (Array[Int], Array[CeMergeInfo]) {
  let info : Array[CeMergeInfo] = Array::new()
  if dendro is None {
    return (labels, info)
  }
  let root = dendro.unwrap()
  // collect leaf order via in-order traversal
  let leaves : Array[Int] = Array::new()
  ce_collect_leaves(root, leaves)
  // map cluster index -> original cluster id
  // (leaves are in dendrogram order)
  // Build merge map: initially each leaf is its own cluster
  let parent : Map[Int, Int] = Map([], capacity=16)
  for i in 0..<leaves.length() {
    parent[i] = leaves[i]
  }
  // walk internal nodes in order of increasing height
  let internals : Array[CeDendroNode] = Array::new()
  ce_collect_internals(root, internals)
  // sort by height ascending (insertion sort, small n)
  for i in 1..<internals.length() {
    let key = internals[i]
    let mut j = i - 1
    while j >= 0 && internals[j].height > key.height {
      internals[j + 1] = internals[j]
      j = j - 1
    }
    internals[j + 1] = key
  }
  let mut node_id = 0
  for node in internals {
    let left_leaves : Array[Int] = Array::new()
    let right_leaves : Array[Int] = Array::new()
    match node.left {
      Some(l) => ce_collect_leaves(l, left_leaves)
      None => ()
    }
    match node.right {
      Some(r) => ce_collect_leaves(r, right_leaves)
      None => ()
    }
    // find current cluster ids for these leaves
    let left_cluster = ce_find_root(parent, left_leaves[0])
    let right_cluster = ce_find_root(parent, right_leaves[0])
    // proportion test across features
    let (prop, qval) = ce_proportion_test(data, labels, left_cluster, right_cluster)
    let merged = qval > merge_cutoff
    if merged {
      // union the two clusters
      parent[right_cluster] = left_cluster
    }
    info.push(
      CeMergeInfo::new(
        node_id=node_id,
        left_cluster=left_cluster,
        right_cluster=right_cluster,
        proportion=prop,
        q_value=qval,
        merged=merged,
      ),
    )
    node_id = node_id + 1
  }
  // relabel: find roots and assign new contiguous ids
  let new_labels = Array::make(labels.length(), -1)
  let root_map : Map[Int, Int] = Map([], capacity=16)
  let mut next_id = 0
  for i in 0..<labels.length() {
    let l = labels[i]
    if l < 0 {
      continue
    }
    let r = ce_find_root(parent, l)
    match root_map.get(r) {
      Some(id) => new_labels[i] = id
      None => {
        root_map[r] = next_id
        new_labels[i] = next_id
        next_id = next_id + 1
      }
    }
  }
  (new_labels, info)
}

///|
/// Union-find: find the root of a cluster id.
fn ce_find_root(parent : Map[Int, Int], x : Int) -> Int {
  let mut cur = x
  while true {
    match parent.get(cur) {
      Some(p) => {
        if p == cur {
          return cur
        }
        cur = p
      }
      None => return cur
    }
  }
  0
}

///|
/// Collect leaf cluster indices from a dendrogram subtree.
fn ce_collect_leaves(node : CeDendroNode, out : Array[Int]) -> Unit {
  if node.cluster_id >= 0 {
    out.push(node.cluster_id)
    return
  }
  match node.left {
    Some(l) => ce_collect_leaves(l, out)
    None => ()
  }
  match node.right {
    Some(r) => ce_collect_leaves(r, out)
    None => ()
  }
}

///|
/// Collect internal nodes from a dendrogram subtree.
fn ce_collect_internals(node : CeDendroNode, out : Array[CeDendroNode]) -> Unit {
  if node.cluster_id >= 0 {
    return
  }
  out.push(node)
  match node.left {
    Some(l) => ce_collect_internals(l, out)
    None => ()
  }
  match node.right {
    Some(r) => ce_collect_internals(r, out)
    None => ()
  }
}

///|
/// Two-sample test per feature: fraction of samples above the pooled median.
/// Returns (max proportion across features, approximate q-value).
/// Uses a simple z-test with Storey-style q-value (pi0 = 1).
fn ce_proportion_test(
  data : Array[Array[Double]],
  labels : Array[Int],
  cluster_a : Int,
  cluster_b : Int,
) -> (Double, Double) {
  let n = data.length()
  let p = if n > 0 { data[0].length() } else { 0 }
  if p == 0 {
    return (0.0, 1.0)
  }
  // gather samples per cluster
  let a_idx : Array[Int] = Array::new()
  let b_idx : Array[Int] = Array::new()
  for i in 0..<n {
    if labels[i] == cluster_a {
      a_idx.push(i)
    } else if labels[i] == cluster_b {
      b_idx.push(i)
    }
  }
  if a_idx.length() == 0 || b_idx.length() == 0 {
    return (0.0, 1.0)
  }
  // per-feature two-sample t-test (unequal variance)
  let mut max_t = 0.0
  let pvals : Array[Double] = Array::new()
  for j in 0..<p {
    let a_vals : Array[Double] = Array::new()
    let b_vals : Array[Double] = Array::new()
    for i in a_idx {
      a_vals.push(data[i][j])
    }
    for i in b_idx {
      b_vals.push(data[i][j])
    }
    let ma = ce_mean(a_vals)
    let mb = ce_mean(b_vals)
    let va = ce_variance(a_vals)
    let vb = ce_variance(b_vals)
    let na = a_vals.length().to_double()
    let nb = b_vals.length().to_double()
    let se = (va / na + vb / nb).sqrt()
    if se <= 0.0 {
      pvals.push(1.0)
      continue
    }
    let t = (ma - mb).abs() / se
    if t > max_t {
      max_t = t
    }
    // two-sided p-value via normal approximation (Welch)
    let z = t
    let pval = 2.0 * (1.0 - ce_normal_cdf(z))
    pvals.push(pval)
  }
  // Storey q-value with pi0 = 1 (BH-style)
  let qvals = ce_bh_fdr(pvals)
  // proportion of features with q <= 0.05
  let mut sig = 0
  for q in qvals {
    if q <= 0.05 {
      sig = sig + 1
    }
  }
  let prop = sig.to_double() / p.to_double()
  // node q-value: minimum q across features (representative)
  let mut min_q = 1.0
  for q in qvals {
    if q < min_q {
      min_q = q
    }
  }
  (prop, min_q)
}

///|
/// Standard normal CDF (Abramowitz & Stegun 26.2.17 approximation).
pub fn ce_normal_cdf(z : Double) -> Double {
  let x = z
  if x < 0.0 {
    return 1.0 - ce_normal_cdf(-x)
  }
  let t = 1.0 / (1.0 + 0.2316419 * x)
  let d = 0.3989423 * @math.exp(-x * x / 2.0)
  let p = t *
    (0.319381530 +
      t *
        (-0.356563782 +
          t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
  1.0 - d * p
}

///|
/// Benjamini-Hochberg FDR correction. Returns q-values.
pub fn ce_bh_fdr(pvals : Array[Double]) -> Array[Double] {
  let n = pvals.length()
  if n == 0 {
    return []
  }
  // sort indices by p-value ascending
  let order : Array[Int] = Array::new()
  for i in 0..<n {
    order.push(i)
  }
  // simple selection sort on indices
  for i in 0..<n {
    let mut min_j = i
    for j in (i + 1)..<n {
      if pvals[order[j]] < pvals[order[min_j]] {
        min_j = j
      }
    }
    let tmp = order[i]
    order[i] = order[min_j]
    order[min_j] = tmp
  }
  let qvals = Array::make(n, 0.0)
  let mut prev = 1.0
  // iterate from largest to smallest p-value (rank = n-1 down to 0)
  let mut rank = n - 1
  while rank >= 0 {
    let idx = order[rank]
    let p = pvals[idx]
    let q = p * n.to_double() / (rank + 1).to_double()
    if q > prev {
      qvals[idx] = prev
    } else {
      qvals[idx] = q
      prev = q
    }
    if qvals[idx] > 1.0 {
      qvals[idx] = 1.0
    }
    rank = rank - 1
  }
  qvals
}

// ============================================================================
// Sequential clustering: iteratively isolate small significant clusters
// ============================================================================

///|
/// Sequential clustering: repeatedly run k-means with k=2 on the current
/// cluster; if the smaller sub-cluster shows feature-wise significance
/// (Welch t-test q-value <= alpha), keep it as a new cluster and continue
/// on the larger sub-cluster. Otherwise stop.
///
/// `beta` is the minimum proportion of the parent that the larger child must
/// retain to continue splitting.
pub fn ce_sequential_cluster(
  data : Array[Array[Double]],
  labels : Array[Int],
  alpha : Double,
  beta : Double,
) -> Array[Int] {
  let n = data.length()
  let out = Array::make(n, -1)
  // collect current cluster ids
  let ids : Array[Int] = Array::new()
  for l in labels {
    let mut found = false
      for id in ids {
        if id == l {
          found = true
          break
        }
      }
      if !found && l >= 0 {
      ids.push(l)
    }
  }
  let mut next_id = 0
  for id in ids {
    // gather samples in this cluster
    let idx : Array[Int] = Array::new()
    for i in 0..<n {
      if labels[i] == id {
        idx.push(i)
      }
    }
    if idx.length() < 2 {
      for i in idx {
        out[i] = next_id
      }
      next_id = next_id + 1
      continue
    }
    // run sequential isolation on this subset
    let subsets = ce_sequential_split(data, idx, alpha, beta)
    for s in subsets {
      for i in s {
        out[i] = next_id
      }
      next_id = next_id + 1
    }
  }
  out
}

///|
/// Recursively split a set of sample indices, isolating small significant
/// clusters.
fn ce_sequential_split(
  data : Array[Array[Double]],
  idx : Array[Int],
  alpha : Double,
  beta : Double,
) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = Array::new()
  let mut current = idx
  while current.length() >= 4 {
    // k-means with k=2
    let pts : Array[Array[Double]] = Array::new()
    for i in current {
      pts.push(data[i])
    }
    let (sub_labels, _) = ce_kmeans(pts, 2, 3, 50)
    let group_a : Array[Int] = Array::new()
    let group_b : Array[Int] = Array::new()
    for i in 0..<current.length() {
      if sub_labels[i] == 0 {
        group_a.push(current[i])
      } else {
        group_b.push(current[i])
      }
    }
    // determine smaller vs larger
    let (small, large) = if group_a.length() <= group_b.length() {
      (group_a, group_b)
    } else {
      (group_b, group_a)
    }
    if small.length() == 0 {
      break
    }
    // significance test: build a subset data matrix with labels
    let subset_data : Array[Array[Double]] = Array::new()
    let subset_labels : Array[Int] = Array::new()
    for idx in small {
      subset_data.push(data[idx])
      subset_labels.push(0)
    }
    for idx in large {
      subset_data.push(data[idx])
      subset_labels.push(1)
    }
    let (_, qval) = ce_proportion_test(subset_data, subset_labels, 0, 1)
    let large_frac = large.length().to_double() / current.length().to_double()
    if qval <= alpha && large_frac >= beta {
      // accept the small cluster
      out.push(small)
      current = large
    } else {
      break
    }
  }
  out.push(current)
  out
}

///|
/// Build a labels array (0 for small, 1 for large) for the proportion test.
fn ce_labels_for(
  all : Array[Int],
  small : Array[Int],
  large : Array[Int],
) -> Array[Int] {
  let n = all.length()
  let out = Array::make(n, -1)
  for i in 0..<n {
    for s in small {
      if all[i] == s {
        out[i] = 0
      }
    }
    for l in large {
      if all[i] == l {
        out[i] = 1
      }
    }
  }
  out
}

// ============================================================================
// RSEC: the full pipeline
// ============================================================================

///|
/// Run the full RSEC workflow on `data`.
///
/// Steps:
///   1. clusterMany over the parameter grid
///   2. makeConsensus to get a single clustering
///   3. makeDendrogram over the consensus clusters
///   4. mergeClusters by significance
///   5. sequential clustering to isolate small clusters
pub fn ce_rsec(
  data : Array[Array[Double]],
  params : CeRsecParams,
) -> ClusterExperiment {
  let n = data.length()
  // 1. clusterMany
  let ensemble = ce_cluster_many(data, params)
  // 2. consensus
  let cocluster = ce_cocluster_matrix(ensemble, n)
  // pick the median k from the grid as the consensus k
  let sorted_ks = Array::new()
  for k in params.ks {
    sorted_ks.push(k)
  }
  for i in 1..<sorted_ks.length() {
    let key = sorted_ks[i]
    let mut j = i - 1
    while j >= 0 && sorted_ks[j] > key {
      sorted_ks[j + 1] = sorted_ks[j]
      j = j - 1
    }
    sorted_ks[j + 1] = key
  }
  let consensus_k = sorted_ks[sorted_ks.length() / 2]
  let consensus_labels = ce_make_consensus(cocluster, consensus_k)
  // 3. dendrogram
  let dendro = ce_make_dendrogram(data, consensus_labels)
  // 4. mergeClusters
  let (merged_labels, merge_info) = ce_merge_clusters(
    data,
    consensus_labels,
    dendro,
    params.merge_cutoff,
  )
  // 5. sequential
  let final_labels = ce_sequential_cluster(
    data,
    merged_labels,
    params.sequential_alpha,
    params.sequential_beta,
  )
  // count clusters
  let ids : Array[Int] = Array::new()
  for l in final_labels {
    if l >= 0 {
      let mut found = false
      for id in ids {
        if id == l {
          found = true
          break
        }
      }
      if !found {
        ids.push(l)
      }
    }
  }
  ClusterExperiment::new(
    labels=final_labels,
    cocluster=cocluster,
    dendro=dendro,
    n_clusters=ids.length(),
    merge_cutoff=params.merge_cutoff,
    merge_info=merge_info,
  )
}

// ============================================================================
// Utilities: accessors and summary
// ============================================================================

///|
/// Count the number of distinct non-negative labels.
pub fn ce_count_clusters(labels : Array[Int]) -> Int {
  let ids : Array[Int] = Array::new()
  for l in labels {
    if l >= 0 {
      let mut found = false
      for id in ids {
        if id == l {
          found = true
          break
        }
      }
      if !found {
        ids.push(l)
      }
    }
  }
  ids.length()
}

///|
/// Compute cluster sizes (map cluster id -> size).
pub fn ce_cluster_sizes(labels : Array[Int]) -> Array[(Int, Int)] {
  let ids : Array[Int] = Array::new()
  let sizes : Array[Int] = Array::new()
  for l in labels {
    if l < 0 {
      continue
    }
    let mut idx = -1
    for i in 0..<ids.length() {
      if ids[i] == l {
        idx = i
        break
      }
    }
    if idx < 0 {
      ids.push(l)
      sizes.push(1)
    } else {
      sizes[idx] = sizes[idx] + 1
    }
  }
  let out : Array[(Int, Int)] = Array::new()
  for i in 0..<ids.length() {
    out.push((ids[i], sizes[i]))
  }
  out
}

///|
/// Format a ClusterExperiment result as a human-readable summary string.
pub fn ce_summary(result : ClusterExperiment) -> String {
  let sizes = ce_cluster_sizes(result.labels)
  let mut s = "=== clusterExperiment RSEC Summary ===\n"
  s = s + "Samples: " + result.labels.length().to_string() + "\n"
  s = s + "Clusters: " + result.n_clusters.to_string() + "\n"
  s = s + "Merge cutoff (q-value): " + result.merge_cutoff.to_string() + "\n"
  s = s + "Merge steps: " + result.merge_info.length().to_string() + "\n"
  let mut merged = 0
  for info in result.merge_info {
    if info.merged {
      merged = merged + 1
    }
  }
  s = s + "Merges applied: " + merged.to_string() + "\n"
  s = s + "Cluster sizes:\n"
  let mut si = 0
  while si < sizes.length() {
    let (id, size) = sizes[si]
    s = s +
      "  cluster " +
      id.to_string() +
      ": " +
      size.to_string() +
      " samples\n"
    si = si + 1
  }
  s
}
