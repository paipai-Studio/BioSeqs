///|
/// rrvgo-style semantic reduction of redundant GO terms.
///
/// Inspired by the R/Bioconductor package rrvgo (Sayols 2020), this module
/// clusters a set of enriched GO terms by their semantic similarity (using
/// the Wang measure from gosemsim) and picks one representative term per
/// cluster, reducing redundancy among GO enrichment results.
///
/// Reference: rrvgo: Reduce and Visualize GO terms.

///|
/// Result of rrvgo term reduction.
pub(all) struct RvgoReduced {
  term_ids : Array[String]        // original term IDs (input order)
  cluster_id : Array[Int]         // cluster assignment (0-based) per term
  representatives : Array[String] // representative term per cluster
  cluster_sizes : Array[Int]      // size of each cluster
}

///|
/// A single merge step in average-linkage hierarchical clustering.
pub(all) struct RvgoMerge {
  first : Int      // index of first merged cluster (negative = leaf -id-1)
  second : Int     // index of second merged cluster
  height : Double  // distance at which merge happens
}

///|
/// Compute the Wang semantic similarity matrix for a list of GO terms.
/// Returns an n x n matrix (symmetric, diagonal = 1.0).
pub fn rvgo_sim_matrix(graph : GOGraph, term_ids : Array[String]) -> Array[Array[Double]] {
  let n = term_ids.length()
  let mat : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i in 0..<n {
    let row = Array::make(n, 0.0)
    mat[i] = row
    for j in 0..<n {
      mat[i][j] = if i == j { 1.0 } else { gosemsim_wang_sim(graph, term_ids[i], term_ids[j]) }
    }
  }
  mat
}

///|
/// Convert a similarity matrix to a distance matrix (dist = 1 - sim).
pub fn rvgo_dist_matrix(sim : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = sim.length()
  let mat : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i in 0..<n {
    let row = Array::make(n, 0.0)
    mat[i] = row
    for j in 0..<n {
      let d = 1.0 - sim[i][j]
      mat[i][j] = if d < 0.0 { 0.0 } else { d }
    }
  }
  mat
}

///|
/// Average-linkage agglomerative hierarchical clustering.
/// Returns n-1 merge steps (for n input points).
///
/// Internally, active clusters are indexed 0..n-1 for leaves; merged
/// clusters get indices n, n+1, ... A negative index in a merge result
/// refers to leaf `-id - 1`.
pub fn rvgo_hclust(dist : Array[Array[Double]]) -> Array[RvgoMerge] {
  let n = dist.length()
  let merges : Array[RvgoMerge] = Array(capacity=0)
  if n <= 1 {
    return merges
  }
  // Distances between current clusters: Map keyed by (min,max) integer pair.
  // Track cluster membership as a list of leaves per cluster.
  let mut members : Array[Array[Int]] = Array::make(n, [])
  for i in 0..<n {
    members[i] = [i]
  }
  // cluster-to-cluster distance map
  let cdist : Map[Int, Double] = Map([], capacity=128)
  let key_of = fn(a : Int, b : Int) -> Int {
    let lo = if a < b { a } else { b }
    let hi = if a > b { a } else { b }
    if lo == hi { lo } else { lo * 100000 + hi }
  }
  for i in 0..<n {
    let j_start = i + 1
    for j in j_start..<n {
      cdist.set(key_of(i, j), dist[i][j])
    }
  }
  let mut active : Array[Int] = Array::make(n, 0)
  for i in 0..<n {
    active[i] = i
  }
  let mut next_id = n
  while active.length() > 1 {
    // Find the closest pair of active clusters
    let mut best_a = -1
    let mut best_b = -1
    let mut best_d = 1.0e30
    for x in 0..<active.length() {
      let y_start = x + 1
      for y in y_start..<active.length() {
        let ca = active[x]
        let cb = active[y]
        let d = cdist.get(key_of(ca, cb)).unwrap_or(1.0e30)
        if d < best_d {
          best_d = d
          best_a = ca
          best_b = cb
        }
      }
    }
    if best_a < 0 {
      break
    }
    // Merge best_a and best_b into next_id
    merges.push(RvgoMerge::{ first: best_a, second: best_b, height: best_d })
    // Build merged membership
    let merged_members : Array[Int] = Array(capacity=0)
    for m in members_index(members, best_a) {
      merged_members.push(m)
    }
    for m in members_index(members, best_b) {
      merged_members.push(m)
    }
    // Update distance from new cluster to all other active clusters (average)
    for oc in active {
      if oc == best_a || oc == best_b {
        continue
      }
      let da = cdist.get(key_of(best_a, oc)).unwrap_or(best_d)
      let db = cdist.get(key_of(best_b, oc)).unwrap_or(best_d)
      let davg = (da + db) / 2.0
      cdist.set(key_of(next_id, oc), davg)
    }
    // Remove best_a, best_b from active; add next_id
    let new_active : Array[Int] = Array(capacity=0)
    for c in active {
      if c != best_a && c != best_b {
        new_active.push(c)
      }
    }
    new_active.push(next_id)
    active = new_active
    members = members_append(members, next_id, merged_members)
    next_id = next_id + 1
  }
  merges
}

///|
/// Retrieve the leaf members for a cluster id (leaf or merged).
fn members_index(members : Array[Array[Int]], cluster_id : Int) -> Array[Int] {
  if cluster_id >= 0 && cluster_id < members.length() {
    members[cluster_id]
  } else {
    []
  }
}

///|
/// Append a new cluster's member list (grow the members table for merged ids).
fn members_append(
  members : Array[Array[Int]],
  new_id : Int,
  merged_members : Array[Int],
) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = Array::make(new_id + 1, [])
  for i in 0..<members.length() {
    out[i] = members[i]
  }
  // fill the rest up to new_id
  let mut i = members.length()
  while i < new_id {
    out[i] = []
    i = i + 1
  }
  out[new_id] = merged_members
  out
}

///|
/// Cut the hierarchical tree at a similarity threshold.
/// Terms with pairwise similarity >= threshold end up in the same cluster
/// (equivalent to cutting the tree at distance 1 - threshold).
///
/// Returns a cluster label (0-based) for each of the n original terms,
/// in input order. Cluster labels are assigned in order of appearance.
pub fn rvgo_cutree(merges : Array[RvgoMerge], n : Int, threshold : Double) -> Array[Int] {
  // union-find over leaves, merging when merge height <= (1 - threshold)
  let parent : Array[Int] = Array::make(n, 0)
  for i in 0..<n {
    parent[i] = i
  }
  fn find(p : Array[Int], x : Int) -> Int {
    let mut root = x
    while p[root] != root {
      root = p[root]
    }
    // path compression
    let mut c = x
    while p[c] != c {
      let nxt = p[c]
      p[c] = root
      c = nxt
    }
    root
  }
  fn union(p : Array[Int], a : Int, b : Int) -> Unit {
    let ra = find(p, a)
    let rb = find(p, b)
    if ra != rb {
      p[rb] = ra
    }
  }
  let cutoff_dist = 1.0 - threshold
  let cluster_of : Array[Int] = Array::make(n, 0)
  for merge in merges {
    if merge.height <= cutoff_dist {
      // merge leaves of first and second
      for lf in cluster_leaves(merge.first, merges, n) {
        for lr in cluster_leaves(merge.second, merges, n) {
          union(parent, lf, lr)
        }
      }
    }
  }
  // Assign cluster labels in order of first appearance
  let label_map : Map[Int, Int] = Map([], capacity=n)
  let mut next_label = 0
  for i in 0..<n {
    let r = find(parent, i)
    match label_map.get(r) {
      Some(l) => cluster_of[i] = l
      None => {
        label_map.set(r, next_label)
        cluster_of[i] = next_label
        next_label = next_label + 1
      }
    }
  }
  cluster_of
}

///|
/// Determine the leaf indices belonging to a cluster id (leaf or merged).
/// Uses the merge list to expand merged ids into leaves.
fn cluster_leaves(cluster_id : Int, merges : Array[RvgoMerge], n : Int) -> Array[Int] {
  // leaf
  if cluster_id >= 0 && cluster_id < n {
    return [cluster_id]
  }
  // merged id >= n
  let midx = cluster_id - n
  if midx < 0 || midx >= merges.length() {
    return []
  }
  let m = merges[midx]
  let result : Array[Int] = Array(capacity=0)
  for l in cluster_leaves(m.first, merges, n) {
    result.push(l)
  }
  for l in cluster_leaves(m.second, merges, n) {
    result.push(l)
  }
  result
}

///|
/// Reduce redundant GO terms to representatives.
///
/// Parameters:
///   graph: the GO graph (from gosemsim)
///   term_ids: enriched GO term IDs (may include names via score map)
///   scores: significance score per term (smaller = more significant),
///           index-aligned with term_ids
///   threshold: similarity cutoff (default 0.7, per rrvgo convention)
///
/// Returns an RvgoReduced with cluster assignment and representatives.
pub fn rvgo_reduce(
  graph : GOGraph,
  term_ids : Array[String],
  scores : Array[Double],
  threshold? : Double = 0.7,
) -> RvgoReduced {
  let n = term_ids.length()
  if n == 0 {
    return RvgoReduced::{
      term_ids: [], cluster_id: [], representatives: [], cluster_sizes: [],
    }
  }
  let sim = rvgo_sim_matrix(graph, term_ids)
  let dist = rvgo_dist_matrix(sim)
  let merges = rvgo_hclust(dist)
  let cluster_of = if n == 1 {
    [0]
  } else {
    rvgo_cutree(merges, n, threshold)
  }

  // Number of clusters
  let mut n_clusters = 0
  for c in cluster_of {
    if c + 1 > n_clusters {
      n_clusters = c + 1
    }
  }
  if n_clusters == 0 {
    n_clusters = 1
  }

  // Pick the most significant (lowest score) term per cluster as representative
  let reps : Array[String] = Array::make(n_clusters, "")
  let rep_score : Array[Double] = Array::make(n_clusters, 1.0e30)
  let sizes : Array[Int] = Array::make(n_clusters, 0)
  for i in 0..<n {
    let c = cluster_of[i]
    sizes[c] = sizes[c] + 1
    if scores[i] < rep_score[c] {
      rep_score[c] = scores[i]
      reps[c] = term_ids[i]
    }
  }
  RvgoReduced::{
    term_ids: term_ids,
    cluster_id: cluster_of,
    representatives: reps,
    cluster_sizes: sizes,
  }
}

///|
/// Number of clusters in a reduction result.
pub fn rvgo_n_clusters(reduced : RvgoReduced) -> Int {
  reduced.representatives.length()
}

///|
/// Get the representative term for a cluster index.
pub fn rvgo_representative(reduced : RvgoReduced, cluster : Int) -> String {
  if cluster >= 0 && cluster < reduced.representatives.length() {
    reduced.representatives[cluster]
  } else {
    ""
  }
}

///|
/// Get the term IDs assigned to a cluster.
pub fn rvgo_cluster_terms(reduced : RvgoReduced, cluster : Int) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  for i in 0..<reduced.term_ids.length() {
    if reduced.cluster_id[i] == cluster {
      result.push(reduced.term_ids[i])
    }
  }
  result
}

///|
/// Export the reduction result as a TSV string.
/// Columns: term_id, cluster, is_representative
pub fn rvgo_export(reduced : RvgoReduced) -> String {
  let sb = StringBuilder()
  sb.write_string("term_id\tcluster\tis_representative\n")
  for i in 0..<reduced.term_ids.length() {
    let c = reduced.cluster_id[i]
    let is_rep = if reduced.representatives[c] == reduced.term_ids[i] { "yes" } else { "no" }
    sb.write_string(reduced.term_ids[i] + "\t" + c.to_string() + "\t" + is_rep + "\n")
  }
  sb.to_string()
}

///|
/// Classic (metric) multidimensional scaling to 2D coordinates for plotting.
/// Implements PCoA: double-center the squared distance matrix and extract the
/// top two eigenpairs.
pub fn rvgo_mds_2d(dist : Array[Array[Double]]) -> Array[(Double, Double)] {
  let n = dist.length()
  let coords : Array[(Double, Double)] = Array(capacity=0)
  if n == 0 {
    return coords
  }
  if n == 1 {
    coords.push((0.0, 0.0))
    return coords
  }
  // Squared distance matrix D^2
  let d2 : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i in 0..<n {
    let row = Array::make(n, 0.0)
    d2[i] = row
    for j in 0..<n {
      d2[i][j] = dist[i][j] * dist[i][j]
    }
  }
  // Double centering: B = -0.5 * (D2 - row_mean - col_mean + grand_mean)
  let row_mean : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let mut s = 0.0
    for j in 0..<n {
      s = s + d2[i][j]
    }
    row_mean[i] = s / n.to_double()
  }
  let mut grand = 0.0
  for i in 0..<n {
    grand = grand + row_mean[i]
  }
  grand = grand / n.to_double()
  let b : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i in 0..<n {
    let row = Array::make(n, 0.0)
    b[i] = row
    for j in 0..<n {
      b[i][j] = -0.5 * (d2[i][j] - row_mean[i] - row_mean[j] + grand)
    }
  }
  let (evals, evecs) = rvgo_topeigen(b, 2)
  for i in 0..<n {
    let x = evecs[0][i] * evals[0].max(0.0).sqrt()
    let y = if 2 <= n { evecs[1][i] * evals[1].max(0.0).sqrt() } else { 0.0 }
    coords.push((x, y))
  }
  coords
}

///|
/// Extract the top `k` eigenpairs of a symmetric matrix using power iteration
/// with deflation. Internal helper.
fn rvgo_topeigen(
  m : Array[Array[Double]],
  k : Int,
) -> (Array[Double], Array[Array[Double]]) {
  let n = m.length()
  let evals : Array[Double] = Array(capacity=0)
  let evecs : Array[Array[Double]] = Array(capacity=0)
  if n == 0 {
    return (evals, evecs)
  }
  // Work on a copy that gets deflated
  let work : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i in 0..<n {
    let row = Array::make(n, 0.0)
    work[i] = row
    for j in 0..<n {
      work[i][j] = m[i][j]
    }
  }
  let nk = if k < n { k } else { n }
  for _t in 0..<nk {
    let (val, vec) = rvgo_power1(work)
    evals.push(val)
    evecs.push(vec)
    // deflate: work -= val * v v^T
    for i in 0..<n {
      for j in 0..<n {
        work[i][j] = work[i][j] - val * vec[i] * vec[j]
      }
    }
  }
  (evals, evecs)
}

///|
/// Power iteration for the dominant eigenpair of a symmetric matrix.
fn rvgo_power1(m : Array[Array[Double]]) -> (Double, Array[Double]) {
  let n = m.length()
  let vec : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    vec[i] = 1.0 / n.to_double().sqrt()
  }
  let mut val = 0.0
  let mut iter = 0
  while iter < 200 {
    // v_next = M * v
    let nxt : Array[Double] = Array::make(n, 0.0)
    for i in 0..<n {
      let mut s = 0.0
      for j in 0..<n {
        s = s + m[i][j] * vec[j]
      }
      nxt[i] = s
    }
    // normalize
    let mut norm = 0.0
    for x in nxt {
      norm = norm + x * x
    }
    norm = norm.sqrt()
    if norm < 1.0e-15 {
      val = 0.0
      break
    }
    for i in 0..<n {
      vec[i] = nxt[i] / norm
    }
    // Rayleigh quotient
    let mut num = 0.0
    let mut den = 0.0
    for i in 0..<n {
      let mut s = 0.0
      for j in 0..<n {
        s = s + m[i][j] * vec[j]
      }
      num = num + vec[i] * s
      den = den + vec[i] * vec[i]
    }
    val = num / den
    iter = iter + 1
  }
  (val, vec)
}
