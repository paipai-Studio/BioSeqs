///|
/// BiocNeighbors - Nearest neighbor search for single-cell analysis.
/// Inspired by Bioconductor BiocNeighbors package.
/// Implements KMKNN (k-means based kNN), Annoy (approximate nearest neighbors),
/// and brute force exact search.

///|
/// Result of a k-nearest neighbor query.
pub struct KNNResult {
  indices : Array[Array[Int]]
  distances : Array[Array[Double]]
  n_query : Int
  k : Int
} derive(Eq, Debug)

///|
/// Neighbor search method enumeration.
pub enum NeighborSearchMethod {
  KMKNN
  Annoy
  BruteForce
} derive(Eq, Debug)

///|
/// Create KMKNN neighbor search method.
pub fn kmknn_method() -> NeighborSearchMethod {
  KMKNN
}

///|
/// Create Annoy neighbor search method.
pub fn annoy_method() -> NeighborSearchMethod {
  Annoy
}

///|
/// Create BruteForce neighbor search method.
pub fn brute_force_method() -> NeighborSearchMethod {
  BruteForce
}

///|
/// Options for k-nearest neighbor search.
pub struct KNNOptions {
  method : NeighborSearchMethod
  k : Int
  distance : String
  n_trees : Int
  search_k : Int
} derive(Eq, Debug)

///|
/// Index result for fast nearest neighbor lookup.
pub struct IndexResult {
  method : NeighborSearchMethod
  n_points : Int
  n_dim : Int
  data : Array[Array[Double]]
  centroids : Array[Array[Double]]
  tree_nodes : Array[AnnoyNode]
} derive(Eq, Debug)

///|
/// A node in the Annoy tree structure.
pub struct AnnoyNode {
  node_id : Int
  is_leaf : Bool
  split_dim : Int
  split_val : Double
  children : Array[Int]
  point_indices : Array[Int]
} derive(Eq, Debug)

///|
/// Create default KNN options with brute force method and k=10.
pub fn knn_options_new() -> KNNOptions {
  KNNOptions::{
    method: BruteForce,
    k: 10,
    distance: "euclidean",
    n_trees: 50,
    search_k: -1,
  }
}

///|
/// Create KNN options with specified method and k value.
pub fn knn_options_with_method(
  method : NeighborSearchMethod,
  k : Int,
) -> KNNOptions {
  let opts = knn_options_new()
  KNNOptions::{ ..opts, method, k }
}

///|
/// Create KNN options with all parameters specified.
pub fn knn_options_full(
  method : NeighborSearchMethod,
  k : Int,
  distance : String,
  n_trees : Int,
  search_k : Int,
) -> KNNOptions {
  KNNOptions::{ method, k, distance, n_trees, search_k }
}

///|
/// Linear congruential generator state for reproducible random numbers.
priv struct LCG {
  state : Int
}

///|
/// Create a new LCG with given seed.
fn lcg_new(seed : Int) -> LCG {
  LCG::{ state: if seed <= 0 { 12345 } else { seed } }
}

///|
/// Generate next random integer from LCG.
fn lcg_next(lcg : LCG) -> Int {
  let next = lcg.state * 1103515245 + 12345
  let _ = LCG::{ state: next }
  next
}

///|
/// Generate random Double in [0.0, 1.0) from LCG.
fn lcg_next_double(lcg : LCG) -> (LCG, Double) {
  let next = lcg_next(lcg)
  let val = (next & 0x7fffffff).to_double() / 2147483648.0
  (LCG::{ state: next }, val)
}

///|
/// Compute distance between two points using the specified method.
/// Dispatches to euclidean_distance, manhattan_distance, or cosine_distance.
pub fn knn_compute_distance(
  p1 : Array[Double],
  p2 : Array[Double],
  method : String,
) -> Double {
  if method == "manhattan" {
    manhattan_distance(p1, p2)
  } else if method == "cosine" {
    cosine_distance(p1, p2)
  } else {
    euclidean_distance(p1, p2)
  }
}

///|
/// Initialize centroids for k-means using random data points.
fn km_init_centroids(
  data : Array[Array[Double]],
  n_points : Int,
  n_dim : Int,
  n_clusters : Int,
) -> Array[Array[Double]] {
  let mut lcg = lcg_new(42)
  let centroids : Array[Array[Double]] = []
  let used_indices = Array::make(n_points, false)
  let mut count = 0
  while count < n_clusters && count < n_points {
    let (new_lcg, rand_val) = lcg_next_double(lcg)
    lcg = new_lcg
    let idx = (rand_val * n_points.to_double()).to_int()
    let actual_idx = if idx >= n_points { n_points - 1 } else { if idx < 0 { 0 } else { idx } }
    if !used_indices[actual_idx] {
      centroids.push(data[actual_idx].copy())
      used_indices[actual_idx] = true
      count = count + 1
    }
  }
  while count < n_clusters {
    centroids.push(Array::make(n_dim, 0.0))
    count = count + 1
  }
  centroids
}

///|
/// Assign each point to the nearest centroid.
fn km_assign_points(
  data : Array[Array[Double]],
  centroids : Array[Array[Double]],
  n_points : Int,
  _n_dim : Int,
  n_clusters : Int,
) -> Array[Int] {
  let assignments = Array::make(n_points, 0)
  for i = 0; i < n_points; i = i + 1 {
    let point = data[i]
    let mut min_dist = 1000000000000000000000000000000000.0
    let mut min_idx = 0
    for c = 0; c < n_clusters; c = c + 1 {
      let dist = euclidean_distance(point, centroids[c])
      if dist < min_dist {
        min_dist = dist
        min_idx = c
      }
    }
    assignments[i] = min_idx
  }
  assignments
}

///|
/// Update centroids as the mean of assigned points.
fn km_update_centroids(
  data : Array[Array[Double]],
  assignments : Array[Int],
  n_points : Int,
  n_dim : Int,
  n_clusters : Int,
) -> Array[Array[Double]] {
  let new_centroids : Array[Array[Double]] = []
  let counts = Array::make(n_clusters, 0)
  for c = 0; c < n_clusters; c = c + 1 {
    new_centroids.push(Array::make(n_dim, 0.0))
  }
  for i = 0; i < n_points; i = i + 1 {
    let cluster = assignments[i]
    counts[cluster] = counts[cluster] + 1
    for d = 0; d < n_dim; d = d + 1 {
      new_centroids[cluster][d] = new_centroids[cluster][d] + data[i][d]
    }
  }
  for c = 0; c < n_clusters; c = c + 1 {
    if counts[c] > 0 {
      for d = 0; d < n_dim; d = d + 1 {
        new_centroids[c][d] = new_centroids[c][d] / counts[c].to_double()
      }
    }
  }
  new_centroids
}

///|
/// Run k-means clustering on data.
fn km_kmeans(
  data : Array[Array[Double]],
  n_points : Int,
  n_dim : Int,
  n_clusters : Int,
  max_iter : Int,
) -> (Array[Array[Double]], Array[Int]) {
  if n_points == 0 {
    return ([], [])
  }
  let actual_clusters = if n_clusters > n_points { n_points } else { n_clusters }
  let mut centroids = km_init_centroids(data, n_points, n_dim, actual_clusters)
  let mut assignments = Array::make(n_points, 0)
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let new_assignments = km_assign_points(data, centroids, n_points, n_dim, actual_clusters)
    let new_centroids = km_update_centroids(data, new_assignments, n_points, n_dim, actual_clusters)
    let mut max_shift = 0.0
    for c = 0; c < actual_clusters; c = c + 1 {
      let shift = euclidean_distance(centroids[c], new_centroids[c])
      if shift > max_shift {
        max_shift = shift
      }
    }
    assignments = new_assignments
    centroids = new_centroids
    if max_shift < 0.0001 {
      break
    }
  }
  (centroids, assignments)
}

///|
/// Build index for fast kNN lookup. Dispatches to the appropriate method.
pub fn build_knn_index(
  data : Array[Array[Double]],
  n_points : Int,
  n_dim : Int,
  options : KNNOptions,
) -> IndexResult {
  match options.method {
    KMKNN => build_kmknn_index(data, n_points, n_dim, options)
    Annoy => build_annoy_index(data, n_points, n_dim, options)
    BruteForce => {
      IndexResult::{
        method: BruteForce,
        n_points,
        n_dim,
        data: data.copy(),
        centroids: [],
        tree_nodes: [],
      }
    }
  }
}

///|
/// Build a KMKNN (k-means based kNN) index.
/// Runs k-means clustering and stores centroids for fast lookup.
pub fn build_kmknn_index(
  data : Array[Array[Double]],
  n_points : Int,
  n_dim : Int,
  _options : KNNOptions,
) -> IndexResult {
  if n_points == 0 {
    return IndexResult::{
      method: KMKNN,
      n_points: 0,
      n_dim,
      data: [],
      centroids: [],
      tree_nodes: [],
    }
  }
  let n_clusters = if n_points <= 10 {
    1
  } else {
    let sqrt_val = ((n_points.to_double()).sqrt()).to_int()
    if sqrt_val < 2 { 2 } else { sqrt_val }
  }
  let (centroids, _) = km_kmeans(data, n_points, n_dim, n_clusters, 20)
  IndexResult::{
    method: KMKNN,
    n_points,
    n_dim,
    data: data.copy(),
    centroids,
    tree_nodes: [],
  }
}

///|
/// Build a single Annoy tree recursively.
fn build_single_annoy_tree(
  data : Array[Array[Double]],
  point_indices : Array[Int],
  n_dim : Int,
  nodes : Array[AnnoyNode],
  node_id_start : Int,
  depth : Int,
  max_depth : Int,
  lcg : LCG,
) -> (Array[AnnoyNode], Int, LCG) {
  let n = point_indices.length()
  let mut nodes_copy = nodes
  let mut current_lcg = lcg
  let mut next_node_id = node_id_start

  if n == 0 {
    return (nodes_copy, -1, current_lcg)
  }

  if n <= 1 || depth >= max_depth {
    let node = AnnoyNode::{
      node_id: next_node_id,
      is_leaf: true,
      split_dim: 0,
      split_val: 0.0,
      children: [],
      point_indices: point_indices.copy(),
    }
    nodes_copy.push(node)
    return (nodes_copy, next_node_id, current_lcg)
  }

  let (new_lcg, rand_dim) = lcg_next_double(current_lcg)
  current_lcg = new_lcg
  let split_dim = (rand_dim * n_dim.to_double()).to_int()
  let actual_dim = if split_dim >= n_dim { n_dim - 1 } else { if split_dim < 0 { 0 } else { split_dim } }

  let (new_lcg2, rand_val) = lcg_next_double(current_lcg)
  current_lcg = new_lcg2

  let mut min_val = data[point_indices[0]][actual_dim]
  let mut max_val = data[point_indices[0]][actual_dim]
  for i = 1; i < n; i = i + 1 {
    let v = data[point_indices[i]][actual_dim]
    if v < min_val { min_val = v }
    if v > max_val { max_val = v }
  }

  let split_val = min_val + rand_val * (max_val - min_val)

  let left_indices : Array[Int] = []
  let right_indices : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    let idx = point_indices[i]
    if data[idx][actual_dim] <= split_val {
      left_indices.push(idx)
    } else {
      right_indices.push(idx)
    }
  }

  if left_indices.length() == 0 || right_indices.length() == 0 {
    let node = AnnoyNode::{
      node_id: next_node_id,
      is_leaf: true,
      split_dim: 0,
      split_val: 0.0,
      children: [],
      point_indices: point_indices.copy(),
    }
    nodes_copy.push(node)
    return (nodes_copy, next_node_id, current_lcg)
  }

  let node_id = next_node_id
  next_node_id = next_node_id + 1

  let (nodes_after_left, left_child_id, lcg_after_left) = build_single_annoy_tree(
    data,
    left_indices,
    n_dim,
    nodes_copy,
    next_node_id,
    depth + 1,
    max_depth,
    current_lcg,
  )
  nodes_copy = nodes_after_left
  current_lcg = lcg_after_left
  next_node_id = if left_child_id >= 0 { left_child_id } else { next_node_id }

  let (nodes_after_right, right_child_id, lcg_after_right) = build_single_annoy_tree(
    data,
    right_indices,
    n_dim,
    nodes_copy,
    next_node_id + 1,
    depth + 1,
    max_depth,
    current_lcg,
  )
  nodes_copy = nodes_after_right
  current_lcg = lcg_after_right

  let node = AnnoyNode::{
    node_id,
    is_leaf: false,
    split_dim: actual_dim,
    split_val,
    children: [left_child_id, right_child_id],
    point_indices: [],
  }

  let idx = nodes_copy.length() - (nodes_after_right.length() - nodes_after_left.length()) - 1
  nodes_copy[idx] = node

  (nodes_copy, node_id, current_lcg)
}

///|
/// Build an Annoy-like approximate nearest neighbor index.
pub fn build_annoy_index(
  data : Array[Array[Double]],
  n_points : Int,
  n_dim : Int,
  options : KNNOptions,
) -> IndexResult {
  if n_points == 0 {
    return IndexResult::{
      method: Annoy,
      n_points: 0,
      n_dim,
      data: [],
      centroids: [],
      tree_nodes: [],
    }
  }

  let n_trees = if options.n_trees <= 0 { 50 } else { options.n_trees }
  let max_depth = if n_points <= 10 {
    5
  } else {
    let d = (@math.log2(n_points.to_double()) + 1.0).to_int()
    if d > 20 { 20 } else { if d < 3 { 3 } else { d } }
  }

  let mut all_nodes : Array[AnnoyNode] = []
  let mut lcg = lcg_new(42)
  let point_indices = Array::make(n_points, 0)
  for i = 0; i < n_points; i = i + 1 {
    point_indices[i] = i
  }

  for t = 0; t < n_trees; t = t + 1 {
    let (new_lcg, _) = lcg_next_double(lcg)
    lcg = new_lcg
    let shuffled : Array[Int] = []
    let mut indices_left = point_indices.copy()
    let mut shuf_lcg = lcg_new(100 + t)
    while indices_left.length() > 0 {
      let (new_lcg2, rand_val) = lcg_next_double(shuf_lcg)
      shuf_lcg = new_lcg2
      let idx = (rand_val * indices_left.length().to_double()).to_int()
      let actual_idx = if idx >= indices_left.length() { indices_left.length() - 1 } else { if idx < 0 { 0 } else { idx } }
      shuffled.push(indices_left[actual_idx])
      let new_left : Array[Int] = []
      for i = 0; i < indices_left.length(); i = i + 1 {
        if i != actual_idx {
          new_left.push(indices_left[i])
        }
      }
      indices_left = new_left
    }

    let tree_start_id = all_nodes.length()
    let (tree_nodes, _, final_lcg) = build_single_annoy_tree(
      data,
      shuffled,
      n_dim,
      all_nodes,
      tree_start_id,
      0,
      max_depth,
      lcg,
    )
    all_nodes = tree_nodes
    lcg = final_lcg
  }

  IndexResult::{
    method: Annoy,
    n_points,
    n_dim,
    data: data.copy(),
    centroids: [],
    tree_nodes: all_nodes,
  }
}

///|
/// Sort indices by their corresponding distances using quicksort.
pub fn knn_quick_sort(
  data : Array[Double],
  indices : Array[Int],
  low : Int,
  high : Int,
) -> (Array[Int], Array[Double]) {
  let mut sorted_indices = indices.copy()
  let mut sorted_distances = data.copy()
  if low < high {
    let mut i = low
    let mut j = high
    let pivot_idx = (low + high) / 2
    let pivot = sorted_distances[pivot_idx]
    let pivot_index = sorted_indices[pivot_idx]
    while i <= j {
      while sorted_distances[i] < pivot || (sorted_distances[i] == pivot && sorted_indices[i] < pivot_index) { i = i + 1 }
      while sorted_distances[j] > pivot || (sorted_distances[j] == pivot && sorted_indices[j] > pivot_index) { j = j - 1 }
      if i <= j {
        let tmp_idx = sorted_indices[i]
        sorted_indices[i] = sorted_indices[j]
        sorted_indices[j] = tmp_idx
        let tmp_dist = sorted_distances[i]
        sorted_distances[i] = sorted_distances[j]
        sorted_distances[j] = tmp_dist
        i = i + 1
        j = j - 1
      }
    }
    let (si1, sd1) = knn_quick_sort(sorted_distances, sorted_indices, low, j)
    sorted_indices = si1
    sorted_distances = sd1
    let (si2, sd2) = knn_quick_sort(sorted_distances, sorted_indices, i, high)
    sorted_indices = si2
    sorted_distances = sd2
  }
  (sorted_indices, sorted_distances)
}

///|
/// Find indices of the k smallest distances.
pub fn knn_find_k_smallest(
  distances : Array[Double],
  k : Int,
) -> Array[Int] {
  let n = distances.length()
  if k >= n {
    let result = Array::make(n, 0)
    for i = 0; i < n; i = i + 1 { result[i] = i }
    return result
  }
  let indices : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    indices.push(i)
  }
  let (sorted_indices, _) = knn_quick_sort(distances, indices, 0, n - 1)
  let result = Array::make(k, 0)
  for i = 0; i < k; i = i + 1 {
    result[i] = sorted_indices[i]
  }
  result
}

///|
/// Brute force exact kNN search. For each query point, compute distance
/// to all data points, sort, and return k nearest.
pub fn knn_brute_force(
  data : Array[Array[Double]],
  n_points : Int,
  _n_dim : Int,
  query_data : Array[Array[Double]],
  n_query : Int,
  k : Int,
  distance : String,
) -> KNNResult {
  if n_query == 0 || n_points == 0 || k == 0 {
    let actual_k = if k > n_points { n_points } else { k }
    return KNNResult::{
      indices: [],
      distances: [],
      n_query,
      k: actual_k,
    }
  }
  let actual_k = if k > n_points { n_points } else { k }
  let all_indices : Array[Array[Int]] = []
  let all_distances : Array[Array[Double]] = []
  for q = 0; q < n_query; q = q + 1 {
    let query_point = query_data[q]
    let dists = Array::make(n_points, 0.0)
    let indices = Array::make(n_points, 0)
    for i = 0; i < n_points; i = i + 1 {
      dists[i] = knn_compute_distance(query_point, data[i], distance)
      indices[i] = i
    }
    let (sorted_indices, sorted_distances) = knn_quick_sort(dists, indices, 0, n_points - 1)
    let top_indices = Array::make(actual_k, 0)
    let top_distances = Array::make(actual_k, 0.0)
    for i = 0; i < actual_k; i = i + 1 {
      top_indices[i] = sorted_indices[i]
      top_distances[i] = sorted_distances[i]
    }
    all_indices.push(top_indices)
    all_distances.push(top_distances)
  }
  KNNResult::{ indices: all_indices, distances: all_distances, n_query, k: actual_k }
}

///|
/// Run kNN query using the specified index. Dispatches to the appropriate method.
pub fn run_knn(
  index : IndexResult,
  query_data : Array[Array[Double]],
  n_query : Int,
  options : KNNOptions,
) -> KNNResult {
  match index.method {
    KMKNN => knn_kmknn(index, query_data, n_query, options)
    Annoy => knn_annoy(index, query_data, n_query, options)
    BruteForce => knn_brute_force(
      index.data,
      index.n_points,
      index.n_dim,
      query_data,
      n_query,
      options.k,
      options.distance,
    )
  }
}

///|
/// KMKNN search: find nearest centroids to query point, search points
/// in those clusters, and expand to neighboring clusters as needed.
pub fn knn_kmknn(
  index : IndexResult,
  query_data : Array[Array[Double]],
  n_query : Int,
  options : KNNOptions,
) -> KNNResult {
  if n_query == 0 || index.n_points == 0 || options.k == 0 {
    let actual_k = if options.k > index.n_points { index.n_points } else { options.k }
    return KNNResult::{
      indices: [],
      distances: [],
      n_query,
      k: actual_k,
    }
  }
  let actual_k = if options.k > index.n_points { index.n_points } else { options.k }
  let n_clusters = index.centroids.length()
  let all_indices : Array[Array[Int]] = []
  let all_distances : Array[Array[Double]] = []

  let point_cluster = Array::make(index.n_points, 0)
  for i = 0; i < index.n_points; i = i + 1 {
    let mut min_c = 0
    let mut min_d = 1000000000000000000000000000000000.0
    for c = 0; c < n_clusters; c = c + 1 {
      let d = euclidean_distance(index.data[i], index.centroids[c])
      if d < min_d {
        min_d = d
        min_c = c
      }
    }
    point_cluster[i] = min_c
  }

  let cluster_points : Array[Array[Int]] = []
  for c = 0; c < n_clusters; c = c + 1 {
    cluster_points.push([])
  }
  for i = 0; i < index.n_points; i = i + 1 {
    cluster_points[point_cluster[i]].push(i)
  }

  for q = 0; q < n_query; q = q + 1 {
    let query = query_data[q]
    let centroid_dists = Array::make(n_clusters, 0.0)
    for c = 0; c < n_clusters; c = c + 1 {
      centroid_dists[c] = euclidean_distance(query, index.centroids[c])
    }
    let sorted_centroid_indices = Array::make(n_clusters, 0)
    for c = 0; c < n_clusters; c = c + 1 { sorted_centroid_indices[c] = c }
    let (sorted_c_idx, _) = knn_quick_sort(centroid_dists, sorted_centroid_indices, 0, n_clusters - 1)

    let candidates : Array[Int] = []
    let mut cluster_idx = 0
    while candidates.length() < actual_k && cluster_idx < n_clusters {
      let c = sorted_c_idx[cluster_idx]
      for p = 0; p < cluster_points[c].length(); p = p + 1 {
        candidates.push(cluster_points[c][p])
      }
      cluster_idx = cluster_idx + 1
    }

    if candidates.length() < actual_k {
      for i = 0; i < index.n_points; i = i + 1 {
        let mut found = false
        for c = 0; c < candidates.length(); c = c + 1 {
          if candidates[c] == i { found = true; break }
        }
        if !found {
          candidates.push(i)
        }
      }
    }

    let cand_dists = Array::make(candidates.length(), 0.0)
    for i = 0; i < candidates.length(); i = i + 1 {
      cand_dists[i] = knn_compute_distance(query, index.data[candidates[i]], options.distance)
    }
    let (sorted_idx, sorted_dist) = knn_quick_sort(cand_dists, candidates, 0, candidates.length() - 1)

    let top_indices = Array::make(actual_k, 0)
    let top_distances = Array::make(actual_k, 0.0)
    for i = 0; i < actual_k; i = i + 1 {
      top_indices[i] = sorted_idx[i]
      top_distances[i] = sorted_dist[i]
    }
    all_indices.push(top_indices)
    all_distances.push(top_distances)
  }

  KNNResult::{ indices: all_indices, distances: all_distances, n_query, k: actual_k }
}

///|
/// Annoy approximate kNN search. Traverses all trees, visits search_k
/// nodes per tree, and returns approximate k nearest neighbors.
pub fn knn_annoy(
  index : IndexResult,
  query_data : Array[Array[Double]],
  n_query : Int,
  options : KNNOptions,
) -> KNNResult {
  if n_query == 0 || index.n_points == 0 || options.k == 0 {
    let actual_k = if options.k > index.n_points { index.n_points } else { options.k }
    return KNNResult::{
      indices: [],
      distances: [],
      n_query,
      k: actual_k,
    }
  }
  let actual_k = if options.k > index.n_points { index.n_points } else { options.k }
  let n_trees = options.n_trees
  let nodes_per_tree = if n_trees > 0 { index.tree_nodes.length() / n_trees } else { 0 }
  let all_indices : Array[Array[Int]] = []
  let all_distances : Array[Array[Double]] = []

  for q = 0; q < n_query; q = q + 1 {
    let query = query_data[q]
    let candidates : Array[Int] = []

    for t = 0; t < n_trees; t = t + 1 {
      let tree_base = t * nodes_per_tree
      let search_k = if options.search_k <= 0 {
        actual_k * n_trees
      } else {
        options.search_k
      }
      let mut nodes_to_visit : Array[(Int, Double)] = []
      let root_node = index.tree_nodes[tree_base]
      let root_dist = if root_node.is_leaf {
        0.0
      } else {
        let dim = root_node.split_dim
        (query[dim] - root_node.split_val).abs()
      }
      nodes_to_visit.push((tree_base, root_dist))
      let mut visited_in_tree = 0

      while nodes_to_visit.length() > 0 && visited_in_tree < search_k {
        let (node_idx, _) = nodes_to_visit[0]
        let node = index.tree_nodes[node_idx]
        let new_visit : Array[(Int, Double)] = []
        for i = 1; i < nodes_to_visit.length(); i = i + 1 {
          new_visit.push(nodes_to_visit[i])
        }
        nodes_to_visit = new_visit

        visited_in_tree = visited_in_tree + 1

        if node.is_leaf {
          for pi = 0; pi < node.point_indices.length(); pi = pi + 1 {
            let point_idx = node.point_indices[pi]
            let mut already = false
            for c = 0; c < candidates.length(); c = c + 1 {
              if candidates[c] == point_idx { already = true; break }
            }
            if !already {
              candidates.push(point_idx)
            }
          }
        } else {
          let dim = node.split_dim
          let val = node.split_val
          let diff = query[dim] - val
          let (near_child, far_child) = if diff <= 0.0 {
            (node.children[0], node.children[1])
          } else {
            (node.children[1], node.children[0])
          }
          if near_child >= 0 && near_child < index.tree_nodes.length() {
            let near_node = index.tree_nodes[near_child]
            let near_dist = if near_node.is_leaf {
              0.0
            } else {
              let nd = near_node.split_dim
              (query[nd] - near_node.split_val).abs()
            }
            nodes_to_visit.push((near_child, near_dist))
          }
          if far_child >= 0 && far_child < index.tree_nodes.length() {
            let far_node = index.tree_nodes[far_child]
            let far_dist = if far_node.is_leaf {
              0.0
            } else {
              let fd = far_node.split_dim
              (query[fd] - far_node.split_val).abs()
            }
            nodes_to_visit.push((far_child, far_dist))
          }
          let mut i = nodes_to_visit.length() - 1
          while i > 0 {
            let key = nodes_to_visit[i]
            let mut j = i - 1
            while j >= 0 {
              if nodes_to_visit[j].1 > key.1 {
                nodes_to_visit[j + 1] = nodes_to_visit[j]
                nodes_to_visit[j] = key
              } else {
                break
              }
              j = j - 1
            }
            i = i - 1
          }
        }
      }
    }

    if candidates.length() > actual_k {
      let cand_dists = Array::make(candidates.length(), 0.0)
      for i = 0; i < candidates.length(); i = i + 1 {
        cand_dists[i] = knn_compute_distance(query, index.data[candidates[i]], options.distance)
      }
      let (sorted_idx, sorted_dist) = knn_quick_sort(cand_dists, candidates, 0, candidates.length() - 1)
      let top_indices = Array::make(actual_k, 0)
      let top_distances = Array::make(actual_k, 0.0)
      for i = 0; i < actual_k; i = i + 1 {
        top_indices[i] = sorted_idx[i]
        top_distances[i] = sorted_dist[i]
      }
      all_indices.push(top_indices)
      all_distances.push(top_distances)
    } else if candidates.length() > 0 {
      let cand_dists = Array::make(candidates.length(), 0.0)
      for i = 0; i < candidates.length(); i = i + 1 {
        cand_dists[i] = knn_compute_distance(query, index.data[candidates[i]], options.distance)
      }
      let (sorted_idx, sorted_dist) = knn_quick_sort(cand_dists, candidates, 0, candidates.length() - 1)
      let top_indices = Array::make(sorted_idx.length(), 0)
      let top_distances = Array::make(sorted_dist.length(), 0.0)
      for i = 0; i < sorted_idx.length(); i = i + 1 {
        top_indices[i] = sorted_idx[i]
        top_distances[i] = sorted_dist[i]
      }
      all_indices.push(top_indices)
      all_distances.push(top_distances)
    } else {
      all_indices.push([])
      all_distances.push([])
    }
  }

  KNNResult::{ indices: all_indices, distances: all_distances, n_query, k: actual_k }
}