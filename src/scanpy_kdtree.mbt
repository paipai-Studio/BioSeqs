///|
/// kd-tree for accelerated exact k-nearest-neighbor search in Euclidean space.
///
/// Replaces the O(n²) brute-force kNN in `scanpy_neighbors`.  For
/// low-dimensional data (d ≲ 20, e.g. PCA-reduced single-cell embeddings)
/// n queries cost O(n·log n) total, versus O(n²) brute-force.  For high
/// dimensions (d > 20) or very small n (≤ 64) the caller falls back to
/// brute-force — kd-tree pruning degrades to O(n) per query when d ≳ 20
/// due to the curse of dimensionality.
///
/// **Correctness guarantee**: the kd-tree kNN produces bit-identical output
/// to brute-force.  Distances are compared as squared distances (monotonic
/// with sqrt), and ties are broken by point index ascending — matching the
/// insertion-sort argsort in the brute-force path.

///|
/// kd-tree node stored in a flat array for cache efficiency.  Children are
/// referenced by index into the node array (-1 = no child).
struct KdNode {
  point_idx : Int       // index into original data
  split_dim : Int       // split dimension (-1 for leaf)
  split_val : Double    // split value at split_dim
  mut left : Int        // left child index (-1 = none)
  mut right : Int       // right child index (-1 = none)
}

///|
/// Build a kd-tree from `n` points in `d` dimensions.
/// Returns a flat array of nodes (root = index 0, or empty if n == 0).
fn scanpy_kdtree_build(
  data : Array[Array[Double]],
  n : Int,
  d : Int,
) -> Array[KdNode] {
  let nodes : Array[KdNode] = Array(capacity=n)
  if n == 0 || d == 0 {
    return nodes
  }
  // Working array of point indices, partitioned in-place during build
  let work : Array[Int] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    work.push(i)
  }
  let _ = scanpy_kdtree_build_rec(data, work, 0, n, d, nodes)
  nodes
}

///|
/// Recursive build: partition `work[lo..hi)` and append nodes to `nodes`.
/// Returns the index of the root node of this subtree in `nodes`, or -1
/// for an empty range.
fn scanpy_kdtree_build_rec(
  data : Array[Array[Double]],
  work : Array[Int],
  lo : Int,
  hi : Int,
  d : Int,
  nodes : Array[KdNode],
) -> Int {
  if lo >= hi {
    return -1
  }
  let count = hi - lo
  // Leaf: single point
  if count == 1 {
    let idx = nodes.length()
    nodes.push(KdNode::{
      point_idx: work[lo],
      split_dim: -1,
      split_val: 0.0,
      left: -1,
      right: -1,
    })
    return idx
  }
  // Find dimension with max spread (best split axis)
  let mut split_dim = 0
  let mut max_spread = 0.0
  for dim = 0; dim < d; dim = dim + 1 {
    let mut min_val = data[work[lo]][dim]
    let mut max_val = min_val
    for t = lo + 1; t < hi; t = t + 1 {
      let v = data[work[t]][dim]
      if v < min_val {
        min_val = v
      }
      if v > max_val {
        max_val = v
      }
    }
    let spread = max_val - min_val
    if spread > max_spread {
      max_spread = spread
      split_dim = dim
    }
  }
  // Sort work[lo..hi) by split_dim to find median
  if max_spread > 0.0 {
    scanpy_kdtree_sort_sub(data, work, lo, hi, split_dim)
  }
  let mid = lo + count / 2
  let median_idx = work[mid]
  let split_val = data[median_idx][split_dim]
  // Create node
  let node_idx = nodes.length()
  nodes.push(KdNode::{
    point_idx: median_idx,
    split_dim,
    split_val,
    left: -1,
    right: -1,
  })
  // Build left subtree [lo, mid) and right subtree [mid+1, hi)
  nodes[node_idx].left = scanpy_kdtree_build_rec(data, work, lo, mid, d, nodes)
  nodes[node_idx].right = scanpy_kdtree_build_rec(data, work, mid + 1, hi, d, nodes)
  node_idx
}

///|
/// Quicksort `work[lo..hi)` by `data[work[i]][dim]` ascending.
/// Uses Lomuto partition with middle-element pivot.
fn scanpy_kdtree_sort_sub(
  data : Array[Array[Double]],
  work : Array[Int],
  lo : Int,
  hi : Int,
  dim : Int,
) -> Unit {
  if lo + 1 >= hi {
    return
  }
  // Lomuto partition with middle pivot
  let pivot_pos = lo + (hi - lo) / 2
  let pivot_val = data[work[pivot_pos]][dim]
  // Swap pivot to hi-1
  let tmp = work[pivot_pos]
  work[pivot_pos] = work[hi - 1]
  work[hi - 1] = tmp
  let mut i = lo
  for j = lo; j < hi - 1; j = j + 1 {
    if data[work[j]][dim] < pivot_val {
      let t = work[i]
      work[i] = work[j]
      work[j] = t
      i = i + 1
    }
  }
  // Swap pivot to final position i
  let t2 = work[i]
  work[i] = work[hi - 1]
  work[hi - 1] = t2
  // Recurse on [lo, i) and [i+1, hi)
  scanpy_kdtree_sort_sub(data, work, lo, i, dim)
  scanpy_kdtree_sort_sub(data, work, i + 1, hi, dim)
}

///|
/// kNN query: find the k nearest neighbors of query point `q` in the
/// kd-tree.  Returns `(indices, squared_distances)` sorted by distance
/// ascending, ties broken by index ascending — matching brute-force
/// argsort semantics.
fn scanpy_kdtree_knn_query(
  nodes : Array[KdNode],
  data : Array[Array[Double]],
  root : Int,
  q : Array[Double],
  k : Int,
  d : Int,
) -> (Array[Int], Array[Double]) {
  // Sorted array of size <= k, ascending by (dist_sq, idx).
  // The last element is the "worst" (largest distance / largest index).
  let heap_dist : Array[Double] = Array(capacity=k)
  let heap_idx : Array[Int] = Array(capacity=k)
  scanpy_kdtree_knn_rec(nodes, data, root, q, k, d, heap_dist, heap_idx)
  (heap_idx, heap_dist)
}

///|
/// Recursive kNN search with branch pruning.
fn scanpy_kdtree_knn_rec(
  nodes : Array[KdNode],
  data : Array[Array[Double]],
  node_idx : Int,
  q : Array[Double],
  k : Int,
  d : Int,
  heap_dist : Array[Double],
  heap_idx : Array[Int],
) -> Unit {
  if node_idx < 0 {
    return
  }
  let node = nodes[node_idx]
  // Compute squared Euclidean distance to this node's point
  let pt = data[node.point_idx]
  let mut dist_sq = 0.0
  for dim = 0; dim < d; dim = dim + 1 {
    let diff = q[dim] - pt[dim]
    dist_sq = dist_sq + diff * diff
  }
  // Try to insert into the candidate heap
  scanpy_kdtree_heap_insert(heap_dist, heap_idx, k, dist_sq, node.point_idx)
  // Leaf node: no children to explore
  if node.split_dim < 0 {
    return
  }
  // Determine which subtree is "near" (contains query point) vs "far"
  let split_dim = node.split_dim
  let diff = q[split_dim] - node.split_val
  let near = if diff < 0.0 { node.left } else { node.right }
  let far = if diff < 0.0 { node.right } else { node.left }
  // Always visit the near subtree first
  scanpy_kdtree_knn_rec(nodes, data, near, q, k, d, heap_dist, heap_idx)
  // Visit far subtree only if it could contain closer points:
  //   - heap not full yet (need k candidates), OR
  //   - distance to splitting plane <= current worst (<= for tie-breaking)
  let heap_size = heap_dist.length()
  if heap_size < k {
    scanpy_kdtree_knn_rec(nodes, data, far, q, k, d, heap_dist, heap_idx)
  } else {
    let diff_sq = diff * diff
    let worst = heap_dist[heap_size - 1]
    if diff_sq <= worst {
      scanpy_kdtree_knn_rec(nodes, data, far, q, k, d, heap_dist, heap_idx)
    }
  }
}

///|
/// Try to insert `(dist_sq, idx)` into the candidate heap.
///
/// The heap is a sorted array of size ≤ k, ordered by `(dist_sq, idx)`
/// ascending.  The last element is the worst (largest distance, or same
/// distance with larger index).  When the heap is full, a new candidate
/// replaces the worst only if it is strictly better — this matches the
/// brute-force argsort tie-breaking (smaller index wins on equal distance).
fn scanpy_kdtree_heap_insert(
  heap_dist : Array[Double],
  heap_idx : Array[Int],
  k : Int,
  dist_sq : Double,
  idx : Int,
) -> Unit {
  let size = heap_dist.length()
  if size < k {
    // Heap not full — insert in sorted position
    // Find insertion position (binary search would be faster, but k is small)
    let mut pos = size
    while pos > 0 {
      let pd = heap_dist[pos - 1]
      let pi = heap_idx[pos - 1]
      // (dist_sq, idx) < (pd, pi) ?
      if dist_sq < pd || (dist_sq == pd && idx < pi) {
        pos = pos - 1
      } else {
        break
      }
    }
    // Shift elements right and insert
    heap_dist.push(0.0)
    heap_idx.push(0)
    let len = heap_dist.length()
    for t = len - 1; t > pos; t = t - 1 {
      heap_dist[t] = heap_dist[t - 1]
      heap_idx[t] = heap_idx[t - 1]
    }
    heap_dist[pos] = dist_sq
    heap_idx[pos] = idx
  } else {
    // Heap full — replace worst only if strictly better
    let worst_dist = heap_dist[size - 1]
    let worst_idx = heap_idx[size - 1]
    // (dist_sq, idx) < (worst_dist, worst_idx) ?
    if dist_sq < worst_dist || (dist_sq == worst_dist && idx < worst_idx) {
      // Find insertion position
      let mut pos = size - 1
      while pos > 0 {
        let pd = heap_dist[pos - 1]
        let pi = heap_idx[pos - 1]
        if dist_sq < pd || (dist_sq == pd && idx < pi) {
          pos = pos - 1
        } else {
          break
        }
      }
      // Shift elements right from pos to size-2, then overwrite pos
      for t = size - 1; t > pos; t = t - 1 {
        heap_dist[t] = heap_dist[t - 1]
        heap_idx[t] = heap_idx[t - 1]
      }
      heap_dist[pos] = dist_sq
      heap_idx[pos] = idx
    }
  }
}

///|
/// kNN distances via kd-tree.  Returns `(indices, distances)` where each
/// row is sorted by distance ascending (index 0 = self at distance 0).
/// Distances are Euclidean (not squared).
fn scanpy_knn_distances_kdtree(
  x : Array[Array[Double]],
  k : Int,
  n : Int,
  d : Int,
) -> (Array[Array[Int]], Array[Array[Double]]) {
  let nodes = scanpy_kdtree_build(x, n, d)
  let idx_all : Array[Array[Int]] = Array(capacity=n)
  let d_all : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let (nbr_idx, nbr_dist_sq) = scanpy_kdtree_knn_query(nodes, x, 0, x[i], k, d)
    // Take sqrt to convert squared distances to Euclidean distances
    let idx_row : Array[Int] = Array(capacity=k)
    let d_row : Array[Double] = Array(capacity=k)
    for t = 0; t < nbr_idx.length(); t = t + 1 {
      idx_row.push(nbr_idx[t])
      d_row.push(nbr_dist_sq[t].sqrt())
    }
    idx_all.push(idx_row)
    d_all.push(d_row)
  }
  (idx_all, d_all)
}
