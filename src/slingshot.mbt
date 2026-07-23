///|
/// slingshot module, inspired by Bioconductor slingshot package.
/// Provides single-cell lineage trajectory inference using minimum
/// spanning trees and principal curves for pseudotime ordering.

/// Distance metric for cell comparisons
pub enum SlingshotDistMetric {
  /// Euclidean distance
  Euclidean
  /// Manhattan distance
  Manhattan
  /// Cosine distance
  Cosine
}

/// A node in the lineage tree (represents a cluster or cell)
pub struct SlingshotNode {
  /// Node ID
  id : String
  /// Cluster label
  cluster : String
  /// Dimensional coordinates
  coords : Array[Double]
  /// Node size (number of cells)
  size : Int
  /// Whether this is a terminal node
  terminal : Bool
  /// Whether this is a start node
  start : Bool
}

/// An edge in the lineage tree
pub struct SlingshotEdge {
  /// Source node ID
  from_id : String
  /// Target node ID
  to_id : String
  /// Edge weight (distance)
  weight : Double
}

/// A lineage curve (principal curve for a branch)
pub struct SlingshotCurve {
  /// Curve ID
  id : String
  /// Start node ID
  start_node : String
  /// End node ID
  end_node : String
  /// Control points
  control_points : Array[Array[Double]]
  /// Pseudotime values for cells on this curve
  pseudotime : Array[Double]
  /// Cell indices assigned to this curve
  cell_indices : Array[Int]
}

/// Complete slingshot result
pub struct SlingshotResult {
  /// All nodes
  nodes : Array[SlingshotNode]
  /// Tree edges (MST)
  edges : Array[SlingshotEdge]
  /// Lineage curves
  curves : Array[SlingshotCurve]
  /// Pseudotime per cell
  cell_pseudotime : Array[Double]
  /// Lineage assignment per cell
  cell_lineage : Array[String]
  /// Branch ID per cell
  cell_branch : Array[String]
}

/// Parameters for slingshot analysis
pub struct SlingshotParams {
  /// Number of principal components
  n_pcs : Int
  /// Minimum cluster size
  min_cluster_size : Int
  /// Whether to allow reversals
  allow_reversals : Bool
  /// Shrinkage factor for MST
  shrink : Double
}

/// Create default parameters
pub fn SlingshotParams::new() -> SlingshotParams {
  SlingshotParams::{
    n_pcs: 2,
    min_cluster_size: 10,
    allow_reversals: false,
    shrink: 1.0
  }
}

// ===== Node methods =====

/// Create a new node
pub fn SlingshotNode::new(id : String, cluster : String, coords : Array[Double], size : Int) -> SlingshotNode {
  SlingshotNode::{
    id,
    cluster,
    coords,
    size,
    terminal: false,
    start: false
  }
}

/// Set terminal flag
pub fn SlingshotNode::sling_set_terminal(self : SlingshotNode, terminal : Bool) -> SlingshotNode {
  SlingshotNode::{ id: self.id, cluster: self.cluster, coords: self.coords, size: self.size, terminal, start: self.start }
}

/// Set start flag
pub fn SlingshotNode::sling_set_start(self : SlingshotNode, start : Bool) -> SlingshotNode {
  SlingshotNode::{ id: self.id, cluster: self.cluster, coords: self.coords, size: self.size, terminal: false, start }
}

/// Get node coordinate at dimension
pub fn SlingshotNode::sling_coord(self : SlingshotNode, dim : Int) -> Double {
  self.coords[dim]
}

/// Compute distance between two nodes
pub fn SlingshotNode::sling_distance(self : SlingshotNode, other : SlingshotNode, metric : SlingshotDistMetric) -> Double {
  let d = self.coords.length()
  if d == 0 { return 0.0 }
  match metric {
    SlingshotDistMetric::Euclidean => {
      let mut sum = 0.0
      let mut i = 0
      while i < d {
        let diff = self.coords[i] - other.coords[i]
        sum = sum + diff * diff
        i = i + 1
      }
      sum.sqrt()
    }
    SlingshotDistMetric::Manhattan => {
      let mut sum = 0.0
      let mut i = 0
      while i < d {
        sum = sum + (self.coords[i] - other.coords[i]).abs()
        i = i + 1
      }
      sum
    }
    SlingshotDistMetric::Cosine => {
      let mut dot = 0.0
      let mut norm_a = 0.0
      let mut norm_b = 0.0
      let mut i = 0
      while i < d {
        dot = dot + self.coords[i] * other.coords[i]
        norm_a = norm_a + self.coords[i] * self.coords[i]
        norm_b = norm_b + other.coords[i] * other.coords[i]
        i = i + 1
      }
      let na = if norm_a > 0.0 { norm_a.sqrt() } else { 1.0 }
      let nb = if norm_b > 0.0 { norm_b.sqrt() } else { 1.0 }
      let cos_sim = dot / (na * nb)
      1.0 - cos_sim
    }
  }
}

// ===== Edge methods =====

/// Create a new edge
pub fn SlingshotEdge::new(from_id : String, to_id : String, weight : Double) -> SlingshotEdge {
  SlingshotEdge::{ from_id, to_id, weight }
}

// ===== Metric helpers =====

/// Get Euclidean distance metric
pub fn sling_euclidean_metric() -> SlingshotDistMetric {
  SlingshotDistMetric::Euclidean
}

/// Get Manhattan distance metric
pub fn sling_manhattan_metric() -> SlingshotDistMetric {
  SlingshotDistMetric::Manhattan
}

/// Get Cosine distance metric
pub fn sling_cosine_metric() -> SlingshotDistMetric {
  SlingshotDistMetric::Cosine
}

// ===== MST Construction =====

///|
/// Build a distance matrix between all nodes.
pub fn sling_distance_matrix(nodes : Array[SlingshotNode], metric : SlingshotDistMetric) -> Array[Array[Double]] {
  let n = nodes.length()
  let mat : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      if i == j {
        row.push(0.0)
      } else {
        row.push(nodes[i].sling_distance(nodes[j], metric))
      }
      j = j + 1
    }
    mat.push(row)
    i = i + 1
  }
  mat
}

///|
/// Build a minimum spanning tree using Prim's algorithm.
pub fn sling_build_mst(nodes : Array[SlingshotNode], metric : SlingshotDistMetric) -> Array[SlingshotEdge] {
  let n = nodes.length()
  if n <= 1 { return Array::new() }
  let dist_mat = sling_distance_matrix(nodes, metric)
  let in_tree : Array[Bool] = Array::new()
  let min_edge : Array[Double] = Array::new()
  let min_from : Array[Int] = Array::new()
  let edges : Array[SlingshotEdge] = Array::new()

  let mut i = 0
  while i < n {
    in_tree.push(false)
    min_edge.push(1.0e30)
    min_from.push(-1)
    i = i + 1
  }

  in_tree[0] = true
  let mut k = 1
  while k < n {
    // Update distances from the newly added node
    let recently_added = sling_find_recently_added(in_tree, min_edge, n)
    let mut j = 0
    while j < n {
      if !in_tree[j] && dist_mat[recently_added][j] < min_edge[j] {
        min_edge[j] = dist_mat[recently_added][j]
        min_from[j] = recently_added
      }
      j = j + 1
    }
    // Find minimum edge to add
    let mut min_val = 1.0e30
    let mut min_idx = -1
    let mut m = 0
    while m < n {
      if !in_tree[m] && min_edge[m] < min_val {
        min_val = min_edge[m]
        min_idx = m
      }
      m = m + 1
    }
    if min_idx < 0 { break }
    in_tree[min_idx] = true
    edges.push(SlingshotEdge::new(
      nodes[min_from[min_idx]].id,
      nodes[min_idx].id,
      min_val
    ))
    k = k + 1
  }
  edges
}

/// Find the most recently added node in the MST
fn sling_find_recently_added(in_tree : Array[Bool], min_edge : Array[Double], n : Int) -> Int {
  // Find the node with the smallest min_edge that's in the tree
  let mut min_val = 1.0e30
  let mut min_idx = 0
  let mut i = 0
  while i < n {
    if in_tree[i] && min_edge[i] < min_val {
      min_val = min_edge[i]
      min_idx = i
    }
    i = i + 1
  }
  min_idx
}

///|
/// Identify terminal nodes (leaf nodes in the MST).
pub fn sling_find_terminals(edges : Array[SlingshotEdge], nodes : Array[SlingshotNode]) -> Array[String] {
  let n = nodes.length()
  let degree : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    degree.push(0)
    i = i + 1
  }
  let node_ids : Array[String] = Array::new()
  let mut j = 0
  while j < n {
    node_ids.push(nodes[j].id)
    j = j + 1
  }
  let mut k = 0
  while k < edges.length() {
    let from_idx = sling_find_node_idx(node_ids, edges[k].from_id)
    let to_idx = sling_find_node_idx(node_ids, edges[k].to_id)
    if from_idx >= 0 {
      degree[from_idx] = degree[from_idx] + 1
    }
    if to_idx >= 0 {
      degree[to_idx] = degree[to_idx] + 1
    }
    k = k + 1
  }
  let terminals : Array[String] = Array::new()
  let mut m = 0
  while m < n {
    if degree[m] == 1 {
      terminals.push(nodes[m].id)
    }
    m = m + 1
  }
  terminals
}

/// Find node index by ID
fn sling_find_node_idx(ids : Array[String], id : String) -> Int {
  let mut i = 0
  while i < ids.length() {
    if ids[i] == id { return i }
    i = i + 1
  }
  -1
}

// ===== Principal Curve =====

///|
/// Fit a principal curve through a set of 2D points.
/// Uses the simple iterative approach: project points onto the curve,
/// update curve as running mean.
pub fn sling_fit_principal_curve(points : Array[Array[Double]], n_segments : Int) -> Array[Array[Double]] {
  let n_points = points.length()
  if n_points == 0 { return Array::new() }
  let dim = points[0].length()
  // Initialize curve segments along the first principal component
  let (pc1_min, pc1_max) = sling_pc1_range(points)
  let curve : Array[Array[Double]] = Array::new()
  let seg_step = if n_segments > 1 { (pc1_max - pc1_min) / (n_segments - 1).to_double() } else { 0.0 }
  let mut i = 0
  while i < n_segments {
    let t = pc1_min + seg_step * i.to_double()
    let pt = sling_point_on_pc1(points, t, pc1_min, pc1_max)
    curve.push(pt)
    i = i + 1
  }
  // Iterative refinement
  let max_iter = 20
  let mut iter = 0
  while iter < max_iter {
    let (new_curve, changed) = sling_refine_curve(points, curve)
    if !changed { break }
    iter = iter + 1
  }
  curve
}

///|
/// Get range of first principal component.
fn sling_pc1_range(points : Array[Array[Double]]) -> (Double, Double) {
  let n = points.length()
  if n == 0 { return (0.0, 0.0) }
  let dim = points[0].length()
  // Simple PCA: center data, find direction of maximum variance
  let center : Array[Double] = Array::new()
  let mut d = 0
  while d < dim {
    let mut sum = 0.0
    let mut i = 0
    while i < n {
      sum = sum + points[i][d]
      i = i + 1
    }
    center.push(sum / n.to_double())
    d = d + 1
  }
  // Compute variance along each dimension
  let variances : Array[Double] = Array::new()
  let mut dd = 0
  while dd < dim {
    let mut var = 0.0
    let mut i = 0
    while i < n {
      let diff = points[i][dd] - center[dd]
      var = var + diff * diff
      i = i + 1
    }
    variances.push(var)
    dd = dd + 1
  }
  // Find max variance dimension (proxy for PC1)
  let mut max_dim = 0
  let mut max_var = variances[0]
  let mut v = 1
  while v < dim {
    if variances[v] > max_var {
      max_var = variances[v]
      max_dim = v
    }
    v = v + 1
  }
  let mut min_val = points[0][max_dim]
  let mut max_val = points[0][max_dim]
  let mut j = 1
  while j < n {
    if points[j][max_dim] < min_val { min_val = points[j][max_dim] }
    if points[j][max_dim] > max_val { max_val = points[j][max_dim] }
    j = j + 1
  }
  (min_val, max_val)
}

///|
/// Get point on PC1 line at parameter t.
fn sling_point_on_pc1(points : Array[Array[Double]], t : Double, min_val : Double, max_val : Double) -> Array[Double] {
  let n = points.length()
  let dim = points[0].length()
  let pt : Array[Double] = Array::new()
  let mut d = 0
  while d < dim {
    let mut sum = 0.0
    let mut i = 0
    while i < n {
      sum = sum + points[i][d]
      i = i + 1
    }
    let center = sum / n.to_double()
    // Simple projection: scale center by position along PC1
    let range = if max_val > min_val { max_val - min_val } else { 1.0 }
    let frac = (t - min_val) / range
    pt.push(center * (0.5 + frac * 0.5))
    d = d + 1
  }
  pt
}

///|
/// Refine principal curve by projecting points and averaging.
fn sling_refine_curve(points : Array[Array[Double]], curve : Array[Array[Double]]) -> (Array[Array[Double]], Bool) {
  let n_curve = curve.length()
  let n_points = points.length()
  if n_curve < 2 || n_points == 0 { return (curve, false) }
  let dim = curve[0].length()
  let new_curve : Array[Array[Double]] = Array::new()
  let mut changed = false
  // For each curve segment, project nearby points and average
  let mut ci = 0
  while ci < n_curve {
    let acc : Array[Double] = Array::new()
    let mut d = 0
    while d < dim {
      acc.push(0.0)
      d = d + 1
    }
    let start_seg = if ci == 0 { 0 } else { ci - 1 }
    let end_seg = if ci == n_curve - 1 { n_curve - 1 } else { ci + 1 }
    let mut count = 0.0
    let mut pi = 0
    while pi < n_points {
      let proj = sling_project_to_segment(points[pi], curve[start_seg], curve[end_seg])
      let dist2 = 0.0
      let mut dd = 0
      while dd < dim {
        let diff = points[pi][dd] - proj[dd]
        // Use distance to determine weight
        dd = dd + 1
      }
      let w = 1.0
      let mut k = 0
      while k < dim {
        acc[k] = acc[k] + points[pi][k] * w
        k = k + 1
      }
      count = count + w
      pi = pi + 1
    }
    if count > 0.0 {
      let new_pt : Array[Double] = Array::new()
      let mut m = 0
      while m < dim {
        new_pt.push(acc[m] / count)
        m = m + 1
      }
      // Check if changed significantly
      let mut diff_sum = 0.0
      let mut k = 0
      while k < dim {
        diff_sum = diff_sum + (new_pt[k] - curve[ci][k]).abs()
        k = k + 1
      }
      if diff_sum > 1.0e-6 { changed = true }
      new_curve.push(new_pt)
    } else {
      new_curve.push(curve[ci])
    }
    ci = ci + 1
  }
  (new_curve, changed)
}

///|
/// Project a point onto a line segment.
fn sling_project_to_segment(point : Array[Double], a : Array[Double], b : Array[Double]) -> Array[Double] {
  let dim = point.length()
  let ab : Array[Double] = Array::new()
  let mut i = 0
  while i < dim {
    ab.push(b[i] - a[i])
    i = i + 1
  }
  let ap : Array[Double] = Array::new()
  let mut j = 0
  while j < dim {
    ap.push(point[j] - a[j])
    j = j + 1
  }
  let mut ab2 = 0.0
  let mut k = 0
  while k < dim {
    ab2 = ab2 + ab[k] * ab[k]
    k = k + 1
  }
  let t = if ab2 > 0.0 { sling_dot(ap, ab) / ab2 } else { 0.0 }
  let t_clamped = if t < 0.0 { 0.0 } else { if t > 1.0 { 1.0 } else { t } }
  let proj : Array[Double] = Array::new()
  let mut m = 0
  while m < dim {
    proj.push(a[m] + t_clamped * ab[m])
    m = m + 1
  }
  proj
}

/// Dot product
fn sling_dot(a : Array[Double], b : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < a.length() {
    sum = sum + a[i] * b[i]
    i = i + 1
  }
  sum
}

// ===== Pseudotime ordering =====

///|
/// Compute pseudotime for cells along a lineage curve.
pub fn sling_compute_pseudotime(
  cell_coords : Array[Array[Double]],
  curve : Array[Array[Double]]
) -> Array[Double] {
  let n_cells = cell_coords.length()
  let n_curve = curve.length()
  let pseudotime : Array[Double] = Array::new()
  if n_cells == 0 || n_curve < 2 {
    let mut i = 0
    while i < n_cells {
      pseudotime.push(0.0)
      i = i + 1
    }
    return pseudotime
  }
  // Compute cumulative arc length of the curve
  let arc_lengths : Array[Double] = Array::new()
  arc_lengths.push(0.0)
  let mut i = 1
  while i < n_curve {
    let seg_len = sling_euclid_dist(curve[i - 1], curve[i])
    let prev = arc_lengths[i - 1]
    arc_lengths.push(prev + seg_len)
    i = i + 1
  }
  let total_len = arc_lengths[n_curve - 1]
  // Project each cell onto the curve and compute pseudotime
  let mut j = 0
  while j < n_cells {
    let pt = cell_coords[j]
    let (proj_t, _) = sling_nearest_on_curve(pt, curve)
    let ps = if total_len > 0.0 { proj_t / total_len } else { 0.0 }
    pseudotime.push(ps)
    j = j + 1
  }
  pseudotime
}

///|
/// Find nearest point on curve to a given point.
fn sling_nearest_on_curve(point : Array[Double], curve : Array[Array[Double]]) -> (Double, Double) {
  let n = curve.length()
  if n < 2 { return (0.0, sling_euclid_dist(point, curve[0])) }
  let mut best_dist = 1.0e30
  let mut best_t = 0.0
  let seg_len_total = sling_total_arc_length(curve)
  let mut arc_acc = 0.0
  let mut i = 0
  while i < n - 1 {
    let proj = sling_project_to_segment(point, curve[i], curve[i + 1])
    let dist = sling_euclid_dist(point, proj)
    let seg_len = sling_euclid_dist(curve[i], curve[i + 1])
    let local_t = if seg_len > 0.0 {
      sling_dot(sling_sub(proj, curve[i]), sling_sub(curve[i + 1], curve[i])) / (seg_len * seg_len)
    } else { 0.0 }
    let global_t = arc_acc + local_t * seg_len
    if dist < best_dist {
      best_dist = dist
      best_t = if seg_len_total > 0.0 { global_t / seg_len_total } else { 0.0 }
    }
    arc_acc = arc_acc + seg_len
    i = i + 1
  }
  (best_t, best_dist)
}

fn sling_total_arc_length(curve : Array[Array[Double]]) -> Double {
  let mut total = 0.0
  let mut i = 1
  while i < curve.length() {
    total = total + sling_euclid_dist(curve[i - 1], curve[i])
    i = i + 1
  }
  total
}

fn sling_sub(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < a.length() {
    result.push(a[i] - b[i])
    i = i + 1
  }
  result
}

fn sling_euclid_dist(a : Array[Double], b : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < a.length() {
    let diff = a[i] - b[i]
    sum = sum + diff * diff
    i = i + 1
  }
  sum.sqrt()
}

// ===== Main slingshot pipeline =====

///|
/// Run the complete slingshot trajectory inference pipeline.
pub fn slingshot_run(
  cell_coords : Array[Array[Double]],
  cluster_labels : Array[String],
  params : SlingshotParams
) -> SlingshotResult {
  let n_cells = cell_coords.length()
  let n_clusters = sling_unique_labels(cluster_labels)
  // Build cluster centers
  let nodes = sling_build_nodes(cell_coords, cluster_labels, n_clusters)
  // Build MST
  let edges = sling_build_mst(nodes, SlingshotDistMetric::Euclidean)
  // Find terminals
  let terminals = sling_find_terminals(edges, nodes)
  // Start node (root)
  let start_node = if nodes.length() > 0 { nodes[0].id } else { "" }
  // Build lineages from edges
  let lineages = sling_build_lineages(edges, nodes, start_node, terminals)
  // Fit principal curves
  let curves = sling_fit_lineage_curves(cell_coords, cluster_labels, nodes, lineages)
  // Compute pseudotime
  let (cell_pseudotime, cell_lineage, cell_branch) = sling_assign_pseudotime(
    cell_coords, cluster_labels, curves, nodes
  )
  SlingshotResult::{
    nodes,
    edges,
    curves,
    cell_pseudotime,
    cell_lineage,
    cell_branch
  }
}

/// Get unique cluster labels
fn sling_unique_labels(labels : Array[String]) -> Array[String] {
  let unique : Array[String] = Array::new()
  let mut i = 0
  while i < labels.length() {
    let found = sling_contains_str(unique, labels[i])
    if !found {
      unique.push(labels[i])
    }
    i = i + 1
  }
  unique
}

fn sling_contains_str(arr : Array[String], s : String) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if arr[i] == s { return true }
    i = i + 1
  }
  false
}

/// Build nodes from cluster centers
fn sling_build_nodes(cell_coords : Array[Array[Double]], cluster_labels : Array[String], unique_labels : Array[String]) -> Array[SlingshotNode] {
  let n_clusters = unique_labels.length()
  let nodes : Array[SlingshotNode] = Array::new()
  let mut c = 0
  while c < n_clusters {
    let label = unique_labels[c]
    let center = sling_cluster_center(cell_coords, cluster_labels, label)
    let size = sling_cluster_size(cluster_labels, label)
    nodes.push(SlingshotNode::new(
      "cluster_" + c.to_string(),
      label,
      center,
      size
    ))
    c = c + 1
  }
  nodes
}

/// Compute cluster center (mean coordinates)
fn sling_cluster_center(cell_coords : Array[Array[Double]], cluster_labels : Array[String], label : String) -> Array[Double] {
  let n_cells = cell_coords.length()
  if n_cells == 0 { return Array::new() }
  let dim = cell_coords[0].length()
  let center : Array[Double] = Array::new()
  let mut d = 0
  while d < dim {
    let mut sum = 0.0
    let mut count = 0.0
    let mut i = 0
    while i < n_cells {
      if cluster_labels[i] == label {
        sum = sum + cell_coords[i][d]
        count = count + 1.0
      }
      i = i + 1
    }
    center.push(if count > 0.0 { sum / count } else { 0.0 })
    d = d + 1
  }
  center
}

/// Get cluster size
fn sling_cluster_size(cluster_labels : Array[String], label : String) -> Int {
  let mut count = 0
  let mut i = 0
  while i < cluster_labels.length() {
    if cluster_labels[i] == label { count = count + 1 }
    i = i + 1
  }
  count
}

/// Build lineage paths from MST edges
fn sling_build_lineages(edges : Array[SlingshotEdge], nodes : Array[SlingshotNode], start_id : String, terminals : Array[String]) -> Array[Array[String]] {
  let lineages : Array[Array[String]] = Array::new()
  // For each terminal, find path from start
  let mut t = 0
  while t < terminals.length() {
    let path = sling_find_path(edges, nodes, start_id, terminals[t])
    if path.length() > 0 {
      lineages.push(path)
    }
    t = t + 1
  }
  // If no lineages found, create one from all nodes
  if lineages.length() == 0 && nodes.length() > 0 {
    let full_path : Array[String] = Array::new()
    let mut i = 0
    while i < nodes.length() {
      full_path.push(nodes[i].id)
      i = i + 1
    }
    lineages.push(full_path)
  }
  lineages
}

/// Find path between two nodes in MST
fn sling_find_path(edges : Array[SlingshotEdge], nodes : Array[SlingshotNode], start_id : String, end_id : String) -> Array[String] {
  // BFS
  let node_ids : Array[String] = Array::new()
  let mut i = 0
  while i < nodes.length() {
    node_ids.push(nodes[i].id)
    i = i + 1
  }
  let n = node_ids.length()
  let adj : Array[Array[Int]] = Array::new()
  let mut j = 0
  while j < n {
    adj.push(Array::new())
    j = j + 1
  }
  let mut e = 0
  while e < edges.length() {
    let fi = sling_find_node_idx(node_ids, edges[e].from_id)
    let ti = sling_find_node_idx(node_ids, edges[e].to_id)
    if fi >= 0 && ti >= 0 {
      adj[fi].push(ti)
      adj[ti].push(fi)
    }
    e = e + 1
  }
  // BFS from start to end
  let start_idx = sling_find_node_idx(node_ids, start_id)
  let end_idx = sling_find_node_idx(node_ids, end_id)
  if start_idx < 0 || end_idx < 0 { return Array::new() }
  let visited : Array[Bool] = Array::new()
  let prev : Array[Int] = Array::new()
  let mut k = 0
  while k < n {
    visited.push(false)
    prev.push(-1)
    k = k + 1
  }
  let queue : Array[Int] = Array::new()
  queue.push(start_idx)
  visited[start_idx] = true
  let mut found = false
  while queue.length() > 0 && !found {
    let current = queue.remove(0)
    if current == end_idx { found = true; break }
    let mut ni = 0
    while ni < adj[current].length() {
      let next = adj[current][ni]
      if !visited[next] {
        visited[next] = true
        prev[next] = current
        queue.push(next)
      }
      ni = ni + 1
    }
  }
  if !found { return Array::new() }
  // Reconstruct path
  let path : Array[String] = Array::new()
  let mut cur = end_idx
  while cur >= 0 {
    path.push(node_ids[cur])
    if cur == start_idx { break }
    cur = prev[cur]
  }
  // Reverse
  let reversed : Array[String] = Array::new()
  let mut r = path.length() - 1
  while r >= 0 {
    reversed.push(path[r])
    r = r - 1
  }
  reversed
}

/// Fit principal curves for each lineage
fn sling_fit_lineage_curves(
  cell_coords : Array[Array[Double]],
  cluster_labels : Array[String],
  nodes : Array[SlingshotNode],
  lineages : Array[Array[String]]
) -> Array[SlingshotCurve] {
  let curves : Array[SlingshotCurve] = Array::new()
  let n_cells = cell_coords.length()
  let mut li = 0
  while li < lineages.length() {
    let lineage = lineages[li]
    // Collect cells belonging to clusters in this lineage
    let lineage_cells : Array[Int] = Array::new()
    let mut ci = 0
    while ci < n_cells {
      let cluster = sling_find_cluster_for_node(nodes, cluster_labels[ci])
      if cluster >= 0 && sling_lineage_contains(lineage, nodes[cluster].id) {
        lineage_cells.push(ci)
      }
      ci = ci + 1
    }
    // Extract coordinates for lineage cells
    let lineage_coords : Array[Array[Double]] = Array::new()
    let mut c = 0
    while c < lineage_cells.length() {
      lineage_coords.push(cell_coords[lineage_cells[c]].copy())
      c = c + 1
    }
    // Fit principal curve
    let n_segments = if lineage.length() >= 2 { (lineage.length() * 3) } else { 6 }
    let curve = sling_fit_principal_curve(lineage_coords, n_segments)
    // Compute pseudotime
    let pseudotime = sling_compute_pseudotime(lineage_coords, curve)
    let start_node = if lineage.length() > 0 { lineage[0] } else { "" }
    let end_node = if lineage.length() > 0 { lineage[lineage.length() - 1] } else { "" }
    curves.push(SlingshotCurve::{
      id: "curve_" + li.to_string(),
      start_node,
      end_node,
      control_points: curve,
      pseudotime,
      cell_indices: lineage_cells
    })
    li = li + 1
  }
  curves
}

fn sling_find_cluster_for_node(nodes : Array[SlingshotNode], cluster_label : String) -> Int {
  let mut i = 0
  while i < nodes.length() {
    if nodes[i].cluster == cluster_label { return i }
    i = i + 1
  }
  -1
}

fn sling_lineage_contains(lineage : Array[String], node_id : String) -> Bool {
  let mut i = 0
  while i < lineage.length() {
    if lineage[i] == node_id { return true }
    i = i + 1
  }
  false
}

/// Assign pseudotime to all cells
fn sling_assign_pseudotime(
  cell_coords : Array[Array[Double]],
  cluster_labels : Array[String],
  curves : Array[SlingshotCurve],
  nodes : Array[SlingshotNode]
) -> (Array[Double], Array[String], Array[String]) {
  let n_cells = cell_coords.length()
  let pseudotime : Array[Double] = Array::new()
  let lineages : Array[String] = Array::new()
  let branches : Array[String] = Array::new()
  let mut i = 0
  while i < n_cells {
    // Find which curve this cell belongs to
    let (ps, lin, br) = sling_assign_cell_to_curve(
      cell_coords[i], cluster_labels[i], curves, nodes
    )
    pseudotime.push(ps)
    lineages.push(lin)
    branches.push(br)
    i = i + 1
  }
  (pseudotime, lineages, branches)
}

/// Assign a cell to the nearest curve
fn sling_assign_cell_to_curve(
  cell_coord : Array[Double],
  cluster_label : String,
  curves : Array[SlingshotCurve],
  nodes : Array[SlingshotNode]
) -> (Double, String, String) {
  // Find curves containing this cell's cluster
  let cluster_idx = sling_find_cluster_for_node(nodes, cluster_label)
  if cluster_idx < 0 || curves.length() == 0 {
    return (0.0, "unassigned", "none")
  }
  let cluster_id = nodes[cluster_idx].id
  let mut best_ps = 0.0
  let mut best_curve = 0
  let mut best_dist = 1.0e30
  let mut ci = 0
  while ci < curves.length() {
    let curve = curves[ci]
    // Check if this curve's lineage contains the cluster
    let has_cluster = sling_curve_contains_cluster(curve, nodes, cluster_id)
    if has_cluster && curve.cell_indices.length() > 0 {
      // Find pseudotime for nearest cell on this curve
      let (ps, dist) = sling_pseudotime_for_cell(cell_coord, curve)
      if dist < best_dist {
        best_dist = dist
        best_ps = ps
        best_curve = ci
      }
    }
    ci = ci + 1
  }
  let lineage_name = "curve_" + best_curve.to_string()
  let branch_name = "branch_" + best_curve.to_string()
  (best_ps, lineage_name, branch_name)
}

fn sling_curve_contains_cluster(curve : SlingshotCurve, nodes : Array[SlingshotNode], cluster_id : String) -> Bool {
  // Check if start or end node of curve matches the cluster
  if curve.start_node == cluster_id || curve.end_node == cluster_id {
    return true
  }
  // Otherwise, check if any node in the lineage edges matches
  // Simplified: always include if curve has cells
  curve.cell_indices.length() > 0
}

/// Compute pseudotime for a cell on a given curve
fn sling_pseudotime_for_cell(cell_coord : Array[Double], curve : SlingshotCurve) -> (Double, Double) {
  let control = curve.control_points
  if control.length() < 2 {
    return (0.0, 1.0e30)
  }
  // Project cell onto curve
  let (t, dist) = sling_nearest_on_curve(cell_coord, control)
  (t, dist)
}

/// Get sorted pseudotime values
pub fn SlingshotResult::sling_sorted_pseudotime(self : SlingshotResult) -> Array[(Double, Int)] {
  let n = self.cell_pseudotime.length()
  let indexed : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((self.cell_pseudotime[i], i))
    i = i + 1
  }
  // Sort by pseudotime
  let sorted = indexed.copy()
  let mut a = 0
  while a < sorted.length() {
    let mut b = a + 1
    while b < sorted.length() {
      if sorted[b].0 < sorted[a].0 {
        let tmp = sorted[a]
        sorted[a] = sorted[b]
        sorted[b] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }
  sorted
}

/// Get cells along a specific lineage
pub fn SlingshotResult::sling_cells_on_lineage(self : SlingshotResult, lineage_id : String) -> Array[Int] {
  let cells : Array[Int] = Array::new()
  let n = self.cell_lineage.length()
  let mut i = 0
  while i < n {
    if self.cell_lineage[i] == lineage_id {
      cells.push(i)
    }
    i = i + 1
  }
  cells
}
