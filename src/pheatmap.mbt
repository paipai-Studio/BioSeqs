///|
/// pheatmap - Pretty heatmaps inspired by Bioconductor pheatmap.
///
/// Provides:
/// - Distance computation (euclidean, manhattan, correlation)
/// - Hierarchical clustering (complete, average, ward)
/// - Heatmap rendering with color schemes
/// - Row/column annotation bars
/// - Cutree-based gaps
/// - Pretty break computation
/// - Font size and display options

///|
pub enum DistanceMethod {
  Euclidean
  Manhattan
  Correlation
  Maximum
}

///|
pub fn distance_method_euclidean() -> DistanceMethod { DistanceMethod::Euclidean }
pub fn distance_method_manhattan() -> DistanceMethod { DistanceMethod::Manhattan }
pub fn distance_method_correlation() -> DistanceMethod { DistanceMethod::Correlation }
pub fn distance_method_maximum() -> DistanceMethod { DistanceMethod::Maximum }

///|
pub enum ClusterMethod {
  Complete
  Average
  Ward
  Single
}

///|
pub fn cluster_method_complete() -> ClusterMethod { ClusterMethod::Complete }
pub fn cluster_method_average() -> ClusterMethod { ClusterMethod::Average }
pub fn cluster_method_ward() -> ClusterMethod { ClusterMethod::Ward }
pub fn cluster_method_single() -> ClusterMethod { ClusterMethod::Single }

///|
pub struct PheatmapInput {
  mat : Array[Array[Double]]
  row_names : Array[String]
  col_names : Array[String]
  color : Array[String]
  breaks : Array[Double]
  legend : Bool
  font_size : Double
  font_size_row : Double
  font_size_col : Double
  show_row_names : Bool
  show_col_names : Bool
  cluster_rows : Bool
  cluster_cols : Bool
  clustering_distance_rows : DistanceMethod
  clustering_distance_cols : DistanceMethod
  clustering_method : ClusterMethod
  cutree_rows : Int
  cutree_cols : Int
  annotation_row : Map[String, Array[String]]
  annotation_col : Map[String, Array[String]]
  annotation_colors : Map[String, Map[String, String]]
  display_numbers : Bool
  number_format : String
  number_color : String
  main_title : String
}

///|
pub struct PheatmapResult {
  text : String
  row_order : Array[Int]
  col_order : Array[Int]
  row_clusters : Array[Int]
  col_clusters : Array[Int]
  row_dendrogram : String
  col_dendrogram : String
  breaks : Array[Double]
  color : Array[String]
  legend : String
}

///|
pub fn distance_euclidean(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    let d = a[i] - b[i]
    sum = sum + d * d
    i = i + 1
  }
  sum.sqrt()
}

///|
pub fn distance_manhattan(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    let d = a[i] - b[i]
    if d < 0.0 { sum = sum - d } else { sum = sum + d }
    i = i + 1
  }
  sum
}

///|
pub fn distance_maximum(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  let mut mx = 0.0
  let mut i = 0
  while i < n {
    let d = a[i] - b[i]
    let ad = if d < 0.0 { -d } else { d }
    if ad > mx { mx = ad }
    i = i + 1
  }
  mx
}

///|
pub fn distance_correlation(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  if n < 2 { return 1.0 }
  let mut ma = 0.0
  let mut mb = 0.0
  let mut i = 0
  while i < n {
    ma = ma + a[i]
    mb = mb + b[i]
    i = i + 1
  }
  ma = ma / n.to_double()
  mb = mb / n.to_double()
  let mut num = 0.0
  let mut da = 0.0
  let mut db = 0.0
  let mut j = 0
  while j < n {
    let va = a[j] - ma
    let vb = b[j] - mb
    num = num + va * vb
    da = da + va * va
    db = db + vb * vb
    j = j + 1
  }
  let denom = (da * db).sqrt()
  if denom == 0.0 { return 1.0 }
  let corr = num / denom
  1.0 - corr
}

///|
pub fn pheatmap_distance_matrix(
  data : Array[Array[Double]],
  method : DistanceMethod,
) -> Array[Array[Double]] {
  let n = data.length()
  let dist : Array[Array[Double]] = Array::new()
  let mut ri = 0
  while ri < n { dist.push(Array::make(n, 0.0)); ri = ri + 1 }
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      if i == j {
        dist[i][j] = 0.0
      } else if i < j {
        let d = match method {
          DistanceMethod::Euclidean => distance_euclidean(data[i], data[j])
          DistanceMethod::Manhattan => distance_manhattan(data[i], data[j])
          DistanceMethod::Maximum => distance_maximum(data[i], data[j])
          DistanceMethod::Correlation => distance_correlation(data[i], data[j])
        }
        dist[i][j] = d
        dist[j][i] = d
      }
      j = j + 1
    }
    i = i + 1
  }
  dist
}

///|
pub fn pheatmap_column_distance_matrix(
  data : Array[Array[Double]],
  method : DistanceMethod,
) -> Array[Array[Double]] {
  let n = data.length()
  if n == 0 { return [] }
  let p = data[0].length()
  let col_data : Array[Array[Double]] = Array::new()
  let mut cj = 0
  while cj < p { col_data.push(Array::make(n, 0.0)); cj = cj + 1 }
  let mut j = 0
  while j < p {
    let mut i = 0
    while i < n {
      col_data[j][i] = data[i][j]
      i = i + 1
    }
    j = j + 1
  }
  pheatmap_distance_matrix(col_data, method)
}

///|
/// Simple hierarchical clustering using complete linkage.
pub struct HclustNode {
  left : Int
  right : Int
  height : Double
  size : Int
}

///|
pub struct HclustResult {
  merge : Array[(Int, Int)]
  height : Array[Double]
  order : Array[Int]
  n : Int
}

///|
pub fn hclust_complete(dist : Array[Array[Double]]) -> HclustResult {
  let n = dist.length()
  if n == 0 {
    return HclustResult::{ merge: [], height: [], order: [], n: 0 }
  }
  if n == 1 {
    return HclustResult::{ merge: [], height: [], order: [0], n: 1 }
  }
  // cluster[i] = current cluster id for point i (initially -i-1 for leaves)
  let cluster : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n { cluster[i] = -i - 1; i = i + 1 }
  // Working distance matrix for merging
  let work_dist : Array[Array[Double]] = Array::new()
  let mut wi = 0
  while wi < n { work_dist.push(Array::make(n, 0.0)); wi = wi + 1 }
  let mut ri = 0
  while ri < n {
    let mut cj = 0
    while cj < n { work_dist[ri][cj] = dist[ri][cj]; cj = cj + 1 }
    ri = ri + 1
  }
  let merge : Array[(Int, Int)] = Array::new()
  let height : Array[Double] = Array::new()
  // Process n-1 merges
  let mut step = 0
  while step < n - 1 {
    // Find minimum distance between different clusters
    let mut min_d = 1.0e18
    let mut min_i = 0
    let mut min_j = 0
    let mut a = 0
    while a < n {
      let mut b = a + 1
      while b < n {
        if cluster[a] != cluster[b] {
          let d = work_dist[a][b]
          if d < min_d { min_d = d; min_i = a; min_j = b }
        }
        b = b + 1
      }
      a = a + 1
    }
    // Merge clusters
    let new_id = step + 1
    merge.push((cluster[min_i], cluster[min_j]))
    height.push(min_d)
    // Update cluster assignments for both merged clusters
    let old_i = cluster[min_i]
    let old_j = cluster[min_j]
    let mut k = 0
    while k < n {
      if cluster[k] == old_i || cluster[k] == old_j { cluster[k] = new_id }
      k = k + 1
    }
    // Update distance matrix for the merged cluster (complete linkage: max distance)
    // Use min_i as the representative; update distances to all other points
    let mut c = 0
    while c < n {
      if cluster[c] != new_id {
        let new_d = if work_dist[min_i][c] > work_dist[min_j][c] {
          work_dist[min_i][c]
        } else {
          work_dist[min_j][c]
        }
        work_dist[min_i][c] = new_d
        work_dist[c][min_i] = new_d
      }
      c = c + 1
    }
    step = step + 1
  }
  // Compute order via inorder traversal of merge tree
  let order = hclust_inorder(merge, n)
  HclustResult::{ merge, height, order, n }
}

fn hclust_inorder(merge : Array[(Int, Int)], n : Int) -> Array[Int] {
  let result : Array[Int] = Array::new()
  if n <= 0 { return result }
  // The last merge (merge[length-1]) contains the final cluster
  let final_idx = merge.length() - 1
  let left_result = hclust_inorder_rec(merge, merge[final_idx].0)
  let right_result = hclust_inorder_rec(merge, merge[final_idx].1)
  // Combine
  let mut i = 0
  while i < left_result.length() { result.push(left_result[i]); i = i + 1 }
  let mut j = 0
  while j < right_result.length() { result.push(right_result[j]); j = j + 1 }
  result
}

fn hclust_inorder_rec(
  merge : Array[(Int, Int)],
  node : Int,
) -> Array[Int] {
  if node < 0 {
    // Leaf node: -node-1 gives the original index
    [-node - 1]
  } else {
    // node is a cluster id (1-indexed), map to merge array index (0-indexed)
    let merge_idx = node - 1
    let left = hclust_inorder_rec(merge, merge[merge_idx].0)
    let right = hclust_inorder_rec(merge, merge[merge_idx].1)
    let result : Array[Int] = Array::new()
    let mut i = 0
    while i < left.length() { result.push(left[i]); i = i + 1 }
    let mut j = 0
    while j < right.length() { result.push(right[j]); j = j + 1 }
    result
  }
}

///|
pub fn hclust_cutree(result : HclustResult, k : Int) -> Array[Int] {
  let n = result.n
  let order = result.order
  let clusters : Array[Int] = Array::make(n, 0)
  if k <= 1 || k >= n {
    let mut i = 0
    while i < n { clusters[order[i]] = 1; i = i + 1 }
    return clusters
  }
  // Find split points from merge tree
  let total_merges = result.merge.length()
  let cut_point = total_merges - k + 1
  // Create cluster assignments by following the merge tree
  let node_cluster : Map[Int, Int] = Map([], capacity = total_merges + n)
  let mut cluster_id = 0
  let mut i = 0
  while i < n {
    node_cluster[-i - 1] = cluster_id
    cluster_id = cluster_id + 1
    i = i + 1
  }
  let mut step = 0
  while step < total_merges && step < cut_point {
    let (left, right) = result.merge[step]
    let left_c = node_cluster[left]
    let right_c = node_cluster[right]
    // Find the minimum non-negative cluster id
    let new_c = if left_c < right_c { left_c } else { right_c }
    node_cluster[left] = new_c
    node_cluster[right] = new_c
    // Update all leaves that have the old cluster ids
    let mut j = 0
    while j < n {
      let leaf_id = -j - 1
      if node_cluster[leaf_id] == left_c || node_cluster[leaf_id] == right_c {
        node_cluster[leaf_id] = new_c
      }
      j = j + 1
    }
    step = step + 1
  }
  // Assign final cluster labels
  let label_map = Map([], capacity = n)
  let mut final_id = 0
  let mut k2 = 0
  while k2 < n {
    let leaf_id = -k2 - 1
    let c = node_cluster[leaf_id]
    if !label_map.contains(c) {
      label_map[c] = final_id
      final_id = final_id + 1
    }
    clusters[k2] = label_map[c] + 1
    k2 = k2 + 1
  }
  clusters
}

///|
/// Compute pretty breaks for color scale.
pub fn pretty_breaks(
  min_val : Double,
  max_val : Double,
  n : Int,
) -> Array[Double] {
  if min_val >= max_val {
    let breaks : Array[Double] = Array::make(1, min_val)
    return breaks
  }
  let range = max_val - min_val
  let step = range / (n.to_double() - 1.0)
  let breaks : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    breaks[i] = min_val + step * i.to_double()
    i = i + 1
  }
  breaks
}

///|
/// Generate a heat color palette (blue to red).
pub fn heat_colors(n : Int) -> Array[String] {
  let colors : Array[String] = Array::new()
  let mut i = 0
  while i < n {
    let t = i.to_double() / (n.to_double() - 1.0)
    let clamped_t = if t < 0.0 { 0.0 } else { if t > 1.0 { 1.0 } else { t } }
    let clamped_inv = if (1.0 - t) < 0.0 { 0.0 } else { if (1.0 - t) > 1.0 { 1.0 } else { 1.0 - t } }
    let r = (255.0 * clamped_t).to_int()
    let b = (255.0 * clamped_inv).to_int()
    let g = if t < 0.5 { (255.0 * t * 2.0).to_int() } else { (255.0 * (1.0 - t) * 2.0).to_int() }
    let clamp_r = if r < 0 { 0 } else if r > 255 { 255 } else { r }
    let clamp_g = if g < 0 { 0 } else if g > 255 { 255 } else { g }
    let clamp_b = if b < 0 { 0 } else if b > 255 { 255 } else { b }
    let hex = pheatmap_hex_color(clamp_r, clamp_g, clamp_b)
    colors.push(hex)
    i = i + 1
  }
  colors
}

fn pheatmap_hex_color(r : Int, g : Int, b : Int) -> String {
  let rh = pheatmap_int_to_hex(r)
  let gh = pheatmap_int_to_hex(g)
  let bh = pheatmap_int_to_hex(b)
  "#\{rh}\{gh}\{bh}"
}

fn pheatmap_int_to_hex(v : Int) -> String {
  let digits = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F"]
  let hi = v / 16
  let lo = v % 16
  "\{digits[hi]}\{digits[lo]}"
}

///|
/// Generate a blues color palette (white to blue).
pub fn blues_colors(n : Int) -> Array[String] {
  let colors : Array[String] = Array::new()
  let mut i = 0
  while i < n {
    let t = i.to_double() / (n.to_double() - 1.0)
    let r = (255.0 - 220.0 * t).to_int()
    let g = (255.0 - 220.0 * t).to_int()
    let b = 255
    let clamp_r = if r < 0 { 0 } else if r > 255 { 255 } else { r }
    let clamp_g = if g < 0 { 0 } else if g > 255 { 255 } else { g }
    colors.push(pheatmap_hex_color(clamp_r, clamp_g, b))
    i = i + 1
  }
  colors
}

///|
/// Generate a greys color palette (black to white).
pub fn greys_colors(n : Int) -> Array[String] {
  let colors : Array[String] = Array::new()
  let mut i = 0
  while i < n {
    let t = i.to_double() / (n.to_double() - 1.0)
    let v = (255.0 * t).to_int()
    let clamp_v = if v < 0 { 0 } else if v > 255 { 255 } else { v }
    colors.push(pheatmap_hex_color(clamp_v, clamp_v, clamp_v))
    i = i + 1
  }
  colors
}

///|
/// Create pheatmap input with defaults.
pub fn pheatmap_input(
  mat : Array[Array[Double]],
  row_names : Array[String],
  col_names : Array[String],
) -> PheatmapInput {
  let default_colors = heat_colors(100)
  let n = mat.length()
  let m = if n > 0 { mat[0].length() } else { 0 }
  let min_val = pheatmap_matrix_min(mat)
  let max_val = pheatmap_matrix_max(mat)
  let default_breaks = pretty_breaks(min_val, max_val, 101)
  PheatmapInput::{
    mat,
    row_names,
    col_names,
    color: default_colors,
    breaks: default_breaks,
    legend: true,
    font_size: 14.0,
    font_size_row: 14.0,
    font_size_col: 14.0,
    show_row_names: true,
    show_col_names: true,
    cluster_rows: true,
    cluster_cols: true,
    clustering_distance_rows: DistanceMethod::Euclidean,
    clustering_distance_cols: DistanceMethod::Euclidean,
    clustering_method: ClusterMethod::Complete,
    cutree_rows: 1,
    cutree_cols: 1,
    annotation_row: Map([], capacity = n),
    annotation_col: Map([], capacity = m),
    annotation_colors: Map([], capacity = 0),
    display_numbers: false,
    number_format: "{:.2f}",
    number_color: "black",
    main_title: "",
  }
}

fn pheatmap_matrix_min(mat : Array[Array[Double]]) -> Double {
  let n = mat.length()
  if n == 0 { return 0.0 }
  let mut mn = mat[0][0]
  let mut i = 0
  while i < n {
    let m = mat[i].length()
    let mut j = 0
    while j < m {
      if mat[i][j] < mn { mn = mat[i][j] }
      j = j + 1
    }
    i = i + 1
  }
  mn
}

fn pheatmap_matrix_max(mat : Array[Array[Double]]) -> Double {
  let n = mat.length()
  if n == 0 { return 0.0 }
  let mut mx = mat[0][0]
  let mut i = 0
  while i < n {
    let m = mat[i].length()
    let mut j = 0
    while j < m {
      if mat[i][j] > mx { mx = mat[i][j] }
      j = j + 1
    }
    i = i + 1
  }
  mx
}

///|
/// Set color scheme.
pub fn pheatmap_set_color(input : PheatmapInput, color : Array[String]) -> PheatmapInput {
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color,
    breaks: input.breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: input.annotation_col,
    annotation_colors: input.annotation_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Set breaks for color scale.
pub fn pheatmap_set_breaks(input : PheatmapInput, breaks : Array[Double]) -> PheatmapInput {
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: input.annotation_col,
    annotation_colors: input.annotation_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Set cutree_rows and cutree_cols for gaps.
pub fn pheatmap_set_cutree(
  input : PheatmapInput,
  cutree_rows : Int,
  cutree_cols : Int,
) -> PheatmapInput {
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks: input.breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows,
    cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: input.annotation_col,
    annotation_colors: input.annotation_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Set font sizes.
pub fn pheatmap_set_fontsize(
  input : PheatmapInput,
  font_size : Double,
) -> PheatmapInput {
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks: input.breaks,
    legend: input.legend,
    font_size,
    font_size_row: font_size,
    font_size_col: font_size,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: input.annotation_col,
    annotation_colors: input.annotation_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Set display options (show row/column names).
pub fn pheatmap_set_display(
  input : PheatmapInput,
  show_row_names : Bool,
  show_col_names : Bool,
) -> PheatmapInput {
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks: input.breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names,
    show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: input.annotation_col,
    annotation_colors: input.annotation_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Set clustering options (enable/disable row and column clustering).
pub fn pheatmap_set_cluster(
  input : PheatmapInput,
  cluster_rows : Bool,
  cluster_cols : Bool,
) -> PheatmapInput {
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks: input.breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows,
    cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: input.annotation_col,
    annotation_colors: input.annotation_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Set main title.
pub fn pheatmap_set_title(input : PheatmapInput, title : String) -> PheatmapInput {
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks: input.breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: input.annotation_col,
    annotation_colors: input.annotation_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: title,
  }
}

///|
/// Add row annotation.
pub fn pheatmap_add_annotation_row(
  input : PheatmapInput,
  name : String,
  values : Array[String],
) -> PheatmapInput {
  let new_annot = Map([], capacity = input.annotation_row.length() + 1)
  for key in input.annotation_row.keys() {
    new_annot[key] = input.annotation_row[key]
  }
  new_annot[name] = values
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks: input.breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: new_annot,
    annotation_col: input.annotation_col,
    annotation_colors: input.annotation_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Add column annotation.
pub fn pheatmap_add_annotation_col(
  input : PheatmapInput,
  name : String,
  values : Array[String],
) -> PheatmapInput {
  let new_annot = Map([], capacity = input.annotation_col.length() + 1)
  for key in input.annotation_col.keys() {
    new_annot[key] = input.annotation_col[key]
  }
  new_annot[name] = values
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks: input.breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: new_annot,
    annotation_colors: input.annotation_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Set annotation colors for a specific annotation.
pub fn pheatmap_set_annotation_color(
  input : PheatmapInput,
  annot_name : String,
  value_color_map : Map[String, String],
) -> PheatmapInput {
  let new_colors = Map([], capacity = input.annotation_colors.length() + 1)
  for key in input.annotation_colors.keys() {
    new_colors[key] = input.annotation_colors[key]
  }
  new_colors[annot_name] = value_color_map
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks: input.breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: input.annotation_col,
    annotation_colors: new_colors,
    display_numbers: input.display_numbers,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Add row annotation colors.
pub fn pheatmap_add_annotation_row_color(
  input : PheatmapInput,
  annot_name : String,
  value_color_map : Map[String, String],
) -> PheatmapInput {
  pheatmap_set_annotation_color(input, annot_name, value_color_map)
}

///|
/// Enable or disable display of numbers.
pub fn pheatmap_set_display_numbers(
  input : PheatmapInput,
  display : Bool,
) -> PheatmapInput {
  PheatmapInput::{
    mat: input.mat,
    row_names: input.row_names,
    col_names: input.col_names,
    color: input.color,
    breaks: input.breaks,
    legend: input.legend,
    font_size: input.font_size,
    font_size_row: input.font_size_row,
    font_size_col: input.font_size_col,
    show_row_names: input.show_row_names,
    show_col_names: input.show_col_names,
    cluster_rows: input.cluster_rows,
    cluster_cols: input.cluster_cols,
    clustering_distance_rows: input.clustering_distance_rows,
    clustering_distance_cols: input.clustering_distance_cols,
    clustering_method: input.clustering_method,
    cutree_rows: input.cutree_rows,
    cutree_cols: input.cutree_cols,
    annotation_row: input.annotation_row,
    annotation_col: input.annotation_col,
    annotation_colors: input.annotation_colors,
    display_numbers: display,
    number_format: input.number_format,
    number_color: input.number_color,
    main_title: input.main_title,
  }
}

///|
/// Render the heatmap.
pub fn pheatmap_render(input : PheatmapInput) -> PheatmapResult {
  let mat = input.mat
  let n = mat.length()
  let m = if n > 0 { mat[0].length() } else { 0 }
  if n == 0 || m == 0 {
    return PheatmapResult::{
      text: "",
      row_order: [],
      col_order: [],
      row_clusters: [],
      col_clusters: [],
      row_dendrogram: "",
      col_dendrogram: "",
      breaks: input.breaks,
      color: input.color,
      legend: "",
    }
  }
  // Compute row and column order via clustering
  let mut row_order : Array[Int] = Array::make(n, 0)
  let mut col_order : Array[Int] = Array::make(m, 0)
  let mut i = 0
  while i < n { row_order[i] = i; i = i + 1 }
  let mut j = 0
  while j < m { col_order[j] = j; j = j + 1 }
  let mut row_clusters : Array[Int] = Array::new()
  let mut col_clusters : Array[Int] = Array::new()
  let mut row_dendrogram = ""
  let mut col_dendrogram = ""
  if input.cluster_rows && n > 1 {
    let row_dist = pheatmap_distance_matrix(mat, input.clustering_distance_rows)
    let row_hc = hclust_complete(row_dist)
    row_order = row_hc.order
    row_dendrogram = "clustered ({row_hc.merge.length()} merges)"
    if input.cutree_rows > 1 && input.cutree_rows <= n {
      row_clusters = hclust_cutree(row_hc, input.cutree_rows)
    }
  }
  if input.cluster_cols && m > 1 {
    let col_dist = pheatmap_column_distance_matrix(mat, input.clustering_distance_cols)
    let col_hc = hclust_complete(col_dist)
    col_order = col_hc.order
    col_dendrogram = "clustered ({col_hc.merge.length()} merges)"
    if input.cutree_cols > 1 && input.cutree_cols <= m {
      col_clusters = hclust_cutree(col_hc, input.cutree_cols)
    }
  }
  // Build text output
  let text = pheatmap_build_text(input, row_order, col_order, row_clusters, col_clusters)
  let legend = pheatmap_build_legend(input.breaks, input.color, input.legend)
  PheatmapResult::{
    text,
    row_order,
    col_order,
    row_clusters,
    col_clusters,
    row_dendrogram,
    col_dendrogram,
    breaks: input.breaks,
    color: input.color,
    legend,
  }
}

fn pheatmap_build_text(
  input : PheatmapInput,
  row_order : Array[Int],
  col_order : Array[Int],
  row_clusters : Array[Int],
  col_clusters : Array[Int],
) -> String {
  let lines : Array[String] = Array::new()
  if input.main_title.length() > 0 {
    lines.push(input.main_title)
    lines.push("")
  }
  let n = row_order.length()
  let m = col_order.length()
  // Header row
  if input.show_col_names {
    let header = "      " + pheatmap_pad_str(" ", 4)
    let mut j = 0
    while j < m {
      let name = input.col_names[col_order[j]]
      let padded = pheatmap_pad_str(name, 8)
      lines.push("      \{padded}")
      j = j + 1
    }
    // Actually, let's build it as a single line
    let header_parts : Array[String] = Array::new()
    let _ = header
    let mut header_line = "      "
    let mut k = 0
    while k < m {
      let name = input.col_names[col_order[k]]
      header_line = header_line + pheatmap_pad_str(name, 8)
      k = k + 1
    }
    lines.push(header_line)
  }
  // Data rows
  let mut i = 0
  while i < n {
    let row_idx = row_order[i]
    let row_name = if input.show_row_names { input.row_names[row_idx] } else { "" }
    let mut line = pheatmap_pad_str(row_name, 8)
    let mut j = 0
    while j < m {
      let col_idx = col_order[j]
      let val = input.mat[row_idx][col_idx]
      let color_idx = pheatmap_value_to_color_idx(val, input.breaks, input.color.length())
      let color = input.color[color_idx]
      let cell = if input.display_numbers {
        pheatmap_format_value(val, input.number_format)
      } else {
        pheatmap_color_cell(color)
      }
      line = line + cell
      j = j + 1
    }
    lines.push(line)
    i = i + 1
  }
  // Cluster info
  if row_clusters.length() > 0 {
    lines.push("")
    lines.push("Row clusters: \{row_clusters.length()} points assigned to \{pheatmap_unique_count(row_clusters)} clusters")
  }
  if col_clusters.length() > 0 {
    lines.push("Column clusters: \{col_clusters.length()} points assigned to \{pheatmap_unique_count(col_clusters)} clusters")
  }
  lines.join("\n")
}

fn pheatmap_pad_str(s : String, width : Int) -> String {
  if s.length() >= width { return s }
  let padded = s
  let mut result = padded
  let mut i = padded.length()
  while i < width { result = result + " "; i = i + 1 }
  result
}

fn pheatmap_value_to_color_idx(
  val : Double,
  breaks : Array[Double],
  n_colors : Int,
) -> Int {
  if n_colors == 0 { return 0 }
  if breaks.length() < 2 { return n_colors / 2 }
  let min_b = breaks[0]
  let max_b = breaks[breaks.length() - 1]
  if val <= min_b { return 0 }
  if val >= max_b { return n_colors - 1 }
  let range = max_b - min_b
  let idx = ((val - min_b) / range * (n_colors - 1).to_double()).to_int()
  if idx < 0 { return 0 }
  if idx >= n_colors { return n_colors - 1 }
  idx
}

fn pheatmap_color_cell(color : String) -> String {
  // Simple text representation using color name
  "[" + color + "]"
}

fn pheatmap_format_value(val : Double, format : String) -> String {
  // Simple format: show 2 decimal places
  let s = val.to_string()
  // Truncate to 5 chars
  if s.length() > 6 { return s[0:6].to_owned() }
  s
}

fn pheatmap_unique_count(arr : Array[Int]) -> Int {
  let seen = Map([], capacity = arr.length())
  let mut count = 0
  for v in arr {
    if !seen.contains(v) { seen[v] = true; count = count + 1 }
  }
  count
}

fn pheatmap_build_legend(
  breaks : Array[Double],
  color : Array[String],
  show : Bool,
) -> String {
  if !show { return "" }
  let n = color.length()
  if n == 0 { return "" }
  let lines : Array[String] = Array::new()
  lines.push("Legend:")
  lines.push("  Min: \{breaks[0]}")
  let mid_idx = n / 2
  let mid_val = if breaks.length() > mid_idx { breaks[mid_idx] } else { breaks[breaks.length() - 1] }
  lines.push("  Mid: \{mid_val}")
  lines.push("  Max: \{breaks[breaks.length() - 1]}")
  lines.join("\n")
}
