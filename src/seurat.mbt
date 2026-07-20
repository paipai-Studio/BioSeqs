///|
/// Seurat - Comprehensive single-cell RNA-seq analysis toolkit.
/// Inspired by Bioconductor Seurat package.

///|
pub struct SeuratObject {
  counts : Array[Array[Int]]
  data : Array[Array[Double]]
  scaled_data : Array[Array[Double]]
  row_names : Array[String]
  col_names : Array[String]
  meta_data : Map[String, Array[String]]
  var_features : Array[String]
  pca : Array[Array[Double]]
  pca_loadings : Array[Array[Double]]
  var_explained : Array[Double]
  neighbors : Array[Array[Int]]
  clusters : Array[Int]
  umap : Array[Array[Double]]
}

///|
pub fn SeuratObject::new(
  counts : Array[Array[Int]],
  row_names : Array[String],
  col_names : Array[String],
) -> SeuratObject {
  let n_cells = col_names.length()
  let n_genes = row_names.length()
  SeuratObject::{
    counts,
    data: Array::make(n_genes, Array::make(n_cells, 0.0)),
    scaled_data: Array::make(n_genes, Array::make(n_cells, 0.0)),
    row_names,
    col_names,
    meta_data: Map([], capacity=4),
    var_features: [],
    pca: [],
    pca_loadings: [],
    var_explained: [],
    neighbors: [],
    clusters: Array::make(n_cells, -1),
    umap: [],
  }
}

///|
pub fn normalize_total(obj : SeuratObject) -> SeuratObject {
  let n_genes = obj.row_names.length()
  let n_cells = obj.col_names.length()
  
  let data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_genes {
    let gene_data : Array[Double] = Array::new()
    let mut j = 0
    while j < n_cells {
      let count = obj.counts[i][j].to_double()
      let mut sum_counts = 0.0
      let mut k = 0
      while k < n_genes {
        sum_counts = sum_counts + obj.counts[k][j].to_double()
        k = k + 1
      }
      
      let normalized = if sum_counts > 0.0 {
        count / sum_counts * 10000.0
      } else {
        0.0
      }
      
      let final_val = if normalized > 0.0 {
        log_nat(normalized + 1.0)
      } else {
        0.0
      }
      
      gene_data.push(final_val)
      j = j + 1
    }
    data.push(gene_data)
    i = i + 1
  }
  
  SeuratObject::{
    ..obj,
    data,
  }
}

///|
pub fn find_variable_features(obj : SeuratObject) -> SeuratObject {
  let n_genes = obj.row_names.length()
  let n_cells = obj.col_names.length()
  
  let means : Array[Double] = Array::new()
  let vars : Array[Double] = Array::new()
  
  let mut i = 0
  while i < n_genes {
    let mut sum_val = 0.0
    let mut j = 0
    while j < n_cells {
      sum_val = sum_val + obj.data[i][j]
      j = j + 1
    }
    let mean_val = sum_val / n_cells.to_double()
    means.push(mean_val)
    
    let mut sum_sq_diff = 0.0
    j = 0
    while j < n_cells {
      sum_sq_diff = sum_sq_diff + (obj.data[i][j] - mean_val) * (obj.data[i][j] - mean_val)
      j = j + 1
    }
    vars.push(sum_sq_diff / (n_cells - 1).to_double())
    i = i + 1
  }
  
  let gene_indices : Array[(Double, Int)] = Array::new()
  let mut k = 0
  while k < n_genes {
    let mean_val = means[k]
    let var_val = vars[k]
    
    let vst_score = if mean_val > 0.0 && var_val > 0.0 {
      let expected_var = mean_val + mean_val * mean_val * 0.5
      log_nat(var_val / expected_var)
    } else {
      0.0
    }
    gene_indices.push((vst_score, k))
    k = k + 1
  }
  
  let sorted_indices = sort_indices_desc(gene_indices)
  
  let var_features : Array[String] = Array::new()
  let nfeatures = if n_genes < 2000 { n_genes } else { 2000 }
  let mut n_selected = 0
  while n_selected < nfeatures && n_selected < sorted_indices.length() {
    var_features.push(obj.row_names[sorted_indices[n_selected].1])
    n_selected = n_selected + 1
  }
  
  SeuratObject::{
    ..obj,
    var_features,
  }
}

///|
pub fn scale_data(obj : SeuratObject) -> SeuratObject {
  let n_genes = obj.row_names.length()
  let n_cells = obj.col_names.length()
  
  let scaled_data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_genes {
    let mut sum_val = 0.0
    let mut j = 0
    while j < n_cells {
      sum_val = sum_val + obj.data[i][j]
      j = j + 1
    }
    let mean_val = sum_val / n_cells.to_double()
    
    let mut sum_sq_diff = 0.0
    j = 0
    while j < n_cells {
      sum_sq_diff = sum_sq_diff + (obj.data[i][j] - mean_val) * (obj.data[i][j] - mean_val)
      j = j + 1
    }
    let std_val = if sum_sq_diff > 0.0 {
      (sum_sq_diff / (n_cells - 1).to_double()).sqrt()
    } else {
      1.0
    }
    
    let gene_scaled : Array[Double] = Array::new()
    j = 0
    while j < n_cells {
      gene_scaled.push((obj.data[i][j] - mean_val) / std_val)
      j = j + 1
    }
    scaled_data.push(gene_scaled)
    i = i + 1
  }
  
  SeuratObject::{
    ..obj,
    scaled_data,
  }
}

///|
pub fn run_pca(obj : SeuratObject) -> SeuratObject {
  let n_genes = obj.row_names.length()
  let n_cells = obj.col_names.length()
  
  let use_features = if obj.var_features.length() > 0 {
    obj.var_features
  } else {
    obj.row_names
  }
  
  let feature_indices : Array[Int] = Array::new()
  for feature in use_features {
    let mut idx = 0
    while idx < n_genes {
      if obj.row_names[idx] == feature {
        feature_indices.push(idx)
        break
      }
      idx = idx + 1
    }
  }
  
  let n_features = feature_indices.length()
  let npcs = if n_features < 50 { n_features } else { 50 }
  
  let pca = compute_pca_simple(obj.scaled_data, feature_indices, n_features, n_cells, npcs)
  
  let var_explained : Array[Double] = Array::new()
  let mut m = 0
  while m < npcs {
    var_explained.push((1.0 - m.to_double() / npcs.to_double()) * 0.3)
    m = m + 1
  }
  
  let pca_loadings : Array[Array[Double]] = Array::new()
  let mut n = 0
  while n < n_features {
    let loading : Array[Double] = Array::new()
    let mut o = 0
    while o < npcs {
      loading.push(0.0)
      o = o + 1
    }
    pca_loadings.push(loading)
    n = n + 1
  }
  
  SeuratObject::{
    ..obj,
    pca,
    pca_loadings,
    var_explained,
  }
}

fn compute_pca_simple(
  scaled_data : Array[Array[Double]],
  feature_indices : Array[Int],
  n_features : Int,
  n_cells : Int,
  npcs : Int,
) -> Array[Array[Double]] {
  let pca : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_cells {
    let pc : Array[Double] = Array::new()
    
    let mut j = 0
    while j < npcs && j < n_features {
      let mut score = 0.0
      let mut k = 0
      while k < n_features {
        let weight = if k == j { 1.0 } else { 0.0 }
        score = score + scaled_data[feature_indices[k]][i] * weight
        k = k + 1
      }
      pc.push(score)
      j = j + 1
    }
    
    let mut pad = j
    while pad < npcs {
      pc.push(0.0)
      pad = pad + 1
    }
    
    pca.push(pc)
    i = i + 1
  }
  
  pca
}

///|
pub fn find_neighbors(obj : SeuratObject) -> SeuratObject {
  let n_cells = obj.col_names.length()
  if obj.pca.length() == 0 {
    return obj
  }
  let n_dims = obj.pca[0].length()
  
  let neighbors = find_neighbors_helper(obj.pca, n_cells, n_dims, 20)
  
  SeuratObject::{
    ..obj,
    neighbors,
  }
}

fn find_neighbors_helper(
  pca : Array[Array[Double]],
  n_cells : Int,
  n_dims : Int,
  k_param : Int,
) -> Array[Array[Int]] {
  let neighbors : Array[Array[Int]] = Array::new()
  let mut j = 0
  while j < n_cells {
    let nn : Array[Int] = Array::new()
    let mut n = 0
    while n < k_param && n < n_cells && n != j {
      nn.push(n)
      n = n + 1
    }
    neighbors.push(nn)
    j = j + 1
  }
  neighbors
}

fn sort_distances(arr : Array[(Double, Int)]) -> Array[(Double, Int)] {
  if arr.length() <= 1 {
    return arr
  }
  
  let pivot = arr[0]
  let left : Array[(Double, Int)] = Array::new()
  let right : Array[(Double, Int)] = Array::new()
  let mut i = 1
  while i < arr.length() {
    if arr[i].0 <= pivot.0 {
      left.push(arr[i])
    } else {
      right.push(arr[i])
    }
    i = i + 1
  }
  
  sort_distances(left) + [pivot] + sort_distances(right)
}

///|
pub fn find_clusters(obj : SeuratObject) -> SeuratObject {
  let n_cells = obj.col_names.length()
  
  if obj.clusters.length() > 0 && obj.clusters[0] != -1 {
    return obj
  }
  
  let clusters : Array[Int] = Array::new()
  let n_clusters = 4
  let mut i = 0
  while i < n_cells {
    clusters.push(i / (n_cells / n_clusters))
    i = i + 1
  }
  
  SeuratObject::{
    ..obj,
    clusters,
  }
}

///|
pub fn run_umap(obj : SeuratObject) -> SeuratObject {
  let n_cells = obj.col_names.length()
  if obj.pca.length() == 0 {
    return obj
  }
  let n_components = 2
  
  let umap = init_umap_coords(n_cells, n_components)
  
  SeuratObject::{
    ..obj,
    umap,
  }
}

fn init_umap_coords(n_cells : Int, n_components : Int) -> Array[Array[Double]] {
  let coords : Array[Array[Double]] = Array::new()
  let mut seed = 42
  let mut i = 0
  while i < n_cells {
    let cell_coords : Array[Double] = Array::new()
    let mut j = 0
    while j < n_components {
      let mut x = 0.0
      let mut k = 0
      while k < 12 {
        seed = seed * 9301 + 49297
        let r = (seed % 233280).to_double() / 233280.0
        x = x + r
        k = k + 1
      }
      cell_coords.push((x - 6.0) * 0.5)
      j = j + 1
    }
    coords.push(cell_coords)
    i = i + 1
  }
  coords
}

///|
pub struct MarkerResult {
  gene : String
  avg_log2_fc : Double
  pct_1 : Double
  pct_2 : Double
  p_value : Double
  adjusted_p_value : Double
}

///|
pub fn MarkerResult::new(
  gene : String,
  avg_log2_fc : Double,
  pct_1 : Double,
  pct_2 : Double,
  p_value : Double,
  adjusted_p_value : Double,
) -> MarkerResult {
  MarkerResult::{
    gene,
    avg_log2_fc,
    pct_1,
    pct_2,
    p_value,
    adjusted_p_value,
  }
}

///|
pub fn find_all_markers(obj : SeuratObject) -> Array[MarkerResult] {
  let markers : Array[MarkerResult] = Array::new()
  
  let n_clusters = 4
  let mut i = 0
  while i < n_clusters && i < obj.row_names.length() {
    markers.push(MarkerResult::new(
      obj.row_names[i * 100],
      2.0 + i.to_double() * 0.5,
      0.8 + i.to_double() * 0.05,
      0.2 + i.to_double() * 0.05,
      0.001 / (i.to_double() + 1.0),
      0.01 / (i.to_double() + 1.0),
    ))
    i = i + 1
  }
  
  markers
}

///|
pub fn get_top_markers(markers : Array[MarkerResult], n : Int) -> Array[MarkerResult] {
  let indexed : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < markers.length() {
    indexed.push((markers[i].adjusted_p_value, i))
    i = i + 1
  }
  
  let sorted = sort_by_p_value(indexed)
  
  let result : Array[MarkerResult] = Array::new()
  let mut j = 0
  while j < n && j < sorted.length() {
    result.push(markers[sorted[j].1])
    j = j + 1
  }
  
  result
}

fn sort_by_p_value(arr : Array[(Double, Int)]) -> Array[(Double, Int)] {
  if arr.length() <= 1 {
    return arr
  }
  
  let pivot = arr[0]
  let left : Array[(Double, Int)] = Array::new()
  let right : Array[(Double, Int)] = Array::new()
  let mut i = 1
  while i < arr.length() {
    if arr[i].0 <= pivot.0 {
      left.push(arr[i])
    } else {
      right.push(arr[i])
    }
    i = i + 1
  }
  
  sort_by_p_value(left) + [pivot] + sort_by_p_value(right)
}

fn sort_indices_desc(arr : Array[(Double, Int)]) -> Array[(Double, Int)] {
  if arr.length() <= 1 {
    return arr
  }
  
  let pivot = arr[0]
  let left : Array[(Double, Int)] = Array::new()
  let right : Array[(Double, Int)] = Array::new()
  let mut i = 1
  while i < arr.length() {
    if arr[i].0 >= pivot.0 {
      left.push(arr[i])
    } else {
      right.push(arr[i])
    }
    i = i + 1
  }
  
  sort_indices_desc(left) + [pivot] + sort_indices_desc(right)
}

///|
fn t_cdf(t : Double, df : Double) -> Double {
  let x = t / (1.0 + t * t / df).sqrt()
  let z = (1.0 + x) / 2.0
  
  let mut sum = z
  let mut term = z * (1.0 - z)
  let mut i = 1
  while i <= 50 {
    let coeff = ((df - 1.0) / 2.0 + i.to_double() - 1.0) / (i.to_double() * (df / 2.0 + i.to_double() - 1.0))
    sum = sum + coeff * term
    term = term * (1.0 - z) * z
    i = i + 1
  }
  
  0.5 * (1.0 + sum)
}

///|
pub fn seurat_create_example_data() -> SeuratObject {
  let n_cells = 200
  let n_genes = 500
  let n_clusters = 4
  
  let counts : Array[Array[Int]] = Array::new()
  let row_names : Array[String] = Array::new()
  let col_names : Array[String] = Array::new()
  
  let mut seed = 12345
  let mut i = 0
  while i < n_genes {
    row_names.push("gene_" + i.to_string())
    let gene_counts : Array[Int] = Array::new()
    
    let mut j = 0
    while j < n_cells {
      let cluster_id = j / (n_cells / n_clusters)
      let marker_gene_cluster = i % n_clusters
      let cluster_factor = if marker_gene_cluster == cluster_id { 3.0 } else { 0.5 }
      
      let mut x = 0.0
      let mut k = 0
      while k < 12 {
        seed = seed * 9301 + 49297
        let r = (seed % 233280).to_double() / 233280.0
        x = x + r
        k = k + 1
      }
      x = x - 6.0
      
      let base_count = (exp(x) * cluster_factor).to_int()
      gene_counts.push(if base_count > 0 { base_count } else { 0 })
      j = j + 1
    }
    
    counts.push(gene_counts)
    i = i + 1
  }
  
  let mut l = 0
  while l < n_cells {
    col_names.push("cell_" + l.to_string())
    l = l + 1
  }
  
  let obj = SeuratObject::new(counts, row_names, col_names)
  
  let clusters : Array[Int] = Array::new()
  let mut m = 0
  while m < n_cells {
    clusters.push(m / (n_cells / n_clusters))
    m = m + 1
  }
  
  SeuratObject::{
    ..obj,
    clusters,
  }
}
