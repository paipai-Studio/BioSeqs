///|
pub struct CellQCStats {
  cell_id : Array[String]
  total_counts : Array[Double]
  total_features : Array[Int]
  percent_mito : Array[Double]
  percent_ribo : Array[Double]
  percent_top50 : Array[Double]
  detected_genes : Array[Int]
}

pub fn CellQCStats::new(n_cells : Int) -> CellQCStats {
  CellQCStats::{
    cell_id: Array::make(n_cells, ""),
    total_counts: Array::make(n_cells, 0.0),
    total_features: Array::make(n_cells, 0),
    percent_mito: Array::make(n_cells, 0.0),
    percent_ribo: Array::make(n_cells, 0.0),
    percent_top50: Array::make(n_cells, 0.0),
    detected_genes: Array::make(n_cells, 0),
  }
}

pub fn CellQCStats::to_string(self : CellQCStats) -> String {
  let mut result = "CellQCStats\n"
  result = result + "===========\n"
  result = result + "Total cells: " + self.cell_id.length().to_string() + "\n"
  if self.cell_id.length() > 0 {
    result = result + "First cell: " + self.cell_id[0] + "\n"
    result = result + "  Total counts: " + self.total_counts[0].to_string() + "\n"
    result = result + "  Total features: " + self.total_features[0].to_string() + "\n"
    result = result + "  Percent mt: " + self.percent_mito[0].to_string() + "%\n"
  }
  result
}

///|
fn scater_is_mito_gene(gene_name : String) -> Bool {
  to_lowercase(gene_name).starts_with("mt-") || 
  to_lowercase(gene_name).starts_with("mt.") ||
  gene_name.starts_with("MT-") ||
  gene_name.starts_with("MT.")
}

///|
fn scater_is_ribo_gene(gene_name : String) -> Bool {
  let lower = to_lowercase(gene_name)
  lower.starts_with("rpl") || lower.starts_with("rps") ||
  lower.starts_with("mrpl") || lower.starts_with("mrps")
}

///|
pub fn scater_calculate_qc_metrics(
  counts : Array[Array[Int]],
  cell_ids : Array[String],
  gene_names : Array[String],
) -> CellQCStats {
  let n_cells = counts.length()
  let n_genes = if n_cells > 0 { counts[0].length() } else { 0 }
  
  let stats = CellQCStats::new(n_cells)
  
  let mut i = 0
  while i < n_cells {
    stats.cell_id[i] = cell_ids[i]
    
    let mut total = 0.0
    let mut n_detected = 0
    let mut mt_counts = 0.0
    let mut ribo_counts = 0.0
    
    let gene_counts : Array[Double] = Array::make(n_genes, 0.0)
    
    let mut j = 0
    while j < n_genes {
      let c = counts[i][j].to_double()
      total = total + c
      gene_counts[j] = c
      
      if counts[i][j] > 0 {
        n_detected = n_detected + 1
      }
      
      if scater_is_mito_gene(gene_names[j]) {
        mt_counts = mt_counts + c
      }
      
      if scater_is_ribo_gene(gene_names[j]) {
        ribo_counts = ribo_counts + c
      }
      
      j = j + 1
    }
    
    stats.total_counts[i] = total
    stats.total_features[i] = n_detected
    stats.detected_genes[i] = n_detected
    
    stats.percent_mito[i] = if total > 0.0 { mt_counts / total * 100.0 } else { 0.0 }
    stats.percent_ribo[i] = if total > 0.0 { ribo_counts / total * 100.0 } else { 0.0 }
    
    let sorted_counts = scater_sort_array_desc(gene_counts)
    let mut top50_sum = 0.0
    let top_n = if 50 < n_genes { 50 } else { n_genes }
    let mut k = 0
    while k < top_n {
      top50_sum = top50_sum + sorted_counts[k]
      k = k + 1
    }
    stats.percent_top50[i] = if total > 0.0 { top50_sum / total * 100.0 } else { 0.0 }
    
    i = i + 1
  }
  
  stats
}

///|
fn scater_sort_array_desc(arr : Array[Double]) -> Array[Double] {
  let result = arr.copy()
  let n = result.length()
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      if result[j] > result[i] {
        let temp = result[i]
        result[i] = result[j]
        result[j] = temp
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
pub struct QCFilterParams {
  min_counts : Double
  max_counts : Double
  min_features : Int
  max_features : Int
  max_percent_mito : Double
  min_detected : Int
}

pub fn QCFilterParams::new() -> QCFilterParams {
  QCFilterParams::{
    min_counts: 500.0,
    max_counts: 50000.0,
    min_features: 200,
    max_features: 10000,
    max_percent_mito: 20.0,
    min_detected: 200,
  }
}

pub fn QCFilterParams::with_params(
  min_counts : Double,
  max_counts : Double,
  min_features : Int,
  max_features : Int,
  max_percent_mito : Double,
  min_detected : Int,
) -> QCFilterParams {
  QCFilterParams::{
    min_counts,
    max_counts,
    min_features,
    max_features,
    max_percent_mito,
    min_detected,
  }
}

///|
pub fn scater_filter_cells(
  stats : CellQCStats,
  params : QCFilterParams,
) -> (Array[Int], Array[Bool]) {
  let n = stats.cell_id.length()
  let keep : Array[Bool] = Array::make(n, false)
  let kept_indices : Array[Int] = Array::new()
  
  let mut i = 0
  while i < n {
    let pass = 
      stats.total_counts[i] >= params.min_counts &&
      stats.total_counts[i] <= params.max_counts &&
      stats.total_features[i] >= params.min_features &&
      stats.total_features[i] <= params.max_features &&
      stats.percent_mito[i] <= params.max_percent_mito &&
      stats.detected_genes[i] >= params.min_detected
    
    if pass {
      keep[i] = true
      kept_indices.push(i)
    }
    
    i = i + 1
  }
  
  (kept_indices, keep)
}

///|
pub fn scater_filter_counts(
  counts : Array[Array[Int]],
  keep_cells : Array[Bool],
) -> Array[Array[Int]] {
  let filtered : Array[Array[Int]] = Array::new()
  
  let mut i = 0
  while i < counts.length() {
    if keep_cells[i] {
      filtered.push(counts[i])
    }
    i = i + 1
  }
  
  filtered
}

///|
pub fn scater_filter_cells_by_name(
  stats : CellQCStats,
  cell_ids_keep : Array[String],
) -> Array[Bool] {
  let n = stats.cell_id.length()
  let keep : Array[Bool] = Array::make(n, false)
  
  let keep_set : Map[String, Bool] = Map([], capacity=cell_ids_keep.length())
  let mut k = 0
  while k < cell_ids_keep.length() {
    keep_set[cell_ids_keep[k]] = true
    k = k + 1
  }
  
  let mut i = 0
  while i < n {
    if keep_set.contains(stats.cell_id[i]) {
      keep[i] = true
    }
    i = i + 1
  }
  
  keep
}

///|
pub struct FeatureQCStats {
  gene_id : Array[String]
  mean_expression : Array[Double]
  n_cells_detected : Array[Int]
  percent_detected : Array[Double]
  total_expression : Array[Double]
}

pub fn FeatureQCStats::new(n_genes : Int) -> FeatureQCStats {
  FeatureQCStats::{
    gene_id: Array::make(n_genes, ""),
    mean_expression: Array::make(n_genes, 0.0),
    n_cells_detected: Array::make(n_genes, 0),
    percent_detected: Array::make(n_genes, 0.0),
    total_expression: Array::make(n_genes, 0.0),
  }
}

///|
pub fn scater_calculate_feature_qc(
  counts : Array[Array[Int]],
  gene_names : Array[String],
) -> FeatureQCStats {
  let n_cells = counts.length()
  let n_genes = gene_names.length()
  
  let stats = FeatureQCStats::new(n_genes)
  
  let mut j = 0
  while j < n_genes {
    stats.gene_id[j] = gene_names[j]
    
    let mut total = 0.0
    let mut detected = 0
    
    let mut i = 0
    while i < n_cells {
      let c = counts[i][j].to_double()
      total = total + c
      if counts[i][j] > 0 {
        detected = detected + 1
      }
      i = i + 1
    }
    
    stats.total_expression[j] = total
    stats.mean_expression[j] = if n_cells > 0 { total / n_cells.to_double() } else { 0.0 }
    stats.n_cells_detected[j] = detected
    stats.percent_detected[j] = if n_cells > 0 { detected.to_double() / n_cells.to_double() * 100.0 } else { 0.0 }
    
    j = j + 1
  }
  
  stats
}

///|
pub fn scater_filter_features(
  stats : FeatureQCStats,
  min_cells_detected : Int,
) -> (Array[Int], Array[Bool]) {
  let n = stats.gene_id.length()
  let keep : Array[Bool] = Array::make(n, false)
  let kept_indices : Array[Int] = Array::new()
  
  let mut i = 0
  while i < n {
    if stats.n_cells_detected[i] >= min_cells_detected {
      keep[i] = true
      kept_indices.push(i)
    }
    i = i + 1
  }
  
  (kept_indices, keep)
}

///|
pub fn scater_normalize_cpm(counts : Array[Array[Int]]) -> Array[Array[Double]] {
  let n_cells = counts.length()
  let n_genes = if n_cells > 0 { counts[0].length() } else { 0 }
  
  let cpm : Array[Array[Double]] = Array::make(n_cells, Array::make(n_genes, 0.0))
  
  let mut i = 0
  while i < n_cells {
    let mut total = 0.0
    let mut j = 0
    while j < n_genes {
      total = total + counts[i][j].to_double()
      j = j + 1
    }
    
    let cell_cpm : Array[Double] = Array::make(n_genes, 0.0)
    let mut k = 0
    while k < n_genes {
      cell_cpm[k] = if total > 0.0 { counts[i][k].to_double() / total * 1000000.0 } else { 0.0 }
      k = k + 1
    }
    cpm[i] = cell_cpm
    
    i = i + 1
  }
  
  cpm
}

///|
pub fn scater_log_normalize_cpm(counts : Array[Array[Int]]) -> Array[Array[Double]] {
  let cpm = scater_normalize_cpm(counts)
  let n_cells = cpm.length()
  let n_genes = if n_cells > 0 { cpm[0].length() } else { 0 }
  
  let log_cpm : Array[Array[Double]] = Array::make(n_cells, Array::make(n_genes, 0.0))
  
  let mut i = 0
  while i < n_cells {
    let cell_log : Array[Double] = Array::make(n_genes, 0.0)
    let mut j = 0
    while j < n_genes {
      cell_log[j] = scater_log2(cpm[i][j] + 1.0)
      j = j + 1
    }
    log_cpm[i] = cell_log
    i = i + 1
  }
  
  log_cpm
}

///|
fn scater_log2(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    scater_ln(x) / scater_ln(2.0)
  }
}

///|
fn scater_ln(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    let mut result = 0.0
    let term = (x - 1.0) / (x + 1.0)
    let mut denom = 1.0
    let mut current = term
    let mut i = 0
    while i < 100 {
      result = result + current / denom
      current = current * term * term
      denom = denom + 2.0
      i = i + 1
    }
    2.0 * result
  }
}

///|
pub fn scater_find_hvgs(
  log_expr : Array[Array[Double]],
  gene_names : Array[String],
  n_hvgs : Int,
) -> Array[String] {
  let n_cells = log_expr.length()
  let n_genes = gene_names.length()
  
  let means : Array[Double] = Array::make(n_genes, 0.0)
  let vars : Array[Double] = Array::make(n_genes, 0.0)
  
  let mut j = 0
  while j < n_genes {
    let mut sum = 0.0
    let mut i = 0
    while i < n_cells {
      sum = sum + log_expr[i][j]
      i = i + 1
    }
    means[j] = sum / n_cells.to_double()
    
    let mut var_sum = 0.0
    let mut k = 0
    while k < n_cells {
      let diff = log_expr[k][j] - means[j]
      var_sum = var_sum + diff * diff
      k = k + 1
    }
    vars[j] = var_sum / (n_cells.to_double() - 1.0)
    
    j = j + 1
  }
  
  let gene_indices : Array[Int] = Array::make(n_genes, 0)
  let mut m = 0
  while m < n_genes {
    gene_indices[m] = m
    m = m + 1
  }
  
  let mut i = 0
  while i < n_genes {
    let mut best_idx = i
    let mut best_var = vars[i]
    let mut j = i + 1
    while j < n_genes {
      if vars[j] > best_var {
        best_var = vars[j]
        best_idx = j
      }
      j = j + 1
    }
    let temp = gene_indices[i]
    gene_indices[i] = gene_indices[best_idx]
    gene_indices[best_idx] = temp
    i = i + 1
  }
  
  let result : Array[String] = Array::new()
  let top_n = if n_hvgs < n_genes { n_hvgs } else { n_genes }
  let mut p = 0
  while p < top_n {
    result.push(gene_names[gene_indices[p]])
    p = p + 1
  }
  
  result
}

///|
pub fn scater_run_pca(
  log_expr : Array[Array[Double]],
  n_components : Int,
) -> Array[Array[Double]] {
  let n_cells = log_expr.length()
  let n_genes = if n_cells > 0 { log_expr[0].length() } else { 0 }
  
  let means : Array[Double] = Array::make(n_genes, 0.0)
  let mut j = 0
  while j < n_genes {
    let mut sum = 0.0
    let mut i = 0
    while i < n_cells {
      sum = sum + log_expr[i][j]
      i = i + 1
    }
    means[j] = sum / n_cells.to_double()
    j = j + 1
  }
  
  let centered : Array[Array[Double]] = Array::make(n_cells, Array::make(n_genes, 0.0))
  let mut k = 0
  while k < n_cells {
    let row : Array[Double] = Array::make(n_genes, 0.0)
    let mut m = 0
    while m < n_genes {
      row[m] = log_expr[k][m] - means[m]
      m = m + 1
    }
    centered[k] = row
    k = k + 1
  }
  
  let components : Array[Array[Double]] = Array::make(n_cells, Array::make(n_components, 0.0))
  
  let comp = scater_min_i(n_components, scater_min_i(n_cells, n_genes))
  let mut c = 0
  while c < comp {
    let mut v : Array[Double] = Array::make(n_genes, 0.0)
    let mut idx = c
    if idx >= n_genes {
      idx = n_genes - 1
    }
    if idx < 0 {
      idx = 0
    }
    v[idx] = 1.0
    
    let mut iter = 0
    while iter < 10 {
      let scores : Array[Double] = Array::make(n_cells, 0.0)
      let mut i = 0
      while i < n_cells {
        let mut dot = 0.0
        let mut j = 0
        while j < n_genes {
          dot = dot + centered[i][j] * v[j]
          j = j + 1
        }
        scores[i] = dot
        i = i + 1
      }
      
      let new_v : Array[Double] = Array::make(n_genes, 0.0)
      let mut norm = 0.0
      let mut j = 0
      while j < n_genes {
        let mut sum = 0.0
        let mut i = 0
        while i < n_cells {
          sum = sum + centered[i][j] * scores[i]
          i = i + 1
        }
        new_v[j] = sum
        norm = norm + sum * sum
        j = j + 1
      }
      
      if norm > 0.0 {
        let n = scater_sqrt_d(norm)
        let mut k2 = 0
        while k2 < n_genes {
          new_v[k2] = new_v[k2] / n
          k2 = k2 + 1
        }
      }
      
      v = new_v
      iter = iter + 1
    }
    
    let scores : Array[Double] = Array::make(n_cells, 0.0)
    let mut i = 0
    while i < n_cells {
      let mut dot = 0.0
      let mut j = 0
      while j < n_genes {
        dot = dot + centered[i][j] * v[j]
        j = j + 1
      }
      components[i][c] = dot
      i = i + 1
    }
    
    let mut i2 = 0
    while i2 < n_cells {
      let mut j2 = 0
      while j2 < n_genes {
        centered[i2][j2] = centered[i2][j2] - components[i2][c] * v[j2]
        j2 = j2 + 1
      }
      i2 = i2 + 1
    }
    
    c = c + 1
  }
  
  components
}

///|
fn scater_min_i(a : Int, b : Int) -> Int {
  if a < b { a } else { b }
}

///|
fn scater_sqrt_d(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    let mut z = x
    let mut prev = z
    let mut i = 0
    while i < 20 {
      z = (z + x / z) / 2.0
      if (z - prev).abs() < 0.0001 {
        break
      }
      prev = z
      i = i + 1
    }
    z
  }
}

///|
pub struct ScaterQCReport {
  cell_stats : CellQCStats
  feature_stats : FeatureQCStats
  n_cells_kept : Int
  n_features_kept : Int
  median_total_counts : Double
  median_n_features : Int
  median_percent_mito : Double
}

pub fn ScaterQCReport::new(
  cell_stats : CellQCStats,
  feature_stats : FeatureQCStats,
  params : QCFilterParams,
) -> ScaterQCReport {
  let (_, keep_cells) = scater_filter_cells(cell_stats, params)
  let (_, keep_features) = scater_filter_features(feature_stats, 10)
  
  let mut n_cells_kept = 0
  let mut i = 0
  while i < keep_cells.length() {
    if keep_cells[i] {
      n_cells_kept = n_cells_kept + 1
    }
    i = i + 1
  }
  
  let mut n_features_kept = 0
  let mut j = 0
  while j < keep_features.length() {
    if keep_features[j] {
      n_features_kept = n_features_kept + 1
    }
    j = j + 1
  }
  
  let sorted_counts = scater_sort_array_desc(cell_stats.total_counts)
  let mid = sorted_counts.length() / 2
  let median_counts = if sorted_counts.length() > 0 { sorted_counts[mid] } else { 0.0 }
  
  let features_arr : Array[Double] = Array::make(cell_stats.total_features.length(), 0.0)
  let mut k = 0
  while k < cell_stats.total_features.length() {
    features_arr[k] = cell_stats.total_features[k].to_double()
    k = k + 1
  }
  let sorted_features = scater_sort_array_desc(features_arr)
  let median_features = if sorted_features.length() > 0 { sorted_features[mid].to_int() } else { 0 }
  
  let sorted_mito = scater_sort_array_desc(cell_stats.percent_mito)
  let median_mito = if sorted_mito.length() > 0 { sorted_mito[mid] } else { 0.0 }
  
  ScaterQCReport::{
    cell_stats,
    feature_stats,
    n_cells_kept,
    n_features_kept,
    median_total_counts: median_counts,
    median_n_features: median_features,
    median_percent_mito: median_mito,
  }
}

pub fn ScaterQCReport::to_string(self : ScaterQCReport) -> String {
  let mut result = "Scater QC Report\n"
  result = result + "=================\n"
  result = result + "Total cells: " + self.cell_stats.cell_id.length().to_string() + "\n"
  result = result + "Cells kept: " + self.n_cells_kept.to_string() + "\n"
  result = result + "Total features: " + self.feature_stats.gene_id.length().to_string() + "\n"
  result = result + "Features kept: " + self.n_features_kept.to_string() + "\n\n"
  result = result + "Median total counts: " + self.median_total_counts.to_string() + "\n"
  result = result + "Median features: " + self.median_n_features.to_string() + "\n"
  result = result + "Median percent mt: " + self.median_percent_mito.to_string() + "%\n"
  result
}
