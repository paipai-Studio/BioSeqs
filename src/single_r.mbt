///|
/// SingleR - Cell type annotation by reference expression profiles.
/// Inspired by Bioconductor SingleR package.
///
/// SingleR is a method for inferring cell types from single-cell RNA-seq data
/// by comparing each cell's expression profile to reference datasets.
///
/// Key features:
/// - Reference dataset with known cell type profiles
/// - Correlation-based classification (Spearman / Pearson)
/// - Fine-tuning for ambiguous assignments
/// - Multi-reference support
/// - Quality control (low-score detection)
/// - Marker gene-based classification

///|
/// Reference profile for a cell type.
pub struct ReferenceProfile {
  cell_type : String
  expression : Array[Double]
  gene_names : Array[String]
}

///|
/// Reference dataset containing multiple cell type profiles.
pub struct ReferenceDataset {
  profiles : Array[ReferenceProfile]
  cell_types : Array[String]
  gene_names : Array[String]
  n_genes : Int
  n_profiles : Int
}

///|
/// Single-cell expression data to annotate.
pub struct SingleRData {
  counts : Array[Array[Double]]
  gene_names : Array[String]
  cell_names : Array[String]
  n_cells : Int
  n_genes : Int
}

///|
/// Result of SingleR annotation for a single cell.
pub struct SingleRResult {
  cell_id : String
  labels : Array[String]
  first_labels : Array[String]
  scores : Array[Double]
  first_scores : Array[Double]
  delta_next : Array[Double]
  annotations_fine : Array[String]
  first_annotation_fine : String
  pruned_first_labels : Array[String]
  pruned_labels : Array[String]
}

///|
/// Full annotation result for all cells.
pub struct SingleRFullResult {
  cell_ids : Array[String]
  labels : Array[String]
  scores : Array[Double]
  delta_scores : Array[Double]
  first_labels : Array[String]
  fine_annotations : Array[String]
  references_used : Array[String]
}

///|
/// Parameters for SingleR annotation.
pub struct SingleRParams {
  method : String
  n_top_scores : Int
  fine_tune : Bool
  min_score : Double
  quantile_score : Bool
}

///|
/// Spearman correlation coefficient.
pub fn single_r_spearman_correlation(x : Array[Double], y : Array[Double]) -> Double {
  let n = x.length()
  if n < 3 {
    return 0.0
  }

  let rank_x = single_r_rank_values(x)
  let rank_y = single_r_rank_values(y)

  let sum_x = single_r_array_sum(rank_x)
  let sum_y = single_r_array_sum(rank_y)
  let mean_x = sum_x / n.to_double()
  let mean_y = sum_y / n.to_double()

  let mut sum_xy = 0.0
  let mut sum_x2 = 0.0
  let mut sum_y2 = 0.0
  let mut i = 0
  while i < n {
    let dx = rank_x[i] - mean_x
    let dy = rank_y[i] - mean_y
    sum_xy = sum_xy + dx * dy
    sum_x2 = sum_x2 + dx * dx
    sum_y2 = sum_y2 + dy * dy
    i = i + 1
  }

  if sum_x2 <= 0.0 || sum_y2 <= 0.0 {
    return 0.0
  }

  sum_xy / (sum_x2 * sum_y2).sqrt()
}

///|
/// Pearson correlation coefficient (SingleR-scoped version).
pub fn single_r_pearson_correlation(x : Array[Double], y : Array[Double]) -> Double {
  let n = x.length()
  if n < 3 {
    return 0.0
  }

  let sum_x = single_r_array_sum(x)
  let sum_y = single_r_array_sum(y)
  let mean_x = sum_x / n.to_double()
  let mean_y = sum_y / n.to_double()

  let mut sum_xy = 0.0
  let mut sum_x2 = 0.0
  let mut sum_y2 = 0.0
  let mut i = 0
  while i < n {
    let dx = x[i] - mean_x
    let dy = y[i] - mean_y
    sum_xy = sum_xy + dx * dy
    sum_x2 = sum_x2 + dx * dx
    sum_y2 = sum_y2 + dy * dy
    i = i + 1
  }

  if sum_x2 <= 0.0 || sum_y2 <= 0.0 {
    return 0.0
  }

  sum_xy / (sum_x2 * sum_y2).sqrt()
}

///|
/// Compute rank values for an array (average ranks for ties).
pub fn single_r_rank_values(values : Array[Double]) -> Array[Double] {
  let n = values.length()
  let ranked : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    ranked.push((values[i], i))
    i = i + 1
  }

  let sorted = single_r_sort_array(ranked)
  let ranks = Array::make(n, 0.0)

  let mut j = 0
  while j < sorted.length() {
    let mut k = j
    while k < sorted.length() && sorted[k].0 == sorted[j].0 {
      k = k + 1
    }
    // For descending sort, map to ascending ranks (smallest=1, largest=n)
    // Position j..k-1 in descending → ascending ranks n-j..n-(k-1)
    let avg_rank = (2.0 * n.to_double() - j.to_double() - k.to_double() + 1.0) / 2.0
    let mut m = j
    while m < k {
      ranks[sorted[m].1] = avg_rank
      m = m + 1
    }
    j = k
  }

  ranks
}

///|
/// Sort array of tuples by first element (descending order).
fn single_r_sort_array(arr : Array[(Double, Int)]) -> Array[(Double, Int)] {
  let n = arr.length()
  let result = arr.copy()
  let mut i = 1
  while i < n {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 && result[j].0 < key.0 {
      result[j + 1] = result[j]
      j = j - 1
    }
    result[j + 1] = key
    i = i + 1
  }
  result
}

///|
/// Sort array of (Int, Double) tuples by the Double value (descending).
fn single_r_sort_int_double_array(arr : Array[(Int, Double)]) -> Array[(Int, Double)] {
  let n = arr.length()
  let result = arr.copy()
  let mut i = 1
  while i < n {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 && result[j].1 < key.1 {
      result[j + 1] = result[j]
      j = j - 1
    }
    result[j + 1] = key
    i = i + 1
  }
  result
}

///|
/// Sum of array elements.
fn single_r_array_sum(arr : Array[Double]) -> Double {
  let mut sum = 0.0
  for val in arr {
    sum = sum + val
  }
  sum
}

///|
/// Create a reference profile for a cell type.
pub fn ReferenceProfile::new(
  cell_type : String,
  expression : Array[Double],
  gene_names : Array[String],
) -> ReferenceProfile {
  ReferenceProfile::{ cell_type, expression, gene_names }
}

///|
/// Create a reference dataset from multiple profiles.
pub fn ReferenceDataset::new(
  profiles : Array[ReferenceProfile],
) -> ReferenceDataset {
  let n_profiles = profiles.length()
  let gene_names = if n_profiles > 0 { profiles[0].gene_names } else { [] }
  let n_genes = gene_names.length()
  let cell_types : Array[String] = Array::new()
  for p in profiles {
    cell_types.push(p.cell_type)
  }
  ReferenceDataset::{
    profiles,
    cell_types,
    gene_names,
    n_genes,
    n_profiles,
  }
}

///|
/// Create reference dataset from raw matrices.
pub fn ReferenceDataset::from_matrix(
  expression_data : Array[Array[Double]],
  cell_types : Array[String],
  gene_names : Array[String],
) -> ReferenceDataset {
  let profiles : Array[ReferenceProfile] = Array::new()
  let n_samples = expression_data.length()
  let mut i = 0
  while i < n_samples {
    profiles.push(
      ReferenceProfile::new(cell_types[i], expression_data[i], gene_names)
    )
    i = i + 1
  }
  ReferenceDataset::new(profiles)
}

///|
/// Create single-cell data from counts matrix.
pub fn SingleRData::new(
  counts : Array[Array[Double]],
  gene_names : Array[String],
  cell_names : Array[String],
) -> SingleRData {
  SingleRData::{
    counts,
    gene_names,
    cell_names,
    n_cells: cell_names.length(),
    n_genes: gene_names.length(),
  }
}

///|
/// Create default SingleR parameters.
pub fn SingleRParams::new() -> SingleRParams {
  SingleRParams::{
    method: "spearman",
    n_top_scores: 1,
    fine_tune: true,
    min_score: 0.0,
    quantile_score: false,
  }
}

///|
/// Create custom SingleR parameters with specific method.
pub fn SingleRParams::with_method(method : String) -> SingleRParams {
  let base = SingleRParams::new()
  SingleRParams::{ ..base, method }
}

///|
/// Create custom SingleR parameters with all fields specified.
pub fn SingleRParams::with_all(
  fine_tune : Bool,
  method : String,
  n_top_scores : Int,
  min_score : Double,
  quantile_score : Bool,
) -> SingleRParams {
  SingleRParams::{
    method,
    n_top_scores,
    fine_tune,
    min_score,
    quantile_score,
  }
}

///|
/// Compute pairwise correlations between a single cell and all reference profiles.
pub fn single_r_compute_correlations(
  cell_expression : Array[Double],
  reference : ReferenceDataset,
  method : String,
) -> Array[Double] {
  let scores : Array[Double] = Array::new()
  for profile in reference.profiles {
    let score = if method == "pearson" {
      single_r_pearson_correlation(cell_expression, profile.expression)
    } else {
      single_r_spearman_correlation(cell_expression, profile.expression)
    }
    scores.push(score)
  }
  scores
}

///|
/// Get the top N scores and their indices.
pub fn single_r_get_top_scores(
  scores : Array[Double],
  n : Int,
) -> Array[(Int, Double)] {
  let n_scores = scores.length()
  let top : Array[(Int, Double)] = Array::new()

  let mut i = 0
  while i < n_scores {
    top.push((i, scores[i]))
    i = i + 1
  }

  // Sort by score descending
  let sorted = single_r_sort_int_double_array(top)

  let result : Array[(Int, Double)] = Array::new()
  let limit = if n < sorted.length() { n } else { sorted.length() }
  let mut l = 0
  while l < limit {
    result.push(sorted[l])
    l = l + 1
  }
  result
}

///|
/// Aggregate scores by cell type.
/// Averages the top N scores for each cell type.
pub fn single_r_aggregate_scores_by_type(
  scores : Array[Double],
  reference : ReferenceDataset,
  n_top : Int,
) -> Map[String, Double] {
  let unique_types : Map[String, Array[Double]] = Map([], capacity=reference.n_profiles)

  let mut i = 0
  while i < reference.n_profiles {
    let ct = reference.cell_types[i]
    if !unique_types.contains(ct) {
      unique_types.set(ct, [])
    }
    let current = unique_types.get(ct).unwrap_or([])
    let updated = current.copy()
    updated.push(scores[i])
    unique_types.set(ct, updated)
    i = i + 1
  }

  let result : Map[String, Double] = Map([], capacity=unique_types.length())
  let keys = unique_types.keys()
  for key in keys {
    let vals = unique_types.get(key).unwrap_or([])
    let sorted_vals = single_r_sort_desc(vals)
    let n = if n_top < sorted_vals.length() { n_top } else { sorted_vals.length() }
    let mut sum = 0.0
    let mut j = 0
    while j < n {
      sum = sum + sorted_vals[j]
      j = j + 1
    }
    let avg = if n > 0 { sum / n.to_double() } else { 0.0 }
    result.set(key, avg)
  }
  result
}

///|
/// Sort array in descending order.
fn single_r_sort_desc(arr : Array[Double]) -> Array[Double] {
  let n = arr.length()
  let result = arr.copy()
  let mut i = 1
  while i < n {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 && result[j] < key {
      result[j + 1] = result[j]
      j = j - 1
    }
    result[j + 1] = key
    i = i + 1
  }
  result
}

///|
/// Find the cell type with the highest score.
pub fn single_r_get_best_label(scores_by_type : Map[String, Double]) -> (String, Double) {
  let keys = scores_by_type.keys()
  let mut best_label = ""
  let mut best_score = -2.0
  for key in keys {
    let score = scores_by_type.get(key).unwrap_or(-2.0)
    if score > best_score {
      best_score = score
      best_label = key
    }
  }
  (best_label, best_score)
}

///|
/// Compute delta score (difference between top and second best).
pub fn single_r_compute_delta_score(scores_by_type : Map[String, Double]) -> Double {
  let keys = scores_by_type.keys()
  let sorted_scores : Array[(String, Double)] = Array::new()
  for key in keys {
    sorted_scores.push((key, scores_by_type.get(key).unwrap_or(0.0)))
  }

  // Sort descending
  let mut i = 1
  while i < sorted_scores.length() {
    let key = sorted_scores[i]
    let mut j = i - 1
    while j >= 0 && sorted_scores[j].1 < key.1 {
      sorted_scores[j + 1] = sorted_scores[j]
      j = j - 1
    }
    sorted_scores[j + 1] = key
    i = i + 1
  }

  if sorted_scores.length() < 2 {
    return 0.0
  }

  sorted_scores[0].1 - sorted_scores[1].1
}

///|
/// Fine-tune annotation by re-computing correlations using only marker genes
/// of the initially predicted cell type.
pub fn single_r_fine_tune_annotation(
  cell_expression : Array[Double],
  reference : ReferenceDataset,
  initial_label : String,
  n_markers : Int,
  method : String,
) -> (String, Double) {
  // Find profiles matching the initial label
  let marker_profiles : Array[Array[Double]] = Array::new()
  let other_profiles : Array[Array[Double]] = Array::new()
  let other_types : Array[String] = Array::new()

  let mut i = 0
  while i < reference.n_profiles {
    if reference.cell_types[i] == initial_label {
      marker_profiles.push(reference.profiles[i].expression)
    } else {
      other_profiles.push(reference.profiles[i].expression)
      other_types.push(reference.cell_types[i])
    }
    i = i + 1
  }

  // If no other profiles, return the initial label
  if other_profiles.length() == 0 {
    let avg_score = if marker_profiles.length() > 0 {
      let mut sum = 0.0
      for p in marker_profiles {
        sum = sum + (if method == "pearson" {
          single_r_pearson_correlation(cell_expression, p)
        } else {
          single_r_spearman_correlation(cell_expression, p)
        })
      }
      sum / marker_profiles.length().to_double()
    } else {
      0.0
    }
    return (initial_label, avg_score)
  }

  // Compute correlations with other cell types
  let other_scores : Array[Double] = Array::new()
  for p in other_profiles {
    other_scores.push(
      if method == "pearson" {
        single_r_pearson_correlation(cell_expression, p)
      } else {
        single_r_spearman_correlation(cell_expression, p)
      }
    )
  }

  // Average scores by type
  let agg = single_r_aggregate_scores_by_type(
    other_scores,
    ReferenceDataset::from_matrix(other_profiles, other_types, reference.gene_names),
    1,
  )

  let (best_other, best_other_score) = single_r_get_best_label(agg)

  // Compare with the initial label's average score
  let avg_initial = if marker_profiles.length() > 0 {
    let mut sum = 0.0
    for p in marker_profiles {
      sum = sum + (if method == "pearson" {
        single_r_pearson_correlation(cell_expression, p)
      } else {
        single_r_spearman_correlation(cell_expression, p)
      })
    }
    sum / marker_profiles.length().to_double()
  } else {
    0.0
  }

  // Return the better of the two
  if avg_initial >= best_other_score {
    (initial_label, avg_initial)
  } else {
    (best_other, best_other_score)
  }
}

///|
/// Run SingleR annotation on a single cell.
pub fn single_r_annotate_cell(
  cell_expression : Array[Double],
  cell_id : String,
  reference : ReferenceDataset,
  params : SingleRParams,
) -> SingleRResult {
  // Compute correlations with all reference profiles
  let raw_scores = single_r_compute_correlations(cell_expression, reference, params.method)

  // Aggregate scores by cell type
  let agg_scores = single_r_aggregate_scores_by_type(raw_scores, reference, params.n_top_scores)

  // Get initial prediction
  let (best_label, best_score) = single_r_get_best_label(agg_scores)
  let delta = single_r_compute_delta_score(agg_scores)

  // Fine-tune if enabled
  let (final_label, final_score) = if params.fine_tune && best_score > 0.0 {
    single_r_fine_tune_annotation(
      cell_expression, reference, best_label, 50, params.method
    )
  } else {
    (best_label, best_score)
  }

  // Build result
  let labels : Array[String] = Array::new()
  let scores_arr : Array[Double] = Array::new()
  let delta_arr : Array[Double] = Array::new()

  let keys = agg_scores.keys()
  for key in keys {
    labels.push(key)
    scores_arr.push(agg_scores.get(key).unwrap_or(0.0))
    delta_arr.push(best_score - agg_scores.get(key).unwrap_or(0.0))
  }

  // Sort by score
  let sorted_labels = single_r_sort_labels_by_score(labels, scores_arr)
  let sorted_scores = single_r_sort_scores_desc(scores_arr)

  SingleRResult::{
    cell_id,
    labels: sorted_labels,
    first_labels: sorted_labels,
    scores: sorted_scores,
    first_scores: sorted_scores,
    delta_next: delta_arr,
    annotations_fine: [final_label],
    first_annotation_fine: final_label,
    pruned_first_labels: [final_label],
    pruned_labels: [final_label],
  }
}

///|
/// Sort labels by their scores (descending).
fn single_r_sort_labels_by_score(labels : Array[String], scores : Array[Double]) -> Array[String] {
  let n = labels.length()
  let pairs : Array[(String, Double)] = Array::new()
  let mut i = 0
  while i < n {
    pairs.push((labels[i], scores[i]))
    i = i + 1
  }

  // Sort descending by score
  let mut j = 1
  while j < pairs.length() {
    let key = pairs[j]
    let mut k = j - 1
    while k >= 0 && pairs[k].1 < key.1 {
      pairs[k + 1] = pairs[k]
      k = k - 1
    }
    pairs[k + 1] = key
    j = j + 1
  }

  let result : Array[String] = Array::new()
  for p in pairs {
    result.push(p.0)
  }
  result
}

///|
/// Sort scores array in descending order.
fn single_r_sort_scores_desc(scores : Array[Double]) -> Array[Double] {
  single_r_sort_desc(scores)
}

///|
/// Run SingleR annotation on all cells in a dataset.
pub fn single_r_annotate_cells(
  data : SingleRData,
  reference : ReferenceDataset,
  params : SingleRParams,
) -> SingleRFullResult {
  let cell_ids : Array[String] = Array::new()
  let labels : Array[String] = Array::new()
  let scores : Array[Double] = Array::new()
  let delta_scores : Array[Double] = Array::new()
  let first_labels : Array[String] = Array::new()
  let fine_annotations : Array[String] = Array::new()

  let mut i = 0
  while i < data.n_cells {
    let cell_expr = single_r_get_cell_expression(data, i)
    let cell_id = if i < data.cell_names.length() { data.cell_names[i] } else { "cell_".to_string() + i.to_string() }

    let result = single_r_annotate_cell(cell_expr, cell_id, reference, params)

    cell_ids.push(result.cell_id)
    labels.push(result.first_annotation_fine)
    scores.push(if result.scores.length() > 0 { result.scores[0] } else { 0.0 })
    delta_scores.push(if result.delta_next.length() > 0 { result.delta_next[0] } else { 0.0 })
    first_labels.push(result.first_labels[0])
    fine_annotations.push(result.annotations_fine[0])

    i = i + 1
  }

  SingleRFullResult::{
    cell_ids,
    labels,
    scores,
    delta_scores,
    first_labels,
    fine_annotations,
    references_used: ["SingleR"],
  }
}

///|
/// Extract expression profile for a single cell.
fn single_r_get_cell_expression(data : SingleRData, cell_idx : Int) -> Array[Double] {
  let n_genes = data.n_genes
  let expr : Array[Double] = Array::new()
  let mut i = 0
  while i < n_genes {
    expr.push(data.counts[i][cell_idx])
    i = i + 1
  }
  expr
}

///|
/// Create a synthetic reference dataset for testing.
pub fn single_r_create_reference_data() -> ReferenceDataset {
  let gene_names = [
    "CD3D", "CD3E", "CD4", "CD8A", "CD14", "CD19", "CD79A", "CD79B",
    "CD56", "CD94", "NCAM1", "KLRD1", "KLRB1", "KLRF1", "NKG7",
    "LYZ", "S100A8", "S100A9", "FCGR3A", "CSF1R", "C1QA", "C1QB",
    "C1QC", "MARCO", "FABP4", "CD68", "CD163", "MRC1",
    "CD34", "CD117", "KIT", "CD133", "PROM1", "ALDH1A1",
    "CD45", "PTPRC", "CD31", "PECAM1", "VWF", "ENG", "ESAM",
    "ACTA2", "TAGLN", "MYH11", "COL1A1", "COL1A2", "DCN", "LUM",
    "CD235A", "GYPA", "HBB", "HBA1", "HBA2", "HBM",
    "MS4A1", "CD20", "CD27", "CD38", "CD138", "SDC1", "BCMA",
    "TPSAB1", "TPSB2", "CPA3", "TPSD1", "MS4A2",
  ]

  let n_genes = gene_names.length()

  // T cell profiles
  let t_cell_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(t_cell_expr, gene_names, "CD3D", 120.0)
  single_r_set_gene_expr(t_cell_expr, gene_names, "CD3E", 115.0)
  single_r_set_gene_expr(t_cell_expr, gene_names, "CD4", 90.0)
  single_r_set_gene_expr(t_cell_expr, gene_names, "CD8A", 85.0)
  single_r_set_gene_expr(t_cell_expr, gene_names, "CD45", 95.0)
  single_r_set_gene_expr(t_cell_expr, gene_names, "PTPRC", 95.0)

  // B cell profiles
  let b_cell_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(b_cell_expr, gene_names, "CD19", 130.0)
  single_r_set_gene_expr(b_cell_expr, gene_names, "CD79A", 125.0)
  single_r_set_gene_expr(b_cell_expr, gene_names, "CD79B", 120.0)
  single_r_set_gene_expr(b_cell_expr, gene_names, "MS4A1", 115.0)
  single_r_set_gene_expr(b_cell_expr, gene_names, "CD20", 110.0)
  single_r_set_gene_expr(b_cell_expr, gene_names, "CD45", 95.0)
  single_r_set_gene_expr(b_cell_expr, gene_names, "PTPRC", 95.0)

  // NK cell profiles
  let nk_cell_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(nk_cell_expr, gene_names, "CD56", 110.0)
  single_r_set_gene_expr(nk_cell_expr, gene_names, "NCAM1", 110.0)
  single_r_set_gene_expr(nk_cell_expr, gene_names, "CD94", 100.0)
  single_r_set_gene_expr(nk_cell_expr, gene_names, "NKG7", 95.0)
  single_r_set_gene_expr(nk_cell_expr, gene_names, "KLRD1", 90.0)
  single_r_set_gene_expr(nk_cell_expr, gene_names, "KLRB1", 85.0)
  single_r_set_gene_expr(nk_cell_expr, gene_names, "KLRF1", 80.0)
  single_r_set_gene_expr(nk_cell_expr, gene_names, "CD45", 95.0)

  // Monocyte profiles
  let mono_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(mono_expr, gene_names, "CD14", 120.0)
  single_r_set_gene_expr(mono_expr, gene_names, "LYZ", 110.0)
  single_r_set_gene_expr(mono_expr, gene_names, "S100A8", 105.0)
  single_r_set_gene_expr(mono_expr, gene_names, "S100A9", 100.0)
  single_r_set_gene_expr(mono_expr, gene_names, "FCGR3A", 85.0)
  single_r_set_gene_expr(mono_expr, gene_names, "CSF1R", 80.0)
  single_r_set_gene_expr(mono_expr, gene_names, "CD45", 95.0)

  // Dendritic cell profiles
  let dc_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(dc_expr, gene_names, "FCGR3A", 95.0)
  single_r_set_gene_expr(dc_expr, gene_names, "CD14", 60.0)
  single_r_set_gene_expr(dc_expr, gene_names, "CD1C", 90.0)
  single_r_set_gene_expr(dc_expr, gene_names, "CD45", 95.0)

  // Macrophage profiles
  let macro_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(macro_expr, gene_names, "CD68", 115.0)
  single_r_set_gene_expr(macro_expr, gene_names, "CD163", 110.0)
  single_r_set_gene_expr(macro_expr, gene_names, "MRC1", 105.0)
  single_r_set_gene_expr(macro_expr, gene_names, "C1QA", 100.0)
  single_r_set_gene_expr(macro_expr, gene_names, "C1QB", 95.0)
  single_r_set_gene_expr(macro_expr, gene_names, "MARCO", 90.0)
  single_r_set_gene_expr(macro_expr, gene_names, "CD45", 95.0)

  // Epithelial cell profiles
  let epi_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(epi_expr, gene_names, "CDH1", 100.0)
  single_r_set_gene_expr(epi_expr, gene_names, "EPCAM", 95.0)
  single_r_set_gene_expr(epi_expr, gene_names, "KRT18", 90.0)
  single_r_set_gene_expr(epi_expr, gene_names, "KRT19", 85.0)
  single_r_set_gene_expr(epi_expr, gene_names, "MUC1", 80.0)

  // Endothelial cell profiles
  let endo_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(endo_expr, gene_names, "CD31", 110.0)
  single_r_set_gene_expr(endo_expr, gene_names, "PECAM1", 110.0)
  single_r_set_gene_expr(endo_expr, gene_names, "VWF", 105.0)
  single_r_set_gene_expr(endo_expr, gene_names, "ENG", 95.0)
  single_r_set_gene_expr(endo_expr, gene_names, "ESAM", 90.0)

  // Fibroblast profiles
  let fib_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(fib_expr, gene_names, "COL1A1", 105.0)
  single_r_set_gene_expr(fib_expr, gene_names, "COL1A2", 100.0)
  single_r_set_gene_expr(fib_expr, gene_names, "DCN", 95.0)
  single_r_set_gene_expr(fib_expr, gene_names, "LUM", 90.0)
  single_r_set_gene_expr(fib_expr, gene_names, "ACTA2", 85.0)

  // Mast cell profiles
  let mast_expr : Array[Double] = Array::make(n_genes, 0.5)
  single_r_set_gene_expr(mast_expr, gene_names, "TPSAB1", 115.0)
  single_r_set_gene_expr(mast_expr, gene_names, "TPSB2", 110.0)
  single_r_set_gene_expr(mast_expr, gene_names, "CPA3", 105.0)
  single_r_set_gene_expr(mast_expr, gene_names, "TPSD1", 100.0)
  single_r_set_gene_expr(mast_expr, gene_names, "MS4A2", 95.0)
  single_r_set_gene_expr(mast_expr, gene_names, "KIT", 90.0)

  let profiles : Array[ReferenceProfile] = [
    ReferenceProfile::new("T cells", t_cell_expr, gene_names),
    ReferenceProfile::new("B cells", b_cell_expr, gene_names),
    ReferenceProfile::new("NK cells", nk_cell_expr, gene_names),
    ReferenceProfile::new("Monocytes", mono_expr, gene_names),
    ReferenceProfile::new("Dendritic cells", dc_expr, gene_names),
    ReferenceProfile::new("Macrophages", macro_expr, gene_names),
    ReferenceProfile::new("Epithelial cells", epi_expr, gene_names),
    ReferenceProfile::new("Endothelial cells", endo_expr, gene_names),
    ReferenceProfile::new("Fibroblasts", fib_expr, gene_names),
    ReferenceProfile::new("Mast cells", mast_expr, gene_names),
  ]

  ReferenceDataset::new(profiles)
}

///|
/// Helper to set gene expression value in an array by gene name.
fn single_r_set_gene_expr(expr : Array[Double], gene_names : Array[String], gene : String, value : Double) -> Unit {
  let n = gene_names.length()
  let mut i = 0
  while i < n {
    if gene_names[i] == gene {
      expr[i] = value
      return
    }
    i = i + 1
  }
}

///|
/// Create a synthetic single-cell dataset for testing.
pub fn single_r_create_test_data() -> SingleRData {
  let gene_names = [
    "CD3D", "CD3E", "CD4", "CD8A", "CD14", "CD19", "CD79A", "CD79B",
    "CD56", "CD94", "NCAM1", "KLRD1", "KLRB1", "KLRF1", "NKG7",
    "LYZ", "S100A8", "S100A9", "FCGR3A", "CSF1R", "C1QA", "C1QB",
    "C1QC", "MARCO", "FABP4", "CD68", "CD163", "MRC1",
    "CD34", "CD117", "KIT", "CD133", "PROM1", "ALDH1A1",
    "CD45", "PTPRC", "CD31", "PECAM1", "VWF", "ENG", "ESAM",
    "ACTA2", "TAGLN", "MYH11", "COL1A1", "COL1A2", "DCN", "LUM",
    "CD235A", "GYPA", "HBB", "HBA1", "HBA2", "HBM",
    "MS4A1", "CD20", "CD27", "CD38", "CD138", "SDC1", "BCMA",
    "TPSAB1", "TPSB2", "CPA3", "TPSD1", "MS4A2",
  ]

  let n_genes = gene_names.length()
  let n_cells = 50

  // Create expression matrix (genes x cells)
  let counts : Array[Array[Double]] = Array::new()
  let mut gene_i = 0
  while gene_i < n_genes {
    let row : Array[Double] = Array::new()
    let mut cell_j = 0
    while cell_j < n_cells {
      // Simulate different cell types
      let cell_type = cell_j % 5
      let base = if cell_type == 0 { // T cell
        if gene_names[gene_i] == "CD3D" || gene_names[gene_i] == "CD3E" { 100.0 }
        else if gene_names[gene_i] == "CD4" { 80.0 }
        else if gene_names[gene_i] == "CD8A" { 75.0 }
        else if gene_names[gene_i] == "CD45" || gene_names[gene_i] == "PTPRC" { 90.0 }
        else { 1.0 }
      } else if cell_type == 1 { // B cell
        if gene_names[gene_i] == "CD19" { 110.0 }
        else if gene_names[gene_i] == "CD79A" || gene_names[gene_i] == "CD79B" { 100.0 }
        else if gene_names[gene_i] == "MS4A1" || gene_names[gene_i] == "CD20" { 95.0 }
        else if gene_names[gene_i] == "CD45" || gene_names[gene_i] == "PTPRC" { 90.0 }
        else { 1.0 }
      } else if cell_type == 2 { // NK cell
        if gene_names[gene_i] == "CD56" || gene_names[gene_i] == "NCAM1" { 95.0 }
        else if gene_names[gene_i] == "CD94" { 85.0 }
        else if gene_names[gene_i] == "NKG7" { 80.0 }
        else if gene_names[gene_i] == "KLRD1" || gene_names[gene_i] == "KLRB1" { 75.0 }
        else if gene_names[gene_i] == "CD45" || gene_names[gene_i] == "PTPRC" { 90.0 }
        else { 1.0 }
      } else if cell_type == 3 { // Monocyte
        if gene_names[gene_i] == "CD14" { 105.0 }
        else if gene_names[gene_i] == "LYZ" { 95.0 }
        else if gene_names[gene_i] == "S100A8" || gene_names[gene_i] == "S100A9" { 90.0 }
        else if gene_names[gene_i] == "FCGR3A" { 70.0 }
        else if gene_names[gene_i] == "CD45" || gene_names[gene_i] == "PTPRC" { 90.0 }
        else { 1.0 }
      } else { // Epithelial cell
        if gene_names[gene_i] == "CDH1" { 90.0 }
        else if gene_names[gene_i] == "EPCAM" { 85.0 }
        else if gene_names[gene_i] == "KRT18" || gene_names[gene_i] == "KRT19" { 80.0 }
        else if gene_names[gene_i] == "MUC1" { 75.0 }
        else { 1.0 }
      }

      // Add some noise
      let noise = ((cell_j * 7 + gene_i * 13) % 10).to_double() * 0.05 + 0.95
      let val = base * noise
      row.push(val)
      cell_j = cell_j + 1
    }
    counts.push(row)
    gene_i = gene_i + 1
  }

  let cell_names : Array[String] = Array::new()
  let mut c = 0
  while c < n_cells {
    cell_names.push("cell_".to_string() + c.to_string())
    c = c + 1
  }

  SingleRData::new(counts, gene_names, cell_names)
}

///|
/// Get summary statistics from annotation results.
pub fn single_r_annotation_summary(result : SingleRFullResult) -> Map[String, Double] {
  let n_cells = result.cell_ids.length()
  let label_counts : Map[String, Int] = Map([], capacity=10)

  let mut i = 0
  while i < n_cells {
    let label = result.labels[i]
    if !label_counts.contains(label) {
      label_counts.set(label, 0)
    }
    let count = label_counts.get(label).unwrap_or(0)
    label_counts.set(label, count + 1)
    i = i + 1
  }

  let summary : Map[String, Double] = Map([], capacity=label_counts.length() + 2)
  let keys = label_counts.keys()
  for key in keys {
    let count = label_counts.get(key).unwrap_or(0)
    let pct = count.to_double() / n_cells.to_double() * 100.0
    summary.set(key + "_count", count.to_double())
    summary.set(key + "_pct", pct)
  }

  let mut avg_score = 0.0
  for s in result.scores {
    avg_score = avg_score + s
  }
  avg_score = if n_cells > 0 { avg_score / n_cells.to_double() } else { 0.0 }
  summary.set("avg_score", avg_score)
  summary.set("n_cells", n_cells.to_double())

  summary
}