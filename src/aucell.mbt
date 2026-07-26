///|
/// AUCell - Gene set scoring for single-cell data.
/// Inspired by Bioconductor AUCell package.
/// Uses AUC (Area Under Curve) to score gene sets per cell.

///|
/// Helper: rank genes by expression for a single cell.
/// Returns indices sorted by expression value (ascending).
pub fn aucell_rank_genes(cell_expression : Array[Double]) -> Array[Int] {
  let n = cell_expression.length()
  let indices : Array[Int] = []
  let mut i = 0
  while i < n {
    indices.push(i)
    i = i + 1
  }
  // Bubble sort by expression value
  let arr = indices.copy()
  let mut swapped = true
  while swapped {
    swapped = false
    let mut j = 0
    while j < n - 1 {
      if cell_expression[arr[j]] > cell_expression[arr[j + 1]] {
        let tmp = arr[j]
        arr[j] = arr[j + 1]
        arr[j + 1] = tmp
        swapped = true
      }
      j = j + 1
    }
  }
  arr
}

///|
/// Helper: compute AUC (Area Under Curve) using the trapezoidal rule.
/// The curve is defined by points (x, y) where x is the rank and y is the
/// cumulative fraction of genes in the gene set.
pub fn aucell_compute_auc(
  total_genes : Int,
  gene_set_indices : Array[Int],
  ranks : Array[Int],
) -> Double {
  if total_genes == 0 {
    return 0.0
  }
  let n_genes = gene_set_indices.length()
  if n_genes == 0 {
    return 0.0
  }

  // Find the rank of each gene in the gene set
  let rank_lookup : Map[Int, Int] = Map([], capacity=total_genes)
  let mut i = 0
  while i < total_genes {
    rank_lookup.set(ranks[i], i)
    i = i + 1
  }

  // Collect ranks of genes in the gene set
  let gene_ranks : Array[Double] = []
  let mut g = 0
  while g < n_genes {
    let idx = gene_set_indices[g]
    if rank_lookup.contains(idx) {
      gene_ranks.push(rank_lookup[idx].to_double())
    }
    g = g + 1
  }

  if gene_ranks.length() == 0 {
    return 0.0
  }

  // Sort the ranks
  let sorted_ranks = gene_ranks.copy()
  let mut swapped = true
  let n = sorted_ranks.length()
  while swapped {
    swapped = false
    let mut j = 0
    while j < n - 1 {
      if sorted_ranks[j] > sorted_ranks[j + 1] {
        let tmp = sorted_ranks[j]
        sorted_ranks[j] = sorted_ranks[j + 1]
        sorted_ranks[j + 1] = tmp
        swapped = true
      }
      j = j + 1
    }
  }

  // Compute AUC using trapezoidal rule
  let max_rank = (total_genes - 1).to_double()
  let mut auc = 0.0
  let mut prev_x = 0.0
  let mut prev_y = 0.0
  let mut idx = 0
  while idx < n {
    let x = sorted_ranks[idx] / max_rank
    let y = (idx + 1).to_double() / n.to_double()
    // Trapezoid area
    auc = auc + (x - prev_x) * (prev_y + y) / 2.0
    prev_x = x
    prev_y = y
    idx = idx + 1
  }

  // Add area from last point to (1, 1)
  let last_x = if n > 0 { sorted_ranks[n - 1] / max_rank } else { 0.0 }
  let last_y = if n > 0 { n.to_double() / n.to_double() } else { 0.0 }
  auc = auc + (1.0 - last_x) * (last_y + 1.0) / 2.0

  auc
}

///|
pub struct AUCellData {
  expression : Array[Array[Double]]
  cell_names : Array[String]
  gene_names : Array[String]
  gene_sets : Array[GeneSet]
}

///|
pub fn AUCellData::new(
  expression : Array[Array[Double]],
  cell_names : Array[String],
  gene_names : Array[String],
  gene_sets : Array[GeneSet],
) -> AUCellData {
  AUCellData::{ expression, cell_names, gene_names, gene_sets }
}

///|
pub struct AUCellResult {
  cell_names : Array[String]
  gene_set_names : Array[String]
  score_matrix : Array[Array[Double]]
}

///|
pub fn AUCellResult::new(
  cell_names : Array[String],
  gene_set_names : Array[String],
  score_matrix : Array[Array[Double]],
) -> AUCellResult {
  AUCellResult::{ cell_names, gene_set_names, score_matrix }
}

///|
/// Map gene names to indices.
fn aucell_gene_index_map(gene_names : Array[String]) -> Map[String, Int] {
  let map : Map[String, Int] = Map([], capacity=gene_names.length())
  let mut i = 0
  while i < gene_names.length() {
    map.set(gene_names[i], i)
    i = i + 1
  }
  map
}

///|
/// Run AUCell scoring on single-cell expression data.
pub fn run_aucell(data : AUCellData) -> AUCellResult {
  let n_cells = data.cell_names.length()
  let n_genes = data.gene_names.length()
  let n_sets = data.gene_sets.length()

  let gene_index_map = aucell_gene_index_map(data.gene_names)

  // Build score matrix (cells × gene sets)
  let score_matrix : Array[Array[Double]] = []
  let mut c = 0
  while c < n_cells {
    let cell_expr = data.expression[c]
    let ranks = aucell_rank_genes(cell_expr)
    let scores : Array[Double] = []

    let mut s = 0
    while s < n_sets {
      let gene_set = data.gene_sets[s]
      // Map gene names to indices
      let gene_indices : Array[Int] = []
      let mut g = 0
      while g < gene_set.genes.length() {
        let gene_name = gene_set.genes[g]
        if gene_index_map.contains(gene_name) {
          gene_indices.push(gene_index_map[gene_name])
        }
        g = g + 1
      }

      // Compute AUC for this gene set
      let auc = if gene_indices.length() > 0 {
        aucell_compute_auc(n_genes, gene_indices, ranks)
      } else {
        0.0
      }
      scores.push(auc)
      s = s + 1
    }
    score_matrix.push(scores)
    c = c + 1
  }

  // Extract gene set names
  let gene_set_names : Array[String] = []
  let mut s = 0
  while s < n_sets {
    gene_set_names.push(data.gene_sets[s].name)
    s = s + 1
  }

  AUCellResult::{
    cell_names: data.cell_names.copy(),
    gene_set_names,
    score_matrix,
  }
}

///|
/// Normalize AUCell scores to [0, 1] range per gene set.
pub fn normalize_aucell_scores(result : AUCellResult) -> AUCellResult {
  let n_cells = result.cell_names.length()
  let n_sets = result.gene_set_names.length()

  let normalized_matrix : Array[Array[Double]] = []
  let mut c = 0
  while c < n_cells {
    let scores : Array[Double] = []
    let mut s = 0
    while s < n_sets {
      scores.push(result.score_matrix[c][s])
      s = s + 1
    }
    normalized_matrix.push(scores)
    c = c + 1
  }

  // Normalize each column
  let mut s = 0
  while s < n_sets {
    let mut min_val = normalized_matrix[0][s]
    let mut max_val = normalized_matrix[0][s]
    let mut c2 = 1
    while c2 < n_cells {
      if normalized_matrix[c2][s] < min_val {
        min_val = normalized_matrix[c2][s]
      }
      if normalized_matrix[c2][s] > max_val {
        max_val = normalized_matrix[c2][s]
      }
      c2 = c2 + 1
    }

    let range = max_val - min_val
    if range > 0.0 {
      let mut c3 = 0
      while c3 < n_cells {
        normalized_matrix[c3][s] = (normalized_matrix[c3][s] - min_val) / range
        c3 = c3 + 1
      }
    } else {
      let mut c3 = 0
      while c3 < n_cells {
        normalized_matrix[c3][s] = 0.5
        c3 = c3 + 1
      }
    }
    s = s + 1
  }

  AUCellResult::{
    cell_names: result.cell_names.copy(),
    gene_set_names: result.gene_set_names.copy(),
    score_matrix: normalized_matrix,
  }
}

///|
/// Get gene set scores for a specific cell.
pub fn get_cell_scores(
  result : AUCellResult,
  cell_name : String,
) -> Array[Double] {
  let mut idx = 0
  while idx < result.cell_names.length() {
    if result.cell_names[idx] == cell_name {
      return result.score_matrix[idx].copy()
    }
    idx = idx + 1
  }
  []
}

///|
/// Get cell scores for a specific gene set.
pub fn get_set_scores(
  result : AUCellResult,
  set_name : String,
) -> Array[Double] {
  let mut idx = 0
  while idx < result.gene_set_names.length() {
    if result.gene_set_names[idx] == set_name {
      let scores : Array[Double] = []
      let mut c = 0
      while c < result.score_matrix.length() {
        scores.push(result.score_matrix[c][idx])
        c = c + 1
      }
      return scores
    }
    idx = idx + 1
  }
  []
}
