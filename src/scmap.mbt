///|
/// scmap - Single-cell cell type annotation by reference projection.
///
/// References:
///   - Bioconductor scmap package
///   - Kiselev VY, Yiu C, Hemberg M (2018) scmap: projection of single-cell
///     RNA-seq data across data sets. Nature Methods 15:359-362.
///
/// scmap assigns cell type labels to query cells by comparing their
/// expression profiles against a reference atlas of known cell types.
/// Two main methods are implemented:
///   1. scmap-cluster: compare query cells to cell type centroids
///      using Spearman correlation
///   2. scmap-cell: find nearest neighbours in the reference and
///      assign by majority vote
///
/// All public identifiers are prefixed with `scmap_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// Reference atlas: gene × cell expression matrix with cell type labels.
pub struct ScmapReference {
  gene_names : Array[String]
  cell_names : Array[String]
  cell_types : Array[String]
  expression : Array[Array[Double]] // rows=genes, cols=cells
}

///|
/// Create a reference atlas.
pub fn ScmapReference::new(
  gene_names : Array[String],
  cell_names : Array[String],
  cell_types : Array[String],
  expression : Array[Array[Double]],
) -> ScmapReference {
  { gene_names, cell_names, cell_types, expression }
}

///|
/// Get gene names.
pub fn ScmapReference::gene_names(self : ScmapReference) -> Array[String] {
  self.gene_names
}

///|
/// Get cell names.
pub fn ScmapReference::cell_names(self : ScmapReference) -> Array[String] {
  self.cell_names
}

///|
/// Get cell type labels.
pub fn ScmapReference::cell_types(self : ScmapReference) -> Array[String] {
  self.cell_types
}

///|
/// Get expression matrix (gene × cell).
pub fn ScmapReference::expression(self : ScmapReference) -> Array[Array[Double]] {
  self.expression
}

///|
/// Get unique cell types in the reference.
pub fn ScmapReference::unique_cell_types(self : ScmapReference) -> Array[String] {
  let types : Array[String] = []
  for t in self.cell_types {
    let mut found = false
    for existing in types {
      if existing == t {
        found = true
        break
      }
    }
    if !found {
      types.push(t)
    }
  }
  types
}

///|
/// Query data: gene × cell expression matrix for cells to be classified.
pub struct ScmapQuery {
  gene_names : Array[String]
  cell_names : Array[String]
  expression : Array[Array[Double]] // rows=genes, cols=cells
}

///|
/// Create query data.
pub fn ScmapQuery::new(
  gene_names : Array[String],
  cell_names : Array[String],
  expression : Array[Array[Double]],
) -> ScmapQuery {
  { gene_names, cell_names, expression }
}

///|
/// Get gene names.
pub fn ScmapQuery::gene_names(self : ScmapQuery) -> Array[String] {
  self.gene_names
}

///|
/// Get cell names.
pub fn ScmapQuery::cell_names(self : ScmapQuery) -> Array[String] {
  self.cell_names
}

///|
/// Get expression matrix.
pub fn ScmapQuery::expression(self : ScmapQuery) -> Array[Array[Double]] {
  self.expression
}

///|
/// Classification result for a single query cell.
pub struct ScmapAssignment {
  cell_name : String
  assigned_type : String // "unassigned" if below threshold
  best_correlation : Double
  second_best_correlation : Double
  method : String // "cluster" or "cell"
}

///|
/// Create an assignment.
pub fn ScmapAssignment::new(
  cell_name : String,
  assigned_type : String,
  best_correlation : Double,
  second_best_correlation : Double,
  method : String,
) -> ScmapAssignment {
  {
    cell_name,
    assigned_type,
    best_correlation,
    second_best_correlation,
    method,
  }
}

///|
/// Get cell name.
pub fn ScmapAssignment::cell_name(self : ScmapAssignment) -> String {
  self.cell_name
}

///|
/// Get assigned cell type.
pub fn ScmapAssignment::assigned_type(self : ScmapAssignment) -> String {
  self.assigned_type
}

///|
/// Get best correlation score.
pub fn ScmapAssignment::best_correlation(self : ScmapAssignment) -> Double {
  self.best_correlation
}

///|
/// Get second-best correlation score.
pub fn ScmapAssignment::second_best_correlation(self : ScmapAssignment) -> Double {
  self.second_best_correlation
}

///|
/// Get method name.
pub fn ScmapAssignment::method(self : ScmapAssignment) -> String {
  self.method
}

// ============================================================================
// Helper: extract a column from gene × cell matrix
// ============================================================================

///|
/// Extract a single cell's expression vector (over common genes).
fn scmap_extract_cell(
  expression : Array[Array[Double]],
  gene_names : Array[String],
  cell_idx : Int,
) -> Array[Double] {
  let n_genes = expression.length()
  let cell_expr : Array[Double] = Array::make(n_genes, 0.0)
  for g in 0..<n_genes {
    let row = expression[g]
    if cell_idx < row.length() {
      cell_expr[g] = row[cell_idx]
    }
  }
  cell_expr
}

///|
/// Compute the centroid (mean expression) for each cell type.
/// Returns: cell_types × genes centroid matrix (as Array of cell-type vectors).
fn scmap_centroids(
  ref : ScmapReference,
) -> Array[(String, Array[Double])] {
  let unique_types = ref.unique_cell_types()
  let n_genes = ref.gene_names.length()
  let centroids : Array[(String, Array[Double])] = []
  for ctype in unique_types {
    let sum : Array[Double] = Array::make(n_genes, 0.0)
    let mut count = 0
    for c in 0..<ref.cell_types.length() {
      if ref.cell_types[c] == ctype {
        for g in 0..<n_genes {
          sum[g] = sum[g] + ref.expression[g][c]
        }
        count = count + 1
      }
    }
    if count > 0 {
      let centroid : Array[Double] = Array::make(n_genes, 0.0)
      for g in 0..<n_genes {
        centroid[g] = sum[g] / count.to_double()
      }
      centroids.push((ctype, centroid))
    }
  }
  centroids
}

// ============================================================================
// scmap-cluster: centroid-based correlation
// ============================================================================

///|
/// Classify query cells using scmap-cluster method.
///
/// For each query cell, compute Spearman correlation with each reference
/// cell type centroid. Assign the cell type with highest correlation
/// if it exceeds the threshold; otherwise mark as "unassigned".
///
/// Parameters:
/// - reference: annotated reference atlas
/// - query: cells to classify
/// - threshold: minimum correlation to assign (default 0.7)
pub fn scmap_cluster(
  reference : ScmapReference,
  query : ScmapQuery,
  threshold : Double,
) -> Array[ScmapAssignment] {
  let centroids = scmap_centroids(reference)
  let n_query = query.cell_names.length()
  let assignments : Array[ScmapAssignment] = []
  for q in 0..<n_query {
    let query_vec = scmap_extract_cell(query.expression, query.gene_names, q)
    // Compute Spearman correlation with each centroid
    let mut best_type = "unassigned"
    let mut best_corr = -2.0
    let mut second_best = -2.0
    for entry in centroids {
      let (ctype, centroid) = entry
      let corr = stat_spearman_correlation(query_vec, centroid)
      if corr > best_corr {
        second_best = best_corr
        best_corr = corr
        best_type = ctype
      } else if corr > second_best {
        second_best = corr
      }
    }
    let assigned = if best_corr >= threshold { best_type } else { "unassigned" }
    assignments.push(
      ScmapAssignment::new(
        query.cell_names[q],
        assigned,
        best_corr,
        second_best,
        "cluster",
      ),
    )
  }
  assignments
}

// ============================================================================
// scmap-cell: nearest-neighbour voting
// ============================================================================

///|
/// Classify query cells using scmap-cell method.
///
/// For each query cell, find the k nearest reference cells by Spearman
/// correlation. If more than `threshold_fraction` of the k neighbours
/// belong to the same cell type, assign that type; otherwise "unassigned".
///
/// Parameters:
/// - reference: annotated reference atlas
/// - query: cells to classify
/// - k_neighbours: number of nearest neighbours (default 10)
/// - threshold_fraction: minimum fraction of agreeing neighbours (default 0.5)
pub fn scmap_cell(
  reference : ScmapReference,
  query : ScmapQuery,
  k_neighbours : Int,
  threshold_fraction : Double,
) -> Array[ScmapAssignment] {
  let n_query = query.cell_names.length()
  let n_ref = reference.cell_names.length()
  let assignments : Array[ScmapAssignment] = []
  for q in 0..<n_query {
    let query_vec = scmap_extract_cell(query.expression, query.gene_names, q)
    // Compute correlation with every reference cell
    let corrs : Array[(Double, Int)] = [] // (correlation, ref_cell_idx)
    for r in 0..<n_ref {
      let ref_vec = scmap_extract_cell(
        reference.expression,
        reference.gene_names,
        r,
      )
      let corr = stat_spearman_correlation(query_vec, ref_vec)
      corrs.push((corr, r))
    }
    // Sort by correlation descending (simple selection sort for top-k)
    let k = if k_neighbours < corrs.length() {
      k_neighbours
    } else {
      corrs.length()
    }
    // Find top-k by partial selection sort
    let top_k : Array[(Double, Int)] = []
    let used : Array[Bool] = Array::make(corrs.length(), false)
    for _ in 0..<k {
      let mut best_idx = -1
      let mut best_val = -2.0
      for i in 0..<corrs.length() {
        if !used[i] {
          let (corr, _) = corrs[i]
          if corr > best_val {
            best_val = corr
            best_idx = i
          }
        }
      }
      if best_idx >= 0 {
        used[best_idx] = true
        top_k.push(corrs[best_idx])
      }
    }
    // Vote among top-k neighbours
    let unique_types = reference.unique_cell_types()
    let votes : Array[Int] = Array::make(unique_types.length(), 0)
    let mut best_corr = -2.0
    let mut second_best = -2.0
    for i in 0..<top_k.length() {
      let (corr, ref_idx) = top_k[i]
      let ref_type = reference.cell_types[ref_idx]
      // Find type index
      for t in 0..<unique_types.length() {
        if unique_types[t] == ref_type {
          votes[t] = votes[t] + 1
          break
        }
      }
      if corr > best_corr {
        second_best = best_corr
        best_corr = corr
      } else if corr > second_best {
        second_best = corr
      }
    }
    // Find majority type
    let mut max_votes = 0
    let mut majority_type = "unassigned"
    for t in 0..<unique_types.length() {
      if votes[t] > max_votes {
        max_votes = votes[t]
        majority_type = unique_types[t]
      }
    }
    let fraction = max_votes.to_double() / k.to_double()
    let assigned = if fraction >= threshold_fraction {
      majority_type
    } else {
      "unassigned"
    }
    assignments.push(
      ScmapAssignment::new(
        query.cell_names[q],
        assigned,
        best_corr,
        second_best,
        "cell",
      ),
    )
  }
  assignments
}

// ============================================================================
// Summary utilities
// ============================================================================

///|
/// Summarise classification results: count cells per assigned type.
pub fn scmap_summary(assignments : Array[ScmapAssignment]) -> Array[(String, Int)] {
  let types : Array[String] = []
  let counts : Array[Int] = []
  for a in assignments {
    let mut found = false
    for i in 0..<types.length() {
      if types[i] == a.assigned_type {
        counts[i] = counts[i] + 1
        found = true
        break
      }
    }
    if !found {
      types.push(a.assigned_type)
      counts.push(1)
    }
  }
  let result : Array[(String, Int)] = []
  for i in 0..<types.length() {
    result.push((types[i], counts[i]))
  }
  result
}

///|
/// Format an assignment as a readable string.
pub fn ScmapAssignment::to_string(self : ScmapAssignment) -> String {
  self.cell_name + " -> " + self.assigned_type + " (corr=" + self.best_correlation.to_string() + ", method=" + self.method + ")"
}

///|
/// Format all assignments as a table.
pub fn scmap_assignments_to_string(
  assignments : Array[ScmapAssignment],
) -> String {
  let mut s = "scmap Assignments (" + assignments.length().to_string() + " cells)\n"
  s = s + "cell\tassigned_type\tbest_corr\tmethod\n"
  for a in assignments {
    s = s + a.cell_name + "\t" + a.assigned_type + "\t" + a.best_correlation.to_string() + "\t" + a.method + "\n"
  }
  s
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate a small reference atlas with 3 cell types and 3 cells each.
/// 10 genes, 9 reference cells.
pub fn scmap_sample_reference() -> ScmapReference {
  let gene_names = [
    "GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8", "GENE9", "GENE10",
  ]
  let cell_names = [
    "T1", "T2", "T3", "B1", "B2", "B3", "N1", "N2", "N3",
  ]
  let cell_types = [
    "T_cell", "T_cell", "T_cell", "B_cell", "B_cell", "B_cell", "NK_cell", "NK_cell", "NK_cell",
  ]
  // gene × cell expression matrix (10 genes × 9 cells)
  // T cells: high GENE1-3, low GENE4-6, GENE7-10 = 0
  // B cells: low GENE1-3, high GENE4-6, GENE7-10 = 0
  // NK cells: low GENE1-3, low GENE4-6, high GENE7-10
  let expression : Array[Array[Double]] = []
  // GENE1-3: T cell markers
  for g in 0..<3 {
    let row : Array[Double] = []
    for c in 0..<9 {
      let val = if c < 3 { 8.0 + (c.to_double() * 0.5) } else { 0.5 }
      row.push(val)
    }
    expression.push(row)
  }
  // GENE4-6: B cell markers
  for g in 0..<3 {
    let row : Array[Double] = []
    for c in 0..<9 {
      let val = if c >= 3 && c < 6 { 8.0 + (c.to_double() * 0.5) } else { 0.5 }
      row.push(val)
    }
    expression.push(row)
  }
  // GENE7-10: NK cell markers
  for g in 0..<4 {
    let row : Array[Double] = []
    for c in 0..<9 {
      let val = if c >= 6 { 8.0 + (c.to_double() * 0.5) } else { 0.5 }
      row.push(val)
    }
    expression.push(row)
  }
  ScmapReference::new(gene_names, cell_names, cell_types, expression)
}

///|
/// Generate 3 query cells: one clearly T, one clearly B, one ambiguous.
pub fn scmap_sample_query() -> ScmapQuery {
  let gene_names = [
    "GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8", "GENE9", "GENE10",
  ]
  let cell_names = ["query_T", "query_B", "query_unknown"]
  // gene × cell (10 genes × 3 query cells)
  let expression : Array[Array[Double]] = []
  // GENE1-3: high in query_T
  for g in 0..<3 {
    let row : Array[Double] = []
    row.push(7.5) // query_T
    row.push(0.6) // query_B
    row.push(1.0) // query_unknown
    expression.push(row)
  }
  // GENE4-6: high in query_B
  for g in 0..<3 {
    let row : Array[Double] = []
    row.push(0.6)
    row.push(7.5)
    row.push(1.0)
    expression.push(row)
  }
  // GENE7-10: low in all query cells
  for g in 0..<4 {
    let row : Array[Double] = []
    row.push(0.5)
    row.push(0.5)
    row.push(0.5)
    expression.push(row)
  }
  ScmapQuery::new(gene_names, cell_names, expression)
}
