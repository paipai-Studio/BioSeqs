///|
/// SCENIC module - Single-Cell rEgulatory Network Inference and Clustering.
/// Reference: Aibar et al. (2017) SCENIC: single-cell regulatory network
/// inference and clustering. Nature Methods 14:1083-1086.
///
/// SCENIC is a three-step pipeline:
///   1. GRN inference: TF-target co-expression modules (GENIE3/GRNBoost2 style)
///   2. Regulon construction: prune modules to direct targets (cisTarget-style)
///   3. AUCell scoring: per-cell regulon activity + binary thresholding
///
/// This implementation provides a self-contained pipeline that can run without
/// external motif databases by using importance-weight pruning as a proxy for
/// motif enrichment. When a motif ranking database is available, the
/// `prune_by_motif_ranking` function can be used for faithful cisTarget-style
/// pruning.

// ============================================================================
// Core data structures
// ============================================================================

///|
/// A TF-target co-expression module from Step 1.
pub struct CoExpressionModule {
  tf_name : String
  targets : Array[String]
  weights : Array[Double]    // importance weights, parallel to targets
} derive(Debug)

///|
/// Create a co-expression module.
pub fn co_expression_module(
  tf_name : String,
  targets : Array[String],
  weights : Array[Double],
) -> CoExpressionModule {
  { tf_name: tf_name, targets: targets, weights: weights }
}

///|
/// A regulon: TF + pruned direct target set with weights (Step 2 output).
pub struct ScenicRegulon {
  tf_name : String
  targets : Array[String]
  weights : Array[Double]    // weight per target (from GRN importance)
  n_targets : Int
} derive(Debug)

///|
/// Create a regulon from a TF name and parallel target/weight arrays.
pub fn scenic_regulon(
  tf_name : String,
  targets : Array[String],
  weights : Array[Double],
) -> ScenicRegulon {
  {
    tf_name: tf_name,
    targets: targets,
    weights: weights,
    n_targets: targets.length(),
  }
}

///|
/// AUC matrix entry: regulon activity score for one cell.
pub struct AUCEntry {
  regulon_name : String
  cell_name : String
  auc : Double               // normalized AUC in [0, 1]
}

///|
/// SCENIC input: expression matrix + gene/cell names + TF list.
pub struct SCENICInput {
  n_genes : Int
  n_cells : Int
  expression : Array[Array[Double]]   // gene x cell (row=gene, col=cell)
  gene_names : Array[String]
  cell_names : Array[String]
  tf_names : Array[String]            // subset of gene_names that are TFs
} derive(Debug)

///|
/// Create SCENIC input from expression matrix and metadata.
pub fn scenic_input(
  expression : Array[Array[Double]],   // gene x cell
  gene_names : Array[String],
  cell_names : Array[String],
  tf_names : Array[String],
) -> SCENICInput {
  {
    n_genes: gene_names.length(),
    n_cells: cell_names.length(),
    expression: expression,
    gene_names: gene_names,
    cell_names: cell_names,
    tf_names: tf_names,
  }
}

///|
/// SCENIC result: regulons + AUC matrix + binary activity + cell states.
pub struct SCENICResult {
  regulons : Array[ScenicRegulon]
  n_regulons : Int
  n_cells : Int
  auc_matrix : Array[Array[Double]]     // regulon x cell
  regulon_names : Array[String]
  cell_names : Array[String]
  binary_matrix : Array[Array[Int]]     // regulon x cell (0/1)
  thresholds : Array[Double]            // per-regulon binarization threshold
  cell_states : Array[Int]              // cluster label per cell
  n_clusters : Int
  master_regulators : Array[(String, Int)]  // (TF name, cluster index)
}

// ============================================================================
// Step 1: GRN inference (co-expression modules)
// ============================================================================

///|
/// Build co-expression modules by computing per-TF importance scores
/// for all target genes using a variance-based feature importance proxy.
///
/// For each TF, we compute the absolute correlation between the TF's expression
/// vector and each target gene's expression vector across all cells.
/// This mimics GENIE3's importance scoring without requiring a full
/// Random Forest implementation (which is available separately in genie3.mbt).
///
/// Parameters:
///   - input: SCENIC input with expression matrix (gene x cell)
///   - min_targets: minimum targets to form a module (default 10)
///   - top_k: keep top-k targets per TF (default 50)
pub fn build_coexpression_modules(
  input : SCENICInput,
  min_targets? : Int = 10,
  top_k? : Int = 50,
) -> Array[CoExpressionModule] {
  let modules : Array[CoExpressionModule] = []
  // Build gene name -> row index map
  let gene_indices_keys : Array[String] = []
  let gene_indices_vals : Array[Int] = []
  for i in 0..<input.n_genes {
    gene_indices_keys.push(input.gene_names[i])
    gene_indices_vals.push(i)
  }
  // For each TF, compute importance (absolute Pearson correlation) to all genes
  for tf in input.tf_names {
    // Find TF row index
    let mut tf_row = -1
    for k in 0..<gene_indices_keys.length() {
      if gene_indices_keys[k] == tf { tf_row = gene_indices_vals[k]; break }
    }
    if tf_row < 0 { continue }
    let tf_expr = input.expression[tf_row]
    // Compute importance for each gene
    let pairs : Array[(String, Double)] = []
    for g in 0..<input.n_genes {
      if g == tf_row { continue }  // skip self
      let gene_expr = input.expression[g]
      let imp = scenic_abs_correlation(tf_expr, gene_expr)
      pairs.push((input.gene_names[g], imp))
    }
    // Sort by importance descending (simple selection sort for top_k)
    scenic_sort_pairs_desc(pairs)
    // Keep top_k targets
    let keep_n = if top_k < pairs.length() { top_k } else { pairs.length() }
    if keep_n < min_targets { continue }
    let targets : Array[String] = []
    let weights : Array[Double] = []
    for i in 0..<keep_n {
      let (gn, w) = pairs[i]
      targets.push(gn)
      weights.push(w)
    }
    modules.push(co_expression_module(tf, targets, weights))
  }
  modules
}

/// Compute absolute Pearson correlation between two vectors.
fn scenic_abs_correlation(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  if n == 0 || n != b.length() { return 0.0 }
  let mut sum_a = 0.0
  let mut sum_b = 0.0
  for i in 0..<n { sum_a += a[i]; sum_b += b[i] }
  let mean_a = sum_a / n.to_double()
  let mean_b = sum_b / n.to_double()
  let mut cov = 0.0
  let mut var_a = 0.0
  let mut var_b = 0.0
  for i in 0..<n {
    let da = a[i] - mean_a
    let db = b[i] - mean_b
    cov += da * db
    var_a += da * da
    var_b += db * db
  }
  let denom = (var_a * var_b).sqrt()
  if denom < 1.0e-12 { 0.0 }
  else { (cov / denom).abs() }
}

/// Sort (String, Double) pairs by Double value descending (selection sort).
fn scenic_sort_pairs_desc(arr : Array[(String, Double)]) -> Unit {
  let n = arr.length()
  for i in 0..<n {
    let mut max_idx = i
    for j in (i + 1)..<n {
      let (_, vj) = arr[j]
      let (_, vmax) = arr[max_idx]
      if vj > vmax { max_idx = j }
    }
    if max_idx != i {
      let tmp = arr[i]
      arr[i] = arr[max_idx]
      arr[max_idx] = tmp
    }
  }
}

// ============================================================================
// Step 2: Regulon construction (cisTarget-style pruning)
// ============================================================================

///|
/// Prune co-expression modules into regulons by applying a weight threshold.
///
/// In the full SCENIC pipeline, this step uses motif enrichment (cisTarget)
/// to keep only direct TF targets. Since motif databases are species-specific
/// and external, this simplified version prunes by importance weight.
///
/// Parameters:
///   - modules: co-expression modules from Step 1
///   - weight_quantile: keep targets above this quantile of weights (default 0.25)
///   - min_targets: minimum targets to keep a regulon (default 10)
pub fn build_regulons(
  modules : Array[CoExpressionModule],
  weight_quantile? : Double = 0.25,
  min_targets? : Int = 10,
) -> Array[ScenicRegulon] {
  let regulons : Array[ScenicRegulon] = []
  for mod in modules {
    // Compute weight threshold (quantile)
    let sorted_weights = mod.weights.copy()
    scenic_sort_doubles_asc(sorted_weights)
    let q_idx = (weight_quantile * sorted_weights.length().to_double()).to_int()
    let threshold = if q_idx >= 0 && q_idx < sorted_weights.length() {
      sorted_weights[q_idx]
    } else {
      0.0
    }
    // Keep targets above threshold
    let keep_targets : Array[String] = []
    let keep_weights : Array[Double] = []
    for i in 0..<mod.targets.length() {
      if mod.weights[i] >= threshold {
        keep_targets.push(mod.targets[i])
        keep_weights.push(mod.weights[i])
      }
    }
    if keep_targets.length() >= min_targets {
      regulons.push(scenic_regulon(mod.tf_name, keep_targets, keep_weights))
    }
  }
  regulons
}

///|
/// Prune a module's targets using a motif ranking vector (cisTarget-style).
///
/// Given a motif ranking (gene names sorted by motif score, best first),
/// keep only module targets that appear in the top `rank_threshold` positions.
///
/// Parameters:
///   - module: the co-expression module
///   - motif_ranked_genes: gene names sorted by motif score (best first)
///   - rank_threshold: keep genes in top N of motif ranking (default: 5% of all genes)
pub fn prune_by_motif_ranking(
  module : CoExpressionModule,
  motif_ranked_genes : Array[String],
  rank_threshold? : Int = 0,
) -> ScenicRegulon {
  let n_ranked = motif_ranked_genes.length()
  let effective_threshold = if rank_threshold > 0 { rank_threshold }
                            else { (n_ranked.to_double() * 0.05).to_int() }
  // Build set of top-ranked genes
  let top_genes_keys : Array[String] = []
  let _top_vals : Array[Int] = []
  for i in 0..<effective_threshold {
    if i < n_ranked {
      top_genes_keys.push(motif_ranked_genes[i])
    }
  }
  // Keep targets that are in top-ranked set
  let keep_targets : Array[String] = []
  let keep_weights : Array[Double] = []
  for i in 0..<module.targets.length() {
    let target = module.targets[i]
    let mut found = false
    for k in 0..<top_genes_keys.length() {
      if top_genes_keys[k] == target { found = true; break }
    }
    if found {
      keep_targets.push(target)
      keep_weights.push(module.weights[i])
    }
  }
  scenic_regulon(module.tf_name, keep_targets, keep_weights)
}

fn scenic_sort_doubles_asc(arr : Array[Double]) -> Unit {
  let n = arr.length()
  for i in 1..<n {
    let key = arr[i]
    let mut j = i - 1
    while j >= 0 && arr[j] > key {
      arr[j + 1] = arr[j]
      j = j - 1
    }
    arr[j + 1] = key
  }
}

// ============================================================================
// Step 3: AUCell scoring (regulon activity per cell)
// ============================================================================

///|
/// Build gene rankings for each cell (genes sorted by expression, highest first).
/// Returns a matrix of gene indices, cell x gene (each row is a cell's ranking).
pub fn build_cell_rankings(
  input : SCENICInput,
) -> Array[Array[Int]] {
  let rankings : Array[Array[Int]] = []
  for c in 0..<input.n_cells {
    // Extract gene expressions for this cell
    let expr : Array[Double] = []
    for g in 0..<input.n_genes {
      expr.push(input.expression[g][c])
    }
    // Create index array and sort by expression descending
    let indices : Array[Int] = []
    for g in 0..<input.n_genes { indices.push(g) }
    // Selection sort by expression descending
    for i in 0..<input.n_genes {
      let mut max_idx = i
      for j in (i + 1)..<input.n_genes {
        if expr[indices[j]] > expr[indices[max_idx]] { max_idx = j }
      }
      if max_idx != i {
        let tmp = indices[i]
        indices[i] = indices[max_idx]
        indices[max_idx] = tmp
      }
    }
    rankings.push(indices)
  }
  rankings
}

///|
/// Compute AUC for a gene set (regulon targets) in one cell's ranking.
///
/// The AUC is the area under the recovery curve, computed over the top
/// `auc_max_rank` genes in the cell's expression ranking.
///
/// Parameters:
///   - gene_ranking: gene indices sorted by expression (highest first)
///   - target_indices: set of gene indices that are regulon targets
///   - auc_max_rank: number of top genes to consider (default: 5% of genes)
fn compute_auc(
  gene_ranking : Array[Int],
  target_indices : Array[Int],
  auc_max_rank : Int,
) -> Double {
  let n_targets = target_indices.length()
  if n_targets == 0 { return 0.0 }
  let n_genes = gene_ranking.length()
  let max_rank = if auc_max_rank > 0 && auc_max_rank < n_genes { auc_max_rank }
                 else { n_genes }
  // Build a set of target gene indices for O(1) lookup
  let target_set : Array[Int] = target_indices.copy()
  // Count how many targets appear in top max_rank positions
  let mut hits = 0
  let mut recovery_area = 0.0
  let prev_recovery = 0
  for i in 0..<max_rank {
    let gene_idx = gene_ranking[i]
    let mut is_target = false
    for t in 0..<n_targets {
      if target_set[t] == gene_idx { is_target = true; break }
    }
    if is_target {
      hits = hits + 1
    }
    // Recovery at position i = hits / n_targets
    // Area = sum of recovery values at each position
    recovery_area += hits.to_double() / n_targets.to_double()
  }
  // Normalize: divide by max_rank (perfect recovery = 1.0)
  recovery_area / max_rank.to_double()
}

///|
/// Compute regulon activity (AUC) matrix: regulon x cell.
///
/// Parameters:
///   - regulons: list of regulons from Step 2
///   - input: SCENIC input with expression matrix
///   - auc_threshold_pct: percentage of top genes for AUC calculation (default 0.05)
pub fn compute_regulon_activity(
  regulons : Array[ScenicRegulon],
  input : SCENICInput,
  auc_threshold_pct? : Double = 0.05,
) -> Array[Array[Double]] {
  // Build cell rankings
  let cell_rankings = build_cell_rankings(input)
  let auc_max_rank = (input.n_genes.to_double() * auc_threshold_pct).to_int()
  let effective_max = if auc_max_rank < 1 { 1 } else { auc_max_rank }
  // Build gene name -> index map
  let gene_keys : Array[String] = []
  let gene_vals : Array[Int] = []
  for g in 0..<input.n_genes {
    gene_keys.push(input.gene_names[g])
    gene_vals.push(g)
  }
  // For each regulon, compute AUC per cell
  let auc_matrix : Array[Array[Double]] = []
  for reg in regulons {
    // Convert target names to indices
    let target_indices : Array[Int] = []
    for tname in reg.targets {
      let mut found = -1
      for k in 0..<gene_keys.length() {
        if gene_keys[k] == tname { found = gene_vals[k]; break }
      }
      if found >= 0 { target_indices.push(found) }
    }
    // Compute AUC for each cell
    let row : Array[Double] = []
    for c in 0..<input.n_cells {
      let auc = compute_auc(cell_rankings[c], target_indices, effective_max)
      row.push(auc)
    }
    auc_matrix.push(row)
  }
  auc_matrix
}

// ============================================================================
// Step 4: Binarization and cell state assignment
// ============================================================================

///|
/// Binarize the AUC matrix using per-regulon thresholds.
///
/// Uses a simple "mean + 0.5 * std" threshold as default. For bimodal
/// distributions, a k-means (k=2) based threshold is also available.
///
/// Returns (binary_matrix, thresholds).
pub fn binarize_activity(
  auc_matrix : Array[Array[Double]],
  n_cells : Int,
  method? : BinarizeMethod = BinarizeMethod::MeanStd,
) -> (Array[Array[Int]], Array[Double]) {
  let n_regulons = auc_matrix.length()
  let binary : Array[Array[Int]] = []
  let thresholds : Array[Double] = []
  for r in 0..<n_regulons {
    let row = auc_matrix[r]
    let threshold = match method {
      BinarizeMethod::MeanStd => {
        let mut sum = 0.0
        for v in row { sum += v }
        let mean = sum / n_cells.to_double()
        let mut sq_sum = 0.0
        for v in row { sq_sum += (v - mean) * (v - mean) }
        let std = (sq_sum / n_cells.to_double()).sqrt()
        mean + 0.5 * std
      }
      BinarizeMethod::KMeans2 => {
        scenic_kmeans_threshold(row)
      }
      BinarizeMethod::Median => {
        let sorted = row.copy()
        scenic_sort_doubles_asc(sorted)
        sorted[n_cells / 2]
      }
    }
    thresholds.push(threshold)
    let bin_row : Array[Int] = []
    for c in 0..<n_cells {
      if row[c] > threshold { bin_row.push(1) } else { bin_row.push(0) }
    }
    binary.push(bin_row)
  }
  (binary, thresholds)
}

///|
/// Binarization method enum.
pub enum BinarizeMethod {
  MeanStd
  KMeans2
  Median
}

///|
/// Mean-std binarization method.
pub fn binarize_mean_std() -> BinarizeMethod {
  BinarizeMethod::MeanStd
}

///|
/// K-means (k=2) binarization method.
pub fn binarize_kmeans2() -> BinarizeMethod {
  BinarizeMethod::KMeans2
}

///|
/// Median binarization method.
pub fn binarize_median() -> BinarizeMethod {
  BinarizeMethod::Median
}

/// Simple 1D k-means with k=2 to find threshold for bimodal distribution.
fn scenic_kmeans_threshold(values : Array[Double]) -> Double {
  let n = values.length()
  if n == 0 { return 0.0 }
  if n == 1 { return values[0] }
  // Initialize centroids: min and max
  let mut c1 = values[0]
  let mut c2 = values[0]
  for v in values {
    if v < c1 { c1 = v }
    if v > c2 { c2 = v }
  }
  if c1 == c2 { return c1 }
  // Iterate k-means
  for _iter in 0..<20 {
    let mut sum1 = 0.0
    let mut count1 = 0
    let mut sum2 = 0.0
    let mut count2 = 0
    for v in values {
      let d1 = (v - c1).abs()
      let d2 = (v - c2).abs()
      if d1 < d2 { sum1 += v; count1 = count1 + 1 }
      else { sum2 += v; count2 = count2 + 1 }
    }
    if count1 == 0 || count2 == 0 { break }
    let new_c1 = sum1 / count1.to_double()
    let new_c2 = sum2 / count2.to_double()
    if (new_c1 - c1).abs() < 1.0e-9 && (new_c2 - c2).abs() < 1.0e-9 { break }
    c1 = new_c1
    c2 = new_c2
  }
  // Threshold = midpoint between two centroids
  (c1 + c2) / 2.0
}

///|
/// Assign cell states based on binary regulon activity matrix.
///
/// Cells are clustered by their binary regulon activity profiles using
/// a simple Hamming-distance-based assignment: each cell is assigned
/// to the regulon with the highest activity (master regulator).
/// Cells with no active regulon are assigned to cluster 0 (inactive).
///
/// Returns (cell_state_labels, n_clusters, master_regulators).
pub fn assign_cell_states(
  binary_matrix : Array[Array[Int]],
  regulon_names : Array[String],
  n_cells : Int,
) -> (Array[Int], Int, Array[(String, Int)]) {
  let n_regulons = binary_matrix.length()
  let labels : Array[Int] = []
  // Count active regulons per cell
  for c in 0..<n_cells {
    let mut best_reg = -1
    let mut best_score = 0
    for r in 0..<n_regulons {
      if binary_matrix[r][c] == 1 {
        // Assign to first active regulon (or could use AUC for tie-breaking)
        if best_reg < 0 { best_reg = r }
        best_score = best_score + 1
      }
    }
    if best_reg < 0 {
      labels.push(0)  // no active regulon → cluster 0
    } else {
      labels.push(best_reg + 1)  // cluster = regulon index + 1
    }
  }
  // Count clusters
  let mut max_cluster = 0
  for l in labels {
    if l > max_cluster { max_cluster = l }
  }
  // Master regulators: (regulon_name, cluster_index)
  let masters : Array[(String, Int)] = []
  for r in 0..<n_regulons {
    let cluster = r + 1
    let mut count = 0
    for c in 0..<n_cells {
      if labels[c] == cluster { count = count + 1 }
    }
    if count > 0 {
      masters.push((regulon_names[r], cluster))
    }
  }
  (labels, max_cluster, masters)
}

// ============================================================================
// Top-level pipeline
// ============================================================================

///|
/// Run the full SCENIC pipeline end-to-end.
///
/// Parameters:
///   - input: expression matrix (gene x cell) + gene/cell/TF names
///   - min_targets: minimum targets per module (default 10)
///   - top_k: top-k targets per TF in modules (default 50)
///   - weight_quantile: pruning quantile for regulon construction (default 0.25)
///   - auc_threshold_pct: AUC calculation top-gene percentage (default 0.05)
///   - binarize_method: threshold method for binarization
pub fn run_scenic(
  input : SCENICInput,
  min_targets? : Int = 10,
  top_k? : Int = 50,
  weight_quantile? : Double = 0.25,
  auc_threshold_pct? : Double = 0.05,
  binarize_method? : BinarizeMethod = BinarizeMethod::MeanStd,
) -> SCENICResult {
  // Step 1: Build co-expression modules
  let modules = build_coexpression_modules(
    input, min_targets=min_targets, top_k=top_k,
  )
  // Step 2: Build regulons
  let regulons = build_regulons(
    modules, weight_quantile=weight_quantile, min_targets=min_targets,
  )
  // Step 3: Compute regulon activity (AUC matrix)
  let auc_matrix = compute_regulon_activity(
    regulons, input, auc_threshold_pct=auc_threshold_pct,
  )
  // Step 4: Binarize and assign cell states
  let (binary_matrix, thresholds) = binarize_activity(
    auc_matrix, input.n_cells, method=binarize_method,
  )
  let regulon_names = regulons.map(fn (r) -> String { r.tf_name })
  let (cell_states, n_clusters, masters) = assign_cell_states(
    binary_matrix, regulon_names, input.n_cells,
  )
  {
    regulons: regulons,
    n_regulons: regulons.length(),
    n_cells: input.n_cells,
    auc_matrix: auc_matrix,
    regulon_names: regulon_names,
    cell_names: input.cell_names,
    binary_matrix: binary_matrix,
    thresholds: thresholds,
    cell_states: cell_states,
    n_clusters: n_clusters,
    master_regulators: masters,
  }
}

// ============================================================================
// Synthetic data generator
// ============================================================================

///|
/// Generate synthetic SCENIC test data with known TF-target relationships.
///
/// Creates `n_genes` genes across `n_cells` cells, with `n_tfs` transcription
/// factors. Each TF controls a set of target genes whose expression correlates
/// with the TF. Tumour-like cells (first half) have elevated TF1 activity,
/// while normal-like cells (second half) have elevated TF2 activity.
///
/// Returns SCENICInput ready for `run_scenic`.
pub fn scenic_sample_data(
  n_genes? : Int = 100,
  n_cells? : Int = 50,
  n_tfs? : Int = 5,
  seed? : Int = 42,
) -> SCENICInput {
  let mut rng = seed
  let next_rand = fn () -> Double {
    rng = (rng * 1103515245 + 12345) & 0x7fffffff
    (rng & 0xffffff).to_double() / 16777215.0
  }
  // Gene names: TF1..TFn, G1..Gm
  let gene_names : Array[String] = []
  let tf_names : Array[String] = []
  for i in 0..<n_tfs {
    let name = "TF" + (i + 1).to_string()
    gene_names.push(name)
    tf_names.push(name)
  }
  for i in 0..<(n_genes - n_tfs) {
    gene_names.push("G" + (i + 1).to_string())
  }
  let cell_names : Array[String] = []
  for i in 0..<n_cells {
    cell_names.push("Cell_" + (i + 1).to_string())
  }
  // Assign targets to each TF (roughly equal split of non-TF genes)
  let n_targets_per_tf = (n_genes - n_tfs) / n_tfs
  // Build TF activity profiles: TF1 high in first half, TF2 high in second half
  let tf_activity : Array[Array[Double]] = []
  for t in 0..<n_tfs {
    let profile : Array[Double] = []
    for c in 0..<n_cells {
      let base = if t == 0 && c < n_cells / 2 { 5.0 }
                 else if t == 1 && c >= n_cells / 2 { 5.0 }
                 else { 1.0 + next_rand() * 0.5 }
      profile.push(base)
    }
    tf_activity.push(profile)
  }
  // Build expression matrix: gene x cell
  let expression : Array[Array[Double]] = []
  // TF expression = activity + noise
  for t in 0..<n_tfs {
    let row : Array[Double] = []
    for c in 0..<n_cells {
      row.push(tf_activity[t][c] + (next_rand() - 0.5) * 0.5)
    }
    expression.push(row)
  }
  // Target gene expression = correlated with controlling TF + noise
  for g in n_tfs..<n_genes {
    let tf_idx = (g - n_tfs) / n_targets_per_tf
    let effective_tf = if tf_idx < n_tfs { tf_idx } else { 0 }
    let row : Array[Double] = []
    for c in 0..<n_cells {
      let val = tf_activity[effective_tf][c] * 0.8 + (next_rand() - 0.5) * 1.0
      row.push(if val < 0.0 { 0.0 } else { val })
    }
    expression.push(row)
  }
  scenic_input(expression, gene_names, cell_names, tf_names)
}

// ============================================================================
// Query methods on SCENICResult
// ============================================================================

///|
/// Get the AUC score for a specific regulon and cell.
pub fn SCENICResult::auc_at(self : SCENICResult, regulon_idx : Int, cell_idx : Int) -> Double {
  if regulon_idx < 0 || regulon_idx >= self.n_regulons { return 0.0 }
  if cell_idx < 0 || cell_idx >= self.n_cells { return 0.0 }
  self.auc_matrix[regulon_idx][cell_idx]
}

///|
/// Get the binary activity (0 or 1) for a specific regulon and cell.
pub fn SCENICResult::binary_at(self : SCENICResult, regulon_idx : Int, cell_idx : Int) -> Int {
  if regulon_idx < 0 || regulon_idx >= self.n_regulons { return 0 }
  if cell_idx < 0 || cell_idx >= self.n_cells { return 0 }
  self.binary_matrix[regulon_idx][cell_idx]
}

///|
/// Get the cell state (cluster label) for a specific cell.
pub fn SCENICResult::cell_state(self : SCENICResult, cell_idx : Int) -> Int {
  if cell_idx < 0 || cell_idx >= self.n_cells { return 0 }
  self.cell_states[cell_idx]
}

///|
/// Return a summary string of the SCENIC result.
pub fn SCENICResult::summary(self : SCENICResult) -> String {
  "SCENICResult(regulons=" + self.n_regulons.to_string() +
    ", cells=" + self.n_cells.to_string() +
    ", clusters=" + self.n_clusters.to_string() +
    ", master_regulators=" + self.master_regulators.length().to_string() + ")"
}

///|
/// Get regulon targets as an array of gene names.
pub fn SCENICResult::regulon_targets(self : SCENICResult, regulon_idx : Int) -> Array[String] {
  if regulon_idx < 0 || regulon_idx >= self.n_regulons { return [] }
  self.regulons[regulon_idx].targets.copy()
}

///|
/// Get the top-N regulons by mean AUC across all cells.
pub fn SCENICResult::top_regulons(self : SCENICResult, n? : Int = 5) -> Array[(String, Double)] {
  let pairs : Array[(String, Double)] = []
  for r in 0..<self.n_regulons {
    let mut sum = 0.0
    for c in 0..<self.n_cells { sum += self.auc_matrix[r][c] }
    let mean = sum / self.n_cells.to_double()
    pairs.push((self.regulon_names[r], mean))
  }
  scenic_sort_pairs_desc(pairs)
  let result : Array[(String, Double)] = []
  let limit = if n < pairs.length() { n } else { pairs.length() }
  for i in 0..<limit {
    result.push(pairs[i])
  }
  result
}
