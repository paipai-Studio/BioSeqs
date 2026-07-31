///|
/// SCENIC module example - Single-Cell Regulatory Network Inference and Clustering.
///
/// Run: moon run examples/scenic_demo
fn main {
  println("=== SCENIC Module Demo ===")
  println("Single-Cell rEgulatory Network Inference and Clustering")
  println("Reference: Aibar et al. (2017) Nature Methods 14:1083-1086")
  println("")

  // ------------------------------------------------------------------
  // 1. Generate synthetic single-cell data with TF-target structure
  // ------------------------------------------------------------------
  println("=== 1. Generate synthetic dataset ===")
  let input = @bio.scenic_sample_data(
    n_genes=80, n_cells=60, n_tfs=5, seed=2025,
  )
  println("Genes: \{input.n_genes}, Cells: \{input.n_cells}, TFs: \{input.tf_names.length()}")
  println("TF names: \{input.tf_names}")
  println("")

  // ------------------------------------------------------------------
  // 2. Step 1: Build co-expression modules
  // ------------------------------------------------------------------
  println("=== 2. Step 1: GRN inference (co-expression modules) ===")
  let modules = @bio.build_coexpression_modules(
    input, min_targets=5, top_k=25,
  )
  println("Modules built: \{modules.length()}")
  let mut mi = 0
  while mi < modules.length() {
    let mod_ = modules[mi]
    println("  \{mod_.tf_name}: \{mod_.targets.length()} targets (top weight: \{mod_.weights[0]})")
    mi = mi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Step 2: Build regulons (cisTarget-style pruning)
  // ------------------------------------------------------------------
  println("=== 3. Step 2: Regulon construction (pruning) ===")
  let regulons = @bio.build_regulons(
    modules, weight_quantile=0.25, min_targets=5,
  )
  println("Regulons after pruning: \{regulons.length()}")
  let mut ri = 0
  while ri < regulons.length() {
    let reg = regulons[ri]
    println("  \{reg.tf_name}: \{reg.n_targets} direct targets")
    ri = ri + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Step 3: AUCell scoring (regulon activity per cell)
  // ------------------------------------------------------------------
  println("=== 4. Step 3: AUCell regulon activity scoring ===")
  let auc_matrix = @bio.compute_regulon_activity(
    regulons, input, auc_threshold_pct=0.1,
  )
  println("AUC matrix: \{auc_matrix.length()} regulons x \{input.n_cells} cells")
  // Show mean AUC per regulon
  let mut ri2 = 0
  while ri2 < regulons.length() {
    let mut sum = 0.0
    let mut ci = 0
    while ci < input.n_cells {
      sum = sum + auc_matrix[ri2][ci]
      ci = ci + 1
    }
    let mean_auc = sum / input.n_cells.to_double()
    println("  \{regulons[ri2].tf_name}: mean AUC = \{mean_auc}")
    ri2 = ri2 + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Step 4: Binarization and cell state assignment
  // ------------------------------------------------------------------
  println("=== 5. Step 4: Binarization + cell state assignment ===")
  let (binary, thresholds) = @bio.binarize_activity(
    auc_matrix, input.n_cells, method=@bio.binarize_kmeans2(),
  )
  let (cell_states, n_clusters, masters) = @bio.assign_cell_states(
    binary, regulons.map(fn (r) -> String { r.tf_name }), input.n_cells,
  )
  println("Clusters identified: \{n_clusters}")
  println("Master regulators:")
  let mut mk = 0
  while mk < masters.length() {
    let (tf, cluster) = masters[mk]
    // Count cells in this cluster
    let mut count = 0
    let mut ci2 = 0
    while ci2 < input.n_cells {
      if cell_states[ci2] == cluster { count = count + 1 }
      ci2 = ci2 + 1
    }
    println("  \{tf} -> cluster \{cluster} (\{count} cells)")
    mk = mk + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Full pipeline in one call
  // ------------------------------------------------------------------
  println("=== 6. Full SCENIC pipeline (run_scenic) ===")
  let result = @bio.run_scenic(
    input, min_targets=5, top_k=25,
    weight_quantile=0.25, auc_threshold_pct=0.1,
    binarize_method=@bio.binarize_mean_std(),
  )
  println(result.summary())
  println("")

  // ------------------------------------------------------------------
  // 7. Top regulons by activity
  // ------------------------------------------------------------------
  println("=== 7. Top regulons by mean AUC ===")
  let top = result.top_regulons(n=5)
  let mut ti = 0
  while ti < top.length() {
    let (name, score) = top[ti]
    println("  \{ti + 1}. \{name}: \{score}")
    ti = ti + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Spot-check: cell state distribution
  // ------------------------------------------------------------------
  println("=== 8. Cell state distribution ===")
  let mut c = 0
  while c <= n_clusters {
    let mut count = 0
    let mut ci3 = 0
    while ci3 < input.n_cells {
      if result.cell_states[ci3] == c { count = count + 1 }
      ci3 = ci3 + 1
    }
    if count > 0 {
      println("  Cluster \{c}: \{count} cells")
    }
    c = c + 1
  }
  println("")

  println("=== Demo complete ===")
}
