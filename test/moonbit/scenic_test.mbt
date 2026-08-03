///|
/// Test file for SCENIC module.

// ---------------------------------------------------------------------------
// Data structure creation
// ---------------------------------------------------------------------------

test "scenic_co_expression_module_creation" {
  let mod_ = @src.co_expression_module(
    "TP53",
    ["G1", "G2", "G3"],
    [0.9, 0.7, 0.5],
  )
  assert_eq(mod_.tf_name, "TP53")
  assert_eq(mod_.targets.length(), 3)
  assert_eq(mod_.weights[0], 0.9)
}

test "scenic_regulon_creation" {
  let reg = @src.scenic_regulon(
    "MYC",
    ["G1", "G2", "G3", "G4"],
    [1.0, 0.8, 0.6, 0.4],
  )
  assert_eq(reg.tf_name, "MYC")
  assert_eq(reg.n_targets, 4)
  assert_eq(reg.targets[2], "G3")
}

test "scenic_input_creation" {
  let expr = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]  // 2 genes x 3 cells
  let input = @src.scenic_input(
    expr, ["TF1", "G1"], ["C1", "C2", "C3"], ["TF1"],
  )
  assert_eq(input.n_genes, 2)
  assert_eq(input.n_cells, 3)
  assert_eq(input.tf_names.length(), 1)
}

// ---------------------------------------------------------------------------
// Binarize method constructors
// ---------------------------------------------------------------------------

test "scenic_binarize_method_constructors" {
  match @src.binarize_mean_std() {
    @src.BinarizeMethod::MeanStd => assert_true(true)
    _ => assert_true(false)
  }
  match @src.binarize_kmeans2() {
    @src.BinarizeMethod::KMeans2 => assert_true(true)
    _ => assert_true(false)
  }
  match @src.binarize_median() {
    @src.BinarizeMethod::Median => assert_true(true)
    _ => assert_true(false)
  }
}

// ---------------------------------------------------------------------------
// Synthetic data generator
// ---------------------------------------------------------------------------

test "scenic_sample_data_shape" {
  let input = @src.scenic_sample_data(
    n_genes=60, n_cells=40, n_tfs=4, seed=10,
  )
  assert_eq(input.n_genes, 60)
  assert_eq(input.n_cells, 40)
  assert_eq(input.tf_names.length(), 4)
  assert_eq(input.gene_names[0], "TF1")
  assert_eq(input.gene_names[4], "G1")
  assert_eq(input.expression.length(), 60)     // gene x cell
  assert_eq(input.expression[0].length(), 40)
}

// ---------------------------------------------------------------------------
// Step 1: Co-expression modules
// ---------------------------------------------------------------------------

test "scenic_build_coexpression_modules_basic" {
  let input = @src.scenic_sample_data(
    n_genes=50, n_cells=30, n_tfs=3, seed=42,
  )
  let modules = @src.build_coexpression_modules(
    input, min_targets=5, top_k=20,
  )
  // Should produce a module for each TF
  assert_eq(modules.length(), 3)
  for mod_ in modules {
    assert_true(mod_.targets.length() >= 5)
    assert_true(mod_.targets.length() <= 20)
    assert_eq(mod_.weights.length(), mod_.targets.length())
  }
}

test "scenic_build_coexpression_modules_skip_unknown_tf" {
  // TF name not in gene_names -> should be skipped
  let expr = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let input = @src.scenic_input(
    expr, ["G1", "G2", "G3"], ["C1", "C2"], ["NONEXISTENT"],
  )
  let modules = @src.build_coexpression_modules(input, min_targets=1, top_k=5)
  assert_eq(modules.length(), 0)
}

// ---------------------------------------------------------------------------
// Step 2: Regulon construction
// ---------------------------------------------------------------------------

test "scenic_build_regulons_prunes_by_weight" {
  let input = @src.scenic_sample_data(
    n_genes=50, n_cells=30, n_tfs=3, seed=7,
  )
  let modules = @src.build_coexpression_modules(
    input, min_targets=5, top_k=20,
  )
  let regulons = @src.build_regulons(
    modules, weight_quantile=0.5, min_targets=3,
  )
  // Each regulon should have fewer or equal targets than its module
  assert_eq(regulons.length(), modules.length())
  for i in 0..<regulons.length() {
    assert_true(regulons[i].n_targets <= modules[i].targets.length())
    assert_true(regulons[i].n_targets >= 3)
  }
}

test "scenic_build_regulons_filters_small_modules" {
  let modules = [
    @src.co_expression_module("TF1", ["G1", "G2"], [0.9, 0.8]),
    @src.co_expression_module("TF2", ["G1", "G2", "G3", "G4", "G5"], [0.9, 0.8, 0.7, 0.6, 0.5]),
  ]
  let regulons = @src.build_regulons(modules, weight_quantile=0.0, min_targets=5)
  // TF1 module has only 2 targets -> filtered out
  assert_eq(regulons.length(), 1)
  assert_eq(regulons[0].tf_name, "TF2")
}

test "scenic_prune_by_motif_ranking" {
  let mod_ = @src.co_expression_module(
    "TF1",
    ["G1", "G2", "G3", "G4", "G5"],
    [0.9, 0.8, 0.7, 0.6, 0.5],
  )
  // Motif ranking: G1, G3, G5 are in top 3
  let motif_ranked = ["G1", "G3", "G5", "G2", "G4", "G6", "G7", "G8"]
  let reg = @src.prune_by_motif_ranking(mod_, motif_ranked, rank_threshold=3)
  assert_eq(reg.tf_name, "TF1")
  // G1, G3, G5 should be kept; G2, G4 should be pruned
  assert_eq(reg.n_targets, 3)
  assert_eq(reg.targets[0], "G1")
  assert_eq(reg.targets[1], "G3")
  assert_eq(reg.targets[2], "G5")
}

// ---------------------------------------------------------------------------
// Step 3: AUCell scoring
// ---------------------------------------------------------------------------

test "scenic_build_cell_rankings" {
  let expr = [
    [5.0, 1.0, 3.0],   // gene 0: high in cell 0
    [2.0, 4.0, 1.0],   // gene 1: high in cell 1
    [1.0, 3.0, 2.0],   // gene 2
  ]
  let input = @src.scenic_input(
    expr, ["G0", "G1", "G2"], ["C0", "C1", "C2"], [],
  )
  let rankings = @src.build_cell_rankings(input)
  // Cell 0: gene 0 has expr 5.0 (highest), gene 1 has 2.0, gene 2 has 1.0
  assert_eq(rankings[0][0], 0)   // gene 0 first
  assert_eq(rankings[0][1], 1)   // gene 1 second
  assert_eq(rankings[0][2], 2)   // gene 2 third
}

test "scenic_compute_regulon_activity_basic" {
  let input = @src.scenic_sample_data(
    n_genes=40, n_cells=20, n_tfs=2, seed=99,
  )
  let modules = @src.build_coexpression_modules(
    input, min_targets=5, top_k=15,
  )
  let regulons = @src.build_regulons(modules, weight_quantile=0.3, min_targets=3)
  let auc_matrix = @src.compute_regulon_activity(
    regulons, input, auc_threshold_pct=0.1,
  )
  assert_eq(auc_matrix.length(), regulons.length())
  assert_eq(auc_matrix[0].length(), input.n_cells)
  // AUC values should be in [0, 1]
  for row in auc_matrix {
    for v in row {
      assert_true(v >= 0.0 && v <= 1.0)
    }
  }
}

// ---------------------------------------------------------------------------
// Step 4: Binarization
// ---------------------------------------------------------------------------

test "scenic_binarize_mean_std" {
  let auc_matrix = [
    [0.1, 0.2, 0.8, 0.9, 0.1, 0.85],  // bimodal: low {0.1,0.2,0.1}, high {0.8,0.9,0.85}
  ]
  let (binary, thresholds) = @src.binarize_activity(
    auc_matrix, 6, method=@src.binarize_mean_std(),
  )
  assert_eq(thresholds.length(), 1)
  assert_true(thresholds[0] > 0.3 && thresholds[0] < 0.7)
  // High-AUC cells should be 1, low-AUC cells should be 0
  assert_eq(binary[0][2], 1)   // 0.8 > threshold
  assert_eq(binary[0][0], 0)   // 0.1 < threshold
}

test "scenic_binarize_kmeans2" {
  let auc_matrix = [
    [0.1, 0.1, 0.1, 0.9, 0.9, 0.9],  // clearly bimodal
  ]
  let (binary, thresholds) = @src.binarize_activity(
    auc_matrix, 6, method=@src.binarize_kmeans2(),
  )
  // K-means should separate at ~0.5
  assert_true(thresholds[0] > 0.3 && thresholds[0] < 0.7)
  assert_eq(binary[0][0], 0)
  assert_eq(binary[0][3], 1)
}

test "scenic_binarize_median" {
  let auc_matrix = [
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.9],
  ]
  let (binary, _thresholds) = @src.binarize_activity(
    auc_matrix, 6, method=@src.binarize_median(),
  )
  // Median of sorted [0.1,0.2,0.3,0.4,0.5,0.9] = 0.35
  // Values > 0.35: 0.4, 0.5, 0.9 -> 1
  assert_eq(binary[0][5], 1)   // 0.9
  assert_eq(binary[0][0], 0)   // 0.1
}

// ---------------------------------------------------------------------------
// Cell state assignment
// ---------------------------------------------------------------------------

test "scenic_assign_cell_states_basic" {
  // 2 regulons, 4 cells
  // Regulon 0 active in cells 0,1; Regulon 1 active in cells 2,3
  let binary = [
    [1, 1, 0, 0],
    [0, 0, 1, 1],
  ]
  let (labels, n_clusters, masters) = @src.assign_cell_states(
    binary, ["TF1", "TF2"], 4,
  )
  assert_eq(labels[0], 1)   // cell 0 -> cluster 1 (TF1)
  assert_eq(labels[1], 1)   // cell 1 -> cluster 1
  assert_eq(labels[2], 2)   // cell 2 -> cluster 2 (TF2)
  assert_eq(labels[3], 2)   // cell 3 -> cluster 2
  assert_eq(n_clusters, 2)
  assert_eq(masters.length(), 2)
}

test "scenic_assign_cell_states_no_active_regulon" {
  let binary = [
    [0, 0, 0],
    [0, 0, 0],
  ]
  let (labels, _n, masters) = @src.assign_cell_states(
    binary, ["TF1", "TF2"], 3,
  )
  // All cells should be in cluster 0 (no active regulon)
  for l in labels { assert_eq(l, 0) }
  assert_eq(masters.length(), 0)
}

// ---------------------------------------------------------------------------
// End-to-end pipeline
// ---------------------------------------------------------------------------

test "scenic_run_pipeline_end_to_end" {
  let input = @src.scenic_sample_data(
    n_genes=60, n_cells=40, n_tfs=4, seed=2025,
  )
  let result = @src.run_scenic(
    input, min_targets=5, top_k=20,
    weight_quantile=0.3, auc_threshold_pct=0.1,
    binarize_method=@src.binarize_mean_std(),
  )
  assert_eq(result.n_cells, 40)
  assert_true(result.n_regulons > 0)
  assert_eq(result.auc_matrix.length(), result.n_regulons)
  assert_eq(result.auc_matrix[0].length(), 40)
  assert_eq(result.binary_matrix.length(), result.n_regulons)
  assert_eq(result.cell_states.length(), 40)
  assert_eq(result.thresholds.length(), result.n_regulons)
}

test "scenic_result_summary" {
  let input = @src.scenic_sample_data(
    n_genes=40, n_cells=20, n_tfs=3, seed=1,
  )
  let result = @src.run_scenic(input)
  let s = result.summary()
  assert_true(s.contains("regulons="))
  assert_true(s.contains("cells="))
  assert_true(s.contains("clusters="))
}

test "scenic_result_auc_at" {
  let input = @src.scenic_sample_data(
    n_genes=40, n_cells=10, n_tfs=2, seed=5,
  )
  let result = @src.run_scenic(input)
  // Valid indices should return a value in [0,1]
  if result.n_regulons > 0 {
    let val = result.auc_at(0, 0)
    assert_true(val >= 0.0 && val <= 1.0)
  }
  // Out of range should return 0.0
  assert_eq(result.auc_at(-1, 0), 0.0)
  assert_eq(result.auc_at(999, 0), 0.0)
}

test "scenic_result_top_regulons" {
  let input = @src.scenic_sample_data(
    n_genes=50, n_cells=30, n_tfs=4, seed=8,
  )
  let result = @src.run_scenic(input)
  let top = result.top_regulons(n=3)
  assert_true(top.length() <= 3)
  if top.length() >= 2 {
    let (_, s0) = top[0]
    let (_, s1) = top[1]
    assert_true(s0 >= s1)   // sorted descending
  }
}

test "scenic_result_regulon_targets" {
  let input = @src.scenic_sample_data(
    n_genes=40, n_cells=20, n_tfs=2, seed=3,
  )
  let result = @src.run_scenic(input)
  if result.n_regulons > 0 {
    let targets = result.regulon_targets(0)
    assert_true(targets.length() > 0)
  }
  // Out of range returns empty
  assert_eq(result.regulon_targets(-1).length(), 0)
}
