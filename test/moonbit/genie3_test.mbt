///|
/// Tests for Bioconductor GENIE3 module - Gene regulatory network inference.

test "genie3_tree_leaf_basic" {
  let leaf = @src.genie3_tree_leaf(3.5)
  assert_true(leaf.is_leaf)
  assert_true((leaf.value - 3.5).abs() < 1.0e-9)
  assert_eq(leaf.feature, -1)
  assert_eq(leaf.left_idx, -1)
  assert_eq(leaf.right_idx, -1)
}

test "genie3_tree_split_basic" {
  let left = @src.genie3_tree_leaf(1.0)
  let right = @src.genie3_tree_leaf(2.0)
  let node = @src.genie3_tree_split(0, 0.5, 0, 1, 10.0)
  let _ = (left, right) // suppress unused warnings
  assert_false(node.is_leaf)
  assert_eq(node.feature, 0)
  assert_true((node.threshold - 0.5).abs() < 1.0e-9)
  assert_eq(node.left_idx, 0)
  assert_eq(node.right_idx, 1)
  assert_true((node.importance_gain - 10.0).abs() < 1.0e-9)
}

test "genie3_build_tree_simple" {
  // Simple dataset: y = 2 * x + noise
  let features = [[1.0], [2.0], [3.0], [4.0], [5.0], [6.0]]
  let target = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
  let tree = @src.genie3_build_tree(features, target, [0], max_depth=3)
  // Tree should make at least one split (non-constant target)
  assert_true(tree.nodes.length() > 1)
  // Predictions should be reasonable
  let p1 = tree.predict([1.5])
  let p2 = tree.predict([5.5])
  assert_true(p1 < p2)
}

test "genie3_build_tree_constant_target" {
  // Constant target should produce a single leaf
  let features = [[1.0], [2.0], [3.0]]
  let target = [5.0, 5.0, 5.0]
  let tree = @src.genie3_build_tree(features, target, [0])
  assert_eq(tree.nodes.length(), 1)
  assert_true(tree.nodes[0].is_leaf)
}

test "genie3_feature_importance" {
  // x0 strongly predicts y, x1 doesn't
  let features = [
    [1.0, 10.0],
    [2.0, 20.0],
    [3.0, 5.0],
    [4.0, 15.0],
    [5.0, 25.0],
    [6.0, 8.0],
  ]
  let target = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0] // y = 2 * x0
  let tree = @src.genie3_build_tree(features, target, [0, 1], max_depth=5)
  let importance = tree.feature_importance(2)
  // x0 should have non-zero importance; x1 should be near zero
  assert_true(importance[0] > 0.0)
}

test "genie3_run_basic" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5)
  assert_eq(result.gene_names.length(), 5)
  assert_eq(result.matrix.length(), 5)
  assert_eq(result.matrix[0].length(), 5)
  // n_trees equals number of target genes
  assert_eq(result.n_trees, 5)
  assert_false(result.symmetrized)
  // Edges list should be non-empty
  assert_true(result.edges.length() > 0)
}

test "genie3_run_symmetrize" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5, symmetrize=true)
  assert_true(result.symmetrized)
  // Symmetrized matrix should be symmetric (matrix[i][j] == matrix[j][i])
  let mut i = 0
  while i < 5 {
    let mut j = 0
    while j < 5 {
      assert_true((result.matrix[i][j] - result.matrix[j][i]).abs() < 1.0e-9)
      j = j + 1
    }
    i = i + 1
  }
}

test "genie3_column_normalization" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5)
  // Each column should sum to 1 (or 0 if no edges)
  let mut j = 0
  while j < 5 {
    let mut sum = 0.0
    let mut i = 0
    while i < 5 {
      sum = sum + result.matrix[i][j]
      i = i + 1
    }
    // Allow tiny tolerance for floating-point
    assert_true((sum - 1.0).abs() < 1.0e-6 || sum < 1.0e-9)
    j = j + 1
  }
}

test "genie3_edges_sorted_descending" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5)
  let mut i = 1
  while i < result.edges.length() {
    assert_true(result.edges[i - 1].weight >= result.edges[i].weight)
    i = i + 1
  }
}

test "genie3_top_edges" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5)
  let top3 = @src.genie3_top_edges(result, 3)
  assert_true(top3.length() <= 3)
  // Should match first 3 of result.edges
  let mut i = 0
  while i < top3.length() {
    assert_eq(top3[i].regulator, result.edges[i].regulator)
    assert_eq(top3[i].target, result.edges[i].target)
    i = i + 1
  }
}

test "genie3_top_edges_more_than_available" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5)
  let top100 = @src.genie3_top_edges(result, 100)
  assert_eq(top100.length(), result.edges.length())
}

test "genie3_regulators_of" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5)
  // G2 should have at least one regulator (likely G1)
  let regs = @src.genie3_regulators_of(result, "G2")
  assert_true(regs.length() > 0)
  let mut found_g1 = false
  let mut i = 0
  while i < regs.length() {
    assert_eq(regs[i].target, "G2")
    if regs[i].regulator == "G1" {
      found_g1 = true
    }
    i = i + 1
  }
  assert_true(found_g1)
}

test "genie3_targets_of" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5)
  // G1 should regulate at least one of G2 or G3
  let tgts = @src.genie3_targets_of(result, "G1")
  assert_true(tgts.length() > 0)
  let mut i = 0
  while i < tgts.length() {
    assert_eq(tgts[i].regulator, "G1")
    i = i + 1
  }
}

test "genie3_regulators_of_unknown_gene" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5)
  let regs = @src.genie3_regulators_of(result, "UNKNOWN_GENE")
  assert_eq(regs.length(), 0)
}

test "genie3_sample_data_shape" {
  let (expr, names) = @src.genie3_sample_data()
  assert_eq(names.length(), 5)
  assert_eq(expr.length(), 30) // 30 samples
  assert_eq(expr[0].length(), 5) // 5 genes
  // Sample data is reproducible (deterministic seed)
  let (expr2, _) = @src.genie3_sample_data()
  assert_true((expr[0][0] - expr2[0][0]).abs() < 1.0e-9)
}

test "genie3_run_no_edges_when_uniform" {
  // All identical samples → no informative splits → no edges
  let expr = [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]
  let names = ["G1", "G2", "G3"]
  let result = @src.genie3_run(expr, names, max_depth=3)
  assert_eq(result.edges.length(), 0)
}

test "genie3_self_loops_excluded" {
  let (expr, names) = @src.genie3_sample_data()
  let result = @src.genie3_run(expr, names, max_depth=5)
  // No edge should have regulator == target
  let mut i = 0
  while i < result.edges.length() {
    assert_not_eq(result.edges[i].regulator, result.edges[i].target)
    i = i + 1
  }
}

test "genie3_predict_after_build" {
  // Build a tree and verify prediction lies within target range
  let features = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 4.0],
    [4.0, 5.0],
    [5.0, 6.0],
    [6.0, 7.0],
  ]
  let target = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let tree = @src.genie3_build_tree(features, target, [0, 1], max_depth=4)
  let p = tree.predict([3.5, 4.5])
  // Prediction should be in range [1, 6]
  assert_true(p >= 0.5 && p <= 6.5)
}
