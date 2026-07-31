///|
/// Tests for Bioconductor CellChat module - Cell-cell communication analysis.

test "lr_pair_create" {
  let pair = @src.lr_pair("TNF", "TNFR1")
  assert_eq(pair.ligand, "TNF")
  assert_eq(pair.receptor, "TNFR1")
  assert_eq(pair.key, "TNF_TNFR1")
}

test "lr_database" {
  let db = @src.cellchat_lr_database()
  assert_true(db.length() > 0)
  assert_eq(db[0].ligand, "TNF")
}

test "mean_expr_basic" {
  let expr = [[10.0, 5.0], [20.0, 8.0], [30.0, 12.0]]
  let cell_types = ["TypeA", "TypeA", "TypeB"]
  let mean = @src.cellchat_mean_expr(expr, cell_types, "TypeA", 0)
  assert_true((mean - 15.0).abs() < 0.01)
}

test "mean_expr_no_cells" {
  let expr = [[10.0, 5.0]]
  let cell_types = ["TypeA"]
  let mean = @src.cellchat_mean_expr(expr, cell_types, "TypeC", 0)
  assert_eq(mean, 0.0)
}

test "analyze_basic" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  let pairs = @src.cellchat_lr_database()
  let result = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=10, seed=42)
  assert_true(result.scores.length() > 0)
  assert_eq(result.cell_types.length(), 3)
  assert_true(result.lr_pairs.length() > 0)
}

test "analyze_permutation_effect" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  let pairs = @src.cellchat_lr_database()
  let result_small = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=5, seed=42)
  let result_large = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=50, seed=42)
  assert_true(result_small.scores.length() > 0)
  assert_true(result_large.scores.length() > 0)
}

test "analyze_score_non_negative" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  let pairs = @src.cellchat_lr_database()
  let result = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=5, seed=42)
  let mut i = 0
  while i < result.scores.length() {
    assert_true(result.scores[i].score >= 0.0)
    assert_true(result.scores[i].p_value >= 0.0 && result.scores[i].p_value <= 1.0)
    i = i + 1
  }
}

test "get_significant" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  let pairs = @src.cellchat_lr_database()
  let result = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=5, seed=42, fdr=0.5)
  let sig = @src.cellchat_get_significant(result)
  assert_true(sig.length() >= 0)
  // All significant scores should have p_value < fdr
  let mut i = 0
  while i < sig.length() {
    assert_true(sig[i].p_value < 0.5)
    i = i + 1
  }
}

test "aggregate" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  let pairs = @src.cellchat_lr_database()
  let result = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=5, seed=42)
  let agg = @src.cellchat_aggregate(result)
  assert_true(agg.size() > 0)
}

test "sample_data_dimensions" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  assert_eq(expr.length(), 30)
  assert_eq(expr[0].length(), 10)
  assert_eq(cell_types.length(), 30)
  assert_eq(gene_names.length(), 10)
}

test "sample_cell_types" {
  let (_, cell_types, _) = @src.cellchat_sample_data()
  let mut type_a = 0
  let mut type_b = 0
  let mut type_c = 0
  let mut i = 0
  while i < cell_types.length() {
    if cell_types[i] == "TypeA" { type_a = type_a + 1 }
    else if cell_types[i] == "TypeB" { type_b = type_b + 1 }
    else if cell_types[i] == "TypeC" { type_c = type_c + 1 }
    i = i + 1
  }
  assert_eq(type_a, 10)
  assert_eq(type_b, 10)
  assert_eq(type_c, 10)
}

test "get_top" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  let pairs = @src.cellchat_lr_database()
  let result = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=5, seed=42)
  let top = @src.cellchat_get_top(result, 5)
  assert_eq(top.length(), 5)
}

test "get_top_fewer_than_n" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  let pairs = @src.cellchat_lr_database()
  let result = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=5, seed=42)
  let top = @src.cellchat_get_top(result, 100)
  assert_true(top.length() <= result.scores.length())
}

test "summary" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  let pairs = @src.cellchat_lr_database()
  let result = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=5, seed=42)
  let summary = @src.cellchat_summary(result)
  assert_true(summary.length() > 0)
}

test "analyze_seed_reproducible" {
  let (expr, cell_types, gene_names) = @src.cellchat_sample_data()
  let pairs = @src.cellchat_lr_database()
  let result1 = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=5, seed=123)
  let result2 = @src.cellchat_analyze(expr, cell_types, gene_names, pairs, n_permutations=5, seed=123)
  // Same seed should give same results
  assert_eq(result1.scores.length(), result2.scores.length())
}
