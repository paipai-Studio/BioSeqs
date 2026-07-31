///|
/// Tests for Bioconductor decoupleR module - Functional activity inference.

test "pkn_edge_create" {
  let e = @src.pkn_edge("TF1", "G1", 1.0)
  assert_eq(e.source, "TF1")
  assert_eq(e.target, "G1")
  assert_true((e.weight - 1.0).abs() < 1.0e-9)
}

test "pkn_from_edges" {
  let edges = [
    @src.pkn_edge("TF1", "G1", 1.0),
    @src.pkn_edge("TF1", "G2", 1.0),
    @src.pkn_edge("TF2", "G3", -1.0),
  ]
  let pkn = @src.pkn_from_edges(edges)
  assert_eq(pkn.edges.length(), 3)
  assert_eq(pkn.regulators.length(), 2)
  assert_eq(pkn.targets.length(), 3)
}

test "decoupler_sample_data_shape" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  assert_eq(samples.length(), 5)
  assert_eq(genes.length(), 11)
  assert_eq(expr.length(), 5)
  assert_eq(expr[0].length(), 11)
  assert_true(pkn.edges.length() > 0)
}

test "decoupler_wsum_basic" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  assert_eq(result.method, "wsum")
  assert_eq(result.samples.length(), 5)
  // 3 TFs in the PKN
  assert_eq(result.regulators.length(), 3)
  // TF1 should have higher activity in S1/S3 than in S4
  let mut tf1_idx = -1
  let mut i = 0
  while i < result.regulators.length() {
    if result.regulators[i] == "TF1" {
      tf1_idx = i
    }
    i = i + 1
  }
  assert_true(tf1_idx >= 0)
  // S1 (idx 0) is TF1-high, S4 (idx 3) is TF2-high
  assert_true(result.matrix[tf1_idx][0] > result.matrix[tf1_idx][3])
}

test "decoupler_wmean_basic" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wmean_method())
  assert_eq(result.method, "wmean")
  assert_eq(result.regulators.length(), 3)
  // wmean normalizes by sum of absolute weights, so magnitude is smaller than wsum
  let wsum_result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  let mut i = 0
  while i < result.regulators.length() {
    let mut j = 0
    while j < result.samples.length() {
      // |wmean| should be <= |wsum| (since abs_w_sum >= 1.0)
      assert_true(result.matrix[i][j].abs() <= wsum_result.matrix[i][j].abs() + 1.0e-9)
      j = j + 1
    }
    i = i + 1
  }
}

test "decoupler_norm_basic" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_norm_method())
  assert_eq(result.method, "norm")
  // Normalized scores should have mean ~0 per regulator across samples
  let mut i = 0
  while i < result.regulators.length() {
    let mut mean = 0.0
    let mut j = 0
    while j < result.samples.length() {
      mean = mean + result.matrix[i][j]
      j = j + 1
    }
    mean = mean / result.samples.length().to_double()
    assert_true(mean.abs() < 1.0e-9)
    i = i + 1
  }
}

test "decoupler_ulm_basic" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_ulm_method())
  assert_eq(result.method, "ulm")
  assert_eq(result.regulators.length(), 3)
  // ULM scores are standardized; mean across samples should be ~0
  let mut i = 0
  while i < result.regulators.length() {
    let mut mean = 0.0
    let mut j = 0
    while j < result.samples.length() {
      mean = mean + result.matrix[i][j]
      j = j + 1
    }
    mean = mean / result.samples.length().to_double()
    assert_true(mean.abs() < 1.0e-6)
    i = i + 1
  }
}

test "decoupler_mlm_basic" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_mlm_method())
  assert_eq(result.method, "mlm")
  assert_eq(result.regulators.length(), 3)
  // MLM slope for TF1 should be positive (TF1 activates its targets)
  let mut tf1_idx = -1
  let mut i = 0
  while i < result.regulators.length() {
    if result.regulators[i] == "TF1" {
      tf1_idx = i
    }
    i = i + 1
  }
  assert_true(tf1_idx >= 0)
  // In S1, TF1 is high → mlm slope (which is roughly y/x) should be positive
  assert_true(result.matrix[tf1_idx][0] > 0.0)
}

test "decoupler_scores_per_cell" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  // Should have one score per regulator × sample combination
  assert_eq(result.scores.length(), result.regulators.length() * result.samples.length())
}

test "decoupler_top_regulators" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  let top = @src.decoupler_top_regulators(result, "S1", 2)
  assert_true(top.length() <= 2)
  // Top regulators for S1 should include TF1 (since S1 is TF1-high)
  let mut found_tf1 = false
  let mut i = 0
  while i < top.length() {
    if top[i].regulator == "TF1" {
      found_tf1 = true
    }
    i = i + 1
  }
  assert_true(found_tf1)
}

test "decoupler_top_regulators_unknown_sample" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  let top = @src.decoupler_top_regulators(result, "UNKNOWN", 3)
  assert_eq(top.length(), 0)
}

test "decoupler_filter_scores" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  let filtered = @src.decoupler_filter_scores(result, 0.5)
  // All filtered scores should have |score| >= 0.5
  let mut i = 0
  while i < filtered.length() {
    assert_true(filtered[i].score.abs() >= 0.5)
    i = i + 1
  }
}

test "decoupler_method_to_string" {
  assert_eq(@src.decoupler_wsum_method().to_string(), "wsum")
  assert_eq(@src.decoupler_wmean_method().to_string(), "wmean")
  assert_eq(@src.decoupler_norm_method().to_string(), "norm")
  assert_eq(@src.decoupler_ulm_method().to_string(), "ulm")
  assert_eq(@src.decoupler_mlm_method().to_string(), "mlm")
}

test "decoupler_default_method_is_ulm" {
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  // Default method argument should be ULM
  let result = @src.decoupler_run(expr, samples, genes, pkn)
  assert_eq(result.method, "ulm")
}

test "decoupler_wsum_symmetric_in_sign" {
  // A PKN with negative weight should produce negative wsum
  let edges = [
    @src.pkn_edge("REP", "G1", -1.0),
  ]
  let pkn = @src.pkn_from_edges(edges)
  let expr = [[0.0, 5.0, 0.0], [0.0, 0.0, 0.0]]
  let samples = ["S1", "S2"]
  let genes = ["REP", "G1", "G2"]
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  // REP's only target is G1 with weight -1.0; in S1, G1=5 → wsum=-5; in S2, G1=0 → wsum=0
  assert_eq(result.regulators.length(), 1)
  assert_true((result.matrix[0][0] - -5.0).abs() < 1.0e-9)
  assert_true((result.matrix[0][1] - 0.0).abs() < 1.0e-9)
}

test "decoupler_pkn_with_unknown_genes" {
  // PKN with regulator/target not in expression matrix should be filtered out
  let edges = [
    @src.pkn_edge("TF1", "G1", 1.0),
    @src.pkn_edge("UNKNOWN_TF", "G1", 1.0),
    @src.pkn_edge("TF1", "UNKNOWN_TARGET", 1.0),
  ]
  let pkn = @src.pkn_from_edges(edges)
  let expr = [[1.0, 2.0], [3.0, 4.0]]
  let samples = ["S1", "S2"]
  let genes = ["TF1", "G1"]
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  // Only TF1 (which appears as a gene) should be in regulators
  assert_eq(result.regulators.length(), 1)
  assert_eq(result.regulators[0], "TF1")
}

test "decoupler_empty_pkn" {
  let pkn = @src.pkn_from_edges([])
  let expr = [[1.0, 2.0], [3.0, 4.0]]
  let samples = ["S1", "S2"]
  let genes = ["TF1", "G1"]
  let result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  assert_eq(result.regulators.length(), 0)
  assert_eq(result.scores.length(), 0)
}
