///|
/// Tests for Bioconductor ESTIMATE module - tumor microenvironment scoring
/// via stromal/immune gene set enrichment.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "est_expression_creation" {
  let expr = @src.EstExpression::new(
    ["G1", "G2"],
    ["S1", "S2"],
    [[1.0, 2.0], [3.0, 4.0]],
  )
  assert_eq(expr.gene_names().length(), 2)
  assert_eq(expr.sample_names().length(), 2)
  assert_eq(expr.matrix()[1][0], 3.0)
}

test "est_signature_creation" {
  let sig = @src.EstSignature::new(
    ["STROMA1", "STROMA2"],
    ["IMMUNE1", "IMMUNE2"],
  )
  assert_eq(sig.stromal_genes().length(), 2)
  assert_eq(sig.immune_genes().length(), 2)
  assert_eq(sig.stromal_genes()[0], "STROMA1")
  assert_eq(sig.immune_genes()[1], "IMMUNE2")
}

test "est_sample_score_accessors" {
  // Use est_run to get real scores and test accessors
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  let s = result.scores()[0]
  assert_eq(s.sample_id(), "HighStromal")
  // Stromal and immune scores should be in [-100, 100]
  assert_true(s.stromal_score() >= -100.1 && s.stromal_score() <= 100.1)
  assert_true(s.immune_score() >= -100.1 && s.immune_score() <= 100.1)
  // ESTIMATE score = stromal + immune
  assert_true(
    (s.estimate_score() - (s.stromal_score() + s.immune_score())).abs() < 0.001,
  )
}

test "est_result_accessors" {
  // Use est_run to get a real result and test accessors
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  assert_eq(result.scores().length(), 4)
  assert_eq(result.sample_names().length(), 4)
  assert_eq(result.method_name(), "estimate")
}

// ---------------------------------------------------------------------------
// Built-in gene signatures
// ---------------------------------------------------------------------------

test "est_stromal_genes_nonempty" {
  let genes = @src.est_stromal_genes()
  assert_true(genes.length() >= 30)
  // Should contain typical stromal markers like COL1A1, DCN
  let mut has_col1a1 = false
  let mut has_dcn = false
  for g in genes {
    if g == "COL1A1" { has_col1a1 = true }
    if g == "DCN" { has_dcn = true }
  }
  assert_true(has_col1a1)
  assert_true(has_dcn)
}

test "est_immune_genes_nonempty" {
  let genes = @src.est_immune_genes()
  assert_true(genes.length() >= 20)
  // Should contain typical immune markers like CD3D, CD3E
  let mut has_cd3d = false
  for g in genes {
    if g == "CD3D" { has_cd3d = true }
  }
  assert_true(has_cd3d)
}

test "est_default_signature" {
  let sig = @src.est_default_signature()
  assert_true(sig.stromal_genes().length() >= 30)
  assert_true(sig.immune_genes().length() >= 20)
}

// ---------------------------------------------------------------------------
// ECDF normalization behavior
// ---------------------------------------------------------------------------

test "est_run_normalizes_across_samples" {
  // With 4 samples, ECDF values are 0.25, 0.5, 0.75, 1.0
  // Scores should be in [-100, 100] range
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  assert_eq(result.scores().length(), 4)
  for s in result.scores() {
    assert_true(s.stromal_score() >= -100.1)
    assert_true(s.stromal_score() <= 100.1)
    assert_true(s.immune_score() >= -100.1)
    assert_true(s.immune_score() <= 100.1)
    // ESTIMATE score = stromal + immune, in [-200, 200]
    assert_true(s.estimate_score() >= -200.1)
    assert_true(s.estimate_score() <= 200.1)
  }
}

// ---------------------------------------------------------------------------
// Sample data correctness
// ---------------------------------------------------------------------------

test "est_sample_data_dimensions" {
  let expr = @src.est_sample_data()
  assert_eq(expr.sample_names().length(), 4)
  assert_true(expr.gene_names().length() >= 50) // 50 stromal + 30 immune (deduped)
  // Matrix rows = genes, cols = samples
  assert_eq(expr.matrix().length(), expr.gene_names().length())
  for row in expr.matrix() {
    assert_eq(row.length(), 4)
  }
}

test "est_sample_data_sample_names" {
  let expr = @src.est_sample_data()
  assert_eq(expr.sample_names()[0], "HighStromal")
  assert_eq(expr.sample_names()[1], "HighImmune")
  assert_eq(expr.sample_names()[2], "HighBoth")
  assert_eq(expr.sample_names()[3], "LowBoth")
}

// ---------------------------------------------------------------------------
// Score semantics
// ---------------------------------------------------------------------------

test "est_run_high_stromal_higher_than_low_stromal" {
  // HighStromal sample should have higher stromal score than LowBoth
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  let high_stromal = result.get_score("HighStromal")
  let low_both = result.get_score("LowBoth")
  assert_true(high_stromal.stromal_score() > low_both.stromal_score())
}

test "est_run_high_immune_higher_than_low_immune" {
  // HighImmune sample should have higher immune score than LowBoth
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  let high_immune = result.get_score("HighImmune")
  let low_both = result.get_score("LowBoth")
  assert_true(high_immune.immune_score() > low_both.immune_score())
}

test "est_run_high_both_highest_estimate_score" {
  // HighBoth should have the highest combined ESTIMATE score
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  let high_both = result.get_score("HighBoth")
  let low_both = result.get_score("LowBoth")
  assert_true(high_both.estimate_score() > low_both.estimate_score())
}

test "est_run_estimate_score_equals_stromal_plus_immune" {
  // ESTIMATEScore = StromalScore + ImmuneScore
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  for s in result.scores() {
    let sum = s.stromal_score() + s.immune_score()
    assert_true((s.estimate_score() - sum).abs() < 0.001)
  }
}

test "est_run_low_both_has_low_scores" {
  // LowBoth (low infiltration) should have negative scores for both
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  let low_both = result.get_score("LowBoth")
  assert_true(low_both.stromal_score() < 0.0)
  assert_true(low_both.immune_score() < 0.0)
}

// ---------------------------------------------------------------------------
// Tumor purity
// ---------------------------------------------------------------------------

test "est_run_tumor_purity_in_valid_range" {
  // Tumor purity should be in [0, 1] or NaN
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  for s in result.scores() {
    let p = s.tumor_purity()
    if p.is_nan() {
      // NaN is acceptable (when arg is out of cosine domain)
      assert_true(p.is_nan())
    } else {
      assert_true(p >= -0.001)
      assert_true(p <= 1.001)
    }
  }
}

test "est_run_low_infiltration_higher_purity" {
  // LowBoth (low infiltration -> higher tumor content) should have higher
  // tumor purity than HighBoth (high infiltration -> lower tumor content)
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  let low_both = result.get_score("LowBoth")
  let high_both = result.get_score("HighBoth")
  if !low_both.tumor_purity().is_nan() &&
    !high_both.tumor_purity().is_nan() {
    assert_true(low_both.tumor_purity() >= high_both.tumor_purity())
  }
}

// ---------------------------------------------------------------------------
// Custom signature
// ---------------------------------------------------------------------------

test "est_run_custom_signature" {
  // Build expression matrix with 4 genes, 2 samples
  // G1, G2 are stromal markers; G3, G4 are immune markers
  let expr = @src.EstExpression::new(
    ["G1", "G2", "G3", "G4"],
    ["S1", "S2"],
    [
      [10.0, 1.0], // G1 - stromal, high in S1
      [8.0, 1.0], // G2 - stromal, high in S1
      [1.0, 10.0], // G3 - immune, high in S2
      [1.0, 8.0], // G4 - immune, high in S2
    ],
  )
  let sig = @src.EstSignature::new(["G1", "G2"], ["G3", "G4"])
  let result = @src.est_run(expr, signature=sig)
  assert_eq(result.scores().length(), 2)
  let s1 = result.get_score("S1")
  let s2 = result.get_score("S2")
  // S1 should have higher stromal score; S2 higher immune score
  assert_true(s1.stromal_score() > s2.stromal_score())
  assert_true(s2.immune_score() > s1.immune_score())
}

test "est_run_empty_signature_returns_zeros" {
  // With empty gene sets, ssGSEA returns 0; ECDF of all-zeros is degenerate
  let expr = @src.EstExpression::new(
    ["G1", "G2"],
    ["S1", "S2"],
    [[1.0, 2.0], [3.0, 4.0]],
  )
  let sig = @src.EstSignature::new([], [])
  let result = @src.est_run(expr, signature=sig)
  for s in result.scores() {
    // All scores should be 0 (or very close to 0) since raw scores are 0
    assert_true(s.stromal_score().abs() < 1.0)
    assert_true(s.immune_score().abs() < 1.0)
  }
}

// ---------------------------------------------------------------------------
// get_score and to_string
// ---------------------------------------------------------------------------

test "est_result_get_score_known_sample" {
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  let s = result.get_score("HighStromal")
  assert_eq(s.sample_id(), "HighStromal")
}

test "est_result_get_score_unknown_sample" {
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  let s = result.get_score("Nonexistent")
  // Unknown sample returns a default score
  assert_eq(s.sample_id(), "Nonexistent")
  assert_true(s.stromal_score().abs() < 1.0e-9)
  assert_true(s.estimate_score().abs() < 1.0e-9)
}

test "est_result_to_string" {
  let expr = @src.est_sample_data()
  let result = @src.est_run(expr)
  let s = result.to_string()
  assert_true(s.contains("ESTIMATE"))
  assert_true(s.contains("StromalScore"))
  assert_true(s.contains("ImmuneScore"))
  assert_true(s.contains("HighStromal"))
  assert_true(s.contains("TumorPurity"))
}

// ---------------------------------------------------------------------------
// Edge cases
// ---------------------------------------------------------------------------

test "est_run_single_sample" {
  // Single sample: ECDF is degenerate (trivially all-equal), scores are 0.0
  let expr = @src.EstExpression::new(
    ["G1", "G2", "G3", "G4"],
    ["Only"],
    [[10.0], [8.0], [5.0], [3.0]],
  )
  let sig = @src.EstSignature::new(["G1", "G2"], ["G3", "G4"])
  let result = @src.est_run(expr, signature=sig)
  assert_eq(result.scores().length(), 1)
  // With only one sample, degenerate case: scores are 0 (neutral)
  let s = result.get_score("Only")
  assert_true(s.stromal_score().abs() < 1.0e-9)
  assert_true(s.immune_score().abs() < 1.0e-9)
}

test "est_run_no_signature_genes_present" {
  // None of the signature genes are in the expression matrix
  let expr = @src.EstExpression::new(
    ["X1", "X2", "X3"],
    ["S1", "S2"],
    [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
  )
  let sig = @src.EstSignature::new(["G1", "G2"], ["G3", "G4"])
  let result = @src.est_run(expr, signature=sig)
  // All scores should be 0 (or near 0) since no genes match
  for s in result.scores() {
    assert_true(s.stromal_score().abs() < 1.0)
    assert_true(s.immune_score().abs() < 1.0)
  }
}

test "est_run_all_samples_identical" {
  // Identical samples should have identical scores
  let expr = @src.EstExpression::new(
    ["G1", "G2", "G3", "G4"],
    ["S1", "S2", "S3"],
    [[5.0, 5.0, 5.0], [5.0, 5.0, 5.0], [5.0, 5.0, 5.0], [5.0, 5.0, 5.0]],
  )
  let sig = @src.EstSignature::new(["G1", "G2"], ["G3", "G4"])
  let result = @src.est_run(expr, signature=sig)
  let s1 = result.get_score("S1")
  let s2 = result.get_score("S2")
  let s3 = result.get_score("S3")
  // All scores should be equal
  assert_true((s1.stromal_score() - s2.stromal_score()).abs() < 1.0e-9)
  assert_true((s1.immune_score() - s3.immune_score()).abs() < 1.0e-9)
}
