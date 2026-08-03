///|
/// Tests for xCell module: immune/stromal cell type scoring using ssGSEA.
/// Tests cover signature creation, ssGSEA scoring, full xCell pipeline,
/// and result accessors.

// ============================================================================
// Signature tests
// ============================================================================

test "xc_signature_new_basic" {
  let sig = @src.XcellSignature::new(
    "CD8+ T cells",
    "immune",
    ["CD8A", "CD8B", "GZMA", "GZMB"],
  )
  assert_eq(sig.cell_type(), "CD8+ T cells")
  assert_eq(sig.category(), "immune")
  assert_eq(sig.genes().length(), 4)
  assert_eq(sig.genes()[0], "CD8A")
}

test "xc_default_signatures_count_and_categories" {
  let sigs = @src.xcell_default_signatures()
  // Should have 63 signatures
  assert_true(sigs.length() == 63)
  // Count categories
  let mut immune_count = 0
  let mut stromal_count = 0
  let mut other_count = 0
  for s in sigs {
    if s.category() == "immune" {
      immune_count = immune_count + 1
    } else if s.category() == "stromal" {
      stromal_count = stromal_count + 1
    } else {
      other_count = other_count + 1
    }
  }
  // At least 30 immune, 10 stromal, 8 other
  assert_true(immune_count >= 30)
  assert_true(stromal_count >= 10)
  assert_true(other_count >= 8)
}

test "xc_default_signatures_each_has_genes" {
  let sigs = @src.xcell_default_signatures()
  for s in sigs {
    assert_true(s.genes().length() >= 4)
    assert_true(s.cell_type().length() > 0)
  }
}

// ============================================================================
// Parameters
// ============================================================================

test "xc_params_default" {
  let p = @src.XcellParams::new()
  assert_eq(p.min_gene_overlap, 3)
  assert_true((p.alpha - 0.25).abs() < 0.0001)
  assert_eq(p.scale_scores, true)
  assert_eq(p.do_spillover, true)
  assert_eq(p.verbose, false)
}

// ============================================================================
// ssGSEA single-sample scoring
// ============================================================================

test "xc_ssgsea_single_basic" {
  let expr = [10.0, 5.0, 20.0, 3.0, 15.0]
  let genes = ["A", "B", "C", "D", "E"]
  // Set genes are all in data, top two highly expressed
  let set = ["A", "C", "E"]
  let score = @src.xcell_ssgsea_single(expr, genes, set, alpha=0.25)
  // Should produce a score between -1 and 1
  assert_true(score >= -1.0 && score <= 1.0)
}

test "xc_ssgsea_single_high_enrichment" {
  // Expression sorted: A=100, B=90, C=80, D=10, E=5 (high A,B,C top ranking)
  // Gene set = [A,B,C] are all top -> high positive enrichment
  let expr = [100.0, 90.0, 80.0, 10.0, 5.0]
  let genes = ["A", "B", "C", "D", "E"]
  let set_top = ["A", "B", "C"]
  let set_bottom = ["D", "E"]
  let score_high = @src.xcell_ssgsea_single(expr, genes, set_top, alpha=0.25)
  let score_low = @src.xcell_ssgsea_single(expr, genes, set_bottom, alpha=0.25)
  // Top ranked set should have higher score than bottom ranked set
  assert_true(score_high > score_low)
}

test "xc_ssgsea_single_empty" {
  let s = @src.xcell_ssgsea_single([], [], ["A"], alpha=0.25)
  assert_eq(s, 0.0)
  let s2 = @src.xcell_ssgsea_single([1.0], ["A"], [], alpha=0.25)
  assert_eq(s2, 0.0)
}

test "xc_ssgsea_single_min_overlap" {
  // Only 1 gene overlap, should return 0 when we need min 3 (but function returns 0 for <3)
  let expr = [1.0, 2.0]
  let genes = ["A", "B"]
  let set = ["A"]
  let s = @src.xcell_ssgsea_single(expr, genes, set, alpha=0.25)
  assert_eq(s, 0.0)
}

// ============================================================================
// Full xCell pipeline
// ============================================================================

test "xc_result_creation" {
  let scores : Array[Array[Double]] = [
    [0.5, 0.8],
    [0.2, 0.1],
  ]
  let cts = ["T cells", "B cells"]
  let cats = ["immune", "immune"]
  let samps = ["s1", "s2"]
  let imm = [0.35, 0.45]
  let str = [0.1, 0.15]
  let menv = [0.45, 0.6]
  let p = @src.XcellParams::new()
  let r = @src.XcellResult::new(
    scores, cts, cats, samps, imm, str, menv, p,
  )
  assert_eq(r.cell_types.length(), 2)
  assert_eq(r.sample_names.length(), 2)
  assert_true((r.immune_scores[0] - 0.35).abs() < 0.0001)
  assert_true((r.microenvironment_scores[1] - 0.6).abs() < 0.0001)
}

test "xc_run_default_simple" {
  // Simple expression matrix with a handful of marker genes
  let gene_names = [
    "CD3D", "CD3E", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7",
    "CD19", "MS4A1", "CD79A", "PECAM1", "VWF", "COL1A1", "FAP", "ALB",
  ]
  let sample_names = ["Tumor1", "Tumor2", "Normal1"]
  // 16 genes x 3 samples: upregulate immune in Normal1, stromal in Tumor2
  let expression = [
    [10.0, 8.0, 100.0],   // CD3D
    [12.0, 9.0, 110.0],   // CD3E
    [15.0, 10.0, 120.0],  // CD8A
    [14.0, 9.0, 115.0],   // CD8B
    [20.0, 15.0, 130.0],  // GZMA
    [18.0, 12.0, 125.0],  // GZMB
    [16.0, 11.0, 110.0],  // PRF1
    [17.0, 13.0, 115.0],  // NKG7
    [5.0, 4.0, 80.0],     // CD19
    [6.0, 5.0, 85.0],     // MS4A1
    [5.0, 3.0, 75.0],     // CD79A
    [20.0, 100.0, 15.0],  // PECAM1
    [18.0, 95.0, 12.0],   // VWF
    [30.0, 150.0, 10.0],  // COL1A1
    [25.0, 140.0, 8.0],   // FAP
    [10.0, 20.0, 200.0],  // ALB
  ]
  let result = @src.xcell_run_default(expression, gene_names, sample_names)
  // Should have returned scores with 63 cell types and 3 samples
  assert_eq(result.scores.length(), 63)
  assert_eq(result.scores[0].length(), 3)
  assert_eq(result.immune_scores.length(), 3)
  assert_eq(result.stromal_scores.length(), 3)
  // Normal1 should have high immune scores (since all immune markers are high)
  // Just check not NaN
  for s in 0..<3 {
    assert_false(result.immune_scores[s].is_nan())
    assert_false(result.stromal_scores[s].is_nan())
  }
}

// ============================================================================
// Result accessors
// ============================================================================

test "xc_result_get_cell_type_scores" {
  let scores : Array[Array[Double]] = [
    [0.3, 0.5],
    [0.8, 0.6],
  ]
  let cts = ["T-cells", "B-cells"]
  let cats = ["immune", "immune"]
  let samps = ["s1", "s2"]
  let imm = [0.55, 0.55]
  let str = [0.1, 0.1]
  let menv = [0.65, 0.65]
  let p = @src.XcellParams::new()
  let r = @src.XcellResult::new(
    scores, cts, cats, samps, imm, str, menv, p,
  )
  let tc = r.get_cell_type_scores("T-cells")
  assert_eq(tc.length(), 2)
  assert_true((tc[0] - 0.3).abs() < 0.0001)
  assert_true((tc[1] - 0.5).abs() < 0.0001)
  let missing = r.get_cell_type_scores("NK-cells")
  assert_eq(missing.length(), 0)
}

test "xc_result_get_sample_scores" {
  let scores : Array[Array[Double]] = [
    [0.3, 0.5],
    [0.8, 0.6],
  ]
  let cts = ["T-cells", "B-cells"]
  let cats = ["immune", "immune"]
  let samps = ["s1", "s2"]
  let imm = [0.55, 0.55]
  let str = [0.1, 0.1]
  let menv = [0.65, 0.65]
  let p = @src.XcellParams::new()
  let r = @src.XcellResult::new(
    scores, cts, cats, samps, imm, str, menv, p,
  )
  let s1 = r.get_sample_scores("s1")
  assert_eq(s1.length(), 2)
  assert_true((s1[0] - 0.3).abs() < 0.0001)
  assert_true((s1[1] - 0.8).abs() < 0.0001)
  let s_missing = r.get_sample_scores("missing")
  assert_eq(s_missing.length(), 0)
}

test "xc_result_get_top_cell_types" {
  let scores : Array[Array[Double]] = [
    [0.1, 0.9],
    [0.5, 0.7],
    [0.9, 0.2],
    [0.3, 0.5],
  ]
  let cts = ["A", "B", "C", "D"]
  let cats = ["immune", "stromal", "immune", "other"]
  let samps = ["s1", "s2"]
  let imm = [0.5, 0.55]
  let str = [0.1, 0.1]
  let menv = [0.6, 0.65]
  let p = @src.XcellParams::new()
  let r = @src.XcellResult::new(
    scores, cts, cats, samps, imm, str, menv, p,
  )
  // For s1, top is C=0.9, then B=0.5
  let top_s1 = r.get_top_cell_types("s1", 2)
  assert_eq(top_s1.length(), 2)
  assert_eq(top_s1[0].0, "C")
  assert_true((top_s1[0].1 - 0.9).abs() < 0.0001)
  assert_eq(top_s1[1].0, "B")
  assert_true((top_s1[1].1 - 0.5).abs() < 0.0001)
  // For s2, top is A=0.9
  let top_s2 = r.get_top_cell_types("s2", 1)
  assert_eq(top_s2.length(), 1)
  assert_eq(top_s2[0].0, "A")
}

test "xc_result_scores_by_category" {
  let scores : Array[Array[Double]] = [
    [0.5, 0.7],
    [0.3, 0.4],
    [0.9, 0.8],
  ]
  let cts = ["T", "B", "Fibro"]
  let cats = ["immune", "immune", "stromal"]
  let samps = ["s1", "s2"]
  let imm = [0.4, 0.55]
  let str = [0.9, 0.8]
  let menv = [1.3, 1.35]
  let p = @src.XcellParams::new()
  let r = @src.XcellResult::new(
    scores, cts, cats, samps, imm, str, menv, p,
  )
  let by_cat = r.scores_by_category()
  // "immune" average of T and B for s1: (0.5+0.3)/2 = 0.4
  let imm_opt = by_cat.get("immune")
  let imm_arr = match imm_opt {
    Some(a) => a
    None => [0.0, 0.0]
  }
  assert_true((imm_arr[0] - 0.4).abs() < 0.001)
  // "stromal" for s1: 0.9
  let str_opt = by_cat.get("stromal")
  let str_arr = match str_opt {
    Some(a) => a
    None => [0.0, 0.0]
  }
  assert_true((str_arr[0] - 0.9).abs() < 0.001)
}
