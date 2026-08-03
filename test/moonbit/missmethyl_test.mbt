///|
/// Tests for the missMethyl module.
/// Covers probe annotations, DMP analysis, GO enrichment, bias correction,
/// value conversions, statistical tests, filtering, summaries, and edge cases.

// ===========================================================================
// MmProbeAnnotation creation and accessors
// ===========================================================================

///|
test "mm_probe_annotation_creation" {
  let ann = @src.MmProbeAnnotation::new(
    "cg001", "chr1", 1000, "TP53", "Island", "I",
  )
  assert_eq(ann.probe_id(), "cg001")
  assert_eq(ann.chr(), "chr1")
  assert_eq(ann.pos(), 1000)
  assert_eq(ann.gene(), "TP53")
  assert_eq(ann.island_relation(), "Island")
  assert_eq(ann.probe_type(), "I")
}

///|
test "mm_probe_annotation_creation_type_II" {
  let ann = @src.MmProbeAnnotation::new(
    "cg999", "chrX", 99999, "BRCA1", "OpenSea", "II",
  )
  assert_eq(ann.probe_id(), "cg999")
  assert_eq(ann.chr(), "chrX")
  assert_eq(ann.pos(), 99999)
  assert_eq(ann.gene(), "BRCA1")
  assert_eq(ann.island_relation(), "OpenSea")
  assert_eq(ann.probe_type(), "II")
}

// ===========================================================================
// MmDmpResult creation, accessors, is_significant
// ===========================================================================

///|
test "mm_dmp_result_creation" {
  let r = @src.MmDmpResult::new(
    "cg001", 1.5, 3.2, 0.001, 0.005, "chr1", 1000, "TP53",
  )
  assert_eq(r.probe_id(), "cg001")
  assert_true((r.log_fc() - 1.5).abs() < 0.001)
  assert_true((r.t_stat() - 3.2).abs() < 0.001)
  assert_true((r.p_value() - 0.001).abs() < 0.001)
  assert_true((r.adj_p_value() - 0.005).abs() < 0.001)
  assert_eq(r.chr(), "chr1")
  assert_eq(r.pos(), 1000)
  assert_eq(r.gene(), "TP53")
}

///|
test "mm_dmp_result_is_significant" {
  let sig = @src.MmDmpResult::new(
    "cg001", 2.0, 5.0, 0.0001, 0.005, "chr1", 1000, "A",
  )
  let nonsig = @src.MmDmpResult::new(
    "cg002", 0.1, 0.5, 0.4, 0.6, "chr1", 2000, "B",
  )
  assert_true(sig.is_significant(0.05))
  assert_true(sig.is_significant(0.01))
  assert_true(!sig.is_significant(0.001))
  assert_true(!nonsig.is_significant(0.05))
  assert_true(!nonsig.is_significant(0.01))
}

// ===========================================================================
// MmGoDatabase creation, add_entry, genes_for_go, go_for_gene
// ===========================================================================

///|
test "mm_go_database_creation_and_add" {
  let db = @src.MmGoDatabase::new()
  assert_eq(db.genes_for_go("GO:0000001").length(), 0)
  db.add_entry("GO:0000001", "TP53", "cell cycle", "BP")
  db.add_entry("GO:0000001", "MDM2", "cell cycle", "BP")
  db.add_entry("GO:0000002", "BAX", "apoptosis", "BP")
  let genes = db.genes_for_go("GO:0000001")
  assert_eq(genes.length(), 2)
  assert_true(genes.contains("TP53"))
  assert_true(genes.contains("MDM2"))
}

///|
test "mm_go_database_go_for_gene" {
  let db = @src.MmGoDatabase::new()
  db.add_entry("GO:0000001", "TP53", "cell cycle", "BP")
  db.add_entry("GO:0000002", "TP53", "apoptosis", "BP")
  db.add_entry("GO:0000003", "BAX", "nucleus", "CC")
  let go_ids = db.go_for_gene("TP53")
  assert_eq(go_ids.length(), 2)
  assert_true(go_ids.contains("GO:0000001"))
  assert_true(go_ids.contains("GO:0000002"))
  let go_bax = db.go_for_gene("BAX")
  assert_eq(go_bax.length(), 1)
  assert_eq(go_bax[0], "GO:0000003")
}

///|
test "mm_go_database_dedup" {
  let db = @src.MmGoDatabase::new()
  db.add_entry("GO:0000001", "A", "term", "BP")
  db.add_entry("GO:0000001", "A", "term", "BP")
  let genes = db.genes_for_go("GO:0000001")
  assert_eq(genes.length(), 1)
}

// ===========================================================================
// mm_beta_to_m and mm_m_to_beta
// ===========================================================================

///|
test "mm_beta_to_m_basic" {
  // beta = 0.5 -> M = log2(0.5/0.5) = log2(1) = 0
  let m = @src.mm_beta_to_m(0.5)
  assert_true(m.abs() < 0.01)
  // beta = 0.25 -> M = log2(0.25/0.75) = log2(1/3) ≈ -1.585
  let m2 = @src.mm_beta_to_m(0.25)
  assert_true((m2 - (-1.585)).abs() < 0.01)
  // beta = 0.75 -> M = log2(0.75/0.25) = log2(3) ≈ 1.585
  let m3 = @src.mm_beta_to_m(0.75)
  assert_true((m3 - 1.585).abs() < 0.01)
}

///|
test "mm_beta_to_m_edge_cases" {
  // beta = 0.0 -> clamped to 1e-6, M ≈ log2(1e-6) ≈ -19.93
  let m0 = @src.mm_beta_to_m(0.0)
  assert_true(m0 < -19.0)
  // beta = 1.0 -> clamped to 1-1e-6, M ≈ log2(1e6) ≈ 19.93
  let m1 = @src.mm_beta_to_m(1.0)
  assert_true(m1 > 19.0)
  // Clamping is symmetric: M(0) ≈ -M(1)
  assert_true((m0 + m1).abs() < 0.1)
}

///|
test "mm_m_to_beta_basic" {
  // M = 0 -> beta = 2^0 / (1 + 2^0) = 1/2 = 0.5
  let b = @src.mm_m_to_beta(0.0)
  assert_true((b - 0.5).abs() < 0.01)
  // M = 1 -> beta = 2 / (1 + 2) = 2/3 ≈ 0.667
  let b2 = @src.mm_m_to_beta(1.0)
  assert_true((b2 - 0.667).abs() < 0.01)
  // M = -1 -> beta = 0.5 / (1 + 0.5) = 1/3 ≈ 0.333
  let b3 = @src.mm_m_to_beta(-1.0)
  assert_true((b3 - 0.333).abs() < 0.01)
}

///|
test "mm_beta_m_round_trip" {
  let betas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
  for b in betas {
    let m = @src.mm_beta_to_m(b)
    let b_back = @src.mm_m_to_beta(m)
    assert_true((b - b_back).abs() < 0.001)
  }
}

// ===========================================================================
// mm_t_test
// ===========================================================================

///|
test "mm_t_test_same_group" {
  let g1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let g2 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let (t, p) = @src.mm_t_test(g1, g2)
  assert_true(t.abs() < 0.01)
  assert_true(p > 0.99)
}

///|
test "mm_t_test_different_groups" {
  let g1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let g2 = [10.0, 11.0, 12.0, 13.0, 14.0]
  let (t, p) = @src.mm_t_test(g1, g2)
  assert_true(t < -5.0)
  assert_true(p < 0.01)
}

// ===========================================================================
// mm_bh_correct
// ===========================================================================

///|
test "mm_bh_correct_basic" {
  let pvalues = [0.01, 0.5, 0.001, 0.2]
  let adj = @src.mm_bh_correct(pvalues)
  assert_eq(adj.length(), 4)
  // All adjusted values are valid probabilities.
  for v in adj {
    assert_true(v >= 0.0 && v <= 1.0)
  }
  // The smallest p (index 2) gets the smallest adjusted value.
  assert_true(adj[2] <= adj[0])
  assert_true(adj[0] <= adj[3])
  assert_true(adj[3] <= adj[1])
}

///|
test "mm_bh_correct_all_ones" {
  let pvalues = [1.0, 1.0, 1.0]
  let adj = @src.mm_bh_correct(pvalues)
  assert_eq(adj.length(), 3)
  for v in adj {
    assert_true((v - 1.0).abs() < 0.001)
  }
}

///|
test "mm_bh_correct_empty" {
  let adj = @src.mm_bh_correct([])
  assert_eq(adj.length(), 0)
}

// ===========================================================================
// mm_fisher_exact
// ===========================================================================

///|
test "mm_fisher_exact_identical_table" {
  // [[5,5],[5,5]] -> no association, p should be 1.0
  let p = @src.mm_fisher_exact(5, 5, 5, 5)
  assert_true(p > 0.99)
}

///|
test "mm_fisher_exact_extreme_table" {
  // [[10,0],[0,10]] -> perfect association, p very small
  let p = @src.mm_fisher_exact(10, 0, 0, 10)
  assert_true(p < 0.001)
}

///|
test "mm_fisher_exact_known_result" {
  // [[3,1],[1,3]] -> two-sided p ≈ 0.4857
  let p = @src.mm_fisher_exact(3, 1, 1, 3)
  assert_true(p > 0.3 && p < 0.7)
}

// ===========================================================================
// mm_map_probes_to_genes
// ===========================================================================

///|
test "mm_map_probes_to_genes" {
  let annotations : Array[@src.MmProbeAnnotation] = [
    @src.MmProbeAnnotation::new("cg001", "chr1", 100, "A", "Island", "I"),
    @src.MmProbeAnnotation::new("cg002", "chr1", 200, "B", "Shore", "II"),
    @src.MmProbeAnnotation::new("cg003", "chr1", 300, "A", "OpenSea", "I"),
    @src.MmProbeAnnotation::new("cg004", "chr1", 400, "C", "Island", "II"),
  ]
  let genes = @src.mm_map_probes_to_genes(
    ["cg001", "cg003", "cg004"],
    annotations,
  )
  assert_eq(genes.length(), 2)
  assert_eq(genes[0], "A")
  assert_eq(genes[1], "C")
}

// ===========================================================================
// mm_filter_probes
// ===========================================================================

///|
test "mm_filter_probes_by_chr" {
  let annotations = @src.mm_sample_annotations()
  let chr1 = @src.mm_filter_probes(annotations, "chr1", "all")
  assert_eq(chr1.length(), 7)
  for ann in chr1 {
    assert_eq(ann.chr(), "chr1")
  }
  let chr2 = @src.mm_filter_probes(annotations, "chr2", "all")
  assert_eq(chr2.length(), 7)
  let chr3 = @src.mm_filter_probes(annotations, "chr3", "all")
  assert_eq(chr3.length(), 6)
}

///|
test "mm_filter_probes_by_island" {
  let annotations = @src.mm_sample_annotations()
  let islands = @src.mm_filter_probes(annotations, "all", "Island")
  for ann in islands {
    assert_eq(ann.island_relation(), "Island")
  }
  assert_true(islands.length() > 0)
  let shores = @src.mm_filter_probes(annotations, "all", "Shore")
  for ann in shores {
    assert_eq(ann.island_relation(), "Shore")
  }
  assert_true(shores.length() > 0)
}

///|
test "mm_filter_probes_combined" {
  let annotations = @src.mm_sample_annotations()
  let chr1_islands = @src.mm_filter_probes(annotations, "chr1", "Island")
  for ann in chr1_islands {
    assert_eq(ann.chr(), "chr1")
    assert_eq(ann.island_relation(), "Island")
  }
  // chr1 has probes cg001(Island) and cg005(Island)
  assert_eq(chr1_islands.length(), 2)
}

///|
test "mm_filter_probes_empty_filter" {
  let annotations = @src.mm_sample_annotations()
  let all = @src.mm_filter_probes(annotations, "", "")
  assert_eq(all.length(), 20)
}

// ===========================================================================
// mm_dmp_analysis
// ===========================================================================

///|
test "mm_dmp_analysis_sample" {
  let beta = @src.mm_sample_beta()
  let annotations = @src.mm_sample_annotations()
  let group = [false, false, false, false, true, true, true, true]
  let results = @src.mm_dmp_analysis(beta, group, annotations)
  assert_eq(results.length(), 20)
  // Each result should have correct annotation fields.
  assert_eq(results[0].probe_id(), "cg001")
  assert_eq(results[0].chr(), "chr1")
  assert_eq(results[0].gene(), "A")
}

///|
test "mm_dmp_analysis_significant" {
  let beta = @src.mm_sample_beta()
  let annotations = @src.mm_sample_annotations()
  let group = [false, false, false, false, true, true, true, true]
  let results = @src.mm_dmp_analysis(beta, group, annotations)
  let mut n_sig = 0
  for r in results {
    if r.is_significant(0.05) {
      n_sig = n_sig + 1
    }
  }
  // The first 5 probes are strongly differential.
  assert_true(n_sig >= 5)
  // The first probe (cg001) should be significant and hypermethylated.
  assert_true(results[0].is_significant(0.05))
  assert_true(results[0].log_fc() > 0.0)
}

// ===========================================================================
// mm_go_enrichment
// ===========================================================================

///|
test "mm_go_enrichment_sample" {
  let annotations = @src.mm_sample_annotations()
  let go_db = @src.mm_sample_go_database()
  let all_probes : Array[String] = [
    "cg001", "cg002", "cg003", "cg004", "cg005",
    "cg006", "cg007", "cg008", "cg009", "cg010",
    "cg011", "cg012", "cg013", "cg014", "cg015",
    "cg016", "cg017", "cg018", "cg019", "cg020",
  ]
  let sig_probes : Array[String] = [
    "cg001", "cg002", "cg003", "cg004", "cg005",
  ]
  let results = @src.mm_go_enrichment(
    sig_probes, all_probes, annotations, go_db,
  )
  // 5 unique GO terms in the database.
  assert_eq(results.length(), 5)
  // All p-values are valid.
  for r in results {
    assert_true(r.p_value() >= 0.0 && r.p_value() <= 1.0)
    assert_true(r.adj_p_value() >= 0.0 && r.adj_p_value() <= 1.0)
    assert_true(r.n_probes() >= 0)
    assert_true(r.n_sig_probes() >= 0)
  }
}

///|
test "mm_go_enrichment_finds_enriched" {
  let annotations = @src.mm_sample_annotations()
  let go_db = @src.mm_sample_go_database()
  let all_probes : Array[String] = [
    "cg001", "cg002", "cg003", "cg004", "cg005",
    "cg006", "cg007", "cg008", "cg009", "cg010",
    "cg011", "cg012", "cg013", "cg014", "cg015",
    "cg016", "cg017", "cg018", "cg019", "cg020",
  ]
  let sig_probes : Array[String] = [
    "cg001", "cg002", "cg003", "cg004", "cg005",
  ]
  let results = @src.mm_go_enrichment(
    sig_probes, all_probes, annotations, go_db,
  )
  // GO:0000001 (genes A, B) should have both A and B in the significant set.
  let mut found_enriched = false
  for r in results {
    if r.go_id() == "GO:0000001" {
      assert_eq(r.n_sig_probes(), 2)
      assert_true(r.p_value() >= 0.0 && r.p_value() <= 1.0)
      assert_true(r.expected() > 0.0)
      found_enriched = true
    }
  }
  assert_true(found_enriched)
}

// ===========================================================================
// mm_probe_bias_correction
// ===========================================================================

///|
test "mm_probe_bias_correction_sample" {
  let annotations = @src.mm_sample_annotations()
  let go_db = @src.mm_sample_go_database()
  let all_probes : Array[String] = [
    "cg001", "cg002", "cg003", "cg004", "cg005",
    "cg006", "cg007", "cg008", "cg009", "cg010",
    "cg011", "cg012", "cg013", "cg014", "cg015",
    "cg016", "cg017", "cg018", "cg019", "cg020",
  ]
  let sig_probes : Array[String] = [
    "cg001", "cg002", "cg003", "cg004", "cg005",
  ]
  let results = @src.mm_probe_bias_correction(
    sig_probes, all_probes, annotations, go_db,
  )
  assert_eq(results.length(), 5)
  for r in results {
    assert_true(r.p_value() >= 0.0 && r.p_value() <= 1.0)
    assert_true(r.adj_p_value() >= 0.0 && r.adj_p_value() <= 1.0)
    assert_true(r.expected() >= 0.0)
  }
}

// ===========================================================================
// Summary functions
// ===========================================================================

///|
test "mm_summary_nonempty" {
  let beta = @src.mm_sample_beta()
  let annotations = @src.mm_sample_annotations()
  let group = [false, false, false, false, true, true, true, true]
  let results = @src.mm_dmp_analysis(beta, group, annotations)
  let s = @src.mm_summary(results)
  assert_true(s.length() > 0)
  assert_true(s.contains("DMP Analysis Summary"))
  assert_true(s.contains("Total probes tested: 20"))
}

///|
test "mm_go_summary_nonempty" {
  let annotations = @src.mm_sample_annotations()
  let go_db = @src.mm_sample_go_database()
  let all_probes : Array[String] = [
    "cg001", "cg002", "cg003", "cg004", "cg005",
    "cg006", "cg007", "cg008", "cg009", "cg010",
    "cg011", "cg012", "cg013", "cg014", "cg015",
    "cg016", "cg017", "cg018", "cg019", "cg020",
  ]
  let sig_probes : Array[String] = [
    "cg001", "cg002", "cg003", "cg004", "cg005",
  ]
  let results = @src.mm_go_enrichment(
    sig_probes, all_probes, annotations, go_db,
  )
  let s = @src.mm_go_summary(results)
  assert_true(s.length() > 0)
  assert_true(s.contains("GO Enrichment Summary"))
  assert_true(s.contains("GO terms tested: 5"))
}

// ===========================================================================
// Sample data validation
// ===========================================================================

///|
test "mm_sample_annotations_validation" {
  let annotations = @src.mm_sample_annotations()
  assert_eq(annotations.length(), 20)
  // Verify first and last entries.
  assert_eq(annotations[0].probe_id(), "cg001")
  assert_eq(annotations[0].chr(), "chr1")
  assert_eq(annotations[0].gene(), "A")
  assert_eq(annotations[19].probe_id(), "cg020")
  assert_eq(annotations[19].chr(), "chr3")
  assert_eq(annotations[19].gene(), "D")
  // Verify gene diversity: genes A-H are present.
  let genes : Array[String] = Array::new()
  for ann in annotations {
    if !genes.contains(ann.gene()) {
      genes.push(ann.gene())
    }
  }
  assert_eq(genes.length(), 8)
}

///|
test "mm_sample_go_database_validation" {
  let db = @src.mm_sample_go_database()
  // 5 unique GO IDs, 10 total entries.
  let go1_genes = db.genes_for_go("GO:0000001")
  assert_eq(go1_genes.length(), 2)
  let go2_genes = db.genes_for_go("GO:0000002")
  assert_eq(go2_genes.length(), 2)
  let go3_genes = db.genes_for_go("GO:0000003")
  assert_eq(go3_genes.length(), 2)
  let go4_genes = db.genes_for_go("GO:0000004")
  assert_eq(go4_genes.length(), 2)
  let go5_genes = db.genes_for_go("GO:0000005")
  assert_eq(go5_genes.length(), 2)
  // Gene A is annotated to 3 GO terms.
  let a_gos = db.go_for_gene("A")
  assert_eq(a_gos.length(), 3)
}

///|
test "mm_sample_beta_validation" {
  let beta = @src.mm_sample_beta()
  // 20 probes x 8 samples.
  assert_eq(beta.length(), 20)
  for row in beta {
    assert_eq(row.length(), 8)
  }
  // First 5 probes: control (s<4) has low beta, treatment (s>=4) has high beta.
  for p in 0..<5 {
    assert_true(beta[p][0] < 0.3)
    assert_true(beta[p][4] > 0.6)
  }
  // Remaining probes: similar beta across groups.
  for p in 5..<20 {
    assert_true(beta[p][0] > 0.4 && beta[p][0] < 0.65)
    assert_true(beta[p][4] > 0.4 && beta[p][4] < 0.65)
  }
}

// ===========================================================================
// Edge cases
// ===========================================================================

///|
test "mm_edge_empty_arrays" {
  // Empty bh_correct.
  let adj = @src.mm_bh_correct([])
  assert_eq(adj.length(), 0)
  // t_test with empty arrays.
  let (t, p) = @src.mm_t_test([], [])
  assert_eq(t, 0.0)
  assert_eq(p, 1.0)
}

///|
test "mm_edge_single_probe" {
  let beta : Array[Array[Double]] = [[0.2, 0.3, 0.7, 0.8, 0.7, 0.8, 0.7, 0.8]]
  let annotations : Array[@src.MmProbeAnnotation] = [
    @src.MmProbeAnnotation::new("cg001", "chr1", 100, "A", "Island", "I"),
  ]
  let group = [false, false, false, false, true, true, true, true]
  let results = @src.mm_dmp_analysis(beta, group, annotations)
  assert_eq(results.length(), 1)
  assert_eq(results[0].probe_id(), "cg001")
}

///|
test "mm_edge_all_same_group" {
  let beta = @src.mm_sample_beta()
  let annotations = @src.mm_sample_annotations()
  // All samples in treatment group (no control).
  let group = [
    true, true, true, true, true, true, true, true,
  ]
  let results = @src.mm_dmp_analysis(beta, group, annotations)
  assert_eq(results.length(), 20)
  // With no control group, all p-values should be 1.0.
  for r in results {
    assert_true(r.p_value() > 0.99)
    assert_true(!r.is_significant(0.05))
  }
}

///|
test "mm_edge_t_test_small_samples" {
  // n < 2 per group -> returns (0.0, 1.0).
  let (t, p) = @src.mm_t_test([1.0], [2.0])
  assert_eq(t, 0.0)
  assert_eq(p, 1.0)
}

///|
test "mm_edge_filter_no_match" {
  let annotations = @src.mm_sample_annotations()
  let filtered = @src.mm_filter_probes(annotations, "chrY", "all")
  assert_eq(filtered.length(), 0)
}

///|
test "mm_edge_fisher_exact_empty" {
  // All zeros -> n=0, returns 1.0.
  let p = @src.mm_fisher_exact(0, 0, 0, 0)
  assert_true(p > 0.99)
}
