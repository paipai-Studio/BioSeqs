///|
/// Tests for the NanoString module.
/// Covers sample/data containers, normalization, QC, differential expression,
/// utility functions, summaries, sample data validation, and edge cases.

// ===========================================================================
// NsSample creation and accessors
// ===========================================================================

///|
test "ns_sample_creation" {
  let raw = [100, 200, 300, 50, 60, 10, 20]
  let pos_indices = [3, 4]
  let neg_indices = [5, 6]
  let hk_indices = [0, 1]
  let endo_indices = [0, 1, 2]
  let sample = @src.NsSample::new(
    "S1", raw, pos_indices, neg_indices, hk_indices, endo_indices,
  )
  assert_eq(sample.sample_id(), "S1")
  assert_eq(sample.raw_counts().length(), 7)
  assert_eq(sample.raw_counts()[0], 100)
}

///|
test "ns_sample_accessors" {
  let raw = [100, 200, 300, 50, 60, 10, 20]
  let sample = @src.NsSample::new(
    "S1", raw, [3, 4], [5, 6], [0, 1], [0, 1, 2],
  )
  let pos = sample.positive_controls()
  assert_eq(pos.length(), 2)
  assert_eq(pos[0], 50)
  assert_eq(pos[1], 60)
  let neg = sample.negative_controls()
  assert_eq(neg.length(), 2)
  assert_eq(neg[0], 10)
  assert_eq(neg[1], 20)
  let hk = sample.housekeeping()
  assert_eq(hk.length(), 2)
  assert_eq(hk[0], 100)
  assert_eq(hk[1], 200)
  let endo = sample.endogenous()
  assert_eq(endo.length(), 3)
  assert_eq(endo[0], 100)
  assert_eq(endo[1], 200)
  assert_eq(endo[2], 300)
}

// ===========================================================================
// NsNanostringData creation, add_sample, n_samples, n_genes
// ===========================================================================

///|
test "ns_nanostring_data_creation" {
  let gene_names = ["G1", "G2", "G3"]
  let data = @src.NsNanostringData::new(gene_names)
  assert_eq(data.n_genes(), 3)
  assert_eq(data.n_samples(), 0)
  assert_eq(data.gene_names().length(), 3)
}

///|
test "ns_nanostring_data_add_sample" {
  let gene_names = ["G1", "G2", "G3"]
  let data = @src.NsNanostringData::new(gene_names)
  let raw = [100, 200, 300, 50, 60, 10, 20]
  let sample = @src.NsSample::new(
    "S1", raw, [3, 4], [5, 6], [0, 1], [0, 1, 2],
  )
  data.add_sample(sample)
  assert_eq(data.n_samples(), 1)
  assert_eq(data.n_genes(), 3)
  let samples = data.samples()
  assert_eq(samples.length(), 1)
  assert_eq(samples[0].sample_id(), "S1")
}

// ===========================================================================
// ns_geometric_mean
// ===========================================================================

///|
test "ns_geometric_mean_known" {
  // geometric mean of [2, 8] = sqrt(16) = 4.
  let gm = @src.ns_geometric_mean([2.0, 8.0])
  assert_true((gm - 4.0).abs() < 0.01)
  // geometric mean of [1, 1, 1] = 1.
  let gm2 = @src.ns_geometric_mean([1.0, 1.0, 1.0])
  assert_true((gm2 - 1.0).abs() < 0.01)
  // geometric mean of [4, 9] = sqrt(36) = 6.
  let gm3 = @src.ns_geometric_mean([4.0, 9.0])
  assert_true((gm3 - 6.0).abs() < 0.01)
}

///|
test "ns_geometric_mean_single" {
  let gm = @src.ns_geometric_mean([5.0])
  assert_true((gm - 5.0).abs() < 0.01)
}

///|
test "ns_geometric_mean_all_same" {
  let gm = @src.ns_geometric_mean([3.0, 3.0, 3.0, 3.0])
  assert_true((gm - 3.0).abs() < 0.01)
}

///|
test "ns_geometric_mean_no_positives" {
  // No positive values -> returns 0.0.
  let gm = @src.ns_geometric_mean([0.0, 0.0, -1.0])
  assert_eq(gm, 0.0)
}

// ===========================================================================
// ns_log2_counts
// ===========================================================================

///|
test "ns_log2_counts_known" {
  let result = @src.ns_log2_counts([3.5])
  // log2(3.5 + 0.5) = log2(4) = 2.
  assert_true((result[0] - 2.0).abs() < 0.01)
}

///|
test "ns_log2_counts_zero_pseudocount" {
  // log2(0 + 0.5) = log2(0.5) = -1.
  let result = @src.ns_log2_counts([0.0])
  assert_true((result[0] - (-1.0)).abs() < 0.01)
  // Multiple values.
  let result2 = @src.ns_log2_counts([0.0, 1.0, 7.0])
  // log2(0.5) = -1, log2(1.5) ≈ 0.585, log2(7.5) ≈ 2.907
  assert_true((result2[0] - (-1.0)).abs() < 0.01)
  assert_true((result2[1] - 0.585).abs() < 0.01)
  assert_true((result2[2] - 2.907).abs() < 0.01)
}

// ===========================================================================
// Positive control normalization
// ===========================================================================

///|
test "ns_positive_control_norm_basic" {
  let data = @src.ns_sample_data()
  let results = @src.ns_positive_control_norm(data)
  assert_eq(results.length(), 6)
  for r in results {
    assert_eq(r.method(), "positive_control")
    assert_eq(r.norm_factors().length(), 1)
    assert_eq(r.normalized_counts().length(), 20)
  }
}

///|
test "ns_positive_control_norm_factors" {
  let data = @src.ns_sample_data()
  let results = @src.ns_positive_control_norm(data)
  // All samples have similar positive control counts, so factors should be
  // close to 1.0.
  for r in results {
    assert_true(r.norm_factors()[0] > 0.9)
    assert_true(r.norm_factors()[0] < 1.1)
  }
  // Normalized counts should be non-negative.
  for r in results {
    for c in r.normalized_counts() {
      assert_true(c >= 0.0)
    }
  }
}

// ===========================================================================
// Negative control normalization
// ===========================================================================

///|
test "ns_negative_control_norm_basic" {
  let data = @src.ns_sample_data()
  let results = @src.ns_negative_control_norm(data)
  assert_eq(results.length(), 6)
  for r in results {
    assert_eq(r.method(), "negative_control")
    assert_eq(r.normalized_counts().length(), 20)
    // Background-subtracted values floored at 0.5.
    for c in r.normalized_counts() {
      assert_true(c >= 0.5)
    }
  }
}

// ===========================================================================
// Housekeeping normalization
// ===========================================================================

///|
test "ns_housekeeping_norm_basic" {
  let data = @src.ns_sample_data()
  let results = @src.ns_housekeeping_norm(data)
  assert_eq(results.length(), 6)
  for r in results {
    assert_eq(r.method(), "housekeeping")
    assert_eq(r.norm_factors().length(), 1)
    assert_eq(r.normalized_counts().length(), 20)
  }
}

// ===========================================================================
// Content normalization
// ===========================================================================

///|
test "ns_content_norm_basic" {
  let data = @src.ns_sample_data()
  let results = @src.ns_content_norm(data)
  assert_eq(results.length(), 6)
  for r in results {
    assert_eq(r.method(), "content")
    // Content norm has two factors: pos and hk.
    assert_eq(r.norm_factors().length(), 2)
    assert_eq(r.normalized_counts().length(), 20)
  }
}

// ===========================================================================
// ns_t_test
// ===========================================================================

///|
test "ns_t_test_same_group" {
  let g1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let g2 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let (t, p) = @src.ns_t_test(g1, g2)
  assert_true(t.abs() < 0.01)
  assert_true(p > 0.99)
}

///|
test "ns_t_test_different_groups" {
  let g1 = [10.0, 20.0, 30.0, 40.0, 50.0]
  let g2 = [100.0, 110.0, 120.0, 130.0, 140.0]
  let (t, p) = @src.ns_t_test(g1, g2)
  assert_true(t < -5.0)
  assert_true(p < 0.01)
}

// ===========================================================================
// ns_bh_correct
// ===========================================================================

///|
test "ns_bh_correct_basic" {
  let pvalues = [0.01, 0.5, 0.001, 0.2]
  let adj = @src.ns_bh_correct(pvalues)
  assert_eq(adj.length(), 4)
  for v in adj {
    assert_true(v >= 0.0 && v <= 1.0)
  }
  // Smallest p (index 2) gets smallest adj.
  assert_true(adj[2] <= adj[0])
  assert_true(adj[0] <= adj[3])
  assert_true(adj[3] <= adj[1])
}

///|
test "ns_bh_correct_empty" {
  let adj = @src.ns_bh_correct([])
  assert_eq(adj.length(), 0)
}

// ===========================================================================
// ns_differential_expression
// ===========================================================================

///|
test "ns_differential_expression_sample" {
  let data = @src.ns_sample_data()
  let norm = @src.ns_positive_control_norm(data)
  let group = [false, false, false, true, true, true]
  let gene_names = data.gene_names()
  let results = @src.ns_differential_expression(norm, group, gene_names)
  assert_eq(results.length(), 20)
  for r in results {
    assert_true(r.p_value() >= 0.0 && r.p_value() <= 1.0)
    assert_true(r.adj_p_value() >= 0.0 && r.adj_p_value() <= 1.0)
  }
}

///|
test "ns_differential_expression_significant" {
  let data = @src.ns_sample_data()
  let norm = @src.ns_positive_control_norm(data)
  let group = [false, false, false, true, true, true]
  let gene_names = data.gene_names()
  let results = @src.ns_differential_expression(norm, group, gene_names)
  let mut n_sig = 0
  let mut n_up = 0
  let mut n_down = 0
  for r in results {
    if r.significant() {
      n_sig = n_sig + 1
      if r.log_fc() > 0.0 {
        n_up = n_up + 1
      } else {
        n_down = n_down + 1
      }
    }
  }
  // Genes 0-4 up-regulated, genes 5-9 down-regulated in treatment.
  assert_true(n_sig >= 5)
  assert_true(n_up > 0)
  assert_true(n_down > 0)
}

// ===========================================================================
// ns_imaging_qc
// ===========================================================================

///|
test "ns_imaging_qc" {
  let data = @src.ns_sample_data()
  // Positive control counts are ~6400 avg, so threshold 100 should pass all.
  let pass = @src.ns_imaging_qc(data, 100)
  assert_eq(pass.length(), 6)
  for p in pass {
    assert_true(p)
  }
  // With a very high threshold, all should fail.
  let fail = @src.ns_imaging_qc(data, 10000)
  for p in fail {
    assert_true(!p)
  }
}

// ===========================================================================
// ns_pos_linearity_qc
// ===========================================================================

///|
test "ns_pos_linearity_qc" {
  let data = @src.ns_sample_data()
  let (r_sq, slope) = @src.ns_pos_linearity_qc(data)
  // Positive controls follow a dilution series with high linearity.
  assert_true(r_sq > 0.95)
  // Slope should be close to 1.0 (counts scale linearly with concentration).
  assert_true(slope > 0.8 && slope < 1.2)
}

// ===========================================================================
// ns_loq_qc
// ===========================================================================

///|
test "ns_loq_qc" {
  let data = @src.ns_sample_data()
  let counts = @src.ns_loq_qc(data, 200)
  assert_eq(counts.length(), 6)
  // With LOQ=200, each sample should have some genes below threshold.
  for c in counts {
    assert_true(c > 0)
    assert_true(c <= 20)
  }
}

// ===========================================================================
// Summary functions
// ===========================================================================

///|
test "ns_qc_summary_nonempty" {
  let data = @src.ns_sample_data()
  let s = @src.ns_qc_summary(data)
  assert_true(s.length() > 0)
  assert_true(s.contains("NanoString QC Summary"))
  assert_true(s.contains("Samples: 6"))
}

///|
test "ns_summary_nonempty" {
  let data = @src.ns_sample_data()
  let s = @src.ns_summary(data)
  assert_true(s.length() > 0)
  assert_true(s.contains("NanoString Data Summary"))
  assert_true(s.contains("Samples: 6"))
  assert_true(s.contains("Genes: 20"))
}

///|
test "ns_de_summary_nonempty" {
  let data = @src.ns_sample_data()
  let norm = @src.ns_positive_control_norm(data)
  let group = [false, false, false, true, true, true]
  let gene_names = data.gene_names()
  let results = @src.ns_differential_expression(norm, group, gene_names)
  let s = @src.ns_de_summary(results)
  assert_true(s.length() > 0)
  assert_true(s.contains("Differential Expression Summary"))
  assert_true(s.contains("Genes tested: 20"))
}

// ===========================================================================
// Sample data validation
// ===========================================================================

///|
test "ns_sample_data_validation" {
  let data = @src.ns_sample_data()
  assert_eq(data.n_samples(), 6)
  assert_eq(data.n_genes(), 20)
  assert_eq(data.gene_names().length(), 20)
  assert_eq(data.gene_names()[0], "GENE01")
  assert_eq(data.gene_names()[19], "GENE20")
}

///|
test "ns_sample_data_sample_ids" {
  let data = @src.ns_sample_data()
  let samples = data.samples()
  assert_eq(samples[0].sample_id(), "C1")
  assert_eq(samples[1].sample_id(), "C2")
  assert_eq(samples[2].sample_id(), "C3")
  assert_eq(samples[3].sample_id(), "T1")
  assert_eq(samples[4].sample_id(), "T2")
  assert_eq(samples[5].sample_id(), "T3")
}

///|
test "ns_sample_data_controls" {
  let data = @src.ns_sample_data()
  for sample in data.samples() {
    // Each sample has 6 positive controls, 6 negative controls,
    // 4 housekeeping, and 20 endogenous.
    assert_eq(sample.positive_controls().length(), 6)
    assert_eq(sample.negative_controls().length(), 6)
    assert_eq(sample.housekeeping().length(), 4)
    assert_eq(sample.endogenous().length(), 20)
    // Raw counts = 20 endo + 6 pos + 6 neg = 32.
    assert_eq(sample.raw_counts().length(), 32)
  }
}

// ===========================================================================
// Edge cases
// ===========================================================================

///|
test "ns_edge_single_sample" {
  let gene_names = ["G1", "G2", "G3"]
  let data = @src.NsNanostringData::new(gene_names)
  let raw = [100, 200, 300, 50, 60, 10, 20]
  let sample = @src.NsSample::new(
    "S1", raw, [3, 4], [5, 6], [0, 1], [0, 1, 2],
  )
  data.add_sample(sample)
  assert_eq(data.n_samples(), 1)
  let results = @src.ns_positive_control_norm(data)
  assert_eq(results.length(), 1)
  assert_eq(results[0].normalized_counts().length(), 3)
}

///|
test "ns_edge_single_gene" {
  let gene_names = ["G1"]
  let data = @src.NsNanostringData::new(gene_names)
  let raw = [100, 50, 10]
  let sample = @src.NsSample::new(
    "S1", raw, [1], [2], [0], [0],
  )
  data.add_sample(sample)
  assert_eq(data.n_genes(), 1)
  let results = @src.ns_positive_control_norm(data)
  assert_eq(results[0].normalized_counts().length(), 1)
}

///|
test "ns_edge_all_zeros" {
  let gene_names = ["G1", "G2", "G3"]
  let data = @src.NsNanostringData::new(gene_names)
  let raw = [0, 0, 0, 0, 0, 0, 0, 0, 0]
  let sample = @src.NsSample::new(
    "S1", raw, [3, 4, 5], [6, 7, 8], [0, 1, 2], [0, 1, 2],
  )
  data.add_sample(sample)
  let results = @src.ns_positive_control_norm(data)
  assert_eq(results.length(), 1)
  // Geometric mean of all-zero positive controls = 0.0, so factor = 1.0.
  assert_true((results[0].norm_factors()[0] - 1.0).abs() < 0.01)
  // Normalized counts should all be 0.0.
  for c in results[0].normalized_counts() {
    assert_eq(c, 0.0)
  }
}

///|
test "ns_edge_t_test_small_samples" {
  let (t, p) = @src.ns_t_test([1.0], [2.0])
  assert_eq(t, 0.0)
  assert_eq(p, 1.0)
}

///|
test "ns_edge_log2_counts_empty" {
  let result = @src.ns_log2_counts([])
  assert_eq(result.length(), 0)
}
