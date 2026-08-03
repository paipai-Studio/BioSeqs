///|
/// Tests for the methylKit module.
/// Mirrors the Bioconductor methylKit package API.

///|
test "methylkit_cytosine_creation" {
  let mc = @src.MethylCytosine::new("chr1", 1000, "+", "CG", 8, 2)
  assert_eq(@src.MethylCytosine::chr(mc), "chr1")
  assert_eq(@src.MethylCytosine::position(mc), 1000)
  assert_eq(@src.MethylCytosine::strand(mc), "+")
  assert_eq(@src.MethylCytosine::context(mc), "CG")
  assert_eq(@src.MethylCytosine::methylated(mc), 8)
  assert_eq(@src.MethylCytosine::unmethylated(mc), 2)
}

///|
test "methylkit_cytosine_coverage_pct" {
  // methylated=8, unmethylated=2 -> coverage=10, pct=80.0
  let mc = @src.MethylCytosine::new("chr1", 100, "+", "CG", 8, 2)
  assert_eq(@src.MethylCytosine::coverage(mc), 10)
  assert_true((@src.MethylCytosine::methylation_pct(mc) - 80.0).abs() < 0.001)
  // 50% methylation
  let mc_half = @src.MethylCytosine::new("chr1", 200, "-", "CG", 5, 5)
  assert_eq(@src.MethylCytosine::coverage(mc_half), 10)
  assert_true((@src.MethylCytosine::methylation_pct(mc_half) - 50.0).abs() < 0.001)
  // Zero coverage -> coverage 0 and pct 0.0
  let mc0 = @src.MethylCytosine::new("chr1", 300, "+", "CG", 0, 0)
  assert_eq(@src.MethylCytosine::coverage(mc0), 0)
  assert_eq(@src.MethylCytosine::methylation_pct(mc0), 0.0)
}

///|
test "methylkit_cytosine_to_string" {
  let mc = @src.MethylCytosine::new("chr1", 100, "+", "CG", 8, 2)
  let s = @src.MethylCytosine::to_string(mc)
  assert_true(s.contains("MethylCytosine"))
  assert_true(s.contains("chr1"))
  assert_true(s.contains("100"))
  assert_true(s.contains("coverage=10"))
}

///|
test "methylkit_sample_new" {
  let s = @src.MethylSample::new("empty")
  assert_eq(@src.MethylSample::sample_id(s), "empty")
  assert_eq(@src.MethylSample::n_cpgs(s), 0)
  assert_eq(@src.MethylSample::coverage(s).length(), 0)
}

///|
test "methylkit_sample_from_coverage" {
  let cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 8, 2),
  ]
  let s = @src.MethylSample::from_coverage("s1", cov)
  assert_eq(@src.MethylSample::sample_id(s), "s1")
  assert_eq(@src.MethylSample::n_cpgs(s), 2)
  assert_eq(@src.MethylSample::coverage(s).length(), 2)
}

///|
test "methylkit_sample_add_cytosine" {
  let s = @src.MethylSample::new("s1")
  assert_eq(@src.MethylSample::n_cpgs(s), 0)
  @src.MethylSample::add_cytosine(
    s,
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5),
  )
  assert_eq(@src.MethylSample::n_cpgs(s), 1)
  @src.MethylSample::add_cytosine(
    s,
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 8, 2),
  )
  assert_eq(@src.MethylSample::n_cpgs(s), 2)
  assert_eq(@src.MethylSample::coverage(s).length(), 2)
  // The added cytosine is retrievable.
  let found = @src.MethylSample::get_cytosine(s, "chr1", 200)
  assert_true(found is Some(_))
}

///|
test "methylkit_sample_get_cytosine_found" {
  let cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 8, 2),
  ]
  let s = @src.MethylSample::from_coverage("s1", cov)
  let found = @src.MethylSample::get_cytosine(s, "chr1", 200)
  assert_true(found is Some(_))
  let mc = found.unwrap()
  assert_eq(@src.MethylCytosine::methylated(mc), 8)
  assert_eq(@src.MethylCytosine::unmethylated(mc), 2)
}

///|
test "methylkit_sample_get_cytosine_not_found" {
  let cov = [@src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5)]
  let s = @src.MethylSample::from_coverage("s1", cov)
  // Wrong chromosome.
  assert_true(@src.MethylSample::get_cytosine(s, "chr2", 100) is None)
  // Wrong position.
  assert_true(@src.MethylSample::get_cytosine(s, "chr1", 999) is None)
}

///|
test "methylkit_sample_methylation_pcts" {
  let cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 8, 2),
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 2, 8),
  ]
  let s = @src.MethylSample::from_coverage("s1", cov)
  let pcts = @src.MethylSample::methylation_pcts(s)
  assert_eq(pcts.length(), 3)
  assert_true((pcts[0] - 80.0).abs() < 0.001)
  assert_true((pcts[1] - 50.0).abs() < 0.001)
  assert_true((pcts[2] - 20.0).abs() < 0.001)
}

///|
test "methylkit_sample_mean_methylation" {
  let cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 8, 2),
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 2, 8),
  ]
  let s = @src.MethylSample::from_coverage("s1", cov)
  // Mean of 80% and 20% = 50%.
  assert_true((@src.MethylSample::mean_methylation(s) - 50.0).abs() < 0.001)
}

///|
test "methylkit_sample_filter_by_coverage" {
  let cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5), // cov=10, kept
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 1, 1), // cov=2, too low
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 500, 500), // cov=1000, too high
  ]
  let s = @src.MethylSample::from_coverage("s1", cov)
  let filtered = @src.MethylSample::filter_by_coverage(s, 5, 100)
  assert_eq(@src.MethylSample::n_cpgs(filtered), 1)
  assert_eq(
    @src.MethylCytosine::position(@src.MethylSample::coverage(filtered)[0]),
    100,
  )
}

///|
test "methylkit_config_defaults" {
  let config = @src.MethylKitConfig::new()
  assert_eq(@src.MethylKitConfig::min_coverage(config), 10)
  assert_eq(@src.MethylKitConfig::max_coverage(config), 999)
  assert_eq(@src.MethylKitConfig::context_filter(config), "CG")
  assert_true((@src.MethylKitConfig::min_perc_samples(config) - 0.6).abs() < 0.001)
  assert_true((@src.MethylKitConfig::min_diff(config) - 25.0).abs() < 0.001)
  assert_true((@src.MethylKitConfig::q_threshold(config) - 0.05).abs() < 0.001)
}

///|
test "methylkit_config_with_params" {
  let config = @src.MethylKitConfig::with_params(5, 500, "CHG", 0.8, 30.0, 0.01)
  assert_eq(@src.MethylKitConfig::min_coverage(config), 5)
  assert_eq(@src.MethylKitConfig::max_coverage(config), 500)
  assert_eq(@src.MethylKitConfig::context_filter(config), "CHG")
  assert_true((@src.MethylKitConfig::min_perc_samples(config) - 0.8).abs() < 0.001)
  assert_true((@src.MethylKitConfig::min_diff(config) - 30.0).abs() < 0.001)
  assert_true((@src.MethylKitConfig::q_threshold(config) - 0.01).abs() < 0.001)
}

///|
test "methylkit_filter_sample_coverage_and_context" {
  let cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5), // cov=10, CG -> kept
    @src.MethylCytosine::new("chr1", 200, "+", "CHG", 5, 5), // wrong context
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 2, 2), // cov=4 < 10
    @src.MethylCytosine::new("chr1", 400, "+", "CG", 500, 500), // cov=1000 > 999
  ]
  let s = @src.MethylSample::from_coverage("s1", cov)
  let config = @src.MethylKitConfig::new()
  let filtered = @src.methylkit_filter_sample(s, config)
  assert_eq(@src.MethylSample::n_cpgs(filtered), 1)
  assert_eq(
    @src.MethylCytosine::position(@src.MethylSample::coverage(filtered)[0]),
    100,
  )
}

///|
test "methylkit_filter_samples_multi" {
  // Two samples that share only position 100.
  let s1_cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 5, 5),
  ]
  let s2_cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 5, 5),
  ]
  let s1 = @src.MethylSample::from_coverage("s1", s1_cov)
  let s2 = @src.MethylSample::from_coverage("s2", s2_cov)
  let config = @src.MethylKitConfig::new()
  let filtered = @src.methylkit_filter_samples([s1, s2], config)
  assert_eq(filtered.length(), 2)
  // min_perc_samples=0.6 * 2 = 1.2 -> a CpG must appear in >=2 samples.
  // Only position 100 is shared, so each sample retains 1 CpG.
  assert_eq(@src.MethylSample::n_cpgs(filtered[0]), 1)
  assert_eq(@src.MethylSample::n_cpgs(filtered[1]), 1)
  assert_eq(
    @src.MethylCytosine::position(@src.MethylSample::coverage(filtered[0])[0]),
    100,
  )
  assert_eq(
    @src.MethylCytosine::position(@src.MethylSample::coverage(filtered[1])[0]),
    100,
  )
}

///|
test "methylkit_normalize_sample" {
  // Coverages 10, 20, 30 -> median = 20.
  let cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5), // cov=10
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 10, 10), // cov=20
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 15, 15), // cov=30
  ]
  let s = @src.MethylSample::from_coverage("s", cov)
  let norm = @src.methylkit_normalize(s)
  assert_eq(@src.MethylSample::n_cpgs(norm), 3)
  // After normalization every cytosine is scaled to the median coverage 20.
  let mut i = 0
  while i < @src.MethylSample::n_cpgs(norm) {
    assert_eq(
      @src.MethylCytosine::coverage(@src.MethylSample::coverage(norm)[i]),
      20,
    )
    i = i + 1
  }
}

///|
test "methylkit_normalize_samples" {
  let cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 10, 10),
  ]
  let s = @src.MethylSample::from_coverage("s", cov)
  let normalized = @src.methylkit_normalize_samples([s])
  assert_eq(normalized.length(), 1)
  assert_eq(@src.MethylSample::n_cpgs(normalized[0]), 2)
}

///|
test "methylkit_fisher_test_same" {
  // Identical methylation (10,10) vs (10,10) -> p-value ~1.0.
  let mc1 = @src.MethylCytosine::new("chr1", 100, "+", "CG", 10, 10)
  let mc2 = @src.MethylCytosine::new("chr1", 200, "+", "CG", 10, 10)
  let p = @src.methylkit_fisher_test(mc1, mc2)
  assert_true(p > 0.99)
}

///|
test "methylkit_fisher_test_different" {
  // Extreme difference (10,0) vs (0,10) -> very small p-value.
  let mc1 = @src.MethylCytosine::new("chr1", 100, "+", "CG", 10, 0)
  let mc2 = @src.MethylCytosine::new("chr1", 200, "+", "CG", 0, 10)
  let p = @src.methylkit_fisher_test(mc1, mc2)
  assert_true(p < 0.01)
}

///|
test "methylkit_bh_correct_basic" {
  let pvalues = [0.01, 0.04, 0.03, 0.20]
  let adj = @src.methylkit_bh_correct(pvalues)
  assert_eq(adj.length(), 4)
  // Adjusted values are valid probabilities in [0, 1].
  let mut i = 0
  while i < adj.length() {
    assert_true(adj[i] >= 0.0 && adj[i] <= 1.0)
    i = i + 1
  }
  // The smallest p-value (index 0, value 0.01) gets the smallest q-value.
  assert_true(adj[0] <= adj[3])
  // 0.01 * 4 / 1 = 0.04.
  assert_true(adj[0] < 0.05)
  // The largest p-value (index 3, value 0.20) is unchanged after correction.
  assert_true(adj[3] >= 0.19 && adj[3] <= 0.21)
}

///|
test "methylkit_bh_correct_empty" {
  let adj = @src.methylkit_bh_correct([])
  assert_eq(adj.length(), 0)
}

///|
test "methylkit_diff_meth_significant" {
  let (control, treatment) = @src.methylkit_sample_data()
  let config = @src.MethylKitConfig::new()
  let results = @src.methylkit_diff_meth(control, treatment, config)
  assert_eq(results.length(), 20)
  let mut n_sig = 0
  let mut n_hyper = 0
  let mut n_hypo = 0
  for r in results {
    if @src.MethylDiffResult::significant(r) {
      n_sig = n_sig + 1
      if @src.MethylDiffResult::meth_diff(r) > 0.0 {
        n_hyper = n_hyper + 1
      } else {
        n_hypo = n_hypo + 1
      }
    }
  }
  // 5 hyper + 5 hypo differential CpGs.
  assert_eq(n_sig, 10)
  assert_eq(n_hyper, 5)
  assert_eq(n_hypo, 5)
}

///|
test "methylkit_tiling" {
  let (control, treatment) = @src.methylkit_sample_data()
  let results = @src.methylkit_tiling([control, treatment], 5000, 5000)
  assert_true(results.length() > 0)
  // Each tiled window is on chr1.
  assert_eq(@src.MethylDiffResult::chr(results[0]), "chr1")
  // At least one window should be significant given the differential CpGs.
  let mut n_sig = 0
  for r in results {
    if @src.MethylDiffResult::significant(r) {
      n_sig = n_sig + 1
    }
  }
  assert_true(n_sig > 0)
}

///|
test "methylkit_get_dmr" {
  let (control, treatment) = @src.methylkit_sample_data()
  let config = @src.MethylKitConfig::new()
  let results = @src.methylkit_diff_meth(control, treatment, config)
  let dmrs = @src.methylkit_get_dmr(results)
  assert_true(dmrs.length() >= 1)
  let (chr, start, end) = dmrs[0]
  assert_eq(chr, "chr1")
  assert_true(end >= start)
}

///|
test "methylkit_sample_correlation_same" {
  // Sample with varying methylation so variance is non-zero.
  let cov = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 8, 2), // 80%
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 5, 5), // 50%
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 2, 8), // 20%
  ]
  let s = @src.MethylSample::from_coverage("s", cov)
  let corr = @src.methylkit_sample_correlation(s, s)
  assert_true((corr - 1.0).abs() < 0.001)
}

///|
test "methylkit_sample_correlation_different" {
  // Two samples that are anti-correlated.
  let cov1 = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 8, 2), // 80%
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 5, 5), // 50%
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 2, 8), // 20%
  ]
  let cov2 = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 2, 8), // 20%
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 5, 5), // 50%
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 8, 2), // 80%
  ]
  let s1 = @src.MethylSample::from_coverage("s1", cov1)
  let s2 = @src.MethylSample::from_coverage("s2", cov2)
  let corr = @src.methylkit_sample_correlation(s1, s2)
  assert_true(corr < 1.0)
}

///|
test "methylkit_correlation_matrix" {
  let cov1 = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 8, 2), // 80%
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 2, 8), // 20%
  ]
  let cov2 = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 6, 4), // 60%
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 4, 6), // 40%
  ]
  let s1 = @src.MethylSample::from_coverage("s1", cov1)
  let s2 = @src.MethylSample::from_coverage("s2", cov2)
  let matrix = @src.methylkit_correlation_matrix([s1, s2])
  assert_eq(matrix.length(), 2)
  assert_eq(matrix[0].length(), 2)
  assert_eq(matrix[1].length(), 2)
  // Diagonal entries are self-correlations -> ~1.0.
  assert_true((matrix[0][0] - 1.0).abs() < 0.001)
  assert_true((matrix[1][1] - 1.0).abs() < 0.001)
}

///|
test "methylkit_cluster_samples" {
  let (control, treatment) = @src.methylkit_sample_data()
  let labels = @src.methylkit_cluster_samples([control, treatment], 2)
  assert_eq(labels.length(), 2)
  // With k >= n, each sample gets its own cluster label.
  assert_true(labels[0] != labels[1])
}

///|
test "methylkit_cluster_samples_k1" {
  let (control, treatment) = @src.methylkit_sample_data()
  let labels = @src.methylkit_cluster_samples([control, treatment], 1)
  assert_eq(labels.length(), 2)
  // k=1 -> all samples in cluster 0.
  assert_eq(labels[0], 0)
  assert_eq(labels[1], 0)
}

///|
test "methylkit_sample_summary" {
  let (control, _) = @src.methylkit_sample_data()
  let s = @src.methylkit_sample_summary(control)
  assert_true(s.contains("MethylSample Summary"))
  assert_true(s.contains("control"))
  assert_true(s.contains("Number of CpGs: 20"))
}

///|
test "methylkit_diff_summary" {
  let (control, treatment) = @src.methylkit_sample_data()
  let config = @src.MethylKitConfig::new()
  let results = @src.methylkit_diff_meth(control, treatment, config)
  let s = @src.methylkit_diff_summary(results)
  assert_true(s.contains("Differential Methylation Summary"))
  assert_true(s.contains("Total CpGs tested: 20"))
  assert_true(s.contains("Significant DMPs: 10"))
  assert_true(s.contains("Hyper-methylated: 5"))
  assert_true(s.contains("Hypo-methylated: 5"))
}

///|
test "methylkit_cpg_to_bed" {
  let mc = @src.MethylCytosine::new("chr1", 100, "+", "CG", 8, 2)
  let bed = @src.methylkit_cpg_to_bed(mc)
  // chr1\t100\t101\t80.0
  assert_true(bed.contains("chr1"))
  assert_true(bed.contains("100"))
  assert_true(bed.contains("101"))
}

///|
test "methylkit_dmr_to_bed" {
  let bed = @src.methylkit_dmr_to_bed(("chr1", 1000, 5000))
  assert_true(bed.contains("chr1"))
  assert_true(bed.contains("1000"))
  assert_true(bed.contains("5000"))
}

///|
test "methylkit_sample_data" {
  let (control, treatment) = @src.methylkit_sample_data()
  assert_eq(@src.MethylSample::sample_id(control), "control")
  assert_eq(@src.MethylSample::sample_id(treatment), "treatment")
  assert_eq(@src.MethylSample::n_cpgs(control), 20)
  assert_eq(@src.MethylSample::n_cpgs(treatment), 20)
  // Control: every CpG at 50% methylation.
  assert_true((@src.MethylSample::mean_methylation(control) - 50.0).abs() < 0.001)
  // Treatment: 5 hyper (90%) + 5 hypo (10%) + 10 same (50%) -> mean 50%.
  assert_true((@src.MethylSample::mean_methylation(treatment) - 50.0).abs() < 0.001)
  // Treatment CpG 0: hyper-methylated (18, 2).
  let t0 = @src.MethylSample::coverage(treatment)[0]
  assert_eq(@src.MethylCytosine::methylated(t0), 18)
  assert_eq(@src.MethylCytosine::unmethylated(t0), 2)
  // Treatment CpG 5: hypo-methylated (2, 18).
  let t5 = @src.MethylSample::coverage(treatment)[5]
  assert_eq(@src.MethylCytosine::methylated(t5), 2)
  assert_eq(@src.MethylCytosine::unmethylated(t5), 18)
  // Treatment CpG 10: same as control (10, 10).
  let t10 = @src.MethylSample::coverage(treatment)[10]
  assert_eq(@src.MethylCytosine::methylated(t10), 10)
  assert_eq(@src.MethylCytosine::unmethylated(t10), 10)
}

///|
test "methylkit_empty_sample" {
  let s = @src.MethylSample::new("empty")
  assert_eq(@src.MethylSample::n_cpgs(s), 0)
  assert_eq(@src.MethylSample::methylation_pcts(s).length(), 0)
  assert_eq(@src.MethylSample::mean_methylation(s), 0.0)
  assert_true(@src.MethylSample::get_cytosine(s, "chr1", 100) is None)
  // Normalizing an empty sample returns an empty sample.
  let norm = @src.methylkit_normalize(s)
  assert_eq(@src.MethylSample::n_cpgs(norm), 0)
  // Filtering an empty sample returns an empty sample.
  let config = @src.MethylKitConfig::new()
  let filtered = @src.methylkit_filter_sample(s, config)
  assert_eq(@src.MethylSample::n_cpgs(filtered), 0)
}

///|
test "methylkit_single_cpg" {
  let cov = [@src.MethylCytosine::new("chr1", 100, "+", "CG", 8, 2)]
  let s = @src.MethylSample::from_coverage("s", cov)
  assert_eq(@src.MethylSample::n_cpgs(s), 1)
  assert_true((@src.MethylSample::mean_methylation(s) - 80.0).abs() < 0.001)
  // Correlation requires >= 2 shared CpGs -> returns 0.0.
  let corr = @src.methylkit_sample_correlation(s, s)
  assert_eq(corr, 0.0)
}

///|
test "methylkit_all_same_methylation" {
  // Two samples with identical methylation -> no significant DMPs.
  let cov1 = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 5, 5),
  ]
  let cov2 = [
    @src.MethylCytosine::new("chr1", 100, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 200, "+", "CG", 5, 5),
    @src.MethylCytosine::new("chr1", 300, "+", "CG", 5, 5),
  ]
  let s1 = @src.MethylSample::from_coverage("s1", cov1)
  let s2 = @src.MethylSample::from_coverage("s2", cov2)
  let config = @src.MethylKitConfig::new()
  let results = @src.methylkit_diff_meth(s1, s2, config)
  assert_eq(results.length(), 3)
  let mut n_sig = 0
  for r in results {
    if @src.MethylDiffResult::significant(r) {
      n_sig = n_sig + 1
    }
  }
  assert_eq(n_sig, 0)
  // Fisher test on identical counts -> p-value ~1.0.
  let p = @src.methylkit_fisher_test(cov1[0], cov2[0])
  assert_true(p > 0.99)
}
