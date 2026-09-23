///|
/// Tests for the DMRcate module.
/// Differentially Methylated Region calling with spatial kernel smoothing.

///|
test "dmrcate_input_new" {
  let beta = [[0.5, 0.6], [0.7, 0.8]]
  let pos = [1000, 2000]
  let chr = ["chr1", "chr1"]
  let design = [[1.0, 1.0], [1.0, 0.0]]
  let names = ["s1", "s2"]
  let input = @src.DMRcateInput::new(beta, pos, chr, design, names)
  assert_eq(@src.DMRcateInput::n_cpgs(input), 2)
  assert_eq(@src.DMRcateInput::n_samples(input), 2)
}

///|
test "dmrcate_spatial_kernel_distance_zero" {
  // At distance 0, the exponential kernel should return 1.0
  let w = @src.dmrcate_spatial_kernel(0.0, 1000.0)
  assert_true((w - 1.0).abs() < 1.0e-10)
}

///|
test "dmrcate_spatial_kernel_decreases_with_distance" {
  let lambda = 1000.0
  let w0 = @src.dmrcate_spatial_kernel(0.0, lambda)
  let w100 = @src.dmrcate_spatial_kernel(100.0, lambda)
  let w500 = @src.dmrcate_spatial_kernel(500.0, lambda)
  let w1000 = @src.dmrcate_spatial_kernel(1000.0, lambda)
  // Kernel should decrease as distance increases
  assert_true(w0 > w100)
  assert_true(w100 > w500)
  assert_true(w500 > w1000)
  // All weights should be positive
  assert_true(w0 > 0.0)
  assert_true(w1000 > 0.0)
}

///|
test "dmrcate_spatial_kernel_lambda_zero" {
  // lambda <= 0 should return 0.0
  let w = @src.dmrcate_spatial_kernel(100.0, 0.0)
  assert_true((w - 0.0).abs() < 1.0e-10)
  let w2 = @src.dmrcate_spatial_kernel(100.0, -1.0)
  assert_true((w2 - 0.0).abs() < 1.0e-10)
}

///|
test "dmrcate_cluster_cpgs_nearby" {
  // Three CpGs on chr1 within 500 bp of each other -> one cluster
  let cpgs = [
    @src.CpGResult::new("CpG1", "chr1", 1000, 0.5, 3.0, 0.01),
    @src.CpGResult::new("CpG2", "chr1", 1300, 0.6, 2.5, 0.02),
    @src.CpGResult::new("CpG3", "chr1", 1600, 0.7, 2.0, 0.03),
  ]
  let clusters = @src.dmrcate_cluster_cpgs(cpgs, 500)
  // All three are within 500 bp of the previous one -> one cluster
  assert_eq(clusters.length(), 1)
  assert_eq(clusters[0].length(), 3)
}

///|
test "dmrcate_cluster_cpgs_separate_chromosomes" {
  // CpGs on different chromosomes -> separate clusters
  let cpgs = [
    @src.CpGResult::new("CpG1", "chr1", 1000, 0.5, 3.0, 0.01),
    @src.CpGResult::new("CpG2", "chr2", 1000, 0.6, 2.5, 0.02),
  ]
  let clusters = @src.dmrcate_cluster_cpgs(cpgs, 500)
  assert_eq(clusters.length(), 2)
  assert_eq(clusters[0].length(), 1)
  assert_eq(clusters[1].length(), 1)
}

///|
test "dmrcate_cluster_cpgs_beyond_lambda" {
  // CpGs beyond lambda bp apart -> separate clusters
  let cpgs = [
    @src.CpGResult::new("CpG1", "chr1", 1000, 0.5, 3.0, 0.01),
    @src.CpGResult::new("CpG2", "chr1", 3000, 0.6, 2.5, 0.02),
  ]
  let clusters = @src.dmrcate_cluster_cpgs(cpgs, 500)
  // Distance is 2000 > 500, so separate clusters
  assert_eq(clusters.length(), 2)
}

///|
test "dmrcate_cluster_cpgs_empty" {
  let clusters = @src.dmrcate_cluster_cpgs([], 500)
  assert_eq(clusters.length(), 0)
}

///|
test "dmrcate_test_returns_results" {
  let input = @src.dmrcate_example_data()
  let result = @src.dmrcate_test(input)
  // Should return CpG results and DMR arrays
  assert_true(@src.DMRcateResult::n_cpgs(result) > 0)
  // CpG results should match input CpG count
  assert_eq(
    @src.DMRcateResult::n_cpgs(result),
    @src.DMRcateInput::n_cpgs(input),
  )
}

///|
test "dmrcate_test_dmr_start_le_end" {
  let input = @src.dmrcate_example_data()
  let result = @src.dmrcate_test(input)
  let dmrs = @src.DMRcateResult::dmrs(result)
  // Every DMR should have start <= end
  let mut i = 0
  while i < dmrs.length() {
    assert_true(@src.DMR::start(dmrs[i]) <= @src.DMR::end_(dmrs[i]))
    i = i + 1
  }
}

///|
test "dmrcate_test_dmr_n_cpgs_at_least_one" {
  let input = @src.dmrcate_example_data()
  let result = @src.dmrcate_test(input)
  let dmrs = @src.DMRcateResult::dmrs(result)
  // Every DMR should contain at least one CpG
  let mut i = 0
  while i < dmrs.length() {
    assert_true(@src.DMR::n_cpgs(dmrs[i]) >= 1)
    i = i + 1
  }
}

///|
test "dmrcate_test_cpg_result_fields" {
  let input = @src.dmrcate_example_data()
  let result = @src.dmrcate_test(input)
  let cpg_results = @src.DMRcateResult::cpg_results(result)
  // Check that CpG results have valid fields
  assert_true(cpg_results.length() > 0)
  let first = cpg_results[0]
  assert_true(@src.CpGResult::cpg_name(first).length() > 0)
  assert_true(@src.CpGResult::chr(first).length() > 0)
  assert_true(@src.CpGResult::pos(first) > 0)
  // P-value should be in [0, 1]
  assert_true(@src.CpGResult::p_value(first) >= 0.0)
  assert_true(@src.CpGResult::p_value(first) <= 1.0)
}

///|
test "dmrcate_test_custom_params" {
  let input = @src.dmrcate_example_data()
  // Run with custom lambda, c_scale, and fdr_threshold
  let result = @src.dmrcate_test(
    input,
    lambda=500,
    c_scale=1.5,
    fdr_threshold=0.1,
  )
  assert_true(@src.DMRcateResult::n_cpgs(result) > 0)
}

///|
test "dmrcate_test_synthetic_data" {
  // Create a small synthetic dataset with a clear DMR signal
  // 6 samples: 3 cases + 3 controls, 10 CpGs on chr1
  let n_cpgs = 10
  let n_samples = 6
  let beta_values : Array[Array[Double]] = Array(capacity=0)
  let positions : Array[Int] = Array(capacity=0)
  let chr : Array[String] = Array(capacity=0)

  let mut i = 0
  while i < n_cpgs {
    chr.push("chr1")
    positions.push(10000 + i * 200)
    let vals : Array[Double] = Array(capacity=0)
    let mut s = 0
    while s < n_samples {
      // Cases (s < 3): elevated methylation for CpGs 3-7
      // Controls (s >= 3): baseline methylation
      let mut beta = 0.4
      if s < 3 && i >= 3 && i <= 7 {
        beta = beta + 0.35
      }
      // Small deterministic noise
      beta = beta + ((i * 3 + s * 7) % 5).to_double() * 0.01 - 0.02
      if beta < 0.0 {
        beta = 0.0
      }
      if beta > 1.0 {
        beta = 1.0
      }
      vals.push(beta)
      s = s + 1
    }
    beta_values.push(vals)
    i = i + 1
  }

  // Design: intercept + group
  let intercept : Array[Double] = Array(capacity=0)
  let group : Array[Double] = Array(capacity=0)
  let mut s = 0
  while s < n_samples {
    intercept.push(1.0)
    group.push(if s < 3 { 1.0 } else { 0.0 })
    s = s + 1
  }
  let design : Array[Array[Double]] = [intercept, group]
  let sample_names = ["c1", "c2", "c3", "n1", "n2", "n3"]

  let input = @src.DMRcateInput::new(
    beta_values, positions, chr, design, sample_names,
  )
  let result = @src.dmrcate_test(input)

  // Should have results for all 10 CpGs
  assert_eq(@src.DMRcateResult::n_cpgs(result), n_cpgs)

  // DMRs should have valid structure
  let dmrs = @src.DMRcateResult::dmrs(result)
  let mut j = 0
  while j < dmrs.length() {
    assert_true(@src.DMR::start(dmrs[j]) <= @src.DMR::end_(dmrs[j]))
    assert_true(@src.DMR::n_cpgs(dmrs[j]) >= 1)
    j = j + 1
  }
}

///|
test "dmrcate_example_data_structure" {
  let input = @src.dmrcate_example_data()
  assert_eq(@src.DMRcateInput::n_cpgs(input), 50)
  assert_eq(@src.DMRcateInput::n_samples(input), 8)
}

///|
test "dmrcate_spatial_autocorrelation" {
  let cpgs = [
    @src.CpGResult::new("CpG1", "chr1", 1000, 0.5, 3.0, 0.01),
    @src.CpGResult::new("CpG2", "chr1", 1200, 0.6, 2.5, 0.02),
  ]
  let corr = @src.dmrcate_spatial_autocorrelation(cpgs, 1000.0)
  // Autocorrelation of positively correlated t-stats should be positive
  assert_true(corr > 0.0)
}

///|
test "dmrcate_filter_dmrs" {
  let dmr1 = @src.DMR::new("chr1", 1000, 2000, 3, 3.5, 0.3, 0.01, [
    "a", "b", "c",
  ])
  let dmr2 = @src.DMR::new("chr1", 3000, 4000, 1, 1.5, 0.1, 0.5, ["d"])
  let filtered = @src.dmrcate_filter_dmrs([dmr1, dmr2], 2, 0.05)
  // Only dmr1 has n_cpgs >= 2 and p_value <= 0.05
  assert_eq(filtered.length(), 1)
  assert_eq(@src.DMR::n_cpgs(filtered[0]), 3)
}

///|
test "dmrcate_dmr_cpg_names" {
  let dmr1 = @src.DMR::new("chr1", 1000, 2000, 2, 3.0, 0.2, 0.01, [
    "CpG1", "CpG2",
  ])
  let dmr2 = @src.DMR::new("chr1", 3000, 4000, 1, 2.0, 0.1, 0.02, ["CpG3"])
  let names = @src.dmrcate_dmr_cpg_names([dmr1, dmr2])
  assert_eq(names.length(), 3)
  assert_eq(names[0], "CpG1")
  assert_eq(names[1], "CpG2")
  assert_eq(names[2], "CpG3")
}

///|
test "dmrcate_summary" {
  let input = @src.dmrcate_example_data()
  let result = @src.dmrcate_test(input)
  let s = @src.dmrcate_summary(result)
  assert_true(s.contains("DMRcate Analysis Summary"))
  assert_true(s.contains("CpG sites tested: 50"))
}

///|
test "dmrcate_results_to_string" {
  let input = @src.dmrcate_example_data()
  let result = @src.dmrcate_test(input)
  let s = @src.dmrcate_results_to_string(result)
  assert_true(s.contains("DMRcate DMRs"))
}

///|
test "dmrcate_dmr_length" {
  let dmr = @src.DMR::new("chr1", 1000, 5000, 5, 3.0, 0.2, 0.01, [
    "a", "b", "c", "d", "e",
  ])
  assert_eq(@src.DMR::length(dmr), 4000)
}

///|
test "dmrcate_dmr_to_string" {
  let dmr = @src.DMR::new("chr1", 1000, 5000, 5, 3.0, 0.2, 0.01, [
    "a", "b", "c", "d", "e",
  ])
  let s = @src.DMR::to_string(dmr)
  assert_true(s.contains("chr1"))
  assert_true(s.contains("1000"))
  assert_true(s.contains("5000"))
  assert_true(s.contains("n_cpgs=5"))
}
