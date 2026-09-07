///|
/// Tests for the LEfSe module (Linear Discriminant Analysis Effect Size).
///
/// Covers Kruskal-Wallis test, Wilcoxon rank-sum test, LDA score,
/// LEfSeInput construction, full LEfSe pipeline (lefse_run),
/// and discriminative feature identification from synthetic data.

// ===========================================================================
// 1. Kruskal-Wallis test
// ===========================================================================

///|
test "lefse_kruskal_wallis_identical_groups_h_near_zero" {
  // Three groups with identical distributions => H should be near 0
  let g1 = [5.0, 6.0, 7.0, 8.0, 9.0]
  let g2 = [5.0, 6.0, 7.0, 8.0, 9.0]
  let g3 = [5.0, 6.0, 7.0, 8.0, 9.0]
  let (h, p) = @src.lefse_kruskal_wallis([g1, g2, g3])
  assert_true(h < 1.0)
  assert_true(p > 0.05)
}

///|
test "lefse_kruskal_wallis_identical_groups_high_pvalue" {
  // Identical groups => p-value should be high (not significant)
  let g1 = [10.0, 12.0, 11.0, 13.0, 10.0]
  let g2 = [11.0, 10.0, 12.0, 11.0, 13.0]
  let (_h, p) = @src.lefse_kruskal_wallis([g1, g2])
  assert_true(p > 0.1)
}

///|
test "lefse_kruskal_wallis_different_groups_high_h" {
  // Groups with clearly different locations => H should be high
  let g1 = [10.0, 12.0, 11.0, 13.0, 10.0]
  let g2 = [50.0, 52.0, 51.0, 53.0, 50.0]
  let g3 = [90.0, 92.0, 91.0, 93.0, 90.0]
  let (h, p) = @src.lefse_kruskal_wallis([g1, g2, g3])
  assert_true(h > 10.0)
  assert_true(p < 0.05)
}

///|
test "lefse_kruskal_wallis_different_groups_low_pvalue" {
  // Clearly different groups => p-value should be very low
  let g1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let g2 = [100.0, 101.0, 102.0, 103.0, 104.0]
  let (_h, p) = @src.lefse_kruskal_wallis([g1, g2])
  assert_true(p < 0.05)
}

///|
test "lefse_kruskal_wallis_single_group_returns_zero" {
  let g1 = [1.0, 2.0, 3.0]
  let (h, p) = @src.lefse_kruskal_wallis([g1])
  assert_true(h == 0.0)
  assert_true(p == 1.0)
}

///|
test "lefse_kruskal_wallis_empty_group_returns_zero" {
  let g1 = [1.0, 2.0, 3.0]
  let g2 : Array[Double] = []
  let (h, p) = @src.lefse_kruskal_wallis([g1, g2])
  assert_true(h == 0.0)
  assert_true(p == 1.0)
}

// ===========================================================================
// 2. Wilcoxon rank-sum test
// ===========================================================================

///|
test "lefse_wilcoxon_identical_samples_high_pvalue" {
  // Identical samples => p-value should be high
  let x = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let (_u, p) = @src.lefse_wilcoxon_rank_sum(x, y)
  assert_true(p > 0.05)
}

///|
test "lefse_wilcoxon_identical_samples_u_near_expected" {
  // For identical samples, U should be near n1*n2/2
  let x = [5.0, 6.0, 7.0, 8.0, 9.0]
  let y = [5.0, 6.0, 7.0, 8.0, 9.0]
  let (u, _p) = @src.lefse_wilcoxon_rank_sum(x, y)
  let expected_u = 5.0 * 5.0 / 2.0 // 12.5
  assert_true((u - expected_u).abs() < 3.0)
}

///|
test "lefse_wilcoxon_different_samples_low_pvalue" {
  // Clearly different samples => p-value should be low
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let y = [50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0]
  let (_u, p) = @src.lefse_wilcoxon_rank_sum(x, y)
  assert_true(p < 0.01)
}

///|
test "lefse_wilcoxon_different_samples_u_small" {
  // When x << y, U should be small (close to 0)
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [100.0, 101.0, 102.0, 103.0, 104.0]
  let (u, _p) = @src.lefse_wilcoxon_rank_sum(x, y)
  assert_true(u < 3.0)
}

///|
test "lefse_wilcoxon_empty_returns_default" {
  let x : Array[Double] = []
  let y = [1.0, 2.0, 3.0]
  let (u, p) = @src.lefse_wilcoxon_rank_sum(x, y)
  assert_true(u == 0.0)
  assert_true(p == 1.0)
}

// ===========================================================================
// 3. LDA score
// ===========================================================================

///|
test "lefse_lda_score_positive_when_best_class_higher" {
  // Best class (index 1) has higher mean => positive LDA score
  let means = [10.0, 100.0, 20.0]
  let lda = @src.lefse_lda_score(means, 1)
  assert_true(lda > 0.0)
}

///|
test "lefse_lda_score_magnitude_for_large_difference" {
  // Large difference => LDA score should be substantially > 0
  let means = [10.0, 1000.0]
  let lda = @src.lefse_lda_score(means, 1)
  assert_true(lda > 1.0)
}

///|
test "lefse_lda_score_zero_when_equal_means" {
  // All means equal => LDA score should be 0 (log10(1) = 0)
  let means = [50.0, 50.0, 50.0]
  let lda = @src.lefse_lda_score(means, 0)
  assert_true(lda.abs() < 0.0001)
}

///|
test "lefse_lda_score_two_classes" {
  // Two classes: best class has mean 200, other has mean 50
  // LDA = log10(200/50) = log10(4) ≈ 0.602
  let means = [50.0, 200.0]
  let lda = @src.lefse_lda_score(means, 1)
  assert_true(lda > 0.5)
  assert_true(lda < 0.7)
}

///|
test "lefse_lda_score_invalid_index_returns_zero" {
  let means = [10.0, 20.0]
  let lda = @src.lefse_lda_score(means, 5)
  assert_true(lda == 0.0)
}

///|
test "lefse_lda_score_single_class_returns_zero" {
  let means = [10.0]
  let lda = @src.lefse_lda_score(means, 0)
  assert_true(lda == 0.0)
}

// ===========================================================================
// 4. LEfSeInput construction
// ===========================================================================

///|
test "lefse_input_new_creates_proper_input" {
  let features = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let feature_names = ["TaxonA", "TaxonB"]
  let sample_names = ["S1", "S2", "S3"]
  let classes = ["A", "A", "B"]
  let subclasses = ["", "", ""]
  let input = @src.LEfSeInput::new(
    features, feature_names, sample_names, classes, subclasses,
  )
  assert_eq(input.feature_names.length(), 2)
  assert_eq(input.sample_names.length(), 3)
  assert_eq(input.classes.length(), 3)
  assert_eq(input.features.length(), 2)
  assert_eq(input.features[0].length(), 3)
}

// ===========================================================================
// 5. LEfSe pipeline (lefse_run)
// ===========================================================================

///|
test "lefse_run_returns_valid_pvalues" {
  let input = @src.lefse_sample_data()
  let output = @src.lefse_run(input)
  for r in output.results() {
    assert_true(r.kw_pvalue() >= 0.0)
    assert_true(r.kw_pvalue() <= 1.0)
  }
}

///|
test "lefse_run_returns_valid_lda_scores" {
  let input = @src.lefse_sample_data()
  let output = @src.lefse_run(input)
  for r in output.results() {
    // LDA scores should be finite (not NaN, not infinite)
    let lda = r.lda_score()
    assert_true(lda == lda) // not NaN
    assert_true(lda < 1000.0)
    assert_true(lda > -1000.0)
  }
}

///|
test "lefse_run_n_significant_consistent" {
  let input = @src.lefse_sample_data()
  let output = @src.lefse_run(input)
  // Count significant results manually and compare
  let mut count = 0
  for r in output.results() {
    if r.is_significant() {
      count = count + 1
    }
  }
  assert_eq(output.n_significant(), count)
}

///|
test "lefse_run_results_count_equals_features" {
  let input = @src.lefse_sample_data()
  let output = @src.lefse_run(input)
  assert_eq(output.results().length(), input.feature_names.length())
}

///|
test "lefse_run_identifies_discriminative_features" {
  // Build synthetic data where features 0 and 1 are clearly differential
  // between two classes, while features 2 and 3 are not.
  let feature_names = ["DiffA", "DiffB", "SameC", "SameD"]
  let sample_names = [
    "H1", "H2", "H3", "H4", "H5", "D1", "D2", "D3", "D4", "D5",
  ]
  let classes = [
    "Healthy", "Healthy", "Healthy", "Healthy", "Healthy",
    "Disease", "Disease", "Disease", "Disease", "Disease",
  ]
  let subclasses : Array[String] = [
    "", "", "", "", "", "", "", "", "", "",
  ]
  let features = [
    // DiffA: high in Disease
    [10.0, 12.0, 11.0, 13.0, 10.0, 100.0, 102.0, 101.0, 103.0, 100.0],
    // DiffB: high in Disease
    [20.0, 22.0, 21.0, 23.0, 20.0, 200.0, 202.0, 201.0, 203.0, 200.0],
    // SameC: similar across classes
    [50.0, 52.0, 51.0, 53.0, 50.0, 51.0, 50.0, 52.0, 49.0, 51.0],
    // SameD: similar across classes
    [30.0, 32.0, 31.0, 33.0, 30.0, 31.0, 30.0, 32.0, 29.0, 31.0],
  ]
  let input = @src.LEfSeInput::new(
    features, feature_names, sample_names, classes, subclasses,
  )
  let output = @src.lefse_run(input, alpha=0.05, lda_threshold=1.0)

  // Differential features should have higher |LDA| than non-differential
  let lda_diff_a = output.results()[0].lda_score().abs()
  let lda_diff_b = output.results()[1].lda_score().abs()
  let lda_same_c = output.results()[2].lda_score().abs()
  let lda_same_d = output.results()[3].lda_score().abs()

  assert_true(lda_diff_a > lda_same_c)
  assert_true(lda_diff_a > lda_same_d)
  assert_true(lda_diff_b > lda_same_c)
  assert_true(lda_diff_b > lda_same_d)
}

///|
test "lefse_run_significant_features_have_low_pvalue_and_high_lda" {
  let input = @src.lefse_sample_data()
  let output = @src.lefse_run(input)
  for r in output.results() {
    if r.is_significant() {
      // Significant => KW p-value < 0.05
      assert_true(r.kw_pvalue() < 0.05)
      // Significant => |LDA| > 2.0 (default threshold)
      assert_true(r.lda_score().abs() > 2.0)
    }
  }
}

///|
test "lefse_run_nonsignificant_features_fail_threshold" {
  let input = @src.lefse_sample_data()
  let output = @src.lefse_run(input)
  for r in output.results() {
    if !r.is_significant() {
      // Non-significant => either p-value >= 0.05 or |LDA| <= 2.0
      let p_fail = r.kw_pvalue() >= 0.05
      let lda_fail = r.lda_score().abs() <= 2.0
      assert_true(p_fail || lda_fail)
    }
  }
}

// ===========================================================================
// 6. LEfSeResult accessors
// ===========================================================================

///|
test "lefse_result_accessors" {
  let r = @src.LEfSeResult::new("Bacteroides", "Disease", 0.001, 3.5, true)
  assert_eq(r.feature_name(), "Bacteroides")
  assert_eq(r.class_name(), "Disease")
  assert_true(r.kw_pvalue() == 0.001)
  assert_true(r.lda_score() == 3.5)
  assert_true(r.is_significant())
}

// ===========================================================================
// 7. LEfSeOutput accessors
// ===========================================================================

///|
test "lefse_output_accessors" {
  let r1 = @src.LEfSeResult::new("A", "C1", 0.01, 2.5, true)
  let r2 = @src.LEfSeResult::new("B", "C2", 0.5, 0.0, false)
  let output = @src.LEfSeOutput::new([r1, r2], 1)
  assert_eq(output.results().length(), 2)
  assert_eq(output.n_significant(), 1)
}

// ===========================================================================
// 8. Sample data validation
// ===========================================================================

///|
test "lefse_sample_data_structure" {
  let input = @src.lefse_sample_data()
  assert_eq(input.feature_names.length(), 6)
  assert_eq(input.sample_names.length(), 10)
  assert_eq(input.classes.length(), 10)
  assert_eq(input.features.length(), 6)
}

///|
test "lefse_sample_data_classes" {
  let input = @src.lefse_sample_data()
  let mut n_healthy = 0
  let mut n_disease = 0
  for c in input.classes {
    if c == "Healthy" {
      n_healthy = n_healthy + 1
    } else if c == "Disease" {
      n_disease = n_disease + 1
    }
  }
  assert_eq(n_healthy, 5)
  assert_eq(n_disease, 5)
}

///|
test "lefse_sample_data_differential_features_elevated" {
  let input = @src.lefse_sample_data()
  // Feature 0 (Bacteroides): Healthy ~120, Disease ~250
  let mut h_sum = 0.0
  let mut d_sum = 0.0
  let mut i = 0
  while i < 10 {
    if input.classes[i] == "Healthy" {
      h_sum = h_sum + input.features[0][i]
    } else {
      d_sum = d_sum + input.features[0][i]
    }
    i = i + 1
  }
  assert_true(d_sum > h_sum)
}

// ===========================================================================
// 9. Result summary
// ===========================================================================

///|
test "lefse_result_summary_contains_info" {
  let input = @src.lefse_sample_data()
  let output = @src.lefse_run(input)
  let summary = @src.lefse_result_summary(output)
  assert_true(summary.contains("LEfSe"))
  assert_true(summary.contains("Bacteroides"))
}

// ===========================================================================
// 10. Edge cases
// ===========================================================================

///|
test "lefse_run_empty_features" {
  let input = @src.LEfSeInput::new([], [], [], [], [])
  let output = @src.lefse_run(input)
  assert_eq(output.results().length(), 0)
  assert_eq(output.n_significant(), 0)
}

///|
test "lefse_run_single_class_returns_empty" {
  let features = [[1.0, 2.0, 3.0]]
  let feature_names = ["A"]
  let sample_names = ["S1", "S2", "S3"]
  let classes = ["X", "X", "X"]
  let subclasses = ["", "", ""]
  let input = @src.LEfSeInput::new(
    features, feature_names, sample_names, classes, subclasses,
  )
  let output = @src.lefse_run(input)
  assert_eq(output.results().length(), 0)
  assert_eq(output.n_significant(), 0)
}
