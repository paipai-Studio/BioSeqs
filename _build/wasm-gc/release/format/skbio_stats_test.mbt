///|
/// Blackbox tests for scikit-bio Stats module.

// ============================================================
// subsample_counts tests
// ============================================================

///|
test "subsample_counts_with_replacement_sum" {
  let counts = [10.0, 20.0, 30.0, 40.0]
  let result = subsample_counts(counts, 50, true, 42)
  let mut total = 0.0
  for x in result {
    total = total + x
  }
  assert_eq(total.to_int(), 50)
}

///|
test "subsample_counts_without_replacement_sum" {
  let counts = [10.0, 20.0, 30.0, 40.0]
  let result = subsample_counts(counts, 50, false, 42)
  let mut total = 0.0
  for x in result {
    total = total + x
  }
  assert_eq(total.to_int(), 50)
}

///|
test "subsample_counts_without_replacement_exceeds_total" {
  let counts = [1.0, 2.0, 3.0]
  let result = subsample_counts(counts, 100, false, 42)
  let mut total = 0.0
  for x in result {
    total = total + x
  }
  assert_eq(total.to_int(), 0)
}

///|
test "subsample_counts_negative_n" {
  let counts = [10.0, 20.0, 30.0]
  let result = subsample_counts(counts, -1, true, 42)
  let mut total = 0.0
  for x in result {
    total = total + x
  }
  assert_eq(total.to_int(), 0)
}

///|
test "subsample_counts_zero_n" {
  let counts = [10.0, 20.0, 30.0]
  let result = subsample_counts(counts, 0, true, 42)
  let mut total = 0.0
  for x in result {
    total = total + x
  }
  assert_eq(total.to_int(), 0)
}

// ============================================================
// p-value adjustment tests
// ============================================================

///|
test "p_adjust_bonferroni_basic" {
  let pvals = [0.01, 0.02, 0.03, 0.5, 0.005, 0.001]
  let adj = p_adjust_bonferroni(pvals)
  assert_true((adj[0] - 0.06).abs() < 0.0001)
  assert_true((adj[1] - 0.12).abs() < 0.0001)
  assert_true((adj[2] - 0.18).abs() < 0.0001)
  assert_true((adj[3] - 1.0).abs() < 0.0001)
  assert_true((adj[4] - 0.03).abs() < 0.0001)
  assert_true((adj[5] - 0.006).abs() < 0.0001)
}

///|
test "p_adjust_bonferroni_empty" {
  let pvals : Array[Double] = Array::new()
  let adj = p_adjust_bonferroni(pvals)
  assert_eq(adj.length(), 0)
}

///|
test "p_adjust_holm_basic" {
  let pvals = [0.01, 0.02, 0.03, 0.5, 0.005, 0.001]
  let adj = p_adjust_holm(pvals)
  assert_true((adj[0] - 0.04).abs() < 0.0001)
  assert_true((adj[1] - 0.06).abs() < 0.0001)
  assert_true((adj[2] - 0.06).abs() < 0.0001)
  assert_true((adj[3] - 0.5).abs() < 0.0001)
  assert_true((adj[4] - 0.025).abs() < 0.0001)
  assert_true((adj[5] - 0.006).abs() < 0.0001)
}

///|
test "p_adjust_fdr_bh_basic" {
  let pvals = [0.01, 0.02, 0.03, 0.5, 0.005, 0.001]
  let adj = p_adjust_fdr_bh(pvals)
  assert_true((adj[0] - 0.02).abs() < 0.0001)
  assert_true((adj[1] - 0.03).abs() < 0.0001)
  assert_true((adj[2] - 0.036).abs() < 0.0001)
  assert_true((adj[3] - 0.5).abs() < 0.0001)
  assert_true((adj[4] - 0.015).abs() < 0.0001)
  assert_true((adj[5] - 0.006).abs() < 0.0001)
}

// ============================================================
// permutation test tests
// ============================================================

///|
test "permutation_test_mean_diff" {
  let s1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let s2 = [2.0, 3.0, 4.0, 5.0, 6.0]
  let result = permutation_test(s1, s2, diff_of_means, 1000, 42)
  assert_eq(result.permutations, 1000)
  assert_true(result.p_value >= 0.0)
  assert_true(result.p_value <= 1.0)
}

///|
test "permutation_test_zero_permutations" {
  let s1 = [1.0, 2.0, 3.0]
  let s2 = [4.0, 5.0, 6.0]
  let result = permutation_test(s1, s2, diff_of_means, 0, 42)
  assert_eq(result.permutations, 0)
}

///|
test "diff_of_means_basic" {
  let s1 = [1.0, 2.0, 3.0]
  let s2 = [4.0, 5.0, 6.0]
  let d = diff_of_means(s1, s2)
  assert_true((d - -3.0).abs() < 0.0001)
}

// ============================================================
// DistanceMatrix tests
// ============================================================

///|
test "distance_matrix_from_2d" {
  let matrix = [[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]]
  let ids = ["A", "B", "C"]
  let dm = DistanceMatrix::from_2d(matrix, ids)
  assert_eq(dm.size, 3)
  assert_true((dm.get(0, 1) - 1.0).abs() < 0.0001)
  assert_true((dm.get(1, 2) - 3.0).abs() < 0.0001)
  assert_true(dm.get(0, 0).abs() < 0.0001)
}

// ============================================================
// ANOSIM tests
// ============================================================

///|
test "anosim_basic" {
  let matrix = [
    [0.0, 1.0, 2.0, 8.0, 9.0, 10.0],
    [1.0, 0.0, 3.0, 7.0, 8.0, 9.0],
    [2.0, 3.0, 0.0, 6.0, 7.0, 8.0],
    [8.0, 7.0, 6.0, 0.0, 1.0, 2.0],
    [9.0, 8.0, 7.0, 1.0, 0.0, 3.0],
    [10.0, 9.0, 8.0, 2.0, 3.0, 0.0],
  ]
  let ids = ["A", "B", "C", "D", "E", "F"]
  let dm = DistanceMatrix::from_2d(matrix, ids)
  let groups = [0, 0, 0, 1, 1, 1]
  let result = anosim(dm, groups, 100, 42)
  assert_eq(result.num_groups, 2)
  assert_eq(result.sample_size, 6)
  assert_eq(result.permutations, 100)
  assert_true(result.r_statistic >= -1.0)
  assert_true(result.r_statistic <= 1.0)
  assert_true(result.p_value >= 0.0)
  assert_true(result.p_value <= 1.0)
}

// ============================================================
// PERMANOVA tests
// ============================================================

///|
test "permanova_basic" {
  let matrix = [
    [0.0, 1.0, 2.0, 8.0, 9.0, 10.0],
    [1.0, 0.0, 3.0, 7.0, 8.0, 9.0],
    [2.0, 3.0, 0.0, 6.0, 7.0, 8.0],
    [8.0, 7.0, 6.0, 0.0, 1.0, 2.0],
    [9.0, 8.0, 7.0, 1.0, 0.0, 3.0],
    [10.0, 9.0, 8.0, 2.0, 3.0, 0.0],
  ]
  let ids = ["A", "B", "C", "D", "E", "F"]
  let dm = DistanceMatrix::from_2d(matrix, ids)
  let groups = [0, 0, 0, 1, 1, 1]
  let result = permanova(dm, groups, 100, 42)
  assert_eq(result.num_groups, 2)
  assert_eq(result.sample_size, 6)
  assert_eq(result.permutations, 100)
  assert_eq(result.df_between, 1)
  assert_eq(result.df_within, 4)
  assert_true(result.f_statistic >= 0.0)
  assert_true(result.p_value >= 0.0)
  assert_true(result.p_value <= 1.0)
  assert_true(
    (result.ss_total - (result.ss_between + result.ss_within)).abs() < 0.0001,
  )
}

// ============================================================
// confidence_bound tests
// ============================================================

///|
test "confidence_bound_basic" {
  let (lower, upper) = confidence_bound(5.0, 2.0, 100, 0.95)
  assert_true((lower - 4.608007).abs() < 0.0001)
  assert_true((upper - 5.391993).abs() < 0.0001)
}

///|
test "confidence_bound_small_n" {
  let (lower, upper) = confidence_bound(5.0, 2.0, 1, 0.95)
  assert_true((lower - 5.0).abs() < 0.0001)
  assert_true((upper - 5.0).abs() < 0.0001)
}

///|
test "confidence_bound_zero_std" {
  let (lower, upper) = confidence_bound(5.0, 0.0, 100, 0.95)
  assert_true((lower - 5.0).abs() < 0.0001)
  assert_true((upper - 5.0).abs() < 0.0001)
}
