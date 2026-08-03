///|
/// Tests for preprocessCore module:
/// Quantile normalization, invariant set normalization,
/// log transform, background correction, and utilities.

// ============================================================================
// Config and basic structure tests
// ============================================================================

test "pc_quantile_config_new" {
  let config = @src.QuantileNormConfig::new()
  assert_eq(config.copy, true)
  assert_eq(config.keep_names, true)
  assert_eq(config.method, "quantile")
}

test "pc_invariant_set_result_creation" {
  let dummy : Array[Array[Double]] = [[1.0, 2.0], [3.0, 4.0]]
  let res = @src.InvariantSetResult::new(dummy, [0, 1], 0, [0.0, 1.0])
  assert_eq(res.normalized_data.length(), 2)
  assert_eq(res.invariant_set_indices.length(), 2)
  assert_eq(res.reference_index, 0)
  assert_eq(res.fit_coefficients.length(), 2)
}

// ============================================================================
// Helper utilities
// ============================================================================

test "pc_interp_linear_basic" {
  let xs = [0.0, 1.0, 2.0, 3.0, 4.0]
  let ys = [0.0, 10.0, 20.0, 30.0, 40.0]
  // Interpolate at 0.5 -> 5.0
  let v = @src.pc_interp_linear(xs, ys, 0.5)
  assert_true((v - 5.0).abs() < 0.0001)
  // Exact points
  let v1 = @src.pc_interp_linear(xs, ys, 2.0)
  assert_true((v1 - 20.0).abs() < 0.0001)
  // Below min
  let v2 = @src.pc_interp_linear(xs, ys, -1.0)
  assert_true((v2 - 0.0).abs() < 0.0001)
  // Above max
  let v3 = @src.pc_interp_linear(xs, ys, 10.0)
  assert_true((v3 - 40.0).abs() < 0.0001)
}

test "pc_interp_linear_empty" {
  let v = @src.pc_interp_linear([], [], 1.0)
  assert_eq(v, 0.0)
}

test "pc_interp_linear_single" {
  let v = @src.pc_interp_linear([5.0], [100.0], 0.0)
  assert_eq(v, 100.0)
}

// ============================================================================
// Quantile normalization core tests
// ============================================================================

test "pc_normalize_quantiles_identical_distributions" {
  // Two columns with same values -> should remain the same
  let matrix = [
    [1.0, 1.0],
    [2.0, 2.0],
    [3.0, 3.0],
    [4.0, 4.0],
    [5.0, 5.0],
  ]
  let res = @src.normalize_quantiles(matrix)
  assert_eq(res.length(), 5)
  assert_eq(res[0].length(), 2)
  // Each pair of columns should be equal after quantile norm
  for row in 0..<5 {
    assert_true((res[row][0] - res[row][1]).abs() < 0.0001)
  }
}

test "pc_normalize_quantiles_swapped_values" {
  // Two columns with reversed order
  // Col1: [1, 2, 3, 4], Col2: [4, 3, 2, 1]
  // After normalization, both should have means of sorted values
  let matrix = [
    [1.0, 4.0],
    [2.0, 3.0],
    [3.0, 2.0],
    [4.0, 1.0],
  ]
  let res = @src.normalize_quantiles(matrix)
  // Check distributions: for each column, after sorting, they should match
  // sorted(col1_normalized) == sorted(col2_normalized) == 2.5, 2.5, 2.5, 2.5
  // Actually: mean(sorted_col1[i], sorted_col2[i])
  // sorted col1 = [1, 2, 3, 4], sorted col2 = [1, 2, 3, 4]
  // means = [1.0, 2.0, 3.0, 4.0]
  // Then put back in original order.
  // Col1 original order is already sorted so: [1.0, 2.0, 3.0, 4.0]
  // Col2 original order is 4,3,2,1 so: [4.0, 3.0, 2.0, 1.0]
  // BUT actually in our case, means should be (1+1)/2, (2+2)/2, etc. = 1,2,3,4
  assert_true((res[0][0] - 1.0).abs() < 0.0001)
  assert_true((res[3][0] - 4.0).abs() < 0.0001)
  // Col 2 after: row 0 had 4 (largest) -> mean of i=3 = 4.0
  assert_true((res[0][1] - 4.0).abs() < 0.0001)
}

test "pc_normalize_quantiles_three_columns" {
  let matrix = [
    [1.0, 4.0, 7.0],
    [2.0, 5.0, 8.0],
    [3.0, 6.0, 9.0],
  ]
  let res = @src.normalize_quantiles(matrix)
  assert_eq(res.length(), 3)
  assert_eq(res[0].length(), 3)
  // Total sum should be preserved approximately
  let mut orig_sum = 0.0
  let mut new_sum = 0.0
  for r in 0..<3 {
    for c in 0..<3 {
      orig_sum = orig_sum + matrix[r][c]
      new_sum = new_sum + res[r][c]
    }
  }
  assert_true((orig_sum - new_sum).abs() < 0.1)
}

test "pc_normalize_quantiles_empty" {
  let res = @src.normalize_quantiles([])
  assert_eq(res.length(), 0)
}

test "pc_normalize_quantiles_one_row" {
  let matrix = [[1.0, 10.0, 100.0]]
  let res = @src.normalize_quantiles(matrix)
  assert_eq(res.length(), 1)
  assert_eq(res[0].length(), 3)
  // mean of [1, 10, 100] = 37, all entries become 37
  let mean_val = (1.0 + 10.0 + 100.0) / 3.0
  for c in 0..<3 {
    assert_true((res[0][c] - mean_val).abs() < 0.0001)
  }
}

// ============================================================================
// Target-based quantile normalization
// ============================================================================

test "pc_normalize_quantiles_target" {
  let matrix = [
    [1.0, 2.0],
    [3.0, 4.0],
    [5.0, 6.0],
  ]
  let target = [10.0, 20.0, 30.0]
  let res = @src.normalize_quantiles_use_target(matrix, target)
  assert_eq(res.length(), 3)
  assert_eq(res[0].length(), 2)
  // Sorted col1 was [1, 3, 5] -> new order uses target [10, 20, 30]
  // Original positions were correct order, so col1 becomes [10, 20, 30]
  assert_true((res[0][0] - 10.0).abs() < 0.0001)
  assert_true((res[1][0] - 20.0).abs() < 0.0001)
  assert_true((res[2][0] - 30.0).abs() < 0.0001)
}

test "pc_normalize_quantiles_determine_target" {
  let matrix = [
    [1.0, 3.0],
    [2.0, 2.0],
    [3.0, 1.0],
  ]
  let target = @src.normalize_quantiles_determine_target(matrix)
  assert_eq(target.length(), 3)
  // Each row of sorted matrix: sorted col1 = [1, 2, 3], col2 = [1, 2, 3]
  // means = [1.0, 2.0, 3.0]
  assert_true((target[0] - 1.0).abs() < 0.0001)
  assert_true((target[1] - 2.0).abs() < 0.0001)
  assert_true((target[2] - 3.0).abs() < 0.0001)
}

// ============================================================================
// Subset quantile normalization
// ============================================================================

test "pc_normalize_quantiles_subset" {
  let matrix = [
    [1.0, 8.0],
    [2.0, 7.0],
    [3.0, 6.0],
    [4.0, 5.0],
    [5.0, 4.0],
    [6.0, 3.0],
    [7.0, 2.0],
    [8.0, 1.0],
  ]
  let res = @src.normalize_quantiles_subset(matrix, [0, 2, 4, 6])
  assert_eq(res.length(), 8)
  assert_eq(res[0].length(), 2)
}

// ============================================================================
// Invariant set normalization
// ============================================================================

test "pc_find_invariant_set_perfect_correlation" {
  // Two arrays with perfect linear relationship
  let reference = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let tgt = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
  let inv = @src.find_invariant_set(reference, tgt, threshold=0.01)
  // Should find large invariant set when relationship is perfect
  assert_true(inv.length() >= 5)
}

test "pc_find_invariant_set_length_mismatch" {
  let inv = @src.find_invariant_set([1.0, 2.0], [1.0])
  assert_eq(inv.length(), 0)
}

test "pc_normalize_invariantset_basic" {
  // 3 samples x 5 rows
  let matrix = [
    [1.0, 2.0, 3.0],
    [2.0, 4.0, 6.0],
    [3.0, 6.0, 9.0],
    [4.0, 8.0, 12.0],
    [5.0, 10.0, 15.0],
  ]
  let res = @src.normalize_invariantset(matrix, reference_index=0, threshold=0.01)
  assert_eq(res.reference_index, 0)
  assert_eq(res.normalized_data.length(), 5)
  // Reference column should be unchanged
  for r in 0..<5 {
    assert_true((res.normalized_data[r][0] - matrix[r][0]).abs() < 0.001)
  }
}

// ============================================================================
// Log transform
// ============================================================================

test "pc_log2_transform_basic" {
  let matrix = [[0.0, 1.0], [3.0, 15.0]]
  let res = @src.log2_transform(matrix, offset=1.0)
  // log2(0 + 1) = 0
  assert_true((res[0][0] - 0.0).abs() < 0.0001)
  // log2(1 + 1) = 1
  assert_true((res[0][1] - 1.0).abs() < 0.0001)
  // log2(3 + 1) = 2
  assert_true((res[1][0] - 2.0).abs() < 0.0001)
  // log2(15 + 1) = 4
  assert_true((res[1][1] - 4.0).abs() < 0.0001)
}

test "pc_log2_transform_empty" {
  let res = @src.log2_transform([])
  assert_eq(res.length(), 0)
}

// ============================================================================
// Background correction
// ============================================================================

test "pc_background_correct_percentile_basic" {
  let matrix = [
    [1.0, 100.0],
    [2.0, 200.0],
    [3.0, 300.0],
    [4.0, 400.0],
    [5.0, 500.0],
    [6.0, 600.0],
    [7.0, 700.0],
    [8.0, 800.0],
    [9.0, 900.0],
    [10.0, 1000.0],
  ]
  let res = @src.background_correct_percentile(matrix, percentile=0.1)
  // No value should be negative
  for r in 0..<10 {
    for c in 0..<2 {
      assert_true(res[r][c] >= 0.0)
    }
  }
  // Corrected values should be <= original
  for r in 0..<10 {
    for c in 0..<2 {
      assert_true(res[r][c] <= matrix[r][c] + 0.0001)
    }
  }
}

// ============================================================================
// Median center
// ============================================================================

test "pc_median_center_columns_basic" {
  let matrix = [
    [1.0, 10.0],
    [2.0, 20.0],
    [3.0, 30.0],
    [4.0, 40.0],
    [5.0, 50.0],
  ]
  let res = @src.median_center_columns(matrix)
  // Column 1 median = 3, centered = [-2, -1, 0, 1, 2]
  assert_true((res[0][0] - (-2.0)).abs() < 0.0001)
  assert_true((res[2][0] - 0.0).abs() < 0.0001)
  assert_true((res[4][0] - 2.0).abs() < 0.0001)
  // Column 2 median = 30, centered = [-20, -10, 0, 10, 20]
  assert_true((res[0][1] - (-20.0)).abs() < 0.0001)
  assert_true((res[2][1] - 0.0).abs() < 0.0001)
  assert_true((res[4][1] - 20.0).abs() < 0.0001)
}

// ============================================================================
// Column summary
// ============================================================================

test "pc_column_summary_basic" {
  let matrix = [
    [1.0, 10.0],
    [2.0, 20.0],
    [3.0, 30.0],
    [4.0, 40.0],
    [5.0, 50.0],
  ]
  let res = @src.column_summary(matrix)
  assert_eq(res.length(), 2)
  // Each column result is [mean, median, sd, min, max]
  // Column 0: mean = 3, median = 3, min = 1, max = 5
  assert_true((res[0][0] - 3.0).abs() < 0.0001)
  assert_true((res[0][1] - 3.0).abs() < 0.0001)
  assert_true((res[0][3] - 1.0).abs() < 0.0001)
  assert_true((res[0][4] - 5.0).abs() < 0.0001)
  // Column 1: mean = 30, median = 30, min = 10, max = 50
  assert_true((res[1][0] - 30.0).abs() < 0.0001)
  assert_true((res[1][1] - 30.0).abs() < 0.0001)
  assert_true((res[1][3] - 10.0).abs() < 0.0001)
  assert_true((res[1][4] - 50.0).abs() < 0.0001)
}

test "pc_column_summary_empty" {
  let res = @src.column_summary([])
  assert_eq(res.length(), 0)
}
