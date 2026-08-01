///|
/// Tests for MatrixGenerics module.

// ---------------------------------------------------------------------------
// Data structure creation
// ---------------------------------------------------------------------------

test "mg_matrix_new_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0]])
  assert_eq(mat.dim_rows(), 2)
  assert_eq(mat.dim_cols(), 2)
  assert_eq(mat.get(0, 0), 1.0)
  assert_eq(mat.get(1, 1), 4.0)
}

test "mg_matrix_new_single_element" {
  let mat = @src.MgMatrix::new([[42.0]])
  assert_eq(mat.dim_rows(), 1)
  assert_eq(mat.dim_cols(), 1)
  assert_eq(mat.get(0, 0), 42.0)
}

test "mg_matrix_new_single_row" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0, 4.0]])
  assert_eq(mat.dim_rows(), 1)
  assert_eq(mat.dim_cols(), 4)
  assert_eq(mat.get(0, 0), 1.0)
  assert_eq(mat.get(0, 3), 4.0)
}

test "mg_matrix_new_single_col" {
  let mat = @src.MgMatrix::new([[1.0], [2.0], [3.0]])
  assert_eq(mat.dim_rows(), 3)
  assert_eq(mat.dim_cols(), 1)
  assert_eq(mat.get(0, 0), 1.0)
  assert_eq(mat.get(2, 0), 3.0)
}

test "mg_matrix_new_empty" {
  let mat = @src.MgMatrix::new([])
  assert_eq(mat.dim_rows(), 0)
  assert_eq(mat.dim_cols(), 0)
}

test "mg_matrix_zeros_basic" {
  let mat = @src.MgMatrix::zeros(3, 4)
  assert_eq(mat.dim_rows(), 3)
  assert_eq(mat.dim_cols(), 4)
  assert_eq(mat.get(0, 0), 0.0)
  assert_eq(mat.get(2, 3), 0.0)
}

test "mg_matrix_zeros_single" {
  let mat = @src.MgMatrix::zeros(1, 1)
  assert_eq(mat.dim_rows(), 1)
  assert_eq(mat.dim_cols(), 1)
  assert_eq(mat.get(0, 0), 0.0)
}

test "mg_matrix_zeros_empty" {
  let mat = @src.MgMatrix::zeros(0, 0)
  assert_eq(mat.dim_rows(), 0)
  assert_eq(mat.dim_cols(), 0)
}

// ---------------------------------------------------------------------------
// Accessor methods
// ---------------------------------------------------------------------------

test "mg_matrix_get_set" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0]])
  mat.set(0, 0, 99.0)
  assert_eq(mat.get(0, 0), 99.0)
  assert_eq(mat.get(0, 1), 2.0)
  assert_eq(mat.get(1, 0), 3.0)
}

test "mg_matrix_get_row" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
  let row0 = mat.get_row(0)
  assert_eq(row0.length(), 3)
  assert_eq(row0[0], 1.0)
  assert_eq(row0[1], 2.0)
  assert_eq(row0[2], 3.0)
  let row1 = mat.get_row(1)
  assert_eq(row1.length(), 3)
  assert_eq(row1[0], 4.0)
  assert_eq(row1[1], 5.0)
  assert_eq(row1[2], 6.0)
}

test "mg_matrix_get_col" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
  let col0 = mat.get_col(0)
  assert_eq(col0.length(), 3)
  assert_eq(col0[0], 1.0)
  assert_eq(col0[1], 3.0)
  assert_eq(col0[2], 5.0)
  let col1 = mat.get_col(1)
  assert_eq(col1.length(), 3)
  assert_eq(col1[0], 2.0)
  assert_eq(col1[1], 4.0)
  assert_eq(col1[2], 6.0)
}

// ---------------------------------------------------------------------------
// Row summary statistics
// ---------------------------------------------------------------------------

test "mg_row_means_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
  let result = @src.mg_row_means(mat)
  assert_eq(result.length(), 3)
  assert_eq(result[0], 2.0)
  assert_eq(result[1], 5.0)
  assert_eq(result[2], 8.0)
}

test "mg_row_means_constant" {
  let mat = @src.MgMatrix::new([[5.0, 5.0, 5.0], [3.0, 3.0, 3.0]])
  let result = @src.mg_row_means(mat)
  assert_eq(result[0], 5.0)
  assert_eq(result[1], 3.0)
}

test "mg_row_means_single_row" {
  let mat = @src.MgMatrix::new([[10.0, 20.0, 30.0]])
  let result = @src.mg_row_means(mat)
  assert_eq(result.length(), 1)
  assert_eq(result[0], 20.0)
}

test "mg_row_sums_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
  let result = @src.mg_row_sums(mat)
  assert_eq(result.length(), 3)
  assert_eq(result[0], 6.0)
  assert_eq(result[1], 15.0)
  assert_eq(result[2], 24.0)
}

test "mg_row_sums_constant" {
  let mat = @src.MgMatrix::new([[2.0, 2.0, 2.0], [0.0, 0.0, 0.0]])
  let result = @src.mg_row_sums(mat)
  assert_eq(result[0], 6.0)
  assert_eq(result[1], 0.0)
}

test "mg_row_vars_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
  let result = @src.mg_row_vars(mat)
  assert_eq(result.length(), 3)
  assert_true(result[0] > 0.99 && result[0] < 1.01)
  assert_true(result[1] > 0.99 && result[1] < 1.01)
  assert_true(result[2] > 0.99 && result[2] < 1.01)
}

test "mg_row_vars_constant" {
  let mat = @src.MgMatrix::new([[5.0, 5.0, 5.0], [2.0, 2.0, 2.0]])
  let result = @src.mg_row_vars(mat)
  assert_eq(result[0], 0.0)
  assert_eq(result[1], 0.0)
}

test "mg_row_vars_single" {
  let mat = @src.MgMatrix::new([[3.0]])
  let result = @src.mg_row_vars(mat)
  assert_eq(result[0], 0.0)
}

test "mg_row_sds_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
  let result = @src.mg_row_sds(mat)
  assert_eq(result.length(), 2)
  assert_true(result[0] > 0.99 && result[0] < 1.01)
  assert_true(result[1] > 0.99 && result[1] < 1.01)
}

test "mg_row_sds_constant" {
  let mat = @src.MgMatrix::new([[7.0, 7.0, 7.0]])
  let result = @src.mg_row_sds(mat)
  assert_eq(result[0], 0.0)
}

test "mg_row_medians_odd" {
  let mat = @src.MgMatrix::new([[1.0, 3.0, 2.0], [6.0, 4.0, 5.0]])
  let result = @src.mg_row_medians(mat)
  assert_eq(result[0], 2.0)
  assert_eq(result[1], 5.0)
}

test "mg_row_medians_even" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0, 4.0]])
  let result = @src.mg_row_medians(mat)
  assert_eq(result[0], 2.5)
}

test "mg_row_medians_constant" {
  let mat = @src.MgMatrix::new([[9.0, 9.0, 9.0]])
  let result = @src.mg_row_medians(mat)
  assert_eq(result[0], 9.0)
}

test "mg_row_mins_basic" {
  let mat = @src.MgMatrix::new([[3.0, 1.0, 2.0], [6.0, 5.0, 4.0], [9.0, 8.0, 7.0]])
  let result = @src.mg_row_mins(mat)
  assert_eq(result[0], 1.0)
  assert_eq(result[1], 4.0)
  assert_eq(result[2], 7.0)
}

test "mg_row_mins_constant" {
  let mat = @src.MgMatrix::new([[4.0, 4.0, 4.0]])
  let result = @src.mg_row_mins(mat)
  assert_eq(result[0], 4.0)
}

test "mg_row_maxs_basic" {
  let mat = @src.MgMatrix::new([[3.0, 1.0, 2.0], [6.0, 5.0, 4.0], [9.0, 8.0, 7.0]])
  let result = @src.mg_row_maxs(mat)
  assert_eq(result[0], 3.0)
  assert_eq(result[1], 6.0)
  assert_eq(result[2], 9.0)
}

test "mg_row_maxs_constant" {
  let mat = @src.MgMatrix::new([[4.0, 4.0, 4.0]])
  let result = @src.mg_row_maxs(mat)
  assert_eq(result[0], 4.0)
}

test "mg_row_ranges_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [10.0, 5.0, 2.0]])
  let result = @src.mg_row_ranges(mat)
  assert_eq(result[0], 2.0)
  assert_eq(result[1], 8.0)
}

test "mg_row_ranges_constant" {
  let mat = @src.MgMatrix::new([[5.0, 5.0, 5.0]])
  let result = @src.mg_row_ranges(mat)
  assert_eq(result[0], 0.0)
}

test "mg_row_counts_basic" {
  let mat = @src.MgMatrix::new([[1.0, 0.0, 3.0], [0.0, 0.0, 0.0], [4.0, 5.0, 6.0]])
  let result = @src.mg_row_counts(mat)
  assert_eq(result.length(), 3)
  assert_eq(result[0], 2)
  assert_eq(result[1], 0)
  assert_eq(result[2], 3)
}

test "mg_row_counts_all_zero" {
  let mat = @src.MgMatrix::new([[0.0, 0.0], [0.0, 0.0]])
  let result = @src.mg_row_counts(mat)
  assert_eq(result[0], 0)
  assert_eq(result[1], 0)
}

test "mg_row_mads_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
  let result = @src.mg_row_mads(mat)
  assert_eq(result.length(), 2)
  assert_true(result[0] > 1.4 && result[0] < 1.6)
  assert_true(result[1] > 1.4 && result[1] < 1.6)
}

test "mg_row_mads_constant" {
  let mat = @src.MgMatrix::new([[5.0, 5.0, 5.0]])
  let result = @src.mg_row_mads(mat)
  assert_eq(result[0], 0.0)
}

// ---------------------------------------------------------------------------
// Row anys and alls
// ---------------------------------------------------------------------------

test "mg_row_anys_basic" {
  let mat = @src.MgMatrix::new([[1.0, 5.0, 3.0], [2.0, 2.0, 2.0], [10.0, 20.0, 30.0]])
  let result = @src.mg_row_anys(mat, 4.0)
  assert_eq(result.length(), 3)
  assert_true(result[0])
  assert_false(result[1])
  assert_true(result[2])
}

test "mg_row_anys_all_below" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0]])
  let result = @src.mg_row_anys(mat, 100.0)
  assert_false(result[0])
}

test "mg_row_anys_threshold_boundary" {
  let mat = @src.MgMatrix::new([[5.0, 5.0, 5.0]])
  let result = @src.mg_row_anys(mat, 5.0)
  assert_false(result[0])
}

test "mg_row_alls_basic" {
  let mat = @src.MgMatrix::new([[10.0, 20.0, 30.0], [1.0, 100.0, 50.0], [2.0, 2.0, 2.0]])
  let result = @src.mg_row_alls(mat, 5.0)
  assert_eq(result.length(), 3)
  assert_true(result[0])
  assert_false(result[1])
  assert_false(result[2])
}

test "mg_row_alls_all_above" {
  let mat = @src.MgMatrix::new([[10.0, 10.0], [20.0, 20.0]])
  let result = @src.mg_row_alls(mat, 5.0)
  assert_true(result[0])
  assert_true(result[1])
}

test "mg_row_alls_threshold_boundary" {
  let mat = @src.MgMatrix::new([[5.0, 5.0, 5.0]])
  let result = @src.mg_row_alls(mat, 5.0)
  assert_false(result[0])
}

// ---------------------------------------------------------------------------
// Column summary statistics
// ---------------------------------------------------------------------------

test "mg_col_means_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
  let result = @src.mg_col_means(mat)
  assert_eq(result.length(), 3)
  assert_eq(result[0], 4.0)
  assert_eq(result[1], 5.0)
  assert_eq(result[2], 6.0)
}

test "mg_col_means_constant" {
  let mat = @src.MgMatrix::new([[3.0, 7.0], [3.0, 7.0], [3.0, 7.0]])
  let result = @src.mg_col_means(mat)
  assert_eq(result[0], 3.0)
  assert_eq(result[1], 7.0)
}

test "mg_col_means_single_col" {
  let mat = @src.MgMatrix::new([[10.0], [20.0], [30.0]])
  let result = @src.mg_col_means(mat)
  assert_eq(result.length(), 1)
  assert_eq(result[0], 20.0)
}

test "mg_col_sums_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
  let result = @src.mg_col_sums(mat)
  assert_eq(result.length(), 3)
  assert_eq(result[0], 12.0)
  assert_eq(result[1], 15.0)
  assert_eq(result[2], 18.0)
}

test "mg_col_sums_constant" {
  let mat = @src.MgMatrix::new([[2.0, 0.0], [2.0, 0.0], [2.0, 0.0]])
  let result = @src.mg_col_sums(mat)
  assert_eq(result[0], 6.0)
  assert_eq(result[1], 0.0)
}

test "mg_col_vars_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
  let result = @src.mg_col_vars(mat)
  assert_eq(result.length(), 3)
  assert_true(result[0] > 8.99 && result[0] < 9.01)
  assert_true(result[1] > 8.99 && result[1] < 9.01)
  assert_true(result[2] > 8.99 && result[2] < 9.01)
}

test "mg_col_vars_constant" {
  let mat = @src.MgMatrix::new([[5.0, 7.0], [5.0, 7.0], [5.0, 7.0]])
  let result = @src.mg_col_vars(mat)
  assert_eq(result[0], 0.0)
  assert_eq(result[1], 0.0)
}

test "mg_col_vars_single" {
  let mat = @src.MgMatrix::new([[3.0], [3.0]])
  let result = @src.mg_col_vars(mat)
  assert_eq(result[0], 0.0)
}

test "mg_col_sds_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [4.0, 5.0], [7.0, 8.0]])
  let result = @src.mg_col_sds(mat)
  assert_eq(result.length(), 2)
  assert_true(result[0] > 2.99 && result[0] < 3.01)
  assert_true(result[1] > 2.99 && result[1] < 3.01)
}

test "mg_col_sds_constant" {
  let mat = @src.MgMatrix::new([[3.0, 5.0], [3.0, 5.0]])
  let result = @src.mg_col_sds(mat)
  assert_eq(result[0], 0.0)
  assert_eq(result[1], 0.0)
}

test "mg_col_medians_odd" {
  let mat = @src.MgMatrix::new([[1.0, 3.0], [4.0, 1.0], [7.0, 2.0]])
  let result = @src.mg_col_medians(mat)
  assert_eq(result[0], 4.0)
  assert_eq(result[1], 2.0)
}

test "mg_col_medians_even" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0]])
  let result = @src.mg_col_medians(mat)
  assert_eq(result[0], 2.0)
  assert_eq(result[1], 3.0)
}

test "mg_col_medians_constant" {
  let mat = @src.MgMatrix::new([[8.0, 2.0], [8.0, 2.0], [8.0, 2.0]])
  let result = @src.mg_col_medians(mat)
  assert_eq(result[0], 8.0)
  assert_eq(result[1], 2.0)
}

test "mg_col_mins_basic" {
  let mat = @src.MgMatrix::new([[3.0, 1.0, 2.0], [6.0, 5.0, 4.0], [9.0, 8.0, 7.0]])
  let result = @src.mg_col_mins(mat)
  assert_eq(result[0], 3.0)
  assert_eq(result[1], 1.0)
  assert_eq(result[2], 2.0)
}

test "mg_col_mins_constant" {
  let mat = @src.MgMatrix::new([[4.0, 6.0], [4.0, 6.0], [4.0, 6.0]])
  let result = @src.mg_col_mins(mat)
  assert_eq(result[0], 4.0)
  assert_eq(result[1], 6.0)
}

test "mg_col_maxs_basic" {
  let mat = @src.MgMatrix::new([[3.0, 1.0, 2.0], [6.0, 5.0, 4.0], [9.0, 8.0, 7.0]])
  let result = @src.mg_col_maxs(mat)
  assert_eq(result[0], 9.0)
  assert_eq(result[1], 8.0)
  assert_eq(result[2], 7.0)
}

test "mg_col_maxs_constant" {
  let mat = @src.MgMatrix::new([[4.0, 6.0], [4.0, 6.0], [4.0, 6.0]])
  let result = @src.mg_col_maxs(mat)
  assert_eq(result[0], 4.0)
  assert_eq(result[1], 6.0)
}

test "mg_col_ranges_basic" {
  let mat = @src.MgMatrix::new([[1.0, 10.0], [5.0, 2.0], [10.0, 1.0]])
  let result = @src.mg_col_ranges(mat)
  assert_eq(result[0], 9.0)
  assert_eq(result[1], 9.0)
}

test "mg_col_ranges_constant" {
  let mat = @src.MgMatrix::new([[5.0, 3.0], [5.0, 3.0], [5.0, 3.0]])
  let result = @src.mg_col_ranges(mat)
  assert_eq(result[0], 0.0)
  assert_eq(result[1], 0.0)
}

// ---------------------------------------------------------------------------
// Col anys and alls
// ---------------------------------------------------------------------------

test "mg_col_anys_basic" {
  let mat = @src.MgMatrix::new([[1.0, 5.0, 3.0], [2.0, 2.0, 2.0], [1.0, 8.0, 1.0]])
  let result = @src.mg_col_anys(mat, 4.0)
  assert_eq(result.length(), 3)
  assert_false(result[0])
  assert_true(result[1])
  assert_false(result[2])
}

test "mg_col_anys_all_below" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0]])
  let result = @src.mg_col_anys(mat, 100.0)
  assert_false(result[0])
  assert_false(result[1])
}

test "mg_col_alls_basic" {
  let mat = @src.MgMatrix::new([[10.0, 5.0], [20.0, 6.0], [30.0, 7.0]])
  let result = @src.mg_col_alls(mat, 4.0)
  assert_true(result[0])
  assert_true(result[1])
}

test "mg_col_alls_some_below" {
  let mat = @src.MgMatrix::new([[10.0, 3.0], [20.0, 6.0], [30.0, 7.0]])
  let result = @src.mg_col_alls(mat, 5.0)
  assert_true(result[0])
  assert_false(result[1])
}

// ---------------------------------------------------------------------------
// Block processing
// ---------------------------------------------------------------------------

test "mg_block_apply_rows_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
  let result = @src.mg_block_apply_rows(mat, 2, (block) => {
    @src.mg_row_means(block)
  })
  assert_eq(result.length(), 3)
  assert_eq(result[0], 1.5)
  assert_eq(result[1], 3.5)
  assert_eq(result[2], 5.5)
}

test "mg_block_apply_rows_single_block" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0]])
  let result = @src.mg_block_apply_rows(mat, 10, (block) => {
    @src.mg_row_sums(block)
  })
  assert_eq(result.length(), 2)
  assert_eq(result[0], 3.0)
  assert_eq(result[1], 7.0)
}

test "mg_block_apply_rows_block_size_exceeds" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0]])
  let result = @src.mg_block_apply_rows(mat, 100, (block) => {
    @src.mg_row_mins(block)
  })
  assert_eq(result.length(), 2)
  assert_eq(result[0], 1.0)
  assert_eq(result[1], 3.0)
}

test "mg_block_apply_cols_basic" {
  let mat = @src.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
  let result = @src.mg_block_apply_cols(mat, 2, (block) => {
    @src.mg_col_means(block)
  })
  assert_eq(result.length(), 3)
  assert_eq(result[0], 2.5)
  assert_eq(result[1], 3.5)
  assert_eq(result[2], 4.5)
}

test "mg_block_apply_cols_single_block" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0]])
  let result = @src.mg_block_apply_cols(mat, 10, (block) => {
    @src.mg_col_sums(block)
  })
  assert_eq(result.length(), 2)
  assert_eq(result[0], 4.0)
  assert_eq(result[1], 6.0)
}

// ---------------------------------------------------------------------------
// Standalone statistical functions
// ---------------------------------------------------------------------------

test "mg_mean_basic" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  assert_eq(@src.mg_mean(data), 3.0)
}

test "mg_mean_single" {
  let data = [42.0]
  assert_eq(@src.mg_mean(data), 42.0)
}

test "mg_mean_empty" {
  let data : Array[Double] = []
  assert_eq(@src.mg_mean(data), 0.0)
}

test "mg_mean_constant" {
  let data = [7.0, 7.0, 7.0, 7.0]
  assert_eq(@src.mg_mean(data), 7.0)
}

test "mg_mean_negative" {
  let data = [-1.0, -2.0, -3.0]
  assert_eq(@src.mg_mean(data), -2.0)
}

test "mg_variance_basic" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.mg_variance(data)
  assert_true(result > 2.49 && result < 2.51)
}

test "mg_variance_single" {
  let data = [5.0]
  assert_eq(@src.mg_variance(data), 0.0)
}

test "mg_variance_empty" {
  let data : Array[Double] = []
  assert_eq(@src.mg_variance(data), 0.0)
}

test "mg_variance_constant" {
  let data = [3.0, 3.0, 3.0]
  assert_eq(@src.mg_variance(data), 0.0)
}

test "mg_variance_two_elements" {
  let data = [2.0, 4.0]
  let result = @src.mg_variance(data)
  assert_true(result > 1.99 && result < 2.01)
}

test "mg_median_odd" {
  let data = [3.0, 1.0, 2.0]
  assert_eq(@src.mg_median(data), 2.0)
}

test "mg_median_even" {
  let data = [1.0, 2.0, 3.0, 4.0]
  assert_eq(@src.mg_median(data), 2.5)
}

test "mg_median_single" {
  let data = [99.0]
  assert_eq(@src.mg_median(data), 99.0)
}

test "mg_median_empty" {
  let data : Array[Double] = []
  assert_eq(@src.mg_median(data), 0.0)
}

test "mg_median_two" {
  let data = [5.0, 10.0]
  assert_eq(@src.mg_median(data), 7.5)
}

test "mg_min_basic" {
  let data = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0]
  assert_eq(@src.mg_min(data), 1.0)
}

test "mg_min_single" {
  let data = [42.0]
  assert_eq(@src.mg_min(data), 42.0)
}

test "mg_min_empty" {
  let data : Array[Double] = []
  assert_eq(@src.mg_min(data), 0.0)
}

test "mg_min_negative" {
  let data = [-5.0, -1.0, -3.0]
  assert_eq(@src.mg_min(data), -5.0)
}

test "mg_max_basic" {
  let data = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0]
  assert_eq(@src.mg_max(data), 9.0)
}

test "mg_max_single" {
  let data = [42.0]
  assert_eq(@src.mg_max(data), 42.0)
}

test "mg_max_empty" {
  let data : Array[Double] = []
  assert_eq(@src.mg_max(data), 0.0)
}

test "mg_max_negative" {
  let data = [-5.0, -1.0, -3.0]
  assert_eq(@src.mg_max(data), -1.0)
}

test "mg_mad_basic" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.mg_mad(data)
  assert_true(result > 1.4 && result < 2.0)
}

test "mg_mad_constant" {
  let data = [5.0, 5.0, 5.0]
  assert_eq(@src.mg_mad(data), 0.0)
}

test "mg_mad_single" {
  let data = [7.0]
  assert_eq(@src.mg_mad(data), 0.0)
}

test "mg_mad_empty" {
  let data : Array[Double] = []
  assert_eq(@src.mg_mad(data), 0.0)
}

test "mg_mad_two_same" {
  let data = [3.0, 3.0]
  assert_eq(@src.mg_mad(data), 0.0)
}

// ---------------------------------------------------------------------------
// Edge cases
// ---------------------------------------------------------------------------

test "mg_empty_matrix_row_stats" {
  let mat = @src.MgMatrix::zeros(0, 5)
  let means = @src.mg_row_means(mat)
  assert_eq(means.length(), 0)
  let sums = @src.mg_row_sums(mat)
  assert_eq(sums.length(), 0)
}

test "mg_empty_matrix_col_stats" {
  let mat = @src.MgMatrix::zeros(5, 0)
  let means = @src.mg_col_means(mat)
  assert_eq(means.length(), 0)
}

test "mg_single_value_matrix" {
  let mat = @src.MgMatrix::new([[7.0]])
  assert_eq(mat.dim_rows(), 1)
  assert_eq(mat.dim_cols(), 1)
  assert_eq(@src.mg_row_means(mat)[0], 7.0)
  assert_eq(@src.mg_col_means(mat)[0], 7.0)
  assert_eq(@src.mg_row_vars(mat)[0], 0.0)
  assert_eq(@src.mg_col_vars(mat)[0], 0.0)
  assert_eq(@src.mg_row_medians(mat)[0], 7.0)
  assert_eq(@src.mg_col_medians(mat)[0], 7.0)
}

test "mg_negative_values_matrix" {
  let mat = @src.MgMatrix::new([[-1.0, -2.0], [-3.0, -4.0]])
  let row_means = @src.mg_row_means(mat)
  assert_eq(row_means[0], -1.5)
  assert_eq(row_means[1], -3.5)
  let col_means = @src.mg_col_means(mat)
  assert_eq(col_means[0], -2.0)
  assert_eq(col_means[1], -3.0)
}

test "mg_large_matrix_stats" {
  let mat = @src.MgMatrix::zeros(100, 50)
  let row_means = @src.mg_row_means(mat)
  assert_eq(row_means.length(), 100)
  let col_means = @src.mg_col_means(mat)
  assert_eq(col_means.length(), 50)
  assert_eq(row_means[0], 0.0)
  assert_eq(col_means[0], 0.0)
}

test "mg_row_anys_empty_matrix" {
  let mat = @src.MgMatrix::zeros(0, 3)
  let result = @src.mg_row_anys(mat, 1.0)
  assert_eq(result.length(), 0)
}

test "mg_row_alls_empty_matrix" {
  let mat = @src.MgMatrix::zeros(0, 3)
  let result = @src.mg_row_alls(mat, 1.0)
  assert_eq(result.length(), 0)
}

test "mg_col_anys_empty_matrix" {
  let mat = @src.MgMatrix::zeros(3, 0)
  let result = @src.mg_col_anys(mat, 1.0)
  assert_eq(result.length(), 0)
}

test "mg_col_alls_empty_matrix" {
  let mat = @src.MgMatrix::zeros(3, 0)
  let result = @src.mg_col_alls(mat, 1.0)
  assert_eq(result.length(), 0)
}

test "mg_block_apply_rows_empty" {
  let mat = @src.MgMatrix::zeros(0, 5)
  let result = @src.mg_block_apply_rows(mat, 10, (block) => {
    @src.mg_row_means(block)
  })
  assert_eq(result.length(), 0)
}

test "mg_block_apply_cols_empty" {
  let mat = @src.MgMatrix::zeros(5, 0)
  let result = @src.mg_block_apply_cols(mat, 10, (block) => {
    @src.mg_col_means(block)
  })
  assert_eq(result.length(), 0)
}

// ---------------------------------------------------------------------------
// Cross-validation: row and col stats should agree on the same data
// ---------------------------------------------------------------------------

test "mg_row_col_stats_consistency" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [3.0, 4.0]])
  let row_means = @src.mg_row_means(mat)
  let col_means = @src.mg_col_means(mat)
  assert_eq(row_means.length(), 2)
  assert_eq(col_means.length(), 2)
  assert_eq(@src.mg_mean(row_means), @src.mg_mean(col_means))
}

test "mg_symmetric_matrix_row_col_stats" {
  let mat = @src.MgMatrix::new([[1.0, 2.0], [2.0, 1.0]])
  let row_means = @src.mg_row_means(mat)
  let col_means = @src.mg_col_means(mat)
  assert_eq(row_means[0], col_means[0])
  assert_eq(row_means[1], col_means[1])
}