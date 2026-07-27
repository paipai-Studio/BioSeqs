///|
/// Tests for DelayedMatrixStats module.

test "create_delayed_matrix_basic" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    ["row1", "row2"],
    ["col1", "col2", "col3"],
  )
  assert_eq(mat.nrow(), 2)
  assert_eq(mat.ncol(), 3)
  assert_eq(mat.row_names.length(), 2)
  assert_eq(mat.col_names.length(), 3)
}

///|
test "create_delayed_matrix_empty" {
  let mat = @src.create_delayed_matrix([], [], [])
  assert_eq(mat.nrow(), 0)
  assert_eq(mat.ncol(), 0)
}

///|
test "delayed_matrix_get_element" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0], [3.0, 4.0]],
    [],
    [],
  )
  assert_eq(mat.get(0, 0), 1.0)
  assert_eq(mat.get(0, 1), 2.0)
  assert_eq(mat.get(1, 0), 3.0)
  assert_eq(mat.get(1, 1), 4.0)
}

///|
test "delayed_matrix_dim" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    [],
    [],
  )
  let d = mat.dim()
  assert_eq(d[0], 2)
  assert_eq(d[1], 3)
}

///|
test "delayed_matrix_subset" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
    ["r1", "r2", "r3"],
    ["c1", "c2", "c3"],
  )
  let sub = @src.delayed_matrix_subset(mat, [0, 2], [1, 2])
  assert_eq(sub.nrow(), 2)
  assert_eq(sub.ncol(), 2)
  assert_eq(sub.get(0, 0), 2.0)
  assert_eq(sub.get(0, 1), 3.0)
  assert_eq(sub.get(1, 0), 8.0)
  assert_eq(sub.get(1, 1), 9.0)
}

///|
test "delayed_matrix_subset_names" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0], [3.0, 4.0]],
    ["r1", "r2"],
    ["c1", "c2"],
  )
  let sub = @src.delayed_matrix_subset(mat, [1], [0])
  assert_eq(sub.row_names.length(), 1)
  assert_eq(sub.col_names.length(), 1)
}

///|
test "row_stats_mean" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    [],
    [],
  )
  let means = @src.row_stats(mat, "mean")
  assert_eq(means.length(), 2)
  assert_eq(means[0], 2.0)
  assert_eq(means[1], 5.0)
}

///|
test "col_stats_mean" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    [],
    [],
  )
  let means = @src.col_stats(mat, "mean")
  assert_eq(means.length(), 3)
  assert_eq(means[0], 2.5)
  assert_eq(means[1], 3.5)
  assert_eq(means[2], 4.5)
}

///|
test "row_medians" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 3.0, 2.0], [4.0, 6.0, 5.0]],
    [],
    [],
  )
  let medians = @src.row_medians(mat)
  assert_eq(medians.length(), 2)
  assert_eq(medians[0], 2.0)
  assert_eq(medians[1], 5.0)
}

///|
test "col_medians" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 5.0], [3.0, 7.0], [2.0, 6.0]],
    [],
    [],
  )
  let medians = @src.col_medians(mat)
  assert_eq(medians.length(), 2)
  assert_eq(medians[0], 2.0)
  assert_eq(medians[1], 6.0)
}

///|
test "row_medians_even" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0, 4.0]],
    [],
    [],
  )
  let medians = @src.row_medians(mat)
  assert_eq(medians[0], 2.5)
}

///|
test "dms_row_means" {
  let mat = @src.create_delayed_matrix(
    [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]],
    [],
    [],
  )
  let means = @src.dms_row_means(mat)
  assert_eq(means.length(), 2)
  assert_eq(means[0], 20.0)
  assert_eq(means[1], 50.0)
}

///|
test "dms_col_means" {
  let mat = @src.create_delayed_matrix(
    [[10.0, 20.0], [30.0, 40.0], [50.0, 60.0]],
    [],
    [],
  )
  let means = @src.dms_col_means(mat)
  assert_eq(means.length(), 2)
  assert_eq(means[0], 30.0)
  assert_eq(means[1], 40.0)
}

///|
test "dms_row_vars" {
  let mat = @src.create_delayed_matrix(
    [[2.0, 4.0, 6.0]],
    [],
    [],
  )
  let vars = @src.dms_row_vars(mat)
  assert_eq(vars.length(), 1)
  assert_eq(vars[0], 4.0)
}

///|
test "col_vars" {
  let mat = @src.create_delayed_matrix(
    [[2.0], [4.0], [6.0]],
    [],
    [],
  )
  let vars = @src.col_vars(mat)
  assert_eq(vars.length(), 1)
  assert_eq(vars[0], 4.0)
}

///|
test "row_allsums" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    [],
    [],
  )
  let sums = @src.row_allsums(mat)
  assert_eq(sums.length(), 2)
  assert_eq(sums[0], 6.0)
  assert_eq(sums[1], 15.0)
}

///|
test "row_n_basic" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    [],
    [],
  )
  let n = @src.row_n(mat)
  assert_eq(n.length(), 2)
  assert_eq(n[0], 3)
  assert_eq(n[1], 3)
}

///|
test "col_n_basic" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
    [],
    [],
  )
  let n = @src.col_n(mat)
  assert_eq(n.length(), 2)
  assert_eq(n[0], 3)
  assert_eq(n[1], 3)
}

///|
test "row_stats_min_max" {
  let mat = @src.create_delayed_matrix(
    [[3.0, 1.0, 4.0, 1.0, 5.0], [9.0, 2.0, 6.0, 5.0, 3.0]],
    [],
    [],
  )
  let mins = @src.row_stats(mat, "min")
  let maxs = @src.row_stats(mat, "max")
  assert_eq(mins[0], 1.0)
  assert_eq(mins[1], 2.0)
  assert_eq(maxs[0], 5.0)
  assert_eq(maxs[1], 9.0)
}

///|
test "col_stats_min_max" {
  let mat = @src.create_delayed_matrix(
    [[3.0, 1.0, 4.0], [9.0, 2.0, 6.0]],
    [],
    [],
  )
  let mins = @src.col_stats(mat, "min")
  let maxs = @src.col_stats(mat, "max")
  assert_eq(mins[0], 3.0)
  assert_eq(mins[1], 1.0)
  assert_eq(mins[2], 4.0)
  assert_eq(maxs[0], 9.0)
  assert_eq(maxs[1], 2.0)
  assert_eq(maxs[2], 6.0)
}

///|
test "row_stats_sd" {
  let mat = @src.create_delayed_matrix(
    [[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]],
    [],
    [],
  )
  let sds = @src.row_stats(mat, "sd")
  assert_eq(sds.length(), 1)
  assert_eq(sds[0], (32.0 / 7.0).sqrt())
}

///|
test "row_stats_with_na" {
  let mat = @src.create_delayed_matrix(
    [[1.0, @double.not_a_number, 3.0], [4.0, 5.0, @double.not_a_number]],
    [],
    [],
  )
  let means = @src.dms_row_means(mat)
  assert_eq(means.length(), 2)
  assert_eq(means[0], 2.0)
  assert_eq(means[1], 4.5)
}

///|
test "row_n_with_na" {
  let mat = @src.create_delayed_matrix(
    [
      [1.0, @double.not_a_number, 3.0],
      [4.0, 5.0, @double.not_a_number],
      [@double.not_a_number, @double.not_a_number, @double.not_a_number],
    ],
    [],
    [],
  )
  let n = @src.row_n(mat)
  assert_eq(n.length(), 3)
  assert_eq(n[0], 2)
  assert_eq(n[1], 2)
  assert_eq(n[2], 0)
}

///|
test "col_n_with_na" {
  let mat = @src.create_delayed_matrix(
    [
      [1.0, 2.0, @double.not_a_number],
      [4.0, @double.not_a_number, 6.0],
      [@double.not_a_number, @double.not_a_number, 9.0],
    ],
    [],
    [],
  )
  let n = @src.col_n(mat)
  assert_eq(n.length(), 3)
  assert_eq(n[0], 2)
  assert_eq(n[1], 1)
  assert_eq(n[2], 2)
}

///|
test "row_medians_with_na" {
  let mat = @src.create_delayed_matrix(
    [
      [1.0, 3.0, @double.not_a_number, 5.0],
      [10.0, @double.not_a_number, 30.0, @double.not_a_number],
    ],
    [],
    [],
  )
  let medians = @src.row_medians(mat)
  assert_eq(medians.length(), 2)
  assert_eq(medians[0], 3.0)
  assert_eq(medians[1], 20.0)
}

///|
test "col_stats_all_na" {
  let mat = @src.create_delayed_matrix(
    [[@double.not_a_number, @double.not_a_number], [@double.not_a_number, @double.not_a_number]],
    [],
    [],
  )
  let means = @src.col_stats(mat, "mean")
  assert_eq(means.length(), 2)
  assert_eq(means[0], 0.0)
  assert_eq(means[1], 0.0)
}

///|
test "row_stats_n_stat" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0], [4.0, @double.not_a_number, 6.0]],
    [],
    [],
  )
  let n_vals = @src.row_stats(mat, "n")
  assert_eq(n_vals.length(), 2)
  assert_eq(n_vals[0], 3.0)
  assert_eq(n_vals[1], 2.0)
}

///|
test "row_stats_nna_stat" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0], [@double.not_a_number, 4.0]],
    [],
    [],
  )
  let nna_vals = @src.row_stats(mat, "nna")
  assert_eq(nna_vals.length(), 2)
  assert_eq(nna_vals[0], 2.0)
  assert_eq(nna_vals[1], 1.0)
}

///|
test "row_stats_nn_stat" {
  let mat = @src.create_delayed_matrix(
    [[1.0, @double.not_a_number, 3.0]],
    [],
    [],
  )
  let nn_vals = @src.row_stats(mat, "nn")
  assert_eq(nn_vals.length(), 1)
  assert_eq(nn_vals[0], 2.0)
}

///|
test "single_element_matrix" {
  let mat = @src.create_delayed_matrix([[42.0]], ["r1"], ["c1"])
  assert_eq(mat.nrow(), 1)
  assert_eq(mat.ncol(), 1)
  let means = @src.dms_row_means(mat)
  assert_eq(means[0], 42.0)
  let medians = @src.row_medians(mat)
  assert_eq(medians[0], 42.0)
  let sums = @src.row_allsums(mat)
  assert_eq(sums[0], 42.0)
}

///|
test "delayed_matrix_subset_single" {
  let mat = @src.create_delayed_matrix(
    [[10.0, 20.0], [30.0, 40.0]],
    ["a", "b"],
    ["x", "y"],
  )
  let sub = @src.delayed_matrix_subset(mat, [0], [1])
  assert_eq(sub.nrow(), 1)
  assert_eq(sub.ncol(), 1)
  assert_eq(sub.get(0, 0), 20.0)
  assert_eq(sub.row_names[0], "a")
  assert_eq(sub.col_names[0], "y")
}

///|
test "row_vars_single_element" {
  let mat = @src.create_delayed_matrix([[5.0]], [], [])
  let vars = @src.dms_row_vars(mat)
  assert_eq(vars.length(), 1)
  assert_eq(vars[0], 0.0)
}

///|
test "col_medians_odd" {
  let mat = @src.create_delayed_matrix(
    [[1.0], [5.0], [3.0]],
    [],
    [],
  )
  let medians = @src.col_medians(mat)
  assert_eq(medians[0], 3.0)
}

///|
test "col_means_with_na" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0], [@double.not_a_number, 4.0], [3.0, 6.0]],
    [],
    [],
  )
  let means = @src.dms_col_means(mat)
  assert_eq(means.length(), 2)
  assert_eq(means[0], 2.0)
  assert_eq(means[1], 4.0)
}

///|
test "row_stats_invalid_stat" {
  let mat = @src.create_delayed_matrix([[1.0, 2.0]], [], [])
  let result = @src.row_stats(mat, "invalid_stat")
  assert_eq(result.length(), 1)
  assert_eq(result[0], 0.0)
}

///|
test "matrix_stats_result_creation" {
  let result = @src.MatrixStatsResult::new(
    [1.0, 2.0, 3.0],
    "mean",
    "row",
    ["g1", "g2", "g3"],
    [],
  )
  assert_eq(result.values.length(), 3)
  assert_eq(result.stat_type, "mean")
  assert_eq(result.dimension, "row")
}

///|
test "row_stats_sum_equals_allsums" {
  let mat = @src.create_delayed_matrix(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    [],
    [],
  )
  let sums1 = @src.row_allsums(mat)
  let sums2 = @src.row_stats(mat, "sum")
  assert_eq(sums1.length(), sums2.length())
  assert_eq(sums1[0], sums2[0])
  assert_eq(sums1[1], sums2[1])
}