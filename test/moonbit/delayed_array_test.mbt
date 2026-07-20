///|
test "da_create_example_data" {
  let darr = @src.da_create_example_data()
  assert_eq(darr.nrow(), 100)
  assert_eq(darr.ncol(), 50)
}

///|
test "da_row_sums" {
  let darr = @src.da_create_example_data()
  let row_sums = @src.da_row_sums(darr)
  assert_eq(row_sums.length(), 100)
}

///|
test "da_col_sums" {
  let darr = @src.da_create_example_data()
  let col_sums = @src.da_col_sums(darr)
  assert_eq(col_sums.length(), 50)
}

///|
test "da_row_means" {
  let darr = @src.da_create_example_data()
  let row_means = @src.da_row_means(darr)
  assert_eq(row_means.length(), 100)
}

///|
test "da_col_means" {
  let darr = @src.da_create_example_data()
  let col_means = @src.da_col_means(darr)
  assert_eq(col_means.length(), 50)
}

///|
test "da_transpose" {
  let darr = @src.da_create_example_data()
  let transposed = @src.da_transpose(darr)
  assert_eq(transposed.nrow(), 50)
  assert_eq(transposed.ncol(), 100)
}

///|
test "da_subset_rows" {
  let darr = @src.da_create_example_data()
  let subset = @src.da_subset_rows(darr, [0, 1, 2])
  assert_eq(subset.nrow(), 3)
  assert_eq(subset.ncol(), 50)
}

///|
test "da_subset_cols" {
  let darr = @src.da_create_example_data()
  let subset = @src.da_subset_cols(darr, [0, 1, 2])
  assert_eq(subset.nrow(), 100)
  assert_eq(subset.ncol(), 3)
}

///|
test "da_compute" {
  let darr = @src.da_create_example_data()
  let data = @src.da_compute(darr)
  assert_eq(data.length(), 100)
}

///|
test "da_summary" {
  let darr = @src.da_create_example_data()
  let summary = @src.da_summary(darr)
  assert_true(summary.contains("DelayedArray"))
}
