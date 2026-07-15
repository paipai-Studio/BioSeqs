///|
test "summarized_experiment_create" {
  let assays = Map([("counts", [[1.0, 2.0], [3.0, 4.0]])])
  let row_ranges = [("chr1", 1, 100), ("chr2", 200, 300)]
  let col_data = [Map([("sample", "s1")]), Map([("sample", "s2")])]
  let metadata = Map([("experiment", "test")])

  let se = @src.summarized_experiment(assays, row_ranges, col_data, metadata)

  assert_eq(@src.se_assay_names(se).length(), 1)
  assert_eq(@src.se_nrow(se), 2)
  assert_eq(@src.se_ncol(se), 2)
}

///|
test "summarized_experiment_assay" {
  let assays = Map([("counts", [[1.0, 2.0], [3.0, 4.0]])])
  let row_ranges = [("chr1", 1, 100), ("chr2", 200, 300)]
  let col_data = [Map([("sample", "s1")]), Map([("sample", "s2")])]
  let metadata = Map([])

  let se = @src.summarized_experiment(assays, row_ranges, col_data, metadata)

  match @src.se_assay(se, "counts") {
    Some(a) => assert_eq(a.length(), 2)
    None => assert_true(false)
  }
}

///|
test "summarized_experiment_row_col_means" {
  let assays = Map([("counts", [[1.0, 3.0], [2.0, 4.0]])])
  let row_ranges = [("chr1", 1, 100), ("chr2", 200, 300)]
  let col_data = [Map([("sample", "s1")]), Map([("sample", "s2")])]
  let metadata = Map([])

  let se = @src.summarized_experiment(assays, row_ranges, col_data, metadata)
  let row_means = @src.se_row_means(se, "counts")
  let col_means = @src.se_col_means(se, "counts")

  assert_eq(row_means.length(), 2)
  assert_eq(col_means.length(), 2)
}

///|
test "summarized_experiment_subset_rows" {
  let assays = Map([("counts", [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])])
  let row_ranges = [("chr1", 1, 100), ("chr2", 200, 300), ("chr3", 400, 500)]
  let col_data = [Map([("sample", "s1")]), Map([("sample", "s2")])]
  let metadata = Map([])

  let se = @src.summarized_experiment(assays, row_ranges, col_data, metadata)
  let se_sub = @src.se_subset_rows(se, [0, 2])

  assert_eq(@src.se_nrow(se_sub), 2)
}

///|
test "summarized_experiment_subset_cols" {
  let assays = Map([("counts", [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])])
  let row_ranges = [("chr1", 1, 100), ("chr2", 200, 300)]
  let col_data = [
    Map([("sample", "s1")]),
    Map([("sample", "s2")]),
    Map([("sample", "s3")]),
  ]
  let metadata = Map([])

  let se = @src.summarized_experiment(assays, row_ranges, col_data, metadata)
  let se_sub = @src.se_subset_cols(se, [0, 2])

  assert_eq(@src.se_ncol(se_sub), 2)
}

///|
test "summarized_experiment_merge" {
  let assays1 = Map([("counts", [[1.0, 2.0], [3.0, 4.0]])])
  let row_ranges1 = [("chr1", 1, 100), ("chr2", 200, 300)]
  let col_data1 = [Map([("sample", "s1")]), Map([("sample", "s2")])]
  let metadata1 = Map([("exp", "test")])

  let assays2 = Map([("counts", [[5.0, 6.0]])])
  let row_ranges2 = [("chr3", 400, 500)]
  let col_data2 = [Map([("sample", "s3")])]
  let metadata2 = Map([])

  let se1 = @src.summarized_experiment(
    assays1, row_ranges1, col_data1, metadata1,
  )
  let se2 = @src.summarized_experiment(
    assays2, row_ranges2, col_data2, metadata2,
  )
  let merged = @src.se_merge(se1, se2)

  assert_eq(@src.se_nrow(merged), 3)
}

///|
test "summarized_experiment_empty" {
  let assays = Map([])
  let row_ranges = []
  let col_data = []
  let metadata = Map([])

  let se = @src.summarized_experiment(assays, row_ranges, col_data, metadata)

  assert_eq(@src.se_nrow(se), 0)
  assert_eq(@src.se_ncol(se), 0)
}
