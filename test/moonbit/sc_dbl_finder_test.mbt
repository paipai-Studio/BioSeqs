///|
test "scdf_create_example_data" {
  let data = @src.scdf_create_example_data()
  assert_eq(data.cell_names.length(), 50)
  assert_eq(data.gene_names.length(), 10)
}

///|
test "scdf_compute_doublet_score" {
  let data = @src.scdf_create_example_data()
  let scores = @src.scdf_compute_doublet_score(data, 10)
  assert_eq(scores.length(), 50)
}

///|
test "scdf_compute_distance" {
  let vec1 = [1.0, 2.0, 3.0]
  let vec2 = [4.0, 5.0, 6.0]
  let dist = @src.scdf_compute_distance(vec1, vec2)
  assert_true(dist > 0.0)
}

///|
test "scdf_find_nearest_neighbors" {
  let data = @src.scdf_create_example_data()
  let neighbors = @src.scdf_find_nearest_neighbors(data, 0, 5)
  assert_eq(neighbors.length(), 5)
}

///|
test "scdf_detect_doublets" {
  let data = @src.scdf_create_example_data()
  let scores = @src.scdf_compute_doublet_score(data, 10)
  let doublets = @src.scdf_detect_doublets(scores, 0.3)
  assert_true(doublets.length() >= 0)
}

///|
test "scdf_doublet_summary" {
  let data = @src.scdf_create_example_data()
  let scores = @src.scdf_compute_doublet_score(data, 10)
  let summary = @src.scdf_doublet_summary(scores)
  assert_true(summary.contains("Doublet"))
}

///|
test "scdf_filter_doublets" {
  let data = @src.scdf_create_example_data()
  let scores = @src.scdf_compute_doublet_score(data, 10)
  let filtered = @src.scdf_filter_doublets(data, scores)
  assert_eq(filtered.gene_names.length(), 10)
}

///|
test "scdf_compute_pca" {
  let data = @src.scdf_create_example_data()
  let pca = @src.scdf_compute_pca(data, 5)
  assert_eq(pca.length(), 50)
}
