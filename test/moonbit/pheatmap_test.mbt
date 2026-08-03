///|
/// Tests for pheatmap - enhanced heatmap with clustering

///|
test "pheatmap_create_input" {
  let mat = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
  ]
  let input = @src.pheatmap_input(
    mat,
    ["row1", "row2", "row3"],
    ["col1", "col2", "col3"],
  )
  assert_eq(input.mat.length(), 3)
  assert_eq(input.row_names.length(), 3)
  assert_eq(input.col_names.length(), 3)
}

///|
test "pheatmap_distance_euclidean" {
  let data = [
    [1.0, 0.0],
    [0.0, 0.0],
  ]
  let dist = @src.pheatmap_distance_matrix(data, @src.distance_method_euclidean())
  assert_eq(dist.length(), 2)
  assert_eq(dist[0][0], 0.0)
  assert_eq(dist[1][1], 0.0)
  assert_true(dist[0][1] > 0.0)
}

///|
test "pheatmap_distance_manhattan" {
  let data = [
    [1.0, 2.0],
    [4.0, 6.0],
  ]
  let dist = @src.pheatmap_distance_matrix(data, @src.distance_method_manhattan())
  assert_eq(dist.length(), 2)
  assert_eq(dist[0][1], 7.0)
}

///|
test "pheatmap_distance_symmetry" {
  let data = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
  ]
  let dist = @src.pheatmap_distance_matrix(data, @src.distance_method_euclidean())
  let n = dist.length()
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      assert_eq(dist[i][j], dist[j][i])
      j = j + 1
    }
    i = i + 1
  }
}

///|
test "pheatmap_column_distance" {
  let data = [
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 0.0],
  ]
  let dist = @src.pheatmap_column_distance_matrix(data, @src.distance_method_euclidean())
  assert_eq(dist.length(), 3)
}

///|
test "pheatmap_render_basic" {
  let mat = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
  ]
  let input = @src.pheatmap_input(
    mat,
    ["r1", "r2", "r3"],
    ["c1", "c2", "c3"],
  )
  let result = @src.pheatmap_render(input)
  assert_eq(result.row_order.length(), 3)
  assert_eq(result.col_order.length(), 3)
  assert_true(result.text.length() > 0)
}

///|
test "pheatmap_set_title" {
  let mat = [[1.0, 2.0], [3.0, 4.0]]
  let input = @src.pheatmap_input(mat, ["r1", "r2"], ["c1", "c2"])
  let titled = @src.pheatmap_set_title(input, "Test Heatmap")
  assert_eq(titled.main_title, "Test Heatmap")
}

///|
test "pheatmap_set_fontsize" {
  let mat = [[1.0, 2.0], [3.0, 4.0]]
  let input = @src.pheatmap_input(mat, ["r1", "r2"], ["c1", "c2"])
  let resized = @src.pheatmap_set_fontsize(input, 18.0)
  assert_eq(resized.font_size, 18.0)
}

///|
test "pheatmap_cutree_rows" {
  let mat = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
  ]
  let input = @src.pheatmap_input(mat, ["r1", "r2", "r3"], ["c1", "c2", "c3"])
  let cut = @src.pheatmap_set_cutree(input, 2, 1)
  assert_eq(cut.cutree_rows, 2)
}

///|
test "pheatmap_no_clustering" {
  let mat = [
    [1.0, 2.0],
    [3.0, 4.0],
  ]
  let input = @src.pheatmap_input(mat, ["r1", "r2"], ["c1", "c2"])
  let unclustered = @src.pheatmap_set_cluster(input, false, false)
  assert_eq(unclustered.cluster_rows, false)
  assert_eq(unclustered.cluster_cols, false)
}

///|
test "pheatmap_empty_matrix" {
  let mat : Array[Array[Double]] = []
  let input = @src.pheatmap_input(mat, [], [])
  let result = @src.pheatmap_render(input)
  assert_eq(result.row_order.length(), 0)
  assert_eq(result.col_order.length(), 0)
}

///|
test "pheatmap_color_scheme" {
  let mat = [[1.0, 2.0], [3.0, 4.0]]
  let input = @src.pheatmap_input(mat, ["r1", "r2"], ["c1", "c2"])
  assert_true(input.color.length() > 0)
}

///|
test "pheatmap_display_numbers" {
  let mat = [[1.5, 2.5], [3.5, 4.5]]
  let input = @src.pheatmap_input(mat, ["r1", "r2"], ["c1", "c2"])
  let with_numbers = @src.pheatmap_set_display_numbers(input, true)
  assert_eq(with_numbers.display_numbers, true)
}
