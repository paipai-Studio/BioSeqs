///|
test "ch_create_example_data" {
  let data = @src.ch_create_example_data()
  assert_eq(data.matrix.length(), 5)
  assert_eq(data.row_names.length(), 5)
  assert_eq(data.col_names.length(), 5)
}

///|
test "ch_create_example_annotation" {
  let annot = @src.ch_create_example_annotation()
  assert_eq(annot.labels.length(), 5)
}

///|
test "ch_create_example_colormap" {
  let cmap = @src.ch_create_example_colormap()
  assert_eq(cmap.colors.length(), 11)
}

///|
test "ch_create_example_options" {
  let options = @src.ch_create_example_options()
  assert_eq(options.show_row_names, true)
}

///|
test "ch_draw_heatmap" {
  let data = @src.ch_create_example_data()
  let options = @src.ch_create_example_options()
  let result = @src.ch_draw_heatmap(data, options)
  assert_true(result.text.length() > 0)
}

///|
test "ch_cluster_rows" {
  let data = @src.ch_create_example_data()
  let row_order = @src.ch_cluster_rows(data.matrix, "complete")
  assert_eq(row_order.length(), 5)
}

///|
test "ch_cluster_cols" {
  let data = @src.ch_create_example_data()
  let col_order = @src.ch_cluster_cols(data.matrix, "complete")
  assert_eq(col_order.length(), 5)
}

///|
test "ch_distance_matrix" {
  let data = @src.ch_create_example_data()
  let dist_matrix = @src.ch_distance_matrix(data.matrix, "euclidean")
  assert_eq(dist_matrix.length(), 5)
}

///|
test "ch_euclidean_distance" {
  let v1 = [1.0, 2.0, 3.0]
  let v2 = [4.0, 5.0, 6.0]
  let dist = @src.ch_euclidean_distance(v1, v2)
  assert_true(dist > 0.0)
}

///|
test "ch_pearson_distance" {
  let v1 = [1.0, 2.0, 3.0]
  let v2 = [4.0, 5.0, 6.0]
  let dist = @src.ch_pearson_distance(v1, v2)
  assert_true(dist > -0.01)
}

///|
test "ch_cosine_distance" {
  let v1 = [1.0, 2.0, 3.0]
  let v2 = [4.0, 5.0, 6.0]
  let dist = @src.ch_cosine_distance(v1, v2)
  assert_true(dist >= 0.0)
}

///|
test "ch_scale_matrix" {
  let data = @src.ch_create_example_data()
  let scaled = @src.ch_scale_matrix(data.matrix, "row", true)
  assert_eq(scaled.length(), 5)
}

///|
test "ch_row_scale" {
  let data = @src.ch_create_example_data()
  let row_scaled = @src.ch_row_scale(data.matrix, true)
  assert_eq(row_scaled.length(), 5)
}

///|
test "ch_col_scale" {
  let data = @src.ch_create_example_data()
  let col_scaled = @src.ch_col_scale(data.matrix, true)
  assert_eq(col_scaled.length(), 5)
}

///|
test "ch_legend" {
  let cmap = @src.ch_create_example_colormap()
  let legend = @src.ch_legend(cmap, "Expression")
  assert_true(legend.length() > 0)
}

///|
test "ch_summary" {
  let data = @src.ch_create_example_data()
  let options = @src.ch_create_example_options()
  let result = @src.ch_draw_heatmap(data, options)
  let summary = @src.ch_summary(result)
  assert_true(summary.length() > 0)
}

///|
test "ch_row_split" {
  let data = @src.ch_create_example_data()
  let groups = ["Group A", "Group A", "Group B", "Group B", "Group C"]
  let row_split = @src.ch_row_split(data, groups)
  assert_eq(row_split.row_order.length(), 5)
}

///|
test "ch_col_split" {
  let data = @src.ch_create_example_data()
  let groups = ["Group A", "Group A", "Group B", "Group B", "Group C"]
  let col_split = @src.ch_col_split(data, groups)
  assert_eq(col_split.col_order.length(), 5)
}

///|
test "ch_merge_heatmaps" {
  let data = @src.ch_create_example_data()
  let options = @src.ch_create_example_options()
  let result = @src.ch_draw_heatmap(data, options)
  let merged = @src.ch_merge_heatmaps([result])
  assert_true(merged.length() > 0)
}

///|
test "ch_utils" {
  let arr = [1.0, 2.0, 3.0, 4.0, 5.0]
  let mean_val = @src.ch_mean(arr)
  assert_eq(mean_val, 3.0)

  let var_val = @src.ch_variance(arr)
  assert_true(var_val > 0.0)

  let std_val = @src.ch_std(arr)
  assert_true(std_val > 0.0)

  let min_val = @src.ch_min(arr)
  assert_eq(min_val, 1.0)

  let max_val = @src.ch_max(arr)
  assert_eq(max_val, 5.0)

  let sqrt_val = @src.ch_sqrt(16.0)
  assert_eq(sqrt_val, 4.0)
}
