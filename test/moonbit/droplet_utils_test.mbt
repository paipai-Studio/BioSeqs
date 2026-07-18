///|
@test "droplet_utils_calculate_stats" {
  let (counts, barcodes) = @src.create_example_droplet_data()
  let stats = @src.calculate_droplet_stats(counts, barcodes)
  assert_eq(stats.total_counts.length(), 100)
  assert_eq(stats.n_genes.length(), 100)
}

@test "droplet_utils_barcode_ranking" {
  let total_counts = [100.0, 50.0, 200.0, 75.0]
  let ranked = @src.barcode_ranking(total_counts)
  assert_eq(ranked[0].0, 200.0)
  assert_eq(ranked[3].0, 50.0)
}

@test "droplet_utils_find_knee_point" {
  let ranked = [(1000.0, 0), (900.0, 1), (800.0, 2), (200.0, 3), (100.0, 4), (50.0, 5)]
  let knee = @src.find_knee_point(ranked)
  assert_true(knee >= 0 && knee < ranked.length())
}

@test "droplet_utils_empty_drops" {
  let (counts, barcodes) = @src.create_example_droplet_data()
  let result = @src.empty_drops(counts, barcodes, 10, 100)
  assert_true(result.barcode.length() > 0)
  assert_eq(result.p_value.length(), result.barcode.length())
}

@test "droplet_utils_filter_cells" {
  let (counts, barcodes) = @src.create_example_droplet_data()
  let result = @src.empty_drops(counts, barcodes, 10, 100)
  let (filtered_counts, filtered_barcodes) = @src.filter_cells_by_empty_drops(counts, barcodes, result)
  assert_eq(filtered_counts.length(), filtered_barcodes.length())
}

@test "droplet_utils_estimate_cell_number" {
  let total_counts = [1000.0, 900.0, 800.0, 200.0, 100.0, 50.0]
  let n_cells = @src.estimate_cell_number(total_counts)
  assert_true(n_cells >= 1 && n_cells <= total_counts.length())
}

@test "droplet_utils_create_example_data" {
  let (counts, barcodes) = @src.create_example_droplet_data()
  assert_eq(counts.length(), 100)
  assert_eq(barcodes.length(), 100)
  assert_eq(counts[0].length(), 500)
}
