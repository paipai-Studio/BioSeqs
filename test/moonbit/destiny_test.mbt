///|
test "destiny_create_cell_data" {
  let cell = @src.CellData::new(
    "cell_001",
    [1.0, 2.0, 3.0],
    "cluster_1"
  )
  
  assert_eq(cell.cell_id, "cell_001")
  assert_eq(cell.expression.length(), 3)
  assert_eq(cell.cluster, "cluster_1")
}

///|
test "destiny_distance_matrix" {
  let cells = @src.create_example_sc_data(5, 3)
  let dist_matrix = @src.compute_distance_matrix(cells, "euclidean")
  
  assert_eq(dist_matrix.n_cells, 5)
  assert_eq(dist_matrix.cells.length(), 5)
  
  // Check symmetry
  assert_true((dist_matrix.distances[0][1] - dist_matrix.distances[1][0]).abs() < 1.0e-10)
  
  // Check diagonal is zero
  assert_true(dist_matrix.distances[0][0].abs() < 1.0e-10)
}

///|
test "destiny_manhattan_distance" {
  let cells = @src.create_example_sc_data(3, 2)
  let dist_matrix = @src.compute_distance_matrix(cells, "manhattan")
  
  assert_eq(dist_matrix.n_cells, 3)
  assert_true(dist_matrix.distances[0][1] >= 0.0)
}

///|
test "destiny_gaussian_kernel" {
  let cells = @src.create_example_sc_data(5, 3)
  let dist_matrix = @src.compute_distance_matrix(cells, "euclidean")
  let kernel = @src.compute_gaussian_kernel(dist_matrix, 1.0)
  
  assert_eq(kernel.cells.length(), 5)
  assert_eq(kernel.bandwidth, 1.0)
  
  // Check kernel values are between 0 and 1
  assert_true(kernel.kernel[0][0] <= 1.0)
  assert_true(kernel.kernel[0][0] >= 0.0)
  
  // Diagonal should be 1 (exp(0))
  assert_true((kernel.kernel[0][0] - 1.0).abs() < 1.0e-10)
}

///|
test "destiny_find_sigma" {
  let cells = @src.create_example_sc_data(10, 3)
  let dist_matrix = @src.compute_distance_matrix(cells, "euclidean")
  let sigma = @src.find_sigma_automatic(dist_matrix)
  
  assert_true(sigma > 0.0)
}

///|
test "destiny_diffusion_map" {
  let cells = @src.create_example_sc_data(20, 5)
  let result = @src.compute_diffusion_map(cells, 3, 1.0)
  
  assert_eq(result.cell_ids.length(), 20)
  assert_eq(result.eigenvalues.length(), 3)
  assert_eq(result.embedding.length(), 20)
  assert_eq(result.explained_variance.length(), 3)
  
  // Check eigenvalues are non-negative
  let mut i = 0
  while i < 3 {
    assert_true(result.eigenvalues[i] >= 0.0)
    i = i + 1
  }
}

///|
test "destiny_auto_embedding" {
  let cells = @src.create_example_sc_data(15, 4)
  let result = @src.destiny_create_embedding(cells, 2)
  
  assert_eq(result.cell_ids.length(), 15)
  assert_eq(result.embedding[0].length(), 2)
}

///|
test "destiny_plot_coordinates" {
  let cells = @src.create_example_sc_data(10, 3)
  let result = @src.destiny_create_embedding(cells, 2)
  let coords = @src.destiny_plot_coordinates(result, 0, 1)
  
  assert_eq(coords.length(), 10)
}

///|
test "destiny_summary" {
  let cells = @src.create_example_sc_data(8, 3)
  let result = @src.destiny_create_embedding(cells, 2)
  let summary = @src.destiny_summary(result)
  
  assert_true(summary.contains("Diffusion Map Summary"))
  assert_true(summary.contains("Number of cells:"))
}

///|
test "destiny_create_sc_data" {
  let cells = @src.create_example_sc_data(10, 3)
  
  assert_eq(cells.length(), 10)
  assert_eq(cells[0].expression.length(), 3)
  assert_true(cells[0].cluster == "cluster_1" || cells[0].cluster == "cluster_2")
}

///|
test "destiny_euclidean_distance_zero" {
  let cells = @src.create_example_sc_data(3, 2)
  let dist_matrix = @src.compute_distance_matrix(cells, "euclidean")
  
  // Same point should have zero distance
  assert_true(dist_matrix.distances[0][0].abs() < 1.0e-10)
  assert_true(dist_matrix.distances[1][1].abs() < 1.0e-10)
}

///|
test "destiny_small_dataset" {
  let cells = @src.create_example_sc_data(2, 2)
  let result = @src.compute_diffusion_map(cells, 1, 0.5)
  
  assert_eq(result.cell_ids.length(), 2)
  assert_eq(result.embedding.length(), 2)
}
