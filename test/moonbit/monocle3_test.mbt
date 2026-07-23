test "monocle3_new_cell_data_set" {
  let counts = [[1, 2], [3, 4]]
  let gene_names = ["Gene_A", "Gene_B"]
  let cell_names = ["Cell_1", "Cell_2"]
  
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  
  assert_eq(cds.counts.length(), 2)
  assert_eq(cds.gene_names.length(), 2)
  assert_eq(cds.cell_names.length(), 2)
  assert_eq(cds.counts[0][0], 1)
}

test "monocle3_new_cell_data_set_empty" {
  let cds = @src.new_cell_data_set([], [], [])
  
  assert_eq(cds.counts.length(), 0)
  assert_eq(cds.normalized.length(), 0)
}

test "monocle3_preprocess_cds" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  
  let preprocessed = @src.preprocess_cds(cds, 10)
  
  assert_eq(preprocessed.counts.length(), 50)
  assert_eq(preprocessed.normalized.length(), 50)
  assert_eq(preprocessed.reduced_dimensions.length(), 1)
}

test "monocle3_reduce_dimension_pca" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  
  let reduced = @src.reduce_dimension(cds, "PCA", 5)
  
  assert_eq(reduced.reduced_dimensions.length(), 1)
  assert_eq(reduced.reduced_dimensions[0].0, "PCA")
}

test "monocle3_reduce_dimension_umap" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  
  assert_true(reduced.reduced_dimensions.length() >= 2)
}

test "monocle3_learn_graph" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  
  let with_graph = @src.learn_graph(reduced)
  
  assert_eq(with_graph.cell_partitions.length(), 50)
  assert_eq(with_graph.principal_graph.length(), 1)
}

test "monocle3_order_cells" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  let with_graph = @src.learn_graph(reduced)
  
  let ordered = @src.order_cells(with_graph, [0])
  
  assert_eq(ordered.pseudotime.length(), 50)
}

test "monocle3_fit_models" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  let with_graph = @src.learn_graph(reduced)
  let ordered = @src.order_cells(with_graph, [0])
  
  let result = @src.fit_models(ordered, ["Gene_0", "Gene_1", "Gene_2"])
  
  assert_eq(result.gene.length(), 3)
  assert_eq(result.q_value.length(), 3)
}

test "monocle3_create_example_data" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  
  assert_eq(counts.length(), 50)
  assert_eq(gene_names.length(), 100)
  assert_eq(cell_names.length(), 50)
}

test "monocle3_full_pipeline" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  let with_graph = @src.learn_graph(reduced)
  let ordered = @src.order_cells(with_graph, [0])
  let result = @src.fit_models(ordered, ["Gene_0"])
  
  assert_eq(cds.counts.length(), 50)
  assert_eq(result.gene.length(), 1)
}

test "monocle3_find_branch_points" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  let with_graph = @src.learn_graph(reduced)
  
  let branch_points = @src.find_branch_points(with_graph)
  
  assert_eq(branch_points.length(), 0)
}

test "monocle3_differential_gene_test_branches" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  let with_graph = @src.learn_graph(reduced)
  
  let branch_points = @src.find_branch_points(with_graph)
  
  if branch_points.length() > 0 {
    let result = @src.differential_gene_test_branches(with_graph, branch_points[0], ["Gene_0", "Gene_1"])
    assert_eq(result.gene.length(), 2)
  }
}