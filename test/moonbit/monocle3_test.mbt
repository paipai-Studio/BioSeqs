///|
@test "monocle3_new_cell_data_set" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  assert_eq(cds.counts.length(), 50)
  assert_eq(cds.gene_names.length(), 100)
  assert_eq(cds.cell_names.length(), 50)
}

@test "monocle3_preprocess_cds" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  assert_true(preprocessed.normalized.length() > 0)
  assert_true(preprocessed.reduced_dimensions.length() > 0)
}

@test "monocle3_reduce_dimension_pca" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let reduced = @src.reduce_dimension(cds, "PCA", 5)
  assert_true(reduced.reduced_dimensions.length() > 0)
}

@test "monocle3_reduce_dimension_umap" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  assert_true(reduced.reduced_dimensions.length() >= 2)
}

@test "monocle3_learn_graph" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  let with_graph = @src.learn_graph(reduced)
  assert_true(with_graph.cell_partitions.length() > 0)
  assert_true(with_graph.principal_graph.length() > 0)
}

@test "monocle3_order_cells" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  let with_graph = @src.learn_graph(reduced)
  let ordered = @src.order_cells(with_graph, [0])
  assert_eq(ordered.pseudotime.length(), 50)
}

@test "monocle3_fit_models" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  let preprocessed = @src.preprocess_cds(cds, 10)
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  let with_graph = @src.learn_graph(reduced)
  let ordered = @src.order_cells(with_graph, [0])
  let result = @src.fit_models(ordered, ["Gene_0", "Gene_1", "Gene_2"])
  assert_true(result.gene.length() > 0)
}

@test "monocle3_create_example_data" {
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  assert_eq(counts.length(), 50)
  assert_eq(gene_names.length(), 100)
  assert_eq(cell_names.length(), 50)
  assert_eq(counts[0].length(), 100)
}
