///|
test "single_cell_experiment_new" {
  let assays = [[[1.0, 2.0], [3.0, 4.0]]]
  let col_data = [["Cell1"], ["Cell2"]]
  let row_data = [["Gene1"], ["Gene2"]]

  let sce = @src.SimpleSCE::new(assays, col_data, row_data)

  assert_eq(sce.assays.length(), 1)
  assert_eq(sce.col_data.length(), 2)
  assert_eq(sce.row_data.length(), 2)
}

///|
test "single_cell_calculate_qc_metrics" {
  let counts = [[10.0, 20.0, 30.0, 5.0], [5.0, 10.0, 15.0, 2.0]]
  let gene_names = ["GeneA", "MT-CO1", "GeneB", "MT-ATP6"]

  let qc = @src.calculate_qc_metrics(counts, gene_names)

  assert_eq(qc.n_umi.length(), 2)
  assert_eq(qc.n_genes.length(), 2)
  assert_eq(qc.mito_percent.length(), 2)

  assert_true(qc.n_umi[0] == 65.0)
  assert_true(qc.n_genes[0] == 4)
}

///|
test "single_cell_filter_cells" {
  let counts = [[100.0, 200.0], [50.0, 100.0], [500.0, 1000.0]]
  let qc = @src.QCMetrics::new([300.0, 150.0, 1500.0], [2, 2, 2], [
    0.0, 0.0, 0.0,
  ])

  let filtered = @src.filter_cells(counts, qc, 100, 1000, 1, 20.0)

  assert_eq(filtered.length(), 2)
}

///|
test "single_cell_log_normalize" {
  let counts = [[10.0, 20.0], [5.0, 10.0]]

  let normalized = @src.log_normalize(counts, None)

  assert_eq(normalized.length(), 2)
  assert_eq(normalized[0].length(), 2)

  assert_true(normalized[0][0] > 0.0)
}

///|
test "single_cell_log_normalize_with_size_factors" {
  let counts = [[10.0, 20.0], [5.0, 10.0]]
  let size_factors = [1.0, 0.5]

  let normalized = @src.log_normalize(counts, Some(size_factors))

  assert_eq(normalized.length(), 2)
}

///|
test "single_cell_pca_power_iteration" {
  let data = [[1.0, 2.0, 3.0], [2.0, 3.0, 4.0], [3.0, 4.0, 5.0]]

  let (pca_result, eigenvalues) = @src.pca_power_iteration(data, 2)

  assert_eq(pca_result.length(), 3)
  assert_eq(pca_result[0].length(), 2)
  assert_eq(eigenvalues.length(), 2)

  assert_true(eigenvalues[0] >= eigenvalues[1])
}

///|
test "single_cell_find_variable_genes" {
  let normalized = [[1.0, 2.0, 10.0], [1.1, 2.1, 20.0], [0.9, 1.9, 15.0]]

  let top_genes = @src.find_variable_genes(normalized, 2)

  assert_eq(top_genes.length(), 2)
}

///|
test "single_cell_select_genes_by_indices" {
  let data = [[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]]
  let indices = [0, 2]

  let selected = @src.select_genes_by_indices(data, indices)

  assert_eq(selected.length(), 2)
  assert_eq(selected[0].length(), 2)
  assert_true(selected[0][0] == 1.0)
  assert_true(selected[0][1] == 3.0)
}

///|
test "single_cell_qc_metrics_new" {
  let n_umi = [1000.0, 2000.0]
  let n_genes = [100, 200]
  let mito_percent = [5.0, 10.0]

  let qc = @src.QCMetrics::new(n_umi, n_genes, mito_percent)

  assert_eq(qc.n_umi.length(), 2)
  assert_eq(qc.n_genes.length(), 2)
  assert_eq(qc.mito_percent.length(), 2)
}
