///|
/// Tests for the single_cell_experiment module.
/// Mirrors the Bioconductor SingleCellExperiment package API.

///|
test "sce_sqrt_basic" {
  let result = @src.sce_sqrt(4.0)
  assert_true(result > 1.99 && result < 2.01)
}

///|
test "sce_sqrt_zero" {
  assert_eq(@src.sce_sqrt(0.0), 0.0)
}

///|
test "sce_sqrt_one" {
  assert_eq(@src.sce_sqrt(1.0), 1.0)
}

///|
test "sce_ln_basic" {
  // ln(1) = 0
  assert_eq(@src.sce_ln(1.0), 0.0)
  // ln(2.718281828459045) ~= 1.0
  let result = @src.sce_ln(2.718281828459045)
  assert_true(result > 0.99 && result < 1.01)
}

///|
test "sce_ln_non_positive" {
  assert_eq(@src.sce_ln(0.0), -10000000000.0)
}

///|
test "sce_exp_basic" {
  // exp(0) = 1
  assert_eq(@src.sce_exp(0.0), 1.0)
  // exp(1) ~= e
  let result = @src.sce_exp(1.0)
  assert_true(result > 2.71 && result < 2.72)
}

///|
test "sce_log2_basic" {
  // log2(2) = 1, log2(8) = 3
  let r1 = @src.sce_log2(2.0)
  assert_true(r1 > 0.99 && r1 < 1.01)
  let r2 = @src.sce_log2(8.0)
  assert_true(r2 > 2.99 && r2 < 3.01)
}

///|
test "sce_mean_basic" {
  assert_eq(@src.sce_mean([1.0, 2.0, 3.0]), 2.0)
  assert_eq(@src.sce_mean([]), 0.0)
  assert_eq(@src.sce_mean([5.0]), 5.0)
}

///|
test "sce_variance_basic" {
  // Sample variance of [1,2,3,4,5] = 2.5
  let v = @src.sce_variance([1.0, 2.0, 3.0, 4.0, 5.0])
  assert_true(v > 2.49 && v < 2.51)
}

///|
test "sce_variance_singleton" {
  assert_eq(@src.sce_variance([5.0]), 0.0)
  assert_eq(@src.sce_variance([]), 0.0)
}

///|
test "sce_median_basic" {
  assert_eq(@src.sce_median([1.0, 2.0, 3.0]), 2.0)
  // Even-length: median of [1,2,3,4] = 2.5
  assert_eq(@src.sce_median([1.0, 2.0, 3.0, 4.0]), 2.5)
  assert_eq(@src.sce_median([]), 0.0)
}

///|
test "sce_geometric_mean_basic" {
  // geometric mean of [1, 4] = 2.
  let gm = @src.sce_geometric_mean([1.0, 4.0])
  assert_true(gm > 1.99 && gm < 2.01)
}

///|
test "sce_geometric_mean_with_zero" {
  // Zero is ignored.
  let gm = @src.sce_geometric_mean([0.0, 4.0, 9.0])
  // geometric mean of [4,9] = 6.
  assert_true(gm > 5.99 && gm < 6.01)
}

///|
test "sce_geometric_mean_all_nonpositive" {
  // All non-positive -> 0.
  assert_eq(@src.sce_geometric_mean([0.0, -1.0]), 0.0)
}

///|
test "sce_create_example_counts" {
  let counts = @src.sce_create_example_counts()
  assert_eq(counts.length(), 5)
  assert_eq(counts[0].length(), 4)
  assert_eq(counts[0][0], 100.0)
}

///|
test "sce_create_example" {
  let sce = @src.sce_create_example()
  // 5 genes x 4 cells.
  assert_eq(sce.row_names.length(), 5)
  assert_eq(sce.col_names.length(), 4)
  assert_eq(sce.row_names[0], "GeneA")
  assert_eq(sce.col_names[0], "Cell1")
}

///|
test "sce_single_cell_experiment_new" {
  let counts = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let row_names = ["GeneA", "GeneB", "GeneC"]
  let col_names = ["Cell1", "Cell2"]
  let sce = @src.SingleCellExperiment::new(counts, row_names, col_names)
  assert_eq(sce.row_names.length(), 3)
  assert_eq(sce.col_names.length(), 2)
}

///|
test "sce_builder_new_and_build" {
  let counts = [[1.0, 2.0], [3.0, 4.0]]
  let sce = @src.SCEBuilder::new()
    |> @src.SCEBuilder::add_assay("counts", counts)
    |> @src.SCEBuilder::add_row_data("gene_type", ["protein", "ncRNA"])
    |> @src.SCEBuilder::add_col_data("batch", ["b1", "b2"])
    |> @src.SCEBuilder::build
  assert_eq(sce.row_names.length(), 2)
  assert_eq(sce.col_names.length(), 2)
  assert_eq(sce.row_names[0], "Gene1")
  assert_eq(sce.col_names[0], "Cell1")
}

///|
test "sce_get_assay_returns_counts" {
  let sce = @src.sce_create_example()
  let counts = @src.sce_get_assay(sce, "counts")
  assert_eq(counts.length(), 5)
  assert_eq(counts[0][0], 100.0)
}

///|
test "sce_get_assay_missing_returns_empty" {
  let sce = @src.sce_create_example()
  let missing = @src.sce_get_assay(sce, "nonexistent")
  assert_eq(missing.length(), 0)
}

///|
test "sce_set_assay_replaces" {
  let sce = @src.sce_create_example()
  let new_data = [[10.0, 20.0], [30.0, 40.0]]
  let _ = @src.sce_set_assay(sce, "custom", new_data)
  let retrieved = @src.sce_get_assay(sce, "custom")
  assert_eq(retrieved.length(), 2)
  assert_eq(retrieved[0][0], 10.0)
}

///|
test "sce_get_counts" {
  let sce = @src.sce_create_example()
  let counts = @src.sce_get_counts(sce)
  assert_eq(counts.length(), 5)
  assert_eq(counts[0].length(), 4)
}

///|
test "sce_get_logcounts_and_set_logcounts" {
  let sce = @src.sce_create_example()
  // Example already has logcounts computed.
  let log1 = @src.sce_get_logcounts(sce)
  assert_eq(log1.length(), 5)
  // Replace and verify.
  let new_log = [[0.5, 0.6], [0.7, 0.8]]
  let _ = @src.sce_set_logcounts(sce, new_log)
  let log2 = @src.sce_get_logcounts(sce)
  assert_eq(log2.length(), 2)
  assert_eq(log2[0][0], 0.5)
}

///|
test "sce_compute_logcounts_basic" {
  let counts = [[10.0, 20.0], [30.0, 40.0]]
  let size_factors = [1.0, 1.0]
  let logcounts = @src.sce_compute_logcounts(counts, size_factors)
  assert_eq(logcounts.length(), 2)
  assert_eq(logcounts[0].length(), 2)
  // log2(10+1) ~= 3.459
  assert_true(logcounts[0][0] > 3.4 && logcounts[0][0] < 3.5)
}

///|
test "sce_compute_logcounts_zero_size_factor" {
  let counts = [[10.0, 20.0]]
  let size_factors = [0.0, 0.0]
  // With zero size factors, function defaults to 1.0.
  let logcounts = @src.sce_compute_logcounts(counts, size_factors)
  assert_eq(logcounts.length(), 1)
  // log2(11) ~= 3.459
  assert_true(logcounts[0][0] > 3.4 && logcounts[0][0] < 3.5)
}

///|
test "sce_get_reduced_dim_missing" {
  let sce = @src.sce_create_example()
  let missing = @src.sce_get_reduced_dim(sce, "PCA")
  assert_eq(missing.length(), 0)
}

///|
test "sce_set_reduced_dim" {
  let sce = @src.sce_create_example()
  let coords = [[0.1, 0.2], [0.3, 0.4], [0.5, 0.6], [0.7, 0.8]]
  let _ = @src.sce_set_reduced_dim(sce, "PCA", coords)
  let retrieved = @src.sce_get_reduced_dim(sce, "PCA")
  assert_eq(retrieved.length(), 4)
  assert_eq(retrieved[0][0], 0.1)
}

///|
test "sce_compute_pca_basic" {
  let logcounts = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [4.0, 3.0, 2.0, 1.0],
  ]
  let result = @src.sce_compute_pca(logcounts, 2)
  assert_eq(result.method, "PCA")
  assert_eq(result.n_components, 2)
  assert_eq(result.coordinates.length(), 4)
  assert_eq(result.variance_explained.length(), 2)
  // Variance explained sums to <= 1.
  let mut total = 0.0
  let mut i = 0
  while i < result.variance_explained.length() {
    total = total + result.variance_explained[i]
    i = i + 1
  }
  assert_true(total >= 0.0 && total <= 1.01)
}

///|
test "sce_compute_pca_empty" {
  let result = @src.sce_compute_pca([], 2)
  assert_eq(result.method, "PCA")
  assert_eq(result.n_components, 0)
  assert_eq(result.coordinates.length(), 0)
}

///|
test "sce_compute_tsne_basic" {
  let logcounts = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [4.0, 3.0, 2.0, 1.0],
  ]
  let result = @src.sce_compute_tsne(logcounts, 5)
  assert_eq(result.method, "tSNE")
  assert_eq(result.n_components, 2)
  assert_eq(result.coordinates.length(), 4)
}

///|
test "sce_compute_tsne_empty" {
  let result = @src.sce_compute_tsne([], 5)
  assert_eq(result.method, "tSNE")
  assert_eq(result.n_components, 0)
}

///|
test "sce_compute_umap_basic" {
  let logcounts = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [4.0, 3.0, 2.0, 1.0],
  ]
  let result = @src.sce_compute_umap(logcounts, 2)
  assert_eq(result.method, "UMAP")
  assert_eq(result.n_components, 2)
  assert_eq(result.coordinates.length(), 4)
}

///|
test "sce_compute_library_size_factors_basic" {
  // Two genes, three cells:
  // col0: 10+20=30, col1: 20+30=50, col2: 30+40=70
  // mean = 50; factors = [0.6, 1.0, 1.4]
  let counts = [[10.0, 20.0, 30.0], [20.0, 30.0, 40.0]]
  let factors = @src.sce_compute_library_size_factors(counts)
  assert_eq(factors.length(), 3)
  assert_true(factors[0] > 0.59 && factors[0] < 0.61)
  assert_true(factors[1] > 0.99 && factors[1] < 1.01)
  assert_true(factors[2] > 1.39 && factors[2] < 1.41)
}

///|
test "sce_compute_deconv_size_factors_basic" {
  let counts = [[10.0, 20.0, 30.0], [20.0, 30.0, 40.0]]
  let factors = @src.sce_compute_deconv_size_factors(counts)
  assert_eq(factors.length(), 3)
  // All factors should be positive.
  let mut i = 0
  while i < factors.length() {
    assert_true(factors[i] > 0.0)
    i = i + 1
  }
}

///|
test "sce_compute_deconv_size_factors_empty" {
  let factors = @src.sce_compute_deconv_size_factors([])
  assert_eq(factors.length(), 0)
}

///|
test "sce_compute_size_factors_library" {
  let counts = [[10.0, 20.0, 30.0], [20.0, 30.0, 40.0]]
  let result = @src.sce_compute_size_factors(counts, "library")
  assert_eq(result.method, "library")
  assert_eq(result.size_factors.length(), 3)
  assert_eq(result.log_means.length(), 3)
}

///|
test "sce_compute_size_factors_deconv" {
  let counts = [[10.0, 20.0, 30.0], [20.0, 30.0, 40.0]]
  let result = @src.sce_compute_size_factors(counts, "deconv")
  assert_eq(result.method, "deconv")
  assert_eq(result.size_factors.length(), 3)
}

///|
test "sce_normalize_library" {
  let sce = @src.sce_create_example()
  let _ = @src.sce_normalize(sce, "library")
  let logcounts = @src.sce_get_logcounts(sce)
  assert_eq(logcounts.length(), 5)
  assert_eq(logcounts[0].length(), 4)
  // logcounts should be non-negative.
  assert_true(logcounts[0][0] >= 0.0)
}

///|
test "sce_normalize_deconv" {
  let sce = @src.sce_create_example()
  let _ = @src.sce_normalize(sce, "deconv")
  let logcounts = @src.sce_get_logcounts(sce)
  assert_eq(logcounts.length(), 5)
  assert_eq(logcounts[0].length(), 4)
}

///|
test "sce_filter_genes_basic" {
  // 3 genes, 3 cells. Gene 0 is expressed in all 3 cells.
  let counts = [[1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 1.0]]
  let sce = @src.SingleCellExperiment::new(counts, ["g1", "g2", "g3"], [
    "c1", "c2", "c3",
  ])
  let filtered = @src.sce_filter_genes(sce, 2)
  // Only gene 0 (3 expressed) and gene 2 (2 expressed) pass min_cells=2.
  assert_eq(filtered.row_names.length(), 2)
  assert_eq(filtered.row_names[0], "g1")
  assert_eq(filtered.row_names[1], "g3")
}

///|
test "sce_filter_cells_basic" {
  // 3 genes, 3 cells. Cell 1 is empty.
  let counts = [[1.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 0.0, 1.0]]
  let sce = @src.SingleCellExperiment::new(counts, ["g1", "g2", "g3"], [
    "c1", "c2", "c3",
  ])
  let filtered = @src.sce_filter_cells(sce, 2)
  // Cells 0 and 2 express 3 genes (>= 2), cell 1 expresses 0.
  assert_eq(filtered.col_names.length(), 2)
  assert_eq(filtered.col_names[0], "c1")
  assert_eq(filtered.col_names[1], "c3")
}

///|
test "sce_get_row_data_missing" {
  let sce = @src.sce_create_example()
  let result = @src.sce_get_row_data(sce, "nonexistent_column")
  assert_eq(result.length(), 0)
}

///|
test "sce_get_col_data_missing" {
  let sce = @src.sce_create_example()
  let result = @src.sce_get_col_data(sce, "nonexistent_column")
  assert_eq(result.length(), 0)
}

///|
test "sce_subset_genes_by_index" {
  let sce = @src.sce_create_example()
  // Keep only genes 0 and 2.
  let subset = @src.sce_subset_genes(sce, [0, 2])
  assert_eq(subset.row_names.length(), 2)
  assert_eq(subset.row_names[0], "GeneA")
  assert_eq(subset.row_names[1], "GeneC")
  let counts = @src.sce_get_counts(subset)
  assert_eq(counts.length(), 2)
}

///|
test "sce_subset_cells_by_index" {
  let sce = @src.sce_create_example()
  // Keep only cells 1 and 3.
  let subset = @src.sce_subset_cells(sce, [1, 3])
  assert_eq(subset.col_names.length(), 2)
  assert_eq(subset.col_names[0], "Cell2")
  assert_eq(subset.col_names[1], "Cell4")
  let counts = @src.sce_get_counts(subset)
  assert_eq(counts[0].length(), 2)
}

///|
test "sce_get_n_genes" {
  let sce = @src.sce_create_example()
  assert_eq(@src.sce_get_n_genes(sce), 5)
}

///|
test "sce_get_n_cells" {
  let sce = @src.sce_create_example()
  assert_eq(@src.sce_get_n_cells(sce), 4)
}

///|
test "sce_get_n_genes_empty" {
  let sce = @src.SingleCellExperiment::new([], [], [])
  assert_eq(@src.sce_get_n_genes(sce), 0)
  assert_eq(@src.sce_get_n_cells(sce), 0)
}

///|
test "sce_merge_basic" {
  let sce1 = @src.SingleCellExperiment::new(
    [[1.0, 2.0], [3.0, 4.0]],
    ["g1", "g2"],
    ["c1", "c2"],
  )
  let sce2 = @src.SingleCellExperiment::new([[5.0, 6.0]], ["g3"], ["c1", "c2"])
  let merged = @src.sce_merge([sce1, sce2])
  assert_eq(merged.row_names.length(), 3)
  assert_eq(merged.row_names[0], "g1")
  assert_eq(merged.row_names[2], "g3")
  let counts = @src.sce_get_counts(merged)
  assert_eq(counts.length(), 3)
}

///|
test "sce_merge_single_returns_same" {
  let sce = @src.sce_create_example()
  let merged = @src.sce_merge([sce])
  assert_eq(merged.row_names.length(), 5)
}

///|
test "sce_merge_empty" {
  let merged = @src.sce_merge([])
  assert_eq(merged.row_names.length(), 0)
}

///|
test "sce_create_example_with_reduced_dims" {
  let sce = @src.sce_create_example_with_reduced_dims()
  let pca = @src.sce_get_reduced_dim(sce, "PCA")
  assert_eq(pca.length(), 4)
  let tsne = @src.sce_get_reduced_dim(sce, "tSNE")
  assert_eq(tsne.length(), 4)
  let umap = @src.sce_get_reduced_dim(sce, "UMAP")
  assert_eq(umap.length(), 4)
}

///|
test "sce_to_string" {
  let sce = @src.sce_create_example()
  let s = @src.SingleCellExperiment::to_string(sce)
  assert_true(s.contains("SingleCellExperiment"))
  assert_true(s.contains("Genes: 5"))
  assert_true(s.contains("Cells: 4"))
  assert_true(s.contains("GeneA"))
  assert_true(s.contains("Cell1"))
}

///|
test "sce_summary" {
  let sce = @src.sce_create_example()
  let s = @src.sce_summary(sce)
  assert_true(s.contains("SingleCellExperiment"))
  assert_true(s.contains("Genes: 5"))
}

///|
test "sce_reduced_dim_result_to_string" {
  let logcounts = [[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0]]
  let result = @src.sce_compute_pca(logcounts, 2)
  let s = @src.ReducedDimResult::to_string(result)
  assert_true(s.contains("ReducedDimResult"))
  assert_true(s.contains("Method: PCA"))
  assert_true(s.contains("Components: 2"))
}

///|
test "sce_size_factor_result_to_string" {
  let counts = [[10.0, 20.0], [30.0, 40.0]]
  let result = @src.sce_compute_size_factors(counts, "library")
  let s = @src.SizeFactorResult::to_string(result)
  assert_true(s.contains("SizeFactorResult"))
  assert_true(s.contains("Method: library"))
}

///|
test "sce_builder_with_reduced_dim" {
  let counts = [[1.0, 2.0], [3.0, 4.0]]
  let rd = [[0.1, 0.2], [0.3, 0.4]]
  let sce = @src.SCEBuilder::new()
    |> @src.SCEBuilder::add_assay("counts", counts)
    |> @src.SCEBuilder::add_reduced_dim("PCA", rd)
    |> @src.SCEBuilder::build
  let pca = @src.sce_get_reduced_dim(sce, "PCA")
  assert_eq(pca.length(), 2)
  assert_eq(pca[0][0], 0.1)
}

///|
test "sce_builder_default_names" {
  // When row_names and col_names are not set, builder generates defaults.
  let counts = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let sce = @src.SCEBuilder::new()
    |> @src.SCEBuilder::add_assay("counts", counts)
    |> @src.SCEBuilder::build
  assert_eq(sce.row_names.length(), 2)
  assert_eq(sce.row_names[0], "Gene1")
  assert_eq(sce.row_names[1], "Gene2")
  assert_eq(sce.col_names.length(), 3)
  assert_eq(sce.col_names[0], "Cell1")
  assert_eq(sce.col_names[2], "Cell3")
}
