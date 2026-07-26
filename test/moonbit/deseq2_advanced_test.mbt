///|
/// DESeq2 Advanced tests - VST (Variance Stabilizing Transformation) and PCA.

///|
test "vst_transformation_correct_dimensions" {
  let dds = @src.sample_deseq_dataset()
  let result = @src.variance_stabilizing_transformation(dds)
  assert_eq(result.transformed.length(), 20)
  assert_eq(result.transformed[0].length(), 6)
  assert_eq(result.row_names.length(), 20)
  assert_eq(result.col_names.length(), 6)
}

///|
test "vst_transformation_blind_true" {
  let dds = @src.sample_deseq_dataset()
  let result = @src.variance_stabilizing_transformation(dds, blind=true)
  assert_eq(result.transformed.length(), 20)
  assert_eq(result.transformed[0].length(), 6)
  assert_true(result.dispersion > 0.0)
}

///|
test "vst_transformation_blind_false" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  let dds_disp = @src.estimate_dispersions(dds_sf)
  let result = @src.variance_stabilizing_transformation(dds_disp, blind=false)
  assert_eq(result.transformed.length(), 20)
  assert_eq(result.transformed[0].length(), 6)
  assert_true(result.dispersion > 0.0)
}

///|
test "fast_vst_correct_dimensions" {
  let dds = @src.sample_deseq_dataset()
  let result = @src.vst(dds)
  assert_eq(result.transformed.length(), 20)
  assert_eq(result.transformed[0].length(), 6)
  assert_eq(result.row_names.length(), 20)
  assert_eq(result.col_names.length(), 6)
}

///|
test "fast_vst_custom_nsub" {
  let dds = @src.sample_deseq_dataset()
  let result = @src.vst(dds, nsub=5)
  assert_eq(result.transformed.length(), 20)
  assert_eq(result.transformed[0].length(), 6)
  assert_true(result.dispersion > 0.0)
}

///|
test "vst_transformed_values_nonnegative_for_positive_counts" {
  let dds = @src.DESeqDataSet::new(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250]],
    ["GeneA", "GeneB", "GeneC"],
    ["S1", "S2", "S3", "S4"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
  )
  let result = @src.variance_stabilizing_transformation(dds)
  for i = 0; i < result.transformed.length(); i = i + 1 {
    for j = 0; j < result.transformed[i].length(); j = j + 1 {
      assert_true(result.transformed[i][j] >= 0.0)
    }
  }
}

///|
test "vst_dispersion_positive" {
  let dds = @src.sample_deseq_dataset()
  let result = @src.variance_stabilizing_transformation(dds)
  assert_true(result.dispersion > 0.0)
}

///|
test "vst_zero_counts_produce_zero_transformed_values" {
  let dds = @src.DESeqDataSet::new(
    [[0, 0, 0], [0, 0, 0]],
    ["GeneA", "GeneB"],
    ["S1", "S2", "S3"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
  )
  let result = @src.variance_stabilizing_transformation(dds)
  for i = 0; i < result.transformed.length(); i = i + 1 {
    for j = 0; j < result.transformed[i].length(); j = j + 1 {
      assert_true(result.transformed[i][j] == 0.0)
    }
  }
}

///|
test "pca_correct_sample_count" {
  let dds = @src.sample_deseq_dataset()
  let vst_result = @src.variance_stabilizing_transformation(dds)
  let pca = @src.plot_pca(vst_result)
  assert_eq(pca.sample_names.length(), 6)
  assert_eq(pca.pc1.length(), 6)
  assert_eq(pca.pc2.length(), 6)
}

///|
test "pca_variance_explained_sums_le_100" {
  let dds = @src.sample_deseq_dataset()
  let vst_result = @src.variance_stabilizing_transformation(dds)
  let pca = @src.plot_pca(vst_result)
  let total = pca.var_explained_pc1 + pca.var_explained_pc2
  assert_true(total <= 100.0)
  assert_true(pca.var_explained_pc1 >= 0.0)
  assert_true(pca.var_explained_pc2 >= 0.0)
}

///|
test "pca_custom_ntop" {
  let dds = @src.sample_deseq_dataset()
  let vst_result = @src.variance_stabilizing_transformation(dds)
  let pca = @src.plot_pca(vst_result, ntop=10)
  assert_eq(pca.sample_names.length(), 6)
  assert_eq(pca.pc1.length(), 6)
  assert_eq(pca.pc2.length(), 6)
  assert_true(pca.var_explained_pc1 >= 0.0)
  assert_true(pca.var_explained_pc2 >= 0.0)
}

///|
test "get_transformed_matrix_returns_correct_dimensions" {
  let dds = @src.sample_deseq_dataset()
  let vst_result = @src.variance_stabilizing_transformation(dds)
  let matrix = @src.get_transformed_matrix(vst_result)
  assert_eq(matrix.length(), 20)
  assert_eq(matrix[0].length(), 6)
}

///|
test "gene_sd_returns_correct_number_of_genes" {
  let dds = @src.sample_deseq_dataset()
  let vst_result = @src.variance_stabilizing_transformation(dds)
  let sds = @src.gene_sd(vst_result)
  assert_eq(sds.length(), 20)
}

///|
test "gene_sd_values_are_nonnegative" {
  let dds = @src.sample_deseq_dataset()
  let vst_result = @src.variance_stabilizing_transformation(dds)
  let sds = @src.gene_sd(vst_result)
  for entry in sds {
    assert_true(entry.1 >= 0.0)
  }
}

///|
test "vst_handles_single_sample_edge_case" {
  let dds = @src.DESeqDataSet::new(
    [[10], [20], [30]],
    ["GeneA", "GeneB", "GeneC"],
    ["S1"],
    [[1.0, 0.0]],
  )
  let result = @src.variance_stabilizing_transformation(dds)
  assert_eq(result.transformed.length(), 3)
  assert_eq(result.transformed[0].length(), 1)
  assert_eq(result.col_names.length(), 1)
}

///|
test "vst_handles_single_gene_edge_case" {
  let dds = @src.DESeqDataSet::new(
    [[10, 20, 30]],
    ["GeneA"],
    ["S1", "S2", "S3"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
  )
  let result = @src.variance_stabilizing_transformation(dds)
  assert_eq(result.transformed.length(), 1)
  assert_eq(result.transformed[0].length(), 3)
  assert_eq(result.row_names.length(), 1)
}
