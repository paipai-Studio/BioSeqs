///|
test "mast_data_new_basic" {
  let expression = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let cell_names = ["Cell1", "Cell2", "Cell3"]
  let groups = [0, 0, 1]
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  assert_eq(data.n_cells, 3)
  assert_eq(data.n_genes, 3)
  assert_eq(data.n_groups, 2)
}

///|
test "mast_data_new_cngeneson" {
  let expression = [[1.0, 0.0, 3.0], [0.0, 0.0, 0.0], [7.0, 8.0, 0.0]]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let cell_names = ["Cell1", "Cell2", "Cell3"]
  let groups = [0, 0, 1]
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  assert_true(data.cngeneson[0] > 0.5)
  assert_eq(data.cngeneson[1], 0.0)
  assert_true(data.cngeneson[2] > 0.5)
}

///|
test "mast_data_new_empty" {
  let expression : Array[Array[Double]] = []
  let gene_names : Array[String] = []
  let cell_names : Array[String] = []
  let groups : Array[Int] = []
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  assert_eq(data.n_cells, 0)
  assert_eq(data.n_genes, 0)
  assert_eq(data.n_groups, 1)
}

///|
test "mast_hurdle_result_new" {
  let result = @src.HurdleResult::new(5)
  assert_eq(result.gene.length(), 5)
  assert_eq(result.discrete_pval.length(), 5)
  assert_eq(result.continuous_pval.length(), 5)
  assert_eq(result.combined_pval.length(), 5)
  assert_eq(result.discrete_fdr.length(), 5)
  assert_eq(result.continuous_fdr.length(), 5)
  assert_eq(result.combined_fdr.length(), 5)
  assert_eq(result.log2_fc.length(), 5)
  assert_eq(result.detection_rate.length(), 5)
}

///|
test "mast_zlm_hurdle_basic" {
  let expression = [
    [10.0, 5.0, 0.0, 8.0],
    [12.0, 6.0, 1.0, 9.0],
    [8.0, 4.0, 0.0, 7.0],
    [100.0, 50.0, 10.0, 80.0],
    [120.0, 60.0, 12.0, 90.0],
    [90.0, 45.0, 9.0, 75.0],
  ]
  let gene_names = ["Gene1", "Gene2", "Gene3", "Gene4"]
  let cell_names = ["C1", "C2", "C3", "C4", "C5", "C6"]
  let groups = [0, 0, 0, 1, 1, 1]
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  let result = @src.zlm_hurdle(data, 0, 1)
  assert_eq(result.gene.length(), 4)
  assert_eq(result.discrete_pval.length(), 4)
  assert_eq(result.combined_pval.length(), 4)
}

///|
test "mast_zlm_hurdle_pval_range" {
  let expression = [
    [10.0, 5.0, 0.0],
    [12.0, 6.0, 1.0],
    [8.0, 4.0, 0.0],
    [100.0, 50.0, 10.0],
    [120.0, 60.0, 12.0],
    [90.0, 45.0, 9.0],
  ]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let cell_names = ["C1", "C2", "C3", "C4", "C5", "C6"]
  let groups = [0, 0, 0, 1, 1, 1]
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  let result = @src.zlm_hurdle(data, 0, 1)
  let mut i = 0
  while i < result.gene.length() {
    assert_true(
      result.discrete_pval[i] >= 0.0 && result.discrete_pval[i] <= 1.0,
    )
    assert_true(
      result.continuous_pval[i] >= 0.0 && result.continuous_pval[i] <= 1.0,
    )
    assert_true(
      result.combined_pval[i] >= 0.0 && result.combined_pval[i] <= 1.0,
    )
    assert_true(result.discrete_fdr[i] >= 0.0 && result.discrete_fdr[i] <= 1.0)
    assert_true(result.combined_fdr[i] >= 0.0 && result.combined_fdr[i] <= 1.0)
    i = i + 1
  }
}

///|
test "mast_zlm_hurdle_detection_rate" {
  let expression = [
    [10.0, 0.0, 0.0],
    [12.0, 0.0, 0.0],
    [8.0, 0.0, 0.0],
    [100.0, 50.0, 10.0],
    [120.0, 60.0, 12.0],
    [90.0, 45.0, 9.0],
  ]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let cell_names = ["C1", "C2", "C3", "C4", "C5", "C6"]
  let groups = [0, 0, 0, 1, 1, 1]
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  let result = @src.zlm_hurdle(data, 0, 1)
  assert_true(result.detection_rate[1] > 0.0 && result.detection_rate[1] < 1.0)
}

///|
test "mast_summary_new" {
  let expression = [
    [10.0, 5.0, 0.0],
    [12.0, 6.0, 1.0],
    [8.0, 4.0, 0.0],
    [100.0, 50.0, 10.0],
    [120.0, 60.0, 12.0],
    [90.0, 45.0, 9.0],
  ]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let cell_names = ["C1", "C2", "C3", "C4", "C5", "C6"]
  let groups = [0, 0, 0, 1, 1, 1]
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  let result = @src.zlm_hurdle(data, 0, 1)
  let summary = @src.MastSummary::new(result, 0.05)
  assert_eq(summary.threshold, 0.05)
  assert_true(summary.n_de_combined >= 0)
  assert_true(summary.n_de_discrete >= 0)
  assert_true(summary.n_de_continuous >= 0)
}

///|
test "mast_summary_top_genes" {
  let expression = [
    [10.0, 5.0, 0.0],
    [12.0, 6.0, 1.0],
    [8.0, 4.0, 0.0],
    [100.0, 50.0, 10.0],
    [120.0, 60.0, 12.0],
    [90.0, 45.0, 9.0],
  ]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let cell_names = ["C1", "C2", "C3", "C4", "C5", "C6"]
  let groups = [0, 0, 0, 1, 1, 1]
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  let result = @src.zlm_hurdle(data, 0, 1)
  let summary = @src.MastSummary::new(result, 0.05)
  assert_eq(summary.top_genes.length(), 3)
}

///|
test "mast_hurdle_result_to_string" {
  let result = @src.HurdleResult::new(3)
  let str = result.to_string()
  assert_true(str.length() > 0)
  assert_true(str.contains("HurdleResult"))
  assert_true(str.contains("Total genes"))
}

///|
test "mast_summary_to_string" {
  let expression = [
    [10.0, 5.0, 0.0],
    [12.0, 6.0, 1.0],
    [8.0, 4.0, 0.0],
    [100.0, 50.0, 10.0],
    [120.0, 60.0, 12.0],
    [90.0, 45.0, 9.0],
  ]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let cell_names = ["C1", "C2", "C3", "C4", "C5", "C6"]
  let groups = [0, 0, 0, 1, 1, 1]
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  let result = @src.zlm_hurdle(data, 0, 1)
  let summary = @src.MastSummary::new(result, 0.05)
  let str = summary.to_string()
  assert_true(str.length() > 0)
  assert_true(str.contains("MAST Summary"))
  assert_true(str.contains("DE genes"))
}

///|
test "mast_single_group" {
  let expression = [[1.0, 2.0], [3.0, 4.0]]
  let gene_names = ["Gene1", "Gene2"]
  let cell_names = ["C1", "C2"]
  let groups = [0, 0]
  let data = @src.MastData::new(expression, gene_names, cell_names, groups)
  assert_eq(data.n_groups, 1)
}
