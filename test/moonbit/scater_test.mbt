///|
test "scater_calculate_qc_metrics_basic" {
  let counts = [[10, 20, 30, 0], [5, 15, 25, 5], [0, 0, 0, 0]]
  let cell_ids = ["cell1", "cell2", "cell3"]
  let gene_names = ["Gene1", "Gene2", "Gene3", "mt-Atp6"]
  let stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  assert_eq(stats.cell_id.length(), 3)
  assert_eq(stats.total_counts.length(), 3)
  assert_eq(stats.total_features.length(), 3)
  assert_eq(stats.percent_mito.length(), 3)
  assert_eq(stats.detected_genes.length(), 3)
}

///|
test "scater_calculate_qc_metrics_total_counts" {
  let counts = [[10, 20, 30], [5, 15, 25]]
  let cell_ids = ["cell1", "cell2"]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  assert_eq(stats.total_counts[0], 60.0)
  assert_eq(stats.total_counts[1], 45.0)
}

///|
test "scater_calculate_qc_metrics_mito" {
  let counts = [[10, 90], [50, 50]]
  let cell_ids = ["cell1", "cell2"]
  let gene_names = ["Gene1", "mt-Atp6"]
  let stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  assert_true(stats.percent_mito[0] > 80.0)
  assert_true(stats.percent_mito[1] > 40.0 && stats.percent_mito[1] < 60.0)
}

///|
test "scater_calculate_qc_metrics_ribo" {
  let counts = [[10, 90], [50, 50]]
  let cell_ids = ["cell1", "cell2"]
  let gene_names = ["Gene1", "Rpl10"]
  let stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  assert_true(stats.percent_ribo[0] > 80.0)
}

///|
test "scater_calculate_qc_metrics_detected" {
  let counts = [[10, 0, 30, 0], [5, 15, 0, 0]]
  let cell_ids = ["cell1", "cell2"]
  let gene_names = ["Gene1", "Gene2", "Gene3", "Gene4"]
  let stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  assert_eq(stats.detected_genes[0], 2)
  assert_eq(stats.detected_genes[1], 2)
}

///|
test "scater_filter_cells_basic" {
  let counts = [[10, 20, 30], [5, 5, 5], [100, 200, 300]]
  let cell_ids = ["cell1", "cell2", "cell3"]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  let params = @src.QCFilterParams::new()
  let (kept_indices, keep) = @src.scater_filter_cells(stats, params)
  assert_eq(keep.length(), 3)
}

///|
test "scater_filter_cells_custom_params" {
  let counts = [[10, 20, 30], [5, 5, 5], [100, 200, 300]]
  let cell_ids = ["cell1", "cell2", "cell3"]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  let params = @src.QCFilterParams::with_params(
    50.0, 
    1000.0, 
    2, 
    100, 
    20.0,
    2
  )
  let (kept_indices, keep) = @src.scater_filter_cells(stats, params)
  assert_true(kept_indices.length() >= 1)
}

///|
test "scater_filter_counts" {
  let counts = [[1, 2], [3, 4], [5, 6]]
  let keep = [true, false, true]
  let filtered = @src.scater_filter_counts(counts, keep)
  assert_eq(filtered.length(), 2)
  assert_eq(filtered[0][0], 1)
  assert_eq(filtered[1][0], 5)
}

///|
test "scater_calculate_feature_qc" {
  let counts = [[10, 20, 0], [5, 15, 0], [8, 0, 0]]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let stats = @src.scater_calculate_feature_qc(counts, gene_names)
  assert_eq(stats.gene_id.length(), 3)
  assert_eq(stats.mean_expression.length(), 3)
  assert_eq(stats.n_cells_detected.length(), 3)
  assert_eq(stats.percent_detected.length(), 3)
}

///|
test "scater_calculate_feature_qc_detected" {
  let counts = [[10, 20, 0], [5, 15, 0], [8, 0, 0]]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let stats = @src.scater_calculate_feature_qc(counts, gene_names)
  assert_eq(stats.n_cells_detected[0], 3)
  assert_eq(stats.n_cells_detected[1], 2)
  assert_eq(stats.n_cells_detected[2], 0)
}

///|
test "scater_filter_features" {
  let counts = [[10, 20, 0], [5, 15, 0], [8, 0, 0]]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  let stats = @src.scater_calculate_feature_qc(counts, gene_names)
  let (kept_indices, keep) = @src.scater_filter_features(stats, 2)
  assert_eq(keep.length(), 3)
  assert_eq(kept_indices.length(), 2)
}

///|
test "scater_normalize_cpm" {
  let counts = [[10, 20, 30], [5, 15, 25]]
  let cpm = @src.scater_normalize_cpm(counts)
  assert_eq(cpm.length(), 2)
  assert_eq(cpm[0].length(), 3)
  let total_cpm = cpm[0][0] + cpm[0][1] + cpm[0][2]
  assert_true(total_cpm > 999999.0 && total_cpm < 1000001.0)
}

///|
test "scater_log_normalize_cpm" {
  let counts = [[10, 20, 30], [5, 15, 25]]
  let log_cpm = @src.scater_log_normalize_cpm(counts)
  assert_eq(log_cpm.length(), 2)
  assert_eq(log_cpm[0].length(), 3)
  assert_true(log_cpm[0][0] > 0.0)
}

///|
test "scater_find_hvgs" {
  let counts : Array[Array[Int]] = Array::make(10, Array::make(20, 0))
  let log_expr = @src.scater_log_normalize_cpm(counts)
  let gene_names : Array[String] = Array::make(20, "")
  let mut i = 0
  while i < 20 {
    gene_names[i] = "Gene" + i.to_string()
    i = i + 1
  }
  let hvgs = @src.scater_find_hvgs(log_expr, gene_names, 10)
  assert_eq(hvgs.length(), 10)
}

///|
test "scater_run_pca" {
  let counts : Array[Array[Int]] = Array::make(5, Array::make(10, 0))
  let mut i = 0
  while i < 5 {
    let mut j = 0
    while j < 10 {
      counts[i][j] = (i * 10 + j)
      j = j + 1
    }
    i = i + 1
  }
  let log_expr = @src.scater_log_normalize_cpm(counts)
  let pca = @src.scater_run_pca(log_expr, 3)
  assert_eq(pca.length(), 5)
  assert_eq(pca[0].length(), 3)
}

///|
test "scater_qc_report" {
  let counts = [[10, 20, 30], [5, 15, 25], [8, 18, 28]]
  let cell_ids = ["cell1", "cell2", "cell3"]
  let gene_names = ["Gene1", "Gene2", "mt-Atp6"]
  let cell_stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  let feature_stats = @src.scater_calculate_feature_qc(counts, gene_names)
  let params = @src.QCFilterParams::new()
  let report = @src.ScaterQCReport::new(cell_stats, feature_stats, params)
  assert_eq(report.cell_stats.cell_id.length(), 3)
  assert_eq(report.feature_stats.gene_id.length(), 3)
  assert_true(report.median_total_counts > 0.0)
}

///|
test "scater_filter_cells_by_name" {
  let counts = [[1, 2], [3, 4], [5, 6]]
  let cell_ids = ["cell1", "cell2", "cell3"]
  let gene_names = ["Gene1", "Gene2"]
  let stats = @src.scater_calculate_qc_metrics(counts, cell_ids, gene_names)
  let keep = @src.scater_filter_cells_by_name(stats, ["cell1", "cell3"])
  assert_eq(keep.length(), 3)
  assert_true(keep[0])
  assert_true(!keep[1])
  assert_true(keep[2])
}
