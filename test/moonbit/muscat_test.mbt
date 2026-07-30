///|
/// Test file for muscat module.

test "muscat_single_cell_creation" {
  let cell = @src.SingleCell::new("cell1", "sample1", "cluster1", "ctrl")
  assert_eq(cell.cell_id, "cell1")
  assert_eq(cell.sample_id, "sample1")
  assert_eq(cell.cluster_id, "cluster1")
  assert_eq(cell.group_id, "ctrl")
}

test "muscat_single_cell_set_get_count" {
  let cell = @src.SingleCell::new("cell1", "sample1", "cluster1", "ctrl")
  cell.set_count("Gene1", 10.0)
  cell.set_count("Gene2", 0.0)
  assert_eq(cell.get_count("Gene1"), 10.0)
  assert_eq(cell.get_count("Gene2"), 0.0)
  assert_eq(cell.get_count("Gene3"), 0.0) // Non-existent
}

test "muscat_single_cell_total_counts" {
  let cell = @src.SingleCell::new("cell1", "sample1", "cluster1", "ctrl")
  cell.set_count("Gene1", 10.0)
  cell.set_count("Gene2", 20.0)
  cell.set_count("Gene3", 5.0)
  assert_eq(cell.total_counts(), 35.0)
}

test "muscat_single_cell_n_expressed" {
  let cell = @src.SingleCell::new("cell1", "sample1", "cluster1", "ctrl")
  cell.set_count("Gene1", 10.0)
  cell.set_count("Gene2", 0.0)
  cell.set_count("Gene3", 5.0)
  assert_eq(cell.n_expressed(), 2)
}

test "muscat_aggregation_sum" {
  let cells : Array[@src.SingleCell] = Array::new()
  let c1 = @src.SingleCell::new("c1", "S1", "C1", "ctrl")
  c1.set_count("Gene1", 10.0)
  c1.set_count("Gene2", 5.0)
  cells.push(c1)
  let c2 = @src.SingleCell::new("c2", "S1", "C1", "ctrl")
  c2.set_count("Gene1", 20.0)
  c2.set_count("Gene2", 3.0)
  cells.push(c2)

  let pb = @src.aggregate_cells(cells, method=@src.aggregation_sum())
  assert_eq(pb.length(), 1)
  assert_eq(pb[0].get_count("Gene1"), 30.0)
  assert_eq(pb[0].get_count("Gene2"), 8.0)
  assert_eq(pb[0].n_cells, 2)
}

test "muscat_aggregation_mean" {
  let cells : Array[@src.SingleCell] = Array::new()
  let c1 = @src.SingleCell::new("c1", "S1", "C1", "ctrl")
  c1.set_count("Gene1", 10.0)
  cells.push(c1)
  let c2 = @src.SingleCell::new("c2", "S1", "C1", "ctrl")
  c2.set_count("Gene1", 20.0)
  cells.push(c2)

  let pb = @src.aggregate_cells(cells, method=@src.aggregation_mean())
  assert_eq(pb.length(), 1)
  assert_eq(pb[0].get_count("Gene1"), 15.0)
}

test "muscat_aggregation_multiple_samples" {
  let cells : Array[@src.SingleCell] = Array::new()
  // Sample S1, Cluster C1
  let c1 = @src.SingleCell::new("c1", "S1", "C1", "ctrl")
  c1.set_count("Gene1", 10.0)
  cells.push(c1)
  // Sample S1, Cluster C2
  let c2 = @src.SingleCell::new("c2", "S1", "C2", "ctrl")
  c2.set_count("Gene1", 20.0)
  cells.push(c2)
  // Sample S2, Cluster C1
  let c3 = @src.SingleCell::new("c3", "S2", "C1", "ctrl")
  c3.set_count("Gene1", 15.0)
  cells.push(c3)

  let pb = @src.aggregate_cells(cells, method=@src.aggregation_sum())
  assert_eq(pb.length(), 3) // 3 (sample, cluster) combinations
}

test "muscat_pseudobulk_total_counts" {
  let cells : Array[@src.SingleCell] = Array::new()
  let c1 = @src.SingleCell::new("c1", "S1", "C1", "ctrl")
  c1.set_count("Gene1", 10.0)
  c1.set_count("Gene2", 20.0)
  cells.push(c1)
  let c2 = @src.SingleCell::new("c2", "S1", "C1", "ctrl")
  c2.set_count("Gene1", 5.0)
  cells.push(c2)

  let pb = @src.aggregate_cells(cells)
  assert_eq(pb[0].total_counts(), 35.0)
}

test "muscat_pseudobulk_n_expressed" {
  let cells : Array[@src.SingleCell] = Array::new()
  let c1 = @src.SingleCell::new("c1", "S1", "C1", "ctrl")
  c1.set_count("Gene1", 10.0)
  c1.set_count("Gene2", 0.0)
  cells.push(c1)
  let c2 = @src.SingleCell::new("c2", "S1", "C1", "ctrl")
  c2.set_count("Gene1", 0.0)
  c2.set_count("Gene3", 5.0)
  cells.push(c2)

  let pb = @src.aggregate_cells(cells, method=@src.aggregation_sum())
  assert_eq(pb[0].n_expressed(), 2) // Gene1 and Gene3 > 0
}

test "muscat_ds_analysis_basic" {
  let cells = @src.muscat_sample_data()
  let pb = @src.aggregate_cells(cells)
  let results = @src.run_ds_analysis(pb)
  assert_true(results.get_n_results() > 0)
  assert_eq(results.n_clusters, 2)
}

test "muscat_ds_analysis_significant" {
  let cells = @src.muscat_sample_data()
  let pb = @src.aggregate_cells(cells)
  let results = @src.run_ds_analysis(pb, fdr_threshold=0.5, log2fc_threshold=0.5)
  let sig = results.get_significant()
  assert_true(sig.length() > 0)
}

test "muscat_ds_analysis_cluster_filter" {
  let cells = @src.muscat_sample_data()
  let pb = @src.aggregate_cells(cells)
  let results = @src.run_ds_analysis(pb)
  let c1_results = results.get_cluster_results("C1")
  assert_true(c1_results.length() > 0)
}

test "muscat_ds_analysis_top_genes" {
  let cells = @src.muscat_sample_data()
  let pb = @src.aggregate_cells(cells)
  let results = @src.run_ds_analysis(pb)
  let top = results.get_top_genes(3)
  assert_eq(top.length(), 3)
}

test "muscat_qc_computation" {
  let cells = @src.muscat_sample_data()
  let pb = @src.aggregate_cells(cells)
  let qc = @src.compute_qc(pb)
  assert_eq(qc.length(), pb.length())
}

test "muscat_summary" {
  let cells = @src.muscat_sample_data()
  let pb = @src.aggregate_cells(cells)
  let results = @src.run_ds_analysis(pb)
  let s = results.summary()
  assert_true(s.contains("muscat"))
}

test "muscat_aggregation_method_to_string" {
  assert_eq(@src.aggregation_sum().to_string(), "sum")
  assert_eq(@src.aggregation_mean().to_string(), "mean")
  assert_eq(@src.aggregation_median().to_string(), "median")
}

test "muscat_ds_method_to_string" {
  assert_eq(@src.ds_method_edger().to_string(), "edgeR")
  assert_eq(@src.ds_method_deseq2().to_string(), "DESeq2")
  assert_eq(@src.ds_method_limma().to_string(), "limma")
}
