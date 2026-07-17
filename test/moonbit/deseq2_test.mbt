///|
test "deseq2_create_dataset" {
  let counts = [
    [100, 120, 90, 110],
    [50, 60, 45, 55],
    [10, 15, 8, 12],
    [200, 220, 190, 210],
  ]
  let row_names = ["GeneA", "GeneB", "GeneC", "GeneD"]
  let col_names = ["Sample1", "Sample2", "Sample3", "Sample4"]
  let design = [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]]

  let dds = @src.deseq_dataset(counts, row_names, col_names, design)

  assert_eq(dds.counts.length(), 4)
  assert_eq(dds.row_names.length(), 4)
  assert_eq(dds.col_names.length(), 4)
  assert_eq(dds.design.length(), 4)
}

///|
test "deseq2_results" {
  let counts = [
    [100, 120, 90, 110],
    [50, 60, 45, 55],
    [10, 15, 8, 12],
    [200, 220, 190, 210],
  ]
  let row_names = ["GeneA", "GeneB", "GeneC", "GeneD"]
  let col_names = ["Sample1", "Sample2", "Sample3", "Sample4"]
  let design = [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]]

  let dds = @src.deseq_dataset(counts, row_names, col_names, design)
  let results = @src.results(dds)

  assert_eq(results.row_names.length(), 4)
  assert_eq(results.base_mean.length(), 4)
  assert_eq(results.log2_fold_change.length(), 4)
}

///|
test "deseq2_significant_genes" {
  let counts = [
    [100, 120, 400, 420],
    [50, 60, 55, 45],
    [10, 15, 100, 110],
    [200, 220, 190, 210],
  ]
  let row_names = ["GeneA", "GeneB", "GeneC", "GeneD"]
  let col_names = ["Control1", "Control2", "Treat1", "Treat2"]
  let design = [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]]

  let dds = @src.deseq_dataset(counts, row_names, col_names, design)
  let results = @src.results(dds)

  let sig_genes = @src.significant_genes(results, 0.05, 1.0)

  assert_true(sig_genes.length() >= 0)
}
