///|
/// DESeq2 tests - Differential Expression Analysis tests.

///|
test "test_sample_deseq_dataset" {
  let dds = @src.sample_deseq_dataset()
  
  assert_eq(dds.counts.length(), 20)
  assert_eq(dds.col_names.length(), 6)
  assert_eq(dds.row_names.length(), 20)
  assert_eq(dds.design.length(), 6)
  assert_eq(dds.design[0].length(), 2)
  assert_eq(dds.size_factors.length(), 6)
  assert_eq(dds.dispersions.length(), 20)
}

///|
test "test_estimate_size_factors" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  
  assert_eq(dds_sf.size_factors.length(), 6)
  
  let mut sum_sf = 0.0
  for sf in dds_sf.size_factors {
    assert_true(sf > 0.0)
    sum_sf = sum_sf + sf
  }
  
  let mean_sf = sum_sf / 6.0
  assert_true(mean_sf > 0.9 && mean_sf < 1.1)
}

///|
test "test_normalize_counts" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  let normalized = @src.normalize_counts(dds_sf)
  
  assert_eq(normalized.length(), 20)
  assert_eq(normalized[0].length(), 6)
}

///|
test "test_log2_cpm" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  let log2_cpm = @src.log2_cpm(dds_sf)
  
  assert_eq(log2_cpm.length(), 20)
  assert_eq(log2_cpm[0].length(), 6)
}

///|
test "test_estimate_dispersions" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  let dds_disp = @src.estimate_dispersions(dds_sf)
  
  assert_eq(dds_disp.dispersions.length(), 20)
  
  for disp in dds_disp.dispersions {
    assert_true(disp >= 0.001)
  }
}

///|
test "test_deseq" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  let dds_disp = @src.estimate_dispersions(dds_sf)
  let res = @src.deseq(dds_disp)
  
  assert_eq(res.row_names.length(), 20)
  assert_eq(res.base_mean.length(), 20)
  assert_eq(res.log2_fold_change.length(), 20)
  assert_eq(res.lfc_se.length(), 20)
  assert_eq(res.stat.length(), 20)
  assert_eq(res.p_value.length(), 20)
  assert_eq(res.padj.length(), 20)
}

///|
test "test_results" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  let dds_disp = @src.estimate_dispersions(dds_sf)
  let res = @src.results(dds_disp)
  
  assert_eq(res.row_names.length(), 20)
  assert_eq(res.base_mean.length(), 20)
  assert_eq(res.log2_fold_change.length(), 20)
  assert_eq(res.lfc_se.length(), 20)
  assert_eq(res.stat.length(), 20)
  assert_eq(res.p_value.length(), 20)
  assert_eq(res.padj.length(), 20)
  
  for p in res.p_value {
    assert_true(p >= 0.0 && p <= 1.0)
  }
  
  for p in res.padj {
    assert_true(p >= 0.0 && p <= 1.0)
  }
}

///|
test "test_lfc_shrink" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  let dds_disp = @src.estimate_dispersions(dds_sf)
  let res = @src.results(dds_disp)
  let res_shrunk = @src.lfc_shrink(dds_disp, res)
  
  assert_eq(res_shrunk.row_names.length(), 20)
  assert_eq(res_shrunk.log2_fold_change.length(), 20)
  
  for i = 0; i < 20; i = i + 1 {
    assert_true(res_shrunk.lfc_se[i] <= res.lfc_se[i])
  }
}

///|
test "test_significant_genes" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  let dds_disp = @src.estimate_dispersions(dds_sf)
  let res = @src.results(dds_disp)
  
  let sig_genes = @src.significant_genes(res)
  assert_true(sig_genes.length() >= 0 && sig_genes.length() <= 20)
  
  let sig_genes_lfc = @src.significant_genes(res, alpha=0.05, lfc_threshold=0.5)
  assert_true(sig_genes_lfc.length() >= 0 && sig_genes_lfc.length() <= sig_genes.length())
}

///|
test "test_top_genes" {
  let dds = @src.sample_deseq_dataset()
  let dds_sf = @src.estimate_size_factors(dds)
  let dds_disp = @src.estimate_dispersions(dds_sf)
  let res = @src.results(dds_disp)
  
  let top = @src.top_genes(res, n=5)
  assert_eq(top.length(), 5)
  
  for i = 0; i < 4; i = i + 1 {
    assert_true(top[i].2 <= top[i + 1].2)
  }
}