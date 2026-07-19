///|
/// DESeq2 tests - Differential Expression Analysis tests.

import "mooncake/test"

///|
fn test_sample_deseq_dataset() {
  let dds = @bio.sample_deseq_dataset()
  
  assert_eq(dds.counts.length(), 20)
  assert_eq(dds.col_names.length(), 6)
  assert_eq(dds.row_names.length(), 20)
  assert_eq(dds.design.length(), 6)
  assert_eq(dds.design[0].length(), 2)
  assert_eq(dds.size_factors.length(), 6)
  assert_eq(dds.dispersions.length(), 20)
  
  println("test_sample_deseq_dataset passed")
}

///|
fn test_estimate_size_factors() {
  let dds = @bio.sample_deseq_dataset()
  let dds_sf = @bio.estimate_size_factors(dds)
  
  assert_eq(dds_sf.size_factors.length(), 6)
  
  let mut sum_sf = 0.0
  for sf in dds_sf.size_factors {
    assert_true(sf > 0.0)
    sum_sf = sum_sf + sf
  }
  
  let mean_sf = sum_sf / 6.0
  assert_true(mean_sf > 0.9 && mean_sf < 1.1)
  
  println("test_estimate_size_factors passed")
}

///|
fn test_normalize_counts() {
  let dds = @bio.sample_deseq_dataset()
  let dds_sf = @bio.estimate_size_factors(dds)
  let normalized = @bio.normalize_counts(dds_sf)
  
  assert_eq(normalized.length(), 20)
  assert_eq(normalized[0].length(), 6)
  
  println("test_normalize_counts passed")
}

///|
fn test_log2_cpm() {
  let dds = @bio.sample_deseq_dataset()
  let dds_sf = @bio.estimate_size_factors(dds)
  let log2_cpm = @bio.log2_cpm(dds_sf)
  
  assert_eq(log2_cpm.length(), 20)
  assert_eq(log2_cpm[0].length(), 6)
  
  println("test_log2_cpm passed")
}

///|
fn test_estimate_dispersions() {
  let dds = @bio.sample_deseq_dataset()
  let dds_sf = @bio.estimate_size_factors(dds)
  let dds_disp = @bio.estimate_dispersions(dds_sf)
  
  assert_eq(dds_disp.dispersions.length(), 20)
  
  for disp in dds_disp.dispersions {
    assert_true(disp >= 0.001)
  }
  
  println("test_estimate_dispersions passed")
}

///|
fn test_deseq() {
  let dds = @bio.sample_deseq_dataset()
  let dds_sf = @bio.estimate_size_factors(dds)
  let dds_disp = @bio.estimate_dispersions(dds_sf)
  let res = @bio.deseq(dds_disp)
  
  assert_eq(res.row_names.length(), 20)
  assert_eq(res.base_mean.length(), 20)
  assert_eq(res.log2_fold_change.length(), 20)
  assert_eq(res.lfc_se.length(), 20)
  assert_eq(res.stat.length(), 20)
  assert_eq(res.p_value.length(), 20)
  assert_eq(res.padj.length(), 20)
  
  println("test_deseq passed")
}

///|
fn test_results() {
  let dds = @bio.sample_deseq_dataset()
  let dds_sf = @bio.estimate_size_factors(dds)
  let dds_disp = @bio.estimate_dispersions(dds_sf)
  let res = @bio.results(dds_disp)
  
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
  
  println("test_results passed")
}

///|
fn test_lfc_shrink() {
  let dds = @bio.sample_deseq_dataset()
  let dds_sf = @bio.estimate_size_factors(dds)
  let dds_disp = @bio.estimate_dispersions(dds_sf)
  let res = @bio.results(dds_disp)
  let res_shrunk = @bio.lfc_shrink(dds_disp, res)
  
  assert_eq(res_shrunk.row_names.length(), 20)
  assert_eq(res_shrunk.log2_fold_change.length(), 20)
  
  for i = 0; i < 20; i = i + 1 {
    assert_true(res_shrunk.lfc_se[i] <= res.lfc_se[i])
  }
  
  println("test_lfc_shrink passed")
}

///|
fn test_significant_genes() {
  let dds = @bio.sample_deseq_dataset()
  let dds_sf = @bio.estimate_size_factors(dds)
  let dds_disp = @bio.estimate_dispersions(dds_sf)
  let res = @bio.results(dds_disp)
  
  let sig_genes = @bio.significant_genes(res)
  assert_true(sig_genes.length() >= 0 && sig_genes.length() <= 20)
  
  let sig_genes_lfc = @bio.significant_genes(res, alpha=0.05, lfc_threshold=0.5)
  assert_true(sig_genes_lfc.length() >= 0 && sig_genes_lfc.length() <= sig_genes.length())
  
  println("test_significant_genes passed")
}

///|
fn test_top_genes() {
  let dds = @bio.sample_deseq_dataset()
  let dds_sf = @bio.estimate_size_factors(dds)
  let dds_disp = @bio.estimate_dispersions(dds_sf)
  let res = @bio.results(dds_disp)
  
  let top = @bio.top_genes(res, n=5)
  assert_eq(top.length(), 5)
  
  for i = 0; i < 4; i = i + 1 {
    assert_true(top[i].2 <= top[i + 1].2)
  }
  
  println("test_top_genes passed")
}

///|
fn main() {
  test_sample_deseq_dataset()
  test_estimate_size_factors()
  test_normalize_counts()
  test_log2_cpm()
  test_estimate_dispersions()
  test_deseq()
  test_results()
  test_lfc_shrink()
  test_significant_genes()
  test_top_genes()
  
  println("\nAll DESeq2 tests passed!")
}
