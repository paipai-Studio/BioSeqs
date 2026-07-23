///|
/// DSS module tests - Differential Expression and Methylation Analysis

///|
test "dss_create_example_counts" {
  let data = @src.dss_create_example_counts()
  assert_eq(data.dss_n_genes(), 5)
  assert_eq(data.dss_n_samples(), 6)
}

///|
test "dss_create_example_methylation" {
  let data = @src.dss_create_example_methylation()
  assert_eq(data.dss_n_cpgs(), 5)
  assert_eq(data.dss_n_samples(), 6)
}

///|
test "dss_estimate_dispersion" {
  let data = @src.dss_create_example_counts()
  let disp = @src.dss_estimate_dispersion(data)
  assert_eq(disp.dispersions.length(), 5)
  assert_eq(disp.shrunken_dispersions.length(), 5)
  assert_true(disp.prior_df > 0.0)
  assert_true(disp.prior_var > 0.0)
  for d in disp.dispersions {
    assert_true(d >= 0.0)
  }
  for sd in disp.shrunken_dispersions {
    assert_true(sd >= 0.0)
  }
}

///|
test "dss_deseq" {
  let data = @src.dss_create_example_counts()
  let disp = @src.dss_estimate_dispersion(data)
  let res = @src.dss_deseq(data, disp)
  assert_eq(res.dss_n_genes(), 5)
  assert_eq(res.gene_ids.length(), 5)
  assert_eq(res.mean1.length(), 5)
  assert_eq(res.mean2.length(), 5)
  assert_eq(res.dispersion.length(), 5)
  assert_eq(res.wald_stat.length(), 5)
  assert_eq(res.p_values.length(), 5)
  assert_eq(res.adj_p_values.length(), 5)
  assert_eq(res.significant.length(), 5)
  for p in res.p_values {
    assert_true(!p.is_nan())
  }
  for p in res.adj_p_values {
    assert_true(!p.is_nan())
  }
}

///|
test "dss_dml" {
  let data = @src.dss_create_example_methylation()
  let res = @src.dss_dml(data)
  assert_eq(res.dss_n_sites(), 5)
  assert_eq(res.cpg_ids.length(), 5)
  assert_eq(res.meth_mean1.length(), 5)
  assert_eq(res.meth_mean2.length(), 5)
  assert_eq(res.delta_meth.length(), 5)
  assert_eq(res.stats.length(), 5)
  assert_eq(res.p_values.length(), 5)
  assert_eq(res.adj_p_values.length(), 5)
  assert_eq(res.significant.length(), 5)
  for p in res.p_values {
    assert_true(!p.is_nan())
  }
  for p in res.adj_p_values {
    assert_true(!p.is_nan())
  }
}

///|
test "dss_dmr" {
  let data = @src.dss_create_example_methylation()
  let dm = @src.dss_dml(data)
  let dmr = @src.dss_dmr(data, dm, 1, 0.1)
  assert_true(dmr.dss_n_dmrs() >= 0)
  let n = dmr.dss_n_dmrs()
  assert_eq(dmr.dmr_ids.length(), n)
  assert_eq(dmr.chrom.length(), n)
  assert_eq(dmr.start.length(), n)
  assert_eq(dmr.end.length(), n)
  assert_eq(dmr.n_cpgs.length(), n)
  assert_eq(dmr.mean_delta.length(), n)
  assert_eq(dmr.stats.length(), n)
  assert_eq(dmr.p_values.length(), n)
}

///|
test "dss_count_data_methods" {
  let data = @src.dss_create_example_counts()
  assert_eq(data.dss_n_genes(), 5)
  assert_eq(data.dss_n_samples(), 6)

  let norm = data.dss_normalize()
  assert_eq(norm.length(), 5)
  assert_eq(norm[0].length(), 6)
  for row in norm {
    for val in row {
      assert_true(!val.is_nan())
    }
  }

  let filtered = data.dss_filter_low_counts(100.0, 2)
  assert_true(filtered.dss_n_genes() <= 5)
  let filtered2 = data.dss_filter_low_counts(1000.0, 10)
  assert_eq(filtered2.dss_n_genes(), 0)
}

///|
test "dss_de_result_methods" {
  let data = @src.dss_create_example_counts()
  let disp = @src.dss_estimate_dispersion(data)
  let res = @src.dss_deseq(data, disp)
  assert_eq(res.dss_n_genes(), 5)

  let n_sig = res.dss_n_significant()
  assert_true(n_sig >= 0 && n_sig <= 5)

  for i = 0; i < res.dss_n_genes(); i = i + 1 {
    let gene_id = res.dss_get_gene_id(i)
    assert_true(gene_id.length() > 0)
    let pval = res.dss_get_pvalue(i)
    assert_true(!pval.is_nan())
    let adj_pval = res.dss_get_adj_pvalue(i)
    assert_true(!adj_pval.is_nan())
    let logfc = res.dss_get_logfc(i)
    assert_true(!logfc.is_nan())
  }
}

///|
test "dss_bh_correction" {
  let data = @src.dss_create_example_counts()
  let disp = @src.dss_estimate_dispersion(data)
  let res = @src.dss_deseq(data, disp)

  for i = 0; i < res.dss_n_genes(); i = i + 1 {
    let pval = res.dss_get_pvalue(i)
    let adj_pval = res.dss_get_adj_pvalue(i)
    assert_true(!pval.is_nan())
    assert_true(!adj_pval.is_nan())
  }
}

///|
test "dss_dm_result_methods" {
  let data = @src.dss_create_example_methylation()
  let res = @src.dss_dml(data)
  assert_eq(res.dss_n_sites(), 5)

  let n_sig = res.dss_n_significant()
  assert_true(n_sig >= 0 && n_sig <= 5)

  for i = 0; i < res.dss_n_sites(); i = i + 1 {
    let delta = res.dss_get_delta_meth(i)
    assert_true(!delta.is_nan())
  }
}

///|
test "dss_dmr_result_methods" {
  let data = @src.dss_create_example_methylation()
  let dm = @src.dss_dml(data)
  let dmr = @src.dss_dmr(data, dm, 1, 0.05)
  let n = dmr.dss_n_dmrs()
  assert_true(n >= 0)
  for i = 0; i < n; i = i + 1 {
    let (id, chrom, start, end, n_cpgs, mean_delta) = dmr.dss_get_dmr(i)
    assert_true(id.length() > 0)
    assert_true(chrom.length() > 0)
    assert_true(end >= start)
    assert_true(n_cpgs >= 1)
    assert_true(!mean_delta.is_nan())
  }
}

///|
test "dss_disp_result_methods" {
  let data = @src.dss_create_example_counts()
  let disp = @src.dss_estimate_dispersion(data)
  let dispersions = disp.dss_get_dispersions()
  let shrunken = disp.dss_get_shrunken_dispersions()
  assert_eq(dispersions.length(), 5)
  assert_eq(shrunken.length(), 5)
  for d in dispersions {
    assert_true(d >= 0.0)
  }
  for sd in shrunken {
    assert_true(sd >= 0.0)
  }
}