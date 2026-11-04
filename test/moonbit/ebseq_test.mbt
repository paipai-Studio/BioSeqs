///|
/// EBSeq tests - Empirical Bayesian differential expression analysis.

///|
test "ebseq_dataset_new_creates_default_size_factors_and_dispersions" {
  let counts = [[100, 120, 110, 200, 210, 220], [50, 60, 55, 50, 55, 60]]
  let gene_names = ["Gene1", "Gene2"]
  let sample_names = ["C1_R1", "C1_R2", "C1_R3", "C2_R1", "C2_R2", "C2_R3"]
  let conditions = [0, 0, 0, 1, 1, 1]
  let data = @src.EBSeqDataSet::new(
    counts, gene_names, sample_names, conditions,
  )

  // Check basic fields
  assert_eq(data.counts.length(), 2)
  assert_eq(data.gene_names.length(), 2)
  assert_eq(data.sample_names.length(), 6)
  assert_eq(data.conditions.length(), 6)

  // Default size_factors should all be 1.0
  assert_eq(data.size_factors.length(), 6)
  for sf in data.size_factors {
    assert_true(sf > 0.99 && sf < 1.01)
  }

  // Default dispersions should all be 0.1
  assert_eq(data.dispersions.length(), 2)
  for d in data.dispersions {
    assert_true(d > 0.09 && d < 0.11)
  }
}

///|
test "ebseq_estimate_size_factors_computes_reasonable_values" {
  let data = @src.ebseq_sample_dataset()
  let data_sf = @src.ebseq_estimate_size_factors(data)

  // Size factors should have same length as samples
  assert_eq(data_sf.size_factors.length(), 6)

  // All size factors should be positive
  for sf in data_sf.size_factors {
    assert_true(sf > 0.0)
  }

  // Geometric mean of size factors should be close to 1.0
  let mut sum_log = 0.0
  for sf in data_sf.size_factors {
    sum_log = sum_log + @math.ln(sf)
  }
  let gm = @math.exp(sum_log / 6.0)
  assert_true(gm > 0.99 && gm < 1.01)
}

///|
test "ebseq_estimate_dispersion_produces_positive_dispersions" {
  let data = @src.ebseq_sample_dataset()
  let data_sf = @src.ebseq_estimate_size_factors(data)
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)

  assert_eq(data_disp.dispersions.length(), 20)

  // All dispersions should be positive (minimum 0.001)
  for d in data_disp.dispersions {
    assert_true(d >= 0.001)
  }
}

///|
test "ebseq_test_returns_ppde_between_0_and_1" {
  let data = @src.ebseq_sample_dataset()
  let data_sf = @src.ebseq_estimate_size_factors(data)
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)
  let res = @src.ebseq_test(data_disp)

  // PPDE values should all be in [0, 1]
  assert_eq(res.ppde.length(), 20)
  for p in res.ppde {
    assert_true(p >= 0.0 && p <= 1.0)
  }
}

///|
test "ebseq_test_identifies_de_genes_from_sample_dataset" {
  let data = @src.ebseq_sample_dataset()
  let data_sf = @src.ebseq_estimate_size_factors(data)
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)
  let res = @src.ebseq_test(data_disp, threshold=0.5)

  // Sample dataset has DE genes (first 10), so at least some should be flagged
  let de_count = @src.ebseq_count_de(res)
  assert_true(de_count > 0)
  assert_true(de_count <= 20)
}

///|
test "ebseq_test_log2_fc_sign_for_up_down_regulated" {
  let data = @src.ebseq_sample_dataset()
  let data_sf = @src.ebseq_estimate_size_factors(data)
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)
  let res = @src.ebseq_test(data_disp)

  // In sample dataset:
  //   Genes 0-4: up-regulated in condition 2 (fold_change=3.0) => log2_fc > 0
  //   Genes 5-9: down-regulated in condition 2 (fold_change=0.3) => log2_fc < 0
  // Check at least the first gene of each group
  assert_true(res.log2_fc[0] > 0.0)
  assert_true(res.log2_fc[5] < 0.0)
}

///|
test "ebseq_count_de_works_correctly" {
  let data = @src.ebseq_sample_dataset()
  let data_sf = @src.ebseq_estimate_size_factors(data)
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)
  let res = @src.ebseq_test(data_disp, threshold=0.95)

  let de_count = @src.ebseq_count_de(res)
  // Count should match the number of true flags
  let mut manual_count = 0
  for flag in res.flag_de {
    if flag {
      manual_count = manual_count + 1
    }
  }
  assert_eq(de_count, manual_count)
}

///|
test "ebseq_get_de_genes_returns_correct_names" {
  let data = @src.ebseq_sample_dataset()
  let data_sf = @src.ebseq_estimate_size_factors(data)
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)
  let res = @src.ebseq_test(data_disp, threshold=0.5)

  let de_genes = @src.ebseq_get_de_genes(res)
  // Number of DE gene names should match count_de
  let de_count = @src.ebseq_count_de(res)
  assert_eq(de_genes.length(), de_count)

  // All returned names should exist in the original gene_names
  for gene in de_genes {
    let mut found = false
    for name in res.gene_names {
      if name == gene {
        found = true
      }
    }
    assert_true(found)
  }
}

///|
test "ebseq_sample_dataset_dimensions" {
  let data = @src.ebseq_sample_dataset()
  assert_eq(data.counts.length(), 20)
  assert_eq(data.counts[0].length(), 6)
  assert_eq(data.gene_names.length(), 20)
  assert_eq(data.sample_names.length(), 6)
  assert_eq(data.conditions.length(), 6)
}

///|
test "ebseq_result_has_consistent_lengths" {
  let data = @src.ebseq_sample_dataset()
  let data_sf = @src.ebseq_estimate_size_factors(data)
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)
  let res = @src.ebseq_test(data_disp)

  let n = res.gene_names.length()
  assert_eq(res.ppde.length(), n)
  assert_eq(res.log2_fc.length(), n)
  assert_eq(res.mean_cond1.length(), n)
  assert_eq(res.mean_cond2.length(), n)
  assert_eq(res.flag_de.length(), n)
}

///|
test "ebseq_top_genes_returns_sorted_by_ppde" {
  let data = @src.ebseq_sample_dataset()
  let data_sf = @src.ebseq_estimate_size_factors(data)
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)
  let res = @src.ebseq_test(data_disp)

  let top = @src.ebseq_top_genes(res, n=5)
  assert_eq(top.length(), 5)

  // PPDE values should be in descending order
  for i = 0; i < 4; i = i + 1 {
    let (_, ppde_i, _) = top[i]
    let (_, ppde_next, _) = top[i + 1]
    assert_true(ppde_i >= ppde_next)
  }
}

///|
test "ebseq_test_with_custom_two_condition_data" {
  // Create a simple 2-condition dataset with clear DE
  let counts = [
    [100, 110, 105, 300, 310, 290], // up in cond2
    [200, 210, 205, 200, 195, 210], // equal expression
    [300, 290, 310, 100, 110, 90], // down in cond2
  ]
  let gene_names = ["UpGene", "EqGene", "DownGene"]
  let sample_names = ["C1_R1", "C1_R2", "C1_R3", "C2_R1", "C2_R2", "C2_R3"]
  let conditions = [0, 0, 0, 1, 1, 1]

  let data = @src.EBSeqDataSet::new(
    counts, gene_names, sample_names, conditions,
  )
  let data_sf = @src.ebseq_estimate_size_factors(data)
  let data_disp = @src.ebseq_estimate_dispersion(data_sf)
  let res = @src.ebseq_test(data_disp)

  // UpGene should have positive log2_fc, DownGene negative
  assert_true(res.log2_fc[0] > 0.0)
  assert_true(res.log2_fc[2] < 0.0)

  // All PPDE should be valid
  for p in res.ppde {
    assert_true(p >= 0.0 && p <= 1.0)
  }
}
