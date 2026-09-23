///|
test "goseq_nullp_basic" {
  // Test basic PWF calculation
  let lengths = [1000.0, 2000.0, 3000.0, 4000.0, 5000.0]
  let de = [1, 0, 1, 0, 0]
  let pwf = @src.goseq_nullp(lengths, de)
  assert_eq(pwf.length(), 5)
  let mut i = 0
  while i < pwf.length() {
    assert_true(pwf[i] >= 0.0 && pwf[i] <= 1.0)
    i = i + 1
  }
}

///|
test "goseq_nullp_empty" {
  let pwf = @src.goseq_nullp([], [])
  assert_eq(pwf.length(), 0)
}

///|
test "goseq_nullp_all_de" {
  let lengths = [1000.0, 2000.0, 3000.0]
  let de = [1, 1, 1]
  let pwf = @src.goseq_nullp(lengths, de)
  assert_eq(pwf.length(), 3)
  let mut i = 0
  while i < pwf.length() {
    assert_true(pwf[i] >= 0.0 && pwf[i] <= 1.0)
    i = i + 1
  }
}

///|
test "goseq_nullp_none_de" {
  let lengths = [1000.0, 2000.0, 3000.0]
  let de = [0, 0, 0]
  let pwf = @src.goseq_nullp(lengths, de)
  assert_eq(pwf.length(), 3)
  let mut i = 0
  while i < pwf.length() {
    assert_true(pwf[i] >= 0.0 && pwf[i] <= 1.0)
    i = i + 1
  }
}

///|
test "goseq_wallenius_pvalue_basic" {
  // Wallenius with odds=1 should match hypergeometric
  let wallenius_p = @src.wallenius_pvalue(3, 20, 10, 100, 1.0)
  let hypergeo_p = @src.hypergeo_pvalue(3, 20, 10, 100)
  assert_true((wallenius_p - hypergeo_p).abs() < 1.0e-10)
}

///|
test "goseq_wallenius_pvalue_biased" {
  // Wallenius with odds > 1 (bias toward DE genes in this category)
  let p = @src.wallenius_pvalue(3, 20, 10, 100, 2.0)
  assert_true(p >= 0.0 && p <= 1.0)
}

///|
test "goseq_wallenius_pvalue_no_bias_low_odds" {
  // Wallenius with odds < 1 (bias against DE genes in this category)
  let p = @src.wallenius_pvalue(3, 20, 10, 100, 0.5)
  assert_true(p >= 0.0 && p <= 1.0)
}

///|
test "goseq_main" {
  let genes = [
    "GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH",
  ]
  let lengths = [1000.0, 2000.0, 1500.0, 3000.0, 2500.0, 1800.0, 4000.0, 3500.0]
  let de = [1, 0, 1, 0, 0, 1, 0, 0]

  let go_sets = [
    @src.GoseqGeneSet::new("GO:0006915", "apoptosis", "BP", [
      "GeneA", "GeneB", "GeneC",
    ]),
    @src.GoseqGeneSet::new("GO:0007049", "cell cycle", "BP", [
      "GeneD", "GeneE", "GeneF",
    ]),
    @src.GoseqGeneSet::new("GO:0008150", "biological process", "BP", [
      "GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH",
    ]),
  ]

  let results = @src.goseq(genes, lengths, de, go_sets)
  assert_true(results.length() > 0)
  let mut i = 0
  while i < results.length() {
    assert_true(results[i].p_value >= 0.0 && results[i].p_value <= 1.0)
    assert_true(
      results[i].adjusted_p_value >= 0.0 && results[i].adjusted_p_value <= 1.0,
    )
    assert_true(results[i].overlap >= 0)
    i = i + 1
  }
  // Results should be sorted by p-value ascending
  i = 0
  while i < results.length() - 1 {
    assert_true(results[i].p_value <= results[i + 1].p_value + 1.0e-10)
    i = i + 1
  }
}

///|
test "goseq_empty_genes" {
  let results = @src.goseq([], [], [], [])
  assert_eq(results.length(), 0)
}

///|
test "goseq_no_go_matches" {
  let genes = ["GeneA", "GeneB", "GeneC"]
  let lengths = [1000.0, 2000.0, 1500.0]
  let de = [1, 0, 1]
  let go_sets = [
    @src.GoseqGeneSet::new("GO:0006915", "apoptosis", "BP", ["GeneX", "GeneY"]),
  ]
  let results = @src.goseq(genes, lengths, de, go_sets)
  assert_eq(results.length(), 0)
}

///|
test "goseq_filter_results" {
  let result1 = @src.GoseqResult::new(
    "GO:0006915",
    "apoptosis",
    "BP",
    0.001,
    0.01,
    3,
    3,
    10,
    20,
    ["GeneA", "GeneB", "GeneC"],
  )
  let result2 = @src.GoseqResult::new(
    "GO:0007049",
    "cell cycle",
    "BP",
    0.05,
    0.15,
    1,
    1,
    5,
    20,
    ["GeneD"],
  )
  let result3 = @src.GoseqResult::new(
    "GO:0008150",
    "biological process",
    "BP",
    0.01,
    0.03,
    2,
    2,
    8,
    20,
    ["GeneA", "GeneD"],
  )

  let results = [result1, result2, result3]
  let filtered = @src.goseq_filter_results(results, 0.05, 2)
  assert_eq(filtered.length(), 2)
  assert_eq(filtered[0].go_id, "GO:0006915")
  assert_eq(filtered[1].go_id, "GO:0008150")
}

///|
test "goseq_result_new" {
  let result = @src.GoseqResult::new(
    "GO:0006915",
    "apoptosis",
    "BP",
    0.001,
    0.01,
    3,
    3,
    10,
    20,
    ["GeneA", "GeneB", "GeneC"],
  )
  assert_eq(result.go_id, "GO:0006915")
  assert_eq(result.term_name, "apoptosis")
  assert_eq(result.ontology, "BP")
  assert_true((result.p_value - 0.001).abs() < 1.0e-10)
  assert_true((result.adjusted_p_value - 0.01).abs() < 1.0e-10)
  assert_eq(result.overlap, 3)
  assert_eq(result.num_de_in_term, 3)
  assert_eq(result.num_in_term, 10)
  assert_eq(result.num_de, 20)
  assert_eq(result.term_genes.length(), 3)
}

///|
test "goseq_set_new" {
  let go_set = @src.GoseqGeneSet::new("GO:0006915", "apoptosis", "BP", [
    "GeneA", "GeneB", "GeneC",
  ])
  assert_eq(go_set.go_id, "GO:0006915")
  assert_eq(go_set.term_name, "apoptosis")
  assert_eq(go_set.ontology, "BP")
  assert_eq(go_set.gene_ids.length(), 3)
}

///|
test "goseq_longer_genes_bias_detection" {
  // Test that goseq can detect length bias:
  // Longer genes are more likely to be DE
  let genes = Array(capacity=0)
  let lengths = Array(capacity=0)
  let de = Array(capacity=0)
  let mut i = 0
  while i < 50 {
    // Short genes: 500-1500 bp, mostly non-DE
    let len = 500.0 + (i * 20).to_double()
    genes.push("ShortGene" + i.to_string())
    lengths.push(len)
    de.push(if i < 10 { 1 } else { 0 })
    i = i + 1
  }
  i = 0
  while i < 50 {
    // Long genes: 3000-5000 bp, mostly DE
    let len = 3000.0 + (i * 40).to_double()
    genes.push("LongGene" + i.to_string())
    lengths.push(len)
    de.push(if i < 40 { 1 } else { 0 })
    i = i + 1
  }

  let pwf = @src.goseq_nullp(lengths, de)
  assert_eq(pwf.length(), 100)
  let mut i = 0
  while i < pwf.length() {
    assert_true(pwf[i] >= 0.0 && pwf[i] <= 1.0)
    i = i + 1
  }
  // Long genes should have higher PWF values on average
  let mut short_pwf_sum = 0.0
  i = 0
  while i < 50 {
    short_pwf_sum = short_pwf_sum + pwf[i]
    i = i + 1
  }
  let mut long_pwf_sum = 0.0
  i = 50
  while i < 100 {
    long_pwf_sum = long_pwf_sum + pwf[i]
    i = i + 1
  }
  let short_mean = short_pwf_sum / 50.0
  let long_mean = long_pwf_sum / 50.0
  // Long genes should have higher PWF (more DE)
  assert_true(long_mean > short_mean)
}
