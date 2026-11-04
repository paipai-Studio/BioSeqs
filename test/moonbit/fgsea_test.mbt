///|
test "fgsea_gene_set_new" {
  let genes = ["GeneA", "GeneB", "GeneC"]
  let gs = @src.GeneSet::new("test_set", genes)

  assert_eq(gs.name, "test_set")
  assert_eq(gs.genes.length(), 3)
}

///|
test "fgsea_rank_genes" {
  let genes = ["GeneA", "GeneB", "GeneC"]
  let stats = [1.0, 2.0, 0.5]
  let ranked = @src.rank_genes(genes, stats)

  assert_eq(ranked.length(), 3)
  assert_eq(ranked[0].0, "GeneB")
  assert_eq(ranked[1].0, "GeneA")
  assert_eq(ranked[2].0, "GeneC")
}

///|
test "fgsea_result_new" {
  let result = @src.FgseaResult::new("pathway1", 0.01, 0.02, 0.5, 1.5, 5, 10, [
    "GeneA", "GeneB",
  ])

  assert_eq(result.pathway, "pathway1")
  assert_eq(result.pval, 0.01)
  assert_eq(result.padj, 0.02)
  assert_eq(result.es, 0.5)
  assert_eq(result.nes, 1.5)
  assert_eq(result.nMoreExtreme, 5)
  assert_eq(result.size, 10)
  assert_eq(result.leadingEdge.length(), 2)
}

///|
test "fgsea_sample_data" {
  let gene_sets = @src.sample_gene_sets()
  let stats = @src.sample_gene_stats()

  assert_eq(gene_sets.length(), 5)
  assert_true(stats.length() > 0)
  assert_eq(stats[0].1 >= stats[stats.length() - 1].1, true)
}

///|
test "fgsea_run_fgsea" {
  let gene_sets = @src.sample_gene_sets()
  let stats = @src.sample_gene_stats()

  let results = @src.fgsea(gene_sets, stats, nperm=50)

  assert_eq(results.length(), 5)
  for r in results {
    assert_true(r.pval >= 0.0 && r.pval <= 1.0)
    assert_true(r.padj >= 0.0 && r.padj <= 1.0)
    assert_true(r.size > 0)
  }
}

///|
test "fgsea_filter" {
  let gene_sets = @src.sample_gene_sets()
  let stats = @src.sample_gene_stats()
  let results = @src.fgsea(gene_sets, stats, nperm=50)

  let filtered = @src.fgsea_filter(results, padj_cutoff=0.5)

  assert_true(filtered.length() <= results.length())
}

///|
test "fgsea_top_up" {
  let gene_sets = @src.sample_gene_sets()
  let stats = @src.sample_gene_stats()
  let results = @src.fgsea(gene_sets, stats, nperm=50)

  let top_up = @src.fgsea_top_up(results, n=3)

  assert_true(top_up.length() <= 3)
  for r in top_up {
    assert_true(r.nes > 0.0)
  }
}

///|
test "fgsea_top_down" {
  let gene_sets = @src.sample_gene_sets()
  let stats = @src.sample_gene_stats()
  let results = @src.fgsea(gene_sets, stats, nperm=50)

  let top_down = @src.fgsea_top_down(results, n=3)

  assert_true(top_down.length() <= 3)
  for r in top_down {
    assert_true(r.nes < 0.0)
  }
}
