///|
test "cluster_profiler_enrich_generic" {
  let gene_list = ["TP53", "BRCA1", "BRCA2"]
  let background_genes = ["TP53", "BRCA1", "BRCA2", "GAPDH", "ACTB"]

  let term_db = Map([], capacity=2)
  term_db["GO:0005515"] = ("protein binding", ["TP53", "BRCA1"])
  term_db["GO:0007283"] = ("spermatogenesis", ["BRCA2"])

  let result = @src.cp_enrich_generic(
    gene_list, term_db, background_genes, "Homo sapiens", "GO",
  )

  assert_true(result.terms.length() >= 1)
}

///|
test "cluster_profiler_enrich_go" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH"]

  let go_db = Map([], capacity=2)
  go_db["GO:0005515"] = ("protein binding", ["TP53", "BRCA1"], "MF")
  go_db["GO:0046982"] = ("protein heterodimerization", ["TP53"], "MF")

  let result = @src.cp_enrich_go(
    gene_list, go_db, background_genes, "Homo sapiens", "MF",
  )

  assert_true(result.terms.length() >= 1)
}

///|
test "cluster_profiler_enrich_kegg" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH"]

  let kegg_db = Map([], capacity=1)
  kegg_db["hsa04115"] = ("p53 signaling pathway", ["TP53", "BRCA1"])

  let result = @src.cp_enrich_kegg(
    gene_list, kegg_db, background_genes, "Homo sapiens",
  )

  assert_eq(result.keytype, "KEGG")
}

///|
test "cluster_profiler_filter_results" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH"]

  let term_db = Map([], capacity=1)
  term_db["GO:0005515"] = ("protein binding", ["TP53", "BRCA1"])

  let result = @src.cp_enrich_generic(
    gene_list, term_db, background_genes, "Homo sapiens", "GO",
  )
  let filtered = @src.cp_filter_results(result, 0.05)

  assert_eq(filtered.query_count, 2)
}

///|
test "cluster_profiler_get_top_terms" {
  let gene_list = ["TP53", "BRCA1", "BRCA2"]
  let background_genes = ["TP53", "BRCA1", "BRCA2", "GAPDH"]

  let term_db = Map([], capacity=3)
  term_db["GO:0005515"] = ("protein binding", ["TP53", "BRCA1"])
  term_db["GO:0007283"] = ("spermatogenesis", ["BRCA2"])
  term_db["GO:0046982"] = ("protein heterodimerization", ["TP53"])

  let result = @src.cp_enrich_generic(
    gene_list, term_db, background_genes, "Homo sapiens", "GO",
  )
  let top = @src.cp_get_top_terms(result, 2)

  assert_eq(top.length(), 2)
}

///|
test "cluster_profiler_dotplot" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH"]

  let term_db = Map([], capacity=1)
  term_db["GO:0005515"] = ("protein binding", ["TP53", "BRCA1"])

  let result = @src.cp_enrich_generic(
    gene_list, term_db, background_genes, "Homo sapiens", "GO",
  )
  let plot = @src.cp_dotplot(result, 5)

  assert_true(plot.length() > 0)
}

///|
test "cluster_profiler_enrich_map" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH"]

  let term_db = Map([], capacity=1)
  term_db["GO:0005515"] = ("protein binding", ["TP53", "BRCA1"])

  let result = @src.cp_enrich_generic(
    gene_list, term_db, background_genes, "Homo sapiens", "GO",
  )
  let map = @src.cp_enrich_map(result)

  assert_true(map.length() > 0)
}

///|
test "cluster_profiler_cnetplot" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH"]

  let term_db = Map([], capacity=1)
  term_db["GO:0005515"] = ("protein binding", ["TP53", "BRCA1"])

  let result = @src.cp_enrich_generic(
    gene_list, term_db, background_genes, "Homo sapiens", "GO",
  )
  let plot = @src.cp_cnetplot(result, 3)

  assert_true(plot.length() > 0)
}

///|
test "cluster_profiler_simplify" {
  let gene_list = ["TP53", "BRCA1", "BRCA2"]
  let background_genes = ["TP53", "BRCA1", "BRCA2", "GAPDH"]

  let term_db = Map([], capacity=2)
  term_db["GO:0005515"] = ("protein binding", ["TP53", "BRCA1"])
  term_db["GO:0046982"] = ("protein heterodimerization", ["TP53"])

  let result = @src.cp_enrich_generic(
    gene_list, term_db, background_genes, "Homo sapiens", "GO",
  )
  let simplified = @src.cp_simplify(result, 0.5)

  assert_eq(simplified.query_count, 3)
}
