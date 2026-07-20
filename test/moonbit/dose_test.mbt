///|
test "dose_enrich_disease" {
  let gene_list = ["TP53", "BRCA1", "BRCA2", "ATM"]
  let background_genes = [
    "TP53", "BRCA1", "BRCA2", "ATM", "GAPDH", "ACTB", "MYC", "EGFR",
  ]

  let disease_db = Map([], capacity=3)
  disease_db["DOID:1234"] = ["TP53", "BRCA1"]
  disease_db["DOID:5678"] = ["BRCA2", "ATM"]
  disease_db["DOID:9012"] = ["GAPDH", "ACTB"]

  let result = @src.enrich_disease(gene_list, disease_db, background_genes)

  assert_true(result.terms.length() >= 1)
}

///|
test "dose_filter_results" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH", "ACTB"]

  let disease_db = Map([], capacity=1)
  disease_db["DOID:1234"] = ["TP53", "BRCA1"]

  let result = @src.enrich_disease(gene_list, disease_db, background_genes)
  let filtered = @src.filter_results(result, 0.05)

  assert_eq(filtered.query_count, 2)
}

///|
test "dose_get_top_terms" {
  let gene_list = ["TP53", "BRCA1", "BRCA2"]
  let background_genes = ["TP53", "BRCA1", "BRCA2", "GAPDH", "ACTB"]

  let disease_db = Map([], capacity=2)
  disease_db["DOID:1234"] = ["TP53", "BRCA1"]
  disease_db["DOID:5678"] = ["BRCA2"]

  let result = @src.enrich_disease(gene_list, disease_db, background_genes)
  let top = @src.get_top_terms(result, 2)

  assert_eq(top.length(), 2)
}

///|
test "dose_summarize_result" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH", "ACTB"]

  let disease_db = Map([], capacity=1)
  disease_db["DOID:1234"] = ["TP53", "BRCA1"]

  let result = @src.enrich_disease(gene_list, disease_db, background_genes)
  let summary = @src.summarize_result(result)

  assert_true(summary.length() > 0)
}

///|
test "dose_empty_gene_list" {
  let gene_list : Array[String] = []
  let background_genes = ["TP53", "BRCA1"]

  let disease_db = Map([], capacity=1)
  disease_db["DOID:1234"] = ["TP53"]

  let result = @src.enrich_disease(gene_list, disease_db, background_genes)

  assert_eq(result.terms.length(), 0)
}
