///|
test "reactome_pa_enrich_reactome" {
  let gene_list = ["TP53", "BRCA1", "ATM"]
  let background_genes = ["TP53", "BRCA1", "ATM", "GAPDH", "ACTB", "MYC"]

  let reactome_db = Map([], capacity=3)
  reactome_db["R-HSA-1234"] = (
    "Cell Cycle",
    ["TP53", "BRCA1"],
    "Cellular Processes",
  )
  reactome_db["R-HSA-5678"] = (
    "DNA Repair",
    ["BRCA1", "ATM"],
    "DNA Replication",
  )
  reactome_db["R-HSA-9012"] = (
    "Metabolism",
    ["GAPDH", "ACTB"],
    "Metabolic Processes",
  )

  let result = @src.enrich_reactome(gene_list, reactome_db, background_genes)

  assert_true(result.pathways.length() >= 1)
}

///|
test "reactome_pa_filter_results" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH", "ACTB"]

  let reactome_db = Map([], capacity=1)
  reactome_db["R-HSA-1234"] = (
    "Cell Cycle",
    ["TP53", "BRCA1"],
    "Cellular Processes",
  )

  let result = @src.enrich_reactome(gene_list, reactome_db, background_genes)
  let filtered = @src.reactome_filter_results(result, 0.05)

  assert_eq(filtered.query_count, 2)
}

///|
test "reactome_pa_get_top_pathways" {
  let gene_list = ["TP53", "BRCA1", "ATM"]
  let background_genes = ["TP53", "BRCA1", "ATM", "GAPDH", "ACTB"]

  let reactome_db = Map([], capacity=2)
  reactome_db["R-HSA-1234"] = (
    "Cell Cycle",
    ["TP53", "BRCA1"],
    "Cellular Processes",
  )
  reactome_db["R-HSA-5678"] = ("DNA Repair", ["ATM"], "DNA Replication")

  let result = @src.enrich_reactome(gene_list, reactome_db, background_genes)
  let top = @src.get_top_pathways(result, 2)

  assert_eq(top.length(), 2)
}

///|
test "reactome_pa_group_by_top_level" {
  let gene_list = ["TP53", "BRCA1", "ATM"]
  let background_genes = ["TP53", "BRCA1", "ATM", "GAPDH", "ACTB"]

  let reactome_db = Map([], capacity=2)
  reactome_db["R-HSA-1234"] = (
    "Cell Cycle",
    ["TP53", "BRCA1"],
    "Cellular Processes",
  )
  reactome_db["R-HSA-5678"] = ("DNA Repair", ["ATM"], "DNA Replication")

  let result = @src.enrich_reactome(gene_list, reactome_db, background_genes)
  let groups = @src.group_by_top_level(result)

  assert_true(groups.size() >= 1)
}

///|
test "reactome_pa_summarize_result" {
  let gene_list = ["TP53", "BRCA1"]
  let background_genes = ["TP53", "BRCA1", "GAPDH", "ACTB"]

  let reactome_db = Map([], capacity=1)
  reactome_db["R-HSA-1234"] = (
    "Cell Cycle",
    ["TP53", "BRCA1"],
    "Cellular Processes",
  )

  let result = @src.enrich_reactome(gene_list, reactome_db, background_genes)
  let summary = @src.reactome_summarize_result(result)

  assert_true(summary.length() > 0)
}
