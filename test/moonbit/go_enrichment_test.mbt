///|
test "go_enrichment_goterm_new" {
  let genes = ["GeneA", "GeneB", "GeneC"]
  let go_term = @src.GOTerm::new("GO:0005515", "protein binding", "BP", genes)

  assert_eq(go_term.go_id, "GO:0005515")
  assert_eq(go_term.term_name, "protein binding")
  assert_eq(go_term.namespace, "BP")
  assert_eq(go_term.genes.length(), 3)
}

///|
test "go_enrichment_hypergeometric_test" {
  let p = @src.hypergeometric_test(3, 10, 90, 15)
  assert_true(p >= 0.0 && p <= 1.0)
}

///|
test "go_enrichment_fishers_exact_test" {
  let p = @src.fishers_exact_test(5, 10, 3, 15)
  assert_true(p >= 0.0 && p <= 1.0)
}

///|
test "go_enrichment_bonferroni_correction" {
  let p_values = [0.01, 0.02, 0.05]
  let corrected = @src.bonferroni_correction(p_values)

  assert_eq(corrected.length(), 3)
  assert_true(corrected[0] == 0.03)
  assert_true(corrected[1] == 0.06)
}

///|
test "go_enrichment_bh_correction" {
  let p_values = [0.01, 0.02, 0.05, 0.1]
  let corrected = @src.bh_correction(p_values)

  assert_eq(corrected.length(), 4)
  assert_true(corrected[0] <= corrected[1])
  assert_true(corrected[1] <= corrected[2])
}

///|
test "go_enrichment_go_enrich" {
  let genes = ["GeneA", "GeneB", "GeneC"]

  let go_term1 = @src.GOTerm::new("GO:0005515", "protein binding", "BP", [
    "GeneA", "GeneB", "GeneD", "GeneE",
  ])
  let go_term2 = @src.GOTerm::new("GO:0008150", "biological process", "BP", [
    "GeneF", "GeneG",
  ])

  let go_terms = [go_term1, go_term2]

  let results = @src.go_enrich(genes, go_terms, 100)

  assert_eq(results.length(), 1)
  assert_eq(results[0].go_id, "GO:0005515")
  assert_eq(results[0].count, 2)
}

///|
test "go_enrichment_filter_enrichment" {
  let result1 = @src.GOEnrichmentResult::new(
    "GO:0005515",
    "protein binding",
    "BP",
    0.01,
    0.02,
    2.0,
    "2/5",
    "10/100",
    2,
    10,
    ["GeneA", "GeneB"],
  )
  let result2 = @src.GOEnrichmentResult::new(
    "GO:0008150",
    "biological process",
    "BP",
    0.5,
    0.6,
    1.0,
    "1/5",
    "20/100",
    1,
    20,
    ["GeneC"],
  )

  let results = [result1, result2]

  let filtered = @src.filter_enrichment(results, 0.05, 2)

  assert_eq(filtered.length(), 1)
  assert_eq(filtered[0].go_id, "GO:0005515")
}

///|
test "go_enrichment_go_terms_by_namespace" {
  let go_term1 = @src.GOTerm::new("GO:0005515", "protein binding", "BP", [
    "GeneA",
  ])
  let go_term2 = @src.GOTerm::new("GO:0005622", "intracellular", "CC", ["GeneB"])
  let go_term3 = @src.GOTerm::new("GO:0003674", "molecular_function", "MF", [
    "GeneC",
  ])

  let go_terms = [go_term1, go_term2, go_term3]

  let bp_terms = @src.go_terms_by_namespace(go_terms, "BP")

  assert_eq(bp_terms.length(), 1)
  assert_eq(bp_terms[0].namespace, "BP")
}
