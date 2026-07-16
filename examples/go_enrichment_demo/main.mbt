///|
fn main {
  println("=== GO Enrichment Analysis Demo ===")

  let gene_list = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"]

  let go_term1 = @bio.GOTerm::new("GO:0005515", "protein binding", "BP", [
    "GeneA", "GeneB", "GeneD", "GeneF", "GeneG",
  ])
  let go_term2 = @bio.GOTerm::new("GO:0008150", "biological process", "BP", [
    "GeneH", "GeneI", "GeneJ",
  ])
  let go_term3 = @bio.GOTerm::new("GO:0005623", "cell", "CC", [
    "GeneA", "GeneC", "GeneK",
  ])
  let go_term4 = @bio.GOTerm::new("GO:0003674", "molecular_function", "MF", [
    "GeneB", "GeneE", "GeneL", "GeneM",
  ])

  let go_terms = [go_term1, go_term2, go_term3, go_term4]

  let results = @bio.go_enrich(gene_list, go_terms, 100)

  println("\nEnrichment Results:")
  println("===================")

  let mut i = 0
  while i < results.length() {
    let r = results[i]
    println("\nGO ID: " + r.go_id)
    println("Term: " + r.term_name)
    println("Namespace: " + r.namespace)
    println("P-value: " + r.p_value.to_string())
    println("Adjusted P-value: " + r.adjusted_p_value.to_string())
    println("Odds Ratio: " + r.odds_ratio.to_string())
    println("Gene Ratio: " + r.gene_ratio)
    println("Background Ratio: " + r.bg_ratio)
    println("Count: " + r.count.to_string())
    i = i + 1
  }

  let filtered = @bio.filter_enrichment(results, 0.1, 2)

  println("\nFiltered Results (padj < 0.1, count >= 2):")
  println("==========================================")

  i = 0
  while i < filtered.length() {
    let r = filtered[i]
    println(
      r.go_id +
      " - " +
      r.term_name +
      " (padj=" +
      r.adjusted_p_value.to_string() +
      ")",
    )
    i = i + 1
  }

  let bp_terms = @bio.go_terms_by_namespace(go_terms, "BP")
  println("\nBP terms: " + bp_terms.length().to_string())

  println("\n=== Demo Complete ===")
}
