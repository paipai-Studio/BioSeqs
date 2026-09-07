///|
/// GOstats Demo - Gene Ontology enrichment analysis using the hypergeometric test.

fn main {
  println("=== GOstats Demo ===")
  println("")

  // 1. Build a GO graph with a few terms and gene mappings
  println("--- Step 1: Build GO Graph ---")

  let go_bp_root = @src.GOStatsTerm::new(
    "GO:0008150",
    "biological_process",
    "BP",
    [],
    [
      "BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "AKT1", "MAPK1", "STAT3",
      "PIK3CA", "PTEN",
    ],
  )

  let go_apoptosis = @src.GOStatsTerm::new(
    "GO:0006915",
    "apoptotic process",
    "BP",
    ["GO:0008150"],
    ["TP53", "PTEN", "AKT1", "PIK3CA", "BRCA1"],
  )

  let go_cell_cycle = @src.GOStatsTerm::new(
    "GO:0007049",
    "cell cycle",
    "BP",
    ["GO:0008150"],
    ["TP53", "BRCA1", "BRCA2", "MYC"],
  )

  let go_signal = @src.GOStatsTerm::new(
    "GO:0023052",
    "signaling",
    "BP",
    ["GO:0008150"],
    ["EGFR", "AKT1", "MAPK1", "STAT3", "PIK3CA"],
  )

  let go_binding = @src.GOStatsTerm::new(
    "GO:0005488",
    "binding",
    "MF",
    [],
    ["BRCA1", "TP53", "EGFR", "STAT3", "MAPK1", "MYC"],
  )

  let graph = @src.GOStatsGraph::new(Map([]))
  let graph = graph.add_term(go_bp_root)
  let graph = graph.add_term(go_apoptosis)
  let graph = graph.add_term(go_cell_cycle)
  let graph = graph.add_term(go_signal)
  let graph = graph.add_term(go_binding)

  println("GO graph built with \{graph.terms.length()} terms")

  // 2. Define universe and selected genes
  println("")
  println("--- Step 2: Define Gene Sets ---")

  let universe = [
    "BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "AKT1", "MAPK1", "STAT3",
    "PIK3CA", "PTEN", "CDK4", "RB1", "APC", "VHL", "NF1", "HRAS", "KRAS",
    "BRAF", "ERBB2", "FGFR1",
  ]

  let selected = ["TP53", "PTEN", "AKT1", "PIK3CA", "BRCA1", "BRCA2", "MYC"]

  println("Universe size: \{universe.length()}")
  println("Selected genes: \{selected.length()}")

  // 3. Create GOHyperGeoTest
  println("")
  println("--- Step 3: Create GOHyperGeoTest ---")

  let hg_test = @src.GOHyperGeoTest::new(universe, selected, graph)
  println("GOHyperGeoTest created")

  // 4. Run hypergeometric test
  println("")
  println("--- Step 4: Hypergeometric Test ---")

  let hyp_results = @src.gostats_hypergeometric_test(hg_test)
  println("Hypergeometric test results: \{hyp_results.length()} GO terms")

  let mut i = 0
  while i < hyp_results.length() {
    let r = hyp_results[i]
    println(
      "  \{r.go_id} (\{r.ontology}) \{r.name}: p=\{r.p_value}, FDR=\{r.fdr}, obs=\{r.observed}, genes=\{r.gene_count}",
    )
    i = i + 1
  }

  // 5. Run conditional test
  println("")
  println("--- Step 5: Conditional Test ---")

  let cond_results = @src.gostats_conditional_test(hg_test)
  println("Conditional test results: \{cond_results.length()} GO terms")

  let mut j = 0
  while j < cond_results.length() {
    let r = cond_results[j]
    println(
      "  \{r.go_id} (\{r.ontology}) \{r.name}: p=\{r.p_value}, FDR=\{r.fdr}, obs=\{r.observed}",
    )
    j = j + 1
  }

  // 6. Filter and print significant GO terms
  println("")
  println("--- Step 6: Significant GO Terms (FDR <= 0.05) ---")

  let sig_results = @src.gostats_filter_results(hyp_results, 0.05, 1)
  if sig_results.length() == 0 {
    println("No significant terms at FDR <= 0.05")
    // Also try a more relaxed threshold
    let relaxed = @src.gostats_filter_results(hyp_results, 0.5, 1)
    println("")
    println("Relaxed threshold (FDR <= 0.5): \{relaxed.length()} terms")
    let mut k = 0
    while k < relaxed.length() {
      let r = relaxed[k]
      println(
        "  \{r.go_id} (\{r.ontology}) \{r.name}: p=\{r.p_value}, FDR=\{r.fdr}, odds=\{r.odds_ratio}",
      )
      k = k + 1
    }
  } else {
    let mut k = 0
    while k < sig_results.length() {
      let r = sig_results[k]
      println(
        "  \{r.go_id} (\{r.ontology}) \{r.name}: p=\{r.p_value}, FDR=\{r.fdr}, odds=\{r.odds_ratio}",
      )
      k = k + 1
    }
  }

  // 7. Print summary
  println("")
  println("--- Summary ---")
  let summary = @src.gostats_summary(hyp_results)
  println(summary)
}
