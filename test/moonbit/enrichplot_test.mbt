///|

test "enrichplot_create_enrich_term" {
  let term = @src.EnrichTerm::new(
    "GO:0005623", "cell", 0.001, 0.01, 1.5, 1.8, 15, 
    ["gene1", "gene2", "gene3"]
  )
  
  assert_eq(term.term_id, "GO:0005623")
  assert_eq(term.term_name, "cell")
  assert_eq(term.pvalue, 0.001)
  assert_eq(term.adjusted_pvalue, 0.01)
  assert_eq(term.enrichment_score, 1.5)
  assert_eq(term.nes, 1.8)
  assert_eq(term.gene_count, 15)
  assert_eq(term.genes.length(), 3)
}

test "enrichplot_dotplot" {
  let terms = [
    @src.EnrichTerm::new("GO:0005623", "cell", 0.001, 0.01, 1.5, 1.8, 15, ["g1", "g2", "g3"]),
    @src.EnrichTerm::new("GO:0005886", "plasma membrane", 0.002, 0.02, 1.2, 1.5, 10, ["g4", "g5"]),
    @src.EnrichTerm::new("GO:0003674", "molecular_function", 0.003, 0.03, 0.8, 1.0, 8, ["g6"]),
  ]
  
  let result = @src.EnrichResult::new(terms, "GO", [])
  let plot = @src.bio_enrichplot_dotplot(result, 10, "Test Dotplot")
  
  assert_true(plot.contains("cell"))
  assert_true(plot.contains("plasma membrane"))
  assert_true(plot.contains("molecular_function"))
}

test "enrichplot_barplot" {
  let terms = [
    @src.EnrichTerm::new("GO:0005623", "cell", 0.001, 0.01, 1.5, 1.8, 15, ["g1"]),
    @src.EnrichTerm::new("GO:0005886", "membrane", 0.002, 0.02, 1.2, 1.5, 10, ["g2"]),
  ]
  
  let result = @src.EnrichResult::new(terms, "GO", [])
  let plot = @src.bio_enrichplot_barplot(result, 5, "padj", "Test Barplot")
  
  assert_true(plot.contains("cell"))
  assert_true(plot.contains("membrane"))
}

test "enrichplot_heatmap" {
  let terms = [
    @src.EnrichTerm::new("GO:0005623", "cell", 0.001, 0.01, 1.5, 1.8, 15, ["g1", "g2", "g3"]),
    @src.EnrichTerm::new("GO:0005886", "membrane", 0.002, 0.02, 1.2, 1.5, 10, ["g2", "g4"]),
  ]
  
  let result = @src.EnrichResult::new(terms, "GO", [])
  let plot = @src.bio_enrichplot_heatmap(result, 5, 5, "Test Heatmap")
  
  assert_true(plot.contains("cell"))
  assert_true(plot.contains("membrane"))
}

test "enrichplot_cnetplot" {
  let terms = [
    @src.EnrichTerm::new("GO:0005623", "cell", 0.001, 0.01, 1.5, 1.8, 15, ["g1", "g2"]),
    @src.EnrichTerm::new("GO:0005886", "membrane", 0.002, 0.02, 1.2, 1.5, 10, ["g2", "g3"]),
  ]
  
  let result = @src.EnrichResult::new(terms, "GO", [])
  let plot = @src.bio_enrichplot_cnetplot(result, 5, "Test Cnetplot")
  
  assert_true(plot.contains("cell"))
  assert_true(plot.contains("membrane"))
  assert_true(plot.contains("g1"))
  assert_true(plot.contains("g2"))
}

test "enrichplot_emapplot" {
  let terms = [
    @src.EnrichTerm::new("GO:0005623", "cell", 0.001, 0.01, 1.5, 1.8, 15, ["g1", "g2", "g3"]),
    @src.EnrichTerm::new("GO:0005886", "membrane", 0.002, 0.02, 1.2, 1.5, 10, ["g2", "g3", "g4"]),
    @src.EnrichTerm::new("GO:0003674", "function", 0.003, 0.03, 0.8, 1.0, 8, ["g5"]),
  ]
  
  let result = @src.EnrichResult::new(terms, "GO", [])
  let plot = @src.bio_enrichplot_emapplot(result, 0.3, "Test Emapplot")
  
  assert_true(plot.contains("cell"))
  assert_true(plot.contains("membrane"))
}

test "enrichplot_empty_result" {
  let result = @src.EnrichResult::new([], "GO", [])
  
  let dotplot = @src.bio_enrichplot_dotplot(result, 10, "Test")
  assert_true(dotplot.contains("Test"))
  
  let barplot = @src.bio_enrichplot_barplot(result, 10, "padj", "Test")
  assert_true(barplot.contains("Test"))
}