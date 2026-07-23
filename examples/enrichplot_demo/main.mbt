///|

fn main {
  println("=== Enrichplot Demo - Enrichment Result Visualization ===\n")
  
  let go_terms = [
    @src.EnrichTerm::new("GO:0005623", "cell", 0.0001, 0.001, 2.5, 2.8, 25, ["TP53", "BRCA1", "ATM", "CHEK2", "PTEN"]),
    @src.EnrichTerm::new("GO:0005886", "plasma membrane", 0.0002, 0.002, 2.2, 2.5, 20, ["EGFR", "HER2", "KRAS", "NRAS"]),
    @src.EnrichTerm::new("GO:0003674", "molecular_function", 0.0005, 0.005, 1.8, 2.0, 15, ["TP53", "BRCA1", "EGFR"]),
    @src.EnrichTerm::new("GO:0007049", "cell cycle", 0.001, 0.01, 1.5, 1.7, 18, ["CDK2", "CDK4", "CCND1", "RB1"]),
    @src.EnrichTerm::new("GO:0016020", "membrane", 0.002, 0.02, 1.2, 1.4, 12, ["EGFR", "HER2", "KRAS"]),
    @src.EnrichTerm::new("GO:0006915", "apoptotic process", 0.005, 0.05, 1.0, 1.2, 10, ["TP53", "BAX", "BCL2"]),
    @src.EnrichTerm::new("GO:0007165", "signal transduction", 0.01, 0.1, 0.8, 1.0, 8, ["EGFR", "KRAS", "BRAF"]),
    @src.EnrichTerm::new("GO:0042254", "ribosome biogenesis", 0.02, 0.2, 0.5, 0.7, 6, ["RPS15", "RPL23"]),
  ]
  
  let result = @src.EnrichResult::new(go_terms, "GO_Biological_Process", [])
  
  println("1. Dot Plot (top 5 terms by adjusted p-value):")
  println("----------------------------------------------")
  println(@src.bio_enrichplot_dotplot(result, 5, "GO Enrichment Dot Plot"))
  
  println("\n2. Bar Plot (top 5 terms by enrichment score):")
  println("-----------------------------------------------")
  println(@src.bio_enrichplot_barplot(result, 5, "es", "GO Enrichment Bar Plot"))
  
  println("\n3. Heatmap (top 5 terms, top 10 genes):")
  println("----------------------------------------")
  println(@src.bio_enrichplot_heatmap(result, 5, 10, "GO Enrichment Heatmap"))
  
  println("\n4. Cnet Plot (top 5 terms):")
  println("----------------------------")
  println(@src.bio_enrichplot_cnetplot(result, 5, "GO Enrichment Cnet Plot"))
  
  println("\n5. Enrichment Map (similarity threshold 0.3):")
  println("-----------------------------------------------")
  println(@src.bio_enrichplot_emapplot(result, 0.3, "GO Enrichment Map"))
  
  println("\n=== Enrichplot Demo Complete ===")
}