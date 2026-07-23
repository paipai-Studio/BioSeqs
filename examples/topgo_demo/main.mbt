///|
/// topGO Demo - Demonstrates topology-based GO enrichment analysis.

fn main {
  println("=== topGO Demo ===")
  
  let graph = @bio.create_example_topgo_graph()
  println("GO graph root: \{graph.root}")
  
  let genes_of_interest = ["gene1", "gene2", "gene3"]
  let all_genes = ["gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9", "gene10"]
  
  println("")
  println("Running elim algorithm...")
  let elim_results = @bio.topgo_elim_algorithm(graph, genes_of_interest, all_genes, "BP")
  println("Elim results: \{elim_results.length()} GO terms")
  
  println("")
  println("Running weight01 algorithm...")
  let weight_results = @bio.topgo_weight01_algorithm(graph, genes_of_interest, all_genes, "BP")
  println("Weight01 results: \{weight_results.length()} GO terms")
  
  println("")
  println("Fisher exact test example:")
  let p = @bio.topgo_fisher_exact(3, 10, 5, 100)
  println("  p-value: \{p}")
}
