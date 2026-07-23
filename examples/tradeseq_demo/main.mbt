///|
/// tradeSeq Demo - Trajectory Differential Expression Analysis.
/// Inspired by Bioconductor tradeSeq package.

///|
fn main {
  println("=== tradeSeq Trajectory DE Demo ===")
  println("Inspired by Bioconductor tradeSeq package")
  println("")
  
  println("1. Creating trajectory data...")
  let data = @src.create_example_trajectory_data()
  println("   Cells: \{data.n_points.to_string()}")
  println("   Genes: \{data.n_genes.to_string()}")
  println("   Conditions: \{data.conditions.length().to_string()}")
  println("")
  
  println("2. Gene names:")
  let genes = data.genes
  let mut g = 0
  while g < genes.length() {
    println("   Gene \{g.to_string()}: \{genes[g]}")
    g = g + 1
  }
  println("")
  
  println("3. Running tradeSeq analysis...")
  let result = @src.run_tradeseq_analysis(
    data,
    data.genes,
    data.conditions,
    4,
    0.05
  )
  println("   Genes tested: \{result.n_genes.to_string()}")
  println("   Significant genes: \{result.n_significant.to_string()}")
  println("")
  
  println("4. Gene-level results:")
  let mut gr = 0
  while gr < result.gene_results.length() {
    let gene_result = result.gene_results[gr]
    println("   \{gene_result.gene_id}: p-value = \{gene_result.p_value.to_string()}, adj p-value = \{gene_result.adjusted_p_value.to_string()}, significant = \{gene_result.significant.to_string()}")
    gr = gr + 1
  }
  println("")
  
  println("5. Gene expression smoothing (GeneA):")
  let smooth_a = @src.calculate_gene_smooth(data, "GeneA", 4, 10)
  let mut s = 0
  while s < smooth_a.fitted.length() {
    println("   Time point \{s.to_string()}: pseudotime = \{smooth_a.pseudotime[s].to_string()}, fitted = \{smooth_a.fitted[s].to_string()}")
    s = s + 1
  }
  println("")
  
  println("6. tradeSeq summary:")
  let summary = @src.tradeseq_summary(result)
  println("   \{summary}")
  
  println("=== tradeSeq demo completed ===")
}
