///|
/// Cyclone Demo - Cell cycle scoring using gene pairs.
/// Inspired by Bioconductor cyclone package.

///|
fn main {
  println("=== Cyclone Cell Cycle Scoring Demo ===")
  println("Inspired by Bioconductor cyclone package")
  println("")
  
  println("1. Loading cell cycle gene pairs...")
  let gene_pairs = @bio.cyclone_get_gene_pairs()
  let phases = @bio.cyclone_get_phases(gene_pairs)
  let counts = @bio.cyclone_count_pairs_by_phase(gene_pairs)
  println("   Total gene pairs: \{gene_pairs.length().to_string()}")
  println("   Phases:")
  for phase in phases {
    let count = counts.get(phase).unwrap_or(0)
    println("     - \{phase}: \{count.to_string()} pairs")
  }
  println("")
  
  println("2. Creating synthetic single-cell data...")
  let test_data = @bio.cyclone_create_test_data()
  let (expr_data, cell_ids, gene_names) = test_data
  println("   Cells: \{cell_ids.length().to_string()}")
  println("   Genes: \{gene_names.length().to_string()}")
  println("")
  
  println("3. Running Cyclone cell cycle scoring...")
  let results = @bio.cyclone_score_cells(expr_data, cell_ids, gene_names)
  println("   Scored \{results.length().to_string()} cells")
  println("")
  
  println("4. Cell cycle phase assignments (first 10 cells):")
  let mut i = 0
  while i < 10 && i < results.length() {
    let result = results[i]
    println("   Cell \{result.cell_id}: G1=\{(result.g1_score * 100.0).to_string()}% S=\{(result.s_score * 100.0).to_string()}% G2=\{(result.g2_score * 100.0).to_string()}% M=\{(result.m_score * 100.0).to_string()}% => \{result.assigned_phase}")
    i = i + 1
  }
  println("")
  
  println("5. Phase distribution:")
  let dist = @bio.cyclone_phase_distribution(results)
  let keys = dist.keys()
  for key in keys {
    let proportion = dist.get(key).unwrap_or(0.0)
    println("   \{key}: \{(proportion * 100.0).to_string()}%")
  }
  println("")
  
  println("6. Average scores per phase:")
  let avg = @bio.cyclone_average_scores(results)
  let avg_keys = avg.keys()
  for key in avg_keys {
    let score = avg.get(key).unwrap_or(0.0)
    println("   \{key}: \{score.to_string()}")
  }
  println("")
  
  println("=== Demo completed successfully! ===")
  println("")
  println("Key takeaways:")
  println("- Cyclone uses gene pairs (up-regulated/down-regulated) to score phases")
  println("- For each cell, the proportion of matching gene pairs determines the score")
  println("- The phase with the highest score is assigned as the cell cycle phase")
  println("- Synthetic data was created with known phase-specific expression patterns")
  println("- Cyclone is robust to different expression levels since it uses comparisons")
}