///|
/// CellChat module example - Cell-cell communication analysis.
///
/// Run: moon run examples/cellchat_demo
fn main {
  println("=== Bioconductor CellChat Module Demo ===")
  println("Cell-Cell Communication via Ligand-Receptor Interactions")
  println("")

  // ---- Section 1: Load sample data ----
  println("=== 1. Load Single-Cell Expression Data ===")
  let (expr, cell_types, gene_names) = @bio.cellchat_sample_data()
  println("Cells: " + expr.length().to_string())
  println("Genes: " + gene_names.length().to_string())
  println("Cell types: " + cell_types.length().to_string() + " assignments")
  let mut a_count = 0
  let mut b_count = 0
  let mut c_count = 0
  let mut i = 0
  while i < cell_types.length() {
    if cell_types[i] == "TypeA" { a_count = a_count + 1 }
    else if cell_types[i] == "TypeB" { b_count = b_count + 1 }
    else { c_count = c_count + 1 }
    i = i + 1
  }
  println("  TypeA: " + a_count.to_string() + " cells")
  println("  TypeB: " + b_count.to_string() + " cells")
  println("  TypeC: " + c_count.to_string() + " cells")
  println("")

  // ---- Section 2: LR interaction database ----
  println("=== 2. Ligand-Receptor Interaction Database ===")
  let db = @bio.cellchat_lr_database()
  println("Total LR pairs in database: " + db.length().to_string())
  println("First 5 pairs:")
  let mut li = 0
  while li < 5 && li < db.length() {
    println("  " + db[li].ligand + " -> " + db[li].receptor)
    li = li + 1
  }
  println("")

  // ---- Section 3: Run CellChat analysis ----
  println("=== 3. Run CellChat Analysis (10 permutations) ===")
  let result = @bio.cellchat_analyze(expr, cell_types, gene_names, db, n_permutations=10, seed=42)
  println("Total scored interactions: " + result.scores.length().to_string())
  println("Cell types analyzed: " + result.cell_types.length().to_string())
  println("Permutations: " + result.n_permutations.to_string())
  println("")

  // ---- Section 4: Top interactions ----
  println("=== 4. Top 10 Cell-Cell Interactions ===")
  let top = @bio.cellchat_get_top(result, 10)
  let mut ti = 0
  while ti < top.length() {
    let s = top[ti]
    let sig = if s.p_value < 0.05 { "*" } else { "" }
    println("  " + (ti + 1).to_string() + ". " + s.source_celltype + " -> " + s.target_celltype +
      " : " + s.ligand + "(" + s.receptor + ")" +
      " | score=" + s.score.to_string() +
      " | p=" + s.p_value.to_string() + sig)
    ti = ti + 1
  }
  println("")

  // ---- Section 5: Significant interactions ----
  println("=== 5. Significant Interactions (p < 0.05) ===")
  let sig = @bio.cellchat_get_significant(result)
  println("Number of significant interactions: " + sig.length().to_string())
  if sig.length() > 0 {
    let mut si = 0
    while si < sig.length() && si < 5 {
      let s = sig[si]
      println("  " + s.ligand + "->" + s.receptor +
        " : " + s.source_celltype + " -> " + s.target_celltype +
        " (score=" + s.score.to_string() + ")")
      si = si + 1
    }
  }
  println("")

  // ---- Section 6: Aggregate communication scores ----
  println("=== 6. Aggregate Communication Scores (per cell-type pair) ===")
  let agg = @bio.cellchat_aggregate(result)
  // Sort and display
  println("Cell-type pair communication strengths:")
  // Since Map doesn't have iteration order, just print size
  println("  Number of cell-type pairs: " + agg.size().to_string())
  println("")

  // ---- Section 7: Summary ----
  println("=== 7. Analysis Summary ===")
  println(@bio.cellchat_summary(result))
}
