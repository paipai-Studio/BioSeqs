///|
/// Demo for pgsea module: Parametric Gene Set Enrichment Analysis.

fn main {
  println("=== pgsea Demo: Parametric Gene Set Enrichment ===")
  println("")

  // Example 1: Score matrix for gene sets
  println("--- Step 1: Score Matrix ---")
  let expr = [
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    [2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
    [5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
    [3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
  ]
  let genes = ["G1", "G2", "G3", "G4"]
  let sets = [("UP", ["G1", "G2"]), ("DOWN", ["G3"]), ("CONST", ["G4"])]
  let (matrix, names) = @src.pgsea_score_matrix(expr, genes, sets)
  println("Score matrix: " + matrix.length().to_string() + " sets x " + matrix[0].length().to_string() + " samples")
  println("Set names: " + names.to_string())
  for i in 0..<names.length() {
    println("  " + names[i] + " scores: " + matrix[i].to_string())
  }
  println("")

  // Example 2: Run PGSEA with two-group comparison
  println("--- Step 2: Run PGSEA ---")
  let labels = ["control", "control", "control", "treat", "treat", "treat"]
  let results = @src.pgsea_run(expr, genes, sets, labels, "control", "treat")
  for r in results {
    println("  " + r.gene_set + ": t=" + r.t_stat.to_string() + ", p=" + r.p_value.to_string())
  }
  println("")

  // Example 3: Top ranking
  println("--- Step 3: Top Ranking ---")
  let top = @src.pgsea_top(results, 2)
  println("Top 2 sets:")
  for r in top {
    println("  " + r.gene_set + " (p=" + r.p_value.to_string() + ")")
  }
  println("")

  // Example 4: Export TSV
  println("--- Step 4: Export TSV ---")
  let tsv = @src.pgsea_export(results)
  println("TSV (first 200 chars): " + tsv.substring(start=0, end=if tsv.length() > 200 { 200 } else { tsv.length() }))
  println("")

  // Example 5: t-distribution p-value
  println("--- Step 5: t Distribution P-value ---")
  let p1 = @src.pgsea_t_pvalue(5.0, 10.0)
  let p2 = @src.pgsea_t_pvalue(0.5, 10.0)
  println("t=5.0, df=10: p=" + p1.to_string())
  println("t=0.5, df=10: p=" + p2.to_string())
  println("")

  println("=== pgsea Demo Complete ===")
}
