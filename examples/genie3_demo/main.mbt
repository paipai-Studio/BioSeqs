///|
/// GENIE3 module example - Gene regulatory network inference.
///
/// Run: moon run examples/genie3_demo
fn main {
  println("=== Bioconductor GENIE3 Module Demo ===")
  println("Gene Regulatory Network Inference via Regression Trees")
  println("")

  // ---- Section 1: Load sample data ----
  println("=== 1. Load Expression Data ===")
  let (expr, gene_names) = @src.genie3_sample_data()
  println("Samples: " + expr.length().to_string())
  println("Genes: " + gene_names.length().to_string())
  println("Gene names: " + gene_names.join(", "))
  println("")

  // ---- Section 2: Build a single regression tree ----
  println("=== 2. Build a Regression Tree ===")
  let features = [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0], [4.0, 8.0], [5.0, 10.0]]
  let target = [2.0, 4.0, 6.0, 8.0, 10.0]
  let tree = @src.genie3_build_tree(features, target, [0, 1], max_depth=3)
  println("Tree nodes: " + tree.nodes.length().to_string())
  let p1 = tree.predict([1.5, 3.0])
  let p2 = tree.predict([4.5, 9.0])
  println("Predict([1.5, 3.0]) = " + p1.to_string())
  println("Predict([4.5, 9.0]) = " + p2.to_string())
  let importance = tree.feature_importance(2)
  println("Feature importance: [" + importance[0].to_string() + ", " + importance[1].to_string() + "]")
  println("")

  // ---- Section 3: Run GENIE3 ----
  println("=== 3. Run GENIE3 (max_depth=5) ===")
  let result = @src.genie3_run(expr, gene_names, max_depth=5)
  println("Trees built: " + result.n_trees.to_string())
  println("Edges found: " + result.edges.length().to_string())
  println("Symmetrized: " + result.symmetrized.to_string())
  println("")

  // ---- Section 4: Top regulatory edges ----
  println("=== 4. Top 10 Regulatory Edges ===")
  let top = @src.genie3_top_edges(result, 10)
  let mut i = 0
  while i < top.length() {
    let e = top[i]
    println("  " + (i + 1).to_string() + ". " + e.regulator + " -> " + e.target + " (weight=" + e.weight.to_string() + ")")
    i = i + 1
  }
  println("")

  // ---- Section 5: Regulators and targets ----
  println("=== 5. Regulators of G2 ===")
  let regs = @src.genie3_regulators_of(result, "G2")
  let mut ri = 0
  while ri < regs.length() {
    println("  " + regs[ri].regulator + " -> " + regs[ri].target + " (weight=" + regs[ri].weight.to_string() + ")")
    ri = ri + 1
  }
  println("")

  println("=== 6. Targets of G1 ===")
  let tgts = @src.genie3_targets_of(result, "G1")
  let mut ti = 0
  while ti < tgts.length() {
    println("  " + tgts[ti].regulator + " -> " + tgts[ti].target + " (weight=" + tgts[ti].weight.to_string() + ")")
    ti = ti + 1
  }
  println("")

  // ---- Section 7: Weight matrix ----
  println("=== 7. Weight Matrix (importance) ===")
  let mut hdr = "          "
  let mut gi = 0
  while gi < gene_names.length() {
    hdr = hdr + gene_names[gi] + "    "
    gi = gi + 1
  }
  println(hdr)
  let mut mi = 0
  while mi < result.matrix.length() {
    let mut row = gene_names[mi] + "  "
    let mut mj = 0
    while mj < result.matrix[mi].length() {
      let v = result.matrix[mi][mj]
      let s = if v < 0.01 { "0.00" } else { v.to_string() }
      row = row + s + "  "
      mj = mj + 1
    }
    println(row)
    mi = mi + 1
  }
  println("")

  // ---- Section 8: Symmetrized network ----
  println("=== 8. Symmetrized Network ===")
  let sym_result = @src.genie3_run(expr, gene_names, max_depth=5, symmetrize=true)
  println("Symmetrized: " + sym_result.symmetrized.to_string())
  println("Edges: " + sym_result.edges.length().to_string())
  println("")

  println("=== Done ===")
}
