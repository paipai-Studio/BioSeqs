///|
/// Phylogenetic Tree Metrics Demo - Advanced phylogenetic analysis.
///
/// Demonstrates tree balance, patristic distances, bipartitions,
/// internal node counting, and Robinson-Foulds distance comparison.
///
/// Run: moon run examples/phylo_metrics_demo/main.mbt
fn format_leaves(leaves : Array[String]) -> String {
  let mut result = "["
  let mut i = 0
  while i < leaves.length() {
    if i > 0 {
      result = result + ", "
    }
    result = result + leaves[i]
    i = i + 1
  }
  result + "]"
}

///|
fn main {
  println("=== Phylogenetic Tree Metrics Demo ===\n")

  // === 1. Create demo trees ===
  println("1. Create Demo Trees:")
  println("=====================")
  let tree1 = @src.create_demo_tree()
  let tree2 = @src.create_demo_tree2()
  println("  Tree 1 (demo): ((A:1, B:1):0.5, (C:1, D:1):0.5)")
  println("  Tree 2 (demo2): ((A:1, C:1):0.5, (B:1, D:1):0.5)")
  let tree1_root = tree1.get_root()
  let tree2_root = tree2.get_root()
  println("  Tree 1 terminals: " + tree1_root.count_terminals().to_string())
  println("  Tree 2 terminals: " + tree2_root.count_terminals().to_string())

  // === 2. Total branch length ===
  println("\n2. Total Branch Length:")
  println("========================")
  let total1 = tree1.total_branch_length()
  let total2 = tree2.total_branch_length()
  println("  Tree 1 total branch length: " + total1.to_string())
  println("  Tree 2 total branch length: " + total2.to_string())

  // === 3. Maximum depth ===
  println("\n3. Maximum Depth (root to deepest leaf):")
  println("=========================================")
  let depth1 = tree1.max_depth()
  let depth2 = tree2.max_depth()
  println("  Tree 1 max depth: " + depth1.to_string())
  println("  Tree 2 max depth: " + depth2.to_string())

  // === 4. Leaf names ===
  println("\n4. Leaf Names:")
  println("==============")
  let leaves1 = tree1_root.get_leaf_names()
  let leaves2 = tree2_root.get_leaf_names()
  println("  Tree 1 leaves: " + format_leaves(leaves1))
  println("  Tree 2 leaves: " + format_leaves(leaves2))

  // === 5. Colless index (tree balance) ===
  println("\n5. Colless Index (tree balance, lower = more balanced):")
  println("========================================================")
  let colless1 = tree1.colless_index()
  let colless2 = tree2.colless_index()
  println("  Tree 1 Colless index: " + colless1.to_string())
  println("  Tree 2 Colless index: " + colless2.to_string())
  if colless1 == 0 {
    println("  Tree 1 is perfectly balanced.")
  }

  // === 6. Patristic distance ===
  println("\n6. Patristic Distance (between pairs of leaves):")
  println("=================================================")
  let pairs = [("A", "B"), ("A", "C"), ("A", "D"), ("B", "C"), ("C", "D")]
  let mut pi = 0
  while pi < pairs.length() {
    let (n1, n2) = pairs[pi]
    let dist = tree1.patristic_distance(n1, n2)
    match dist {
      Some(d) => println("  Tree1 " + n1 + " - " + n2 + ": " + d.to_string())
      None => println("  Tree1 " + n1 + " - " + n2 + ": N/A")
    }
    pi = pi + 1
  }

  // === 7. Patristic distance matrix ===
  println("\n7. Patristic Distance Matrix:")
  println("==============================")
  let matrix = tree1.patristic_distance_matrix()
  println("  Number of pairwise distances: " + matrix.length().to_string())
  let mut mi = 0
  while mi < matrix.length() {
    let (n1, n2, d) = matrix[mi]
    println("  " + n1 + " - " + n2 + ": " + d.to_string())
    mi = mi + 1
  }

  // === 8. Internal node count ===
  println("\n8. Internal Node Count:")
  println("=======================")
  let internal1 = tree1.count_internal()
  let internal2 = tree2.count_internal()
  println("  Tree 1 internal nodes: " + internal1.to_string())
  println("  Tree 2 internal nodes: " + internal2.to_string())

  // === 9. Bipartitions ===
  println("\n9. Bipartitions (tree splits):")
  println("===============================")
  let biparts1 = tree1_root.get_bipartitions()
  let biparts2 = tree2_root.get_bipartitions()
  println("  Tree 1 bipartitions: " + biparts1.length().to_string())
  let mut bi = 0
  while bi < biparts1.length() {
    let (left, right) = biparts1[bi]
    println(
      "  [" +
      bi.to_string() +
      "] " +
      format_leaves(left) +
      " | " +
      format_leaves(right),
    )
    bi = bi + 1
  }
  println("  Tree 2 bipartitions: " + biparts2.length().to_string())
  let mut bi2 = 0
  while bi2 < biparts2.length() {
    let (left, right) = biparts2[bi2]
    println(
      "  [" +
      bi2.to_string() +
      "] " +
      format_leaves(left) +
      " | " +
      format_leaves(right),
    )
    bi2 = bi2 + 1
  }

  // === 10. Robinson-Foulds distance ===
  println("\n10. Robinson-Foulds Distance:")
  println("=============================")
  // Same tree comparison
  let rf_same = @src.robinson_foulds_distance(tree1, tree1)
  println("  Tree1 vs Tree1 (same): " + rf_same.to_string())
  // Different tree comparison
  let rf_diff = @src.robinson_foulds_distance(tree1, tree2)
  println("  Tree1 vs Tree2 (different topology): " + rf_diff.to_string())
  if rf_same == 0 {
    println("  RF=0 confirms identical tree topologies.")
  }
  if rf_diff > 0 {
    println("  RF>0 confirms different tree topologies.")
  }

  println("\n=== Phylogenetic Tree Metrics Demo Complete ===")
}
