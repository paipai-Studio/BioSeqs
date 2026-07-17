///|
/// Tree Construction example: build phylogenetic trees from distance matrices.
///
/// Run: moon run examples/tree_construction_demo/main.mbt
fn main {
  println("=== Phylogenetic Tree Construction Demo ===")
  println("")

  // === Example 1: Distance calculation ===
  println("--- 1. Sequence Distance Calculation ---")

  let calc = @bio.DistanceCalculator::new(model="jukes_cantor")

  let seq1 = "ACGTACGTACGTACGTACGTACGTACGT"
  let seq2 = "ACGTACGTACGTACGTACGTACGTGCGA"
  let seq3 = "ACGTACGTACGTAAAGCTTGCATGCCTA"
  let seq4 = "TTTACGTACGTACGTACGTACGTACGTA"

  println("  Seq1: \{seq1}")
  println("  Seq2: \{seq2}")
  println("  Distance (1 vs 2): \{calc.pairwise(seq1, seq2)}")
  println("  Distance (1 vs 3): \{calc.pairwise(seq1, seq3)}")
  println("  Distance (2 vs 3): \{calc.pairwise(seq2, seq3)}")
  println("")

  // === Example 2: Distance matrix from sequences ===
  println("--- 2. Distance Matrix from Sequences ---")

  let sequences = [
    ("seq1", "ACGTACGTACGTACGT"),
    ("seq2", "ACGTACGTACGTAATT"),
    ("seq3", "ACGTACGTAAAAAATT"),
    ("seq4", "TTTTACGTACGTACGT"),
  ]

  let matrix = calc.get_distance(sequences)
  println("  Distance matrix:")
  for i = 0; i < matrix.length(); i = i + 1 {
    let mut row_str = "    "
    for j = 0; j < matrix[i].length(); j = j + 1 {
      row_str = row_str + "\{matrix[i][j]}  "
    }
    println(row_str)
  }
  println("")

  // === Example 3: UPGMA tree ===
  println("--- 3. UPGMA Tree Construction ---")

  let names = ["seq1", "seq2", "seq3", "seq4"]
  let upgma_cons = @bio.TreeConstructor::new(tree_method="upgma")
  let upgma_tree = upgma_cons.build_tree(matrix, names)

  println("  Tree built successfully")
  match upgma_tree.name {
    Some(n) => println("  Root name: \{n}")
    None => println("  Root name: (none)")
  }
  println("  Root clades: \{upgma_tree.clades.length()}")
  println("")

  // === Example 4: NJ tree ===
  println("--- 4. Neighbor-Joining (NJ) Tree Construction ---")

  let nj_cons = @bio.TreeConstructor::new(tree_method="nj")
  let nj_tree = nj_cons.build_tree(matrix, names)

  println("  Tree built successfully")
  println("  Root clades: \{nj_tree.clades.length()}")
  println("")

  // === Example 5: WPGMA tree ===
  println("--- 5. WPGMA Tree Construction ---")

  let wpgma_cons = @bio.TreeConstructor::new(tree_method="wpgma")
  let wpgma_tree = wpgma_cons.build_tree(matrix, names)

  println("  Tree built successfully")
  println("  Root clades: \{wpgma_tree.clades.length()}")
  println("")

  // === Example 6: Kimura 2-parameter model ===
  println("--- 6. Kimura 2-Parameter Distance ---")

  let kimura_calc = @bio.DistanceCalculator::new(model="kimura")
  println("  Kimura distance (1 vs 2): \{kimura_calc.pairwise(seq1, seq2)}")
  println("  Kimura distance (1 vs 3): \{kimura_calc.pairwise(seq1, seq3)}")
  println("")

  println("=== Demo Complete ===")
}
