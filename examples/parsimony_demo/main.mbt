///|
/// Parsimony scoring demo - Fitch and Sankoff algorithms.
///
/// Demonstrates phylogenetic parsimony scoring using:
///   - Fitch parsimony (equal cost for all transitions)
///   - Sankoff parsimony (variable cost via cost matrix)
///
/// Inspired by Biopython's Bio.Phylo.TreeConstruction.ParsimonyScorer.
/// Run: moon run examples/parsimony_demo
fn main {
  println("=== Parsimony Module Demo ===")
  println("Fitch and Sankoff parsimony scoring")
  println("Reference: Biopython Bio.Phylo.TreeConstruction.ParsimonyScorer")
  println("")

  // ------------------------------------------------------------------
  // 1. Sample alignment and tree
  // ------------------------------------------------------------------
  println("=== 1. Sample Data ===")
  let aln = @bio.parsimony_sample_alignment()
  println("Alignment (4 taxa, 5 characters):")
  let aln_keys = aln.keys()
  let aln_arr : Array[String] = aln_keys.collect()
  let mut i = 0
  while i < aln_arr.length() {
    let name = aln_arr[i]
    println("  " + name + ": " + aln[name])
    i = i + 1
  }
  println("")
  println("Tree topology: ((A, B), (C, D))")
  println("")

  let tree = @bio.parsimony_sample_tree()

  // ------------------------------------------------------------------
  // 2. Fitch parsimony
  // ------------------------------------------------------------------
  println("=== 2. Fitch Parsimony ===")
  println("Algorithm: equal cost (1) for all state transitions")
  println("")
  let fitch_score = @bio.fitch_parsimony_score(tree, aln)
  println("Fitch parsimony score: " + fitch_score.to_string())
  println("")

  // Explain per-position scores
  println("Per-position analysis:")
  println("  Position 0: A/A/A/A -> 0 changes (all same)")
  println("  Position 1: C/C/T/T -> 1 change (C<->T at root)")
  println("  Position 2: G/G/C/C -> 1 change (G<->C at root)")
  println("  Position 3: G/G/T/T -> 1 change (G<->T at root)")
  println("  Position 4: A/A/A/A -> 0 changes (all same)")
  println("  Total: 3 changes")
  println("")

  // ------------------------------------------------------------------
  // 3. Sankoff parsimony with Fitch matrix
  // ------------------------------------------------------------------
  println("=== 3. Sankoff Parsimony (Fitch matrix) ===")
  let fitch_matrix = @bio.parsimony_fitch_dna_matrix()
  println("Cost matrix (Fitch - all changes cost 1):")
  println("     A    C    G    T")
  println("  A  0.0  1.0  1.0  1.0")
  println("  C  1.0  0.0  1.0  1.0")
  println("  G  1.0  1.0  0.0  1.0")
  println("  T  1.0  1.0  1.0  0.0")
  println("")
  let sankoff_fitch_score = @bio.sankoff_parsimony_score(tree, aln, fitch_matrix)
  println("Sankoff score (Fitch matrix): " + sankoff_fitch_score.to_string())
  println("  (Should equal Fitch score: " + fitch_score.to_string() + ")")
  println("")

  // ------------------------------------------------------------------
  // 4. Sankoff parsimony with transition/transversion bias
  // ------------------------------------------------------------------
  println("=== 4. Sankoff Parsimony (Transition/Transversion bias) ===")
  let sankoff_matrix = @bio.parsimony_sankoff_dna_matrix()
  println("Cost matrix (transitions=1, transversions=2):")
  println("     A    C    G    T")
  println("  A  0.0  2.0  1.0  2.0  (A<->G = transition)")
  println("  C  2.0  0.0  2.0  1.0  (C<->T = transition)")
  println("  G  1.0  2.0  0.0  2.0")
  println("  T  2.0  1.0  2.0  0.0")
  println("")
  let sankoff_biased_score = @bio.sankoff_parsimony_score(tree, aln, sankoff_matrix)
  println("Sankoff score (biased matrix): " + sankoff_biased_score.to_string())
  println("")

  // ------------------------------------------------------------------
  // 5. Comparison with varied alignment
  // ------------------------------------------------------------------
  println("=== 5. Comparison with Varied Alignment ===")
  let aln_varied = @bio.parsimony_sample_alignment_varied()
  println("Varied alignment (5 taxa, 8 characters):")
  let var_keys = aln_varied.keys()
  let var_arr : Array[String] = var_keys.collect()
  let mut v = 0
  while v < var_arr.length() {
    let name = var_arr[v]
    println("  " + name + ": " + aln_varied[name])
    v = v + 1
  }
  println("")

  // Build a tree for 5 taxa: ((A, B), (C, (D, E)))
  let leaf_a = @bio.Clade::new(name=Some("A"))
  let leaf_b = @bio.Clade::new(name=Some("B"))
  let leaf_c = @bio.Clade::new(name=Some("C"))
  let leaf_d = @bio.Clade::new(name=Some("D"))
  let leaf_e = @bio.Clade::new(name=Some("E"))
  let inner_de = @bio.Clade::new(clades=[leaf_d, leaf_e])
  let inner_left = @bio.Clade::new(clades=[leaf_a, leaf_b])
  let inner_right = @bio.Clade::new(clades=[leaf_c, inner_de])
  let root = @bio.Clade::new(clades=[inner_left, inner_right])
  let tree5 = @bio.Tree::new(root, rooted=true)

  let fitch5 = @bio.fitch_parsimony_score(tree5, aln_varied)
  let sankoff5 = @bio.sankoff_parsimony_score(tree5, aln_varied, sankoff_matrix)
  println("Fitch score: " + fitch5.to_string())
  println("Sankoff score (biased): " + sankoff5.to_string())
  println("")

  // ------------------------------------------------------------------
  // 6. Unified interface
  // ------------------------------------------------------------------
  println("=== 6. Unified parsimony_get_score Interface ===")
  println("Without matrix (Fitch): " + @bio.parsimony_get_score(tree, aln).to_string())
  println("With Fitch matrix: " + @bio.parsimony_get_score(tree, aln, matrix=Some(fitch_matrix)).to_string())
  println("With biased matrix: " + @bio.parsimony_get_score(tree, aln, matrix=Some(sankoff_matrix)).to_string())
  println("")

  println("=== Parsimony Demo Complete ===")
}
