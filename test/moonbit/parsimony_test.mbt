///|
/// Tests for the parsimony scoring module (parsimony.mbt).
/// Tests Fitch and Sankoff parsimony algorithms.

// ============================================================================
// ParsimonyMatrix Tests
// ============================================================================

///|
test "parsimony_matrix_fitch_dna_construction" {
  let m = @src.parsimony_fitch_dna_matrix()
  assert_eq(m.states.length(), 4)
  assert_eq(m.states[0], "A")
  assert_eq(m.states[1], "C")
  assert_eq(m.states[2], "G")
  assert_eq(m.states[3], "T")
}

///|
test "parsimony_matrix_fitch_dna_costs" {
  let m = @src.parsimony_fitch_dna_matrix()
  // Same state = 0 cost
  assert_eq(@src.parsimony_matrix_cost(m, 0, 0), 0.0)
  // Different state = 1 cost
  assert_eq(@src.parsimony_matrix_cost(m, 0, 1), 1.0)
  assert_eq(@src.parsimony_matrix_cost(m, 0, 2), 1.0)
  assert_eq(@src.parsimony_matrix_cost(m, 0, 3), 1.0)
}

///|
test "parsimony_matrix_sankoff_dna_construction" {
  let m = @src.parsimony_sankoff_dna_matrix()
  assert_eq(m.states.length(), 4)
  // A->G (transition) = 1
  assert_eq(@src.parsimony_matrix_cost(m, 0, 2), 1.0)
  // A->C (transversion) = 2
  assert_eq(@src.parsimony_matrix_cost(m, 0, 1), 2.0)
  // A->T (transversion) = 2
  assert_eq(@src.parsimony_matrix_cost(m, 0, 3), 2.0)
}

///|
test "parsimony_matrix_sankoff_custom_costs" {
  let m = @src.parsimony_sankoff_dna_matrix(transition_cost=0.5, transversion_cost=2.5)
  // A->G (transition) = 0.5
  assert_eq(@src.parsimony_matrix_cost(m, 0, 2), 0.5)
  // A->C (transversion) = 2.5
  assert_eq(@src.parsimony_matrix_cost(m, 0, 1), 2.5)
}

///|
test "parsimony_matrix_index_lookup" {
  let m = @src.parsimony_fitch_dna_matrix()
  assert_eq(@src.parsimony_matrix_index(m, "A"), 0)
  assert_eq(@src.parsimony_matrix_index(m, "C"), 1)
  assert_eq(@src.parsimony_matrix_index(m, "G"), 2)
  assert_eq(@src.parsimony_matrix_index(m, "T"), 3)
  // Not found
  assert_eq(@src.parsimony_matrix_index(m, "X"), -1)
}

///|
test "parsimony_matrix_new_custom" {
  let states = ["0", "1"]
  let costs = [[0.0, 1.0], [1.0, 0.0]]
  let m = @src.ParsimonyMatrix::new(states, costs)
  assert_eq(m.states.length(), 2)
  assert_eq(@src.parsimony_matrix_cost(m, 0, 1), 1.0)
  assert_eq(@src.parsimony_matrix_cost(m, 1, 0), 1.0)
  assert_eq(@src.parsimony_matrix_cost(m, 0, 0), 0.0)
}

// ============================================================================
// Fitch Parsimony Tests
// ============================================================================

///|
test "fitch_parsimony_identical_sequences" {
  // All sequences identical => score = 0
  let aln : Map[String, String] = Map([], capacity=4)
  aln["A"] = "AAAA"
  aln["B"] = "AAAA"
  aln["C"] = "AAAA"
  aln["D"] = "AAAA"
  let leaf_a = @src.Clade::new(name=Some("A"))
  let leaf_b = @src.Clade::new(name=Some("B"))
  let leaf_c = @src.Clade::new(name=Some("C"))
  let leaf_d = @src.Clade::new(name=Some("D"))
  let inner_left = @src.Clade::new(clades=[leaf_a, leaf_b])
  let inner_right = @src.Clade::new(clades=[leaf_c, leaf_d])
  let root = @src.Clade::new(clades=[inner_left, inner_right])
  let tree = @src.Tree::new(root, rooted=true)
  let score = @src.fitch_parsimony_score(tree, aln)
  assert_eq(score, 0.0)
}

///|
test "fitch_parsimony_two_groups" {
  // A=B and C=D, but A!=C => 1 change per position
  let aln : Map[String, String] = Map([], capacity=4)
  aln["A"] = "ACGTA"
  aln["B"] = "ACGTA"
  aln["C"] = "ATCTA"
  aln["D"] = "ATCTA"
  let leaf_a = @src.Clade::new(name=Some("A"))
  let leaf_b = @src.Clade::new(name=Some("B"))
  let leaf_c = @src.Clade::new(name=Some("C"))
  let leaf_d = @src.Clade::new(name=Some("D"))
  let inner_left = @src.Clade::new(clades=[leaf_a, leaf_b])
  let inner_right = @src.Clade::new(clades=[leaf_c, leaf_d])
  let root = @src.Clade::new(clades=[inner_left, inner_right])
  let tree = @src.Tree::new(root, rooted=true)
  let score = @src.fitch_parsimony_score(tree, aln)
  // Positions: 0=A/A (same), 1=C/T (diff), 2=G/C (diff), 3=T/T (same), 4=A/A (same)
  // 2 changes expected
  assert_eq(score, 2.0)
}

///|
test "fitch_parsimony_sample_alignment" {
  let aln = @src.parsimony_sample_alignment()
  let tree = @src.parsimony_sample_tree()
  let score = @src.fitch_parsimony_score(tree, aln)
  // A=B=ACGTA, C=D=ATCTA
  // Position 0: A,A,A,A -> 0
  // Position 1: C,C,T,T -> 1
  // Position 2: G,G,C,C -> 1
  // Position 3: T,T,T,T -> 0
  // Position 4: A,A,A,A -> 0
  // Total: 2
  assert_eq(score, 2.0)
}

///|
test "fitch_parsimony_empty_alignment" {
  let aln : Map[String, String] = Map([], capacity=0)
  let tree = @src.parsimony_sample_tree()
  let score = @src.fitch_parsimony_score(tree, aln)
  assert_eq(score, -1.0)
}

///|
test "fitch_parsimony_all_different" {
  let aln : Map[String, String] = Map([], capacity=4)
  aln["A"] = "AAAA"
  aln["B"] = "CCCC"
  aln["C"] = "GGGG"
  aln["D"] = "TTTT"
  let leaf_a = @src.Clade::new(name=Some("A"))
  let leaf_b = @src.Clade::new(name=Some("B"))
  let leaf_c = @src.Clade::new(name=Some("C"))
  let leaf_d = @src.Clade::new(name=Some("D"))
  let inner_left = @src.Clade::new(clades=[leaf_a, leaf_b])
  let inner_right = @src.Clade::new(clades=[leaf_c, leaf_d])
  let root = @src.Clade::new(clades=[inner_left, inner_right])
  let tree = @src.Tree::new(root, rooted=true)
  let score = @src.fitch_parsimony_score(tree, aln)
  // Each position: A,C,G,T -> 2 changes (intersection empty at both internal nodes and root)
  // Total: 4 * 2 = 8... actually let me think again.
  // For each position with A,C,G,T:
  //   inner_left: {A} ∩ {C} = empty -> union {A,C}, cost+1
  //   inner_right: {G} ∩ {T} = empty -> union {G,T}, cost+1
  //   root: {A,C} ∩ {G,T} = empty -> union {A,C,G,T}, cost+1
  //   Total per position: 3
  // 4 positions: 12
  assert_eq(score, 12.0)
}

// ============================================================================
// Sankoff Parsimony Tests
// ============================================================================

///|
test "sankoff_parsimony_identical_sequences" {
  let aln : Map[String, String] = Map([], capacity=4)
  aln["A"] = "AAAA"
  aln["B"] = "AAAA"
  aln["C"] = "AAAA"
  aln["D"] = "AAAA"
  let leaf_a = @src.Clade::new(name=Some("A"))
  let leaf_b = @src.Clade::new(name=Some("B"))
  let leaf_c = @src.Clade::new(name=Some("C"))
  let leaf_d = @src.Clade::new(name=Some("D"))
  let inner_left = @src.Clade::new(clades=[leaf_a, leaf_b])
  let inner_right = @src.Clade::new(clades=[leaf_c, leaf_d])
  let root = @src.Clade::new(clades=[inner_left, inner_right])
  let tree = @src.Tree::new(root, rooted=true)
  let matrix = @src.parsimony_fitch_dna_matrix()
  let score = @src.sankoff_parsimony_score(tree, aln, matrix)
  assert_eq(score, 0.0)
}

///|
test "sankoff_parsimony_with_fitch_matrix" {
  // Sankoff with Fitch matrix should give same results as Fitch
  let aln = @src.parsimony_sample_alignment()
  let tree = @src.parsimony_sample_tree()
  let matrix = @src.parsimony_fitch_dna_matrix()
  let score = @src.sankoff_parsimony_score(tree, aln, matrix)
  let fitch_score = @src.fitch_parsimony_score(tree, aln)
  assert_eq(score, fitch_score)
}

///|
test "sankoff_parsimony_with_transition_matrix" {
  let aln = @src.parsimony_sample_alignment()
  let tree = @src.parsimony_sample_tree()
  let matrix = @src.parsimony_sankoff_dna_matrix()
  let score = @src.sankoff_parsimony_score(tree, aln, matrix)
  // With transition/transversion bias, the score should be > 0
  assert_true(score > 0.0)
}

///|
test "sankoff_parsimony_empty_alignment" {
  let aln : Map[String, String] = Map([], capacity=0)
  let tree = @src.parsimony_sample_tree()
  let matrix = @src.parsimony_fitch_dna_matrix()
  let score = @src.sankoff_parsimony_score(tree, aln, matrix)
  assert_eq(score, -1.0)
}

///|
test "sankoff_parsimony_lower_than_fitch_for_transitions" {
  // When transitions are cheaper, Sankoff should give lower or equal score
  // compared to Fitch (which counts all changes as 1)
  let aln : Map[String, String] = Map([], capacity=4)
  aln["A"] = "AGAG"
  aln["B"] = "AGAG"
  aln["C"] = "GAGA"
  aln["D"] = "GAGA"
  let leaf_a = @src.Clade::new(name=Some("A"))
  let leaf_b = @src.Clade::new(name=Some("B"))
  let leaf_c = @src.Clade::new(name=Some("C"))
  let leaf_d = @src.Clade::new(name=Some("D"))
  let inner_left = @src.Clade::new(clades=[leaf_a, leaf_b])
  let inner_right = @src.Clade::new(clades=[leaf_c, leaf_d])
  let root = @src.Clade::new(clades=[inner_left, inner_right])
  let tree = @src.Tree::new(root, rooted=true)
  let fitch_score = @src.fitch_parsimony_score(tree, aln)
  let matrix = @src.parsimony_sankoff_dna_matrix(transition_cost=1.0, transversion_cost=2.0)
  let sankoff_score = @src.sankoff_parsimony_score(tree, aln, matrix)
  // All changes are A<->G (transitions), so Sankoff = Fitch here
  assert_eq(sankoff_score, fitch_score)
}

// ============================================================================
// parsimony_get_score (unified interface)
// ============================================================================

///|
test "parsimony_get_score_fitch_default" {
  let aln = @src.parsimony_sample_alignment()
  let tree = @src.parsimony_sample_tree()
  let score = @src.parsimony_get_score(tree, aln)
  assert_eq(score, 2.0)
}

///|
test "parsimony_get_score_sankoff_with_matrix" {
  let aln = @src.parsimony_sample_alignment()
  let tree = @src.parsimony_sample_tree()
  let matrix = @src.parsimony_fitch_dna_matrix()
  let score = @src.parsimony_get_score(tree, aln, matrix=Some(matrix))
  assert_eq(score, 2.0)
}

// ============================================================================
// Sample Data Tests
// ============================================================================

///|
test "parsimony_sample_alignment_structure" {
  let aln = @src.parsimony_sample_alignment()
  assert_eq(aln.size(), 4)
  assert_eq(aln["A"], "ACGTA")
  assert_eq(aln["B"], "ACGTA")
  assert_eq(aln["C"], "ATCTA")
  assert_eq(aln["D"], "ATCTA")
}

///|
test "parsimony_sample_tree_structure" {
  let tree = @src.parsimony_sample_tree()
  assert_eq(tree.rooted, true)
  // Root has 2 children
  assert_eq(tree.root.clades.length(), 2)
  // Each internal node has 2 leaves
  assert_eq(tree.root.clades[0].clades.length(), 2)
  assert_eq(tree.root.clades[1].clades.length(), 2)
  // Leaves are terminal
  assert_eq(tree.root.clades[0].clades[0].is_terminal(), true)
  assert_eq(tree.root.clades[0].clades[1].is_terminal(), true)
}

///|
test "parsimony_sample_alignment_varied_structure" {
  let aln = @src.parsimony_sample_alignment_varied()
  assert_eq(aln.size(), 5)
  assert_eq(aln["A"].length(), 8)
}

///|
test "fitch_parsimony_single_position" {
  // Single position, single change
  let aln : Map[String, String] = Map([], capacity=2)
  aln["A"] = "A"
  aln["B"] = "G"
  let leaf_a = @src.Clade::new(name=Some("A"))
  let leaf_b = @src.Clade::new(name=Some("B"))
  let root = @src.Clade::new(clades=[leaf_a, leaf_b])
  let tree = @src.Tree::new(root, rooted=true)
  let score = @src.fitch_parsimony_score(tree, aln)
  assert_eq(score, 1.0)
}

///|
test "sankoff_parsimony_single_position" {
  let aln : Map[String, String] = Map([], capacity=2)
  aln["A"] = "A"
  aln["B"] = "G"
  let leaf_a = @src.Clade::new(name=Some("A"))
  let leaf_b = @src.Clade::new(name=Some("B"))
  let root = @src.Clade::new(clades=[leaf_a, leaf_b])
  let tree = @src.Tree::new(root, rooted=true)
  let matrix = @src.parsimony_fitch_dna_matrix()
  let score = @src.sankoff_parsimony_score(tree, aln, matrix)
  assert_eq(score, 1.0)
}
